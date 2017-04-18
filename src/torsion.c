// Copyright 2015 Astex Therapeutics Ltd.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.



#include "pli.h"

static TORSION_TYPE chemscore_torsion_types[] = { { "sp3-sp3 bond", 0, 1, 3, 3, 0.18750, 3.0, 3.1415936, 1 },
						  { "sp2-sp3 bond", 0, 1, 3, 2, 0.09375, 6.0, 0.0000000, 1 },
						  { "sp2-sp2 bond", 0, 1, 2, 2, 0.18750, 2.0, 0.0000000, 1 },
						  { "unknown bond", 0, 1,-1,-1, 0.00000, 0.0, 0.0000000, 1 } };


// new torsion terms to reflect the fact we're using united atoms here
// should make this a lot more sophisticated!

static TORSION_TYPE pliff_torsion_types[] = { { "sp3-sp3 bond", 1, 0, 3, 3, 0.40000, 3.0, 3.1415936, 1 },
					      { "sp2-sp3 bond", 1, 0, 3, 2, 0.13000, 6.0, 0.0000000, 2 },
					      { "sp2-sp2 bond", 1, 0, 2, 2, 0.18000, 2.0, 0.0000000, 1 },
					      { "unknown bond", 1, 0,-1,-1, 0.00000, 0.0, 0.0000000, 1 } };



static void set_bond_torsion(BOND*,RING_LIST*,MOLECULE*);
static int ring_bond(BOND*,RING_LIST*);
static TORSION_TYPE* get_torsion_type(ATOM*,ATOM*);
static void init_molecule_torsions(MOLECULE*);
static void free_molecule_torsions(MOLECULE*);
static void init_torsion(TORSION*);
static void free_torsion(TORSION*);
static int has_other_hybridization(ATOMLIST*, int);



void set_molecule_torsions(MOLECULE *molecule) {

  int i;
  BOND *bond;
  RING_LIST *ring_list;

  init_molecule_torsions(molecule);

  if ((!molecule->selection) || (!molecule->selection->natoms)) {

    return;
  }

  ring_list = atomlist2ringlist(molecule->selection,molecule->selection->natoms);

  for (i=0,bond=molecule->bond;i<molecule->nbonds;i++,bond++) {

    set_bond_torsion(bond,ring_list,molecule);
  }

  free_ring_list(ring_list);
}



TORSION* get_molecule_torsion(MOLECULE *molecule,ATOM *atom1,ATOM *atom2) {

  int i;
  TORSION *torsion;

  for (i=0,torsion=molecule->torsions;i<molecule->n_torsions;i++,torsion++) {

    if (((torsion->atom1 == atom1) && (torsion->atom2 == atom2)) ||
	((torsion->atom2 == atom1) && (torsion->atom1 == atom2))) {

      return(torsion);
    }
  }

  return(NULL);
}



void set_torsion_angle(TORSION *torsion,double angle) {

  rotate_torsion(torsion,angle - torsion->angle);
}



void rotate_torsion(TORSION *torsion,double angle) {

  int i;
  double rm[4][4];
  ATOM **atom;
  ATOMLIST *rotlist;

  rotation_matrix(torsion->atom1->position,torsion->atom2->position,angle,rm);

  rotlist = torsion->rotlist;

  for (i=0,atom=rotlist->atom;i<rotlist->natoms;i++,atom++) {
    
    transform_atom(*atom,rm);
  }

  calc_torsion_angle(torsion);
}



void calc_torsion_angle(TORSION *torsion) {

  torsion->angle = dihedral_angle(torsion->batom1->position,
				  torsion->atom1->position,
				  torsion->atom2->position,
				  torsion->batom2->position);
}



static void set_bond_torsion(BOND *bond,RING_LIST *ring_list,MOLECULE *molecule) {

  int i;
  ATOM *atom1,*atom2,**batomp,*batom,**batoms1,**batoms2;
  ATOMLIST *list1,*list2;
  TORSION *torsion;

  bond->torsion = NULL;

  if ((bond->type != 1) || (ring_bond(bond,ring_list))) {
 
    return;
  }


  atom1 = bond->atom1;
  atom2 = bond->atom2;


  list1 = alloc_atomlist();
  list2 = alloc_atomlist();

  if ((list1 == NULL) || (list2 == NULL)) {

    error_fn("set_bond_torsion: out of memory allocating lists");
  }

  init_atomlist(list1);
  init_atomlist(list2);

  get_connected_atoms(list1,atom1,atom2,0,0,molecule->natoms);
  get_connected_atoms(list2,atom2,atom1,0,0,molecule->natoms);

  if ((list1->natoms < 2) || (list2->natoms < 2) || (!has_other_hybridization(list1, 1)) || (!has_other_hybridization(list2, 1))) {

    free_atomlist(list1);
    free_atomlist(list2);

    return;  
  }

  TORSION_TYPE *torsion_type = get_torsion_type(atom1, atom2);
  if (!torsion_type) {
    warning_fn("Torsion type can not be defined for the bond between atoms %d and %d in molecule %s", atom1->id, atom2->id, atom1->molecule->name);
    free_atomlist(list1);
    free_atomlist(list2);
    return;
  }

  torsion = molecule->torsions + molecule->n_torsions;

  init_torsion(torsion);

  torsion->bond = bond;

  torsion->type = torsion_type;


  if (list1->natoms > list2->natoms) {

    torsion->atom1 = atom1;
    torsion->atom2 = atom2;

    torsion->alist1 = list1;
    torsion->alist2 = list2;

  } else {

    torsion->atom1 = atom2;
    torsion->atom2 = atom1;

    torsion->alist1 = list2;
    torsion->alist2 = list1;
  }

  torsion->rotlist = torsion->alist2;

  batoms1 = torsion->atom1->connections->atom;
  batoms2 = torsion->atom2->connections->atom;

  torsion->batom1 = (batoms1[0] != torsion->atom2) ? batoms1[0] : batoms1[1];
  torsion->batom2 = (batoms2[0] != torsion->atom1) ? batoms2[0] : batoms2[1];

  calc_torsion_angle(torsion);

  bond->torsion = torsion;

  molecule->n_torsions++;
}


static int has_other_hybridization(ATOMLIST *list, int this) {
  ATOM *atom;
  for (int i=0; i<list->natoms; i++) {
    atom = list->atom[i];
    if ((atom->type == NULL) || (atom->type->united_atom == NULL)) {
      return 1;
    }

    if (atom->type->united_atom->hybridisation != this) {
      return 1;
    }
  }
  return 0;
}

static int ring_bond(BOND *bond,RING_LIST *ring_list) {

  int i,j,count;
  RING **ringp,*ring;
  ATOM **atomp,*atom,*atom1,*atom2;

  atom1 = bond->atom1;
  atom2 = bond->atom2;

  for (i=0,ringp=ring_list->rings;i<ring_list->n_rings;i++,ringp++) {

    count = 0;

    ring = *ringp;

    for (j=0,atomp=ring->atom;j<ring->size;j++,atomp++) {
      
      atom = *atomp;

      if ((atom == atom1) || (atom == atom2)) {

	count++;
      }

      if (count == 2) {

	return(1);
      }
    }
  }

  return(0);
}



static TORSION_TYPE* get_torsion_type(ATOM *atom1,ATOM *atom2) {

  int hyb1,hyb2;
  char *scheme;
  ATOM_TYPE *type1,*type2;
  UNITED_ATOM *uatom1,*uatom2;
  TORSION_TYPE *type;

  scheme = (params_get_parameter("torsion_scheme"))->value.s;

  type1 = atom1->type;
  type2 = atom2->type;

  if ((!type1) || (!type2)) {

    return(NULL);
  }

  uatom1 = type1->united_atom;
  uatom2 = type2->united_atom;

  if ((!uatom1) || (!uatom2)) {

    return(NULL);
  }

  hyb1 = uatom1->hybridisation;
  hyb2 = uatom2->hybridisation;

  type = (!strcmp(scheme,"pliff")) ? pliff_torsion_types : chemscore_torsion_types;

  while ((strcmp(type->name,"unknown bond"))) {

    if (type->n_terms != 0) {

      if (((type->hyb1 == hyb1) && (type->hyb2 == hyb2)) ||
	  ((type->hyb1 == hyb2) && (type->hyb2 == hyb1))) {

	return(type);
      }
    }

    type++;
  }

  return(NULL);
}



static void init_molecule_torsions(MOLECULE *molecule) {

  int i;
  ATOM *atom;

  free_molecule_torsions(molecule);

  molecule->torsions = calloc(molecule->nbonds,sizeof(TORSION));

  if (molecule->torsions == NULL) {

    error_fn("init_molecule_torsions: out of memory allocating torsions");
  }

  for (i=0,atom=molecule->atom;i<molecule->natoms;i++,atom++) {

    find_ring(atom,molecule->natoms);
  }
}



static void free_molecule_torsions(MOLECULE *molecule) {

  int i;
  TORSION *torsion;

  if (molecule->torsions) {

    for (i=0,torsion=molecule->torsions;i<molecule->n_torsions;i++,torsion++) {

      free_torsion(torsion);
    }

    free(molecule->torsions);
  }
}



static void init_torsion(TORSION *torsion) {

  torsion->atom1 = NULL;
  torsion->atom2 = NULL;
  torsion->batom1 = NULL;
  torsion->batom2 = NULL;
  torsion->bond = NULL;
  torsion->alist1 = NULL;
  torsion->alist2 = NULL;
  torsion->rotlist = NULL;
  torsion->type = NULL;
}



static void free_torsion(TORSION *torsion) {

   free_atomlist(torsion->alist1);
   free_atomlist(torsion->alist2);
   free_atomlist(torsion->rotlist);
 }

