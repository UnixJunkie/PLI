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



static void set_molecule_bonds(MOLECULE*);
static void set_molecule_bond_angles(MOLECULE*);
static void set_molecule_bond_angle(MOLECULE*,ATOM*,ATOM*,ATOM*);
static void calc_bond_angle(BOND_ANGLE*);
static void set_molecule_nb_atom_pairs(MOLECULE*);
static void set_molecule_nb_atom_pairs(MOLECULE*);
static void init_molecule_nb_atom_pairs(MOLECULE*);



double molecule_internal_energy(MOLECULE *molecule,ATOM_TYPING_SCHEME *scheme,FORCE_FIELD *ff) {

  int intra_vdw;
  double E;

  intra_vdw = (params_get_parameter("pliff_intra_vdw"))->value.i;

  E = molecule_torsion_energy(molecule) + molecule_vdw_energy(molecule,scheme,ff,intra_vdw);

  return(E);
}



void set_molecule_internal(MOLECULE *molecule) {

  set_molecule_bonds(molecule);
  set_molecule_bond_angles(molecule);
  set_molecule_torsions(molecule);

  set_molecule_nb_atom_pairs(molecule);
}



double molecule_bond_energy(MOLECULE *molecule) {

  int i;
  double E;
  BOND *bond;

  E = 0.0;

  for (i=0,bond=molecule->bond;i<molecule->nbonds;i++,bond++) {

    if (!(bond->flags & SKIP_BOND)) {

      bond->distance = distance(bond->atom1->position,bond->atom2->position);
            
      E += (bond->k)*sqr(bond->distance - bond->Do);
    }
  }

  return(E);
}



double molecule_bond_angle_energy(MOLECULE *molecule) {

  int i;
  double E;
  BOND_ANGLE *bond_angle;

  E = 0.0;

  for (i=0,bond_angle=molecule->bond_angles;i<molecule->n_bond_angles;i++,bond_angle++) {

    if (!(bond_angle->flags & SKIP_BOND_ANGLE)) {

      calc_bond_angle(bond_angle);

      E += (bond_angle->k)*sqr(bond_angle->angle - bond_angle->Ao);
    }
  }

  return(E);
}



double molecule_torsion_energy(MOLECULE *molecule) {

  int i,j,k,l;
  double score,tscore,*pos1,*pos2,*bpos1,*bpos2,X;
  ATOM *atom1,*atom2,**batom1p,**batom2p,*batom1,*batom2;
  BOND *bond;
  ATOMLIST *list1,*list2;
  TORSION *torsions,*torsion;
  TORSION_TYPE *type,*t;

  torsions = molecule->torsions;

  if ((!torsions) || (molecule->n_torsions == 0)) {

    return(0.0);
  }

  score = 0.0;

  for (i=0,torsion=torsions;i<molecule->n_torsions;i++,torsion++) {

    type = torsion->type;

    if (type) {

      calc_torsion_angle(torsion);

      if (type->united) {

	// use only one united torsion angle:

	X = (PI/180.0)*(torsion->angle);

	for (l=0,t=type;l<type->n_terms;l++,t++) {

	  score += (t->A)*(1.0 - cos((t->n)*X - (t->Xo)));
	}

      } else {

	// loop over all 1-4 combinations:

	bond = torsion->bond;

	if (bond) {
	  
	  atom1 = bond->atom1;
	  atom2 = bond->atom2;
	  
	  list1 = atom1->connections;
	  list2 = atom2->connections;
	  
	  if ((list1) && (list2)) {
	    
	    tscore = 0.0;
	    
	    pos1 = atom1->position;
	    pos2 = atom2->position;
	    
	    for (j=0,batom1p=list1->atom;j<list1->natoms;j++,batom1p++) {
	      
	      batom1 = *batom1p;
	      
	      if (batom1 != atom2) {
		
		bpos1 = batom1->position;
		
		for (k=0,batom2p=list2->atom;k<list2->natoms;k++,batom2p++) {
		  
		  batom2 = *batom2p;
		  
		  if (batom2 != atom1) {
		    
		    bpos2 = batom2->position;
		    
		    X = (PI/180.0)*dihedral_angle(bpos1,pos1,pos2,bpos2);
		    
		    for (l=0,t=type;l<type->n_terms;l++,t++) {
		      
		      tscore += (t->A)*(1.0 - cos((t->n)*X - (t->Xo)));
		    }
		  }
		}
	      }
	    }
	    
	    score += tscore;
	  }
	}
      }
    }
  }
	
  return(score);
}



double molecule_vdw_energy(MOLECULE *molecule,ATOM_TYPING_SCHEME *scheme,FORCE_FIELD *ff,int include_attractive) {

  int i,j,**pairs,sid1,sid2;
  double E;
  unsigned int flags;
  ATOM *atom1,*atom2;
  CONTACTLIST *contactlist;
  CONTACT *contact;
  NONBONDED_FF *nonbonded_ff;
  HISTOGRAM *hist_R;

  pairs = molecule->nb_atom_pairs;

  if (pairs == NULL) {

    return(0.0);
  }

  E = 0.0;

  for (i=0,atom1=molecule->atom;i<molecule->natoms;i++,atom1++) {

    sid1 = atom1->seqid;

    contactlist = atom1->contactlist;

    if (contactlist) {

      for (j=0,contact=contactlist->contacts;j<contactlist->ncontacts;j++,contact++) {

	flags = contact->flags;

	if ((!(flags & COVALENT_CONTACT)) && (flags & INTRAMOLECULAR_CONTACT)) {

	  atom2 = contact->atom2;

	  sid2 = atom2->seqid;

	  if ((sid1 < sid2) && (pairs[sid1][sid2])) {

	    calc_contact_geometry(contact);

	    pliff_score_contact(contact,scheme,ff);

	    nonbonded_ff = get_nonbonded_ff(ff->nonbonded,atom1->type,atom2->type);

	    hist_R = get_histogram(nonbonded_ff->R_histogram_id);

	    if (((contact->distance < nonbonded_ff->Dc) && (contact->score > 0.0)) || 
		((include_attractive) && (contact->score < 0.0))) {

	      E += contact->score;
	    }
	  }
	}
      }
    }
  }

  return(E);
}



static void set_molecule_bonds(MOLECULE *molecule) {

  int i;
  BOND *bond;

  for (i=0,bond=molecule->bond;i<molecule->nbonds;i++,bond++) {

    if ((bond->atom1->element->id != HYDROGEN) && (bond->atom2->element->id != HYDROGEN)) {

      bond->distance = distance(bond->atom1->position,bond->atom2->position);
      
      bond->Do = bond->distance;
      bond->k = 300.0;

    } else {

      bond->flags |= SKIP_BOND;
    }
  }
}



static void set_molecule_bond_angles(MOLECULE *molecule) {

  int i,j,k;
  ATOM *atom1,*atom2,*atom3,**atom1p,**atom3p;
  ATOMLIST *conns;

  for (i=0,atom2=molecule->atom;i<molecule->natoms;i++,atom2++) {

    conns = atom2->connections;

    if ((conns) && (conns->natoms > 1)) {

      for (j=0,atom1p=conns->atom;j<conns->natoms-1;j++,atom1p++) {

	atom1 = *atom1p;

	for (k=j+1,atom3p=conns->atom+j+1;k<conns->natoms;k++,atom3p++) {

	  atom3 = *atom3p;
	  
	  if ((atom1->element->id != HYDROGEN) && (atom2->element->id != HYDROGEN) && (atom3->element->id != HYDROGEN)) {

	    set_molecule_bond_angle(molecule,atom1,atom2,atom3);
	  }
	}
      }
    }
  }
}



static void set_molecule_bond_angle(MOLECULE *molecule,ATOM *atom1,ATOM *atom2,ATOM *atom3) {

  BOND_ANGLE *bond_angle;

  if (molecule->bond_angles == NULL) {

    molecule->n_alloc_bond_angles = 100;

    molecule->bond_angles = (BOND_ANGLE*) calloc(molecule->n_alloc_bond_angles,sizeof(BOND_ANGLE));

  } else if (molecule->n_alloc_bond_angles == molecule->n_bond_angles) {

    molecule->n_alloc_bond_angles *= 2;

    molecule->bond_angles = (BOND_ANGLE*) realloc(molecule->bond_angles,(molecule->n_alloc_bond_angles)*sizeof(BOND_ANGLE));
  }

  if (molecule->bond_angles == NULL) {

    error_fn("%s: out of memory allocating bond angles",__func__);
  }

  bond_angle = molecule->bond_angles + molecule->n_bond_angles;

  bond_angle->atom1 = atom1;
  bond_angle->atom2 = atom2;
  bond_angle->atom3 = atom3;

  calc_bond_angle(bond_angle);

  bond_angle->Ao = bond_angle->angle;
  bond_angle->k = 60.0;

  bond_angle->flags = 0;

  molecule->n_bond_angles++;
}



static void calc_bond_angle(BOND_ANGLE *bond_angle) {

  bond_angle->angle = three_point_angle(bond_angle->atom1->position,
					bond_angle->atom2->position,
					bond_angle->atom3->position);
}



static void set_molecule_nb_atom_pairs(MOLECULE *molecule) {

  int i,j,k,**pairs,lid1,lid2;
  int atom1in1,atom1in2,atom2in1,atom2in2;
  ATOM *atom1,*atom2,**latom1p,**latom2p,*latom1,*latom2;
  ATOMLIST *conns1,*conns2,*list1,*list2;
  BOND *bond;
  TORSION *torsions,*torsion;
  TORSION_TYPE *type;

  init_molecule_nb_atom_pairs(molecule);

  // TODO: currently no internal clashes are considered
  //return;

  torsions = molecule->torsions;

  if ((!torsions) || (molecule->n_torsions == 0)) {

    return;
  }

  pairs = molecule->nb_atom_pairs;

  for (i=0,torsion=torsions;i<molecule->n_torsions;i++,torsion++) {

    type = torsion->type;

    bond = torsion->bond;

    atom1 = bond->atom1;
    atom2 = bond->atom2;

    conns1 = atom1->connections;
    conns2 = atom2->connections;

    list1 = torsion->alist1;
    list2 = torsion->alist2;

    if ((list1) && (list2) && (list1->natoms) && (list2->natoms)) {

      for (j=0,latom1p=list1->atom;j<list1->natoms;j++,latom1p++) {

	latom1 = *latom1p;

	if ((latom1 != atom1) && (latom1 != atom2)) {

	  atom1in1 = atom_in_list(conns1,latom1);
	  atom1in2 = atom_in_list(conns2,latom1);

	  lid1 = latom1->seqid;

	  for (k=0,latom2p=list2->atom;k<list2->natoms;k++,latom2p++) {

	    latom2 = *latom2p;

	    if ((latom2 != atom1) && (latom2 != atom2)) {

	      atom2in1 = atom_in_list(conns1,latom2);
	      atom2in2 = atom_in_list(conns2,latom2);

	      lid2 = latom2->seqid;

	      if (type->implicit14) {

		// 1-4 contacts are incorporated in torsion term, so don't include here:

		if (!((atom1in1 && atom2in2) || (atom2in1 && atom1in2))) {
		  
		  pairs[lid1][lid2] = pairs[lid2][lid1] = 1;
		}

	      } else {

		// 1-4 contacts need to be included here:

		pairs[lid1][lid2] = pairs[lid2][lid1] = 1;
	      }

	    }
	  }
	}
      }	
    }
  }
}



static void init_molecule_nb_atom_pairs(MOLECULE *molecule) {

  int i,j,**pairs;

  if (molecule->nb_atom_pairs != NULL) {

    free_2d_imatrix(molecule->nb_atom_pairs);
  }

  molecule->nb_atom_pairs = alloc_2d_imatrix(molecule->natoms,molecule->natoms);

  if (!molecule->nb_atom_pairs) {

    error_fn("init_molecule_nb_atom_pairs: out of memory allocating atom pair matrix");
  }

  pairs = molecule->nb_atom_pairs;

  for (i=0;i<molecule->natoms;i++) {

    for (j=0;j<molecule->natoms;j++) {

      pairs[i][j] = 0;
    }
  }
}
