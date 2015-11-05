// Copyright 2015 Astex Therapautics Ltd.
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



static double molecule_torsion_energy(MOLECULE*);
static double molecule_clash_energy(MOLECULE*,ATOM_TYPING_SCHEME*,FORCE_FIELD*);
static void set_molecule_bond_distances(MOLECULE*);
static void set_molecule_bond_angles(MOLECULE*);
static void set_molecule_nb_atom_pairs(MOLECULE*);
static void set_molecule_nb_atom_pairs(MOLECULE*);
static void init_molecule_nb_atom_pairs(MOLECULE*);



double molecule_internal_energy(MOLECULE *molecule,ATOM_TYPING_SCHEME *scheme,FORCE_FIELD *ff) {

  double E;

  E = molecule_torsion_energy(molecule) + molecule_clash_energy(molecule,scheme,ff);

  return(E);
}



void set_molecule_internal(MOLECULE *molecule) {

  set_molecule_bond_distances(molecule);
  set_molecule_bond_angles(molecule);
  set_molecule_torsions(molecule);

  set_molecule_nb_atom_pairs(molecule);
}



static double molecule_torsion_energy(MOLECULE *molecule) {

  int i,j,k;
  double score,tscore,*pos1,*pos2,*bpos1,*bpos2,X;
  ATOM *atom1,*atom2,**batom1p,**batom2p,*batom1,*batom2;
  BOND *bond;
  ATOMLIST *list1,*list2;
  TORSION *torsions,*torsion;
  TORSION_TYPE *type;

  torsions = molecule->torsions;

  if ((!torsions) || (molecule->n_torsions == 0)) {

    return(0.0);
  }

  score = 0.0;

  for (i=0,torsion=torsions;i<molecule->n_torsions;i++,torsion++) {

    type = torsion->type;

    if (type) {

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

		  X = (PI/180.0)*torsion_angle(bpos1,pos1,pos2,bpos2);

		  tscore += (type->A)*(1.0 - cos((type->n)*X - (type->Xo)));
		}
	      }
	    }
	  }

	  score += tscore;
	}
      }
    }
  }

  return(score);
}



static double molecule_clash_energy(MOLECULE *molecule,ATOM_TYPING_SCHEME *scheme,FORCE_FIELD *ff) {

  int i,j,**pairs,sid1,sid2;
  double E;
  unsigned int flags;
  ATOM *atom1,*atom2;
  CONTACTLIST *contactlist;
  CONTACT *contact;
  NONBONDED_FF *nonbonded_ff;
  HISTOGRAM *hist_R;

  pairs = molecule->nb_atom_pairs;

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

	    if ((contact->distance < nonbonded_ff->Dc) && (contact->score > 0.0)) {

	      E += contact->score;
	    }
	  }
	}
      }
    }
  }

  return(E);
}



static void set_molecule_bond_distances(MOLECULE *molecule) {

}



static void set_molecule_bond_angles(MOLECULE *molecule) {
}



static void set_molecule_nb_atom_pairs(MOLECULE *molecule) {

  int i,j,k,**pairs,lid1,lid2;
  ATOM *atom1,*atom2,**latom1p,**latom2p,*latom1,*latom2;
  ATOMLIST *conns1,*conns2,*list1,*list2;
  BOND *bond;
  TORSION *torsions,*torsion;

  init_molecule_nb_atom_pairs(molecule);

  return;

  torsions = molecule->torsions;

  if ((!torsions) || (molecule->n_torsions == 0)) {

    return;
  }

  pairs = molecule->nb_atom_pairs;

  for (i=0,torsion=torsions;i<molecule->n_torsions;i++,torsion++) {

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

	if ((latom1 != atom1) && (latom1 != atom2) && (!atom_in_list(conns1,latom1)) && (!atom_in_list(conns2,latom1))) {

	  lid1 = latom1->seqid;

	  for (k=0,latom2p=list2->atom;k<list2->natoms;k++,latom2p++) {

	    latom2 = *latom2p;

	    if ((latom2 != atom1) && (latom2 != atom2) && (!atom_in_list(conns1,latom2)) && (!atom_in_list(conns2,latom2))) {

	      lid2 = latom2->seqid;

	      pairs[lid1][lid2] = pairs[lid2][lid1] = 1;
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
