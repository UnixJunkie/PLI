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


static double calc_atom_simb(ATOM*,SYSTEM*,int);




void normalise_bfactors(SYSTEM *system) {

  int i,natoms;
  double sumb,avgb,sumd2,stdb;
  ATOM *atom;
  MOLECULE *ligand,*protein,*symmetry;

  ligand = system->ligand;
  protein = system->protein;
  symmetry = system->symmetry;

  natoms = 0;
  sumb = 0.0;

  if (ligand != NULL) {

    for (i=0,atom=ligand->atom;i<ligand->natoms;i++,atom++) {

      if ((atom->element != NULL) && (atom->element->id != HYDROGEN)) {

	sumb += atom->bfactor;

	natoms++;
      }
    }
  }

  if (protein != NULL) {

    for (i=0,atom=protein->atom;i<protein->natoms;i++,atom++) {

      if ((atom->element != NULL) && (atom->element->id != HYDROGEN)) {

	sumb += atom->bfactor;

	natoms++;
      }
    }
  }

  if (natoms < 2) {

    warning_fn("normalise_bfactors: not enough atoms");

    return;
  }

  avgb = sumb/natoms;

  sumd2 = 0.0;

  if (ligand != NULL) {

    for (i=0,atom=ligand->atom;i<ligand->natoms;i++,atom++) {

      if ((atom->element != NULL) && (atom->element->id != HYDROGEN)) {

	sumd2 += sqr(atom->bfactor - avgb);
      }
    }
  }

  if (protein != NULL) {

    for (i=0,atom=protein->atom;i<protein->natoms;i++,atom++) {

      if ((atom->element != NULL) && (atom->element->id != HYDROGEN)) {

	sumd2 += sqr(atom->bfactor - avgb);
      }
    }
  }

  stdb = sqrt(sumd2/(natoms-1));

  if (ligand != NULL) {

    for (i=0,atom=ligand->atom;i<ligand->natoms;i++,atom++) {

      atom->bfactor = (stdb > 1E-10) ? (atom->bfactor - avgb)/stdb : 0.0;
    }
  }

  if (protein != NULL) {

    for (i=0,atom=protein->atom;i<protein->natoms;i++,atom++) {

      atom->bfactor = (stdb > 1E-10) ? (atom->bfactor - avgb)/stdb : 0.0;
    }
  }

  if (symmetry != NULL) {

    for (i=0,atom=symmetry->atom;i<symmetry->natoms;i++,atom++) {

      atom->bfactor = (stdb > 1E-10) ? (atom->bfactor - avgb)/stdb : 0.0;
    }
  }
}



void calc_system_simb(SYSTEM *system,int free) {

  int i;
  ATOMLIST *selection;
  ATOM **atomp,*atom;
  SETTINGS *settings;
  FORCE_FIELD *ff;

  selection = system->selection;

  if (selection) {

    for (i=0,atomp=selection->atom;i<selection->natoms;i++,atomp++) {

      atom = *atomp;

      if (free) {

	atom->tdf = calc_atom_simb(atom,system,10);

      } else {

	atom->td = calc_atom_simb(atom,system,10);
      }
    }
  }
}



static double calc_atom_simb(ATOM *atom1,SYSTEM *system,int n_bonds) {

  int i;
  double sumA,sumDA,D,minD,td;
  ATOM *atom2;
  CONTACTLIST *contactlist;
  CONTACT *contact;
  FORCE_FIELD *ff;
  NONBONDED_FF *nonbonded_ff;

  // NEED TO CHANGE THIS - this wont work if contactlist was allocated, but needs to be calculated!

  if (!(atom1->cstatus & ATOM_CONTACTS_CALCULATED)) {

    set_contacts_atom(atom1,system,0);
  }

  contactlist = atom1->contactlist;

  sumA = atom1->exposed_area;
  sumDA = (atom1->exposed_area)*BULK_WATER_TD;

  if ((contactlist) && (contactlist->ncontacts > 0)) {

    ff = system->settings->force_field;

    for (i=0,contact=contactlist->contacts;i<contactlist->ncontacts;i++,contact++) {

      atom2 = contact->atom2;

      if (contact->flags & COVALENT_CONTACT) {

	if (n_bonds > 0) {

	  sumA += contact->area;
	  sumDA += (contact->area)*calc_atom_simb(atom2,system,n_bonds-1);
	}

      } else {

	nonbonded_ff = get_nonbonded_ff(ff->nonbonded,atom1->type,atom2->type);

	minD = nonbonded_ff->Dc + atom1->vdw_radius + atom2->vdw_radius;

	D = (contact->distance < minD) ? 0.0 : contact->distance - minD;
	
	sumA += contact->area;
	sumDA += (contact->area)*D;
      }
    }
  }

  td = (sumA > 1.0E-10) ? sumDA/sumA : 0.0;

  return(td);
}
