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



#define LOOSE_WATER_SCORE_CUTOFF -0.2
#define MAX_WATER_ENERGY 1.0



static int flag_all_waters(SYSTEM*);
static int flag_loose_waters(SYSTEM*);
static int flag_non_mediating_waters(SYSTEM*);
static ATOM* next_loose_water(ATOMLIST*);
static void add_waters_to_selection(SYSTEM*);
static void remove_waters_from_selection(SYSTEM*);
static ATOMLIST* system2waters(SYSTEM*);


int flag_skipped_waters(SYSTEM *system) {

  int n_flagged;
  SETTINGS *settings;

  settings = system->settings;

  if (!strcmp(settings->keep_waters,"all")) {

    return(0);

  } else if (!strcmp(settings->keep_waters,"none")) {

    n_flagged = flag_all_waters(system);

  } else if (!strcmp(settings->keep_waters,"tight")) {

    n_flagged = flag_loose_waters(system);

  } else if (!strcmp(settings->keep_waters, "mediating")) {

    n_flagged = flag_non_mediating_waters(system);

  } else {

    error_fn("flag_skipped_waters: unknown option '%s'",settings->keep_waters);
  }

  return(n_flagged);
}


static int flag_non_mediating_waters(SYSTEM *system) {
  MOLECULE *ligand = system->ligand;
  MOLECULE *protein = system->protein;

  if (!ligand || !protein) {
    error_fn("Mediating waters can't be identified without both protein and ligand being specified");
    return(0);
  }

  set_contacts_system(system, 0);

  MOLECULE *molecule;
  ATOM *atom;
  int n_skip = 0;
  for (int i=0; i<system->molecule_list->n_molecules; i++) {
    molecule = system->molecule_list->molecules[i];
    for (int j=0; j<molecule->natoms; j++) {
      atom = molecule->atom + j;
      if (atom->flags & WATER_OXYGEN) {
        int protein_contact = 0;
        int ligand_contact = 0;
        atom->flags |= SKIP_ATOM;
        n_skip++;
        for (int k=0; k<atom->contactlist->ncontacts; k++) {
          if ((atom->contactlist->contacts+k)->atom2->molecule == ligand) {
            ligand_contact = 1;
          }
          if ((atom->contactlist->contacts+k)->atom2->molecule == protein) {
            protein_contact = 1;
          }
          if (protein_contact && ligand_contact) {
            atom->flags ^= SKIP_ATOM;
            n_skip--;
            break;
          }
        }
      }
    }
  }
  reset_system(system, 0);
  return(n_skip);
}


static int flag_all_waters(SYSTEM *system) {

  int i,j,n_flagged;
  MOLECULE **moleculep,*molecule;
  ATOM *atom;
  MOLECULE_LIST *list;

  list = system->molecule_list;

  if (list == NULL) {

    return(0);
  }

  n_flagged = 0;

  for (i=0,moleculep=list->molecules;i<list->n_molecules;i++,moleculep++) {

    molecule = *moleculep;

    for (j=0,atom=molecule->atom;j<molecule->natoms;j++,atom++) {

      if (atom->flags & WATER_OXYGEN) {

	atom->flags |= SKIP_ATOM;

	n_flagged++;
      }
    }
  }

  return(n_flagged);
}



static int flag_loose_waters(SYSTEM *system) {

  int i,j,n_flagged;
  ATOM *water;
  ATOMLIST *selection,*waters;
  PLI_SFUNC *sfunc;

  n_flagged = 0;

  selection = system->selection;

  waters = system2waters(system);

  //system->selection = waters;
  
  //reset_system(system,1);

  //minimise(system,10);

  sfunc = get_sfunc("pliff");

  if (waters) {

    selection = system->selection;
  
    reset_system(system,1);
    
    system->selection = waters;

    do {

      score_system(system,sfunc);

      water = next_loose_water(waters);

      if (water) {

	remove_atom_from_list(waters,water);
	remove_atom_from_list(selection,water);

	water->flags |= SKIP_ATOM;
	water->status = ATOM_CHANGED;
	water->cstatus = ATOM_NOTHING_CALCULATED;

	n_flagged++;
      }

      //reset_system(system,0);

    } while (water);

    system->selection = selection;

    reset_system(system,0);

    free_atomlist(waters);
  }

  printf("n_flagged = %d\n",n_flagged);

  return(n_flagged);
}



void flag_active_waters(SYSTEM *system) {

  int i,j;
  ATOM *water,*atom;
  MOLECULE *ligand,*protein;

  ligand = system->ligand;
  protein = system->protein;

  if ((ligand) && (protein) && (ligand->n_active_waters)) {

    for (i=0,water=ligand->active_waters;i<ligand->n_active_waters;i++,water++) {

      if (water->flags & SKIP_ATOM) {

	for (j=0,atom=protein->atom;j<protein->natoms;j++,atom++) {

	  if (atom->id == water->id) {

	    if (atom->type == water->type) {

	      atom->flags |= SKIP_ATOM;

	    } else {
	      warning_fn("flag_active_waters: atom %d is not a water molecule");
	    }
	  }
	}
      }
    }
  }
}



static ATOM* next_loose_water(ATOMLIST *waters) {

  int i,j,set;
  double E,Ep,Ew,maxE;
  ATOM **waterp,*water,*loose_water;
  CONTACT *contact;
  CONTACTLIST *contactlist;

  loose_water = NULL;

  maxE = -0.70;
 
  for (i=0,waterp=waters->atom;i<waters->natoms;i++,waterp++) {

    water = *waterp;

    if (!(water->flags & SKIP_ATOM)) {

      Ep = Ew = 0.0;

      contactlist = water->contactlist;

      if (contactlist) {

	for (j=0,contact=contactlist->contacts;j<contactlist->ncontacts;j++,contact++) {

	  set = contact->atom2->set;

	  if (set == WATER_ATOM_SET) {

	    Ew += contact->score;

	  } else if (set != LIGAND_ATOM_SET) {

	    Ep += contact->score;
	  }
	}
      }

      E = Ep + Ew;

      if (E > maxE) {

	maxE = E;

	loose_water = water;

      } else if ((Ep > maxE) && (Ep >= Ew)) {

	//maxE = Ep;

	//loose_water = water;
      }
    }
  }

  if (loose_water) {
    
    //printf("water %5d %5d %10.4lf\n",loose_water->id,loose_water->subid,maxE);
  }

  return(loose_water);
}



static ATOMLIST* system2waters(SYSTEM *system) {

  int i,j;
  ATOMLIST *list;
  ATOM *atom;
  MOLECULE **moleculep,*molecule;

  if (system->molecule_list == NULL) {

    return(NULL);
  }

  list = alloc_atomlist();

  if (list == NULL) {

    error_fn("system2waters: out of memory allocating list");
  }

  init_atomlist(list);

  for (i=0,moleculep=system->molecule_list->molecules;i<system->molecule_list->n_molecules;i++,moleculep++) {

    molecule = *moleculep;

    for (j=0,atom=molecule->atom;j<molecule->natoms;j++,atom++) {

      if (atom->flags & WATER_OXYGEN) {

	add_atom_to_list(list,atom,0);
      }
    }
  }

  return(list);
}



static void add_waters_to_selection(SYSTEM *system) {

  int i,j;
  ATOMLIST *list,*selection;
  ATOM **atom1p,*atom1,*atom2;
  CONTACTLIST *contactlist;
  CONTACT *contact;

  set_contacts_system(system,0);

  selection = system->selection;

  list = alloc_atomlist();

  if (list == NULL) {

    error_fn("add_waters_to_selection: out of memory allocating list");
  }

  init_atomlist(list);

  for (i=0,atom1p=selection->atom;i<selection->natoms;i++,atom1p++) {

    atom1 = *atom1p;

    contactlist = atom1->contactlist;

    if (contactlist) {

      for (j=0,contact=contactlist->contacts;j<contactlist->ncontacts;j++,contact++) {

	atom2 = contact->atom2;

	if ((atom2->flags & WATER_OXYGEN) && (!atom_in_list(list,atom2))) {

	  add_atom_to_list(list,atom2,0);
	}
      }
    }
  }

  for (i=0,atom1p=list->atom;i<list->natoms;i++,atom1p++) {

    atom1 = *atom1p;

    if (!atom_in_list(selection,atom1)) {

      add_atom_to_list(selection,atom1,0);
    }
  }

  free_atomlist(list);
}



static void remove_waters_from_selection(SYSTEM *system) {

  int i,n_waters;
  ATOMLIST *selection;
  ATOM **atomp,*atom;

  selection = system->selection;

  n_waters = 0;

  for (i=0,atomp=selection->atom;i<selection->natoms;i++,atomp++) {

    atom = *atomp;

    if (atom->flags & WATER_OXYGEN) {

      n_waters++;
    }
  }

  selection->natoms -= n_waters;
}

