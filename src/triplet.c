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


void set_hbond_triplets_atom(ATOM *atom,ATOMLIST *list);
static void init_tripletlist(TRIPLETLIST*);
static int triplet_in_list(TRIPLETLIST*,ATOM*,ATOM*,ATOM*);
static void add_triplet_to_list(TRIPLETLIST*,ATOM*,ATOM*,ATOM*,double,double,double);
static void init_triplet(TRIPLET *triplet);
static void set_hbond_flags_system(SYSTEM *system);
static void set_hbond_flags_contact(CONTACT *contact,SETTINGS *settings);

#define PLIFF_HBOND_GEOMETRY_SCORE 7


void set_hbond_triplets_system(SYSTEM *system) {

  int i,j;
  ATOM **atom1p,*atom1,*atom2;
  ATOMLIST *selection;
  CONTACT *contact1;
  CONTACTLIST *clist1;

  selection = system->selection;

  if (selection == NULL) {

    error_fn("set_hbond_triplets_system: no atoms selected");
  }

  set_hbond_flags_system(system);


  for (i=0,atom1p=selection->atom;i<selection->natoms;i++,atom1p++) {

    atom1 = *atom1p;

    set_hbond_triplets_atom(atom1,selection);

    clist1 = atom1->contactlist;

    if (clist1 != NULL) {

      for (j=0,contact1=clist1->contacts;j<clist1->ncontacts;j++,contact1++) {

	atom2 = contact1->atom2;

	set_hbond_triplets_atom(atom2,selection);
      }
    }
  }
}



void set_hbond_triplets_atom(ATOM *atom,ATOMLIST *list) {

  int i,j,t1_ok,t2_ok,t3_ok;
  double v1[4],v2[4],angle;
  CONTACT *contact1,*contact2;
  CONTACTLIST *clist;


  clist = atom->contactlist;

  if (clist == NULL) {

    return;
  }

  t1_ok = (atom_in_list(list,atom)) ? 1 : 0;

  for (i=0,contact1=clist->contacts;i<clist->ncontacts;i++,contact1++) {

    if ((contact1->flags & HBOND_CONTACT) || (contact1->flags & METAL_COORDINATING_CONTACT)) {

      t2_ok = ((t1_ok) || (atom_in_list(list,contact1->atom2))) ? 1 : 0;

      for (j=0,contact2=clist->contacts;j<clist->ncontacts;j++,contact2++) {

	if ((contact1 != contact2) && ((contact2->flags & HBOND_CONTACT) || (contact2->flags & METAL_COORDINATING_CONTACT))) {

	  if (contact1->atom2 != contact2->atom2) {

	    t3_ok = ((t2_ok) || (atom_in_list(list,contact2->atom2))) ? 1 : 0;

	    // alloc tripletlist
	    if (atom->tripletlist == NULL) {

	      atom->tripletlist = (TRIPLETLIST*) malloc(sizeof(TRIPLETLIST));

	      if (atom->tripletlist == NULL) {

	        error_fn("set_atom_hbond_triplets: out of memory allocating triplet_list");
	      }

	      init_tripletlist(atom->tripletlist);
	    }


	    if ((t3_ok) && (!triplet_in_list(atom->tripletlist,atom,contact1->atom2,contact2->atom2))) {

	      calc_vector(atom->position,contact1->atom2->position,v1);
	      calc_vector(atom->position,contact2->atom2->position,v2);

	      angle = vector_angle(v1,v2);

	      add_triplet_to_list(atom->tripletlist,atom,contact1->atom2,contact2->atom2,angle,
		  contact1->scores[PLIFF_HBOND_GEOMETRY_SCORE], contact2->scores[PLIFF_HBOND_GEOMETRY_SCORE]);
	    }
	  }
	}
      }
    }
  }
}


void set_hbond_flags_system(SYSTEM *system) {

  int i,j;
  ATOM **atomp,*atom;
  ATOMLIST *list;
  CONTACT *contact;
  CONTACTLIST *clist;
  SETTINGS *settings;

  list = system->selection;
  settings = system->settings;

  if (list == NULL) {

    error_fn("set_hbond_flags_system: no atoms");
  }

  for (i=0,atomp=list->atom;i<list->natoms;i++,atomp++) {

    atom = *atomp;

    clist = atom->contactlist;

    if (clist != NULL) {

      for (j=0,contact=clist->contacts;j<clist->ncontacts;j++,contact++) {

	set_hbond_flags_contact(contact,settings);
      }
    }
  }
}



void set_hbond_flags_contact(CONTACT *contact,SETTINGS *settings) {

  ATOM *atom1,*atom2;
  double hbond_geom_score;

  if ((contact->flags & COVALENT_CONTACT) || (contact->flags & SECONDARY_CONTACT)) {

    return;
  }

  atom1 = contact->atom1;
  atom2 = contact->atom2;

  if (!connected_atoms(atom1,atom2,NULL,4)) {
    if (hbond_type_match(atom1,atom2)) {
      hbond_geom_score = hbond_geometry_score(contact);
      if (hbond_geom_score > 0.1) {
	contact->flags |= HBOND_CONTACT;
	contact->scores[PLIFF_HBOND_GEOMETRY_SCORE] = hbond_geom_score;
      }
    }
  }

  if (metal_coordinating_contact(contact,0.4)) {

    contact->flags |= METAL_COORDINATING_CONTACT;
    // -1.0 to distinguish from hbond contact where the value will be between 0 and 1
    contact->scores[PLIFF_HBOND_GEOMETRY_SCORE] = -1.0;
  }
}


static void init_triplet(TRIPLET *triplet) {

  triplet->atom1 = NULL;
  triplet->atom2 = NULL;
  triplet->atom3 = NULL;

  triplet->angle = 0.0;
}


static void init_tripletlist(TRIPLETLIST *list) {

  if (list == NULL) {

    error_fn("init_tripletlis: list undefined");
  }

  list->ntriplets = 0;
  list->n_alloc_triplets = 0;
  list->triplets = NULL;
}



static int triplet_in_list(TRIPLETLIST *list,ATOM *atom1,ATOM *atom2,ATOM *atom3) {

  int i;
  TRIPLET *triplet;

  if (list == NULL) {

    return(0);
  }

  for (i=0,triplet=list->triplets;i<list->ntriplets;i++,triplet++) {

    if ((atom1 == triplet->atom1) &&
	(((atom2 == triplet->atom2) && (atom3 == triplet->atom3)) ||
	 ((atom3 == triplet->atom2) && (atom2 == triplet->atom3)))) {

      return(1);
    }
  }

  return(0);
}



static void add_triplet_to_list(TRIPLETLIST *list,ATOM *atom1,ATOM *atom2,ATOM *atom3,double angle, double hbond_geom_score1, double hbond_geom_score2) {

  TRIPLET *triplet;

  if (list == NULL) {

    error_fn("add_triplet_to_list: list undefined");

    return;
  }

  if (list->n_alloc_triplets == 0) {

    list->n_alloc_triplets = 2;

    list->triplets = (TRIPLET*) calloc(list->n_alloc_triplets,sizeof(TRIPLET));

  } else if (list->ntriplets == list->n_alloc_triplets) {

    list->n_alloc_triplets *= 2;

    list->triplets = (TRIPLET*) realloc(list->triplets,(list->n_alloc_triplets)*sizeof(TRIPLET));
  }

  if (list->triplets == NULL) {

    error_fn("add_triplet_to_list: out of memory allocating list");
  }

  triplet = list->triplets + list->ntriplets;

  init_triplet(triplet);

  triplet->atom1 = atom1;
  triplet->atom2 = atom2;
  triplet->atom3 = atom3;
  triplet->hbond_geom_score1 = hbond_geom_score1;
  triplet->hbond_geom_score2 = hbond_geom_score2;


  triplet->angle = angle;

  list->ntriplets++;
}


void write_hbond_triplets_system(PLI_FILE *file,SYSTEM *system, enum OUTPUT_FORMAT oformat) {

  int i;
  ATOM **atomp;
  ATOMLIST *list;

  list = system->selection;

  for (i=0,atomp=list->atom;i<list->natoms;i++,atomp++) {

    write_hbond_triplets_atom(file,*atomp, oformat);
  }
}



void write_hbond_triplets_atom(PLI_FILE *file,ATOM *atom, enum OUTPUT_FORMAT oformat) {

  int i;
  TRIPLET *triplet;
  TRIPLETLIST *list;

  list = atom->tripletlist;

  if (list == NULL) {

    return;
  }

  for (i=0,triplet=list->triplets;i<list->ntriplets;i++,triplet++) {

    write_hbond_triplet(file,triplet, oformat);
  }
}



void write_hbond_triplet(PLI_FILE *file,TRIPLET *triplet, enum OUTPUT_FORMAT oformat) {

  int sym_atom1, sym_atom2, sym_atom3, json;
  ATOM *atom1, *atom2, *atom3;

  atom1 = triplet->atom1;
  atom2 = triplet->atom2;
  atom3 = triplet->atom3;

  sym_atom1 = (atom1->molecule->flags & SYMMETRY_MOLECULE) ? 1 : 0;
  sym_atom2 = (atom2->molecule->flags & SYMMETRY_MOLECULE) ? 1 : 0;
  sym_atom3 = (atom3->molecule->flags & SYMMETRY_MOLECULE) ? 1 : 0;

  json = (oformat == JSON) ? 1 : 0;
  write_line(file,(json) ?
	   "{\"triplet\":{\"aid1\":%d,\"aid2\":%d,\"aid3\":%d,\"geom_score1\":%.3f,\"geom_score2\":%.3f,\"angle\":%.4f,\"sym_flag1\":%d,\"sym_flag2\":%d,\"sym_flag3\":%d}}\n":
	   "TRIPLET %6d %6d %6d %10.3f %10.3f %10.4f\n",
	   triplet->atom1->id, triplet->atom2->id, triplet->atom3->id,
	   triplet->hbond_geom_score1, triplet->hbond_geom_score2,
	   triplet->angle, sym_atom1, sym_atom2, sym_atom3);
}

