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



static struct ContactSettings {
  int voronoi_contacts;
  int exact_voronoi_areas;
  int allow_bad_clashes;
  int allow_covalent_bonds;
} *contact_settings = NULL;



static void init_contactlist(CONTACTLIST*);
static void add_contact_to_list(CONTACTLIST*,ATOM*,ATOM*,double);
static CONTACTLIST* alloc_contactlist(CONTACTLIST*);
static void flag_bond_mediated_contacts(ATOM*,ATOM*,ATOMLIST*);
static int bond_mediated_contact(CONTACT*,ATOMLIST*,ATOMLIST*);
static int secondary_contact(CONTACT*,CONTACTLIST*);
static void swap_contacts(CONTACT*,CONTACT*);
static int allow_contact(ATOM*,ATOM*);
static void set_contacts_profile_atom(ATOM*);
static void write_contact_profile_list(FILE*,ATOMLIST*);
static void molecule2contacts(ATOM*,MOLECULE*);
static void map2contacts(ATOM*,MAP*);
static void pair2contacts(ATOM*,ATOM*);
static void flag_atom_list_secondary_contacts(ATOMLIST*);
static double sum_sas_atomlist(ATOMLIST*,ATOM_TYPE*,int);


void run_contacts(SETTINGS *settings) {

  int flag;
  SYSTEM *system;

  system = settings2system(settings);

  if ((settings->oflags & OUTPUT_SCORES) || (settings->oflags & OUTPUT_HBONDS)) {

    prep_score(system,settings->sfunc);

    score_system(system,settings->sfunc);

  } else {

    prep_system(system,ANY_MOLECULE);

    set_contacts_system(system,1);
  }

  flag_atom_list_secondary_contacts(system->selection);

  flag = flag_clashing_contacts(system->selection);

  if ((flag == 2) && (!(contact_settings->allow_bad_clashes))) {

    error_fn("contacts: severely clashing contacts observed");
  }

  if ((flag == 1) && (!(contact_settings->allow_covalent_bonds))) {

    error_fn("contacts: covalent bonds observed");
  }

  normalise_bfactors(system);

  output_system(PLI_STDOUT,system);

  unprep_system(system,ANY_MOLECULE);

  free_system(system);
}



void init_contact_settings(void) {

  if (contact_settings == NULL) {

    contact_settings = (struct ContactSettings*) malloc(sizeof(struct ContactSettings));
  }

  if (contact_settings == NULL) {

    error_fn("init_contact_settings: out of memory allocating settings");
  }

  contact_settings->voronoi_contacts = (params_get_parameter("voronoi_contacts"))->value.i;
  contact_settings->exact_voronoi_areas = (params_get_parameter("exact_voronoi_areas"))->value.i;
  contact_settings->allow_bad_clashes = (params_get_parameter("allow_bad_clashes"))->value.i;
  contact_settings->allow_covalent_bonds = (params_get_parameter("allow_covalent_bonds"))->value.i;
}



void set_contacts_system(SYSTEM *system,unsigned int calc_geometries) {

  int i,j;
  unsigned int flag;
  ATOM **atomp,*atom,*catom;
  ATOMLIST *selection,*updated;
  CONTACT *contact,*icontact;
  CONTACTLIST *contactlist;

  selection = system->selection;

  if (selection == NULL) {

    error_fn("set_contacts_system: no atoms selected");
  }

  // find the list of atoms that need their contacts recalculated
  // moved or changed atoms have their contacts recalculated in this first step

  updated = (ATOMLIST*) alloc_atomlist();

  init_atomlist(updated);

  for (i=0,atomp=selection->atom;i<selection->natoms;i++,atomp++) {

    atom = *atomp;

    if (atom->status & ATOM_MOVED_OR_CHANGED) {

      if (!atom_in_list(updated,atom)) {

	// add atom that has moved or changed:

	add_atom_to_list(updated,atom,1);

	// add atoms that were previously contacting moved or changed atom:

	contactlist = atom->contactlist;

	if (contactlist != NULL) {

	  for (j=0,contact=contactlist->contacts;j<contactlist->ncontacts;j++,contact++) {

	    catom = contact->atom2;

	    flag = (atom_in_list(selection,catom)) ? 2 : 3;

	    if ((flag == 3) || (!(catom->status & ATOM_MOVED_OR_CHANGED))) {

	      if (!atom_in_list(updated,catom)) {
	      
		add_atom_to_list(updated,catom,flag);
	      }
	    }
	  }
	}

	// set new contacts for moved or changed atom:

	set_contacts_atom(atom,system,calc_geometries);

	// add atoms that are now contacting moved or changed atom:

	contactlist = atom->contactlist;

	if (contactlist != NULL) {

	  for (j=0,contact=contactlist->contacts;j<contactlist->ncontacts;j++,contact++) {

	    catom = contact->atom2;
	    flag = (atom_in_list(selection,catom)) ? 2 : 3;

	    if ((flag == 3) || (!(catom->status & ATOM_MOVED_OR_CHANGED))) {

	      if (!atom_in_list(updated,catom)) {  

		add_atom_to_list(updated,catom,flag);
	      }
	    }
	  }
	}
      }
    }
  }

  // recalculate contacts for updated atoms:

  for (i=0,atomp=updated->atom;i<updated->natoms;i++,atomp++) {

    if ((updated->flags[i] == 2) || ((updated->flags[i] == 3) && (contact_settings->exact_voronoi_areas))) {

      atom = *atomp;

      set_contacts_atom(atom,system,calc_geometries);
    }
  }

  // update iarea values for contacts within selection:

  for (i=0,atomp=updated->atom;i<updated->natoms;i++,atomp++) {

    if (updated->flags[i] < 3) {

      atom = *atomp;

      contactlist = atom->contactlist;

      if (contactlist != NULL) {

	for (j=0,contact=contactlist->contacts;j<contactlist->ncontacts;j++,contact++) {

	  catom = contact->atom2;

	  if ((contact_settings->exact_voronoi_areas) || (atom_in_list(selection,catom))) {

	    icontact = find_contact(catom,atom);

	    if (icontact != NULL) {

	      contact->iarea = icontact->area;

	    } else {

	      contact->iarea = estimate_contact_iarea(contact);
	    }

	  } else {

	    contact->iarea = estimate_contact_iarea(contact);
	  }
	}
      }
    }
  }

  // set status for atoms in selection for which contacts changed:

  for (i=0,atomp=updated->atom;i<updated->natoms;i++,atomp++) {

    if (updated->flags[i] == 2) {

      atom = *atomp;
 
      atom->status |= ATOM_CHANGED;
    }
  }
  
  free_atomlist(updated);

  // set metal coordination mediated contacts:

  for (i=0,atomp=selection->atom;i<selection->natoms;i++,atomp++) {

    atom = *atomp;

    contactlist = atom->contactlist;

    if (contactlist != NULL) {

      for (j=0,contact=contactlist->contacts;j<contactlist->ncontacts;j++,contact++) {

	if (contact->flags & METAL_COORDINATING_CONTACT) {

	  flag_bond_mediated_contacts(contact->atom1,contact->atom2,selection);
	}
      }
    }
  }
}



void init_contactlist(CONTACTLIST *list) {

  list->ncontacts = 0;
  list->n_alloc_contacts = 0;
  list->contacts = NULL;
  list->poly = NULL;
}



static init_contact(CONTACT *contact) {

  contact->atom1 = NULL;
  contact->atom2 = NULL;
  contact->distance = 0.0;
  contact->alpha1 = 0.0;
  contact->beta1 = 0.0;
  contact->alpha2 = 0.0;
  contact->beta2 = 0.0;
  contact->area = 0.0;
  contact->iarea = 0.0;
  contact->score = 0.0;
  contact->poly = NULL;
  contact->flags = 0;
}



void add_contact_to_list(CONTACTLIST *list,ATOM *atom1,ATOM *atom2,double dist) {

  int i;
  ATOM **atom;
  ATOMLIST *conns;
  CONTACT *contact;

  if (list->n_alloc_contacts == 0) {

    list->n_alloc_contacts = 10;

    list->contacts = (CONTACT*) calloc(list->n_alloc_contacts,sizeof(CONTACT));

  } else if (list->ncontacts == list->n_alloc_contacts) {

    list->n_alloc_contacts *= 2;

    list->contacts = (CONTACT*) realloc(list->contacts,(list->n_alloc_contacts)*sizeof(CONTACT));
  }

  if (list->contacts == NULL) {

    error_fn("add_contact_to_list: out of memory allocating contactlist (%d contacts)",list->n_alloc_contacts);
  }

  contact = list->contacts + list->ncontacts;

  init_contact(contact);

  contact->atom1 = atom1;
  contact->atom2 = atom2;
  contact->distance = dist;

  if (contact->atom1->connections != NULL) {

    conns = atom1->connections;

    for (i=0,atom=conns->atom;i<conns->natoms;i++,atom++) {

      if (*atom == contact->atom2) {

	contact->flags |= COVALENT_CONTACT;

	break;
      }
    }
  }

  if ((contact->flags & COVALENT_CONTACT) || (same_molecule_atoms(atom1,atom2))) {

    contact->flags |= INTRAMOLECULAR_CONTACT;
  }

  if (metal_coordinating_contact(contact,0.4)) {

    contact->flags |= METAL_COORDINATING_CONTACT;
  }

  list->ncontacts++;
}



int contact_in_list(CONTACTLIST *list,CONTACT *contact1) {

  int i;
  CONTACT *contact2;

  for (i=0,contact2=list->contacts;i<list->ncontacts;i++,contact2++) {

    if (contact2 == contact1) {

      return(1);
    }
  }

  return(0);
}



static CONTACTLIST* alloc_contactlist(CONTACTLIST *list) {

  if (list == NULL) {

    list = (CONTACTLIST*) malloc(sizeof(CONTACTLIST));

    if (list == NULL) {

      error_fn("alloc_contact_list: out of memory allocating contactlist");
    }

    init_contactlist(list);

  } else {

    list->ncontacts = 0;
  }

  return(list);
}



void free_contactlist(CONTACTLIST *list) {

  if (list) {

    if (list->contacts) {

      free(list->contacts);
    }

    free(list);
  }
}



CONTACT* find_contact(ATOM* atom1,ATOM* atom2) {

  int i;
  CONTACT *contact;
  CONTACTLIST *contactlist;

  contactlist = atom1->contactlist;

  if (atom1->contactlist == NULL) {

    return(NULL);
  }

  for (i=0,contact=contactlist->contacts;i<contactlist->ncontacts;i++,contact++) {

    if (contact->atom2 == atom2) {

      return(contact);
    }
  }

  return(NULL);
}



void set_contacts_atom(ATOM *atom,SYSTEM *system,unsigned int calc_geometries) {

  int i;
  MOLECULE **molecule;
  CONTACT *contact;

  atom->contactlist = alloc_contactlist(atom->contactlist);

  for (i=0,molecule=system->molecule_list->molecules;i<system->molecule_list->n_molecules;i++,molecule++) {

    molecule2contacts(atom,*molecule);
  }

  atom->cstatus |= ATOM_CONTACTS_CALCULATED;

  if (contact_settings->voronoi_contacts == 1) {

    contacts2voronoi(atom);

    atom->cstatus |= ATOM_AREAS_CALCULATED;
  }

  if (calc_geometries) {

    for (i=0,contact=atom->contactlist->contacts;i<atom->contactlist->ncontacts;i++,contact++) {

      calc_contact_geometry(contact);
    }

    atom->cstatus |= ATOM_GEOMETRIES_CALCULATED;
  }
}



static void set_contacts_profile_atom(ATOM *atom) {

  int i;
  unsigned int flags;
  double area;
  ATOM_TYPE *type;
  CONTACT *contact;
  CONTACTLIST *contactlist;

  contactlist = atom->contactlist;

  if (atom->contactlist == NULL) {

    return;
  }

  atom->da_area = 0.0;
  atom->don_area = 0.0;
  atom->acc_area = 0.0;
  atom->ch_area = 0.0;
  atom->met_area = 0.0;
  atom->lip_area = 0.0;

  for (i=0,contact=contactlist->contacts;i<contactlist->ncontacts;i++,contact++) {

    if (!(contact->flags & COVALENT_CONTACT)) {

      type = contact->atom2->type;

      if (type == NULL) {

	error_fn("set_contacts_profile_atom: type undefined");
      }

      flags = type->flags;
      area = contact->area;

      if ((flags & HBOND_DA_ATOM_TYPE) == HBOND_DA_ATOM_TYPE) {

	atom->da_area += area;

      } else if (flags & HBOND_DONOR_ATOM_TYPE) {

	atom->don_area += area;

      } else if (flags & CHBOND_DONOR_ATOM_TYPE) {

	atom->ch_area += area;

      } else if ((flags & HBOND_ACCEPTOR_ATOM_TYPE) || (flags & METAL_ACCEPTOR_ATOM_TYPE)) {

	atom->acc_area += area;

      } else if (flags & METAL_ATOM_TYPE) {

	atom->met_area += area;

      } else {

	atom->lip_area += area;
      }
    }
  }

  atom->pol_area = atom->da_area + atom->don_area + atom->ch_area + atom->acc_area + atom->met_area;
}



static void write_contact_profile_list(FILE *file,ATOMLIST *list) {

  int i,sym_atom;
  ATOM **atomp,*atom;

  for (i=0,atomp=list->atom;i<list->natoms;i++,atomp++) {

    atom = *atomp;

    sym_atom = (atom->molecule->flags & SYMMETRY_MOLECULE) ? 1 : 0;

    fprintf(file,"PROFILE %6d %3d %10.4lf %10.4lf %10.4lf %10.4lf %10.4lf %10.4lf %10.4lf %10.4lf\n",
	    atom->id,sym_atom,
	    atom->pol_area,atom->lip_area,atom->exposed_area,
	    atom->da_area,atom->don_area,atom->ch_area,atom->acc_area,atom->met_area);
  }
}



static void molecule2contacts(ATOM *atom1,MOLECULE *molecule) {

  int i;
  ATOM *atom2;
  MAP *contacts_map;

  contacts_map = molecule->contacts_map;

  if (contacts_map) {

    map2contacts(atom1,contacts_map);

  } else {

    for (i=0,atom2=molecule->atom;i<molecule->natoms;i++,atom2++) {

      pair2contacts(atom1,atom2);
    }
  }
}



static void map2contacts(ATOM *atom1,MAP *contacts_map) {

  int i,j,iv[3],flag,irange[3][2],ix,iy,iz;
  GRID *grid;
  ATOM **atom2;
  ATOMLIST ***matrix,*list;

  grid = contacts_map->grid;
  matrix = (ATOMLIST***) contacts_map->matrix;

  flag = pos2grid(atom1->position,iv,grid);

  if (flag != 0) {

    return;
  }

  for (i=0;i<3;i++) {

    irange[i][0] = (iv[i] == 0) ? 0 : iv[i] - 1;
    irange[i][1] = (iv[i]+1 == grid->npoints[i]) ? iv[i] : iv[i] + 1;
  }

  for (ix=irange[0][0];ix<=irange[0][1];ix++) {

    for (iy=irange[1][0];iy<=irange[1][1];iy++) {

      for (iz=irange[2][0];iz<=irange[2][1];iz++) {

	list = &(matrix[ix][iy][iz]);

	for (i=0,atom2=list->atom;i<list->natoms;i++,atom2++) {

	  pair2contacts(atom1,*atom2);
	}
      }
    }
  }
}



static void pair2contacts(ATOM *atom1,ATOM *atom2) {

  double R1,R2,D;

  if ((atom1 != atom2) && (atom2->element != NULL) && (atom2->element->id != HYDROGEN)) {

    if (allow_contact(atom1,atom2)) {

      R1 = atom1->vdw_radius_H2O;
      R2 = atom2->vdw_radius_H2O;

      D = distance(atom1->position,atom2->position);

      if (D < R1 + R2) {

	add_contact_to_list(atom1->contactlist,atom1,atom2,D);
      }
    }
  }
}



int metal_coordinating_contact(CONTACT *contact,double tolerance) {

  double ideal_distance;
  ATOM *atom1,*atom2,*metal,*acceptor;

  atom1 = contact->atom1;
  atom2 = contact->atom2;

  if ((atom1->element == NULL) || (atom2->element == NULL)) {

    return(0);
  }

  metal = (atom1->element->flags & METAL_ELEMENT) ? atom1 : (atom2->element->flags & METAL_ELEMENT) ? atom2 : NULL;

  if (metal == NULL) {

    return(0);
  }

  acceptor = (contact->atom1 == metal) ? contact->atom2 : contact->atom1;

  if ((!(acceptor->type->flags & HBOND_ACCEPTOR_ATOM_TYPE)) && (!(acceptor->type->flags & METAL_ACCEPTOR_ATOM_TYPE))) {

    return(0);
  }

  ideal_distance = acceptor->vdw_radius + metal->vdw_radius - 1.0;

  if (contact->distance < ideal_distance + tolerance) {

    return(1);
  }

  return(0);
}



static void flag_atom_list_secondary_contacts(ATOMLIST *list) {

  int i,j;
  ATOM **atomp,*atom,**batomp,*batom;
  ATOMLIST *alist;
  CONTACT *contact;
  CONTACTLIST *clist;
  SETTINGS *settings;

  for (i=0,atomp=list->atom;i<list->natoms;i++,atomp++) {

    // loop contacts:

    atom = *atomp;
    clist = atom->contactlist;

    if (clist != NULL) {

      for (j=0,contact=clist->contacts;j<clist->ncontacts;j++,contact++) {

	if (secondary_contact(contact,clist)) {

	  contact->flags |= SECONDARY_CONTACT;
	}

	if (contact->flags & METAL_COORDINATING_CONTACT) {

	  flag_bond_mediated_contacts(contact->atom1,contact->atom2,list);
	}
      }
    }

    alist = atom->connections;

    // find covalent bonds between molecules

    if (alist != NULL) {

      for (j=0,batomp=alist->atom;j<alist->natoms;j++,batomp++) {

	batom = *batomp;

	if (!same_molecule_atoms(atom,batom)) {

	  //flag_bond_mediated_contacts(atom,batom,list);
	}
      }
    }
  }
}



int flag_clashing_contacts(ATOMLIST *list) {

  int i,j,clashing,covalent;
  double Do,D;
  ATOM **atomp,*atom1,*atom2;
  CONTACTLIST *clist;
  CONTACT *contact;

  clashing = covalent = 0;

  for (i=0,atomp=list->atom;i<list->natoms;i++,atomp++) {

    atom1 = *atomp;
    clist = atom1->contactlist;

    if (clist != NULL) {

      for (j=0,contact=clist->contacts;j<clist->ncontacts;j++,contact++) {

	if ((!(contact->flags & INTRAMOLECULAR_CONTACT)) && (!(contact->flags & SECONDARY_CONTACT))) {

	  atom2 = contact->atom2;

	  Do = atom1->element->cov_radius + atom2->element->cov_radius;

	  D = contact->distance;

	  if ((D < Do + COVALENT_TOLERANCE) && (allow_covalent_bond(atom1,atom2)))  {
	 
	    contact->flags |= COVALENT_CONTACT;

	    if (D < Do - COVALENT_TOLERANCE) {

	      atom1->error_flags |= ATOM_GEOMETRY_ERROR;
	      atom2->error_flags |= ATOM_GEOMETRY_ERROR;

	      clashing = 1;

	    } else {

	      covalent = 1;
	    }

	  } else if (!metal_coordinating_contact(contact,0.4)) {

	    Do = atom1->vdw_radius + atom2->vdw_radius;

	    if (D < Do - 2.0) {

	      atom1->error_flags |= ATOM_GEOMETRY_ERROR;
	      atom2->error_flags |= ATOM_GEOMETRY_ERROR;

	      clashing = 1;
	    }
	  }
	}
      }
    }
  }

  if (clashing) {

    return(2);

  } else if (covalent) {

    return(1);
  }

  return(0);
}



static int allow_contact(ATOM *atom1,ATOM *atom2) {

  if ((atom1->flags & SKIP_ATOM) || (atom2->flags & SKIP_ATOM)) {

    return(0);
  }

  if (((atom1->altloc != ' ') && (atom2->altloc != ' ')) && (atom1->altloc != atom2->altloc))
    return(0);

  return(1);
}



static void flag_bond_mediated_contacts(ATOM *atom1,ATOM *atom2,ATOMLIST *list) {

  int i,j;
  ATOM **atomp,*atom;
  CONTACT *contact;
  ATOMLIST *list1,*list2;

  list1 = (ATOMLIST*) alloc_atomlist();
  list2 = (ATOMLIST*) alloc_atomlist();

  init_atomlist(list1);
  init_atomlist(list2);

  get_connected_atoms(list1,atom1,NULL,1,0,3);
  get_connected_atoms(list2,atom2,NULL,1,0,3);

  for (i=0,atomp=list->atom;i<list->natoms;i++,atomp++) {

    atom = *atomp;

    if (atom->contactlist != NULL) {

      for (j=0,contact=atom->contactlist->contacts;j<atom->contactlist->ncontacts;j++,contact++) {

	if ((!(contact->flags & SECONDARY_CONTACT)) && (!(contact->flags & METAL_COORDINATING_CONTACT))) {

	  if (bond_mediated_contact(contact,list1,list2)) {

	    contact->flags |= SECONDARY_CONTACT;
	  }
	}
      }
    }
  }

  free_atomlist(list1);
  free_atomlist(list2);
}



static int bond_mediated_contact(CONTACT *contact,ATOMLIST *list1,ATOMLIST *list2) {

  int i,j,nbonds;
  unsigned int *dp1,*dp2;
  ATOM *atom1,*atom2,**ap1,**ap2,*a1,*a2;

  atom1 = contact->atom1;
  atom2 = contact->atom2;

  for (i=0,ap1=list1->atom,dp1=list1->flags;i<list1->natoms;i++,ap1++,dp1++) {

    a1 = *ap1;

    if ((a1 == atom1) || (a1 == atom2)) {

      for (j=0,ap2=list2->atom,dp2=list2->flags;j<list2->natoms;j++,ap2++,dp2++) {

	a2 = *ap2;

	if ((a2 != a1) && ((a2 == atom1) || (a2 == atom2))) {
      
	  nbonds = *dp1 + *dp2 + 1;

	  if ((nbonds > 1) && (nbonds < 5)) {

	    return(1);
	  }
	}
      }
    }
  }

  return(0);
}



void remove_flagged_contacts(CONTACTLIST *list,unsigned int flag) {

  int i,j,imax,nsec=0;
  CONTACT *contact1,*contact2;

  if (list->ncontacts == 0)
    return;

  // re-order list so secondary contacts are at the end:

  imax = list->ncontacts;

  for (i=0,contact1=list->contacts;i<imax;i++,contact1++) {

    if (contact1->flags & flag) {

      for (j=imax-1,contact2=list->contacts+imax-1;j>i;j--,contact2--) {

	if (!(contact2->flags & flag)) {

	  imax = j;

	  swap_contacts(contact1,contact2);

	  break;
	}
      }
    }
  }

  // update ncontacts:

  for (i=0,contact1=list->contacts;i<list->ncontacts;i++,contact1++) {

    if (contact1->flags & flag) {

      nsec++;
    }
  }

  list->ncontacts -= nsec;
}



static int secondary_contact(CONTACT *contact1,CONTACTLIST *list) {

  int i;
  double v1[4],v2[4],a,d;
  CONTACT *contact2;

  calc_vector(contact1->atom1->position,contact1->atom2->position,v1);

  for (i=0,contact2=list->contacts;i<list->ncontacts;i++,contact2++) {

    if (contact1 != contact2) {

      if (contact2->distance < contact1->distance) {

	calc_vector(contact2->atom1->position,contact2->atom2->position,v2);

	a = vector_angle(v1,v2);

	if (a < 90.0) {

	  d = (contact2->distance)*sin((PI/180.0)*a);

	  if (d < 1.00) {

	    return(1);
	  }
	}
      }
    }
  }

  return(0);
}



static void swap_contacts(CONTACT *contact1,CONTACT *contact2) {

  CONTACT contact;

  memcpy(&contact,contact1,sizeof(CONTACT));
  memcpy(contact1,contact2,sizeof(CONTACT));
  memcpy(contact2,&contact,sizeof(CONTACT));
}



int count_clashes_atom_list(ATOMLIST *list) {

  int i,j,uid1,uid2,n_clashes;
  double Dc,Do;
  ATOM **atomp,*atom1,*atom2;
  CONTACTLIST *contactlist;
  CONTACT *contact;

  flag_atom_list_secondary_contacts(list);

  n_clashes = 0;

  for (i=0,atomp=list->atom;i<list->natoms;i++,atomp++) {

    atom1 = *atomp;
    uid1 = atom1->unique_id;
    contactlist = atom1->contactlist;

    if ((contactlist != NULL) && (contactlist->ncontacts > 0)) {
	
      for (j=0,contact=contactlist->contacts;j<contactlist->ncontacts;j++,contact++) {
	
	atom2 = contact->atom2;
	uid2 = atom2->unique_id;

	if ((uid1 < uid2) || (!atom_in_list(list,atom2))) {

	  if ((!(contact->flags & COVALENT_CONTACT)) && (!(contact->flags & INTRAMOLECULAR_CONTACT)) &&
	      (!(contact->flags & SECONDARY_CONTACT))) {
	    
	    if (clash_contact(contact,&Dc,&Do)) {

	      n_clashes++;
	    }
	  }
	}
      }
    }
  }

  return(n_clashes);
}



void write_contacts_atom_list(PLI_FILE *file,ATOMLIST *list,enum OUTPUT_FORMAT oformat,unsigned long int oflags) {

  int i,j,uid1,uid2;
  ATOM **atomp,*atom1,*atom2;
  CONTACTLIST *contactlist;
  CONTACT *contact;

  flag_atom_list_secondary_contacts(list);

  for (i=0,atomp=list->atom;i<list->natoms;i++,atomp++) {

    atom1 = *atomp;
    uid1 = atom1->unique_id;
    contactlist = atom1->contactlist;

    if ((contactlist != NULL) && (contactlist->ncontacts > 0)) {
	
      for (j=0,contact=contactlist->contacts;j<contactlist->ncontacts;j++,contact++) {
	
	atom2 = contact->atom2;
	uid2 = atom2->unique_id;

	if ((uid1 < uid2) || (!atom_in_list(list,atom2))) {

	  if (((oflags & OUTPUT_COVALENT) || (!(contact->flags & COVALENT_CONTACT))) &&
	      ((oflags & OUTPUT_INTRAMOLECULAR) || (!(contact->flags & INTRAMOLECULAR_CONTACT)))) {
	    
	    write_contact(file,contact,oformat,oflags);
	  }
	}
      }
    }
  }
}



void write_clashes_atom_list(PLI_FILE *file,ATOMLIST *list,enum OUTPUT_FORMAT oformat,unsigned long int oflags) {

  int i,j,uid1,uid2;
  double Dc,Do;
  ATOM **atomp,*atom1,*atom2;
  CONTACTLIST *contactlist;
  CONTACT *contact;

  flag_atom_list_secondary_contacts(list);

  for (i=0,atomp=list->atom;i<list->natoms;i++,atomp++) {

    atom1 = *atomp;
    uid1 = atom1->unique_id;
    contactlist = atom1->contactlist;

    if ((contactlist != NULL) && (contactlist->ncontacts > 0)) {
	
      for (j=0,contact=contactlist->contacts;j<contactlist->ncontacts;j++,contact++) {
	
	atom2 = contact->atom2;
	uid2 = atom2->unique_id;

	if ((uid1 < uid2) || (!atom_in_list(list,atom2))) {

	  if ((!(contact->flags & COVALENT_CONTACT)) && (!(contact->flags & INTRAMOLECULAR_CONTACT)) &&
	      (!(contact->flags & SECONDARY_CONTACT))) {
	    
	    if (clash_contact(contact,&Dc,&Do)) {
      
	      write_clash(file,contact,Dc,Do,oformat,oflags);
	    }
	  }
	}
      }
    }
  }
}



void write_hbonds_atom_list(PLI_FILE *file,ATOMLIST *list,enum OUTPUT_FORMAT oformat,unsigned long int oflags) {

  int i,j,uid1,uid2;
  double score;
  ATOM **atomp,*atom1,*atom2;
  CONTACTLIST *contactlist;
  CONTACT *contact;

  flag_atom_list_secondary_contacts(list);

  for (i=0,atomp=list->atom;i<list->natoms;i++,atomp++) {

    atom1 = *atomp;
    uid1 = atom1->unique_id;
    contactlist = atom1->contactlist;

    if ((contactlist != NULL) && (contactlist->ncontacts > 0)) {
	
      for (j=0,contact=contactlist->contacts;j<contactlist->ncontacts;j++,contact++) {
	
	atom2 = contact->atom2;
	uid2 = atom2->unique_id;

	if ((uid1 < uid2) || (!atom_in_list(list,atom2))) {

	  if ((!(contact->flags & COVALENT_CONTACT)) && (!(contact->flags & INTRAMOLECULAR_CONTACT)) &&
	      (!(contact->flags & SECONDARY_CONTACT))) {

	    if (hbond_flags_match(atom1->type->flags,atom2->type->flags)) {

	      score = hbond_geometry_score(contact);

	      if ((score > 0.01) && (contact->score < 0.0)) {

		write_hbond(file,contact,score,oformat,oflags);
	      }
	    }
	  }
	}
      }
    }
  }
}



void write_contacts_atom(PLI_FILE *file,ATOM *atom,enum OUTPUT_FORMAT oformat,unsigned long int oflags) {

  int i;
  CONTACTLIST *contactlist;
  CONTACT *contact;

  contactlist = atom->contactlist;

  if ((contactlist != NULL) && (contactlist->ncontacts > 0)) {

    for (i=0,contact=contactlist->contacts;i<contactlist->ncontacts;i++,contact++) {

      if (((oflags & OUTPUT_COVALENT) || (!(contact->flags & COVALENT_CONTACT))) &&
	  ((oflags & OUTPUT_INTRAMOLECULAR) || (!(contact->flags & INTRAMOLECULAR_CONTACT)))) {
	
	calc_contact_geometry(contact);
	
	write_contact(file,contact,oformat,oflags);
      }
    }
  }
}



void write_clashes_atom(PLI_FILE *file,ATOM *atom,enum OUTPUT_FORMAT oformat,unsigned long int oflags) {

  int i;
  double Dc,Do;
  CONTACTLIST *contactlist;
  CONTACT *contact;

  contactlist = atom->contactlist;

  if ((contactlist != NULL) && (contactlist->ncontacts > 0)) {

    for (i=0,contact=contactlist->contacts;i<contactlist->ncontacts;i++,contact++) {
      
      if ((!(contact->flags & COVALENT_CONTACT)) && (!(contact->flags & INTRAMOLECULAR_CONTACT)) &&
	  (!(contact->flags & SECONDARY_CONTACT))) {
	
	if (clash_contact(contact,&Dc,&Do)) {
	  
	  write_clash(file,contact,Dc,Do,oformat,oflags);
	}
      }
    }
  }
}



void write_contact(PLI_FILE *file,CONTACT *contact,enum OUTPUT_FORMAT oformat,unsigned long int oflags) {

  int sym_flag1,sym_flag2,hbond,lipo;
  double hbond_score;
  ATOM *atom1,*atom2;
  MOLECULE *molecule1,*molecule2;

  atom1 = contact->atom1;
  atom2 = contact->atom2;

  molecule1 = atom1->molecule;
  molecule2 = atom2->molecule;

  sym_flag1 = (molecule1->flags & SYMMETRY_MOLECULE) ? 1 : 0;
  sym_flag2 = (molecule2->flags & SYMMETRY_MOLECULE) ? 1 : 0;

  lipo = hbond = 0;

  hbond_score = 0.0;

  if (hbond_flags_match(atom1->type->flags,atom2->type->flags)) {

    hbond_score = hbond_geometry_score(contact);
  }

  if ((contact->score < 0.0) || (!(oflags & OUTPUT_SCORES))) {

    if (hbond_score > 0.01) {

      hbond = 1;

    } else {

      if ((atom1->type->flags & LIPOPHILIC_ATOM_TYPE) && (atom2->type->flags & LIPOPHILIC_ATOM_TYPE)) {

	lipo = 1;
      }
    }
  }

  write_line(file,(oformat == JSON) ?
	     "{\"contact\":{\"aid1\":%d,\"set1\":%d,\"aid2\":%d,\"set2\":%d,\"distance\":%lf,\"area1\":%lf,\"area2\":%lf,\"secondary\":%d,\"hbond\":%d,\"lipo\":%d" :
	     "CONTACT %6d %6d %6d %6d      %9.4f %9.4f %9.4f     %1d     %1d     %1d",
	     atom1->id,atom1->set,
	     atom2->id,atom2->set,
	     contact->distance,
	     contact->area,
	     contact->iarea,
	     (contact->flags & SECONDARY_CONTACT) ? 1 : 0,
	     hbond,lipo);

  if (oflags & OUTPUT_GEOMETRIES) {

    if (oformat == JSON) {

      write_line(file,",");
    }

    write_contact_geometry(file,contact,oformat,oflags);
  }

  write_line(file,(oformat == JSON) ? ",\"sym_flag1\":%d,\"sym_flag2\":%d" : " %5d %5d",sym_flag1,sym_flag2);

  if (oflags & OUTPUT_SCORES) {

    if (oformat == JSON) {

      write_line(file,",");
    }

    write_contact_scores(file,contact,oformat,oflags);
  }

  write_line(file,(oformat == JSON) ? "}}\n" : "\n");
}



void write_clash(PLI_FILE *file,CONTACT *contact,double Dc,double Do,enum OUTPUT_FORMAT oformat,unsigned long int oflags) {

  int sym_flag1,sym_flag2;
  ATOM *atom1,*atom2;
  MOLECULE *molecule1,*molecule2;

  atom1 = contact->atom1;
  atom2 = contact->atom2;

  molecule1 = atom1->molecule;
  molecule2 = atom2->molecule;

  sym_flag1 = (molecule1->flags & SYMMETRY_MOLECULE) ? 1 : 0;
  sym_flag2 = (molecule2->flags & SYMMETRY_MOLECULE) ? 1 : 0;

  write_line(file,(oformat == JSON) ?
	     "{\"clash\":{\"aid1\":%d,\"set1\":%d,\"aid2\":%d,\"set2\":%d,\"distance\":%lf,\"Dc\":%lf,\"Do\":%lf,\"secondary\":%d" :
	     "CLASH %6d %6d %6d %6d      %9.4f %9.4f %9.4f     %1d",
	     atom1->id,atom1->set,
	     atom2->id,atom2->set,
	     contact->distance,
	     Dc,
	     Do,
	     (contact->flags & SECONDARY_CONTACT) ? 1 : 0);

  if (oflags & OUTPUT_GEOMETRIES) {

    if (oformat == JSON) {

      write_line(file,",");
    }

    write_contact_geometry(file,contact,oformat,oflags);
  }

  write_line(file,(oformat == JSON) ? ",\"sym_flag1\":%d,\"sym_flag2\":%d" : " %5d %5d",sym_flag1,sym_flag2);

  if (oflags & OUTPUT_SCORES) {

    if (oformat == JSON) {

      write_line(file,",");
    }

    write_contact_scores(file,contact,oformat,oflags);
  }

  write_line(file,(oformat == JSON) ? "}}\n" : "\n");
}



void write_hbond(PLI_FILE *file,CONTACT *contact,double score,enum OUTPUT_FORMAT oformat,unsigned long int oflags) {

  int sym_flag1,sym_flag2,chbond;
  unsigned int flags1,flags2;
  char direction;
  ATOM *atom1,*atom2;
  MOLECULE *molecule1,*molecule2;

  atom1 = contact->atom1;
  atom2 = contact->atom2;

  molecule1 = atom1->molecule;
  molecule2 = atom2->molecule;

  sym_flag1 = (molecule1->flags & SYMMETRY_MOLECULE) ? 1 : 0;
  sym_flag2 = (molecule2->flags & SYMMETRY_MOLECULE) ? 1 : 0;

  flags1 = (atom1->type) ? atom1->type->flags : 00;
  flags2 = (atom2->type) ? atom2->type->flags : 00;

  direction = '?';

  if ((atom1->element->id == CARBON) || ((flags1 & HBOND_DONOR_ATOM_TYPE) && (!(flags1 & HBOND_ACCEPTOR_ATOM_TYPE)))) {

    direction = 'D';

  } else if ((atom2->element->id == CARBON) || ((flags2 & HBOND_DONOR_ATOM_TYPE) && (!(flags2 & HBOND_ACCEPTOR_ATOM_TYPE)))) {

    direction = 'A';

  } else if ((flags1 & HBOND_ACCEPTOR_ATOM_TYPE) && (!(flags1 & HBOND_DONOR_ATOM_TYPE))) {

    direction = 'A';

  } else if ((flags2 & HBOND_ACCEPTOR_ATOM_TYPE) && (!(flags2 & HBOND_DONOR_ATOM_TYPE))) {

    direction = 'D';
  }

  write_line(file,(oformat == JSON) ?
	     "{\"hbond\":{\"aid1\":%d,\"set1\":%d,\"aid2\":%d,\"set2\":%d,\"distance\":%lf,\"chbond\":%d,\"secondary\":%d,\"direction\":\"%c\"" :
	     "HBOND %6d %6d %6d %6d      %9.4f     %1d    %1d    %c",
	     atom1->id,atom1->set,
	     atom2->id,atom2->set,
	     contact->distance,
	     (((atom1->element->id == CARBON) || (atom2->element->id == CARBON)) ? 1 : 0),
	     (contact->flags & SECONDARY_CONTACT) ? 1 : 0,
	     direction);

  if (oflags & OUTPUT_GEOMETRIES) {

    if (oformat == JSON) {

      write_line(file,",");
    }

    write_contact_geometry(file,contact,oformat,oflags);
  }

  write_line(file,(oformat == JSON) ? ",\"sym_flag1\":%d,\"sym_flag2\":%d" : " %5d %5d",sym_flag1,sym_flag2);

  if (oflags & OUTPUT_SCORES) {

    write_line(file,(oformat == JSON) ? ",\"hbond_score\":%.4lf," : "%10.4lf",score);

    write_contact_scores(file,contact,oformat,oflags);
  }

  write_line(file,(oformat == JSON) ? "}}\n" : "\n");
}



void write_sas_stats_atom_list(PLI_FILE *file,ATOMLIST *list,enum OUTPUT_FORMAT oformat) {

  int i,set,eid,uid;
  double sum_sas;
  UNITED_ATOM *uatom;
  ATOM_TYPE *type;
  ATOM_TYPING_SCHEME *scheme;

  scheme = get_atom_typing_scheme();

  for (i=0,type=scheme->atom_types;i<scheme->n_atom_types;i++,type++) {

    for (set=1;set<5;set++) {

      sum_sas = sum_sas_atomlist(list,type,set);

      if (sum_sas > 1.0E-6) {

	uatom = (type->united_atom) ? type->united_atom : NULL;
	uid = (uatom) ? uatom->id : -1;
	eid = (uatom) ? uatom->atom_node->element->id : (type->element) ? type->element->id : -1;

	write_line(file,(oformat == JSON) ? 
		   "{\"sumsas\":{\"set\":%d,\"eid\":%d,\"uid\":%d,\"tid\":%d,\"sas\":%.6lf}}\n" : "SUMSAS %5d %5d %5d %5d %10.4f\n",
		   set,eid,uid,type->id,sum_sas);
      }
    }
  }

  for (set=1;set<5;set++) {

    sum_sas = sum_sas_atomlist(list,NULL,set);

    if (sum_sas > 1.0E-6) {

      write_line(file,(oformat == JSON) ? 
		 "{\"sumsas\":{\"set\":%d,\"eid\":%d,\"uid\":%d,\"tid\":%d,\"sas\":%.6lf}}\n" : "SUMSAS %5d %5d %5d %5d %10.4f\n",
		 set,-1,-1,-1,sum_sas);
    }
  }  
}



static double sum_sas_atomlist(ATOMLIST *list,ATOM_TYPE *type,int set) {

  int i;
  double sum_sas;
  ATOM **atomp,*atom;

  sum_sas = 0.0;

  for (i=0,atomp=list->atom;i<list->natoms;i++,atomp++) {

    atom = *atomp;

    if ((atom->error_flags == 0) && (atom->type == type) && (atom->set == set)) {

      sum_sas += atom->exposed_area;
    }
  }

  return(sum_sas);
}



int clash_contact(CONTACT *contact,double *Dc,double *Do) {

  ATOM_TYPE *type1,*type2;
  FORCE_FIELD *ff;
  NONBONDED_FF *nonbonded_ff;

  ff = get_ff();

  type1 = contact->atom1->type;
  type2 = contact->atom2->type;

  nonbonded_ff = get_nonbonded_ff(ff->nonbonded,type1,type2);
  
  *Dc = nonbonded_ff->Dc;
  *Do = nonbonded_ff->Do;

  if (contact->distance < nonbonded_ff->Dc - 0.1) {

    return(1);
  }

  return(0);
}
