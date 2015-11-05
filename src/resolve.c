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



typedef struct RulesCoord {
  int n_metal;
  int n_acc;
  int n_don;
  int n_ch_don;
  int n_ch_acc;
} RULES_COORD;



static int report_alt_types = 0;
static int random_type_assignment = 0;



static int resolve_side_chains(SYSTEM*);
static int resolve_metals(SYSTEM*);
static int resolve_ambiguous_types(SYSTEM*);
static void check_hbond_mismatches(SYSTEM*);
static ATOMLIST* get_flippable_side_chain_atoms(SYSTEM*);
static int resolve_next_side_chain(ATOMLIST*,SYSTEM*,double);
static double flip_primary_amide(ATOM*,SYSTEM*,unsigned int);
static double flip_histidine(ATOM*,SYSTEM*,unsigned int);
static int resolve_next_type_from_contacts(SYSTEM*,ATOMLIST*,double,int);
static int resolve_unambiguous_types(ATOMLIST*,ATOM_TYPING_SCHEME*,int);
static int resolve_next_type_from_probabilities(ATOMLIST*,ATOM_TYPING_SCHEME*);
static ALT_ATOM_TYPE* alt_type_from_rules(SYSTEM*,ATOM*,double*);
static void count_rules_coord(ATOM*,RULES_COORD*,ATOM_TYPING_SCHEME*);
static ALT_ATOM_TYPE* alt_type_from_contacts(SYSTEM*,ATOM*,double*,int*);
static ATOMLIST *get_ring_unsubstituted_nitrogens(RING*);



void resolve_system(SYSTEM *system) {

  int n_resolved;

  // resolve side chains:

  n_resolved = resolve_side_chains(system);

  // check if metals look more like waters:

  n_resolved += resolve_metals(system);

  // resolve ambiguous types (tautomers and charge states):

  n_resolved += resolve_ambiguous_types(system);

  if (n_resolved) {
    
    // check for problematic donors/acceptors:
    
    check_hbond_mismatches(system);
  }
}



void unresolve_system(SYSTEM *system) {

  int i,j;
  ATOM *atom;
  MOLECULE **moleculep,*molecule;
  MOLECULE_LIST *list;
 
  list = system->molecule_list;

  if (list == NULL) {

    return;
  }

  for (i=0,moleculep=list->molecules;i<list->n_molecules;i++,moleculep++) {

    molecule = *moleculep;

    for (j=0,atom=molecule->atom;j<molecule->natoms;j++,atom++) {

      if (atom->ambiguous_type) {

	atom->type = atom->ambiguous_type;
      }
    }
  }
}



static int resolve_side_chains(SYSTEM *system) {

  int n_resolved;
  double score;
  ATOMLIST *atomlist;

  atomlist = get_flippable_side_chain_atoms(system);

  n_resolved = 0;

  while (resolve_next_side_chain(atomlist,system,-1.0)) {

    n_resolved++;
  }

  free_atomlist(atomlist);

  return(n_resolved);
}



static int resolve_metals(SYSTEM *system) {

  int i,j,sym_mol,n_resolved;
  double metal_energy,water_energy;
  ATOM *atom;
  ELEMENT *element;
  ATOM_TYPE *type,*water;
  ATOM_TYPING_SCHEME *scheme;
  MOLECULE **moleculep,*molecule;
  MOLECULE_LIST *list;
 
  list = system->molecule_list;

  if (list == NULL) {

    return(0);
  }

  if (system->settings->resolve_metals == 0) {

    return(0);
  }

  scheme = system->settings->atom_typing_scheme;

  water = get_atom_type("H2O",scheme);

  n_resolved = 0;

  for (i=0,moleculep=list->molecules;i<list->n_molecules;i++,moleculep++) {

    molecule = *moleculep;

    for (j=0,atom=molecule->atom;j<molecule->natoms;j++,atom++) {

      element = atom->element;
      type = atom->type;

      if ((type) && (type->flags & METAL_ATOM_TYPE)) {

	// calculate energy for metal:

	set_contacts_atom(atom,system,1);

	metal_energy = atom_hbond_energy(atom,0);

	// now calculate energy if this atom were a water:

	change_atom_type(atom,water);

	set_contacts_atom(atom,system,1);

	water_energy = atom_hbond_energy(atom,0);

	sym_mol = (molecule != NULL) ? ((molecule->flags & SYMMETRY_MOLECULE) ? 1 : 0) : -1;
	
	if (report_alt_types) {

	  printf("METAL %10d %5d %10d %10.4lf %10.4lf\n",atom->id,sym_mol,element->id,metal_energy,water_energy);
	}

	if (metal_energy - water_energy > 2.0) {

	  // leave as water:

	  warning_fn("changing metal atom %d to a water molecule",atom->id);

	  atom->error_flags |= ATOM_TYPING_ERROR;

	  atom->ambiguous_type = type;

	  n_resolved++;

	} else {

	  // change atom back:

	  change_atom_type(atom,type);

	  if (metal_energy - water_energy > 0.0) {

	    warning_fn("metal atom %d might be a water molecule",atom->id);

	    atom->error_flags |= ATOM_TYPING_ERROR;
	  }
	}
      }
    }
  }

  return(n_resolved);
}



static int resolve_ambiguous_types(SYSTEM *system) {

  int i,j,nc,nu,np,n_resolved;
  ATOM **atomp,*atom;
  ATOMLIST *atomlist;
  CONTACTLIST *contactlist;
  MOLECULE **moleculep,*molecule;
  MOLECULE_LIST *list;
  ATOM_TYPE *type;
  ATOM_TYPING_SCHEME *scheme;
  SETTINGS *settings;

  list = system->molecule_list;

  if (list == NULL) {

    return(0);
  }

  settings = system->settings;

  report_alt_types = settings->report_alt_types;

  random_type_assignment = (!strcmp(settings->alt_type_assignment,"random")) ? 1 : 0;

  // create list with ambiguous atoms:

  atomlist = alloc_atomlist();

  init_atomlist(atomlist);

  for (i=0,moleculep=list->molecules;i<list->n_molecules;i++,moleculep++) {

    molecule = *moleculep;

    for (j=0,atom=molecule->atom;j<molecule->natoms;j++,atom++) {

      type = atom->type;

      if ((type) && (type->n_alt_types)) {

	add_atom_to_list(atomlist,atom,0);

	set_contacts_atom(atom,system,1);

	calc_atom_coordination(atom);

	atom->ambiguous_type = type;
      }
    }
  }

  // iteratively resolve ambigous types:

  scheme = system->settings->atom_typing_scheme;

  np = 0;

  n_resolved = 0;

  do {

    do {

      nc = resolve_next_type_from_contacts(system,atomlist,0.5,np);

      nu = resolve_unambiguous_types(atomlist,scheme,np);

      n_resolved += nc + nu;

    } while ((nc) || (nu));

    np = resolve_next_type_from_probabilities(atomlist,scheme);

    n_resolved += np;

  } while (np);

  // write out unresolved atoms:

  if (atomlist->natoms) {

    warning_fn("resolve_ambiguous_types: unresolved atom types for the following atoms");

    for (i=0,atomp=atomlist->atom;i<atomlist->natoms;i++,atomp++) {

      write_atom(PLI_STDERR,*atomp,BASIC_ASTYLE,FORMATTED,0);
    }
  }

  free_atomlist(atomlist);

  return(n_resolved);
}



static void check_hbond_mismatches(SYSTEM *system) {

  int i,j,k;
  double Ea,score;
  ATOM *atom;
  ATOM_TYPE *type1,*type2;
  MOLECULE **moleculep,*molecule;
  MOLECULE_LIST *list;
  CONTACT *contact;
 
  list = system->molecule_list;

  if (list == NULL) {

    return;
  }

  for (i=0,moleculep=list->molecules;i<list->n_molecules;i++,moleculep++) {

    molecule = *moleculep;

    for (j=0,atom=molecule->atom;j<molecule->natoms;j++,atom++) {

      if (atom_residue_in_list(atom,system->selection)) {

	type1 = atom->type;

	if ((type1->flags & HBOND_DONOR_ATOM_TYPE) && (!(type1->flags & HBOND_ACCEPTOR_ATOM_TYPE))) {

	  set_contacts_atom(atom,system,1);

	  Ea = atom_hbond_energy(atom,0);

	  if (Ea > -1.0) {

	    for (k=0,contact=atom->contactlist->contacts;k<atom->contactlist->ncontacts;k++,contact++) {

	      type2 = contact->atom2->type;

	      if ((type2->flags & HBOND_DONOR_ATOM_TYPE) && (!(type2->flags & HBOND_ACCEPTOR_ATOM_TYPE))) {

		score = hbond_geometry_score(contact);

		if (score > 0.5) {

		  warning_fn("donor-donor hbond mismatch for atom %5d %c:%3s%d in molecule %s",
			     atom->id,atom->chain,atom->subname,atom->subid,molecule->name);

		  atom->error_flags |= ATOM_HBOND_MISMATCH;
		}
	      }
	    }
	  }

	} else if ((type1->flags & HBOND_ACCEPTOR_ATOM_TYPE) && (!(type1->flags & HBOND_DONOR_ATOM_TYPE))) {

	  set_contacts_atom(atom,system,1);

	  Ea = atom_hbond_energy(atom,0);

	  if (Ea > -1.0) {

	    for (k=0,contact=atom->contactlist->contacts;k<atom->contactlist->ncontacts;k++,contact++) {

	      type2 = contact->atom2->type;

	      if ((type2->flags & HBOND_ACCEPTOR_ATOM_TYPE) && (!(type2->flags & HBOND_DONOR_ATOM_TYPE))) {
	      
		score = hbond_geometry_score(contact);
	      
		if (score > 0.5) {
		
		  warning_fn("acceptor-acceptor hbond mismatch for atom %5d %c:%3s%d in molecule %s",
			   atom->id,atom->chain,atom->subname,atom->subid,molecule->name);
		
		  atom->error_flags |= ATOM_HBOND_MISMATCH;
		}
	      }
	    }
	  }
	}
      }
    }
  }
}



static ATOMLIST* get_flippable_side_chain_atoms(SYSTEM *system) {

  int i,j,flip_all;
  char *flip;
  ATOM *atom;
  ATOMLIST *atomlist,*selection;
  MOLECULE **moleculep,*molecule;
  MOLECULE_LIST *list;
  ATOM_TYPE *type,*primary_amide_o;
  ATOM_TYPING_SCHEME *scheme;

  scheme = system->settings->atom_typing_scheme;

  primary_amide_o = get_atom_type("Primary amide =O",scheme);

  atomlist = alloc_atomlist();

  init_atomlist(atomlist);

  list = system->molecule_list;

  if (list == NULL) {

    return(atomlist);
  }

  flip = system->settings->resolve_side_chains;

  if (!strcmp(flip,"none")) {

    return(atomlist);
  }

  flip_all = (!strcmp(flip,"all")) ? 1 : 0;

  for (i=0,moleculep=list->molecules;i<list->n_molecules;i++,moleculep++) {

    molecule = *moleculep;

    if (molecule->flags & PROTEIN_MOLECULE) {

      selection = (!strcmp(flip,"site")) ? molecule2site(molecule) : molecule->selection;

      for (j=0,atom=molecule->atom;j<molecule->natoms;j++,atom++) {

	if ((atom_residue_in_list(atom,selection)) || (flip_all)) {

	  type = atom->type;

	  if (type == primary_amide_o) {

	    add_atom_to_list(atomlist,atom,1);

	  } else if ((!strcmp(atom->subname,"HIS")) && (!strcmp(atom->name," NE2"))) {

	    add_atom_to_list(atomlist,atom,2);
	  }
	}
      }

      if (!strcmp(flip,"site")) {

	free_atomlist(selection);
      }
    }
  }

  return(atomlist);
}



static int resolve_next_side_chain(ATOMLIST *atomlist,SYSTEM *system,double max_energy) {

  int i;
  unsigned int *flagp,flag;
  double energy,min_energy;
  ATOM **atomp,*atom,*best_atom;
  ATOM_TYPE *type,*primary_amide_o;
  ATOM_TYPING_SCHEME *scheme;
  double (*flip)(ATOM*,SYSTEM*,unsigned int);

  scheme = system->settings->atom_typing_scheme;

  primary_amide_o = get_atom_type("Primary amide =O",scheme);

  min_energy = 0.0;

  for (i=0,atomp=atomlist->atom,flagp=atomlist->flags;i<atomlist->natoms;i++,atomp++,flagp++) {

    atom = *atomp;
    flag = *flagp;

    if (flag == 1) {

      energy = flip_primary_amide(atom,system,0);

      if (energy < min_energy) {

	min_energy = energy;
	best_atom = atom;
	flip = flip_primary_amide;
      }

    } else if (flag == 2) {

      energy = flip_histidine(atom,system,0);

      if (energy < min_energy) {

	min_energy = energy;
	best_atom = atom;
	flip = flip_histidine;
      }
    }
  }

  if (min_energy < max_energy) {

    flip(best_atom,system,1);

    warning_fn("flipping atom %d in residue %s%d (energy=%.2lf)",best_atom->id,best_atom->subname,best_atom->subid,min_energy);

    remove_atom_from_list(atomlist,best_atom);

    return(1);
  }

  return(0);
}



static double flip_primary_amide(ATOM *oxygen,SYSTEM *system,unsigned int keep) {

  int i;
  double m[4][4],energy1,energy2;
  ATOM *carbon,*nitrogen,*other_atom,**atomp,*atom;
  ATOM_TYPE *primary_amide_n;
  ATOM_TYPING_SCHEME *scheme;

  if ((!oxygen) || (!oxygen->connections) || (oxygen->connections->natoms != 1)) {

    return(0.0);
  }

  carbon = oxygen->connections->atom[0];

  if ((!carbon) || (!carbon->connections) || (carbon->connections->natoms != 3)) {

    return(0.0);
  }

  scheme = system->settings->atom_typing_scheme;

  primary_amide_n = get_atom_type("Amide -NH2",scheme);

  nitrogen = other_atom = NULL;

  for (i=0,atomp=carbon->connections->atom;i<3;i++,atomp++) {

    atom = *atomp;

    if (atom != oxygen) {

      if (atom->type == primary_amide_n) { 

	nitrogen = atom;

      } else {

	other_atom = atom;
      }
    }
  }

  if ((!nitrogen) || (!other_atom) || (!nitrogen->connections) || (nitrogen->connections->natoms != 1)) {

    return(0.0);
  }

  set_contacts_atom(oxygen,system,1);
  set_contacts_atom(nitrogen,system,1);

  energy1 = atom_hbond_energy(oxygen,0) + atom_hbond_energy(nitrogen,0);

  rotation_matrix(other_atom->position,carbon->position,180.0,m);

  transform_atom(oxygen,m);
  transform_atom(nitrogen,m);

  set_contacts_atom(oxygen,system,1);
  set_contacts_atom(nitrogen,system,1);

  energy2 = atom_hbond_energy(oxygen,0) + atom_hbond_energy(nitrogen,0);

  if (!keep) {

    // flip back:

    transform_atom(oxygen,m);
    transform_atom(nitrogen,m);
  }

  return(energy2 - energy1);
}



static double flip_histidine(ATOM *NE2,SYSTEM *system,unsigned int keep) {

  int i;
  double m[4][4],energy1,energy2;
  ATOM **atomp,*atom,*CB,*CG,*ND1,*CD2,*CE1;
  RING *ring;

  if ((!NE2) || (!NE2->connections) || (NE2->connections->natoms != 2)) {

    return(0.0);
  }

  if (!NE2->ring) {

    find_ring(NE2,5);
  }

  ring = NE2->ring;

  if ((!ring) || (ring->size != 5)) {

    return(0.0);
  }

  CG = ND1 = CD2 = CE1 = NULL;

  for (i=0,atomp=ring->atom;i<5;i++,atomp++) {

    atom = *atomp;

    if (!strcmp(atom->name," CG ")) {

      CG = atom;

    } else if (!strcmp(atom->name," ND1")) {

      ND1 = atom;

    } else if (!strcmp(atom->name," CD2")) {

      CD2 = atom;

    } else if (!strcmp(atom->name," CE1")) {

      CE1 = atom;
    }
  }

  if ((!CG) || (!ND1) || (!CD2) || (!CE1) || 
      (!CG->connections) || (!ND1->connections) || (!CD2->connections) || (!CE1->connections) ||
      (CG->connections->natoms != 3) || (ND1->connections->natoms != 2) || 
      (CD2->connections->natoms != 2) || (CE1->connections->natoms != 2)) {

    return(0.0);
  }

  CB = NULL;

  for (i=0,atomp=CG->connections->atom;i<3;i++,atomp++) {

    atom = *atomp;

    if ((atom != CD2) && (atom != ND1)) {

      CB = atom;
    }
  }

  if (!CB) {

    return(0.0);
  }

  set_contacts_atom(ND1,system,1);
  set_contacts_atom(NE2,system,1);

  energy1 = atom_hbond_energy(ND1,1) + atom_hbond_energy(NE2,1);

  rotation_matrix(CB->position,CG->position,180.0,m);

  for (i=0,atomp=ring->atom;i<5;i++,atomp++) {

    transform_atom(*atomp,m);
  }

  set_contacts_atom(ND1,system,1);
  set_contacts_atom(NE2,system,1);

  energy2 = atom_hbond_energy(ND1,1) + atom_hbond_energy(NE2,1);

  if (!keep) {

    // flip back:

    for (i=0,atomp=ring->atom;i<5;i++,atomp++) {

      transform_atom(*atomp,m);
    }
  }

  return(energy2 - energy1);
}



static int resolve_next_type_from_contacts(SYSTEM *system,ATOMLIST *atomlist,double min_p,int flag) {

  int i,sym_mol,rule,rule_based;
  double p,p_max;
  ATOM **atomp,*atom,*best_atom;
  ATOM_TYPE *type;
  ALT_ATOM_TYPE *alt_type,*best_alt_type;
  MOLECULE **moleculep,*molecule;
  ATOM_TYPING_SCHEME *scheme;

  scheme = system->settings->atom_typing_scheme;

  // find the ambiguous atom that is most unambiguous:

  best_atom = NULL;
  best_alt_type = NULL;
  p_max = 0.0;
  rule_based = 0;

  for (i=0,atomp=atomlist->atom;i<atomlist->natoms;i++,atomp++) {

    atom = *atomp;

    alt_type = alt_type_from_contacts(system,atom,&p,&rule);

    if ((alt_type) && ((p > p_max) || ((rule) && (!rule_based)))) {

      best_atom = atom;
      best_alt_type = alt_type;
      p_max = p;

      rule_based = rule;
    }
  }
  
  if ((best_atom) && (best_alt_type) && (p_max > min_p)) {

    molecule = best_atom->molecule;

    sym_mol = (molecule != NULL) ? ((molecule->flags & SYMMETRY_MOLECULE) ? 1 : 0) : -1;
    
    if (report_alt_types) {

      printf("ALT_TYPE %10d %5d %-5s %-5s %5d %10d %10d %10d %10.6lf %10d\n",
	     best_atom->id,sym_mol,best_atom->name,best_atom->subname,best_atom->subid,best_atom->type->id,
	     best_alt_type->type->id,flag,p_max,best_atom->group_status);
    }

    update_group_statuses(best_atom,best_alt_type,scheme);

    best_atom->type = best_alt_type->type;
    best_atom->type_probability = p_max;

    remove_atom_from_list(atomlist,best_atom);

    return(1);
  }

  return(0);
}



static int resolve_unambiguous_types(ATOMLIST *atomlist,ATOM_TYPING_SCHEME *scheme,int flag) {

  int i,j,n,n_alt_types,sym_mol;
  ATOM **atomp,*atom;
  ATOM_TYPE *type;
  ALT_ATOM_TYPE *alt_type,*best_alt_type;
  MOLECULE *molecule;

  if (atomlist->natoms == 0) {

    return(0);
  }

  n = 0;

  for (i=0,atomp=atomlist->atom;i<atomlist->natoms;i++,atomp++) {

    atom = *atomp;

    type = atom->type;

    if ((type) && (type->n_alt_types)) {

      n_alt_types = 0;
      best_alt_type = NULL;

      for (j=0,alt_type=type->alt_types;j<type->n_alt_types;j++,alt_type++) {

	if (alt_type->group_status == atom->group_status) {

	  best_alt_type = alt_type;

	  n_alt_types++;
	}
      }

      if ((n_alt_types == 1) && (best_alt_type)) {

	molecule = atom->molecule;

	sym_mol = (molecule != NULL) ? ((molecule->flags & SYMMETRY_MOLECULE) ? 1 : 0) : -1;

	if (report_alt_types) {

	  printf("ALT_TYPE %10d %5d %-5s %-5s %5d %10d %10d %10d %10.6lf %10d\n",
		 atom->id,sym_mol,atom->name,atom->subname,atom->subid,atom->type->id,
		 best_alt_type->type->id,flag,best_alt_type->probability,atom->group_status);
	}

	update_group_statuses(atom,best_alt_type,scheme);

	atom->type = best_alt_type->type;
	atom->type_probability = best_alt_type->probability;

	remove_atom_from_list(atomlist,atom);

	n++;
      }
    }
  }

  return(n);
}



static ALT_ATOM_TYPE* alt_type_from_contacts(SYSTEM *system,ATOM *atom,double *probability,int *rule) {

  int i;
  double p,p0,p1,maxp,sump,k0,k1;
  ATOM_TYPE *type;
  ALT_ATOM_TYPE *best_alt_type,*alt_type;

  // first check if atom falls under any of the rules for H-bonding atom types:

  best_alt_type = alt_type_from_rules(system,atom,probability);

  if (best_alt_type) {

    //*rule = 1;

    //return(best_alt_type);
  }

  // now do general, probability-based assessment:

  *rule = 0;

  maxp = sump = 0.0;

  type = atom->type;
  best_alt_type = NULL;

  for (i=0,alt_type=type->alt_types;i<type->n_alt_types;i++,alt_type++) {

    if (alt_type->group_status == atom->group_status) {

      atom->type = alt_type->type;

      p = (alt_type->probability)*exp(-atom_hbond_energy(atom,0)/RT);

      if (p > maxp) {

	maxp = p;

	best_alt_type = alt_type;
      }

      sump += p;
    }
  }

  atom->type = type;

  if (sump > 1.0E-6) {

    p0 = best_alt_type->probability;

    p1 = maxp/sump;

    k0 = (1.0-p0)/p0;

    k1 = (1.0-p1)/p1;

    if ((k1/k0) < 0.8) {

      *probability = p1;

      return(best_alt_type);
    }
  }

  *probability = 0.0;

  return(NULL);
}



static ALT_ATOM_TYPE* alt_type_from_rules(SYSTEM *system,ATOM *atom,double *p) {

  int i;
  ATOM *oxygen;
  ATOM_TYPE *type;
  ALT_ATOM_TYPE *alt_type;
  ATOM_TYPING_SCHEME *scheme;
  RULES_COORD rcoord;

  scheme = system->settings->atom_typing_scheme;

  type = atom->type;

  if (!type) {

    return(NULL);
  }

  if (type == get_atom_type("Carboxyl O?",scheme)) {

    // carboxylic acid oxygens:

    count_rules_coord(atom,&rcoord,scheme);

    if ((rcoord.n_acc) && (!(rcoord.n_metal))) {

      // if the oxygen H-bonds to an acceptor, it becomes a -OH
      // this has higher p, because it needs to happen first

      *p = 2.0;

      return(get_alt_type(type,"Carboxyl -OH",atom->group_status));

    } else if ((rcoord.n_ch_don) || (rcoord.n_metal)) {

      // if the oxygen H-bonds to a charged donor, it becomes a =O

      *p = 1.0;

      return(get_alt_type(type,"Carboxyl =O",atom->group_status));

    } else if (rcoord.n_don > 1) {

      // if the oxygen H-bonds to at least two donors, it becomes a =O

      *p = 1.0;

      return(get_alt_type(type,"Carboxyl =O",atom->group_status));
      
    } else {

      // if the other oxygen on this carboxyl group doesn't H-bond to an acceptor
      // and does H-bond to a charged donor, then this atom becomes a =O

      oxygen = carboxyl_o_from_o(atom);

      if (oxygen) {

	count_rules_coord(oxygen,&rcoord,scheme);

	if ((!rcoord.n_acc) && (rcoord.n_ch_don)) {

	  *p = 1.0;

	  return(get_alt_type(type,"Carboxyl =O",atom->group_status));
	}
      }
    }

  } else if (type == get_atom_type("-N?",scheme)) {

    // primary amine:

    count_rules_coord(atom,&rcoord,scheme);

    if (rcoord.n_don) {

      // if the nitrogen H-bonds to a donor, it becomes a -NH2
      // this has higher p, because it needs to happen first

      *p = 2.0;

      return(get_alt_type(type,"-NH2",atom->group_status));

    } else if (rcoord.n_ch_acc) {

      // if the nitrogen H-bonds to a charged acceptor it becomes a -NH3

      *p = 1.0;

      return(get_alt_type(type,"-NH3",atom->group_status));

    } else if (rcoord.n_acc > 2) {

      // if the nitrogen H-bonds to at least three acceptors it becomes a -NH3

      *p = 1.0;

      return(get_alt_type(type,"-NH3",atom->group_status));
    }

  } else if (type == get_atom_type("-CH2-N?",scheme)) {

    // primary amine:

    count_rules_coord(atom,&rcoord,scheme);

    if (rcoord.n_don) {

      // if the nitrogen H-bonds to a donor, it becomes a -NH2
      // this has higher p, because it needs to happen first

      *p = 2.0;

      return(get_alt_type(type,"-NH2",atom->group_status));

    } else if (rcoord.n_ch_acc) {

      // if the nitrogen H-bonds to a charged acceptor it becomes a -CH2-NH3

      *p = 1.0;

      return(get_alt_type(type,"-CH2-NH3",atom->group_status));

    } else if (rcoord.n_acc > 2) {

      // if the nitrogen H-bonds to at least three acceptors it becomes a -CH2-NH3

      *p = 1.0;

      return(get_alt_type(type,"-CH2-NH3",atom->group_status));
    }
  }

  return(NULL);
}



static void count_rules_coord(ATOM *atom,RULES_COORD *rcoord,ATOM_TYPING_SCHEME *scheme) {

  int i;
  unsigned int flags2;
  ATOM *atom2;
  ATOM_TYPE *type2;
  UNITED_ATOM *uatom2;
  ATOM_TYPE *carboxyl_o,*carboxyl_o2,*primary_amine_n,*primary_amine_nh3,*primary_amine_ch2_n,*primary_amine_ch2_nh3;
  ATOM_COORDINATION *coord;
  COORDINATION *item;

  carboxyl_o = get_atom_type("Carboxyl O?",scheme);
  carboxyl_o2 = get_atom_type("Carboxyl =O",scheme);

  primary_amine_n = get_atom_type("-N?",scheme);
  primary_amine_nh3 = get_atom_type("-NH3",scheme);
  primary_amine_ch2_n = get_atom_type("-CH2-N?",scheme);
  primary_amine_ch2_nh3 = get_atom_type("-CH2-NH3",scheme);

  rcoord->n_metal = rcoord->n_acc = rcoord->n_don = rcoord->n_ch_don = rcoord->n_ch_acc = 0;

  coord = atom->coordination;

  if ((!coord) || (coord->max_z == 0)) {

    return;
  }

  for (i=1,item=coord->list+1;i<=coord->max_z;i++,item++) {
    
    if (item->s > 0.5) {

      atom2 = item->contact->atom2;
      type2 = atom2->type;
      uatom2 = type2->united_atom;
      flags2 = type2->flags;

      if (type2) {

	if (flags2 & METAL_ATOM_TYPE) {

	  (rcoord->n_metal)++;
	
	} else if ((uatom2) && ((type2 == carboxyl_o) || (type2 == carboxyl_o2) ||
				((flags2 & HBOND_ACCEPTOR_ATOM_TYPE) && (!(flags2 & HBOND_DONOR_ATOM_TYPE))))) {

	  (rcoord->n_acc)++;

	  if ((uatom2->formal_charge < 0) || (type2 == carboxyl_o) || (type2 == carboxyl_o2)) {

	    (rcoord->n_ch_acc)++;
	  }
	      
	} else if ((uatom2) && ((type2 == primary_amine_n) || (type2 == primary_amine_nh3) ||
				(type2 == primary_amine_ch2_n) || (type2 == primary_amine_ch2_nh3) ||
				((flags2 & HBOND_DONOR_ATOM_TYPE) && (!(flags2 & HBOND_ACCEPTOR_ATOM_TYPE))))) {

	  (rcoord->n_don)++;

	  if ((uatom2->formal_charge > 0) ||
	      (type2 == primary_amine_n) || (type2 == primary_amine_nh3) ||
	      (type2 == primary_amine_ch2_n) || (type2 == primary_amine_ch2_nh3)) {

	    (rcoord->n_ch_don)++;
	  }
	}
      }
    }
  }
}



static int resolve_next_type_from_probabilities(ATOMLIST *atomlist,ATOM_TYPING_SCHEME *scheme) {

  int atom_index,alt_type_index,i,n_alt_types,sym_mol;
  double p;
  ATOM_TYPE *type;
  ALT_ATOM_TYPE **alt_types,*alt_type,*atype;
  ATOM *atom;
  MOLECULE *molecule;

  if (atomlist->natoms == 0) {

    return(0);
  }

  // pick next ambiguous atom:

  atom_index = (random_type_assignment) ? irand() % (atomlist->natoms) : 0;

  atom = *(atomlist->atom + atom_index);

  type = atom->type;

  if ((type) && (type->n_alt_types)) {

    // construct a list of potential alternative atome types:

    alt_types = (ALT_ATOM_TYPE**) calloc(type->n_alt_types,sizeof(ALT_ATOM_TYPE*));

    if (!alt_types) {

      error_fn("resolve_next_type_from_probabilities: out of memory allocating alt types");
    }

    n_alt_types = 0;

    for (i=0,alt_type=type->alt_types;i<type->n_alt_types;i++,alt_type++) {

      if (alt_type->group_status == atom->group_status) {

	alt_types[n_alt_types++] = (type->alt_types+i);
      }
    }

    if (n_alt_types) {

      alt_type = NULL;

      if (random_type_assignment) {

	// pick an alternative atom type at random:

	do {
	
	  alt_type_index = irand() % (n_alt_types);
	
	  alt_type = *(alt_types + alt_type_index);

	  // check against probability:

	  p = frand();

	} while (p > alt_type->probability);

      } else {

	// pick the most probable atom type:

	p = 0.0;

	for (i=0;i<n_alt_types;i++) {

	  atype = *(alt_types + i);

	  if (atype->probability > p) {

	    p = atype->probability + 1.0E-10; // adding a small number ensures the first type is always selected if the probabilities are the same

	    alt_type = atype;
	  }
	}
      }

      if (alt_type) {

	molecule = atom->molecule;

	sym_mol = (molecule != NULL) ? ((molecule->flags & SYMMETRY_MOLECULE) ? 1 : 0) : -1;

	if (report_alt_types) {

	  printf("ALT_TYPE %10d %5d %-5s %-5s %5d %10d %10d %10d %10.6lf %10d\n",
		 atom->id,sym_mol,atom->name,atom->subname,atom->subid,atom->type->id,
		 alt_type->type->id,1,alt_type->probability,atom->group_status);
	}

	update_group_statuses(atom,alt_type,scheme);

	atom->type = alt_type->type;
	atom->type_probability = alt_type->probability;
      }
    }

    free(alt_types);
  }

  remove_atom_from_list(atomlist,atom);

  return(1);
}



void update_group_statuses(ATOM *atom,ALT_ATOM_TYPE *alt_type,ATOM_TYPING_SCHEME *scheme) {

  int i;
  ATOM **nitrogenp,*nitrogen,*oxygen;
  ATOM_TYPE *type;
  ATOMLIST *alist;

  type = atom->type;

  if ((!strcmp(type->name,"Imidazole N?")) ||
      (!strcmp(type->name,"Pyrazole N?")) ||
      (!strcmp(type->name,"Triazole N?")) ||
      (!strcmp(type->name,">N HIS-ND1")) ||
      (!strcmp(type->name,">N HIS-NE2"))) {

    find_ring(atom,5);

    alist = get_ring_unsubstituted_nitrogens(atom->ring);
    
    for (i=0,nitrogenp=alist->atom;i<alist->natoms;i++,nitrogenp++) {

      (*nitrogenp)->group_status = alt_type->assigned_group_status;
    }

    free_atomlist(alist);

  } else if (!strcmp(type->name,"Pyridone N?")) {

    atom->group_status = alt_type->assigned_group_status;

    oxygen = pyridone_o_from_n(atom);

    if (oxygen) {

      oxygen->group_status = alt_type->assigned_group_status;
    }

  } else if (!strcmp(type->name,"Pyridone O?")) {

    atom->group_status = alt_type->assigned_group_status;

    nitrogen = pyridone_n_from_o(atom);

    if (nitrogen) {

      nitrogen->group_status = alt_type->assigned_group_status;
    }

  } else if (!strcmp(type->name,"Carboxyl O?")) {

    atom->group_status = alt_type->assigned_group_status;

    oxygen = carboxyl_o_from_o(atom);

    if (oxygen) {

      oxygen->group_status = alt_type->assigned_group_status;
    }
  }
}



static ATOMLIST *get_ring_unsubstituted_nitrogens(RING *ring) {

  int i;
  ATOM **atomp,*atom;
  ATOM_NODE *node;
  ATOMLIST *list;

  list = alloc_atomlist();

  init_atomlist(list);

  if (ring) {

    for (i=0,atomp=ring->atom;i<ring->size;i++,atomp++) {

      atom = *atomp;

      node = atom->node;

      if ((!strcmp(node->name,"=N-")) || (!strcmp(node->name,">N"))) {

	add_atom_to_list(list,atom,0);
      }
    }
  }

  return(list);
}


