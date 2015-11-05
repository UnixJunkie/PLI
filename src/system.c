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



static SYSTEM_LIST *ligands2syslist(SETTINGS*);
static SYSTEM_LIST *complexes2syslist(SETTINGS*);
static SYSTEM_LIST *protein2syslist(SETTINGS*);
static void init_system(SYSTEM*);
static void setup_system(SYSTEM*,MOLECULE*,MOLECULE*,MOLECULE*,MOLECULE*,SETTINGS*);
static void init_syslist(SYSTEM_LIST*);
static void select_system_atoms(SYSTEM*);
static void check_system_ligand(PLI_FILE *file,SYSTEM*,enum OUTPUT_FORMAT oformat);
static void write_system(PLI_FILE*,SYSTEM*,enum OUTPUT_FORMAT,unsigned long int);
static ATOMLIST* extended_system_atomlist(SYSTEM*);



SYSTEM* settings2system(SETTINGS *settings) {

  SYSTEM *system;
  SYSMOLS *sysmols;

  system = (SYSTEM*) malloc(sizeof(SYSTEM));

  if (system == NULL) {

    error_fn("settings2system: out of memory allocating system");
  }

  init_system(system);

  system->id = 0;

  sysmols = settings->sysmols;

  setup_system(system,sysmols->protein,sysmols->water,sysmols->ligand,sysmols->symmetry,settings);

  return(system);
}



SYSTEM_LIST* settings2syslist(SETTINGS *settings) {

  SYSTEM_LIST *syslist;

  if (settings->sysmols->ligand_list) {

    syslist = ligands2syslist(settings);

  } else if (settings->sysmols->complex_list) {

    syslist = complexes2syslist(settings);

  } else {

    syslist = protein2syslist(settings);
  }

  return(syslist);
}


static SYSTEM_LIST *ligands2syslist(SETTINGS *settings) {

  int i;
  MOLECULE **ligandp;
  SYSTEM_LIST *syslist;
  SYSMOLS *sysmols;
  SYSTEM *system;

  syslist = (SYSTEM_LIST*) malloc(sizeof(SYSTEM_LIST));

  if (syslist == NULL) {

    error_fn("ligands2syslist: out of memory allocating syslist");
  }

  init_syslist(syslist);

  sysmols = settings->sysmols;

  syslist->n_systems = sysmols->ligand_list->n_molecules;

  syslist->systems = (SYSTEM*) calloc(syslist->n_systems,sizeof(SYSTEM));

  if (!syslist->systems) {

    error_fn("ligands2syslist: out of memory allocating systems");
  }

  for (i=0,system=syslist->systems,ligandp=sysmols->ligand_list->molecules;i<syslist->n_systems;i++,system++,ligandp++) {

    init_system(system);

    system->id = i;

    setup_system(system,sysmols->protein,sysmols->water,*ligandp,sysmols->symmetry,settings);
  }

  return(syslist);
}


static SYSTEM_LIST *complexes2syslist(SETTINGS *settings) {

  int i;
  COMPLEX *complex;
  SYSTEM_LIST *syslist;
  SYSMOLS *sysmols;
  SYSTEM *system;

  syslist = (SYSTEM_LIST*) malloc(sizeof(SYSTEM_LIST));

  if (syslist == NULL) {

    error_fn("complexes2syslist: out of memory allocating syslist");
  }

  init_syslist(syslist);

  sysmols = settings->sysmols;

  syslist->n_systems = sysmols->complex_list->n_complexes;

  syslist->systems = (SYSTEM*) calloc(syslist->n_systems,sizeof(SYSTEM));

  if (!syslist->systems) {

    error_fn("complexes2syslist: out of memory allocating systems");
  }

  for (i=0,system=syslist->systems,complex=sysmols->complex_list->complexes;i<syslist->n_systems;i++,system++,complex++) {

    init_system(system);

    system->id = i;

    setup_system(system,complex->protein,NULL,complex->ligand,NULL,settings);
  }

  return(syslist);
}



static SYSTEM_LIST *protein2syslist(SETTINGS *settings) {

  int i;
  MOLECULE **ligandp;
  SYSTEM_LIST *syslist;
  SYSMOLS *sysmols;
  SYSTEM *system;

  syslist = (SYSTEM_LIST*) malloc(sizeof(SYSTEM_LIST));

  if (syslist == NULL) {

    error_fn("protein2syslist: out of memory allocating syslist");
  }

  init_syslist(syslist);

  sysmols = settings->sysmols;

  syslist->n_systems = 1;

  syslist->systems = (SYSTEM*) calloc(syslist->n_systems,sizeof(SYSTEM));

  if (!syslist->systems) {

    error_fn("protein2syslist: out of memory allocating systems");
  }

  system = syslist->systems;

  init_system(system);

  system->id = i;

  setup_system(system,sysmols->protein,sysmols->water,NULL,sysmols->symmetry,settings);

  return(syslist);
}



void output_system(PLI_FILE *file,SYSTEM *system) {

  int free_alist;
  unsigned long int oflags;
  char *oatoms,filename[MAX_LINE_LEN];
  enum OUTPUT_FORMAT oformat;
  enum ATOM_STYLE astyle;
  ATOMLIST *alist,*clist;
  PLI_FILE *ofile;
  SETTINGS *settings;

  settings = system->settings;

  oformat = settings->oformat;
  oflags = settings->oflags;
  astyle = settings->astyle;
  oatoms = settings->oatoms;

  // set atom list:

  oatoms = system->settings->oatoms;

  free_alist = 0;

  if (!strcmp(oatoms,"selection")) {

    alist = system->selection;

    clist = alist;

  } else if (!strcmp(oatoms,"extended")) {

    alist = extended_system_atomlist(system);

    clist = system->selection;

    free_alist = 1;

  } else if (!strcmp(oatoms,"ligand")) {

    alist = system->ligand->selection;

    clist = alist;

  } else if (!strcmp(oatoms,"water")) {

    alist = atomlist2waters(system->selection);

    clist = alist;

    free_alist = 1;

  } else {

    error_fn("output_system: unknown atom selection '%s'",oatoms);
  }

  if (oflags & OUTPUT_SIMB) {

    calc_system_simb(system,0);
  }

  // output system:

  if (oflags & OUTPUT_SYSTEM) {

    write_system(file,system,oformat,oflags);
  }

  // output atoms:

  if (oflags & OUTPUT_ATOMS) {

    write_atom_list(file,alist,astyle,oformat,oflags);
  }

  // output contacts:

  if (oflags & OUTPUT_CONTACTS) {

    write_contacts_atom_list(file,clist,oformat,oflags);
  }

  // output hbonds:

  if (oflags & OUTPUT_HBONDS) {

    write_hbonds_atom_list(file,clist,oformat,oflags);
  }

  // output clashes:

  if (oflags & OUTPUT_CLASHES) {

    write_clashes_atom_list(file,clist,oformat,oflags);
  }

  if (oflags & OUTPUT_SITE_STATS) {

    //write_burial_system(file,system);    
  }

  // output SAS statistics:

  if (oflags & OUTPUT_SAS_STATS) {

    write_sas_stats_atom_list(file,alist,oformat);
  }

  // output files:

  if (oflags & OUTPUT_PDB) {

    sprintf(filename,"%s_%d.pdb",settings->jobname,system->id);

    ofile = open_file(filename,"w");

    write_pdb_atom_list(ofile,alist,oflags);

    close_file(ofile);
  }

  if (oflags & OUTPUT_MOL) {

    sprintf(filename,"%s_%d.mol",settings->jobname,system->id);

    ofile = open_file(filename,"w");

    write_mdl_atom_list(ofile,alist);

    close_file(ofile);
  }

  if (oflags & OUTPUT_LIGAND) {

    check_system_ligand(file,system,oformat);
  }

  if (free_alist) {

    free_atomlist(alist);
  }
}



SYSTEM* molecule2system(MOLECULE *molecule,SETTINGS *settings) {

  SYSTEM *system;

  system = (SYSTEM*) malloc(sizeof(SYSTEM));

  if (system == NULL) {

    error_fn("molecule2system: out of memory allocating system");
  }

  init_system(system);

  system->settings = settings;

  system->molecule_list = alloc_molecule_list();

  init_molecule_list(system->molecule_list);

  add_molecule_to_list(system->molecule_list,molecule);

  system->selection = molecule->selection;

  strcpy(system->name,molecule->name);

  return(system);
}



void prep_system(SYSTEM *system,unsigned int flags) {

  int i,all_molecules_prepped;
  MOLECULE **moleculep,*molecule,*protein,*symmetry;
  MOLECULE_LIST *list;
  SETTINGS *settings;

  list = system->molecule_list;

  if (list == NULL) {

   return;
  }

  settings = system->settings;

  all_molecules_prepped = 1;

  for (i=0,moleculep=list->molecules;i<list->n_molecules;i++,moleculep++) {

    molecule = *moleculep;

    if (molecule->flags & flags) {

      molecule->system = system;

      prep_molecule(molecule,settings); 
    }

    if (!(molecule->prepped)) {

      all_molecules_prepped = 0;
    }
  }

  protein = system->protein;
  symmetry = system->symmetry;

  if ((protein) && (symmetry) && (protein->selection) && (symmetry->selection)) {

    add_atomlist_to_atomlist(protein->selection,symmetry->selection);

    free_atomlist(symmetry->selection);

    symmetry->selection = NULL;
    symmetry->molsystem->selection = NULL;

    protein->molsystem->selection = protein->selection;
  }

  if ((!strcmp(system->name,"system")) && (system->ligand)) {

    strcpy(system->name,system->ligand->name);
  }

  select_system_atoms(system);

  if (all_molecules_prepped) {

    resolve_system(system);

    flag_skipped_waters(system);

    flag_active_waters(system);

  } else {

    for (i=0,moleculep=list->molecules;i<list->n_molecules;i++,moleculep++) {

      molecule = *moleculep;

      if (molecule->prepped) {

	if (molecule->molsystem == NULL) {

	  error_fn("no molsystem for molecule %s\n",molecule->name);
	}

	resolve_system(molecule->molsystem);

	flag_skipped_waters(molecule->molsystem);

	flag_active_waters(molecule->molsystem);
      }
    }
  }

  system->coords = get_list_coords(system->selection);

  if ((all_molecules_prepped) && (system->settings->minimise)) {

    minimise_system(system,settings->sfunc);
  }
}



void unprep_system(SYSTEM *system,unsigned int flags) {

  int i;
  MOLECULE **moleculep,*molecule;
  MOLECULE_LIST *list;
  SETTINGS *settings;

  reset_system(system,1);

  set_list_coords(system->selection,system->coords);

  list = system->molecule_list;

  if (list) {

    settings = system->settings;

    for (i=0,moleculep=list->molecules;i<list->n_molecules;i++,moleculep++) {

      molecule = *moleculep;

      if (molecule->flags & flags) {

	unprep_molecule(molecule,settings);
      }
    }
  }

  unresolve_system(system);

  free_atomlist(system->selection);

  free(system->coords);

  system->selection = NULL;
  system->coords = NULL;
}



void reset_system(SYSTEM *system,int reset_skipped_waters) {

  int i,j;
  ATOM **atomp,*atom;
  ATOMLIST *selection;
  MOLECULE **moleculep,*molecule;
  MOLECULE_LIST *list;
  CONTACTLIST *contactlist;

  list = system->molecule_list;

  if (list == NULL) {

    return;
  }

  // reset atoms:

  for (i=0,moleculep=list->molecules;i<list->n_molecules;i++,moleculep++) {

    molecule = *moleculep;

    for (j=0,atom=molecule->atom;j<molecule->natoms;j++,atom++) {

      contactlist = atom->contactlist;

      if (contactlist) {

	contactlist->ncontacts = 0;
      }

      atom->status = ATOM_NEW;

      atom->cstatus = ATOM_NOTHING_CALCULATED;

      if ((reset_skipped_waters) && (atom->flags & WATER_OXYGEN) && (atom->flags & SKIP_ATOM)) {

	atom->flags ^= SKIP_ATOM;
      }
    }
  }
}



void free_molsystem(SYSTEM *system) {

  if (system) {

    free_molecule_list(system->molecule_list);

    free(system);
  }
}



static void init_system(SYSTEM *system) {

  int i;

  system->id = 0;

  strcpy(system->name,"system");

  system->settings = NULL;

  system->protein = NULL;
  system->water = NULL;
  system->ligand = NULL;
  system->symmetry = NULL;

  system->molecule_list = NULL;

  system->selection = NULL;

  system->score = 0.0;
  system->ref_score = 0.0;
  system->constraint_score = 0.0;

  for (i=0;i<MAX_SYSTEM_SCORES;i++) {

    system->scores[i] = 0.0;
    system->ref_scores[i] = 0.0;
  }

  system->min_dx = 0.0;
  system->min_ds = 0.0;

  system->coords = NULL;
}



static void setup_system(SYSTEM *system,MOLECULE *protein,MOLECULE *water,MOLECULE *ligand,MOLECULE *symmetry,SETTINGS *settings) {
  
  system->settings = settings;

  system->protein = protein;
  system->water = water;
  system->ligand = ligand;
  system->symmetry = symmetry;
    
  system->molecule_list = alloc_molecule_list();

  init_molecule_list(system->molecule_list);

  add_molecule_to_list(system->molecule_list,system->protein);
  add_molecule_to_list(system->molecule_list,system->water);
  add_molecule_to_list(system->molecule_list,system->ligand);
  add_molecule_to_list(system->molecule_list,system->symmetry);
}



static void init_syslist(SYSTEM_LIST *syslist) {

  syslist->n_systems = 0;

  syslist->systems = NULL;
}



void free_system(SYSTEM *system) {

}



static void select_system_atoms(SYSTEM *system) {

  int i;
  MOLECULE **moleculep,*molecule;

  if (system->selection == NULL) {

    system->selection = (ATOMLIST*) alloc_atomlist();

    if (system->selection == NULL) {

      error_fn("select_system_atoms: out of memory allocating selection");
    }
  }

  init_atomlist(system->selection);

  for (i=0,moleculep=system->molecule_list->molecules;i<system->molecule_list->n_molecules;i++,moleculep++) {

    molecule = *moleculep;

    add_atomlist_to_atomlist(system->selection,molecule->selection);
  }
}



static void check_system_ligand(PLI_FILE *file,SYSTEM *system,enum OUTPUT_FORMAT oformat) {

  int i,j,n_atoms,n_amino_atoms,n_total_atoms,n_ligand_chain_atoms,n_protein_chain_atoms;
  char chain;
  ATOM *atom1,*atom2;
  MOLECULE *ligand,*protein;
  ligand = system->ligand;
  protein = system->protein;

  if ((ligand == NULL) || (protein == NULL)) {

    return;
  }

  // identify ligand chain:

  chain = ' ';

  n_amino_atoms = n_total_atoms = n_protein_chain_atoms = n_ligand_chain_atoms = 0;

  for (i=0,atom1=ligand->atom;i<ligand->natoms;i++,atom1++) {

    if (atom1->flags & AMINO_ACID_ATOM) {

      n_atoms = 0;

      for (j=0,atom2=ligand->atom;j<ligand->natoms;j++,atom2++) {

	if (atom2->chain == atom1->chain) {

	  n_atoms++;
	}
      }

      if (n_atoms > n_ligand_chain_atoms) {

	n_ligand_chain_atoms = n_atoms;
	chain = atom1->chain;
      }

      n_amino_atoms++;
    }

    n_total_atoms++;
  }

  if (n_ligand_chain_atoms) {

    n_protein_chain_atoms = 0;

    for (i=0,atom1=protein->atom;i<protein->natoms;i++,atom1++) {

      if ((atom1->flags & AMINO_ACID_ATOM) && (atom1->chain == chain)) {

	n_protein_chain_atoms++;
      }
    }
  }

  printf((oformat == JSON) ? "{\"ligand\":{\"n_aa_ligand\":%d,\"n_aa_chain\":%d}}\n" : "LIGAND %5d %5d\n",n_ligand_chain_atoms,n_protein_chain_atoms); 
}



static void write_system(PLI_FILE *file,SYSTEM *system,enum OUTPUT_FORMAT oformat,unsigned long int oflags) {

  write_line(file,(oformat == JSON) ? 
	     "{\"system\":{\"name\":\"%s\",\"natoms\":%d,\"min_dx\":%.4lf,\"min_ds\":%.4lf" : "system %-20s %5d %10.4lf %10.4lf",
	     system->name,system->selection->natoms,system->min_dx,system->min_ds);

  if (oflags & OUTPUT_SCORES) {

    if (oformat == JSON) {

      write_line(file,",");
    }

    write_system_scores(file,system,oformat,oflags);
  }

  write_line(file,(oformat == JSON) ? "}}\n" : "\n");
}



static ATOMLIST* extended_system_atomlist(SYSTEM *system) {

  int i,j,covalent_contacts,intramolecular_contacts;
  ATOM **atomp,*atom;
  ATOMLIST *list,*selection;
  CONTACT *contact;
  CONTACTLIST *contactlist;
  SETTINGS *settings;

  list = alloc_atomlist();

  init_atomlist(list);

  selection = system->selection;

  settings = system->settings;

  covalent_contacts = (settings->oflags & OUTPUT_COVALENT) ? 1 : 0;
  intramolecular_contacts = (settings->oflags & OUTPUT_INTRAMOLECULAR) ? 1 : 0;

  for (i=0,atomp=selection->atom;i<selection->natoms;i++,atomp++) {

    atom = *atomp;

    if (!atom_in_list(list,atom)) {

      add_atom_to_list(list,atom,0);
    }

    contactlist = atom->contactlist;

    if (contactlist){

      for (j=0,contact=contactlist->contacts;j<contactlist->ncontacts;j++,contact++) {

	if (((covalent_contacts) || (!(contact->flags & COVALENT_CONTACT))) &&
	    ((intramolecular_contacts) || (!(contact->flags & INTRAMOLECULAR_CONTACT)))) {

	  if (!atom_in_list(list,contact->atom2)) {

	    add_atom_to_list(list,contact->atom2,0);
	  }
	}
      }
    }
  }

  return(list);
}
