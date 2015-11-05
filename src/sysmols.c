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


static void sysmols_init_complex(COMPLEX*);
static void init_sysmols(SYSMOLS*);
static void setup_sysmol(MOLECULE*,char*,unsigned int,SETTINGS*);
static void fix_symmetry_atoms(MOLECULE*,MOLECULE*);
static COMPLEX_LIST* sysmols_read_complex_list(char*);
static sysmols_read_complex_list_molecules(COMPLEX_LIST*,SETTINGS*);



SYSMOLS* load_sysmols(SETTINGS *settings) {

  int i;
  unsigned int flags;
  ATOM *atom;
  SYSMOLS *sysmols;
  MOLECULE *protein,**ligandp,*symmetry;

  // allocate memory for molecules:

  sysmols = (SYSMOLS*) malloc(sizeof(SYSMOLS));

  if (sysmols == NULL) {

    error_fn("load_sysmols: out of memory allocating sysmols");
  }

  init_sysmols(sysmols);

  // read protein:

  if (strcmp(settings->protein_file,"undefined")) {

    sysmols->protein = read_molecule(settings->protein_file);

    setup_sysmol(sysmols->protein,"protein",PROTEIN_MOLECULE,settings);
  }

  // read ligand(s):

  if (strcmp(settings->ligand_file,"undefined")) {

    sysmols->ligand_list = read_molecule_list(settings->ligand_file);

    for (i=0,ligandp=sysmols->ligand_list->molecules;i<sysmols->ligand_list->n_molecules;i++,ligandp++) {

      setup_sysmol(*ligandp,"ligand",LIGAND_MOLECULE,settings);
    }

    sysmols->ligand = *(sysmols->ligand_list->molecules);
  }

  // read symmetry:

  if (strcmp(settings->symmetry_file,"undefined")) {

    flags = PROTEIN_MOLECULE;
    flags |= SYMMETRY_MOLECULE;

    sysmols->symmetry = read_molecule(settings->symmetry_file);
  }

  if ((sysmols->protein) && (sysmols->symmetry)) {

    fix_symmetry_atoms(sysmols->symmetry,sysmols->protein);
  }

  if (sysmols->symmetry) {

    setup_sysmol(sysmols->symmetry,"symmetry",flags,settings);
  }

  // read complex list:

  if (strcmp((params_get_parameter("complexes"))->value.s,"undefined")) {

    if ((sysmols->ligand_list) || (sysmols->protein) || (sysmols->symmetry)) {

      error_fn("load_sysmols: cannot combine complexes with other molecule files");
    }

    sysmols->complex_list = sysmols_read_complex_list((params_get_parameter("complexes"))->value.s);

    sysmols_read_complex_list_molecules(sysmols->complex_list,settings);
  }

  return(sysmols);
}



static sysmols_read_complex_list_molecules(COMPLEX_LIST *complex_list,SETTINGS *settings) {

  int i;
  COMPLEX *complex;

  if (complex_list == NULL) {

    return;
  }

  for (i=0,complex=complex_list->complexes;i<complex_list->n_complexes;i++,complex++) {

    complex->protein = read_molecule(complex->protein_file);
    complex->ligand = read_molecule(complex->ligand_file);

    setup_sysmol(complex->protein,"protein",PROTEIN_MOLECULE,settings);
    setup_sysmol(complex->ligand,"ligand",LIGAND_MOLECULE,settings);
  }
}



static COMPLEX_LIST* sysmols_read_complex_list(char *filename) {

  int n_read,n_alloc_complexes,n_complexes;
  char line[MAX_LINE_LEN];
  PLI_FILE *file;
  COMPLEX *complex,*complexes;
  COMPLEX_LIST *complex_list;

  file = open_file(filename,"r");

  if (file == NULL) {

    error_fn("sysmols_read_complex_list: could not open complex list file '%s'",filename); 
  }

  n_alloc_complexes = 10;

  complexes = (COMPLEX*) calloc(n_alloc_complexes,sizeof(COMPLEX));

  if (complexes == NULL) {

    error_fn("sysmols_read_complex_list: out of memory allocating complexes");
  }

  n_complexes = 0;

  complex = complexes;

  while (!end_of_file(file)) {
	      
    if (read_line(line,MAX_LINE_LEN,file) == NULL)
      break;
 
    if (line[0] != '#') {

      sysmols_init_complex(complex);

      n_read = sscanf(line,"%s %s",complex->protein_file,complex->ligand_file);

      if (n_read == 2) {

	n_complexes++;

	if (n_complexes == n_alloc_complexes) {

	  n_alloc_complexes *= 2;

	  complexes = (COMPLEX*) realloc(complexes,n_alloc_complexes*sizeof(COMPLEX));

	  if (complexes == NULL) {

	    error_fn("sysmols_read_complex_list: out of memory reallocating complexes");
	  }
	}

	complex = complexes + n_complexes;
      }
    }
  }

  close_file(file);

  complexes = (COMPLEX*) realloc(complexes,n_complexes*sizeof(COMPLEX));

  complex_list = (COMPLEX_LIST*) malloc(sizeof(COMPLEX_LIST));

  if (complex_list == NULL) {

    error_fn("sysmols_read_complex_list: out of memory allocating complex list");
  }

  complex_list->n_alloc_complexes = n_alloc_complexes;
  complex_list->n_complexes = n_complexes;

  complex_list->complexes = complexes;

  return(complex_list);
}



static void sysmols_init_complex(COMPLEX *complex) {

  strcpy(complex->protein_file,"undefined");
  strcpy(complex->ligand_file,"undefined");

  complex->protein = NULL;
  complex->ligand = NULL;
}



static void init_sysmols(SYSMOLS *sysmols) {

  sysmols->protein = NULL;
  sysmols->water = NULL;
  sysmols->ligand = NULL;
  sysmols->symmetry = NULL;
  sysmols->ligand_list = NULL;
  sysmols->complex_list = NULL;
}



static void setup_sysmol(MOLECULE *molecule,char *name,unsigned int flags,SETTINGS *settings) {

  strcpy(molecule->name,name);

  molecule->flags = flags;

  if (flags & PROTEIN_MOLECULE) {

    molecule->use_pdb_dict = settings->protein_use_pdb_dict;
    molecule->use_hydrogens = settings->protein_use_hydrogens;
    molecule->use_grids = 1;

  } else if (flags & LIGAND_MOLECULE) {

    molecule->use_pdb_dict = settings->ligand_use_pdb_dict;
    molecule->use_hydrogens = settings->ligand_use_hydrogens;
    molecule->use_grids = 0;

  } else {

    warning_fn("setup_sysmol: unkown molecule type for molecule '%s'",molecule->name);
  }
}



static void fix_symmetry_atoms(MOLECULE *symmetry,MOLECULE *protein) {

  int i,j,duplicate;
  ATOM *satom,*patom;

  i = 0;

  while (i < symmetry->natoms) {

    satom = symmetry->atom + i;

    j = 0;
    duplicate = 0;

    while ((j < protein->natoms) && (!duplicate)) {

      patom = protein->atom + j;

      if (points_within_distance(satom->position,patom->position,1.0E-6)) {

	duplicate = 1;
      }

      j++;
    }

    if (duplicate) {

      warning_fn("fix_symmetry_atoms: deleting duplicate symmetry atom");

      delete_atom(symmetry,i);

    } else {

      i++;
    }
  }
}
