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



static void read_sybyl_header(PLI_FILE*,MOLECULE*);
static void read_sybyl_atoms(PLI_FILE*,MOLECULE*);
static void read_sybyl_bonds(PLI_FILE*,MOLECULE*);
static void read_sybyl_atom(char*,MOLECULE*,ATOM_TYPING_SCHEME*);
static void read_sybyl_bond(char*,MOLECULE*);
static int find_sybyl_header(PLI_FILE*,char*);



MOLECULE* read_sybyl_molecule(char *filename) {

  char header[MAX_LINE_LEN],line[MAX_LINE_LEN];
  MOLECULE *molecule;
  PLI_FILE *sybylfile;

  sybylfile = open_file(filename,"r");

  if (sybylfile == NULL) {

    error_fn("read_sybyl_molecule: could not open file '%s'",filename);
  }

  molecule = (MOLECULE*) malloc(sizeof(MOLECULE));

  if (molecule == NULL) {

    error_fn("read_pdb_molecule: out of memory allocating molecule");
  }

  init_molecule(molecule,5000);

  strcpy(molecule->filename,filename);

  read_sybyl_header(sybylfile,molecule);

  read_sybyl_atoms(sybylfile,molecule);

  // not sure why this is commented out... should check!

  //read_sybyl_bonds(sybylfile,molecule);

  molecule->loaded = 1;

  return(molecule);
}



void write_sybyl_molecule(MOLECULE *molecule,char *filename) {

  error_fn("%s: not implemented yet",__func__);
}



static void read_sybyl_header(PLI_FILE *sybylfile,MOLECULE *molecule) {

  int n_read;
  char line[MAX_LINE_LEN];

  if (!find_sybyl_header(sybylfile,"@<TRIPOS>MOLECULE")) {

    error_fn("read_sybyl_header: cannot find molecule header");
  }

  if (!read_line(line,MAX_LINE_LEN,sybylfile)) {

    error_fn("read_sybyl_header: corrupt molecule header (1)");
  }

  sscanf(line,"%[^\n]",molecule->name);

  if (!read_line(line,MAX_LINE_LEN,sybylfile)) {

    error_fn("read_sybyl_header: corrupt molecule header (2)");
  }

  n_read = sscanf(line,"%d",&(molecule->natoms));

  if (n_read != 1) {

    error_fn("read_sybyl_header: corrupt molecule header (3)");
  }
}



static void read_sybyl_atoms(PLI_FILE *sybylfile,MOLECULE *molecule) {

  int i,natoms;
  char line[MAX_LINE_LEN];
  ATOM_TYPING_SCHEME *scheme;

  if (!find_sybyl_header(sybylfile,"@<TRIPOS>ATOM")) {

    error_fn("read_sybyl_atoms: cannot find atom header");
  }

  natoms = molecule->natoms;

  molecule->natoms = 0;

  scheme = get_atom_typing_scheme();

  for (i=0;i<natoms;i++) {

    if (!read_line(line,MAX_LINE_LEN,sybylfile)) {

      error_fn("read_sybyl_atoms: not enough atoms");
    }

    read_sybyl_atom(line,molecule,scheme);
  }
}



static void read_sybyl_bonds(PLI_FILE *sybylfile,MOLECULE *molecule) {

  int i,nbonds;
  char line[MAX_LINE_LEN];

  if (!find_sybyl_header(sybylfile,"@<TRIPOS>BOND")) {

    error_fn("read_sybyl_bonds: cannot find bond header");
  }

  nbonds = molecule->nbonds;

  molecule->nbonds = 0;

  for (i=0;i<nbonds;i++) {

    if (!read_line(line,MAX_LINE_LEN,sybylfile)) {

      error_fn("read_sybyl_bonds: not enough bonds");
    }

    read_sybyl_bond(line,molecule);
  }
}



static void read_sybyl_atom(char *line,MOLECULE *molecule,ATOM_TYPING_SCHEME *scheme) {

  int n_words,len,start,end;
  char sybyl_name[10],sybyl_type[10],sybyl_subname[10],elname[5];
  ATOM *atom;

  realloc_atoms(molecule);

  atom = molecule->atom + molecule->natoms;

  init_atom(atom);

  n_words = sscanf(line,"%d %s %lf %lf %lf %s %*d %s",
		   &(atom->id),sybyl_name,
		   &(atom->position[0]),&(atom->position[1]),&(atom->position[2]),
		   sybyl_type,sybyl_subname);

  if (n_words != 7) {

    error_fn("read_sybyl_atom: corrupt atom line\n%s",line);
  }

  // attempt to convert atom name to pdb-style name, so residue/atom-name matching stands a chance:

  strcpy(atom->name,"    ");

  len = strlen(sybyl_name);
  start = (len < 4) ? 1 : 0;
  end = (len < 4) ? 4 : len;

  strncpy(atom->name+start,sybyl_name,len);

  atom->name[end+1] = '\0';

  // similar for the substructure name and number:

  n_words = sscanf(sybyl_subname,"%3s%d",atom->subname,&(atom->subid));

  if (n_words != 2) {

    n_words = sscanf(sybyl_subname,"%[^0-9]%d",atom->subname,&(atom->subid));

    if (n_words != 2) {

      error_fn("read_sybyl_atom: corrupt atom line\n%s",line);
    }
  }

  if (!strcmp(atom->subname,"HID")) {

    strcpy(atom->subname,"HIS");

  } else if (!strcmp(atom->subname,"CYH")) {

    strcpy(atom->subname,"CYS");

  } else if (!strcmp(atom->subname,"LY+")) {

    strcpy(atom->subname,"LYS");
  }

  // get element from sybyl atom type:

  sscanf(sybyl_type,"%[^.]",elname);

  upper_case(elname);

  // skip lone pairs and dummy atoms:

  if ((!strcmp(elname,"LP")) || (!(strcmp(elname,"DU")))) {

    return;
  }

  atom->element = get_element(elname,scheme);

  if (atom->element == NULL) {

    atom->element = pdb_atom_element(atom,scheme);

    if (atom->element == NULL) {

      warning_fn("read_sybyl_atom: skipping ATOM line (cannot derive element '%s'):\n%s",elname,line);

      return;
    }
  }

  // set flag for amino acid atoms:

  if (get_amino_acid(atom->subname) != NULL) {

    atom->flags |= AMINO_ACID_ATOM;
  }

  atom->molecule = molecule;

  molecule->natoms++;
}



static void read_sybyl_bond(char *line,MOLECULE *molecule) {

  int id1,id2,n_words,btype;
  char btype_name[5];
  ATOM *atom1,*atom2;
  BOND *bond;

  n_words = sscanf(line,"%*d %d %d %s",&id1,&id2,&btype_name);

  if (n_words != 3) {

    error_fn("read_sybyl_bond: corrupt bond line\n%s",line);
  }

  if (!strcmp(btype_name,"ar")) {

    btype = 4;

  } else {

    if (!sscanf(btype_name,"%d",&btype)) {

      error_fn("read_sybyl_bond: corrupt bond line\n%s",line);
    }
  }

  if ((btype < 1) || (btype > 4)) {

    error_fn("read_sybyl_bond: bond order is %d",btype);
  }

  atom1 = get_atom(molecule,id1);
  atom2 = get_atom(molecule,id2);

  if ((atom1 == NULL) || (atom2 == NULL)) {

    warning_fn("read_sybyl_bond: skipping bond for line\n%s",line);

    return;
  }

  bond = get_bond(molecule,atom1,atom2);

  if (bond == NULL) {

    add_bond(molecule,atom1,atom2,btype);
  }
}



static int find_sybyl_header(PLI_FILE *sybylfile,char *header) {

  int len;
  char line[MAX_LINE_LEN];

  len = strlen(header);

  do {

    if (!read_line(line,MAX_LINE_LEN,sybylfile)) {

      return(0);
    }

  } while ((strncmp(line,header,len)) && (!end_of_file(sybylfile)));

  return(1);
}
