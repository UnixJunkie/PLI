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



MOLECULE* read_molecule(char *filename) {

  int len,end;
  MOLECULE *molecule;

  len = strlen(filename);

  end = (!strcmp(filename+len-3,".gz")) ? len - 3 : len;

  if ((!strncmp(filename+end-4,".sdf",4)) || (!strncmp(filename+end-4,".mol",4))) {

    molecule = read_mdl_molecule(filename);

  } else if (!strncmp(filename+end-5,".mol2",5)) {

    molecule = read_sybyl_molecule(filename);

  } else if ((!strncmp(filename+end-4,".pdb",4)) || (!strncmp(filename+end-4,".ent",4))) {

    molecule = read_pdb_molecule(filename);

  } else {

    error_fn("read_molecule: unknown file format for file '%s'",filename);
  }

  return(molecule);
}



void write_molecule(MOLECULE *molecule,char *filename) {

  int len,end;

  len = strlen(filename);

  end = (!strcmp(filename+len-3,".gz")) ? len - 3 : len;

  if ((!strncmp(filename+end-4,".sdf",4)) || (!strncmp(filename+end-4,".mol",4))) {

    write_mdl_molecule(molecule,filename);

  } else if (!strncmp(filename+end-5,".mol2",5)) {

    write_sybyl_molecule(molecule,filename);

  } else if ((!strncmp(filename+end-4,".pdb",4)) || (!strncmp(filename+end-4,".ent",4))) {

    write_pdb_molecule(molecule,filename);

  } else {

    error_fn("%s: unknown file format for file '%s'",filename);
  } 
}



void read_molecule_fp(MOLECULE *molecule) {

  int len,end,flag;
  char *filename;

 if (molecule->loaded) {

    return;
  }

  if ((molecule->file == NULL) || (molecule->file_pos == -1)) {

    error_fn("read_molecule_fp: file or file_pos undefined");
  }

  filename = molecule->filename;

  len = strlen(filename);

  end = (!strcmp(filename+len-3,".gz")) ? len - 3 : len;

  flag = 0;

  if ((!strncmp(filename+end-4,".sdf",4)) || (!strncmp(filename+end-4,".mol",4))) {

    flag = read_mdl_molecule_fp(molecule);

  } else {

    error_fn("read_molecule_fp: unknown file format");
  }

  if (flag) {

    error_fn("read_molecule_fp: corrupt molecule read from '%s'",filename);
  }
}



MOLECULE_LIST* read_molecule_list(char *filename) {

  int len,end;
  MOLECULE *molecule;
  MOLECULE_LIST *molecule_list;

  len = strlen(filename);

  end = (!strcmp(filename+len-3,".gz")) ? len - 3 : len;

  if (!strncmp(filename+end-4,".sdf",4)) {

    molecule_list = read_mdl_molecule_list(filename);

  } else {

    molecule_list = alloc_molecule_list();

    init_molecule_list(molecule_list);

    molecule = read_molecule(filename);

    add_molecule_to_list(molecule_list,molecule);
  }

  return(molecule_list);
}



void write_atom_list(PLI_FILE *file,ATOMLIST *list) {

  int i;
  ATOM **atomp;

  if (list == NULL) {

    error_fn("%s: no atoms",__func__);
  }

  for (i=0,atomp=list->atom;i<list->natoms;i++,atomp++) {

    write_atom(file,*atomp);
  }
}



void write_residue(PLI_FILE *file,RESIDUE *residue) {

  char text[MAX_LINE_LEN];

  residue2text(residue,text);

  if (!strcmp((params_get_parameter("residue_oformat"))->value.s,"json")) {

    write_line(file,"{\"residue\":%s}\n",text);

  } else {

    write_line(file,"%s\n",text);
  }
}



void write_atom(PLI_FILE *file,ATOM *atom) {

  char text[MAX_LINE_LEN];

  atom2text(atom,text);

  if (!strcmp((params_get_parameter("atom_oformat"))->value.s,"json")) {

    write_line(file,"{\"atom\":%s}\n",text);

  } else {

    write_line(file,"%s\n",text);
  }
}



void write_bond(PLI_FILE *file,BOND *bond) {

  char text[MAX_LINE_LEN];

  bond2text(bond,text);

  if (!strcmp((params_get_parameter("bond_oformat"))->value.s,"json")) {

    write_line(file,"{\"bond\":%s}\n",text);

  } else {

    write_line(file,"%s\n",text);
  }
}



int duplicate_atom(ATOM *atom1,MOLECULE *molecule) {

  int i;
  double dist;
  ATOM *atom2;

  for (i=0,atom2=molecule->atom;i<molecule->natoms;i++,atom2++) {

    if (points_within_distance(atom1->position,atom2->position,1.0E-6)) {

      return(1);
    }
  }

  return(0);
}
