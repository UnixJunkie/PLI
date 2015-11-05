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



static int get_next_mdl_molecule(MOLECULE*,PLI_FILE*,int,int);
static long int locate_next_mdl_molecule(PLI_FILE*,int);
static void read_mdl_atom(char*,MOLECULE*,ATOM_TYPING_SCHEME*,int);
static void read_mdl_bond(char*,MOLECULE*);
static void read_mdl_water(char*,MOLECULE*,ATOM_TYPING_SCHEME*);
static void write_mdl_atom(PLI_FILE*,ATOM*);



MOLECULE* read_mdl_molecule(char *filename) {

  int flag;
  PLI_FILE *mdlfile;
  MOLECULE *molecule;
  
  mdlfile = open_file(filename,"r");

  if (mdlfile == NULL) {

    error_fn("read_mdl_molecule: could not open file '%s'",filename);
  }

  // allocate memory for molecule and initialise:

  molecule = (MOLECULE*) malloc(sizeof(MOLECULE));

  if (molecule == NULL) {

    error_fn("read_mdl_molecule: out of memory allocating molecule");
  }

  molecule->file = mdlfile;
  molecule->file_pos = 0;
  molecule->lazy_load = 0;

  flag = read_mdl_molecule_fp(molecule);

  if (flag) {

    error_fn("read_mdl_molecule: corrupt file '%s' (error=%d)",filename,flag);
  }

  strcpy(molecule->filename,filename);  

  close_file(mdlfile);

  return(molecule);
}



MOLECULE_LIST* read_mdl_molecule_list(char *filename) {

  int i,j,n_molecules,lazy_load;
  int n_alloc_molecules = 100;
  PLI_FILE *mdlfile;
  ATOM *atom;
  MOLECULE *molecule,**moleculep,*molecules,*test;
  MOLECULE_LIST *list;
  
  lazy_load = 1;

  mdlfile = open_file(filename,"r");

  if (mdlfile == NULL) {

    error_fn("read_mdl_molecule: could not open file '%s'",filename);
  }

  list = alloc_molecule_list();

  init_molecule_list(list);

  molecules = (MOLECULE*) calloc(n_alloc_molecules,sizeof(MOLECULE));

  if (!molecules) {

    error_fn("read_mdl_molecule_list: out of memory allocating molecules");
  }

  n_molecules = 0;

  while (!get_next_mdl_molecule(molecules+n_molecules,mdlfile,lazy_load,n_molecules)) {

    if (n_molecules == n_alloc_molecules-1) {

      n_alloc_molecules *= 2;

      molecules = (MOLECULE*) realloc(molecules,(n_alloc_molecules)*sizeof(MOLECULE));
    }

    n_molecules++;

    if (!molecules) {

      error_fn("read_mdl_molecule_list: out of memory reallocating (%d) molecules",n_alloc_molecules);
    }
  }

  list->n_molecules = n_molecules;

  if (n_molecules) {

    if (n_molecules < n_alloc_molecules) {

      molecules = (MOLECULE*) realloc(molecules,(n_molecules)*sizeof(MOLECULE));
    }

    list->molecules = (MOLECULE**) calloc(n_molecules,sizeof(MOLECULE*));

    if (!list->molecules) {

      error_fn("read_mdl_molecule_list: out of memory allocating molecule list");
    }

    for (i=0,moleculep=list->molecules,molecule=molecules;i<n_molecules;i++,moleculep++,molecule++) {
      
      *moleculep = molecule;
    }
    
  } else {

    free(molecules);
  }

  if (lazy_load) {

    list->file = mdlfile;

  } else {

    close_file(mdlfile);
  }

  return(list);
}



static int get_next_mdl_molecule(MOLECULE *molecule,PLI_FILE *mdlfile,int lazy_load,int id) {

  int flags = 0;

  if (lazy_load) {

    init_molecule(molecule,0);

    molecule->file = mdlfile;

    molecule->lazy_load = 1;

    molecule->file_pos = locate_next_mdl_molecule(mdlfile,id);

    if (molecule->file_pos == -1) {

      flags = 1;
    }

  } else {

    molecule->file = mdlfile;

    flags = read_mdl_molecule_fp(molecule);
  }

  strcpy(molecule->filename,mdlfile->filename);

  // always rely on the bonds specified in the file:

   molecule->rely_on_bonds = 1;

  return(flags);
}



static long int locate_next_mdl_molecule(PLI_FILE *mdlfile,int id) {

  long int file_pos;
  char line[MAX_LINE_LEN];

  if (id == 0) {

    return(0);
  }

  do {

    if (!read_line(line,MAX_LINE_LEN,mdlfile)) {

      return(-1);
    }
    
  } while (strncmp(line,"$$$$",4));

  file_pos = pli_ftell(mdlfile);

  if (!read_line(line,MAX_LINE_LEN,mdlfile)) {

    return(-1);
  }
  
  return(file_pos);
}



int read_mdl_molecule_fp(MOLECULE *molecule) {

  int i,natoms,nbonds;
  char line[MAX_LINE_LEN],header[MAX_LINE_LEN];
  ATOM_TYPING_SCHEME *scheme;
  PLI_FILE *mdlfile;
  
  mdlfile = molecule->file;

  if (molecule->file_pos != -1) {

    pli_fseek(mdlfile,molecule->file_pos,SEEK_SET);
  }

  if (!(molecule->lazy_load)) {

    init_molecule(molecule,50);
  }

  // read molecule name:

  if (!read_line(line,MAX_LINE_LEN,mdlfile)) {

    return(1);
  }
  
  sscanf(line,"%[^\n]",molecule->name);

  // skip two next lines:

  for (i=0;i<2;i++) {

    if (!read_line(line,MAX_LINE_LEN,mdlfile)) {

      return(2);
    }
  }
  
  // read number of atoms and bonds:

  if (!read_line(line,MAX_LINE_LEN,mdlfile)) {

    return(3);
  }

  read_word(line,"%3d",&natoms);
  read_word(line+3,"%3d",&nbonds);

  // read atoms:

  scheme = get_atom_typing_scheme();

  for (i=0;i<natoms;i++) {

    if (!read_line(line,MAX_LINE_LEN,mdlfile)) {

      return(4);
    }

    read_mdl_atom(line,molecule,scheme,i+1);
  }

  // read bonds:

  for (i=0;i<nbonds;i++) {

    if (!read_line(line,MAX_LINE_LEN,mdlfile)) {

      return(5);
    }

    read_mdl_bond(line,molecule);
  }

  // read meta data:

  strcpy(header,"unknown");

  do {

    if (!read_line(line,MAX_LINE_LEN,mdlfile)) {

      return(0);
    }

    if (line[0] == '>') {

      sscanf(line,"%*s <%[^>]",header);

    } else if (!strcmp(header,"Gold.Protein.RotatedWaterAtoms")) {

      read_mdl_water(line,molecule,scheme);
    }

  } while (strncmp(line,"$$$$",4));
  
  if (molecule->natoms) {

    if (molecule->natoms < molecule->n_alloc_atoms) {

      molecule->atom = (ATOM*) realloc(molecule->atom,(molecule->natoms)*sizeof(ATOM));

      molecule->n_alloc_atoms = molecule->natoms;
    }
  }
  
  if (molecule->nbonds) {

    if (molecule->nbonds < molecule->n_alloc_bonds) {

      molecule->bond = (BOND*) realloc(molecule->bond,(molecule->nbonds)*sizeof(BOND));

      molecule->n_alloc_bonds = molecule->nbonds;
    }
  }

  molecule->loaded = 1;

  return(0);
}



static void read_mdl_atom(char *line,MOLECULE *molecule,ATOM_TYPING_SCHEME *scheme,int id) {

  int n_words;
  char elname[3];
  ATOM *atom;

  realloc_atoms(molecule);

  atom = molecule->atom + molecule->natoms;

  init_atom(atom);

  n_words = sscanf(line,"%lf %lf %lf %s",&(atom->position[0]),&(atom->position[1]),&(atom->position[2]),elname);

  if (n_words != 4) {

    error_fn("read_mdl_atom: corrupt atom line");
  }

  upper_case(elname);

  atom->id = id;

  atom->position[3] = 1.0;

  copy_vector(atom->position,atom->original_position);

  sprintf(atom->name,"%s%d",elname,atom->id);

  atom->element = get_element(elname,scheme);

  if (atom->element == NULL) {

    error_fn("read_mdl_atom: unknown element '%s' in molecule '%s'\nfilename: %s\nline: %s",elname,molecule->name,molecule->filename,line);
  }

  atom->molecule = molecule;

  atom->seqid = molecule->natoms;
  molecule->natoms++;
}



static void read_mdl_bond(char *line,MOLECULE *molecule) {

  int id1,id2,btype,n_words;
  ATOM *atom1,*atom2;
  BOND *bond;

  n_words = sscanf(line,"%d %d %d",&id1,&id2,&btype);

  if (n_words != 3) {

    error_fn("read_mdl_bond: corrupt bond line");
  }

  if ((btype < 1) || (btype > 3)) {

    error_fn("read_mdl_bond: bond order is %d",btype);
  }

  atom1 = get_atom(molecule,id1);
  atom2 = get_atom(molecule,id2);

  if ((atom1 == NULL) || (atom2 == NULL)) {

    error_fn("read_mdl_bond: atom(s) not found");
  }

  bond = get_bond(molecule,atom1,atom2);

  if (bond == NULL) {

    add_bond(molecule,atom1,atom2,btype);
  }
}



static void read_mdl_water(char *line,MOLECULE *molecule,ATOM_TYPING_SCHEME *scheme) {

  int n_read,id;
  double x,y,z;
  char type[20],toggle[5];
  ATOM *water;

  n_read = sscanf(line,"%lf %lf %lf %s %*[^#]# atno %d %*[^O]%3s",&x,&y,&z,type,&id,&toggle); 

  if (n_read == 6) {

    if (!strcmp(type,"O.3")) {

      if (molecule->active_waters == NULL) {

	molecule->n_alloc_active_waters = 5;

	molecule->active_waters = (ATOM*) calloc(5,sizeof(ATOM));

      } else if (molecule->n_active_waters == molecule->n_alloc_active_waters) {

	molecule->n_alloc_active_waters *= 2;

	molecule->active_waters = (ATOM*) realloc(molecule->active_waters,(molecule->n_alloc_active_waters)*sizeof(ATOM));
      }

      if (molecule->active_waters == NULL) {

	error_fn("read_mdl_water: out of memory (re)allocating active waters");
      }

      water = molecule->active_waters + molecule->n_active_waters;

      init_atom(water);

      water->id = id;

      water->type = get_atom_type("H2O",scheme);

      water->position[0] = x;
      water->position[1] = y;
      water->position[2] = z;

      if (!strcmp(toggle,"OFF")) {

	water->flags |= SKIP_ATOM;
      }

      molecule->n_active_waters++;
    }
  }
}



void write_mdl_atom_list(PLI_FILE *file,ATOMLIST *atomlist) {

  int i,id1,id2;
  ATOM **atom;
  BOND **bond;
  BONDLIST *bondlist;

  bondlist = atomlist2bondlist(atomlist);

  write_line(file,"%s\n","pli mdl mol file");
  write_line(file,"\n");
  write_line(file,"\n");
  write_line(file,"%3d%3d  0     0  0              0 V2000\n",atomlist->natoms,bondlist->nbonds);

  for (i=0,atom=atomlist->atom;i<atomlist->natoms;i++,atom++) {

    write_mdl_atom(file,*atom);
  }

  for (i=0,bond=bondlist->bond;i<bondlist->nbonds;i++,bond++) {

    id1 = atomlist_index(atomlist,(*bond)->atom1);
    id2 = atomlist_index(atomlist,(*bond)->atom2);

    if ((id1 != -1) && (id2 != -2)) {

      write_line(file,"%3d%3d%3d  0  0  0  0\n",id1+1,id2+1,(*bond)->type);

    } else {

      warning_fn("write_mdl_atom_list: cannot find atom(s) in list");
    }
  }

  write_line(file,"M  END\n");
  write_line(file,"$$$$\n");

  free_bondlist(bondlist);
}



static void write_mdl_atom(PLI_FILE *file,ATOM *atom) {

  int charge = 0;

  if ((atom->type) && (atom->type->united_atom)) {
      
    if (atom->type->united_atom->formal_charge != 0) {

      charge = 4 - atom->type->united_atom->formal_charge;
    }

  } else {

    warning_fn("write_mdl_atom: atom type not set for atom '%s'",atom->name);
  }

  write_line(file,"%10.4lf%10.4lf%10.4lf %-2s  0%3d  0  0  0  0  0  0  0  0  0  0\n",
	     atom->position[0],atom->position[1],atom->position[2],
	     atom->element->name,charge);
}
