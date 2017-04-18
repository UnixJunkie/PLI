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



static int get_next_mdl_molecule(MOLECULE*,PLI_FILE*,int,int);
static long int locate_next_mdl_molecule(PLI_FILE*,int);
static void read_mdl_atom(char*,MOLECULE*,ATOM_TYPING_SCHEME*,int);
static void read_mdl_bond(char*,MOLECULE*);
static void read_mdl_water(char*,MOLECULE*,ATOM_TYPING_SCHEME*);
static void write_mdl_atom(PLI_FILE*,ATOM*,ATOM_TYPING_SCHEME*);



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



void write_mdl_molecule(MOLECULE *molecule,char *filename) {

  error_fn("%s: not implemented yet",__func__);
}



void write_mdl_atom_list(PLI_FILE *file,ATOMLIST *atomlist, char *mol_name, int write_hydrogens, int n_properties, ...) {

  int i,j,id1,id2;
  ATOM **atom;
  BOND **bond;
  BONDLIST *bondlist;
  MOLECULE_PROPERTY *mp;
  ATOM_TYPING_SCHEME *scheme;

  // create a new list of output atoms as connected hydrogens will be added to it
  ATOMLIST *output_atoms;

  if (write_hydrogens) {
    output_atoms = alloc_atomlist();
    init_atomlist(output_atoms);
    copy_atomlist(atomlist, output_atoms);
    for (i=0; i<atomlist->natoms; i++){
      ATOM *a = atomlist->atom[i];
      for (j=0; j<a->h_connections->natoms; j++) {
        ATOM *ha = a->h_connections->atom[j];
        if (!atom_in_list(output_atoms, ha)) {
          add_atom_to_list(output_atoms, ha, ha->flags);
        }
      }
    }
  } else {
    output_atoms = atomlist;
  }

  bondlist = atomlist2bondlist(output_atoms);

  //write_line(file,"%s\n","pli mdl mol file");
  write_line(file,"%s\n", mol_name);
  write_line(file,"\n");
  write_line(file,"\n");
  write_line(file,"%3d%3d  0     0  0              0 V2000\n",output_atoms->natoms,bondlist->nbonds);

  scheme = get_atom_typing_scheme();

  for (i=0,atom=output_atoms->atom;i<output_atoms->natoms;i++,atom++) {
    write_mdl_atom(file,*atom,scheme);
  }

  //bonds
  for (i=0,bond=bondlist->bond;i<bondlist->nbonds;i++,bond++) {

    id1 = atomlist_index(output_atoms,(*bond)->atom1);
    id2 = atomlist_index(output_atoms,(*bond)->atom2);

    if ((id1 != -1) && (id2 != -2)) {

      write_line(file,"%3d%3d%3d  0  0  0  0\n",id1+1,id2+1,(*bond)->type);

    } else {

      warning_fn("write_mdl_atom_list: cannot find atom(s) in list");
    }
  }


  write_line(file,"M  END\n");

  if (n_properties > 0) {
    va_list arguments;
    va_start(arguments, n_properties);
    for (i=0; i<n_properties; i++) {
      write_line(file, "\n");
      mp = va_arg(arguments, MOLECULE_PROPERTY*);
      write_line(file, ">  <%s>\n", mp->name);
      if (mp->type == 'i') {
	write_line(file, "%d\n", mp->value.i);
      }
      if (mp->type == 'd') {
	write_line(file, "%.2f\n", mp->value.d);
      }
      if (mp->type == 's') {
	write_line(file, "%s\n", mp->value.s);
      }
      write_line(file, "\n");
    }
  }

  write_line(file,"$$$$\n");

  free_bondlist(bondlist);
  if (output_atoms != atomlist) {
    free_atomlist(output_atoms);
  }
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
  
  read_word(line,"%[^\n]",molecule->name);

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

  if (!read_word(line,"%3d",&natoms)) {

    return(3);
  }

  if (!read_word(line+3,"%3d",&nbonds)) {

    return(3);
  }

  // read atoms:

  scheme = get_atom_typing_scheme();

  for (i=0;i<natoms;i++) {

    if (!read_line(line,MAX_LINE_LEN,mdlfile)) {

      return(4);
    }

    read_mdl_atom(line,molecule,scheme,i+1);
  }
  
  if (molecule->natoms) {

    if (molecule->natoms < molecule->n_alloc_atoms) {

      molecule->atom = (ATOM*) realloc(molecule->atom,(molecule->natoms)*sizeof(ATOM));

      molecule->n_alloc_atoms = molecule->natoms;
    }
  }

  // read bonds:

  for (i=0;i<nbonds;i++) {

    if (!read_line(line,MAX_LINE_LEN,mdlfile)) {

      return(5);
    }

    read_mdl_bond(line,molecule);
  }
  
  if (molecule->nbonds) {

    if (molecule->nbonds < molecule->n_alloc_bonds) {

      molecule->bond = (BOND*) realloc(molecule->bond,(molecule->nbonds)*sizeof(BOND));

      molecule->n_alloc_bonds = molecule->nbonds;
    }
  }

  // read meta data:

  strcpy(header,"unknown");

  do {

    if (!read_line(line,MAX_LINE_LEN,mdlfile)) {

      //
      molecule->loaded = 1;
      //
      return(0);
    }

    if (line[0] == '>') {

      sscanf(line,"%*s <%[^>]",header);

    } else if (!strcmp(header,"Gold.Protein.RotatedWaterAtoms")) {

      read_mdl_water(line,molecule,scheme);
    }

  } while (strncmp(line,"$$$$",4));

  molecule->loaded = 1;
  molecule->rely_on_bonds = 1;

  return(0);
}



static void read_mdl_atom(char *line,MOLECULE *molecule,ATOM_TYPING_SCHEME *scheme,int id) {

  int n_words;
  char elname[3];
  ATOM *atom;

  realloc_atoms(molecule);

  atom = molecule->atom + molecule->natoms;

  init_atom(atom);
  // "%10.4lf%10.4lf%10.4lf %-2s  0%3d  0  0  0  0  0  %3d
  n_words = sscanf(line,"%10lf%10lf%10lf%*[ ]%s%*d%*3c%*3c%*3c%*3c%*3c%*3c%3d",
                   &(atom->position[0]),&(atom->position[1]),&(atom->position[2]),elname, &(atom->type_id_file));
  if (n_words < 4) {

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

  n_words = read_word(line,"%3d",&id1);
  n_words += read_word(line+3,"%3d",&id2);
  n_words += read_word(line+6,"%3d",&btype);

  if (n_words != 3) {

    error_fn("read_mdl_bond: corrupt bond line");
  }

  if ((btype < 1) || (btype > 3)) {

    error_fn("read_mdl_bond: bond order is %d",btype);
  }

  atom1 = get_atom(molecule,id1);
  atom2 = get_atom(molecule,id2);

  if ((atom1 == NULL) || (atom2 == NULL)) {

    error_fn("read_mdl_bond: atom(s) not found (%d,%d)\n%s\n",id1,id2,line);
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



static void write_mdl_atom(PLI_FILE *file,ATOM *atom,ATOM_TYPING_SCHEME *scheme) {

  int charge = 0;
  char elem_name[10];
  ATOM *carbon;
  ATOMLIST *conns;
  BOND *bond;
  int atom_type = 0;

  if ((atom->type) && (atom->type->united_atom)) {
      
    atom_type = atom->type->id;
    if (atom->type->united_atom->formal_charge != 0) {

      charge = 4 - atom->type->united_atom->formal_charge;
    }

    if (atom->type == get_atom_type("Carboxyl =O",scheme)) {

      // sort out carboxyl groups:

      if (atom->connections) {

	conns = atom->connections;

	if (conns->natoms == 1) {

	  carbon = conns->atom[0];

	  bond = get_bond(atom->molecule,atom,carbon);

	  if (bond->type == 1) {

	    charge = 5;
	  }
	}
      }

    } else if ((atom->type == get_atom_type("Amidine =NH2",scheme)) ||
	       (atom->type == get_atom_type("=NH2 ARG-NH",scheme))) {

      // sort out amidines:

      if (atom->connections) {

	conns = atom->connections;

	if (conns->natoms == 1) {

	  carbon = conns->atom[0];

	  bond = get_bond(atom->molecule,atom,carbon);

	  if (bond->type == 1) {

	    charge = 0;
	  }
	}
      }
    }
  } else if (atom->type->element->id != HYDROGEN) {

    warning_fn("write_mdl_atom: atom type not set for atom '%s'",atom->name);
  }

  strcpy(elem_name, atom->element->name);
  title_case(elem_name);

  write_line(file,"%10.4lf%10.4lf%10.4lf %-2s  0%3d  0  0  0  0  0%3d  0  0  0  0\n",
	     atom->position[0],atom->position[1],atom->position[2],
	     elem_name,charge,atom_type);
}
