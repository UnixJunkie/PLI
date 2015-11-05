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



static int N_PDB_DICTS = 0;
static MOLECULE **PDB_DICTS = NULL;

static int N_PDB_DICT_FAILURES = 0;
static MOLECULE *PDB_DICT_FAILURES = NULL;



static void read_pdb_atom(char*,MOLECULE*,ATOM_TYPING_SCHEME*);
static void read_pdb_bonds(char*,MOLECULE*);
static void read_cif_atom(char*,MOLECULE*,ATOM_TYPING_SCHEME*);
static void read_cif_bond(char*,MOLECULE*);
static void write_pdb_atom_axes(PLI_FILE*,ATOM*,ATOM*);
static void write_pdb_atom_vpts(PLI_FILE*,ATOM*,ATOM*);



MOLECULE* read_pdb_molecule(char *filename) {

  int i,n_models;
  char line[MAX_LINE_LEN];
  PLI_FILE *pdbfile;
  ATOM_TYPING_SCHEME *scheme;
  ATOM *atom;
  MOLECULE *molecule;

  pdbfile = open_file(filename,"r");

  if (pdbfile == NULL) {

    error_fn("read_pdb_molecule: could not open file '%s'",filename);
  }

  scheme = get_atom_typing_scheme();

  molecule = (MOLECULE*) malloc(sizeof(MOLECULE));

  if (molecule == NULL) {

    error_fn("read_pdb_molecule: out of memory allocating molecule");
  }

  init_molecule(molecule,5000);

  strcpy(molecule->filename,filename);

  n_models = 0;

  while (!end_of_file(pdbfile)) {

    if (read_line(line,MAX_LINE_LEN,pdbfile) == NULL)
      break;

    if (!strncmp(line,"MODEL",5)) {

      if (n_models == 1) {

	warning_fn("read_pdb_molecule: multi-model structure; only reading first model");

	break;
      }

      n_models++;
    }

    if ((!strncmp(line,"ATOM",4)) || (!strncmp(line,"HETATM",6)))
      read_pdb_atom(line,molecule,scheme);
    
    if (!strncmp(line,"CONECT",6))  
      read_pdb_bonds(line,molecule);
  }

  for (i=0,atom=molecule->atom;i<molecule->natoms;i++,atom++) {

    if (atom->element == NULL) {

      atom->element = pdb_atom_element(atom,scheme);
    }

    if ((atom->element != NULL) && (atom->element->id == OXYGEN)) {

      if ((!strcmp(atom->subname,"HOH")) || (!strcmp(atom->subname,"WAT"))) {

	atom->flags |= WATER_OXYGEN;
      }
    }
  }

  close_file(pdbfile);

  molecule->loaded = 1;

  return(molecule);
}



static void read_pdb_atom(char *line,MOLECULE *molecule,ATOM_TYPING_SCHEME *scheme) {

  int len;
  char elname[3];
  ATOM *atom;

  len = strlen(line);

  if (len < 66) {

    warning_fn("read_pdb_atom: skipping corrupt ATOM line:\n%s",line);

    return;
  }

  realloc_atoms(molecule);

  atom = molecule->atom + molecule->natoms;

  init_atom(atom);

  sscanf(line+6,"%5d",&atom->id);
  sscanf(line+16,"%c",&atom->altloc);
  sscanf(line+21,"%c",&atom->chain);
  sscanf(line+22,"%4d",&atom->subid);
  sscanf(line+26,"%c",&atom->icode);

  sscanf(line+30,"%8lf%8lf%8lf",&atom->position[0],&atom->position[1],&atom->position[2]);

  atom->position[3] = 1.0;

  copy_vector(atom->position,atom->original_position);

  sscanf(line+54,"%6lf",&atom->occupancy);
  sscanf(line+60,"%6lf",&atom->bfactor);

  substring(line,12,4,atom->name);
  substring(line,17,3,atom->subname);

  if (get_amino_acid(atom->subname) != NULL) {

    atom->flags |= AMINO_ACID_ATOM;
  }

  if (len >= 78) {

    if (line[76] == ' ') {

      substring(line,77,1,elname);

    } else {

      substring(line,76,2,elname);
    }

    atom->element = get_element(elname,scheme);
  }

  if (atom->element == NULL) {

    atom->element = pdb_atom_element(atom,scheme);
  }

  if (atom->element == NULL) {

    warning_fn("read_pdb_atom: skipping ATOM line (cannot derive element):\n%s",line);

    return;
  }

  if (duplicate_atom(atom,molecule)) {

    warning_fn("read_pdb_atom: skipping ATOM line (duplicate atom):\n%s",line);

    return;
  }

  atom->molecule = molecule;

  atom->seqid = molecule->natoms;

  molecule->natoms++;
}



static void read_pdb_bonds(char *line,MOLECULE *molecule) {

  int i,id1,id2;
  int len;
  char idstr[6];
  ATOM *atom,*atom1,*atom2;
  BOND *bond;

  idstr[5] = '\0';

  strncpy(idstr,line+6,5);

  id1 = atoi(idstr);

  atom1 = get_atom(molecule,id1);

  if (atom1 == NULL) {

    warning_fn("read_pdb_bonds: could not find atom %d (%s)",id1,molecule->filename);

    return;
  }

  len = strlen(line) - 12;

  i = 0;

  while ((i+1)*5 <= len) {

    strncpy(idstr,line+11+5*i,5);

    if (strcmp(idstr,"     ")) {

      id2 = atoi(idstr);

      atom2 = get_atom(molecule,id2);

      if (atom2 != NULL) {

	bond = get_bond(molecule,atom1,atom2);

	if (bond == NULL) {

	  add_bond(molecule,atom1,atom2,1);

	} else if ((bond->atom1 == atom1) && (bond->atom2 == atom2)) {

	  bond->type += 1;
	}

      } else {

	warning_fn("read_pdb_bonds: could not find atom %d (%s)",id2,molecule->filename);
      }
    }

    i++;
  }
}



MOLECULE* get_pdb_dict(char *name) {

  int i,dirlen;
  char *pli_dir,code[4],dir[3],filename[MAX_LINE_LEN];
  MOLECULE **dict,*failure,*new_dict;

  if (!strcmp(name,"unknown")) {

    return(NULL);
  }

  // find dictionary in list:

  for (i=0,dict=PDB_DICTS;i<N_PDB_DICTS;i++,dict++) {

    if (!strcmp((*dict)->name,name)) {

      return(*dict);
    }
  }

  // check previous failures:

  for (i=0,failure=PDB_DICT_FAILURES;i<N_PDB_DICT_FAILURES;i++,failure++) {
  
    if (!strcmp(failure->name,name)) {

      return(NULL);
    }
  }

  // try to open a new dictionary:

  pli_dir = get_pli_dir();

  if (strlen(name) > 3) {

    warning_fn("get_pdb_dict: substructure name '%s' skipped",name);

    return(NULL);
  }

  remove_spaces(name,code);

  dirlen = (strlen(name) < 2) ? 1 : 2;

  substring(code,0,dirlen,dir);

  sprintf(filename,"%s/pdbdict/cif/%s/%s.cif",pli_dir,dir,code);

  new_dict = read_cif_molecule(filename);

  if (new_dict == NULL) {

    warning_fn("get_pdb_dict: no dictionary file for '%s' in '%s'",code,filename);

    if (N_PDB_DICT_FAILURES == 0) {

      PDB_DICT_FAILURES = (MOLECULE*) calloc(1,sizeof(MOLECULE));

    } else {

      PDB_DICT_FAILURES = realloc(PDB_DICT_FAILURES,(N_PDB_DICT_FAILURES+1)*sizeof(MOLECULE));
    }

    if (PDB_DICT_FAILURES == NULL) {

      error_fn("get_pdb_dict: out of memory allocation %d pdb dictionary failures",N_PDB_DICT_FAILURES+1);
    }

    strcpy(PDB_DICT_FAILURES[N_PDB_DICT_FAILURES].name,code);

    N_PDB_DICT_FAILURES++;

    return(NULL);
  }

  strcpy(new_dict->name,name);

  // add new dictionary to list:

  if (N_PDB_DICTS == 0) {

    PDB_DICTS = (MOLECULE**) calloc(1,sizeof(MOLECULE*));

  } else {

    PDB_DICTS = (MOLECULE**) realloc(PDB_DICTS,(N_PDB_DICTS+1)*sizeof(MOLECULE*));
  }

  if (PDB_DICTS == NULL) {

    error_fn("get_pdb_dict: out of memory allocation %d pdb dictionaries",N_PDB_DICTS+1);
  }

  PDB_DICTS[N_PDB_DICTS] = new_dict;

  N_PDB_DICTS++;

  // return new dictionary:

  return(new_dict);
}



MOLECULE* read_cif_molecule(char *filename) {

  char line[MAX_LINE_LEN];
  char data_type[MAX_LINE_LEN],data_name[MAX_LINE_LEN];
  char field_name[MAX_LINE_LEN],field_value[MAX_LINE_LEN];
  PLI_FILE *ciffile;
  ATOM_TYPING_SCHEME *scheme;
  MOLECULE *molecule;

  ciffile = open_file(filename,"r");

  if (ciffile == NULL) {

    return(NULL);
  }

  scheme = get_atom_typing_scheme();

  molecule = (MOLECULE*) malloc(sizeof(MOLECULE));

  if (molecule == NULL) {

    error_fn("read_cif_molecule: out of memory allocating molecule");
  }

  init_molecule(molecule,50);

  strcpy(molecule->filename,filename);

  strcpy(data_type,"unknown");
  strcpy(data_name,"unknown");
  strcpy(field_name,"unknown");
  strcpy(field_value,"unknown");

  while (!end_of_file(ciffile)) {

    if (read_line(line,MAX_LINE_LEN,ciffile) == NULL)
      break;

    if (!strncmp(line,"data_",5)) {

      strcpy(data_type,"data_");

    } else if (!strncmp(line,"loop_",4)) {

      strcpy(data_type,"loop_");

    } else if (!strncmp(line,"_",1)) {

      if (!strcmp(data_type,"data_")) {

	sscanf(line,"%[^.].%s %s",data_name,field_name,field_value);

      } else if (!strcmp(data_type,"loop_")) {

	sscanf(line,"%[^.].%*s",data_name);
      }

    } else if (strncmp(line,"#",1)) {

      if (!strcmp(data_type,"loop_")) {

	if (!strcmp(data_name,"_chem_comp_atom")) {

	  read_cif_atom(line,molecule,scheme);

	} else if (!strcmp(data_name,"_chem_comp_bond")) {

	  read_cif_bond(line,molecule);
	}
      }
    }
  }

  close_file(ciffile);

  return(molecule);
}



static void read_cif_atom(char *line,MOLECULE *molecule,ATOM_TYPING_SCHEME *scheme) {

  char name[MAX_LINE_LEN],alt_name[MAX_LINE_LEN],elname[MAX_LINE_LEN];
  char format1[20],format2[20],format[MAX_LINE_LEN];
  ATOM *atom;

  // first work out whether atom names have double quotes around them:

  if (sscanf(line,"%*s \"%s",name) == 1) {

    strcpy(format1,"\"%[^\"]\"");

  } else {

    strcpy(format1,"%s");
  }

  sprintf(format,"%%*s %s \"%%s",format1);

  if (sscanf(line,format,name,alt_name) == 2) {

    strcpy(format2,"\"%[^\"]\"");

  } else {

    strcpy(format2,"%s");
  }

  // read the atom line:

  sprintf(format,"%%*s %s %s %%s",format1,format2);

  if (sscanf(line,format,name,alt_name,elname) != 3) {

    warning_fn("read_cif_atom: corrupt atom line:\n%s",line);

    return;
  }

  if ((strlen(name)>10) || (strlen(alt_name) > 10) || (strlen(elname) > 2)) {

    warning_fn("read_cif_atom: atom or element name too long:\n%s",line);
 
    return;
  }

  // add the atom to the molecule:

  realloc_atoms(molecule);

  atom = molecule->atom + molecule->natoms;

  init_atom(atom);

  atom->id = molecule->natoms + 1;

  strcpy(atom->name,name);
  strcpy(atom->alt_name,alt_name);

  atom->element = get_element(elname,scheme);

  atom->molecule = molecule;

  molecule->natoms++;

  return;
}



static void read_cif_bond(char *line,MOLECULE *molecule) {

  char name1[MAX_LINE_LEN],name2[MAX_LINE_LEN],type[MAX_LINE_LEN];
  char format1[20],format2[20],format[MAX_LINE_LEN];
  ATOM *atom1,*atom2;
  BOND *bond;

  // first work out whether atom names have double quotes around them:

  if (sscanf(line,"%*s \"%s",name1) == 1) {

    strcpy(format1,"\"%[^\"]\"");

  } else {

    strcpy(format1,"%s");
  }

  sprintf(format,"%%*s %s \"%%s",format1);

  if (sscanf(line,format,name1,name2) == 2) {

    strcpy(format2,"\"%[^\"]\"");

  } else {

    strcpy(format2,"%s");
  }

  // read the bond line:

  sprintf(format,"%%*s %s %s %%s",format1,format2);

  if (sscanf(line,format,name1,name2,type) != 3) {

    warning_fn("read_cif_bond: corrupt bond line:\n%s",line);

    return;
  }

  atom1 = get_cif_atom(name1,molecule);
  atom2 = get_cif_atom(name2,molecule);

  if ((atom1 == NULL) || (atom2 == NULL)) {

    warning_fn("read_cif_bond: cannot find bond atom(s):\n",line);

    return;
  }

  if (((atom1->element != NULL) && (atom1->element->flags & METAL_ELEMENT)) ||
      ((atom2->element != NULL) && (atom2->element->flags & METAL_ELEMENT))) {

    return;
  }

  // add the bond to the molecule:

  realloc_bonds(molecule);

  bond = get_bond(molecule,atom1,atom2);

  if (bond == NULL) {

    bond = molecule->bond + molecule->nbonds;
      
    bond->atom1 = atom1;
    bond->atom2 = atom2;

    if (!strcmp(type,"SING")) {

      bond->type = 1; 	

    } else if (!strcmp(type,"DOUB")) {

      bond->type = 2;

    } else if (!strcmp(type,"TRIP")) {

      bond->type = 3;

    } else {

      warning_fn("read_cif_bond: unknown bond type '%s' - set to single",type);

      bond->type = 1;
    }

    molecule->nbonds++;
  }

  return;
}



ATOM* get_cif_atom(char *name,MOLECULE *molecule) {

  int i;
  char name1[10];
  ATOM *atom;

  if (strlen(name) > 10) {

    warning_fn("get_cif_atom: atom name '%s' too long",name);

    return(NULL);
  }

  remove_outer_spaces(name,name1);

  for (i=0,atom=molecule->atom;i<molecule->natoms;i++,atom++) {

    if ((!strcmp(atom->name,name1)) || (!strcmp(atom->alt_name,name1))) {

      return(atom);
    }
  }

  return(NULL);
}



ELEMENT* pdb_atom_element(ATOM *atom,ATOM_TYPING_SCHEME *scheme) {

  int i,len,start;
  char c;
  ELEMENT *element;
  ATOM_TYPE *type;

  type = resatoms2type(atom,scheme);

  if (type != NULL) {

    if ((type->united_atom != NULL) &&
	(type->united_atom->atom_node != NULL) &&
	(type->united_atom->atom_node->element != NULL)) {

      return(type->united_atom->atom_node->element);

    } else if (type->element != NULL) {
      
      return(type->element);
    }
  }

  c = atom->name[0];

  start = ((c == ' ') || (c < 'A') || (c > 'Z')) ? 1 : 0;

  len = (start == 0) ? 2 : 1;

  for (i=0,element=scheme->elements;i<scheme->n_elements;i++,element++) {

    if (strlen(element->name) == len) {

      if (!strncmp(atom->name+start,element->name,len)) {

	return(element);
      }
    }
  }

  return(NULL);
}



int write_pdb_molecule(MOLECULE *molecule,char *filename) {

  int i,j,id,maxid,maxsubid;
  PLI_FILE *pdbfile;
  ELEMENT *elements;
  ATOM *atom,vatom;

  pdbfile = open_file(filename,"w");

  if (pdbfile == NULL) {

    return(1);
  }

  maxid = 0;
  maxsubid = 0;

  for (i=0,atom=molecule->atom;i<molecule->natoms;i++,atom++) {

    if (!(atom->flags & SKIP_ATOM)) {

      write_pdb_atom(pdbfile,atom);

      if (atom->id > maxid)
	maxid = atom->id;

      if ((atom->chain == 'X') && (atom->subid > maxsubid))
	maxsubid = atom->subid;
    }
  }

  init_atom(&vatom);

  vatom.element = get_element("H",NULL);
  strcpy(vatom.subname,"UVW");
  vatom.chain = 'X';
  vatom.subid = maxsubid + 1;

  id = maxid + 1;

  for (i=0,atom=molecule->atom;i<molecule->natoms;i++,atom++) {

    if (!(atom->flags & SKIP_ATOM)) {

      if (atom->geometry != NULL) {

	if (atom->geometry->u_axis) {

	  vatom.id = id;

	  strcpy(vatom.name," XU ");

	  for (j=0;j<3;j++) {

	    vatom.position[j] = atom->position[j] + atom->u[j];
	  }
	  
	  //write_pdb_atom(pdbfile,&vatom);

	  //id++;
	}

	if (atom->geometry->v_axis) {

	  vatom.id = id;

	  strcpy(vatom.name," XV ");

	  for (j=0;j<3;j++) {

	    vatom.position[j] = atom->position[j] + atom->v[j];
	  }

	  //write_pdb_atom(pdbfile,&vatom);

	  //id++;
	}
      }
    }
  }

  close_file(pdbfile);
}



void write_pdb_atom_list(PLI_FILE *pli_file,ATOMLIST *atomlist,unsigned long int oflags) {

  int i,j,k,nbonds;
  ATOM **atomp,*atom,*batom,vatom;
  BOND **bond;
  BONDLIST *bondlist;

  bondlist = atomlist2bondlist(atomlist);

  init_atom(&vatom);

  for (i=0,atomp=atomlist->atom;i<atomlist->natoms;i++,atomp++) {

    atom = *atomp;

    if (!(atom->flags & SKIP_ATOM)) {

      write_pdb_atom(pli_file,atom);

      if (atom->id > vatom.id)
	vatom.id = atom->id;

      if ((atom->chain == 'X') && (atom->subid > vatom.subid))
	vatom.subid = atom->subid;
    }
  }

  if ((oflags & OUTPUT_AXES) || (oflags & OUTPUT_VPTS)) {

    vatom.element = get_element("H",NULL);
    strcpy(vatom.subname,"UVW");
    vatom.chain = 'X';
    vatom.id++;
    vatom.subid++;
  }

  if (oflags & OUTPUT_AXES) {

    for (i=0,atomp=atomlist->atom;i<atomlist->natoms;i++,atomp++) {

      if (!(atom->flags & SKIP_ATOM)) {

	write_pdb_atom_axes(pli_file,*atomp,&vatom);
      }
    }
  }

  if (oflags & OUTPUT_VPTS) {

    for (i=0,atomp=atomlist->atom;i<atomlist->natoms;i++,atomp++) {

      if (!(atom->flags & SKIP_ATOM)) {

	write_pdb_atom_vpts(pli_file,*atomp,&vatom);
      }
    }
  }

  for (i=0,atomp=atomlist->atom;i<atomlist->natoms;i++,atomp++) {

    atom = *atomp;

    if (!(atom->flags & SKIP_ATOM)) {

      nbonds = 0;

      for (j=0,bond=bondlist->bond;j<bondlist->nbonds;j++,bond++) {

	batom = ((*bond)->atom1 == atom) ? (*bond)->atom2 : ((*bond)->atom2 == atom) ? (*bond)->atom1 : NULL;

	if ((batom != NULL) && (!(batom->flags & SKIP_ATOM))) {

	  if (nbonds == 0) {

	    write_line(pli_file,"CONECT%5d",atom->id);
	  }

	  for (k=0;k<(*bond)->type;k++) {

	    write_line(pli_file,"%5d",batom->id);
	  }
	  
	  nbonds++;
	}
      }

      if (nbonds > 0) {

	write_line(pli_file,"\n");
      }
    }
  }

  write_line(pli_file,"END\n");

  free_bondlist(bondlist);
}



void write_pdb_atom(PLI_FILE *pli_file,ATOM *atom) {

  char subname[4],formal_charge[3];
  UNITED_ATOM *uatom;

  substring(atom->subname,0,3,subname);

  upper_case(subname);

  strcpy(formal_charge,"");

  if ((atom->type) && (atom->type->united_atom)) {

    uatom = atom->type->united_atom;

    if (uatom->formal_charge < 0) {

      sprintf(formal_charge,"%1d-",uatom->formal_charge);

    } else if (uatom->formal_charge > 0) {

      sprintf(formal_charge,"%1d+",uatom->formal_charge);
    }
  }

  write_line(pli_file,"%-6s%5d %4s%c%3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n","ATOM ",
	     atom->id,
	     atom->name,
	     atom->altloc,
	     subname,
	     atom->chain,
	     atom->subid,
	     atom->icode,
	     atom->position[0],atom->position[1],atom->position[2],
	     atom->occupancy,
	     atom->bfactor,
	     atom->element->name,
	     formal_charge);
}



static void write_pdb_atom_axes(PLI_FILE *file,ATOM *atom,ATOM *vatom) {

  if (atom->geometry != NULL) {

    if (atom->geometry->u_axis) {

      strcpy(vatom->name," XU ");

      sum_vector(atom->position,atom->u,vatom->position);

      write_pdb_atom(file,vatom);

      vatom->id++;
    }

    if (atom->geometry->v_axis) {

      strcpy(vatom->name," XV ");

      sum_vector(atom->position,atom->v,vatom->position);

      write_pdb_atom(file,vatom);

      vatom->id++;
    }
  }
}



static void write_pdb_atom_vpts(PLI_FILE *file,ATOM *atom,ATOM *vatom) {

  int i;
  HBOND_GEOMETRY *hbond_geometry;

  hbond_geometry = atom->hbond_geometry;

  if ((hbond_geometry) && (hbond_geometry->npts)) {

    for (i=0;i<hbond_geometry->npts;i++) {

      strcpy(vatom->name," VP ");

      sum_vector(atom->position,atom->vpts[i],vatom->position);

      write_pdb_atom(file,vatom);

      vatom->id++;
    }
  }
}



void atom2pdbstr(ATOM *atom,char *str) {

  sprintf(str,"%-6s%5d %4s%c%3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n","ATOM ",
	  atom->id,
	  atom->name,
	  atom->altloc,
	  atom->subname,
	  atom->chain,
	  atom->subid,
	  atom->icode,
	  atom->position[0],atom->position[1],atom->position[2],
	  atom->occupancy,
	  atom->bfactor,
	  atom->element->name);
}
