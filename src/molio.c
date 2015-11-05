#include "pli.h"



static void write_atom_basic(PLI_FILE*,ATOM*,enum OUTPUT_FORMAT,unsigned long int);
static void write_atom_pliff(PLI_FILE*,ATOM*,enum OUTPUT_FORMAT,unsigned long int);
static void write_atom_type(PLI_FILE*,ATOM*,enum OUTPUT_FORMAT,unsigned long int);
static void write_atom_stype(PLI_FILE*,ATOM*,enum OUTPUT_FORMAT,unsigned long int);




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



void write_atom_list(PLI_FILE *file,ATOMLIST *list,enum ATOM_STYLE astyle,enum OUTPUT_FORMAT oformat,unsigned long int oflags) {

  int i;
  ATOM **atomp;

  if (list == NULL) {

    error_fn("write_atoms_list: no atoms");
  }

  for (i=0,atomp=list->atom;i<list->natoms;i++,atomp++) {

    write_atom(file,*atomp,astyle,oformat,oflags);
  }
}



enum ATOM_STYLE get_atom_style(char *name) {

  if (!strcmp(name,"basic")) {

    return(BASIC_ASTYLE);

  } else if (!strcmp(name,"pdb")) {

    return(PDB_ASTYLE);

  } else if (!strcmp(name,"type")) {

    return(TYPE_ASTYLE);

  } else if (!strcmp(name,"pliff")) {

    return(PLIFF_ASTYLE);
  }

  error_fn("get_atom_style: no such atom format '%s'",name);
}



void write_atom(PLI_FILE *file,ATOM *atom,enum ATOM_STYLE astyle,enum OUTPUT_FORMAT oformat,unsigned long int oflags) {

  if (astyle == BASIC_ASTYLE) {

    write_atom_basic(file,atom,oformat,oflags);

  } else if (astyle == PDB_ASTYLE) {

    write_pdb_atom(file,atom);

  } else if (astyle == PLIFF_ASTYLE) {

    write_atom_pliff(file,atom,oformat,oflags);

  } else if (astyle == TYPE_ASTYLE) {

    write_atom_type(file,atom,oformat,oflags);

  } else {

    error_fn("write_atom: unknown atom format");
  }
}



static void write_atom_basic(PLI_FILE *file,ATOM *atom,enum OUTPUT_FORMAT oformat,unsigned long int oflags) {

  int sym_mol,json;
  char subname[4];
  MOLECULE *molecule;

  molecule = atom->molecule;

  sym_mol = (molecule != NULL) ? ((molecule->flags & SYMMETRY_MOLECULE) ? 1 : 0) : -1;

  subname[3] = '\0';

  (strlen(atom->subname) > 3) ? (strncpy(subname,atom->subname,3)) : (strcpy(subname,atom->subname));

  json = (oformat == JSON) ? 1 : 0;

  write_line(file,(json) ? 
	     "{\"atom\":{\"id\":%d,\"name\":\"%s\",\"altloc\":\"%c\",\"subname\":\"%s\",\"chain\":\"%c\",\"subid\":%d,\"icode\":\"%c\",\"set\":%d,\"sym_mol\":%d" : 
	     "ATOM %5d %4s %c %3s %c %4d %c %5d %5d",
	     atom->id,
	     atom->name,
	     ((!json) && (atom->altloc == ' ') ? '_' : atom->altloc),
	     subname,
	     ((!json) && (atom->chain == ' ') ? '_' : atom->chain),
	     atom->subid,
	     ((!json) && (atom->icode == ' ') ? '_' : atom->icode),
	     atom->set,
             sym_mol);

  if (oflags & OUTPUT_STYPES) {

    if (json) {

      write_line(file,",");
    }

    write_atom_stype(file,atom,oformat,oflags);
  }

  if (oflags & OUTPUT_SCORES) {

    if (json) {

      write_line(file,",");
    }

    write_atom_scores(file,atom,oformat,oflags);
  }

  if (oflags & OUTPUT_SIMB) {

    write_line(file,(json) ? ",\"tdf\"=%.4lf,\"td\"=%.4lf" : "%10.4lf%10.4lf",atom->tdf,atom->td);
  }

  write_line(file,(json) ? "}}\n" : "\n");
}



static void write_atom_pliff(PLI_FILE *file,ATOM *atom,enum OUTPUT_FORMAT oformat,unsigned long int oflags) {

  int eid,tid,uid,sym_mol,json;
  char subname[4];
  MOLECULE *molecule;

  molecule = atom->molecule;

  eid = (atom->element != NULL) ? atom->element->id : -1;
  tid = (atom->type != NULL) ? atom->type->id : -1;
  uid = (atom->type != NULL) ? ((atom->type->united_atom != NULL) ? atom->type->united_atom->id : -1) : -1;

  sym_mol = (molecule != NULL) ? ((molecule->flags & SYMMETRY_MOLECULE) ? 1 : 0) : -1;

  subname[3] = '\0';

 (strlen(atom->subname) > 3) ? (strncpy(subname,atom->subname,3)) : (strcpy(subname,atom->subname));

  json = (oformat == JSON) ? 1 : 0;

  write_line(file,(json) ?
	     "{\"atom\":{\"id\":%d,\"name\":\"%s\",\"altloc\":\"%c\",\"subname\":\"%s\",\"chain\":\"%c\",\"subid\":%d,\"icode\":\"%c\",\"x\":%lf,\"y\":%lf,\"z\":%lf,\"vdw_radius\":%lf,\"occupancy\":%lf,\"bfactor\":%lf,\"intra_area\":%lf,\"contact_area\":%lf,\"exposed_area\":%lf,\"set\":%d,\"eid\":%d,\"uid\":%d,\"tid\":%d,\"error_flags\":%d,\"sym_mol\":%d" :
	     "ATOM %5d %4s %c %3s %c %4d %c   %8.3f%8.3f%8.3f   %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f    %3d %3d %3d %3d %3d %3d",
	     atom->id,
	     atom->name,
	     ((atom->altloc == ' ') ? '_' : atom->altloc),
	     subname,
	     ((atom->chain == ' ') ? '_' : atom->chain),
	     atom->subid,
	     ((atom->icode == ' ') ? '_' : atom->icode),
	     atom->position[0],atom->position[1],atom->position[2],
	     atom->vdw_radius,
	     atom->occupancy,
	     atom->bfactor,
	     (atom->covalent_area)+(atom->intra_area),
	     atom->contact_area,
	     atom->exposed_area,
	     atom->set,
	     eid,uid,tid,
	     atom->error_flags,
	     sym_mol);

  if (oflags & OUTPUT_STYPES) {

    if (json) {

      write_line(file,",");
    }

    write_atom_stype(file,atom,oformat,oflags);
  }

  if (oflags & OUTPUT_SCORES) {

    if (json) {

      write_line(file,",");
    }

    write_atom_scores(file,atom,oformat,oflags);
  }

  if (oflags & OUTPUT_SIMB) {

    write_line(file,(json) ? ",\"tdf\"=%.4lf,\"td\"=%.4lf" : "%10.4lf%10.4lf",atom->tdf,atom->td);
  }

  write_line(file,(json) ? "}}\n" : "\n");
}



static void write_atom_type(PLI_FILE *file,ATOM *atom,enum OUTPUT_FORMAT oformat,unsigned long int oflags) {

  int tid,uid,eid,json,sym_mol;
  char type_name[MAX_LINE_LEN],uatom_name[MAX_LINE_LEN],element_name[MAX_LINE_LEN],subname[4];
  ATOM_TYPE *type;
  UNITED_ATOM *uatom;
  ELEMENT *element;
  MOLECULE *molecule;

  molecule = atom->molecule;

  sym_mol = (molecule != NULL) ? ((molecule->flags & SYMMETRY_MOLECULE) ? 1 : 0) : -1;

  subname[3] = '\0';

  (strlen(atom->subname) > 3) ? (strncpy(subname,atom->subname,3)) : (strcpy(subname,atom->subname));

  type = atom->type;
  uatom = type->united_atom;
  element = atom->element;

  tid = (type == NULL) ? -1 : type->id;
  uid = (uatom == NULL) ? -1 : uatom->id;
  eid = (element == NULL) ? -1 : element->id;

  sprintf(type_name,"\"%s\"",((tid == -1) ? "Unknown" : type->name));
  sprintf(uatom_name,"\"%s\"",((uid == -1) ? "Unknown" : uatom->name));
  sprintf(element_name,"\"%s\"",((eid == -1) ? "Unknown" : element->name));

  json = (oformat == JSON) ? 1 : 0;

  write_line(file,(json) ? 
	     "{\"atom\":{\"id\":%d,\"name\":\"%s\",\"altloc\":\"%c\",\"subname\":\"%s\",\"chain\":\"%c\",\"subid\":%d,\"icode\":\"%c\",\"set\":%d,\"tid\":%d,\"tname\":%s,\"uid\":%d,\"uname\":%s,\"eid\":%d,\"ename\":%s,\"sym_mol\":%d" : 
	     "ATOM %5d %4s %c %3s %c %4d %c %5d %5d %-30s %5d %-20s %5d %-5s %5d",
	     atom->id,
	     atom->name,
	     ((!json) && (atom->altloc == ' ') ? '_' : atom->altloc),
	     subname,
	     ((!json) && (atom->chain == ' ') ? '_' : atom->chain),
	     atom->subid,
	     ((!json) && (atom->icode == ' ') ? '_' : atom->icode),
	     atom->set,
	     tid,type_name,uid,uatom_name,eid,element_name,sym_mol);

  if (oflags & OUTPUT_STYPES) {

    if (json) {

      write_line(file,",");
    }

    write_atom_stype(file,atom,oformat,oflags);
  }

  if (oflags & OUTPUT_SCORES) {

    if (json) {

      write_line(file,",");
    }

    write_atom_scores(file,atom,oformat,oflags);
  }

  if (oflags & OUTPUT_SIMB) {

    write_line(file,(json) ? ",\"tdf\"=%.4lf,\"td\"=%.4lf" : "%10.4lf%10.4lf",atom->tdf,atom->td);
  }

  write_line(file,(json) ? "}}\n" : "\n");
}



static void write_atom_stype(PLI_FILE *file,ATOM *atom,enum OUTPUT_FORMAT oformat,unsigned long int oflags) {

  int don,acc,met,lip,aro;
  ATOM_TYPE *type;
  UNITED_ATOM * uatom;
  ELEMENT *element;

  type = atom->type;

  don = acc = met = lip = aro = 0;

  if (type->id != -1) {

    element = atom->element;
    uatom = type->united_atom;

    don = (type->flags & HBOND_DONOR_ATOM_TYPE) ? 1 : 0;
    acc = (type->flags & HBOND_ACCEPTOR_ATOM_TYPE) ? 1 : 0;
    met = (type->flags & METAL_ATOM_TYPE) ? 1 : 0;

    aro = (type->flags & AROMATIC_ATOM_TYPE) ? 1 : 0;
    lip = (type->flags & LIPOPHILIC_ATOM_TYPE) ? 1 : 0;
  }

  write_line(file,(oformat == JSON) ? 
	     "\"don\":%d,\"acc\":%d,\"met\":%d,\"lip\":%d,\"aro\":%d" : 
	     " %10d %10d %10d %10d %10d",don,acc,met,lip,aro);
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
