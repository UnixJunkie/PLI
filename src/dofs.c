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



#define FD_DX  0.000010
#define MIN_DX 0.01



static void add_dof_to_list(DOF_LIST*,enum DOF_TYPE,int,TORSION*,ATOMLIST*);
static void add_ligand_rigid_body_to_dof_list(SYSTEM*,enum DOF_TYPE,DOF_LIST*);
static void add_ligand_torsions_to_dof_list(SYSTEM*,DOF_LIST*);
static void add_ligand_atoms_to_dof_list(SYSTEM*,DOF_LIST*,unsigned int);
static void add_water_rigid_body_to_dof_list(SYSTEM*,DOF_LIST*);
static DOF_LIST* alloc_dof_list(void);
static void init_dof_list(DOF_LIST*);
static int setup_dof(DOF*,enum DOF_TYPE,int,TORSION*,ATOMLIST*);
static double dof_shift_dx(ATOMLIST*,double**);



unsigned int dofs_get_space(char *flags) {

  unsigned int space;
  char *flag,cflags[MAX_LINE_LEN];

  space = 00;

  strcpy(cflags,flags);

  flag = strtok(cflags,",");

  while (flag) {

    if (!strcmp(flag,"ligand_translation")) {

      space |= DOF_LIGAND_TRANSLATION;
      
    } else if (!strcmp(flag,"ligand_rotation")) {
      
      space |= DOF_LIGAND_ROTATION;
      
    } else if (!strcmp(flag,"ligand_torsions")) {
      
      space |= DOF_LIGAND_TORSIONS;
      
    } else if (!strcmp(flag,"ligand_rigid_body")) {
      
      space |= DOF_LIGAND_RIGID_BODY;

    } else if (!strcmp(flag,"ligand_atoms")) {
     
      space |= DOF_LIGAND_ATOMS;
     
    } else if (!strcmp(flag,"ligand")) {
      
      space |= DOF_LIGAND;
      
    } else if (!strcmp(flag,"water_translation")) {
      
      space |= DOF_WATER_TRANSLATION;
      
    } else if (!strcmp(flag,"water")) {
      
      space |= DOF_WATER;
      
    } else {

      error_fn("update_space: unknown output option '%s'",flag);
    }

    flag = strtok(NULL,",");
  }

  return(space);
}



unsigned int dofs_molecule_space(unsigned int dof,unsigned int flags) {

  unsigned int space;

  space = 0;

  if (flags & LIGAND_MOLECULE) {

    if (dof & DOF_LIGAND_TRANSLATION) {

      space |= DOF_LIGAND_TRANSLATION;
    }

    if (dof & DOF_LIGAND_ROTATION) {
      
      space |= DOF_LIGAND_ROTATION;
    }

    if (dof & DOF_LIGAND_TORSIONS) {
      
      space |= DOF_LIGAND_TORSIONS;
    }

    if (dof & DOF_LIGAND_RIGID_BODY) {

      space |= DOF_LIGAND_RIGID_BODY;
    }

    if (space & DOF_LIGAND_ATOMS) {
     
      space |= DOF_LIGAND_ATOMS;
    }

    if (dof & DOF_LIGAND) {

      space |= DOF_LIGAND;
    }

  } else if (flags & PROTEIN_MOLECULE) {
      
    if (space & DOF_WATER_TRANSLATION) {

      space |= DOF_WATER_TRANSLATION;
    }

    if (space & DOF_WATER) {

      space |= DOF_WATER;
    }
  }

  return(space);
}



void set_dof_list_shifts(DOF_LIST *list,double shift) {

  int i;
  double sum,scale;
  DOF *variable;

  // scale dof shifts:

  sum = 0.0;

  for (i=0,variable=list->variables;i<list->n_variables;i++,variable++) {

    sum += (list->cg_initialised) ? sqr((variable->cg)) : sqr((variable->fd));
  }

  scale = shift/sqrt(sum);

  for (i=0,variable=list->variables;i<list->n_variables;i++,variable++) {

    variable->shift = (list->cg_initialised) ? -scale*((variable->cg)) : -scale*((variable->fd));
  }
}



void apply_dof_list_shifts(DOF_LIST *list) {
  
  int i;
  DOF *variable;

  for (i=0,variable=list->variables;i<list->n_variables;i++,variable++) {

    apply_dof_shift(variable,variable->shift);
  }
}



void apply_dof_shift(DOF *variable,double shift) {
  
  int i;
  double v[4],c[4],ci[4],tm[4][4],tmi[4][4],rm[4][4],m1[4][4],m2[4][4];
  ATOMLIST *atomlist;
  ATOM **atomp,*atom;

  atomlist = variable->atomlist;

  null_vector(v);

  if (variable->type == RB_TRANSLATION) {

    v[variable->axis_id] = shift;

    for (i=0,atomp=atomlist->atom;i<atomlist->natoms;i++,atomp++) {

      atom = *atomp;

      shift_atom(atom,v);

      if (atomlist->natoms == 1) {

	set_atom_axes(atom);
      }
    }

  } else if (variable->type == RB_ROTATION) {

    v[variable->axis_id] = shift;

    atomlist2centroid(atomlist,c);

    invert_vector(c,ci);

    translation_matrix(c,tm);
    translation_matrix(ci,tmi);

    euler_matrix(v[0],v[1],v[2],rm);

    calc_matrix_product(tmi,rm,m1);
    calc_matrix_product(m1,tm,m2);

    for (i=0,atomp=atomlist->atom;i<atomlist->natoms;i++,atomp++) {
      
      transform_atom(*atomp,m2);
    }

  } else if (variable->type == BOND_ROTATION) {

    rotate_torsion(variable->torsion,shift);

    //rotation_matrix(variable->torsion->atom1->position,variable->torsion->atom2->position,shift,rm);

    //for (i=0,atomp=atomlist->atom;i<atomlist->natoms;i++,atomp++) {
      
    //  transform_atom(*atomp,rm);
    //}
  }
}



DOF_LIST* setup_dof_list(SYSTEM *system) {

  int i,j;
  unsigned int space;
  DOF *variable;
  DOF_LIST *list = NULL;
  ATOM **atomp,*atom;
  ATOMLIST *atomlist;

  space = system->dof;

  list = (DOF_LIST*) alloc_dof_list();

  init_dof_list(list);

  if (space & DOF_LIGAND_TRANSLATION) {

    add_ligand_rigid_body_to_dof_list(system,RB_TRANSLATION,list);
  }

  if (space & DOF_LIGAND_ROTATION) {

    add_ligand_rigid_body_to_dof_list(system,RB_ROTATION,list);
  }

  if (space & DOF_LIGAND_TORSIONS) {

    add_ligand_torsions_to_dof_list(system,list);
  }

  if (space & DOF_LIGAND_ATOMS) {

    add_ligand_atoms_to_dof_list(system,list,0);
  }

  if (space & DOF_LIGAND_MATCHED_ATOMS) {

    add_ligand_atoms_to_dof_list(system,list,MATCHED_ATOM);
  }

  if (space & DOF_WATER_TRANSLATION) {

    add_water_rigid_body_to_dof_list(system,list);
  }

  // setup list of moving atoms:

  list->atomlist = alloc_atomlist();

  init_atomlist(list->atomlist);
 
  for (i=0,variable=list->variables;i<list->n_variables;i++,variable++) {

    atomlist = variable->atomlist;

    if (atomlist == NULL) {

      error_fn("no atoms\n");
    }

    for (j=0,atomp=atomlist->atom;j<atomlist->natoms;j++,atomp++) {

      atom = *atomp;

      if (!atom_in_list(list->atomlist,atom)) {

	add_atom_to_list(list->atomlist,atom,0);
      }
    }
  }

  if (list->atomlist->natoms) {

    list->pos = (double**) alloc_2d_fmatrix(list->atomlist->natoms,4);

    if (list->pos == NULL) {

      error_fn("setup_dof_list: out of memory allocating positions");
    }

    list->u = (double**) alloc_2d_fmatrix(list->atomlist->natoms,4);

    if (list->u == NULL) {

      error_fn("setup_dof_list: out of memory allocating axes");
    }

    list->v = (double**) alloc_2d_fmatrix(list->atomlist->natoms,4);

    if (list->v == NULL) {

      error_fn("setup_dof_list: out of memory allocating axes");
    }

    list->w = (double**) alloc_2d_fmatrix(list->atomlist->natoms,4);

    if (list->w == NULL) {

      error_fn("setup_dof_list: out of memory allocating axes");
    }
  }

  return(list);
}



static void add_ligand_rigid_body_to_dof_list(SYSTEM *system,enum DOF_TYPE type,DOF_LIST *list) {

  int i,j;
  ATOMLIST *selection;
  MOLECULE **moleculep,*molecule;

  for (i=0,moleculep=system->molecule_list->molecules;i<system->molecule_list->n_molecules;i++,moleculep++) {

    molecule = *moleculep;

    selection = molecule->selection;

    if ((selection) && (selection->natoms)) {

      if (molecule->flags & LIGAND_MOLECULE) {

	for (j=0;j<3;j++) {

	  add_dof_to_list(list,type,j,NULL,molecule->selection);
	}
      }
    }
  }
}



static void add_ligand_torsions_to_dof_list(SYSTEM *system,DOF_LIST *list) {

  int i,j;
  MOLECULE **moleculep,*molecule;
  TORSION *torsion;

  for (i=0,moleculep=system->molecule_list->molecules;i<system->molecule_list->n_molecules;i++,moleculep++) {

    molecule = *moleculep;

    if (molecule->flags & LIGAND_MOLECULE) {

      if (molecule->torsions) {

	for (j=0,torsion=molecule->torsions;j<molecule->n_torsions;j++,torsion++) {
	  
	  add_dof_to_list(list,BOND_ROTATION,-1,torsion,torsion->rotlist);
	}
      }
    }
  }
}



static void add_ligand_atoms_to_dof_list(SYSTEM *system,DOF_LIST *list,unsigned int flags) {

  int i,j,k;
  MOLECULE **moleculep,*molecule;
  ATOM *atom;
  ATOMLIST *alist;

  for (i=0,moleculep=system->molecule_list->molecules;i<system->molecule_list->n_molecules;i++,moleculep++) {

    molecule = *moleculep;

    if (molecule->flags & LIGAND_MOLECULE) {

      if (molecule->atom) {

	for (j=0,atom=molecule->atom;j<molecule->natoms;j++,atom++) {

	  if (atom->element->id != HYDROGEN) {

	    if ((!flags) || (atom->flags & flags)) {

	      alist = atom2atomlist(atom,0);
	    
	      for (k=0;k<3;k++) {
	      
		add_dof_to_list(list,RB_TRANSLATION,k,NULL,alist);
	      }
	    }
	  }
	}
      }
    }
  }
}



static void add_water_rigid_body_to_dof_list(SYSTEM *system,DOF_LIST *list) {

  int i,j;
  ATOMLIST *selection,*alist;
  ATOM **atomp,*atom;

  selection = system->selection;

  if ((selection) && (selection->natoms)) {

    for (i=0,atomp=selection->atom;i<selection->natoms;i++,atomp++) {

      atom = *atomp;

      if ((atom->flags & WATER_OXYGEN) && (!(atom->flags & SKIP_ATOM))) {

	alist = atom2atomlist(atom,0);

	for (j=0;j<3;j++) {

	  add_dof_to_list(list,RB_TRANSLATION,j,NULL,alist);
	}

	free_atomlist(alist);
      }
    }
  }
}



static void add_dof_to_list(DOF_LIST *list,enum DOF_TYPE type,int axis_id,TORSION *torsion,ATOMLIST *atomlist) {

  DOF *variable;

  if (list->n_alloc_variables == 0) {

    list->n_alloc_variables = 6;

    list->variables = (DOF*) calloc(sizeof(DOF),6);

    if (list->variables == NULL) {

      error_fn("add_dof_to_list: out of memory allocating variables");
    }

  } else if (list->n_alloc_variables == list->n_variables) {

    list->n_alloc_variables += 6;

    list->variables = (DOF*) realloc(list->variables,(list->n_alloc_variables)*sizeof(DOF));

    if (list->variables == NULL) {

      error_fn("add_dof_to_list: out of memory allocating variables");
    }
  }

  variable = list->variables + (list->n_variables);

  if (setup_dof(variable,type,axis_id,torsion,atomlist)) {

    list->n_variables += 1;

  } else {

    warning_fn("add_dof_to_list: skipping variable - insufficient effect on coordinates");

    free_dof(variable);
  }
}



static int setup_dof(DOF *dof,enum DOF_TYPE type,int axis_id,TORSION *torsion,ATOMLIST *list) {

  double f,dx;

  dof->type = type;
  dof->axis_id = axis_id;
  dof->shift = 0.0;

  dof->fd = 0.0;
  dof->fd_prev = 0.0;
  dof->cg = 0.0;

  dof->torsion = torsion;

  dof->atomlist = NULL;

  dof->pos = NULL;

  dof->u = NULL;
  dof->v = NULL;
  dof->w = NULL;

  if (list->natoms) {

    // copy the list of moving atoms so it can be freed later:

    dof->atomlist = alloc_atomlist();
    
    init_atomlist(dof->atomlist);

    copy_atomlist(list,dof->atomlist);

    // allocate space for positions of moving atoms:

    dof->pos = (double**) alloc_2d_fmatrix(list->natoms,4);

    if (dof->pos == NULL) {

      error_fn("setup_dof: out of memory allocating positions");
    }

    dof->u = (double**) alloc_2d_fmatrix(list->natoms,4);

    if (dof->u == NULL) {

      error_fn("setup_dof: out of memory allocating axes");
    }

    dof->v = (double**) alloc_2d_fmatrix(list->natoms,4);

    if (dof->v == NULL) {

      error_fn("setup_dof: out of memory allocating axes");
    }

    dof->w = (double**) alloc_2d_fmatrix(list->natoms,4);

    if (dof->w == NULL) {

      error_fn("setup_dof: out of memory allocating axes");
    }
  }

  if (type == RB_TRANSLATION) {

    dof->fd_shift = FD_DX;

  } else {

    //dof->fd_shift = 0.0005;

    store_positions(dof->atomlist,dof->pos,dof->u,dof->v,dof->w);

    apply_dof_shift(dof,1.0);

    dx = dof_shift_dx(dof->atomlist,dof->pos);

    restore_positions(dof->atomlist,dof->pos,dof->u,dof->v,dof->w);

    if (dx > MIN_DX) {

      dof->fd_shift = FD_DX/dx;

    } else {

      return(0);
    }
  }

  return(1);
}



void free_dof(DOF *dof) {

  if (dof->atomlist) {

    free_atomlist(dof->atomlist);
  }

  if (dof->pos) {

    free_2d_fmatrix(dof->pos);
    free_2d_fmatrix(dof->u);
    free_2d_fmatrix(dof->v);
    free_2d_fmatrix(dof->w);
  }
}



static DOF_LIST* alloc_dof_list(void) {

  DOF_LIST *list;

  list = (DOF_LIST*) malloc(sizeof(DOF_LIST));

  if (list == NULL) {

    error_fn("alloc_dof_list: out of memory allocating dof list");
  }

  return(list);
}



void free_dof_list(DOF_LIST *list) {

  int i;
  DOF *variable;

  if (list) {

    if (list->n_alloc_variables) {
      
      for (i=0,variable=list->variables;i<list->n_variables;i++,variable++) {

	free_dof(variable);
      }

      free(list->variables);
    }

    if (list->atomlist) {

      free_atomlist(list->atomlist);
    }

    if (list->pos) {

      free_2d_fmatrix(list->pos);
    }

    free(list);
  }
}



static void init_dof_list(DOF_LIST *list) {

  list->n_alloc_variables = 0;
  list->n_variables = 0;
  list->variables = NULL;
  list->cg_initialised = 0;
  list->atomlist = NULL;
  list->pos = NULL;
  list->u = NULL;
  list->v = NULL;
  list->w = NULL;
}



void store_positions(ATOMLIST *atomlist,double **pos,double **u,double **v,double **w) {

  int i;
  ATOM **atomp,*atom;

  for (i=0,atomp=atomlist->atom;i<atomlist->natoms;i++,atomp++) {

    atom = *atomp;

    copy_vector(atom->position,pos[i]);

    copy_vector(atom->u,u[i]);
    copy_vector(atom->v,v[i]);
    copy_vector(atom->w,w[i]);
  }
}



void restore_positions(ATOMLIST *atomlist,double **pos,double **u,double **v,double **w) {

  int i;
  ATOM **atomp,*atom;

  for (i=0,atomp=atomlist->atom;i<atomlist->natoms;i++,atomp++) {

    atom = *atomp;

    move_atom(atom,pos[i]);

    copy_vector(u[i],atom->u);
    copy_vector(v[i],atom->v);
    copy_vector(w[i],atom->w);
  }
}


ATOM_COORDS* get_list_coords(ATOMLIST *atomlist) {

  int i;
  ATOM **atomp,*atom;
  ATOM_COORDS *coords,*coord;

  if ((atomlist == NULL) || (atomlist->natoms == 0)) {

    return(NULL);
  }

  coords = (ATOM_COORDS*) calloc(atomlist->natoms,sizeof(ATOM_COORDS));

  if (coords == NULL) {

    error_fn("get_list_coords: out of memory allocating coords");
  }

  for (i=0,atomp=atomlist->atom,coord=coords;i<atomlist->natoms;i++,atomp++,coord++) {

    atom = *atomp;

    copy_vector(atom->position,coord->position);

    copy_vector(atom->u,coord->u);
    copy_vector(atom->v,coord->v);
    copy_vector(atom->w,coord->w);
  }

  return(coords);
}



void set_list_coords(ATOMLIST *atomlist,ATOM_COORDS *coords) {

  int i;
  ATOM **atomp,*atom;
  ATOM_COORDS *coord;

  if (atomlist == NULL) {

    return;
  }

  if (coords == NULL) {

    error_fn("set_list_coords: coords undefined");
  }

  for (i=0,atomp=atomlist->atom,coord=coords;i<atomlist->natoms;i++,atomp++,coord++) {

    atom = *atomp;

    move_atom(atom,coord->position);

    copy_vector(coord->u,atom->u);
    copy_vector(coord->v,atom->v);
    copy_vector(coord->w,atom->w);
  }
}



static double dof_shift_dx(ATOMLIST *atomlist,double **pos) {

  int i;
  double sumdx;
  ATOM **atomp,*atom;

  if ((atomlist) && (atomlist->natoms)) {

    sumdx = 0.0;

    for (i=0,atomp=atomlist->atom;i<atomlist->natoms;i++,atomp++) {

      atom = *atomp;

      sumdx += distance(atom->position,pos[i]);
    }

    return(sumdx/atomlist->natoms);
  }

  return(0.0);
}
