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



static int UNIQUE_ATOM_ID = 0;



static int linear_atom(ATOM*);
static void check_molecule_dict(MOLECULE*);
static void set_atom_connections(ATOM*,int);
static void molecule2connections(ATOM*,MOLECULE*,MOLECULE*,ATOM*,int,int);
static void bonds2connections(ATOM*,MOLECULE*);
static void map2connections(ATOM*,MAP*,MOLECULE*,ATOM*,int,int);
static void pair2connections(ATOM*,ATOM*,MOLECULE*,ATOM*,int,int);
static void update_atom_gridlist(ATOM*);



void atom2name(ATOM *atom,char *name) {

  char chain,icode;

  chain = (atom->chain == ' ') ? '_' : atom->chain;
  icode = (atom->icode == ' ') ? '_' : atom->icode;

  sprintf(name,"%6d %-6s %-6s %c %5d %c",atom->id,atom->name,atom->subname,chain,atom->subid,icode);
}



void init_atom(ATOM *atom) {

  int i;

  atom->unique_id = UNIQUE_ATOM_ID;

  UNIQUE_ATOM_ID++;

  atom->id = 0;
  atom->seqid = 0;
  atom->subid = 0;

  atom->element = NULL;

  strcpy(atom->name,"unknown");
  strcpy(atom->alt_name,"unknown");
  strcpy(atom->subname,"unknown");

  atom->chain = ' ';
  atom->icode = ' ';
  atom->altloc = ' ';

  null_vector(atom->position);
  null_vector(atom->original_position);

  null_vector(atom->u);
  null_vector(atom->v);
  null_vector(atom->w);

  for (i=0;i<MAX_VIRTUAL_POINTS;i++) {

    null_vector(atom->vpts[i]);
  }

  atom->tether_position = NULL;
  atom->tether_k = 1.0;

  atom->occupancy = 0.0;
  atom->bfactor = 0.0;
  atom->tdf = 0.0;
  atom->td = 0.0;
  atom->vdw_radius = DEFAULT_VDW_RADIUS;
  atom->vdw_radius_H2O = DEFAULT_VDW_RADIUS + 1.4;

  atom->node = NULL;
  atom->type = NULL;
  atom->ambiguous_type = NULL;
  atom->gridmap = NULL;
  atom->gridlist = NULL;
  atom->connections = NULL;
  atom->contactlist = NULL;
  atom->geometry = NULL;
  atom->hbond_geometry = NULL;
  atom->coordination = NULL;
  atom->ring = NULL;
  atom->molecule = NULL;
  atom->set = 0;
  atom->planar = -1;
  atom->n_hydrogens = 0;
  atom->score = 0.0;
  atom->ref_score = 0.0;
  atom->constraint_score = 0.0;
  atom->type_probability = 1.0;

  atom->exposed_area = 0.0;
  atom->contact_area = 0.0;
  atom->intra_area = 0.0;
  atom->covalent_area = 0.0;

  atom->da_area = 0.0;
  atom->don_area = 0.0;
  atom->ch_area = 0.0;
  atom->acc_area = 0.0;
  atom->met_area = 0.0;
  atom->lip_area = 0.0;
  atom->pol_area = 0.0;

  atom->group_status = 0;

  atom->cstatus = ATOM_NOTHING_CALCULATED;
  atom->status = ATOM_NEW;

  atom->flags = 0;
  atom->error_flags = 0;
}



void unprep_atom(ATOM *atom) {

  free_atomlist(atom->connections);
  free_contactlist(atom->contactlist);

  if (atom->coordination) {

    free(atom->coordination);
  }

  if (atom->ring) {

    free(atom->ring);
  }

  atom->node = NULL;
  atom->type = NULL;
  atom->geometry = NULL;
  atom->hbond_geometry = NULL;

  atom->connections = NULL;
  atom->contactlist = NULL;
  atom->coordination = NULL;
  atom->ring = NULL;
}



void init_bond(BOND *bond) {

  bond->type = -1;
  bond->atom1 = NULL;
  bond->atom2 = NULL;
  bond->torsion = NULL;
}



void init_molecule(MOLECULE *molecule,int n) {

  molecule->loaded = 0;
  molecule->lazy_load = 0;
  molecule->file_pos = -1;

  strcpy(molecule->name,"unknown");

  molecule->natoms = 0;
  molecule->nbonds = 0;

  molecule->n_alloc_atoms = n;
  molecule->n_alloc_bonds = n;

  if (n) {

    molecule->atom = (ATOM*) calloc(molecule->n_alloc_atoms,sizeof(ATOM));
    molecule->bond = (BOND*) calloc(molecule->n_alloc_bonds,sizeof(BOND));

  } else {

    molecule->atom = NULL;
    molecule->bond = NULL;
  }

  molecule->torsions = NULL;
  molecule->n_torsions = 0;

  molecule->nb_atom_pairs = NULL;

  molecule->use_pdb_dict = 0;
  molecule->use_hydrogens = 0;
  molecule->use_grids = 0;
  molecule->rely_on_bonds = 0;

  molecule->prepped = 0;

  molecule->flags = 0;

  molecule->covalent_map = NULL;
  molecule->contacts_map = NULL;

  molecule->selection = NULL;
  molecule->molsystem = NULL;

  molecule->n_active_waters = 0;
  molecule->n_alloc_active_waters = 0;
  molecule->active_waters = NULL;

  molecule->score = 0.0;

  molecule->system = NULL;
}



void prep_molecule(MOLECULE *molecule,SETTINGS *settings) {

  int i;
  double f;
  ATOM *atom;

  if (molecule->prepped) {

    return;
  }

  read_molecule_fp(molecule);

  for (i=0,atom=molecule->atom;i<molecule->natoms;i++,atom++) {

    atom->molecule = molecule;
  }

  // molecule atom selection:

  select_molecule_atoms(molecule,settings->selection);

  if (molecule->use_grids) {

    molecule->covalent_map = molecule2atommap(molecule,settings->max_covalent_dist,NULL,0);
    molecule->contacts_map = molecule2atommap(molecule,settings->max_contact_dist,NULL,1);
  }

  // create molecule system (for molecule in isolation):

  molecule->molsystem = molecule2system(molecule,settings);

  // set molecule connections:

  if ((molecule->use_pdb_dict) && (molecule->flags & LIGAND_MOLECULE)) {

    set_molecule_connections(molecule,1);

    f = molecule_dict_match(molecule);

    if (f < 0.8) {

      warning_fn("only %.1f%% of molecule atoms could be matched to dictionary; switching off dictionary for molecule '%s'",f,molecule->name);

      molecule->use_pdb_dict = 0;

      for (i=0,atom=molecule->atom;i<molecule->natoms;i++,atom++) {

	atom->error_flags |= ATOM_SKIP_DICTIONARY;
      }
    }

    reset_molecule_connections(molecule);
  }

  set_molecule_connections(molecule,0);

  // set nodes, atom types and radii:

  set_molecule_atom_nodes(molecule,settings->atom_typing_scheme);
  set_molecule_atom_types(molecule,settings->atom_typing_scheme);
  set_molecule_vdw_radii(molecule,settings->water_vdw_radius);
  set_molecule_atom_geometries(molecule,settings->atom_typing_scheme);
  
  // set special atom types (heme nitrogens for now):

  atom_type_heme_nitrogens(molecule->molsystem);

  // molecule torsions:

  if (molecule->flags & LIGAND_MOLECULE) {

    set_molecule_internal(molecule);
  }
  
  // simulated bfactors:

  //calc_system_simb(molecule->molsystem,1);

  reset_system(molecule->molsystem,0);

  molecule->prepped = 1;
}



void unprep_molecule(MOLECULE *molecule,SETTINGS *settings) {

  int i;
  ATOM *atom;
  TORSION *torsion;

  if (!(molecule->prepped)) {

    return;
  }

  free_map(molecule->covalent_map);
  free_map(molecule->contacts_map);

  for (i=0,atom=molecule->atom;i<molecule->natoms;i++,atom++) {

    unprep_atom(atom);
  }

  if (molecule->torsions) {

    for (i=0,torsion=molecule->torsions;i<molecule->n_torsions;i++,torsion++) {

      free_atomlist(torsion->alist1);
      free_atomlist(torsion->alist2);
      free_atomlist(torsion->rotlist);
    }

    free(molecule->torsions);
  }

  free_2d_imatrix(molecule->nb_atom_pairs);

  molecule->nb_atom_pairs = NULL;

  free_atomlist(molecule->selection);

  molecule->selection = NULL;

  free_molsystem(molecule->molsystem);

  molecule->molsystem = NULL;

  molecule->torsions = NULL;
  molecule->n_torsions = 0;
  molecule->nb_atom_pairs = NULL;
  molecule->selection = NULL;
  molecule->molsystem = NULL;

  if (molecule->lazy_load) {

    if (molecule->atom) {

      free(molecule->atom);
    }

    if (molecule->bond) {

      free(molecule->bond);
    }

    molecule->n_alloc_atoms = 0;
    molecule->natoms = 0;
    molecule->atom = NULL;

    molecule->n_alloc_bonds = 0;
    molecule->nbonds = 0;
    molecule->bond = NULL;

    molecule->loaded = 0;
  }

  molecule->prepped = 0;
}



void realloc_atoms(MOLECULE *molecule) {

  if (molecule->n_alloc_atoms == 0) {

    molecule->n_alloc_atoms = 50;

    molecule->atom = (ATOM*) calloc(molecule->n_alloc_atoms,sizeof(ATOM));

    if (molecule->atom == NULL) {

      error_fn("realloc_atoms: out of memory allocating %d atoms",molecule->n_alloc_atoms);
    }

  } else if (molecule->natoms == molecule->n_alloc_atoms) {

    molecule->n_alloc_atoms *= 2;

    molecule->atom = (ATOM*) realloc(molecule->atom,molecule->n_alloc_atoms*sizeof(ATOM));

    if (molecule->atom == NULL) {

      error_fn("realloc_atoms: out of memory reallocating %d atoms",molecule->n_alloc_atoms);
    }
  }
}



void realloc_bonds(MOLECULE *molecule) {

  if (molecule->n_alloc_bonds == 0) {

    molecule->n_alloc_bonds = 50;

    molecule->bond = (BOND*) calloc(molecule->n_alloc_bonds,sizeof(BOND));

    if (molecule->bond == NULL) {

      error_fn("realloc_bonds: out of memory allocating %d bonds",molecule->n_alloc_bonds);
    }

  } else if (molecule->nbonds == molecule->n_alloc_bonds) {

    molecule->n_alloc_bonds *= 2;

    molecule->bond = (BOND*) realloc(molecule->bond,molecule->n_alloc_bonds*sizeof(BOND));

    if (molecule->bond == NULL) {

      error_fn("realloc_bonds: out of memory reallocating %d bonds",molecule->n_alloc_bonds);
    }
  }
}



ATOM* get_atom(MOLECULE *molecule,int id) {

  int i;
  ATOM *atom;

  for (i=0,atom=molecule->atom;i<molecule->natoms;i++,atom++) {

    if (atom->id == id) {

      return(atom);
    }
  }

  return(NULL);
}



BOND* get_bond(MOLECULE *molecule,ATOM *atom1,ATOM *atom2) {

  int i;
  BOND *bond;

  for (i=0,bond=molecule->bond;i<molecule->nbonds;i++,bond++) {

    if (((bond->atom1 == atom1) && (bond->atom2 == atom2)) ||
	((bond->atom2 == atom1) && (bond->atom1 == atom2))) {

      return(bond);
    }
  }

  return(NULL);
}



BOND* add_bond(MOLECULE *molecule,ATOM *atom1,ATOM *atom2,int type) {

  ELEMENT *element1,*element2;
  BOND *bond;

  realloc_bonds(molecule);

  bond = molecule->bond + molecule->nbonds;
      
  bond->atom1 = atom1;
  bond->atom2 = atom2;
  bond->type = type; 	

  molecule->nbonds++;

  element1 = atom1->element;
  element2 = atom2->element;

  return(bond);
}



void delete_atom(MOLECULE *molecule,int id) {

  int i;
  ATOM *atom;

  for (i=id,atom=molecule->atom+id;i<molecule->natoms-1;i++,atom++) {

    memcpy(atom,atom+1,sizeof(ATOM));
  }

  molecule->natoms--;
}



void delete_bond(MOLECULE *molecule,int id) {

  int i;
  BOND *bond;

  for (i=id,bond=molecule->bond+id;i<molecule->nbonds-1;i++,bond++) {

    memcpy(bond,bond+1,sizeof(BOND));
  }

  molecule->nbonds--;
}



MAP* system2atommap(SYSTEM *system,double spacing) {

  GRID *grid;
  MAP *map;

  grid = system2grid(system,spacing,0.0);

  if (grid == NULL) {

    return(NULL);
  }

  map = new_map("PLI system atom map");

  map->grid = grid;

  map->type = ATOMLIST_MAP;

  alloc_map_matrix(map);

  if (system->protein != NULL) {

    map = molecule2atommap(system->protein,spacing,map,0);
  }

  if (system->ligand != NULL) {

    map = molecule2atommap(system->ligand,spacing,map,0);
  }

  if (system->symmetry != NULL) {

    map = molecule2atommap(system->symmetry,spacing,map,0);
  }

  return(map);
}



GRID* system2grid(SYSTEM *system,double spacing,double padding) {

  GRID *grid = NULL;
  
  if (system->protein != NULL) {

    grid = molecule2grid(system->protein,spacing,padding,grid);
  }

  if (system->ligand != NULL) {

    grid = molecule2grid(system->ligand,spacing,padding,grid);
  }

  if (system->symmetry != NULL) {

    grid = molecule2grid(system->symmetry,spacing,padding,grid);
  }

  return(grid);
}



MAP* molecule2atommap(MOLECULE *molecule,double spacing,MAP *map,int set_gridlists) {

  int i,j,iv[3],flag;
  double *v;
  ATOM *atom;
  GRID *grid;
  ATOMLIST ***matrix,*al;

  if (map == NULL) {

    map = new_map("PLI molecule atom map");

    map->grid = molecule2grid(molecule,spacing,3.0,NULL);

    map->type = ATOMLIST_MAP;

    alloc_map_matrix(map);
  }

  grid = map->grid;

  matrix = (ATOMLIST***) map->matrix;

  for (i=0,atom=molecule->atom;i<molecule->natoms;i++,atom++) {

    v = atom->position;

    flag = pos2grid(v,iv,grid);

    if (flag != 0) {

      write_pdb_atom(PLI_STDERR,atom);

      error_fn("molecule2atommap: atom oustide grid");
    }

    al = (ATOMLIST*) &(matrix[iv[0]][iv[1]][iv[2]]);

    add_atom_to_list(al,atom,0);

    if (set_gridlists) {

      atom->gridmap = map;
      atom->gridlist = al;
    }
  }

  return(map);
}



GRID* molecule2grid(MOLECULE *molecule,double spacing,double padding,GRID *grid) {

  int i;
  ATOM *atom;

  if (grid == NULL) {

    grid = (GRID*) malloc(sizeof(GRID));

    if (grid == NULL) {

      error_fn("molecule2grid: out of memory allocating grid");
    }

    init_grid(grid,spacing,padding);
  }

  for (i=0,atom=molecule->atom;i<molecule->natoms;i++,atom++) {

    update_gridlimits(grid,atom->position);
  }

  grid->npoints[0] = grid->limit[0][1] - grid->limit[0][0] + 1;
  grid->npoints[1] = grid->limit[1][1] - grid->limit[1][0] + 1;
  grid->npoints[2] = grid->limit[2][1] - grid->limit[2][0] + 1;

  return(grid);
}



GRID* atomlist2grid(ATOMLIST *list,double spacing,double padding,GRID *gridin) {

  int i;
  ATOM **atomp,*atom;
  GRID *grid;

  if (gridin == NULL) {

    grid = (GRID*) malloc(sizeof(GRID));

    if (grid == NULL) {

      error_fn("atomlist2grid: out of memory allocating grid");
    }

    init_grid(grid,spacing,padding);

  } else {

    grid = gridin;
  }

  for (i=0,atomp=list->atom;i<list->natoms;i++,atomp++) {

    atom = *atomp;

    if ((!(atom->flags & SKIP_ATOM)) || (atom->element->id == HYDROGEN)) {

      update_gridlimits(grid,(*atomp)->position);
    }
  }

  grid->npoints[0] = grid->limit[0][1] - grid->limit[0][0] + 1;
  grid->npoints[1] = grid->limit[1][1] - grid->limit[1][0] + 1;
  grid->npoints[2] = grid->limit[2][1] - grid->limit[2][0] + 1;

  return(grid);
}




void move_atom(ATOM *atom,double *position) {

  int i;

  // move the atom:

  for (i=0;i<3;i++) {

    atom->position[i] = position[i];
  }

  // update grid:

  update_atom_gridlist(atom);

  atom->status |= ATOM_MOVED;

  atom->cstatus = ATOM_NOTHING_CALCULATED;
}



void shift_atom(ATOM *atom,double *shift) {

  int i;
 
  // shift the atom:

  for (i=0;i<3;i++) {

    atom->position[i] += shift[i];
  }

  // update grid:

  update_atom_gridlist(atom);

  atom->status |= ATOM_MOVED; 

  atom->cstatus = ATOM_NOTHING_CALCULATED;
}



void transform_atom(ATOM *atom,double m[4][4]) {

  int i;
  double *pos,old_pos[4],u[4],v[4],w[4],*vpt;
  ATOM_GEOMETRY *geometry;
  HBOND_GEOMETRY *hbond_geometry;

  pos = atom->position;

  // save original coordinates for transformation of u and v:
 
  for (i=0;i<4;i++) {

    old_pos[i] = atom->position[i];
  }

  // transform the atom:

  transform_vector(atom->position,m);

  // transform u and v:

  geometry = atom->geometry;

  if (geometry) {

    if (geometry->u_axis) {

      sum_vector(old_pos,atom->u,u);

      transform_vector(u,m);

      calc_vector(pos,u,atom->u);
    } 

    if (geometry->v_axis) {

      sum_vector(old_pos,atom->v,v);

      transform_vector(v,m);

      calc_vector(pos,v,atom->v);

      sum_vector(old_pos,atom->w,w);

      transform_vector(w,m);

      calc_vector(pos,w,atom->w);
    }
  }

  // transform virtual points:

  hbond_geometry = atom->hbond_geometry;

  if (hbond_geometry) {

    for (i=0;i<hbond_geometry->npts;i++) {

      vpt = atom->vpts[i];

      sum_vector(old_pos,vpt,v);

      transform_vector(v,m);

      calc_vector(pos,v,vpt);
    }
  }

  // update grid:

  update_atom_gridlist(atom);

  atom->status |= ATOM_MOVED; 

  atom->cstatus = ATOM_NOTHING_CALCULATED;
}



void select_molecule_atoms(MOLECULE *molecule,char *selstr) {

  if (molecule->flags & PROTEIN_MOLECULE) {

    if ((!strcmp(selstr,"protein")) || (!strcmp(selstr,"all"))) {

      molecule->selection = molecule2atoms(molecule);

    } else if (!strcmp(selstr,"water")) {

      molecule->selection = molecule2waters(molecule);

    } else if ((!strcmp(selstr,"complex")) || (!strcmp(selstr,"site"))) {

      molecule->selection = molecule2site(molecule);
    }

  } else if (molecule->flags & LIGAND_MOLECULE) {

    if ((!strcmp(selstr,"ligand")) || (!strcmp(selstr,"complex")) || (!strcmp(selstr,"all"))) {

      molecule->selection = molecule2atoms(molecule);
    }
  }
}



static void update_atom_gridlist(ATOM *atom) {

  int i,iv[3],flag,n1,n2,j,k,nx,ny,nz;
  GRID *grid;
  MAP *gridmap;
  ATOMLIST ***matrix,*al,*l,***am;

  // remove atom from gridlist:

  if (atom->gridlist) {

    flag = remove_atom_from_list(atom->gridlist,atom);

    atom->gridlist = NULL;
  }

  // add atom to gridlist at new position:

  gridmap = atom->gridmap;

  if (gridmap) {

    grid = gridmap->grid;
    matrix = (ATOMLIST***) gridmap->matrix;

    flag = pos2grid(atom->position,iv,grid);

    if (flag == 0) {

      al = (ATOMLIST*) &(matrix[iv[0]][iv[1]][iv[2]]);

      add_atom_to_list(al,atom,0);

      atom->gridlist = al;

    } else {

      error_fn("update_atom_gridlist: cannot find atom grid position for atom %d (%s) in molecule %s\n",atom->id,atom->name,atom->molecule->name);
    }
  }

  if ((!atom->gridlist) || ((atom->gridlist) && (!atom_in_list(atom->gridlist,atom)))) {

  }
}



void molecule_center(MOLECULE *molecule,double *center) {

  int i,j,n;
  double *v;
  ATOM *atom;

  null_vector(center);

  n = 0;

  for (i=0,atom=molecule->atom;i<molecule->natoms;i++,atom++) {

    if (atom->element->id != HYDROGEN) {

      v = atom->position;

      for (j=0;j<3;j++) {

        center[j] += v[j];
      }

      n++;
    }
  }

  for (i=0;i<3;i++) {

    center[i] /= n;
  }

  center[3] = 1.0;
}



ATOMLIST* molecule_sphere_to_atom_list(MOLECULE *molecule1,MOLECULE *molecule2,double max_dist) {

  int i,j,contact;
  ATOM *atom1,*atom2;
  ATOMLIST *list;

  list = alloc_atomlist();

  init_atomlist(list);

  for (i=0,atom1=molecule1->atom;i<molecule1->natoms;i++,atom1++) {

    if ((atom1->element != NULL) && (atom1->element->id != HYDROGEN)) {

      contact = 0;

      for (j=0,atom2=molecule2->atom;j<molecule2->natoms;j++,atom2++) {

	if (!contact) {

	  if ((atom2->element != NULL) && (atom2->element->id != HYDROGEN)) {

	    if (points_within_distance(atom1->position,atom2->position,max_dist)) {

	      contact = 1;
	    }
	  }
	}
      }

      if (contact) {

	add_atom_to_list(list,atom1,0);
      }
    }
  }

  return(list);
}



ATOMLIST* molecule_atom_ids_to_atom_list(MOLECULE *molecule,int *atom_ids,int n_atom_ids) {

  int i,j,*atom_id;
  ATOM *atom;
  ATOMLIST *list;

  list = alloc_atomlist();

  init_atomlist(list);

  for (i=0,atom=molecule->atom;i<molecule->natoms;i++,atom++) {

    if ((atom->element != NULL) && (atom->element->id != HYDROGEN)) {

      for (j=0,atom_id=atom_ids;j<n_atom_ids;j++,atom_id++) {

	if (*atom_id == atom->id) {

	  add_atom_to_list(list,atom,0);
	}
      }
    }
  }

  return(list);
}



MOLECULE_LIST* alloc_molecule_list(void) {

  MOLECULE_LIST *list;

  list = (MOLECULE_LIST*) malloc(sizeof(MOLECULE_LIST));

  if (list == NULL) {

    error_fn("alloc_molecule_list: out of memory allocating list");
  }

  return(list);
}



void init_molecule_list(MOLECULE_LIST *list) {

  list->n_molecules = 0;
  list->n_alloc_molecules = 0;
  list->molecules = NULL;
  list->file = NULL;
}



void add_molecule_to_list(MOLECULE_LIST *list,MOLECULE* molecule) {

  if (molecule == NULL) {

    return;
  }

  if (list->n_alloc_molecules == 0) {

    list->n_alloc_molecules = 10;

    list->molecules = (MOLECULE**) calloc(list->n_alloc_molecules,sizeof(MOLECULE*));

  } else if (list->n_molecules == list->n_alloc_molecules) {

    list->n_alloc_molecules *= 2;

    list->molecules = (MOLECULE**) realloc(list->molecules,(list->n_alloc_molecules)*sizeof(MOLECULE*));
  }

  if (list->molecules == NULL) {

    error_fn("add_molecule_to_list: out of memory allocating molecule_list");
  }

  list->molecules[list->n_molecules] = molecule;

  list->n_molecules++;
}



int remove_molecule_from_list(MOLECULE_LIST *list,MOLECULE *molecule) {

  int i,j;
  MOLECULE **lmol;

  if (list->n_molecules == 0) {

    return(1);
  }

  i = 0;
  lmol = list->molecules;

  while ((*lmol != molecule) && (i < list->n_molecules)) {

    i++;

    if (i < list->n_molecules) {

      lmol++;
    }
  }

  if (i < list->n_molecules) {

    for (j=i;j<list->n_molecules-1;j++,lmol++) {

      *lmol = *(lmol+1);
    }

    list->n_molecules--;

    return(0);
  }

  return(1);
}



void free_molecule_list(MOLECULE_LIST *list) {

  if (list) {

    if (list->molecules) {

      free(list->molecules);
    }

    if (list->file) {

      close_file(list->file);
    }

    free(list);
  }
}



ATOMLIST* alloc_atomlist(void) {

  ATOMLIST *list;

  list = (ATOMLIST*) malloc(sizeof(ATOMLIST));

  if (list == NULL) {

    error_fn("alloc_atomlist: out of memory allocating list");
  }

  return(list);
}



void init_atomlist(ATOMLIST *list) {

  list->natoms = 0;
  list->n_alloc_atoms = 0;
  list->atom = NULL;
  list->flags = NULL;
}



void add_atom_to_list(ATOMLIST *list,ATOM* atom,unsigned int flag) {

  if (list->n_alloc_atoms == 0) {

    list->n_alloc_atoms = 10;

    list->atom = (ATOM**) calloc(list->n_alloc_atoms,sizeof(ATOM*));

    if (list->atom != NULL)
      list->flags = (unsigned int*) calloc(list->n_alloc_atoms,sizeof(unsigned int));

  } else if (list->natoms == list->n_alloc_atoms) {

    list->n_alloc_atoms *= 2;

    list->atom = (ATOM**) realloc(list->atom,(list->n_alloc_atoms)*sizeof(ATOM*));

    if (list->atom != NULL)
      list->flags = (unsigned int*) realloc(list->flags,(list->n_alloc_atoms)*sizeof(unsigned int));
  }

  if ((list->atom == NULL) || (list->flags == NULL)) {

    error_fn("add_atom_to_list: out of memory allocating atomlist");
  }

  list->atom[list->natoms] = atom;
  list->flags[list->natoms] = flag;

  list->natoms++;
}



int remove_atom_from_list(ATOMLIST *list,ATOM *atom) {

  int i,j;
  unsigned int *flag;
  ATOM **latom;

  if (list->natoms == 0) {

    return(1);
  }

  i = 0;
  latom = list->atom;
  flag = list->flags;

  while ((*latom != atom) && (i < list->natoms)) {

    i++;

    if (i < list->natoms) {

      latom++;
      flag++;
    }
  }

  if (i < list->natoms) {

    for (j=i;j<list->natoms-1;j++,latom++,flag++) {

      *latom = *(latom+1);
      *flag = *(flag+1);
    }

    list->natoms--;

    return(0);

  }

  return(1);
}



int atom_in_list(ATOMLIST *list,ATOM *atom) {

  int i;
  ATOM **latom;

  for (i=0,latom=list->atom;i<list->natoms;i++,latom++) {

    if (*latom == atom) {

      return(1);
    }
  }

  return(0);
}



int atomlist_index(ATOMLIST *list,ATOM *atom) {

  int i;
  ATOM **latom;

  for (i=0,latom=list->atom;i<list->natoms;i++,latom++) {

    if (*latom == atom) {

      return(i);
    }
  }

  return(-1);
}



ATOMLIST* atom2atomlist(ATOM *atom,unsigned int flag) {

  ATOMLIST *list;

  list = alloc_atomlist();

  init_atomlist(list);

  add_atom_to_list(list,atom,flag);

  return(list);
}



ATOMLIST* molecule2atoms(MOLECULE *molecule) {

  ATOMLIST *list;

  list = (ATOMLIST*) alloc_atomlist();

  if (list == NULL) {

    error_fn("molecule2atoms: out of memory allocating list");
  }

  init_atomlist(list);

  add_molecule_atoms_to_list(list,molecule);

  return(list);
}



ATOMLIST* molecule2waters(MOLECULE *molecule) {

  int i;
  ATOM *atom;
  ATOMLIST *list;

  list = (ATOMLIST*) alloc_atomlist();

  if (list == NULL) {

    error_fn("molecule2waters: out of memory allocating list");
  }

  init_atomlist(list);

  if ((molecule == NULL) || (molecule->atom == NULL)) {

    error_fn("molecule2waters: no molecule (or atoms)");
  }

  for (i=0,atom=molecule->atom;i<molecule->natoms;i++,atom++) {

    if (atom->flags & WATER_OXYGEN) {

      add_atom_to_list(list,atom,0);
    }
  }

  return(list);
}



ATOMLIST* atomlist2waters(ATOMLIST *list) {

  int i;
  ATOM **atomp,*atom;
  ATOMLIST *wlist;

  wlist = (ATOMLIST*) alloc_atomlist();

  if (wlist == NULL) {

    error_fn("atomlist2waters: out of memory allocating wlist");
  }

  init_atomlist(wlist);

  if (list) {

    for (i=0,atomp=list->atom;i<list->natoms;i++,atomp++) {

      atom = *atomp;

      if (atom->flags & WATER_OXYGEN) {

	add_atom_to_list(wlist,atom,0);
      }
    }
  }

  return(wlist);
}



void add_molecule_atoms_to_list(ATOMLIST *list,MOLECULE *molecule) {

  int i;
  ATOM *atom;

  if (list == NULL) {

    error_fn("add_molecule_atoms_to_list: no list");
  }

  if ((molecule == NULL) || (molecule->atom == NULL)) {

    error_fn("add_molecule_atoms_to_list: no molecule (or atoms)");
  }

  for (i=0,atom=molecule->atom;i<molecule->natoms;i++,atom++) {

    if ((atom->element != NULL) && (atom->element->id != HYDROGEN)) {

      add_atom_to_list(list,atom,0);
    }
  }
}



void add_atomlist_to_atomlist(ATOMLIST *list1,ATOMLIST *list2) {

  int i;
  ATOM **atom2p;

  if (list1 == NULL) {

    error_fn("add_atomlist_to_atomlist: no list");
  }

  if (list2 != NULL) {

    for (i=0,atom2p=list2->atom;i<list2->natoms;i++,atom2p++) {
      
      add_atom_to_list(list1,*atom2p,0);
    }
  }
}



ATOMLIST* list2hydrogen_free_list(ATOMLIST* list) {

  int i;
  ATOM **atom;
  ATOMLIST *hf_list;

  hf_list = (ATOMLIST*) alloc_atomlist();

  init_atomlist(hf_list);

  for (i=0,atom=list->atom;i<list->natoms;i++,atom++) {

    if (((*atom)->element != NULL) && ((*atom)->element->id != HYDROGEN)) {

      add_atom_to_list(hf_list,*atom,0);
    }
  }

  return(hf_list);
}



void reset_atom_status_list(ATOMLIST *list) {

  int i;
  ATOM **atomp;

  if (list == NULL) {

    return;
  }

  for (i=0,atomp=list->atom;i<list->natoms;i++,atomp++) {

    (*atomp)->status = ATOM_MOVED;

    (*atomp)->cstatus = ATOM_NOTHING_CALCULATED;
  }
}



void free_atomlist(ATOMLIST* list) {

  if (list != NULL) {

    if (list->n_alloc_atoms) {

      free(list->atom);
      free(list->flags);
    }

    free(list);
  }
}



void init_bondlist(BONDLIST *list) {

  list->nbonds = 0;
  list->n_alloc_bonds = 0;
  list->bond = NULL;
}



void add_bond_to_list(BONDLIST *list,BOND* bond) {

  if (list->n_alloc_bonds == 0) {

    list->n_alloc_bonds = 10;

    list->bond = (BOND**) calloc(list->n_alloc_bonds,sizeof(BOND*));

  } else if (list->nbonds == list->n_alloc_bonds) {

    list->n_alloc_bonds *= 2;

    list->bond = (BOND**) realloc(list->bond,(list->n_alloc_bonds)*sizeof(BOND*));
  }

  if (list->bond == NULL) {

    error_fn("add_bond_to_list: out of memory allocating bondlist");
  }

  list->bond[list->nbonds] = bond;

  list->nbonds++;
}



int bond_in_list(BONDLIST *list,BOND *bond) {

  int i;
  BOND **lbond;

  for (i=0,lbond=list->bond;i<list->nbonds;i++,lbond++) {

    if (*lbond == bond) {

      return(1);
    }

    if ((((*lbond)->atom1 == bond->atom1) && ((*lbond)->atom2 == bond->atom2)) ||
	(((*lbond)->atom1 == bond->atom2) && ((*lbond)->atom2 == bond->atom1))) {

      return(1);
    }
  }

  return(0);
}



void free_bondlist(BONDLIST* list) {

  if (list != NULL) {

    if (list->n_alloc_bonds) {

      free(list->bond);
    }

    free(list);
  }
}



BONDLIST* atomlist2bondlist(ATOMLIST *atomlist) {

  int i,j;
  ATOM **atom,**batom;
  ATOMLIST *conns;
  BOND *bond;
  MOLECULE *mol,*bmol;
  BONDLIST *bondlist;

  bondlist = (BONDLIST*) malloc(sizeof(BONDLIST));

  if (bondlist == NULL) {

    error_fn("atomlist2bondlist: out of memory allocating bond list");
  }

  init_bondlist(bondlist);

  for (i=0,atom=atomlist->atom;i<atomlist->natoms;i++,atom++) {
    
    mol = (*atom)->molecule;
    conns = (*atom)->connections;

    if (conns != NULL) {

      for (j=0,batom=conns->atom;j<conns->natoms;j++,batom++) {

	bmol = (*batom)->molecule;

	bond = get_bond(mol,*atom,*batom);

	if ((bond == NULL) && (mol != bmol)) {

	  bond = get_bond(bmol,*atom,*batom);
	}

	if ((bond != NULL) && (!bond_in_list(bondlist,bond)) && (atom_in_list(atomlist,*batom))) {

	  add_bond_to_list(bondlist,bond);
	}
      }
    }
  }

  return(bondlist);
}



void set_system_connections(SYSTEM *system) {

  int i,j;
  ATOM *atom;
  MOLECULE **moleculep,*molecule;
  MOLECULE_LIST *list;

  list = system->molecule_list;

  if (list == NULL) {

    return;
  }

  if (system->ligand) {

    check_molecule_dict(system->ligand);
  }

  for (i=0,moleculep=list->molecules;i<list->n_molecules;i++,moleculep++) {

    molecule = *moleculep;

    for (j=0,atom=molecule->atom;j<molecule->natoms;j++,atom++) {

      set_atom_connections(atom,0);
    }
  }
}



void set_molecule_connections(MOLECULE *molecule,int silent) {

  int i;
  ATOM *atom;

  for (i=0,atom=molecule->atom;i<molecule->natoms;i++,atom++) {

    set_atom_connections(atom,silent);
  }
}



double molecule_dict_match(MOLECULE *molecule) {

  int i,n_atoms,n_dict_failures;
  ATOM *atom;

  n_atoms = n_dict_failures = 0;

  for (i=0,atom=molecule->atom;i<molecule->natoms;i++,atom++) {

    if (atom->element->id != HYDROGEN) {

      if (atom->error_flags & ATOM_DICTIONARY_ERROR) {

	n_dict_failures++;
      }

      n_atoms++;
    }
  }

  if (n_atoms == 0) {

    return(1.0);
  }

  return(((double) (n_atoms - n_dict_failures))/((double) n_atoms));
}



void reset_molecule_connections(MOLECULE *molecule) {

  int i;
  ATOM *atom;

  for (i=0,atom=molecule->atom;i<molecule->natoms;i++,atom++) {

    if (atom->element->id != HYDROGEN) {  

      if (atom->connections != NULL) {

	free_atomlist(atom->connections);

	atom->connections = NULL;

	atom->n_hydrogens = 0;
      }
    }
  }
}



static void check_molecule_dict(MOLECULE *molecule) {

  int i,n_atoms,n_dict_failures;
  ATOM *atom;

  if ((molecule == NULL) || (!molecule->use_pdb_dict)) {

    return;
  }

  n_atoms = n_dict_failures = 0;

  for (i=0,atom=molecule->atom;i<molecule->natoms;i++,atom++) {

    if (atom->element->id != HYDROGEN) {

      set_atom_connections(atom,1);

      if (atom->error_flags & ATOM_DICTIONARY_ERROR) {

	n_dict_failures++;
      }

      n_atoms++;
    }
  }

  if (n_atoms == 0) {

    return;
  }

  if (((double) n_dict_failures)/((double) n_atoms) > 0.20) {

    warning_fn("%.1f%% of molecule atoms could not be matched to dictionary; switching dictionary off for molecule",
	       100*((double) n_dict_failures)/((double) n_atoms));


    molecule->use_pdb_dict = 0;
  }

  for (i=0,atom=molecule->atom;i<molecule->natoms;i++,atom++) {

    if (atom->element->id != HYDROGEN) {  

      if (atom->connections != NULL) {

	free_atomlist(atom->connections);

	atom->connections = NULL;

	if (molecule->use_pdb_dict == 0) {

	  atom->error_flags = ATOM_SKIP_DICTIONARY;
	}
      }
    }
  }
}



static void set_atom_connections(ATOM *atom1,int silent) {

  int i,j,dict_atom_ok;
  int nbonds,dnbonds;
  ATOM **atom2,*datom1,*batom;
  BOND *bond;
  MOLECULE *molecule1,*dict;

  if (atom1->connections == NULL) {

    atom1->connections = (ATOMLIST*) alloc_atomlist();

    init_atomlist(atom1->connections);
  }

  if (atom1->element->id == HYDROGEN) {

    return;
  }

  molecule1 = atom1->molecule;

  if ((molecule1 != NULL) && (molecule1->use_pdb_dict) && (!(atom1->flags & AMINO_ACID_ATOM)) &&
      (!(atom1->flags & WATER_OXYGEN)) && (!atom1->element->flags & METAL_ELEMENT)) {

    dict = get_pdb_dict(atom1->subname);

    if ((dict != NULL) && (dict->natoms == 0)) {

      dict = NULL;
    }

  } else {

    dict = NULL;
  }

  datom1 = NULL;

  dict_atom_ok = 1;

  if (dict != NULL) {

    datom1 = get_cif_atom(atom1->name,dict);

    if (datom1 == NULL) {

      dict_atom_ok = 0;

      atom1->error_flags |= ATOM_DICTIONARY_ERROR;

      if (!silent) {

	warning_fn("set_atom_connections: atom1 '%s' (%d) not found in dictionary '%s'",atom1->name,atom1->id,dict->name);
      }

    } else if (datom1->element != atom1->element) {
      
      dict_atom_ok = 0;

      atom1->error_flags |= ATOM_DICTIONARY_ERROR;

      if (!silent) {

	warning_fn("set_atom_connections: element mismatch for atom1 '%s' (%d) against dictionary '%s'",atom1->name,atom1->id,dict->name);
      }
    }
  }

  molecule2connections(atom1,molecule1,dict,datom1,dict_atom_ok,silent);

  // now remove any bonds that are not in list of connected atoms:
  // these are bonds that were defined in, e.g. the CONECT records, but werent found in the structure

  if (!silent) {

    if (molecule1 != NULL) {

      i = 0;

      while (i < molecule1->nbonds) {

	bond = molecule1->bond + i;

	batom = (bond->atom1 == atom1) ? bond->atom2 : (bond->atom2 == atom1) ? bond->atom1 : NULL;

	if ((batom != NULL) && (!atom_in_list(atom1->connections,batom))) {

	  delete_bond(molecule1,i);

	} else {

	  i++;
	}
      }
    }
  }

  // now check for any bonds that are in the dictionary, but missing here:

  if ((dict != NULL) && (datom1 != NULL)) {

    dnbonds = 0;

    for (i=0,bond=dict->bond;i<dict->nbonds;i++,bond++) {

      batom = (bond->atom1 == datom1) ? bond->atom2 : (bond->atom2 == datom1) ? bond->atom1 : NULL;

      if ((batom != NULL) && (batom->element->id != HYDROGEN)) {

	dnbonds++;
      }
    }

    nbonds = 0;

    for (i=0,atom2=atom1->connections->atom;i<atom1->connections->natoms;i++,atom2++) {

      if (((*atom2)->element->id != HYDROGEN) && (same_residue_atoms(atom1,*atom2))) {

	nbonds++;
      }
    }

    if (nbonds != dnbonds) {

      atom1->error_flags |= ATOM_DICTIONARY_ERROR;

      if (!silent) {

	warning_fn("set_atom_connections: atom '%s' (%d) connectivity differs from dictionary '%s' (%d)",
		   atom1->name,nbonds,dict->name,dnbonds);
      }
    }
  }
}



static void molecule2connections(ATOM *atom1,MOLECULE *molecule,MOLECULE *dict,ATOM *datom1,int dict_atom_ok,int silent) {

  int i;
  ATOM *atom2;
  MAP *covalent_map;

  if ((atom1->molecule == molecule) && (molecule->rely_on_bonds)) {

    bonds2connections(atom1,molecule);

  } else {

    covalent_map = molecule->covalent_map;

    if (covalent_map) {

      map2connections(atom1,covalent_map,dict,datom1,dict_atom_ok,silent);

    } else {

      for (i=0,atom2=molecule->atom;i<molecule->natoms;i++,atom2++) {

	pair2connections(atom1,atom2,dict,datom1,dict_atom_ok,silent);
      }
    }
  }
}



static void bonds2connections(ATOM *atom,MOLECULE *molecule) {

  int i;
  ATOM *atom1,*atom2;
  ATOMLIST *connections;
  BOND *bond;

  connections = atom->connections;

  for (i=0,bond=molecule->bond;i<molecule->nbonds;i++,bond++) {

    atom1 = bond->atom1;
    atom2 = bond->atom2;

    if (atom1 == atom) {

      if (atom2->element->id == HYDROGEN) {

	atom->n_hydrogens++;

      } else if (!atom_in_list(connections,atom2)) {

	add_atom_to_list(connections,atom2,bond->type);
      }

    } else if (atom2 == atom) {

      if (atom1->element->id == HYDROGEN) {

	atom->n_hydrogens++;

      } else if (!atom_in_list(connections,atom1)) {

	add_atom_to_list(connections,atom1,bond->type);
      }
    }
  }
}



static void map2connections(ATOM *atom1,MAP *map,MOLECULE *dict,ATOM *datom1,int dict_atom_ok,int silent) {

  int i,flag,iv[3],irange[3][2],ix,iy,iz;
  GRID *grid;
  ATOM **atom2;
  ATOMLIST ***matrix,*list;

  grid = map->grid;
  matrix = (ATOMLIST***) map->matrix;

  flag = pos2grid(atom1->position,iv,grid);

  if (flag != 0) {

    return;
  }

  for (i=0;i<3;i++) {

    irange[i][0] = (iv[i] == 0) ? 0 : iv[i] - 1;
    irange[i][1] = (iv[i]+1 == grid->npoints[i]) ? iv[i] : iv[i] + 1;
  }

  for (ix=irange[0][0];ix<=irange[0][1];ix++) {

    for (iy=irange[1][0];iy<=irange[1][1];iy++) {

      for (iz=irange[2][0];iz<=irange[2][1];iz++) {

	list = &(matrix[ix][iy][iz]);

	for (i=0,atom2=list->atom;i<list->natoms;i++,atom2++) {

	  pair2connections(atom1,*atom2,dict,datom1,dict_atom_ok,silent);
	}
      }
    }
  }
}



void pair2connections(ATOM *atom1,ATOM *atom2,MOLECULE *dict,ATOM *datom1,int dict_atom_ok,int silent) {

  int dict_bond_ok,bb_bond_type;
  double Do,D;
  ATOM *datom2;
  BOND *dbond,*bond;
  MOLECULE *molecule;

  molecule = atom1->molecule;

  if ((atom2)->element == NULL) {

    error_fn("pair2connections: no element for atom %d '%s'\n",atom2->id,atom2->name);
  }

  if (!allow_covalent_bond(atom1,atom2)) {

    return;
  }

  Do = atom1->element->cov_radius + atom2->element->cov_radius;

  D = distance(atom1->position,atom2->position);

  if (D < Do + COVALENT_TOLERANCE) {

    if (atom2->element->id == HYDROGEN) {

      atom1->n_hydrogens++;

    } else {

      dbond = NULL;

      dict_bond_ok = dict_atom_ok;

      if ((datom1 != NULL) && (same_residue_atoms(atom1,atom2))) {
	
	datom2 = get_cif_atom(atom2->name,dict);

	if (datom2 != NULL) {

	  dbond = get_bond(dict,datom1,datom2);

	  if (dbond == NULL) {

	    dict_bond_ok = 0;

	    if (!silent) {

	      warning_fn("pair2connections: set_atom_connections: no bond between atoms '%s' and '%s' in dictionary '%s'",
			 datom1->name,datom2->name,dict->name);
	    }
	  }
	  
	} else {

	  dict_bond_ok = 0;

	  if (!silent) {

	    warning_fn("pair2connections: set_atom_connections: atom2 '%s' (%d,%d) not found in dictionary '%s'",
		       atom2->name,atom1->id,atom2->id,dict->name);
	  }
	}
      }

      // add bond to list and update molecule bonds:

      if ((atom2->element != NULL) && (atom2->element->id != HYDROGEN)) {

	add_atom_to_list(atom1->connections,atom2,0);

	if (!silent) {

	  bond = get_bond(molecule,atom1,atom2);

	  if (bond == NULL) {

	    bond = add_bond(molecule,atom1,atom2,1);
	  }

	  bb_bond_type = backbone_bond_type(bond);

	  if (bb_bond_type) {

	    bond->type = bb_bond_type;
		  
	  } else if (dbond != NULL) {
		   
	    bond->type = dbond->type;
	  }
	}
		  
	if (D < Do - COVALENT_TOLERANCE) {

	  atom1->error_flags |= ATOM_GEOMETRY_ERROR;

	  if (!silent) {

	    warning_fn("pair2connections: very short bond:\n%5d %-6s%c %-3s %6d - %5d %-6s%c %-3s %6d : %6.4fA",
		       atom1->id,atom1->name,atom1->altloc,atom1->subname,atom1->subid,
		       atom2->id,atom2->name,atom2->altloc,atom2->subname,atom2->subid,D);
	  }
	}
      }

      if (!dict_bond_ok) {

	atom1->error_flags |= ATOM_DICTIONARY_ERROR;
      }
    }
  }
}



void get_connected_atoms(ATOMLIST *list,ATOM *atom,ATOM *banned_atom,int metal_bonds,unsigned int depth,unsigned int max_depth) {

  int i,j;
  ATOM **batomp,*batom,*catom;
  CONTACT *contact;

  if ((depth == max_depth) || (atom == banned_atom)) {

    return;
  }

  if (!atom_in_list(list,atom)) {

    add_atom_to_list(list,atom,depth);
  }

  if (atom->connections == NULL) {

    return;
  }

  for (i=0,batomp=atom->connections->atom;i<atom->connections->natoms;i++,batomp++) {

    batom = *batomp;

    if (!atom_in_list(list,batom)) {

      get_connected_atoms(list,batom,banned_atom,metal_bonds,depth+1,max_depth);
    }
  }

  if ((!metal_bonds) || (atom->contactlist == NULL)) {

    return;
  }

  for (i=0,contact=atom->contactlist->contacts;i<atom->contactlist->ncontacts;i++,contact++) {

    catom = (contact->atom1 == atom) ? contact->atom2 : contact->atom1;

    if ((!atom_in_list(list,catom)) && (contact->flags & METAL_COORDINATING_CONTACT)) {

      get_connected_atoms(list,catom,banned_atom,metal_bonds,depth+1,max_depth);
    }
  }
}



int same_molecule_atoms(ATOM *atom1,ATOM *atom2) {

  if (atom1->molecule != atom2->molecule) {

    return(0);
  }

  if (atom1->chain != atom2->chain) {

    return(0);
  }

  if ((atom1->flags & WATER_OXYGEN) || (atom2->flags & WATER_OXYGEN)) {

    return(0);
  }

  return(1);
}



int atom_residue_in_list(ATOM *atom1,ATOMLIST *list) {

  int i;
  ATOM **atom2p,*atom2;

  if (list == NULL) {

    return(0);
  }

  for (i=0,atom2p=list->atom;i<list->natoms;i++,atom2p++) {

    atom2 = *atom2p;

    if (same_residue_atoms(atom1,atom2)) {

      return(1);
    }
  }

  return(0);
}



int same_residue_atoms(ATOM *atom1,ATOM *atom2) {

  if ((!strcmp(atom1->subname,atom2->subname)) && (atom1->subid == atom2->subid) &&
      (atom1->chain == atom2->chain) && (atom1->icode == atom2->icode) && (atom1->altloc == atom2->altloc)) {

    return(1);
  }

  return(0);
}



ATOM* atomlist2atoms(ATOMLIST *list) {

  int i;
  ATOM *atoms,**atomp;
  
  if ((list == NULL) || (list->natoms == 0)) {

    error_fn("atomlist2atoms: no atoms");
  }

  atoms = (ATOM*) calloc(list->natoms,sizeof(ATOM));

  if (atoms == NULL) {

    error_fn("atomlist2atoms: out of memory allocating atoms");
  }

  for (i=0,atomp=list->atom;i<list->natoms;i++,atomp++) {

    memcpy(atoms+i,*atomp,sizeof(ATOM));
  }

  return(atoms);
}



void atoms2atomlist(ATOMLIST *list,ATOM *atoms) {

  int i;
  ATOM **atomp,*atom;

  for (i=0,atomp=list->atom;i<list->natoms;i++,atomp++) {

    atom = *atomp;

    memcpy(atom,atoms+i,sizeof(ATOM));

    update_atom_gridlist(atom);
    
    atom->status |= ATOM_MOVED;

    atom->cstatus = ATOM_NOTHING_CALCULATED;
  }
}



void copy_atomlist(ATOMLIST *list1,ATOMLIST *list2) {

  int i;
  unsigned int *flag;
  ATOM **atomp;

  for (i=0,atomp=list1->atom,flag=list1->flags;i<list1->natoms;i++,atomp++,flag++) {

    add_atom_to_list(list2,*atomp,*flag);
  }
}



double** atomlist2coordinates(ATOMLIST *list) {

  int i,j;
  double **v;
  ATOM **atomp,*atom;

  if ((list == NULL) || (list->natoms == 0)) {

    error_fn("atomlist2coordinates: no atoms");
  }

  v = alloc_2d_fmatrix(list->natoms,4);

  for (i=0,atomp=list->atom;i<list->natoms;i++,atomp++) {

    atom = *atomp;

    for (j=0;j<4;j++) {

      v[i][j] = atom->position[j];
    }
  }

  return(v);
}



void coordinates2atomlist(ATOMLIST *list,double **v) {

  int i;
  ATOM **atomp,*atom;

  for (i=0,atomp=list->atom;i<list->natoms;i++,atomp++) {

    atom = *atomp;

    move_atom(atom,v[i]);
  }
}



double** molecule2coordinates(MOLECULE *molecule) {

  int i,j;
  double **v;
  ATOM *atom;

  if ((molecule == NULL) || (molecule->natoms == 0)) {

    error_fn("molecule2coordinates: no atoms");
  }

  v = alloc_2d_fmatrix(molecule->natoms,4);

  for (i=0,atom=molecule->atom;i<molecule->natoms;i++,atom++) {

    for (j=0;j<4;j++) {

      v[i][j] = atom->position[j];
    }
  }

  return(v);
}



void coordinates2molecule(MOLECULE *molecule,double **v) {

  int i;
  ATOM *atom;

  for (i=0,atom=molecule->atom;i<molecule->natoms;i++,atom++) {

    move_atom(atom,v[i]);
  }
}



static int linear_atom(ATOM *atom) {

  double v1[4],v2[4];
  ATOMLIST *conns;

  conns = atom->connections;

  if (conns == NULL) {

    warning_fn("linear_atom: atom connections undefined");

    return(0);
  }

  if (conns->natoms != 2)
    return(0);
  
  calc_vector(atom->position,(*(conns->atom))->position,v1);
  calc_vector(atom->position,(*(conns->atom+1))->position,v2);

  if (vector_angle(v1,v2) > 160.0) {

    return(1);
  }

  return(0);
}



int planar_atom(ATOM *atom) {

  double v1[4],v2[4],v3[4],angle_sum;
  ATOMLIST *conns;

  if (atom->planar != -1) {

    return(atom->planar);
  }

  conns = atom->connections;

  if (conns == NULL) {

    warning_fn("planar_atom: atom connections undefined");

    return(0);
  }

  if (conns->natoms != 3) {

    atom->planar = 0;

    return(0);
  }

  calc_vector(atom->position,(*(conns->atom+0))->position,v1);
  calc_vector(atom->position,(*(conns->atom+1))->position,v2);
  calc_vector(atom->position,(*(conns->atom+2))->position,v3);

  angle_sum = vector_angle(v1,v2) + vector_angle(v1,v3) + vector_angle(v2,v3);

  if (angle_sum > 350.0) {

    atom->planar = 1;

    return(1);
  }

  atom->planar = 0;

  return(0);
}



void set_system_vdw_radii(SYSTEM *system) {

  int i,j;
  double water_vdw_radius;
  ATOM *atom;
  MOLECULE **moleculep,*molecule;
  MOLECULE_LIST *list;

  list = system->molecule_list;

  if (list == NULL) {

    return;
  }

  water_vdw_radius = system->settings->water_vdw_radius;

  for (i=0,moleculep=list->molecules;i<list->n_molecules;i++,moleculep++) {

    molecule = *moleculep;

    for (j=0,atom=molecule->atom;j<molecule->natoms;j++,atom++) {

      set_atom_vdw_radii(atom,water_vdw_radius);
    }
  }
}



void set_molecule_vdw_radii(MOLECULE *molecule,double water_vdw_radius) {

  int i;
  ATOM *atom;

  for (i=0,atom=molecule->atom;i<molecule->natoms;i++,atom++) {
    
    set_atom_vdw_radii(atom,water_vdw_radius);
  }
}



void set_atom_vdw_radii(ATOM *atom,double water_vdw_radius) {

  ELEMENT *element;
  UNITED_ATOM *uatom;
  ATOM_TYPE *type;

  if (atom->flags & ATOM_VDW_RADIUS_SET) {

    return;
  }

  atom->flags |= ATOM_VDW_RADIUS_SET;

  type = atom->type;

  if (type != NULL) {

    uatom = type->united_atom;

    if (uatom != NULL) {

      atom->vdw_radius = uatom->vdw_radius;
      atom->vdw_radius_H2O = uatom->vdw_radius + water_vdw_radius;

      return;
    }
  }

  element = atom->element;

  if (element != NULL) {

    atom->vdw_radius = element->vdw_radius;
    atom->vdw_radius_H2O = element->vdw_radius + water_vdw_radius;

    return;
  }

  warning_fn("get_atom_vdw_radius: using default vdw radius for atom %d in molecule %s\n",atom->id,atom->molecule->name);

  atom->vdw_radius = DEFAULT_VDW_RADIUS;
  atom->vdw_radius_H2O = DEFAULT_VDW_RADIUS + water_vdw_radius;
}



int allow_covalent_bond(ATOM *atom1,ATOM *atom2) {

  if ((atom1->flags & WATER_OXYGEN) || (atom2->flags & WATER_OXYGEN))
    return(0);

  if ((atom1->element->flags & METAL_ELEMENT) || (atom2->element->flags & METAL_ELEMENT))
    return(0);

  if ((atom1->altloc != ' ') && (atom2->altloc != ' ') && (atom1->altloc != atom2->altloc))
    return(0);

  if (atom1 == atom2)
    return(0);

  return(1);
}



void atomlist2centroid(ATOMLIST *list,double *c) {

  int i,j;
  double *v;
  ATOM **atomp;

  null_vector(c);

  for (i=0,atomp=list->atom;i<list->natoms;i++,atomp++) {

    v = (*atomp)->position;

    for (j=0;j<3;j++) {

      c[j] += v[j];
    }
  }

  for (i=0;i<3;i++) {

    c[i] /= list->natoms;
  }
}

