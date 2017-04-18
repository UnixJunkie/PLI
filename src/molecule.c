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



static int UNIQUE_ATOM_ID = 0;



static void set_molecule_hydrogens(MOLECULE*);
static void set_atom_hydrogens(ATOM*);
static void set_atom_hydrogens_from_bonds(ATOM*);
static void set_atom_hydrogens_from_map(ATOM*,MAP*);
static void add_atom_hydrogen(ATOM*,ATOM*);
static void set_molecule_substructures(MOLECULE*);
static void set_atom_substructure(ATOM*);

static void check_molecule_dict(MOLECULE*);
static void set_atom_connections(ATOM*,int);
static void molecule2connections(ATOM*,MOLECULE*,MOLECULE*,ATOM*,int,int);
static void bonds2connections(ATOM*,MOLECULE*);
static void map2connections(ATOM*,MAP*,MOLECULE*,ATOM*,int,int);
static void pair2connections(ATOM*,ATOM*,MOLECULE*,ATOM*,int,int);
static void update_atom_gridlist(ATOM*);
static void add_atom_type_to_list(ATOM_TYPE*, ATOM_TYPE_LIST*);
static void set_molecule_edges(MOLECULE*);
static void set_atom_edges(ATOM*);
static ATOMLIST* molecule_atom_ids2atoms(MOLECULE*);
static ATOMLIST* molecule2shell(MOLECULE*,MOLECULE*,double);
static void select_molecule_atoms(MOLECULE*,char*);
static MOLECULE_PROBE* create_mol_probe_from_file(char*, SETTINGS*, INT_LIST*);



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
  atom->lnPu = 0.0;
  atom->vdw_radius = DEFAULT_VDW_RADIUS;
  atom->vdw_radius_H2O = DEFAULT_VDW_RADIUS + 1.4;

  atom->node = NULL;
  atom->type = NULL;
  atom->ambiguous_type = NULL;
  atom->gridmap = NULL;
  atom->gridlist = NULL;
  atom->connections = NULL;
  atom->h_connections = NULL;
  atom->contactlist = NULL;
  atom->geometry = NULL;
  atom->hbond_geometry = NULL;
  atom->coordination = NULL;
  atom->ring = NULL;
  atom->molecule = NULL;
  atom->substructure = NULL;
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
  //
  atom->tripletlist = NULL;
  atom->type_id_file = -1;

  atom->hybridisation = 0;
  atom->query = NULL;
  atom->fgroup = NULL;

  atom->depth_map = NULL;

  atom->hydrogens = NULL;
}



void unprep_atom(ATOM *atom) {

  free_atomlist(atom->connections);
  free_atomlist(atom->h_connections);
  free_contactlist(atom->contactlist);

  if (atom->coordination) {

    free(atom->coordination);
  }

  if (atom->ring) {
    free_ring(atom->ring);
    atom->ring = NULL;
  }

  if (atom->query) {

    free(atom->query);
  }

  atom->node = NULL;
  atom->type = NULL;
  atom->geometry = NULL;
  atom->hbond_geometry = NULL;

  atom->connections = NULL;
  atom->h_connections = NULL;
  atom->contactlist = NULL;
  atom->coordination = NULL;
  atom->ring = NULL;

  atom->query = NULL;

  atom->n_edges = 0;

  atom->hydrogens = NULL;
}



void init_bond(BOND *bond) {

  bond->type = -1;
  bond->atom1 = NULL;
  bond->atom2 = NULL;
  bond->torsion = NULL;

  bond->flags = 0;
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

  molecule->n_alloc_bond_angles = 0;

  if (n) {

    molecule->atom = (ATOM*) calloc(molecule->n_alloc_atoms,sizeof(ATOM));
    molecule->bond = (BOND*) calloc(molecule->n_alloc_bonds,sizeof(BOND));

  } else {

    molecule->atom = NULL;
    molecule->bond = NULL;
  }

  molecule->bond_angles = NULL;
  molecule->n_bond_angles = 0;

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
  molecule->conmap = NULL;
  molecule->depth_map = NULL;

  molecule->atom_ids = NULL;
  molecule->selection = NULL;
  molecule->active_atoms = NULL;
  molecule->molsystem = NULL;

  molecule->n_active_waters = 0;
  molecule->n_alloc_active_waters = 0;
  molecule->active_waters = NULL;

  molecule->score = 0.0;

  molecule->system = NULL;

  molecule->fgroups = NULL;
  molecule->feature_lists = NULL;
  molecule->features = NULL;

  molecule->nonh_atoms = NULL;
  molecule->hydrogens = NULL;
  molecule->hydrogen_lists = NULL;

  molecule->substructures = NULL;
  molecule->subatoms = NULL;
  molecule->atom_queries = NULL;
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

  // molecule active atoms:

  molecule->active_atoms = molecule_active_atoms(molecule);

  if (molecule->use_grids) {

    molecule->covalent_map = molecule2atommap(molecule,settings->max_covalent_dist,NULL,0);
    molecule->contacts_map = molecule2atommap(molecule,settings->max_contact_dist,NULL,1);
  }

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

  set_molecule_edges(molecule);

  set_molecule_hydrogens(molecule);

  set_molecule_substructures(molecule);

  // set nodes, atom types and radii:

  set_molecule_atom_nodes(molecule,settings->atom_typing_scheme);

  // set rings for ligands
  if (molecule->flags & LIGAND_MOLECULE) {

    set_molecule_rings(molecule);
  }

  set_molecule_atom_types(molecule,settings->atom_typing_scheme);
  set_molecule_vdw_radii(molecule,settings->water_vdw_radius);
  set_molecule_atom_geometries(molecule,settings->atom_typing_scheme);
  
  // create molecule system (for molecule in isolation):

  molecule->molsystem = molecule2system(molecule,settings);

  // set special atom types (heme nitrogens for now):

  atom_type_heme_nitrogens(molecule->molsystem);

  // molecule torsions:

  if (molecule->flags & LIGAND_MOLECULE) {

    set_molecule_internal(molecule);
  }
  
  // atomic mobility factors:

  if ((molecule->flags & PROTEIN_MOLECULE) && ((params_get_parameter("ff_ami_scaling"))->value.i)) {

    calc_molecule_ami(molecule,settings);
  }

  // apply shake if requested:

  if (((molecule->flags & PROTEIN_MOLECULE) && ((params_get_parameter("shake_protein"))->value.i)) ||
      ((molecule->flags & LIGAND_MOLECULE) && ((params_get_parameter("shake_ligand"))->value.i))) {

    shake_molecule(molecule);
  }

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
  free_map(molecule->conmap);

  molecule->covalent_map = NULL;
  molecule->contacts_map = NULL;
  molecule->conmap = NULL;

  for (i=0,atom=molecule->atom;i<molecule->natoms;i++,atom++) {

    unprep_atom(atom);
  }

  if (molecule->bond_angles) {

    free(molecule->bond_angles);
  }

  if (molecule->torsions) {

    for (i=0,torsion=molecule->torsions;i<molecule->n_torsions;i++,torsion++) {

      free_atomlist(torsion->alist1);
      free_atomlist(torsion->alist2);

      torsion->rotlist = NULL;
    }

    free(molecule->torsions);
  }

  if (molecule->fgroups) {

    free_list(molecule->fgroups);

    molecule->fgroups = NULL;
  }

  free_2d_imatrix(molecule->nb_atom_pairs);

  molecule->nb_atom_pairs = NULL;

  free_atomlist(molecule->selection);

  molecule->selection = NULL;

  free_molsystem(molecule->molsystem);

  molecule->molsystem = NULL;

  molecule->bond_angles = NULL;
  molecule->n_bond_angles = 0;
  molecule->n_alloc_bond_angles = 0;

  molecule->torsions = NULL;
  molecule->n_torsions = 0;
  molecule->nb_atom_pairs = NULL;
  molecule->selection = NULL;
  molecule->molsystem = NULL;

  free_list(molecule->active_atoms);
  free_list(molecule->nonh_atoms);
  free_list(molecule->hydrogens);
  free_list(molecule->hydrogen_lists);
  free_list(molecule->fgroups);
  free_list(molecule->feature_lists);
  free_list(molecule->features);
  free_list(molecule->substructures);
  free_list(molecule->subatoms);
  free_list(molecule->atom_queries);

  molecule->active_atoms = NULL;
  molecule->nonh_atoms = NULL;
  molecule->hydrogens = NULL;
  molecule->hydrogen_lists = NULL;
  molecule->fgroups = NULL;
  molecule->feature_lists = NULL;
  molecule->features = NULL;
  molecule->substructures = NULL;
  molecule->subatoms = NULL;
  molecule->atom_queries = NULL;

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



void free_molecule(MOLECULE *molecule) {

  SETTINGS *settings;

  settings = get_settings();

  unprep_molecule(molecule,settings);

  if (molecule->atom) {

    free(molecule->atom);
  }

  if (molecule->bond) {

    free(molecule->bond);
  }

  free(molecule);
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

    map = new_map("PLI molecule atom map","atomlist",NULL);

    map->grid = molecule2grid(molecule,spacing,3.0,NULL);

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
  ATOM *atom,**atomp;
  LIST *active_atoms;

  grid = (grid) ? grid : new_grid(spacing,padding);

  // all molecule atoms are used to define grid:

  for (i=0,atom=molecule->atom;i<molecule->natoms;i++,atom++) {

    update_gridlimits(grid,atom->position);
  }

  grid->npoints[0] = grid->limit[0][1] - grid->limit[0][0] + 1;
  grid->npoints[1] = grid->limit[1][1] - grid->limit[1][0] + 1;
  grid->npoints[2] = grid->limit[2][1] - grid->limit[2][0] + 1;

  return(grid);
}



GRID* alist2grid(LIST *list,double spacing,double padding,GRID *grid) {

  int i;
  ATOM **atomp,*atom;

  grid = (grid) ? grid : new_grid(spacing,padding);

  for (i=0,atomp=(ATOM**) list->items;i<list->n_items;i++,atomp++) {

    atom = *atomp;

    // TODO: Marcel, please check
    if ((!(atom->flags & SKIP_ATOM)) && (atom->element->id != HYDROGEN)) {
    //if ((!(atom->flags & SKIP_ATOM)) || (atom->element->id == HYDROGEN)) {

      update_gridlimits(grid,(*atomp)->position);
    }
  }

  grid->npoints[0] = grid->limit[0][1] - grid->limit[0][0] + 1;
  grid->npoints[1] = grid->limit[1][1] - grid->limit[1][0] + 1;
  grid->npoints[2] = grid->limit[2][1] - grid->limit[2][0] + 1;

  return(grid);
}



GRID* atomlist2grid(ATOMLIST *list,double spacing,double padding,GRID *grid) {

  int i;
  ATOM **atomp,*atom;

  grid = (grid) ? grid : new_grid(spacing,padding);

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



void move_molecule_centre(MOLECULE *molecule, double *position) {
  int i;
  double shift[3];
  ATOM *atom;
  for (i=0; i<3; i++) {
    shift[i] = position[i] - molecule->geometric_center[i];
  }
  for (i=0,atom=molecule->atom;i<molecule->natoms;i++,atom++) {
    shift_atom(atom, shift);
  }
  copy_vector(position, molecule->geometric_center);
}

double molecule_radius(MOLECULE *molecule) {
  int i;
  ATOM *atom;
  double dist, molecule_radius;
  int atom_found;

  atom_found = 0;
  add_molecule_center(molecule, 0);
  molecule_radius = 0.0;

  for (i=0, atom=molecule->atom; i<molecule->natoms; i++, atom++) {
    if ((atom->element->id == HYDROGEN) || (atom->flags & SKIP_ATOM)){
      continue;
    }
    dist = distance(molecule->geometric_center, atom->position);
    if (dist > molecule_radius) {
      molecule_radius = dist;
      atom_found = 1;
    }
  }
  if (! atom_found) {
    error_fn("%s: There is no selected atom in the molecule", __func__);
  }
  return molecule_radius;
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

void shift_molecule(MOLECULE *molecule, double *shift) {
  int i, x;
  ATOM *atom;
  ATOM_GEOMETRY *geometry;
  HBOND_GEOMETRY *hbond_geometry;
  double *vpt;

  for (i=0,atom=molecule->atom;i<molecule->natoms;i++,atom++) {
    shift_atom(atom, shift);
  }
}

void skip_atom_ids(MOLECULE *molecule, INT_LIST *skip_atom_ids) {
  int i, j;
  ATOM *atom;
  for (i=0, atom=molecule->atom; i<molecule->natoms; i++, atom++) {
    for (j=0; j<skip_atom_ids->n_items; j++) {
      if (skip_atom_ids->items[j] == atom->id) {
	//printf("skiping atom %d\n", atom->id);
	atom->flags |= SKIP_ATOM;
	break;
      }
    }
  }
}

void orient_molecule(MOLECULE *molecule, double *lebedev_point, double rotation_angle) {
  int i;
  ATOM *atom;
  double shift[4], angle, rot_mat[4][4], cp[4];
  add_molecule_center(molecule, 0);

  double coords_origin[4] = {0.0, 0.0, 0.0, 1.0};
  move_molecule_centre(molecule, coords_origin);

  // angle b/w lebedev point and the axis (center->atom1)
  if (distance(molecule->atom->position, molecule->geometric_center) > 1e-3) {
    angle = vector_angle(molecule->atom->position, lebedev_point);
  // If the first atom is on the center
  } else {
    angle = vector_angle((molecule->atom+1)->position, lebedev_point);
  }

  //TODO: move to rotation_matrix function
  if (angle > 0.0) {
    if (abs(angle-180.0) < 0.0001) {
      scalar_matrix(-1.0, rot_mat);
    } else {
      calc_crossproduct(molecule->atom->position, lebedev_point, cp);
      rotation_matrix(molecule->geometric_center, cp, angle, rot_mat);
    }
    for (i=0,atom=molecule->atom; i<molecule->natoms; i++,atom++) {
      transform_atom(atom, rot_mat);
    }
  }

  // rotate about the lebedev axis
  if (rotation_angle > 0.0) {
    if (abs(rotation_angle-180.0) < 0.0001) {
      scalar_matrix(-1.0, rot_mat);
    } else {
      rotation_matrix(molecule->geometric_center, lebedev_point, rotation_angle, rot_mat);
    }
    for (i=0,atom=molecule->atom; i<molecule->natoms; i++,atom++) {
      transform_atom(atom, rot_mat);
    }
  }

}


ATOM_COORDS** get_molecule_orientations(MOLECULE *molecule, int n_lebedev_points, double **lebdev_axes, int n_rot_lebedev_axis, int n_dihedral_steps, int *n_orientations) {
  //PLI_FILE *tempfile;
  //tempfile = open_file("probe_orientations.sdf", "w");
  //ATOMLIST *atomlist;

  ATOM_COORDS **ppcoords, *pinitial_coords;
  int i, j, k, l, index;
  //int n_orientations;
  double lebedev_point[4];
  double lebedev_rot_angle;
  ATOM **atomp;
  double rm[4][4];
  double origin[4] = {0.0, 0.0, 0.0, 1.0};

  pinitial_coords = molecule2coords(molecule);
  if (n_dihedral_steps) {
    (*n_orientations) = n_lebedev_points * n_rot_lebedev_axis * n_dihedral_steps;
  } else {
    (*n_orientations) = n_lebedev_points * n_rot_lebedev_axis;
  }

  ppcoords = (ATOM_COORDS**) calloc((*n_orientations), sizeof(ATOM_COORDS*));

  if (ppcoords == NULL) {
    error_fn("%s: out of memory allocating coords");
  }

  for (i=0; i<n_lebedev_points; i++) {
    for (j=0; j<n_rot_lebedev_axis; j++) {
      index = i*n_rot_lebedev_axis+j;
      memcpy(lebedev_point,*(lebdev_axes + i), 3*sizeof(double));
      lebedev_rot_angle = (360.0/(double) n_rot_lebedev_axis) * (double) j;
      orient_molecule(molecule, lebedev_point, lebedev_rot_angle); // also centers on the origin

      //
      if (n_dihedral_steps) {
        for (k=0; k<n_dihedral_steps; k++) {
          index = i*n_rot_lebedev_axis*n_dihedral_steps + j*n_dihedral_steps + k;
          rotation_matrix(molecule->torsions->atom1->position, molecule->torsions->atom2->position, 360.0 / (double) n_dihedral_steps, rm);
          for (l=0,atomp=molecule->torsions->rotlist->atom;l<molecule->torsions->rotlist->natoms;l++,atomp++) {
            transform_atom(*atomp,rm);
          }
          // center on the origin
          add_molecule_center(molecule, 0);
          move_molecule_centre(molecule, origin);
          *(ppcoords + index) = molecule2coords(molecule);
        }
      } else {
        *(ppcoords + index) = molecule2coords(molecule);
      }
      coords2molecule(pinitial_coords, molecule);
    }
  }


  // reset the center
  add_molecule_center(molecule, 0);
  free(pinitial_coords);
  //close_file(tempfile);
  return(ppcoords);
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



static void select_molecule_atoms(MOLECULE *molecule,char *selstr) {

  // symmetry molecules are part of the calculation but are not selected

  if (molecule->flags & SYMMETRY_MOLECULE) {

    return;
  }

  if (molecule->flags & PROTEIN_MOLECULE) {

    if ((!strcmp(selstr,"protein")) || (!strcmp(selstr,"all"))) {

      molecule->selection = (molecule->atom_ids) ? molecule_atom_ids2atoms(molecule) : molecule2atoms(molecule);

    } else if (!strcmp(selstr,"water")) {

      molecule->selection = (molecule->atom_ids) ? molecule_atom_ids2atoms(molecule) : molecule2waters(molecule);

    } else if ((!strcmp(selstr,"complex")) || (!strcmp(selstr,"site"))) {

      molecule->selection = (molecule->atom_ids) ? molecule_atom_ids2atoms(molecule) : molecule2site(molecule);

      if (molecule->selection == NULL) {

	error_fn("%s: failed to define site - check site ligand or site atoms",__func__);
      }
    }

  } else if (molecule->flags & LIGAND_MOLECULE) {

    if ((!strcmp(selstr,"ligand")) || (!strcmp(selstr,"complex")) || (!strcmp(selstr,"all"))) {

      molecule->selection = (molecule->atom_ids) ? molecule_atom_ids2atoms(molecule) : molecule2atoms(molecule);
    }
  }
}



LIST* molecule_active_atoms(MOLECULE *molecule) {

  ATOMLIST *alist;
  LIST *list;

  alist = NULL;

  if (molecule->atom_ids) {

    alist = molecule_atom_ids2atoms(molecule);
    
  } else {
    
    alist = molecule2site(molecule);
    
    if (!alist) {
      
      alist = molecule2atoms(molecule);
    }
  }

  list = (alist) ? atomlist2list(alist,"molecule active atoms") : NULL;

  return(list);
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



void add_molecule_center(MOLECULE *molecule, int mass_weighted) {

  int i,j,n;
  double *v;
  ATOM *atom;
  double center[4];
  null_vector(center);
  int total_mass = 0;
  n = 0;

  for (i=0,atom=molecule->atom;i<molecule->natoms;i++,atom++) {

    if ((atom->element->id != HYDROGEN) && (!(atom->flags & SKIP_ATOM))) {
      total_mass += atom->element->id;
      n++;
      v = atom->position;

      for (j=0;j<3;j++) {
	if (mass_weighted == 1) {
	  center[j] += v[j] * atom->element->id;
	} else {
	  center[j] += v[j];
	}
      }
    }
  }
  for (i=0;i<3;i++) {

    if (mass_weighted == 1) {
      center[i] /= (double) total_mass;
    } else {
      center[i] /= (double) n;
    }
  }

  center[3] = 1.0;
  if (mass_weighted == 1) {
    copy_vector(center, molecule->center_of_mass);
  } else {
    copy_vector(center, molecule->geometric_center);
  }
}


void add_molecule_moi(MOLECULE *molecule) {
  int i, j;
  ATOM *atom;
  // TODO: recalculate center of mass (com) just in case... may be add a way to check
  add_molecule_center(molecule, 1);
  double com2[4];
  null_vector(com2);
  double pos2[4];

  double moi[4];
  null_vector(moi);

  for (i=0; i<3; i++) {
    com2[i] = molecule->center_of_mass[i] * molecule->center_of_mass[i];
  }

  for (i=0,atom=molecule->atom;i<molecule->natoms;i++,atom++) {
    if (atom->element->id != HYDROGEN) {
      for (j=0; j<3; j++){
	pos2[i] = atom->position[i] * atom->position[i];
      }
      moi[0] += ((atom->id * (pos2[1] + pos2[2])) - com2[1] - com2[2]);
      moi[1] += ((atom->id * (pos2[0] + pos2[2])) - com2[0] - com2[2]);
      moi[2] += ((atom->id * (pos2[0] + pos2[1])) - com2[0] - com2[1]);
    }
  }
  copy_vector(moi, molecule->moment_of_inertia);
  //printf("moi from fn: %.2f %.2f %.2f", moi[0], moi[1], moi[2]);
}


static ATOMLIST* molecule2shell(MOLECULE *molecule1,MOLECULE *molecule2,double max_dist) {

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



static ATOMLIST* molecule_atom_ids2atoms(MOLECULE *molecule) {

  int i,*atom_id;
  ATOM *atom;
  ATOMLIST *list;

  if (!molecule->atom_ids) {

    return(NULL);
  }

  list = alloc_atomlist();

  init_atomlist(list);

  for (i=0,atom_id=(int*) molecule->atom_ids->items;i<molecule->atom_ids->n_items;i++,atom_id++) {

    atom = get_atom(molecule,*atom_id);

    if (!atom) {

      error_fn("%s: could not find atom %d in molecule '%s'",__func__,*atom_id,molecule->name);
    }

    add_atom_to_list(list,atom,0);
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



ATOMLIST* molecule2site(MOLECULE *molecule) {

  char *def_type;
  ATOMLIST *site_atoms;

  site_atoms = NULL;

  if (molecule->flags & PROTEIN_MOLECULE) {

    def_type = (params_get_parameter("site_def_type"))->value.s;
    
    if (!strcmp(def_type,"ligand")) {
      
      if ((molecule->system) && (molecule->system->site_ligand)) {
	
	site_atoms = molecule2shell(molecule,molecule->system->site_ligand,(params_get_parameter("site_max_dist"))->value.d);
      }
    }
  }

  return(site_atoms);
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


// TODO: atom edges should replace atom connections
// as the current connections are just atom lists
// and cannot contain information on the nature of the edge/bond

static void set_molecule_edges(MOLECULE *molecule) {

  int i;
  ATOM *atom;

  for (i=0,atom=molecule->atom;i<molecule->natoms;i++,atom++) {

    set_atom_edges(atom);
  }
}



static void set_molecule_hydrogens(MOLECULE *molecule) {

  int i,n_hydrogens;
  LIST *nonh_atoms;
  ATOM **atomp,*atom;

  if (molecule->hydrogens != NULL) {

    error_fn("%s: can only set hydrogens once for molecule %s",__func__,molecule->name);
  }

  // create list of non-hydrogen atoms:

  molecule->nonh_atoms = new_list("non-hydrogen atoms",sizeof(ATOM*),0);

  nonh_atoms = molecule->nonh_atoms;

  for (i=0,atom=molecule->atom;i<molecule->natoms;i++,atom++) {
 
    if (atom->element->id != HYDROGEN) {
  
      atomp = (ATOM**) add_list_item(nonh_atoms);

      *atomp = atom;
    }
  }

  // create lists of hydrogen atoms:

  n_hydrogens = molecule->natoms - nonh_atoms->n_items;

  if (n_hydrogens) {

    molecule->hydrogens = new_list("hydrogen atoms",sizeof(ATOM*),n_hydrogens);

    molecule->hydrogen_lists = new_list("hydrogen lists",sizeof(LIST),nonh_atoms->n_items);

    for (i=0,atomp=(ATOM**) nonh_atoms->items;i<nonh_atoms->n_items;i++,atomp++) {
      
      set_atom_hydrogens(*atomp);
    }
  }
}



static void set_atom_hydrogens(ATOM *atom) {

  int i;
  double Do,D;
  ATOM *hydrogen;
  MOLECULE *molecule;
  MAP *covalent_map;

  if (atom->hydrogens) {

    error_fn("%s: can only set atom hydrogens once",__func__);
  }

  molecule = atom->molecule;

  atom->hydrogens = (LIST*) add_list_item(molecule->hydrogen_lists);
  
  atom->hydrogens->items = get_list_item(molecule->hydrogens,molecule->hydrogens->n_items);
  
  atom->hydrogens->n_items = 0;
 
  if (molecule->rely_on_bonds) {

    set_atom_hydrogens_from_bonds(atom);

  } else {

    covalent_map = molecule->covalent_map;

    if (covalent_map) {

      set_atom_hydrogens_from_map(atom,covalent_map);

    } else {

      Do = atom->element->cov_radius + (get_element("H",NULL))->cov_radius;

      for (i=0,hydrogen=molecule->atom;i<molecule->natoms;i++,hydrogen++) {

	if (hydrogen->element->id == HYDROGEN) {

	  D = distance(atom->position,hydrogen->position);

	  if (D < Do + COVALENT_TOLERANCE) {

	    add_atom_hydrogen(atom,hydrogen);
	  }
	}
      }
    }
  }
}



static void set_atom_hydrogens_from_bonds(ATOM *atom) {

  int i;
  MOLECULE *molecule;
  ATOM *atom1,*atom2;
  BOND *bond;

  molecule = atom->molecule;

  for (i=0,bond=molecule->bond;i<molecule->nbonds;i++,bond++) {

    atom1 = bond->atom1;
    atom2 = bond->atom2;

    if ((atom1 == atom) && (atom2->element->id == HYDROGEN)) {

      add_atom_hydrogen(atom,atom2);

    } else if ((atom2 == atom) && (atom1->element->id == HYDROGEN)) {

      add_atom_hydrogen(atom,atom1);
    }
  }
}



static void set_atom_hydrogens_from_map(ATOM *atom,MAP *map) {

  int i,flag,iv[3],irange[3][2],ix,iy,iz;
  double Do,D;
  GRID *grid;
  ATOM **atomp,*hydrogen;
  ATOMLIST ***matrix,*list;

  grid = map->grid;
  matrix = (ATOMLIST***) map->matrix;

  flag = pos2grid(atom->position,iv,grid);

  if (flag != 0) {

    return;
  }

  for (i=0;i<3;i++) {

    irange[i][0] = (iv[i] == 0) ? 0 : iv[i] - 1;
    irange[i][1] = (iv[i]+1 == grid->npoints[i]) ? iv[i] : iv[i] + 1;
  }

  Do = atom->element->cov_radius + (get_element("H",NULL))->cov_radius;

  for (ix=irange[0][0];ix<=irange[0][1];ix++) {

    for (iy=irange[1][0];iy<=irange[1][1];iy++) {

      for (iz=irange[2][0];iz<=irange[2][1];iz++) {

	list = &(matrix[ix][iy][iz]);

	for (i=0,atomp=list->atom;i<list->natoms;i++,atomp++) {

	  hydrogen = *atomp;

	  if (hydrogen->element->id == HYDROGEN) {

	    D = distance(atom->position,hydrogen->position);
	    
	    if (D < Do + COVALENT_TOLERANCE) {
	      
	      add_atom_hydrogen(atom,hydrogen);
	    }
	  }
	}
      }
    }
  }
}



static void add_atom_hydrogen(ATOM *atom,ATOM *hydrogen) {

  int i;
  MOLECULE *molecule;
  ATOM **hp;

  molecule = atom->molecule;

  if (atom->hydrogens == NULL) {

    error_fn("%s: atom hydrogens need to be initialised for atom %d in molecule %s",__func__,atom->id,molecule->name);
  }

  for (i=0,hp=(ATOM**) atom->hydrogens->items;i<atom->hydrogens->n_items;i++,hp++) {

    if (*hp == hydrogen) {

      return;
    }
  }

  if (molecule->hydrogens == NULL) {

    error_fn("%s: molecule hydrogens need to be initialised for atom %d in molecule %s",__func__,atom->id,molecule->name);
  }

  hp = (ATOM**) add_list_item(molecule->hydrogens);

  *hp = hydrogen;

  atom->hydrogens->n_items++;
}



static void set_molecule_substructures(MOLECULE *molecule) {

  int i,j;
  ATOM **atom;
  LIST *substructure;

  if ((molecule->subatoms) || (molecule->substructures)) {

    error_fn("%s: cannot set molecule substructures twice",__func__);
  }

  molecule->subatoms = new_list("molecule substructure atoms",sizeof(ATOM*),molecule->nonh_atoms->n_items);
  molecule->substructures = new_list("molecule substructures",sizeof(LIST),molecule->nonh_atoms->n_items);

  for (i=0,atom=(ATOM**)molecule->nonh_atoms->items;i<molecule->nonh_atoms->n_items;i++,atom++) {

    set_atom_substructure(*atom);
  }

  // some tedious book-keeping:

  condense_list(molecule->substructures);

  for (i=0,substructure=(LIST*)molecule->substructures->items;i<molecule->substructures->n_items;i++,substructure++) {

    for (j=0,atom=(ATOM**)substructure->items;j<substructure->n_items;j++,atom++) {

      (*atom)->substructure = substructure;
    }
  }
}



static void set_atom_substructure(ATOM *atom1) {

  int i,j;
  MOLECULE *molecule;
  ATOM **atom2p,*atom2,**subatom;
  LIST *substructure,*subatoms;

  if (atom1->substructure) {

    return;
  }

  molecule = atom1->molecule;
  subatoms = molecule->subatoms;

  substructure = add_list_item(molecule->substructures);

  substructure->n_items = 0;

  substructure->items = get_list_item(subatoms,molecule->subatoms->n_items);

  for (i=0,atom2p=(ATOM**)molecule->nonh_atoms->items;i<molecule->nonh_atoms->n_items;i++,atom2p++) {

    atom2 = *atom2p;

    if (atom2->substructure == NULL) {

      if (same_residue_atoms(atom1,atom2)) {

	subatom = (ATOM**) add_list_item(subatoms);

	*subatom = atom2;

	atom2->substructure = substructure;

	substructure->n_items++;
      }
    }
  } 
}



LIST* atoms2substructures(LIST *atoms) {

  int i,j;
  ATOM **atom;
  LIST **substructure,*substructures,**item;

  substructures = new_list("substructures",sizeof(LIST*),0);

  for (i=0,atom=(ATOM**)atoms->items;i<atoms->n_items;i++,atom++) {

    substructure = &((*atom)->substructure);

    if (!item_in_list(substructures,(void*) substructure)) {

      item = (LIST**) add_list_item(substructures);

      *item = *substructure;
    }
  }

  return(substructures);
}



// TODO: these are currently calculated from the atom connections
// but in future it would be better to remove the atom connections
// and use the edges instead; then we will have to calculate the edges
// from scratch

void set_atom_edges(ATOM *atom) {

  int i;
  ATOM **batomp,*batom;
  BOND *bond;
  ATOM_EDGE *edge;
  ATOMLIST *conns;
  MOLECULE *molecule;

  conns = atom->connections;

  if (conns == NULL) {

    return;
  }

  if (conns->natoms > MAX_ATOM_EDGES) {

    atom->n_edges = MAX_ATOM_EDGES;

    warning_fn("%s: exceeded maxiumum number of edges (%d) for atom %d",__func__,MAX_ATOM_EDGES,atom->id);

  } else {

    atom->n_edges = conns->natoms;
  }

  molecule = atom->molecule;

  for (i=0,batomp=conns->atom,edge=atom->edges;i<atom->n_edges;i++,batomp++,edge++) {

    batom = *batomp;

    edge->atom = batom;
    edge->bond = get_bond(molecule,atom,batom);
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

	free_atomlist(atom->h_connections);

	atom->connections = NULL;

	atom->h_connections = NULL;

	atom->n_hydrogens = 0;

	atom->n_edges = 0;

	atom->hydrogens = NULL;
      }
    }
  }

  free_list(molecule->hydrogens);
  free_list(molecule->hydrogen_lists);

  molecule->hydrogens = NULL;
  molecule->hydrogen_lists = NULL;
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
    atom1->h_connections = (ATOMLIST*) alloc_atomlist();
    init_atomlist(atom1->h_connections);
  }

  if (atom1->element->id == HYDROGEN) {
    molecule2connections(atom1, atom1->molecule, NULL, NULL, 0, 0);
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
	// don't delete bonds to hydrogens as the connections are stored in hydrogen atom
	if ((batom != NULL) && (batom->element->id != HYDROGEN) && (!atom_in_list(atom1->connections,batom))) {

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
  ATOMLIST *connections, *h_connections;
  BOND *bond;

  connections = atom->connections;
  h_connections = atom->h_connections;

  for (i=0,bond=molecule->bond;i<molecule->nbonds;i++,bond++) {

    atom1 = bond->atom1;
    atom2 = bond->atom2;

    if (atom1 == atom) {

      if (atom2->element->id == HYDROGEN) {

	atom->n_hydrogens++;
	if (!atom_in_list(h_connections, atom2)){
	  add_atom_to_list(h_connections, atom2, bond->type);
	}

      } else if (!atom_in_list(connections,atom2)) {

	add_atom_to_list(connections,atom2,bond->type);
      }

    } else if (atom2 == atom) {

      if (atom1->element->id == HYDROGEN) {

	atom->n_hydrogens++;
	if (!atom_in_list(h_connections, atom1)){
          add_atom_to_list(h_connections, atom1, bond->type);
        }

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

      if (!atom_in_list(atom1->h_connections, atom2)){

        add_atom_to_list(atom1->h_connections, atom2, 1);
      }

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



int connected_atoms(ATOM *atom1,ATOM *atom2,ATOM *last_atom,int max_bonds) {

  int i;
  ATOM **batomp,*batom;
  ATOMLIST *list;

  if (max_bonds == 0) {

    return(0);
  }

  list = atom1->connections;

  if (list == NULL) {

    return(0);
  }

  for (i=0,batomp=list->atom;i<list->natoms;i++,batomp++) {

    batom = *batomp;

    if (batom == atom2) {

      return(1);
    }

    if ((batom != last_atom) && (connected_atoms(batom,atom2,atom1,max_bonds-1))) {

      return(1);
    }
  }

  return(0);
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


RESIDUE* atomlist2residues(ATOMLIST *list,int *n_residues) {

  int i,j,n_res,new_residue;
  ATOM **atomp,*atom;
  RESIDUE *residues,*residue;

  *n_residues = 0;

  if (list == NULL) {

    return(NULL);
  }

  residues = (RESIDUE*) calloc(list->natoms,sizeof(RESIDUE));

  if (residues == NULL) {

    error_fn("%s: out of memory allocating residues",__func__);
  }

  n_res = 0;

  for (i=0,atomp=list->atom;i<list->natoms;i++,atomp++) {

    atom = *atomp;

    new_residue = 1;

    for (j=0,residue=residues;j<n_res;j++,residue++) {

      if ((!strcmp(atom->subname,residue->name)) && (atom->subid == residue->id) &&
      	  (atom->chain == residue->chain) && (atom->icode == residue->icode)) {

      	new_residue = 0;
      }
    }

    if (new_residue) {

      residue = residues + n_res;

      residue->id = atom->subid;
      strcpy(residue->name,atom->subname);
      residue->chain = atom->chain;
      residue->icode = atom->icode;

      n_res++;
    }
  }

  residues = (RESIDUE*) realloc(residues,n_res*sizeof(RESIDUE));

  *n_residues = n_res;

  return(residues);
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

    return(TRUE);
  }

  return(FALSE);
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



void copy_atom_coords(ATOM *atom1,ATOM *atom2) {

  copy_vector(atom1->position,atom2->position);
  copy_vector(atom1->u,atom2->u);
  copy_vector(atom1->v,atom2->v);
  copy_vector(atom1->w,atom2->w);

  atom2->status |= ATOM_MOVED;
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

ATOM_COORDS* molecule2coords(MOLECULE *molecule) {
  int i;
  ATOM *atom;
  ATOM_COORDS *coords,*coord;

  if (molecule->natoms == 0) {
    return(NULL);
  }

  coords = (ATOM_COORDS*) calloc(molecule->natoms,sizeof(ATOM_COORDS));

  if (coords == NULL) {
    error_fn("get_list_coords: out of memory allocating coords");
  }

  for (i=0, atom=molecule->atom, coord=coords; i<molecule->natoms; i++,atom++,coord++) {
    copy_vector(atom->position,coord->position);
    copy_vector(atom->u,coord->u);
    copy_vector(atom->v,coord->v);
    copy_vector(atom->w,coord->w);
  }

  return(coords);
}



void coords2molecule(ATOM_COORDS *coords, MOLECULE *molecule) {
  int i, natoms_coords;
  ATOM *atom;
  ATOM_COORDS *coord;


  for (i=0, atom=molecule->atom, coord=coords; i<molecule->natoms; i++,atom++,coord++) {
    copy_vector(coord->position, atom->position);
    copy_vector(coord->u, atom->u);
    copy_vector(coord->v, atom->v);
    copy_vector(coord->w, atom->w);
  }
}




void coordinates2molecule(MOLECULE *molecule,double **v) {

  int i;
  ATOM *atom;

  for (i=0,atom=molecule->atom;i<molecule->natoms;i++,atom++) {

    move_atom(atom,v[i]);
  }
}



int linear_atom(ATOM *atom) {

  double v1[4],v2[4];
  ATOMLIST *conns;

  conns = atom->connections;

  if (conns == NULL) {

    return(UNDEFINED);
  }

  if (conns->natoms != 2) {

    return(UNDEFINED);
  }
  
  calc_vector(atom->position,(*(conns->atom))->position,v1);
  calc_vector(atom->position,(*(conns->atom+1))->position,v2);

  if (vector_angle(v1,v2) > 160.0) {

    return(TRUE);
  }

  return(FALSE);
}



int planar_atom(ATOM *atom) {

  int i,n_edges;
  double v[3][4],angle_sum;
  ATOM **hydrogen;
  ATOM_EDGE *edge;

  if (atom->planar != UNDEFINED) {

    return(atom->planar);
  }

  n_edges = (atom->molecule->use_hydrogens) ? atom->n_edges + atom->hydrogens->n_items : atom->n_edges;

  if (n_edges != 3) {

    return(UNDEFINED);
  }

  n_edges = 0;

  for (i=0,edge=(ATOM_EDGE*) atom->edges;i<atom->n_edges;i++,edge++,n_edges++) {

    calc_vector(atom->position,edge->atom->position,v[n_edges]);
  }

  if (atom->molecule->use_hydrogens) {

    for (i=0,hydrogen=(ATOM**) atom->hydrogens->items;i<atom->hydrogens->n_items;i++,hydrogen++,n_edges++) {
    
      calc_vector(atom->position,(*hydrogen)->position,v[n_edges]);
    }
  }

  angle_sum = vector_angle(v[0],v[1]) + vector_angle(v[0],v[2]) + vector_angle(v[1],v[2]);

  if (angle_sum > 355.0) {

    atom->planar = TRUE;

    return(TRUE);
  }
  
  atom->planar = FALSE;
  
  return(FALSE);
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

    c[i] /= (double) list->natoms;
  }
}

int heavy_atom_count(MOLECULE *molecule) {
  int i;
  int n = 0;
  ATOM *atom;
  for (i=0,atom=molecule->atom;i<molecule->natoms;i++,atom++) {
    if ((atom->element->id != HYDROGEN) && (!(atom->flags & SKIP_ATOM))) {
      n++;
    }
  }
  return(n);
}


ATOM_TYPE_LIST* molecule2atomtypelist(MOLECULE *molecule) {
  int i, j, present;
  ATOM *atom;

  ATOM_TYPE_LIST *list = malloc(sizeof(ATOM_TYPE_LIST));
  if (list == NULL) {
    error_fn("%s: error allocating ATOM_TYP_LIST\n", __func__);
  }
  list->n_alloc_types = 0;
  list->n_types = 0;

  for (i=0; i<molecule->natoms; i++) {
    atom = molecule->atom+i;
    present = 0;

    if (atom->element->id == HYDROGEN) {
      continue;
    }

    for (j=0; j<list->n_types; j++) {
      if (atom->type == list->types[j]) {
	present = 1;
	break;
      }
    }

    if (! present) {
      add_atom_type_to_list(atom->type, list);
    }
  }

  return(list);
}


static void add_atom_type_to_list(ATOM_TYPE *type, ATOM_TYPE_LIST *list) {
  if (list->n_alloc_types == 0) {
    list->n_alloc_types = 10;
    list->types = calloc(list->n_alloc_types, sizeof(ATOM_TYPE*));
  } else if (list->n_types == list->n_alloc_types) {
    list->n_alloc_types += 10;
    list->types = realloc(list->types,(list->n_alloc_types)*sizeof(ATOM_TYPE*));
  }

  if (list->types == NULL) {
    error_fn("%s: out of memory allocating atom_type_list");
  }

  list->types[list->n_types] = type;
  list->n_types++;
}


void set_molecule_min_coords(MOLECULE *molecule, double *min_coords) {
  double orig_min_coords[3] = {1e10, 1e10, 1e10};
  int i, j;
  ATOM *atom;

  for (i=0; i<molecule->natoms; i++) {
    atom = molecule->atom+i;
    for (j=0; j<3; j++) {
      if (atom->position[j] < orig_min_coords[j]) {
	orig_min_coords[j] = atom->position[j];
      }
    }
  }

  for (i=0; i<molecule->natoms; i++) {
    atom = molecule->atom+i;
    for (j=0; j<3; j++) {
      atom->position[j] = atom->position[j] - orig_min_coords[j] + min_coords[j];
    }
  }
}

double molecule_max_size(MOLECULE *molecule) {
  int i, j;
  ATOM *atom;
  double min_coords[3], max_coords[3];
  int extremes_set = 0;
  double max_size;

  for (i=0, atom=molecule->atom; i<molecule->natoms; i++,atom++) {
    if ((atom->element->id == HYDROGEN) || (atom->flags & SKIP_ATOM)) {
      continue;
    }

    for (j=0; j<3; j++) {
      if (! extremes_set) {
	min_coords[j] = max_coords[j] = atom->position[j];
      } else {
	if (atom->position[j] < min_coords[j]) {
	  min_coords[j] = atom->position[j];
	}
	if (atom->position[j] > max_coords[j]) {
	  max_coords[j] = atom->position[j];
	}
      }
    }
    extremes_set = 1;
  }

  if (! extremes_set) {
    error_fn("%s: can not get maximum size of the molecule\n", __func__);
  }
  max_size = max_coords[0] - min_coords[0];
  for (i=1; i<3; i++) {
    if ((max_coords[i] - min_coords[i]) > max_size) {
      max_size = max_coords[i] - min_coords[i];
    }
  }
  return max_size;
}

int molecule_accessible(MAP *field, MOLECULE *molecule) {
  int i;
  int gp[3];
  int index;
  ATOM *atom;

  for (i=0, atom=molecule->atom; i<molecule->natoms; i++,atom++) {
    if ((atom->element->id != HYDROGEN) && (!(atom->flags & SKIP_ATOM))) {
      if (pos2grid_round(atom->position, gp, field->grid)) {
        error_fn("%s: atom lies outside the grid\n", __func__);
      }
      index = field->grid->npoints[1] * field->grid->npoints[2] * gp[0] + \
              field->grid->npoints[2] * gp[1] + \
              gp[2];
      if (is_mask_bit_set(field->mask, index)) {
        return(0);
      }
    }
  }
  return(1);
}


void shake_molecule(MOLECULE *molecule) {
  int i, j;
  ATOM *atom;
  double shake;

  for (i=0, atom=molecule->atom; i<molecule->natoms; i++, atom++) {

    shake_atom(atom,MIN_COORD_SHAKE, MAX_COORD_SHAKE);
  }
}



void shake_atom(ATOM *atom,double min_shift,double max_shift) {

  int i;
  double shift;

  for (i=0;i<3;i++) {

    shift = uniform_rand(min_shift,max_shift);

    if (uniform_rand_int(0,1)) {

      shift *= -1.0;
    }

    atom->position[i] += shift;
  }
}



BOOLEAN **create_adjacency_matrix (MOLECULE *molecule)
{
  ATOMLIST *list;
  BOOLEAN **matrix;
  int i, j;
  ATOM *atom1, *atom2;

  list = molecule->selection;
  if (list == NULL) {
    error_fn("%s: atom selection list is not defined", __func__);
  }
  if (list->natoms == 0) {
    error_fn("%s: no atom is selected", __func__);
  }

  matrix = alloc_2d_bool_matrix(list->natoms, list->natoms);

  for (i=0; i<list->natoms; i++) {
    atom1 = list->atom[i];
    for (j=i; j<list->natoms; j++) {
      atom2 = list->atom[j];
      if (atom_in_list(atom1->connections, atom2)) {
        matrix [i][j] = matrix [j][i] = 1;
      } else {
        matrix [i][j] = matrix [j][i] = 0;
      }
    }
  }

  return matrix;
}


void get_selected_atom_indices(MOLECULE *molecule, int *ids) {
  int i;
  for (i=0; i<molecule->selection->natoms; i++) {
    ids[i] = molecule->selection->atom[i] - molecule->atom;
  }
}

void coords_bounds(ATOM_COORDS **coords, int n_coords, MOLECULE *molecule, double bounds[3][2]) {
  int i, j, k;
  int selected_indices[molecule->natoms];

  if (! molecule->selection) {
    error_fn("%s: no atom selected in the molecule", __func__);
  }

  get_selected_atom_indices(molecule, selected_indices);

  for (i=0; i<3; i++) {
    bounds[i][0] =  1.0E30;
    bounds[i][1] = -1.0E30;
  }

  for (i=0; i<n_coords; i++) {
    for (j=0; j<molecule->selection->natoms; j++) {
      for (k=0; k<3; k++) {
        if ((coords[i] + selected_indices[j])->position[k] < bounds[k][0]) {
          bounds[k][0] = (coords[i] + selected_indices[j])->position[k];
        }
        if ((coords[i] + selected_indices[j])->position[k] > bounds[k][1]) {
          bounds[k][1] = (coords[i] + selected_indices[j])->position[k];
        }
      }
    }
  }
}


GRID* coords2grid(ATOM_COORDS **coords, int n_coords, MOLECULE *molecule, double spacing, double padding, double *center) {

  int i, j;
  double bounds[3][2];
  GRID *grid;


  coords_bounds(coords, n_coords, molecule, bounds);

  grid = (GRID*) malloc(sizeof(GRID));


  if (grid == NULL) {

    error_fn("%s: out of memory allocating grid", __func__);
  }
  grid->padding = padding;
  grid->spacing = spacing;

  int n_points;
  int limit1, limit2;

  for (i=0; i<3; i++) {
    if (center[i] < bounds[i][0] || center[i] > bounds[i][1]) {
      error_fn("%s: center not within the bounds", __func__);
    }
    // symmtric grid around the center
    if ((center[i] - bounds[i][0]) >= (bounds[i][1] - center[i])) {
      n_points = (int) ceil((center[i] - bounds[i][0] + grid->padding)/grid->spacing);
    } else {
      n_points = (int) ceil((-center[i] + bounds[i][1] + grid->padding)/grid->spacing);
    }
    grid->npoints[i] = n_points * 2 + 1;
    grid->limit[i][0] = (int) round(center[i] - (double) n_points);
    grid->limit[i][1] = grid->limit[i][0] + grid->npoints[i] -1;

    grid->flimit[i][0] = ((double) grid->limit[i][0])*(grid->spacing);
    grid->flimit[i][1] = ((double) grid->limit[i][1])*(grid->spacing);
  }

  return(grid);
}

ATOMLIST* get_hbond_atomlist(MOLECULE *molecule){
  int i;

  ATOMLIST *hbond_atomlist;
  hbond_atomlist = alloc_atomlist();
  init_atomlist(hbond_atomlist);

  if (! molecule->selection) {
    error_fn("%s: selection not defined", __func__);
  }

  if (! molecule->selection->natoms) {
    error_fn("%s: no atom selected");
  }


  for (i=0; i<molecule->selection->natoms; i++) {
    if (molecule->selection->atom[i]->type->flags & HBOND_DA_ATOM_TYPE) {
      add_atom_to_list(hbond_atomlist, molecule->selection->atom[i], molecule->selection->atom[i]->flags);
    }
  }
  return hbond_atomlist;
}

ATOMLIST* molecule2hbonding_atoms(MOLECULE *molecule) {
  ATOMLIST *list = alloc_atomlist();
  init_atomlist(list);

  int i;
  ATOM *atom;
  for (i=0, atom=molecule->atom; i<molecule->natoms; i++, atom++) {
    if ((atom->flags & SKIP_ATOM) || (atom->element->id == HYDROGEN)) {
      continue;
    }
    if (((atom->type->flags & HBOND_DONOR_ATOM_TYPE) == HBOND_DONOR_ATOM_TYPE) || ((atom->type->flags & HBOND_ACCEPTOR_ATOM_TYPE) == HBOND_ACCEPTOR_ATOM_TYPE)) {
      add_atom_to_list(list, atom, atom->flags);
    }
  }
  return (list);
}

int has_hbonding_atom(MOLECULE *molecule) {
  int i;
  ATOM *atom;
  for (i=0, atom=molecule->atom; i<molecule->natoms; i++, atom++) {
    if ((atom->flags & SKIP_ATOM) || (atom->element->id == HYDROGEN)) {
      continue;
    }
    if (((atom->type->flags & HBOND_DONOR_ATOM_TYPE) == HBOND_DONOR_ATOM_TYPE) || ((atom->type->flags & HBOND_ACCEPTOR_ATOM_TYPE) == HBOND_ACCEPTOR_ATOM_TYPE)) {
      return(1);
    }
  }
  return (0);
}

MOLECULE* create_atom_type_molecules(ATOM_TYPE *type, int n_molecule, char *molecule_name, SETTINGS *settings) {
  MOLECULE *molecules = calloc(n_molecule, sizeof(MOLECULE));

  if (molecules == NULL) {
    error_fn("realloc_probe_atoms: out of memory allocating probe molecule");
  }

  int i;
  MOLECULE *molecule;
  for (i=0, molecule = molecules; i<n_molecule; i++, molecule++) {
    init_molecule(molecule,1);
    ATOM *atom = molecule->atom + molecule->natoms;
    init_atom(atom);
    null_vector(atom->position);
    atom->type = type;
    atom->node = type->united_atom->atom_node;
    atom->element = (atom->type->id == -1) ? atom->type->element : atom->type->united_atom->atom_node->element;
    atom->hbond_geometry = type->hbond_geometry;
    atom->molecule = molecule;
    atom->geometry = type->geometry;
    //atom->flags |= SINGLE_ATOM_PROBE;

    molecule->natoms++;
    strcpy(molecule->name, molecule_name);
    molecule->selection = molecule2atoms(molecule);
    molecule->flags |= LIGAND_MOLECULE;
    molecule->molsystem = molecule2system(molecule, settings);
    set_molecule_vdw_radii(molecule,settings->water_vdw_radius);
  }
  return(molecules);
}

MOLECULE_PROBE* get_molecule_probe(char *probe_name, SETTINGS *settings, char *skip_atoms_csv) {

  INT_LIST *skip_atoms = NULL;
  if ((skip_atoms_csv) && (strcmp(skip_atoms_csv, "undefined"))) {
    skip_atoms = parse_int_csv(skip_atoms_csv, ",");
    if (skip_atoms->n_items == 0) {
      error_fn("%s: No atom ids could be parsed from \"%s\"", __func__, skip_atoms_csv);
    }
  }
  if (is_readable_file(probe_name)){
    MOLECULE_PROBE *probe = create_mol_probe_from_file(probe_name, settings, skip_atoms);
    return probe;
  }

  // look for the probe_name in the probes directory
  char *pli_dir = get_pli_dir();
  char file_name[MAX_LINE_LEN];
  sprintf(file_name, "%s/params/probes/%s.sdf", pli_dir, probe_name);
  if (is_readable_file(file_name)){
    MOLECULE_PROBE *probe = create_mol_probe_from_file(file_name, settings, skip_atoms);
    return probe;
  }

  error_fn("%s: Error initializing molecule probe %s", __func__, probe_name);
}


static MOLECULE_PROBE* create_mol_probe_from_file(char *filename, SETTINGS *settings, INT_LIST *skip_aids) {

  int i;
  ATOM *atom;
  MOLECULE_PROBE *probe;
  MOLECULE *molecule;

  probe = alloc_molecule_probe();

  molecule = read_molecule(filename);
  molecule->flags |= LIGAND_MOLECULE;

  // to resolve the tautomers
  molecule->use_hydrogens = 1;

  molecule->selection = molecule2atoms(molecule);
  prep_molecule(molecule, settings);

  //
  shake_molecule(molecule);

  if (skip_aids) {
    skip_atom_ids(molecule, skip_aids);
  }

  // remove flagged atoms from the selection
  for (i=0, atom=molecule->atom; i<molecule->natoms; i++, atom++) {
    if ((atom->flags & SKIP_ATOM) && (atom_in_list(molecule->selection, atom))) {
      remove_atom_from_list(molecule->selection, atom);
    }
  }

  if (heavy_atom_count(molecule) == 1) {
    error_fn("%s: molecular probe should have at least two atoms, only one atom found in \"%s\". Consider running mode \"field\"", __func__, filename);
  }

  molecule->molsystem = molecule2system(molecule,settings);
  probe->molecule = molecule;
  probe->selected_ids = calloc(molecule->selection->natoms, sizeof(int));
  get_selected_atom_indices(molecule, probe->selected_ids);
  return(probe);
}

MOLECULE_PROBE* alloc_molecule_probe(void) {
  MOLECULE_PROBE *probe;

  probe = (MOLECULE_PROBE*) malloc(sizeof(MOLECULE_PROBE));
  if (probe == NULL) {
    error_fn("%s: out of memory allocating molecule probe", __func__);
  }

  probe->atom_maps = NULL;
  probe->molecule = NULL;
  probe->n_orientations = 0;
  probe->orientations = NULL;
  probe->score = 0.0;
  probe->orientation_id = -1;
  probe->auto_mappings = NULL;
  probe->n_alloc_mappings;
  probe->n_mappings = 0;
  return(probe);
}


double sqr_distance_offset(double *v1,double *v2, double *offset1, double *offset2) {

  return(sqr((v2[0]+offset2[0])-(v1[0]+offset1[0])) + sqr((v2[1]+offset2[1])-(v1[1]+offset1[1])) + sqr((v2[2]+offset2[2])-(v1[2]+offset1[2])));
}

double rms_coords_atleat(ATOM_COORDS *coords1, ATOM_COORDS *coords2, int *selected1, int *selected2, int n_selected, BOOLEAN ***isomorphs, int n_isomprphs, double upper_bound) {

  double min_rms, temp, upper_bound_sqr;
  double rms;


  min_rms = 1E10;
  int i,j,k;
  int coords_ind1, coords_ind2;


  upper_bound_sqr = upper_bound * upper_bound * n_selected;

  for (i=0; i<n_isomprphs; i++) {
    rms = 0.0;
    for (j=0; j<n_selected; j++) {
      coords_ind1 = (selected1 == NULL) ? j: selected1[j];
      for (k=0; k<n_selected; k++) {
        if (isomorphs[i][j][k]) {
          coords_ind2 = (selected2 == NULL) ? k: selected2[k];
          rms+= sqr_distance((coords1+coords_ind1)->position, (coords2+coords_ind2)->position);
          if (rms > upper_bound_sqr) {
            j = n_selected;
            break;
          }
        }
      }
    }
    rms = sqrt(rms/(double)n_selected);
    // early termination if close to zero rms is alreay found
    if (rms < 1E-10) {
      return(rms);
    }
    if (rms < min_rms) {
      min_rms = rms;
    }
  }
  return (min_rms);
}

double rms_coords(ATOM_COORDS *coords1, ATOM_COORDS *coords2, int *selected1, int *selected2, int n_selected, BOOLEAN ***isomorphs, int n_isomprphs, double *offset1, double *offset2) {

  double min_rms;
  double rms;


  min_rms = 1E100;
  int i,j,k;
  int coords_ind1, coords_ind2;



  for (i=0; i<n_isomprphs; i++) {
    rms = 0.0;
    for (j=0; j<n_selected; j++) {
      coords_ind1 = (selected1 == NULL) ? j: selected1[j];
      for (k=0; k<n_selected; k++) {
        if (isomorphs[i][j][k]) {
          coords_ind2 = (selected2 == NULL) ? k: selected2[k];
          if ((!offset1) || (!offset2)) {
            rms += sqr_distance((coords1+coords_ind1)->position, (coords2+coords_ind2)->position);
          } else {
            rms += sqr_distance_offset((coords1+coords_ind1)->position, (coords2+coords_ind2)->position, offset1, offset2);
          }
        }
      }
    }
    rms = sqrt(rms/(double)n_selected);
    // early termination if close to zero rms is alreay found
    if (rms < 1E-10) {
      return(rms);
    }
    if (rms < min_rms) {
      min_rms = rms;
    }
  }
  return (min_rms);
}
