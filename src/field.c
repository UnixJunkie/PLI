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



static struct FieldSettings {
  int n_alloc_probes;
  int n_probes;
  FIELD_PROBE *probes;
  FIELD_PROBE *probe;
  char probe_name[MAX_LINE_LEN];
  double grid_spacing;
  double grid_padding;
  char grid_file[MAX_LINE_LEN];
  int optimise_u_axis;
  int n_lebedev_points;
} *field_settings = NULL;



static void add_atom_type_probe(ATOM_TYPE*,SETTINGS*);
static void alloc_field_probes(void);
static void init_field_probe(FIELD_PROBE*);
static void realloc_probe_atoms(FIELD_PROBE*);
static double calc_field_score(FIELD_PROBE*,MAP*,SYSTEM*,PLI_SFUNC*);
static GRID* field_settings2grid(SYSTEM*);
static FIELD_PROBE* alloc_field_probe(void);
static double calc_field_score(FIELD_PROBE*,MAP*,SYSTEM*,PLI_SFUNC*);
static MAP* new_field(char*,SYSTEM*);



void run_field(SETTINGS *settings) {

  SYSTEM *system;
  FIELD_PROBE *probe;
  MAP *field;

  system = settings2system(settings);

  prep_system(system,ANY_MOLECULE);

  probe = get_field_probe(field_settings->probe_name);

  if (probe == NULL) {

    error_fn("unknown probe '%s'",field_settings->probe_name);
  }

  field = calculate_field(system,probe,NULL, settings->sfunc);

  write_insight_map(field_settings->grid_file,field,0);
  //write_insight_map("mask.grd",field,1);

  free_map(field);

  output_system(PLI_STDOUT,system);

  unprep_system(system,ANY_MOLECULE);

  free_system(system);
}



MAP* calculate_field(SYSTEM *system,FIELD_PROBE *probe, MAP *field, PLI_SFUNC *sfunc) {

  int ix,iy,iz,nx,ny,nz,index,score_inter,score_intra;
  double pos[3],***matrix,score,bad_score,radius;
  GRID *grid;
  ATOMLIST *selection;
  BYTE *volume_mask,*clash_mask;

  // only calculate inter-molecular terms:

  score_inter = (params_get_parameter("score_inter"))->value.i;
  score_intra = (params_get_parameter("score_intra"))->value.i;

  (params_get_parameter("score_inter"))->value.i = 1;
  (params_get_parameter("score_intra"))->value.i = 0;

  reset_system(system,0);

  if (field == NULL) {

    field = new_field(probe->name,system);
  }

  grid = field->grid;

  radius = 2.0*(probe->atom->vdw_radius_H2O) - (probe->atom->vdw_radius);

  volume_mask = mask_system_volume(system,NULL,grid,PROBE_MASK,radius);
  clash_mask = mask_system_atoms(system,NULL,grid,CONTACT_MASK,0.0);

  selection = system->selection;

  init_field_system(system,probe);

  prep_score_system(system,sfunc);

  // TODO: To account for edge cases
  params_get_parameter("grid_padding")->value.d += params_get_parameter("grid_spacing")->value.d;

  matrix = (double***) field->matrix;

  nx = grid->npoints[0];
  ny = grid->npoints[1];
  nz = grid->npoints[2];

  bad_score = sfunc->bad_score;

  for (ix=0;ix<nx;ix++) {

    pos[0] = grid->flimit[0][0] + ((double) ix)*(grid->spacing);

    for (iy=0;iy<ny;iy++) {

     pos[1] = grid->flimit[1][0] + ((double) iy)*(grid->spacing);

      for (iz=0;iz<nz;iz++) {

	index = ny*nz*ix + nz*iy + iz;

	if (is_mask_bit_set(clash_mask,index)) {

	  matrix[ix][iy][iz] = bad_score;

	} else if (is_mask_bit_set(volume_mask,index)) {

	  pos[2] = grid->flimit[2][0] + ((double) iz)*(grid->spacing);

	  move_probe_to_new_position(probe,pos);

	  score = calc_field_score(probe,field,system,sfunc);

	  matrix[ix][iy][iz] = score;
	}
      }
    }
  }

  free(volume_mask);
  free(clash_mask);

  finish_field_system(system,probe);

  system->selection = selection;

  reset_system(system,0);

  (params_get_parameter("score_inter"))->value.i = score_inter;
  (params_get_parameter("score_intra"))->value.i = score_intra;

  unprep_score_system(system,ANY_MOLECULE,sfunc);

  return(field);
}



FIELD_PROBE* get_field_probe(char *name) {

  int i;
  FIELD_PROBE *probe;

  for (i=0,probe=field_settings->probes;i<field_settings->n_probes;i++,probe++) {

    if (!strcmp(probe->name,name)) {

      return(probe);
    }
  }

  return(NULL);
}



void init_field_settings(void) {

  char *pli_dir;

  if (field_settings == NULL) {

    field_settings = (struct FieldSettings*) malloc(sizeof(struct FieldSettings));
  }

  if (field_settings == NULL) {

    error_fn("init_field_settings: out of memory allocating settings");
  }

  pli_dir = get_pli_dir();

  field_settings->n_alloc_probes = 0;
  field_settings->n_probes = 0;

  field_settings->probes = NULL;
  field_settings->probe = NULL;

  // grid settings:

  field_settings->grid_spacing = (params_get_parameter("grid_spacing"))->value.d;
  field_settings->grid_padding = (params_get_parameter("grid_padding"))->value.d;

  strcpy(field_settings->grid_file,(params_get_parameter("grid_file"))->value.s);

  // field probe:

  strcpy(field_settings->probe_name,(params_get_parameter("probe_name"))->value.s);

  // rotational parameters:

  field_settings->optimise_u_axis = (params_get_parameter("optimise_u_axis"))->value.i;

  field_settings->n_lebedev_points = (params_get_parameter("n_lebedev_points"))->value.i;
}



void init_field_probes(SETTINGS *settings) {

  int i;
  ATOM_TYPE *type;
  ATOM_TYPING_SCHEME *scheme;

  scheme = settings->atom_typing_scheme;

  for (i=0,type=scheme->atom_types;i<scheme->n_atom_types;i++,type++) {

    add_atom_type_probe(type,settings);
  }
}



static void add_atom_type_probe(ATOM_TYPE *type,SETTINGS *settings) {

  ATOM *atom;
  MOLECULE *molecule;
  FIELD_PROBE *probe;

  alloc_field_probes();

  probe = field_settings->probes + (field_settings->n_probes);

  init_field_probe(probe);

  strcpy(probe->name,type->name);

  realloc_probe_atoms(probe);

  probe->type = ATOM_PROBE;

  probe->atom = probe->molecule->atom;

  atom = probe->atom;

  molecule = probe->molecule;

  null_vector(atom->position);

  atom->type = type;

  //
  //atom->node = type->united_atom->atom_node;

  atom->element = (atom->type->id == -1) ? atom->type->element : atom->type->united_atom->atom_node->element;

  //
  atom->hbond_geometry = type->hbond_geometry;

  atom->molecule = molecule;

  atom->flags |= SINGLE_ATOM_PROBE;

  molecule->natoms++;

  strcpy(molecule->name,probe->name);

  molecule->selection = molecule2atoms(molecule);

  molecule->flags |= LIGAND_MOLECULE;

  molecule->molsystem = molecule2system(molecule,settings);

  set_molecule_vdw_radii(molecule,settings->water_vdw_radius);

  molecule->active_atoms = molecule_active_atoms(molecule);

  field_settings->n_probes++;
}


// TODO: use this for adding atom type probe. Or actually just create one probe, dont add all
FIELD_PROBE* create_atom_type_probe(ATOM_TYPE *type, SETTINGS *settings, FIELD_PROBE *probe_in) {
  // TODO: use create_atom_type_molecule function
  ATOM *atom;
  MOLECULE *molecule;
  FIELD_PROBE *probe;

  if (probe_in == NULL) {
    probe = alloc_field_probe();
    init_field_probe(probe);
  } else {
    probe = probe_in;
  }

  probe->type = ATOM_PROBE;

  strcpy(probe->name,type->name);

  realloc_probe_atoms(probe);

  probe->atom = probe->molecule->atom;

  atom = probe->atom;

  molecule = probe->molecule;

  null_vector(atom->position);

  atom->type = type;

  //printf("type: %s\n", type->name);

  atom->node = type->united_atom->atom_node;

  atom->element = (atom->type->id == -1) ? atom->type->element : atom->type->united_atom->atom_node->element;

  atom->hbond_geometry = type->hbond_geometry;

  atom->molecule = molecule;

  atom->flags |= SINGLE_ATOM_PROBE;

  molecule->natoms++;

  strcpy(molecule->name,probe->name);

  molecule->selection = molecule2atoms(molecule);

  molecule->active_atoms = molecule_active_atoms(molecule);

  molecule->flags |= LIGAND_MOLECULE;

  molecule->molsystem = molecule2system(molecule,settings);

  set_molecule_vdw_radii(molecule,settings->water_vdw_radius);

  return(probe);
}



static void alloc_field_probes(void) {

  if (field_settings->probes == NULL) {

    field_settings->n_alloc_probes = 10;

    field_settings->probes = (FIELD_PROBE*) calloc(field_settings->n_alloc_probes,sizeof(FIELD_PROBE));

  } else if (field_settings->n_probes == field_settings->n_alloc_probes) {

    field_settings->n_alloc_probes += 10;

    field_settings->probes = (FIELD_PROBE*) realloc(field_settings->probes,(field_settings->n_alloc_probes)*sizeof(FIELD_PROBE));
  }

  if (field_settings->probes == NULL) {
    
    error_fn("alloc_field_probes: out of memory (re)allocating probes");
  }
}



static void init_field_probe(FIELD_PROBE *probe) {

  probe->id = -1;

  strcpy(probe->name,"");

  probe->atom = NULL;
  probe->molecule = NULL;
}



static void realloc_probe_atoms(FIELD_PROBE *probe) {

  MOLECULE *molecule;

  if (probe->molecule == NULL) {

    probe->molecule = (MOLECULE*) malloc(sizeof(MOLECULE));

    if (probe->molecule == NULL) {

      error_fn("realloc_probe_atoms: out of memory allocating probe molecule");
    }

    init_molecule(probe->molecule,1);
  }

  molecule = probe->molecule;

  realloc_atoms(molecule);
  
  init_atom(molecule->atom + (molecule->natoms));
}



static MAP* new_field(char *probe_name,SYSTEM *system) {

  char field_name[MAX_LINE_LEN];
  MAP *field;
  GRID *grid;

  (probe_name != NULL) ? sprintf(field_name,"PLI field for %s",probe_name) : sprintf(field_name,"PLI field");

  grid = system_grid(system,(params_get_parameter("grid_spacing"))->value.d,(params_get_parameter("grid_padding"))->value.d);

  field = new_map(field_name,"double",grid);

  return(field);
}



void init_field_system(SYSTEM *system,FIELD_PROBE *probe) {

  ATOM *atom;
  MOLECULE *molecule;
  FORCE_FIELD *ff;

  ff = system->settings->force_field;

  molecule = probe->molecule;

  atom = probe->atom;
  atom->geometry = get_atom_geometry("X",system->settings->atom_typing_scheme);

  if (field_settings->optimise_u_axis) {

    atom->geometry = get_atom_geometry("-X",system->settings->atom_typing_scheme);

    probe->u_axes = read_lebedev_sphere_points(field_settings->n_lebedev_points);
  }
    
  add_molecule_to_list(system->molecule_list,probe->molecule);

  system->selection = probe->molecule->selection;
}



void finish_field_system(SYSTEM *system,FIELD_PROBE *probe) {

  if (field_settings->optimise_u_axis) {

    free_2d_fmatrix(probe->u_axes);
  }

  remove_molecule_from_list(system->molecule_list,probe->molecule);
}



void move_probe_to_new_position(FIELD_PROBE* probe,double *position) {

  if (probe->type == ATOM_PROBE) {

    move_atom(probe->atom,position);

  } else {

    error_fn("move_probe_to_new_position: unknown probe type");
  }
}



double calc_field_score(FIELD_PROBE *probe,MAP *field,SYSTEM *system,PLI_SFUNC *sfunc) {

  int i;
  double score,**u_axis;
  ATOM *atom;

  if (probe->type == ATOM_PROBE) {

    if (field_settings->optimise_u_axis) {

      atom = probe->atom;

      score = 999.9;

      for (i=0,u_axis=probe->u_axes;i<field_settings->n_lebedev_points;i++,u_axis++) {

	memcpy(atom->u,*u_axis,3*sizeof(double));
	
	set_atom_virtual_points(atom);

	atom->status |= ATOM_ROTATED;

	atom->cstatus &= ~ATOM_GEOMETRIES_CALCULATED;
	atom->cstatus &= ~ATOM_SCORES_CALCULATED;

	score_system(system,sfunc);

	if (system->score < score) {

	  score = system->score;
	}
      }

    } else {

      score_system(system,sfunc);

      score = system->score;
    }

  } else {

    error_fn("calc_field_score: unknown probe type");
  }

  return(score);
}



static GRID* field_settings2grid(SYSTEM *system) {

  GRID *grid = NULL;

  grid = atomlist2grid(system->selection,field_settings->grid_spacing,field_settings->grid_padding,grid);

  if (grid == NULL) {

    error_fn("settings2grid: grid undefined");
  }

  return(grid);
}

static FIELD_PROBE* alloc_field_probe(void) {
  FIELD_PROBE *probe;

  probe = (FIELD_PROBE*) malloc(sizeof(FIELD_PROBE));
  if (probe == NULL) {
    error_fn("alloc_field_probe: out of memory allocating probe");
  }
  return(probe);
}

