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

#include <libgen.h>

static void init_fragmap_settings(SYSTEM*);
static int get_optimal_n_rot_lebedev(double **lebedev_points, int n_points);
static double calculate_score(MOLECULE_PROBE *probe, double pos[4], SYSTEM *system, MAP *buffer_field);
static void update_molecule_score_matrix(MOLECULE_PROBE*, MAP*, MAP*, SYSTEM*);
static double score_probe_using_fields(MOLECULE_PROBE*, SETTINGS*);
static void write_probe_orientations (MOLECULE*, ATOM_COORDS**, int, char*);
static GRID* copy_grid_extend(GRID *to_grid, GRID *from_grid, double extend_by);
void calculate_atom_fields_for_probe(MOLECULE_PROBE *probe, SETTINGS *settings, SYSTEM *system);
void copy_mask_from_field(unsigned char* mask_to, MAP *field_from);
static void mask_ref_ligand(MAP *field, MOLECULE *ref_ligand, double max_dist, int mask_islands);
static MAP_ISLAND* get_ligand_island(MAP_ISLAND*, int, MAP*, MOLECULE*);
static void output_fragmap(MAP *field, MAP *poselist_map);
static void remove_redundant_orientations(ATOM_COORDS **probe_orientations, int *n_orientations, MOLECULE *molecule, BOOLEAN ***mappings, int n_mappings, double rms_cutoff);
static MAP* fragmap(SYSTEM *system, MOLECULE_PROBE *probe);
static void alloc_fragmap_data();
static void add_automorph(BOOLEAN**);
MAP* get_atom_field(MAP *map_in, MOLECULE *molecule, FIELD_PROBE *, int);
double solvate_pose(SYSTEM *system, int ixc, int iyc, int izc, DOCKED_WATERS*);

# define MAX_AUTOMORPHS 100
# define SOLVATE_LIG_MAX_SCORE -0.15
# define SOLVATE_PROTEIN_MAX_SCORE -0.15
# define WATER_PENALTY -0.3


double COORDS_ORIGIN[3] = {0.0, 0.0, 0.0};

// some variables for analyzing the speed of the code
static int N_SCORE_GRID, N_SCORE_FULL, N_SAVE;
static void init_debug();
//

static struct FragmapSettings {
  MOLECULE_PROBE *probe;
  INT_LIST *skip_atoms;
  char probe_name[MAX_LINE_LEN];
  char probe_file[MAX_LINE_LEN];
  double grid_spacing;
  double grid_padding;
  char grid_file[MAX_LINE_LEN];
  int n_lebedev_points;
  int precalculate;
  int save_intermediate_grids;
  int n_rot_lebedev_axis;
  int rescore;
  double **lebedev_points;
  double rescore_cutoff;
  double save_cutoff;
  int scan_dihedral;
  int n_dihedral_steps;
  MOLECULE *ref_ligand;
  double ref_ligand_max_dist;
  CLUSTER_SETTINGS *cluster_settings;
  int solvate;
  FIELD_PROBE *water_probe;
  int n_max_waters;
  double ligand_grid_padding;
  DOCKED_WATERS *docked_waters;
  MAP *ligand_water_field;
  MAP *ligand_water_field_backup;
  MAP *protein_water_field;
  MOLECULE *mask_mol;
} *fragmap_settings = NULL;


static struct FragmapData {
  MAP *field;
  MAP *buffer_field;
  POSE_LIST *saved_pose_list;
  PROBE_POSE *saved_allocated_poses;
  POSE_LIST *orientation_pose_list;
  PROBE_POSE *orientation_allocated_poses;
  POSE_LIST *gridpoint_pose_list;
  PROBE_POSE *gridpoint_allocated_poses;
  CLUSTER_SET *div_cluster_set;
  MAP *poselist_map;
  CLUSTER_SET *peak_cluster_set;
  PROBE_POSE *peak_allocated_poses;
} *fragmap_data = NULL;


static void init_debug() {
  N_SCORE_FULL = 0;
  N_SCORE_GRID = 0;
  N_SAVE = 0;
}

void run_fragmap(SETTINGS *settings) {
  init_debug();

  int i;
  SYSTEM *system;
  MOLECULE_PROBE *probe;
  MAP *poselist_map;

  // get system and prepare
  system = settings2system(settings);
  prep_system(system,ANY_MOLECULE);
  prep_score_system(system,system->settings->sfunc);


  // init fragment settings and probe
  init_fragmap_settings(system);
  probe = fragmap_settings->probe;

  // calculate atom field if required
  if (fragmap_settings->precalculate) {
    calculate_atom_fields_for_probe(probe, settings, system);
  }

  // allocate probe poses for saving
  alloc_fragmap_data();

  // calculate field
  fragmap_data->field = fragmap(system, probe);

  // TODO:
  // free_saved_poses_lists();
  // free_map(field);
  // free mol probe
  // free fragmap_settings
  // free orientations
  // free atom type grids
  output_fragmap(fragmap_data->field, fragmap_data->poselist_map);
  printf("N_SCORE GRID: %d\nN_SCORE_FULL: %d\nN_SAVE: %d\n", N_SCORE_GRID, N_SCORE_FULL, N_SAVE);

  output_system(PLI_STDOUT,system);

  unprep_score_system(system, ANY_MOLECULE,system->settings->sfunc);
  unprep_system(system, ANY_MOLECULE);

  free_system(system);
}


static MAP* fragmap(SYSTEM *system, MOLECULE_PROBE *probe) {

  ATOMLIST *selection = system->selection;
  reset_system(system,0);

  // fragmap_settings
  CLUSTER_SETTINGS *cluster_settings = fragmap_settings->cluster_settings;
  int n_orientations = probe->n_orientations;
  ATOM_COORDS **probe_orientations = probe->orientations;
  //write_probe_orientations(probe->molecule, probe_orientations, n_orientations, "probe_ori");

  // get map
  GRID *grid = system_grid(system,(params_get_parameter("grid_spacing"))->value.d,(params_get_parameter("grid_padding"))->value.d);
  MAP *field = new_map("Fragment map","double",grid);

  mask_molecule(system->protein,field->mask,field->grid,CONTACT_MASK);
  if (fragmap_settings->ref_ligand) {
    mask_ref_ligand(field, fragmap_settings->ref_ligand, fragmap_settings->ref_ligand_max_dist, 0);
  }
  if (fragmap_settings->mask_mol) {
    mask_molecule(fragmap_settings->mask_mol,field->mask,field->grid,ATOM_MASK);
  }


  int nx = field->grid->npoints[0];
  int ny = field->grid->npoints[1];
  int nz = field->grid->npoints[2];

  unsigned char *mask = field->mask;
  double ***matrix = (double***) field->matrix;


  // initialize the matrix
  int ix, iy, iz;
  for (ix=0;ix<nx;ix++) {
    for (iy=0;iy<ny;iy++) {
      for (iz=0;iz<nz;iz++) {
        matrix[ix][iy][iz] = system->settings->sfunc->bad_score * (system->settings->sfunc->order == DESCENDING) ? -1.0: 1.0;
      }
    }
  }


  // set grid for assessing molecule accessibility; it should be larger than the field grid
  // because molecule center is placed on each grid point and atoms will lie outside
  // the field grid as well

  // A better way is to use another buffer in the grid settings to be used only for mask
  // and for calculation and writing out grid, this extra buffer is not considered

  double bounds[3][2];
  MOLECULE *molecule = probe->molecule;
  coords_bounds(probe->orientations, probe->n_orientations, molecule, bounds);
  double molecule_size = 0.0;
  for (int i=0; i<3; i++) {
    if ((bounds[i][1]-bounds[i][0]) > molecule_size) {
      molecule_size = bounds[i][1]-bounds[i][0];
    }
  }


  MAP *buffer_field = new_map("buffer","null",NULL);

  buffer_field->grid = copy_grid_extend(buffer_field->grid, field->grid, molecule_size + fragmap_settings->ligand_grid_padding);
  alloc_map_matrix(buffer_field);
  mask_molecule(system->protein,buffer_field->mask,buffer_field->grid,CONTACT_MASK);
  if (system->ligand) {
    mask_molecule(system->ligand,buffer_field->mask,buffer_field->grid,CONTACT_MASK);
  }
  if (fragmap_settings->mask_mol) {
    mask_molecule(fragmap_settings->mask_mol,buffer_field->mask,buffer_field->grid,ATOM_MASK);
  }

  fragmap_data->buffer_field = buffer_field;

  // protein water field
  if (fragmap_settings->solvate) {
    fragmap_settings->protein_water_field = new_map("protein water map","double",NULL);
    fragmap_settings->protein_water_field->grid = copy_grid_extend(fragmap_settings->protein_water_field->grid, buffer_field->grid, 0.0);
    alloc_map_matrix(fragmap_settings->protein_water_field);
    fragmap_settings->protein_water_field = get_atom_field(fragmap_settings->protein_water_field, system->protein, fragmap_settings->water_probe, 0);
    //write_insight_map("protein_water.grd", fragmap_settings->protein_water_field, 0);
    //write_insight_map("protein_water_mask.grd", fragmap_settings->protein_water_field, 1);
  }

  // poselist map saving all poses on each grid point
  if (cluster_settings->save_peak_clusters) {

    fragmap_data->poselist_map = new_map("pose list map","null",copy_grid(NULL,grid));

    fragmap_data->poselist_map->matrix = (void***) alloc_3d_poselist_matrix(nx, ny, nz);
    init_3d_poselist_matrix(fragmap_data->poselist_map);
  }

  // add probe to the system
  add_molecule_to_list(system->molecule_list, molecule);
  // this is needed for calculating internal energy: may we can calculate energy once for each conformation and use that instead!
  system->ligand = molecule;
  system->selection = molecule->selection;


  // save n orientations per grid point
  int index;
  double pos[4];
  double score;
  if (cluster_settings->n_poses_per_gridpoint) {
    for (ix=0;ix<nx;ix++) {
      pos[0] = grid->flimit[0][0] + ((double) ix)*(grid->spacing);
      for (iy=0;iy<ny;iy++) {
        pos[1] = grid->flimit[1][0] + ((double) iy)*(grid->spacing);
        for (iz=0;iz<nz;iz++) {
          index = ny*nz*ix + nz*iy + iz;
          if (!(is_mask_bit_set(mask,index))) {
            pos[2] = grid->flimit[2][0] + ((double) iz)*(grid->spacing);
            score = calculate_score(probe, pos, system, buffer_field);
            matrix[ix][iy][iz] = score;
          }
        }
      }
    }
  } else {  // save n positions for each orientation or not save poses
    for (int i=0; i<n_orientations; i++) {
      printf("orientation %d out of %d\n", i+1, n_orientations);
      coords2molecule(probe_orientations[i], molecule);
      // set the center to the origin as the orientations have their centers on the origin
      null_vector(molecule->geometric_center);
      probe->orientation_id = i;

      if (fragmap_settings->solvate) {
        move_grid_center(fragmap_settings->ligand_water_field->grid, COORDS_ORIGIN);
        fragmap_settings->ligand_water_field = get_atom_field(fragmap_settings->ligand_water_field, molecule, fragmap_settings->water_probe, 1);
        write_insight_map("ligand_water.grd", fragmap_settings->ligand_water_field, 0);
        write_pdb_molecule(molecule, "ligand.pdb");
      }
      update_molecule_score_matrix(probe, field, fragmap_data->buffer_field, system);
    }
  }

  // TODO:
  // free memory?
  system->selection = selection;
  reset_system(system,0);
  return(field);
}


void update_molecule_score_matrix(MOLECULE_PROBE *probe, MAP *field, MAP *buffer_field, SYSTEM *system) {

  double ***matrix = (double***) field->matrix;
  unsigned char *mask = field->mask;
  GRID *grid = field->grid;
  MOLECULE *molecule = probe->molecule;
  int use_precalculated_grids = fragmap_settings->precalculate;
  PLI_SFUNC *sfunc = system->settings->sfunc;
  CLUSTER_SETTINGS *cluster_settings = fragmap_settings->cluster_settings;

  int nx = grid->npoints[0];
  int ny = grid->npoints[1];
  int nz = grid->npoints[2];

  int n_waters = 0;

  // reset scores for the positions for each orientation
  if (fragmap_data->orientation_pose_list != NULL) {
    fragmap_data->orientation_pose_list->n_poses = 0;
    fragmap_data->orientation_pose_list->worst_score_pose = NULL;
  }

  POSELIST_LITE ***poselist_matrix;
  if (cluster_settings->save_peak_clusters) {
    poselist_matrix = (POSELIST_LITE***) fragmap_data->poselist_map->matrix;
  }

  //int ixc, iyc, izc;
  double pos[4];
  int index;
  //double score;
  double water_score = 0.0;
  for (int ix=0; ix<nx; ix++) {
    pos[0] = grid->flimit[0][0] + ((double) ix)*(grid->spacing);
    for (int iy=0; iy<ny; iy++) {
      pos[1] = grid->flimit[1][0] + ((double) iy)*(grid->spacing);
      for (int iz=0; iz<nz; iz++) {
        index = ny*nz*ix + nz*iy + iz;
        if (is_mask_bit_set(mask, index)) {
          continue;
        }
        pos[2] = grid->flimit[2][0] + ((double) iz)*(grid->spacing);
        probe->score = sfunc->bad_score * (system->settings->sfunc->order == DESCENDING)? -1.0: 1.0;
        move_molecule_centre(molecule, pos);

        if (molecule_accessible(buffer_field, probe->molecule)) {
          // score
          if (use_precalculated_grids) {
            probe->score = score_probe_using_fields(probe, system->settings);
            N_SCORE_GRID++;
            if ((fragmap_settings->rescore) && (probe->score < fragmap_settings->rescore_cutoff)) {
              score_system(system, sfunc);
              probe->score = system->score;
              N_SCORE_FULL++;
            }
          } else {
            N_SCORE_FULL++;
            score_system(system, sfunc);
            probe->score = system->score;
          }

          // solvate
          if (fragmap_settings->solvate && probe->score < fragmap_settings->save_cutoff) {
            move_grid_center(fragmap_settings->ligand_water_field->grid, pos);
            // adjust ixc, etc... according to protein water field
            int ixc = ix + (fragmap_settings->protein_water_field->grid->npoints[0] - nx)/2;
            int iyc = iy + (fragmap_settings->protein_water_field->grid->npoints[1] - ny)/2;
            int izc = iz + (fragmap_settings->protein_water_field->grid->npoints[2] - nz)/2;

            double score_solvated = solvate_pose(system, ixc, iyc, izc, fragmap_settings->docked_waters);

            if (n_waters) {
              water_score = score_solvated - probe->score;
              probe->score = score_solvated;
            }
          } /* end solvate */

          // update score

          if (system->settings->sfunc->order == DESCENDING) {
            probe->score *= -1.0;
          }

          if (probe->score < matrix[ix][iy][iz]) {
            matrix[ix][iy][iz] = probe->score;
          }

          // save the probe
          if (cluster_settings->save && (probe->score < fragmap_settings->save_cutoff)) {
            N_SAVE++;
            if (cluster_settings->n_poses_per_orientation) {
              save_pose(fragmap_data->orientation_pose_list, fragmap_data->orientation_allocated_poses, probe);
            } else if (cluster_settings->save_poses) {
              save_pose(fragmap_data->saved_pose_list, fragmap_data->saved_allocated_poses, probe);
            }
            if (cluster_settings->save_diverse_clusters) {
              save_diverse_pose(fragmap_data->div_cluster_set, probe, fragmap_settings->docked_waters);
            }
            if (cluster_settings->save_peak_clusters) {
              add_pose_to_poselist(&(poselist_matrix[ix][iy][iz]), probe->orientation_id, probe->score);
            }
          }
        } /* end molecule accessible */
      } /* end iz*/
    } /* end iy*/
  } /* end ix*/

  if (cluster_settings->n_poses_per_orientation && cluster_settings->save_poses) {
    PROBE_POSE *pose;
    for (int i=0; i<fragmap_data->orientation_pose_list->n_poses; i++) {
      pose = fragmap_data->orientation_pose_list->poses[i];
      probe->orientation_id = pose->orientation_index;
      probe->score = pose->score;
      save_pose(fragmap_data->saved_pose_list, fragmap_data->saved_allocated_poses, probe);
    }
    probe->orientation_id = -1;
    probe->score = 0.0;
  }
}



MAP* get_atom_field(MAP *map_in, MOLECULE *molecule, FIELD_PROBE *atom_probe, int da_atoms_only) {

  if ((! (molecule->flags & LIGAND_MOLECULE)) && (! (molecule->flags & PROTEIN_MOLECULE))) {
    error_fn("%s: only protein or ligand molecule can be passed to the function", __func__);
  }

  SYSTEM *system = molecule->molsystem;

  if ((! system->ligand) && (molecule->flags & LIGAND_MOLECULE)) {
    system->ligand = molecule;
  } else if ((! system->protein) && (molecule->flags & PROTEIN_MOLECULE)) {
    system->protein = molecule;
  }

  // reset mask
  init_mask(map_in->mask,map_in->grid);
  mask_molecule(molecule,map_in->mask,map_in->grid,CONTACT_MASK);

  if (da_atoms_only) {
    ATOMLIST *hbond_atoms = molecule2hbonding_atoms(molecule);
    mask_distant_points(map_in, NULL, NULL, hbond_atoms, 4.0, 45.0);

    free(hbond_atoms);
  }

  map_in = calculate_field(system, atom_probe, map_in, system->settings->sfunc);

  double max_score = (molecule->flags & LIGAND_MOLECULE) ? SOLVATE_LIG_MAX_SCORE: SOLVATE_PROTEIN_MAX_SCORE;
  mask_map_cutoff(map_in, max_score, "above");

  return(map_in);
}


double solvate_pose(SYSTEM *system, int ixc, int iyc, int izc, DOCKED_WATERS *docked_waters) {
  double old_score = system->score;

  fragmap_settings->ligand_water_field_backup = copy_map(fragmap_settings->ligand_water_field_backup,fragmap_settings->ligand_water_field);

  MAP *protein_field = fragmap_settings->protein_water_field;
  MAP *ligand_field = fragmap_settings->ligand_water_field_backup;

  int nxp = protein_field->grid->npoints[0];
  int nyp = protein_field->grid->npoints[1];
  int nzp = protein_field->grid->npoints[2];
  int nxl = ligand_field->grid->npoints[0];
  int nyl = ligand_field->grid->npoints[1];
  int nzl = ligand_field->grid->npoints[2];

  double ***protein_matrix = (double***) protein_field->matrix;
  double ***ligand_matrix = (double***) ligand_field->matrix;

  // combine scores
  for (int ixp = ixc - (nxl -1) / 2, ixl=0; ixl < nxl; ixp++, ixl++) {
    for (int iyp = iyc - (nyl -1) / 2, iyl=0; iyl < nyl; iyp++, iyl++) {
      for (int izp = izc - (nzl -1) / 2, izl=0; izl < nzl; izp++, izl++) {
        int indexp = nyp*nzp*ixp + nzp*iyp + izp;
        int indexl = nyl*nzl*ixl + nzl*iyl + izl;
        // can this be done in one go?
        if (! is_mask_bit_set(ligand_field->mask, indexl)) {
          if (is_mask_bit_set(protein_field->mask, indexp)) {
            set_mask_bit(ligand_field->mask, indexl);
          } else {
            ligand_matrix[ixl][iyl][izl] += protein_matrix[ixp][iyp][izp];
          }
        }
      }
    }
  }

  // place waters
  // TODO: can we do this only for the affected atoms to make it faster?
  reset_system(system, 0);
  docked_waters->n_waters = 0;
  double score_grid = 0.0;
  int water_point[3];
  for (int i=0; i<fragmap_settings->n_max_waters; i++) {
    int water_index = map_extreme_index(ligand_field, 1, "min", water_point, water_point+1, water_point+2);
    if (water_index == -1) {
      break;
    }
    score_grid += ligand_matrix[water_point[0]][water_point[1]][water_point[2]];

    MOLECULE *water = docked_waters->mols+ docked_waters->n_waters;
    water->atom->position[0] = ligand_field->grid->flimit[0][0] + ((double) water_point[0])*(ligand_field->grid->spacing);
    water->atom->position[1] = ligand_field->grid->flimit[1][0] + ((double) water_point[1])*(ligand_field->grid->spacing);
    water->atom->position[2] = ligand_field->grid->flimit[2][0] + ((double) water_point[2])*(ligand_field->grid->spacing);


    water->atom->flags |= ATOM_NEW;
    water->atom->cstatus |= ATOM_NOTHING_CALCULATED;
    add_molecule_to_list(system->molecule_list, water);
    add_atom_to_list(system->selection, water->atom, water->atom->flags);
    mask_molecule(water,ligand_field->mask,ligand_field->grid,CONTACT_MASK);
    docked_waters->n_waters++;
  }

  // score system
  double new_score;

  if (docked_waters->n_waters) {
    if (fragmap_settings->precalculate) {
      new_score = score_grid;
    } else {
      score_system(system, system->settings->sfunc);
      new_score = system->score - WATER_PENALTY * docked_waters->n_waters;
    }
  }

  // remove waters from the system
  for (int i=0; i<docked_waters->n_waters; i++) {
    remove_atom_from_list(system->selection, (docked_waters->mols +i)->atom);
    remove_molecule_from_list(system->molecule_list, docked_waters->mols + i);
  }

  if (docked_waters->n_waters) {
    if (new_score < old_score) {
      return (new_score);
    } else {
      docked_waters->n_waters = 0;
      return(old_score);
    }
  } else {
    return(old_score);
  }
}


void alloc_fragmap_data() {
  int i;

  CLUSTER_SETTINGS *cluster_settings;

  cluster_settings = fragmap_settings->cluster_settings;

  fragmap_data = malloc(sizeof(struct FragmapData));

  if (fragmap_data == NULL) {
    error_fn("%s: Error allocating fragmap data");
  }


  if (cluster_settings->n_poses_per_gridpoint > 0) {
    fragmap_data->gridpoint_pose_list = alloc_pose_list(cluster_settings->n_poses_per_gridpoint);
    fragmap_data->gridpoint_allocated_poses = alloc_poses(cluster_settings->n_poses_per_gridpoint, 0, 0);
  } else {
    fragmap_data->gridpoint_pose_list = NULL;
    fragmap_data->gridpoint_allocated_poses = NULL;
  }

  if (cluster_settings->n_poses_per_orientation > 0) {
    fragmap_data->orientation_pose_list = alloc_pose_list(cluster_settings->n_poses_per_orientation);
    fragmap_data->orientation_allocated_poses = alloc_poses(cluster_settings->n_poses_per_orientation, 0, 0);
  } else {
    fragmap_data->orientation_pose_list = NULL;
    fragmap_data->orientation_allocated_poses = NULL;
  }

  if (cluster_settings->save_poses) {
    fragmap_data->saved_pose_list = alloc_pose_list(cluster_settings->n_saved_poses);
    fragmap_data->saved_allocated_poses = alloc_poses(cluster_settings->n_saved_poses, 0, 0);
  } else {
    fragmap_data->saved_pose_list = NULL;
    fragmap_data->saved_allocated_poses = NULL;
  }

  if (cluster_settings->save_diverse_clusters) {
    int n_waters = fragmap_settings->solvate ? fragmap_settings->n_max_waters : 0;
    //int n_atoms = (cluster_settings->dist_type == RMSD_DIST) ? fragmap_settings->probe->molecule->selection->natoms: 0;
    fragmap_data->div_cluster_set = alloc_div_cluster_set(cluster_settings->n_diverse_poses, 0, n_waters);
    //fragmap_data->diverse_allocated_poses = alloc_poses(cluster_settings->n_diverse_poses, n_atoms, n_waters);
  } else {
    fragmap_data->div_cluster_set = NULL;
    //fragmap_data->diverse_allocated_poses = NULL;
  }

  if (cluster_settings->save_peak_clusters) {
    fragmap_data->peak_cluster_set = NULL;
    fragmap_data->peak_allocated_poses = alloc_poses(cluster_settings->n_peaks * cluster_settings->peak_cluster_size, 0, 0);
  }


  fragmap_data->field = NULL;
  fragmap_data->buffer_field = NULL;
  fragmap_data->poselist_map = NULL;
}


static void output_fragmap(MAP *field, MAP *poselist_map) {
  MOLECULE_PROBE *probe;
  CLUSTER_SETTINGS *cluster_settings;

  probe = fragmap_settings->probe;
  cluster_settings = fragmap_settings->cluster_settings;

  write_insight_map(fragmap_settings->grid_file,field,0);

  // save probe poses
  if (cluster_settings->save_poses) {
    write_pose_list(fragmap_data->saved_pose_list, probe->molecule, probe->orientations);
  }

  if (cluster_settings->save_diverse_clusters) {
    if (fragmap_settings->solvate) {
      write_cluster_set(fragmap_data->div_cluster_set, "diverse_cluster", probe, fragmap_settings->docked_waters->mols);
    } else {
      write_cluster_set(fragmap_data->div_cluster_set, "diverse_cluster", probe, NULL);
    }
    //test_div_cluster_score(div_cluster_set, system, probe_orientations, probe);
  }

  if (cluster_settings->save_peak_clusters) {
    fragmap_data->peak_cluster_set = get_peak_clusters(fragmap_data->poselist_map,
                                                       fragmap_data->field, probe,
                                                       probe->orientations,
                                                       cluster_settings,
                                                       fragmap_data->peak_allocated_poses,
                                                       fragmap_settings->probe->auto_mappings,
                                                       fragmap_settings->probe->n_mappings);
    if (fragmap_settings->solvate) {
      write_cluster_set(fragmap_data->peak_cluster_set, "peak_cluster", probe, fragmap_settings->docked_waters->mols);
    } else {
      write_cluster_set(fragmap_data->peak_cluster_set, "peak_cluster", probe, NULL);
    }

  }
}




MAP_ISLAND* get_ligand_island(MAP_ISLAND *islands, int n_islands, MAP *field, MOLECULE *ligand) {
  int i, nx, ny, nz, index;
  int ipos[3];
  ATOM *atom;
  MAP_ISLAND *island;

  nx = field->grid->npoints[0];
  ny = field->grid->npoints[1];
  nz = field->grid->npoints[2];

  for (i=0, atom=ligand->atom; i<ligand->natoms; i++, atom++) {
    if (pos2grid(atom->position, ipos, field->grid)) {
      error_fn("%s: reference ligand is not within the grid bound", __func__);
    }
    index = ipos[0]*ny*nz + ipos[1]*nz+ipos[2];
    island = index2island(index, islands, n_islands);
    if (island != NULL) {
      return (island);
    }
  }
  return NULL;
}


void mask_ref_ligand(MAP *field, MOLECULE *ref_ligand, double max_dist, int mask_islands) {
  int nx, ny, nz, ix, iy, iz, i, index, j;
  GRID *grid;
  unsigned char *mask;
  double ***matrix;
  int n_islands;
  MAP_ISLAND *islands;
  MAP_ISLAND *ligand_island,*island;
  int ipos[3];


  grid = field->grid;
  mask = field->mask;
  matrix = (double***) field->matrix;
  nx = grid->npoints[0];
  ny = grid->npoints[1];
  nz = grid->npoints[2];

  // mask other islands
  mask_distant_points(field, NULL, ref_ligand, NULL, max_dist, -1.0);


  if (mask_islands) {
    // set accessible and nearby points to 1.0 so that with that level, islands can be found
    for (ix=0;ix<nx;ix++) {
      for (iy=0;iy<ny;iy++) {
        for (iz=0;iz<nz;iz++) {
          index = ny*nz*ix + nz*iy + iz;
          if (!(is_mask_bit_set(mask,index))) {
            matrix[ix][iy][iz] = 1.0;
          }
        }
      }
    }


    islands = map2islands(field, 1.0, &n_islands);

    ligand_island = get_ligand_island(islands, n_islands, field, ref_ligand);
    if (ligand_island == NULL) {
      error_fn("%s: reference ligand is masked by protein", __func__);
    }

    memset(mask, -1, (nx*ny*nz)/8 + 1);

    for (i=0; i<ligand_island->n_points; i++) {
      unset_mask_bit(mask, ligand_island->points[i]);
    }


    // write island map
    init_3d_matrix(field);
    for (i=0; i<n_islands; i++) {
      island = islands+i;
      for (j=0; j<island->n_points; j++) {
        map_index2grid(island->points[j], &ix, &iy, &iz, field->grid);
        matrix[ix][iy][iz] = (double) island->id;
      }
    }
    //write_insight_map("islands.grd", field, 0);

    init_3d_matrix(field);
  }
}




static double calculate_score(MOLECULE_PROBE *probe, double pos[4], SYSTEM *system, MAP *buffer_field) {
  int i;
  double min_score, **u_axis;
  PROBE_POSE probe_pose;
  int n_orientations;
  PROBE_POSE *pose;
  PLI_SFUNC *sfunc;
  CLUSTER_SETTINGS *cluster_settings;

  sfunc = system->settings->sfunc;
  cluster_settings = fragmap_settings->cluster_settings;

  min_score = sfunc->bad_score * (sfunc->order == DESCENDING)? -1.0: 1.0;
  // reset saved poses for grid point
  if (cluster_settings->n_poses_per_gridpoint) {
    fragmap_data->gridpoint_pose_list->n_poses = 0;
    fragmap_data->gridpoint_pose_list->worst_score_pose = NULL;
  }

  n_orientations = probe->n_orientations;

  for (i=0; i<n_orientations; i++) {
    coords2molecule(probe->orientations[i], probe->molecule);
    probe->orientation_id = i;
    // set the center to the origin as the orientations have their centers on the origin
    null_vector(probe->molecule->geometric_center);
    move_molecule_centre(probe->molecule, pos);
    probe->score = sfunc->bad_score;
    if (molecule_accessible(buffer_field, probe->molecule)) {
      if (fragmap_settings->precalculate) {
        probe->score = score_probe_using_fields(probe, system->settings);
        if (fragmap_settings->rescore && (probe->score < fragmap_settings->rescore_cutoff)) {
          score_system(system, sfunc);
          probe->score = system->score;
        }
      } else {
        score_system(system, sfunc);
        probe->score = system->score;
      }
    }

    if (sfunc->order == DESCENDING) {
      probe->score *= -1;
    }

    if (probe->score < min_score) {
      min_score = probe->score;
    }

    // save pose
    if ((probe->score < fragmap_settings->save_cutoff) && (cluster_settings->save)) {
      if (cluster_settings->n_poses_per_gridpoint) {
        save_pose(fragmap_data->gridpoint_pose_list, fragmap_data->gridpoint_allocated_poses, probe);
      } else if (cluster_settings->save_poses) {
        save_pose(fragmap_data->saved_pose_list, fragmap_data->saved_allocated_poses, probe);
      }

      if (cluster_settings->save_diverse_clusters) {
        save_diverse_pose(fragmap_data->div_cluster_set, probe, NULL);
      }
    }

    // save poses for this grid point to the saved poses
    if (cluster_settings->n_poses_per_gridpoint && cluster_settings->save_poses) {
      for (i=0; i<fragmap_data->gridpoint_pose_list->n_poses; i++) {
        pose = fragmap_data->gridpoint_pose_list->poses[i];
        probe->score = pose->score;
        probe->orientation_id = pose->orientation_index;
        save_pose(fragmap_data->saved_pose_list, fragmap_data->saved_allocated_poses, probe);
      }
    }

  }

  return(min_score);
}




static int get_optimal_n_rot_lebedev(double **lebedev_points, int n_points) {
  double min_angle = 360.0;
  int i;
  for (i=0; i<n_points; i++) {
    if (fabs(lebedev_points[i][4]) > 1.0E-6 && fabs(lebedev_points[i][4] < min_angle)) {
      min_angle = fabs(lebedev_points[i][4]);
    }
  }
  return((int)floor(360.0/min_angle));
}


static void init_fragmap_settings(SYSTEM *system) {

  SETTINGS *settings = system->settings;

  // alloc
  if (fragmap_settings == NULL) {
    fragmap_settings = (struct FragmapSettings*) malloc(sizeof(struct FragmapSettings));
  }
  if (fragmap_settings == NULL) {
      error_fn("%s: out of memory allocating fragmap settings", __func__);
  }

  // cluster settings
  fragmap_settings->cluster_settings = get_cluster_settings();
  CLUSTER_SETTINGS *cluster_settings = fragmap_settings->cluster_settings;

  // grid settings:
  fragmap_settings->grid_spacing = (params_get_parameter("grid_spacing"))->value.d;
  fragmap_settings->grid_padding = (params_get_parameter("grid_padding"))->value.d;
  strcpy(fragmap_settings->grid_file,(params_get_parameter("grid_file"))->value.s);

  // molecule probe:
  if (strcmp((params_get_parameter("probe"))->value.s, "undefined") == 0) {
    error_fn("%s: probe is required for %s", __func__, (params_get_parameter("mode"))->value.s);
  }
  strcpy(fragmap_settings->probe_name,  params_get_parameter("probe")->value.s);
  fragmap_settings->probe = get_molecule_probe(fragmap_settings->probe_name, settings, params_get_parameter("skip_probe_atoms")->value.s);
  MOLECULE *molecule = fragmap_settings->probe->molecule;

  // automorphs
  do_ullman(molecule, molecule, add_automorph);

  // rotational parameters:
  fragmap_settings->n_lebedev_points = (params_get_parameter("n_lebedev_points"))->value.i;
  fragmap_settings->n_dihedral_steps = (params_get_parameter("n_dihedral_steps"))->value.i;
  fragmap_settings->scan_dihedral = (params_get_parameter("scan_dihedral"))->value.i;
  fragmap_settings->lebedev_points = read_lebedev_sphere_points(fragmap_settings->n_lebedev_points);
  fragmap_settings->n_rot_lebedev_axis = (params_get_parameter("n_rot_lebedev_axis"))->value.i;
  if (heavy_atom_count(molecule) == 2) {
    fragmap_settings->n_rot_lebedev_axis = 1;
  }

  if (fragmap_settings->n_rot_lebedev_axis == -1) {
    fragmap_settings->n_rot_lebedev_axis = get_optimal_n_rot_lebedev(fragmap_settings->lebedev_points, fragmap_settings->n_lebedev_points);
  }

  // get the probe orientation
  if ((fragmap_settings->scan_dihedral) && (molecule->n_torsions)){
    if (molecule->n_torsions > 1) {
      warning_fn("More than one rotatable bonds found, only the first one being considered for rotation");
    }
    fragmap_settings->probe->orientations = get_molecule_orientations(molecule,
                                                                     fragmap_settings->n_lebedev_points,
                                                                     fragmap_settings->lebedev_points,
                                                                     fragmap_settings->n_rot_lebedev_axis,
                                                                     fragmap_settings->n_dihedral_steps,
                                                                     &(fragmap_settings->probe->n_orientations));
  } else {
    fragmap_settings->probe->orientations = get_molecule_orientations(molecule,
                                                                     fragmap_settings->n_lebedev_points,
                                                                     fragmap_settings->lebedev_points,
                                                                     fragmap_settings->n_rot_lebedev_axis,
                                                                     0,
                                                                     &(fragmap_settings->probe->n_orientations));
  }
  //write_probe_orientations(fragmap_settings->probe->molecule, fragmap_settings->probe->orientations, fragmap_settings->probe->n_orientations, "probe_orientation");

  // move molecule to origin as the orientation are centered on the origin
  null_vector(molecule->geometric_center);

  // remove similar orientations due to symmetry
  remove_redundant_orientations(fragmap_settings->probe->orientations, &(fragmap_settings->probe->n_orientations),
                                molecule, fragmap_settings->probe->auto_mappings, fragmap_settings->probe->n_mappings, 0.01);
  // score and probe save settings
  fragmap_settings->precalculate = (params_get_parameter("precalculate"))->value.i;
  fragmap_settings->save_intermediate_grids = (params_get_parameter("save_intermediate_grids"))->value.i;
  fragmap_settings->rescore = (params_get_parameter("fragmap_rescore"))->value.i;

  // reference ligand
  fragmap_settings->ref_ligand_max_dist = (params_get_parameter("ref_ligand_max_dist"))->value.d;
  // ref_ligand
  if (strcmp((params_get_parameter("ref_ligand"))->value.s, "undefined") != 0) {
    if (is_readable_file((params_get_parameter("ref_ligand"))->value.s)) {
      fragmap_settings->ref_ligand = read_molecule((params_get_parameter("ref_ligand"))->value.s);
    } else {
      error_fn("%s: Molecule could not be read from the ref_ligand (%s)", __func__, (params_get_parameter("ref_ligand"))->value.s);
    }
  } else {
    fragmap_settings->ref_ligand = NULL;
  }

  // mask molecule
  fragmap_settings->mask_mol = NULL;
  if (strcmp((params_get_parameter("mask"))->value.s, "undefined") != 0) {
    if (is_readable_file((params_get_parameter("mask"))->value.s)) {
      fragmap_settings->mask_mol = read_molecule((params_get_parameter("mask"))->value.s);
      fragmap_settings->mask_mol->flags |= LIGAND_MOLECULE;
      fragmap_settings->mask_mol->use_hydrogens = 1;
      prep_molecule(fragmap_settings->mask_mol, system->settings);
    } else {
      error_fn("%s: Molecule could not be read from \"%s\" file", __func__, (params_get_parameter("mask"))->value.s);
    }
  }

  // cluster_settings
  if (cluster_settings->n_poses_per_orientation > 0) {
    cluster_settings->n_saved_poses = fragmap_settings->n_lebedev_points * fragmap_settings->n_rot_lebedev_axis * cluster_settings->n_poses_per_orientation;
  }
  fragmap_settings->save_cutoff = cluster_settings->save_cutoff_atom * fragmap_settings->probe->molecule->selection->natoms;
  fragmap_settings->rescore_cutoff = params_get_parameter("fragmap_rescore_cutoff")->value.d * fragmap_settings->probe->molecule->selection->natoms;

  // fragmap solvate
  MAP *map;
  fragmap_settings->water_probe = NULL;
  map = NULL;

  fragmap_settings->solvate = params_get_parameter("fragmap_solvate")->value.i;

  // check if the probe has hbonding atoms
  if (fragmap_settings->solvate) {
    if (! has_hbonding_atom(molecule)) {
      warning_fn("%s: No hbonding atoms found in the probe; fragmap_solvate option is being unset", __func__);
      fragmap_settings->solvate = 0;
    }
  }

  fragmap_settings->ligand_grid_padding = 0.0;
  fragmap_settings->docked_waters = NULL;
  // alloc water map for the probe
  if (fragmap_settings->solvate) {
    fragmap_settings->ligand_grid_padding = 6.2;
    fragmap_settings->n_max_waters = 3;
    ATOM_TYPE *type = get_atom_type("H2O", settings->atom_typing_scheme);
    fragmap_settings->water_probe = create_atom_type_probe(type, settings, NULL);
    map = new_map("ligand water field","double",NULL);
    map->grid = coords2grid(fragmap_settings->probe->orientations, fragmap_settings->probe->n_orientations,
                            molecule, fragmap_settings->grid_spacing, fragmap_settings->ligand_grid_padding, molecule->geometric_center);
    alloc_map_matrix(map);
    fragmap_settings->ligand_water_field_backup = NULL;
    fragmap_settings->docked_waters = malloc(sizeof(DOCKED_WATERS));
    fragmap_settings->docked_waters->n_waters = 0;
    fragmap_settings->docked_waters->mols = create_atom_type_molecules(type, fragmap_settings->n_max_waters, "water", settings);
    fragmap_settings->docked_waters->score = 0.0;
  }
  fragmap_settings->ligand_water_field = map;
}


static void add_automorph(BOOLEAN **m) {
  MOLECULE *molecule;
  int n_atoms;
  BOOLEAN *automorph_matrix;

  if (fragmap_settings->probe->n_mappings == MAX_AUTOMORPHS) {
    error_fn("%s: the probe molecule has more that 100 automorphs", __func__);
  }

  molecule = fragmap_settings->probe->molecule;
  n_atoms = molecule->selection->natoms;

  if (fragmap_settings->probe->n_mappings == 0) {
    fragmap_settings->probe->n_alloc_mappings = 2;
    fragmap_settings->probe->auto_mappings = (BOOLEAN***) malloc(fragmap_settings->probe->n_alloc_mappings * sizeof(BOOLEAN **));
  }

  if (fragmap_settings->probe->n_mappings == fragmap_settings->probe->n_alloc_mappings) {
    fragmap_settings->probe->n_alloc_mappings += 2;
    fragmap_settings->probe->auto_mappings = (BOOLEAN***) realloc(fragmap_settings->probe->auto_mappings,
                                                           fragmap_settings->probe->n_alloc_mappings * sizeof(BOOLEAN **));
  }
  fragmap_settings->probe->auto_mappings[fragmap_settings->probe->n_mappings] = alloc_2d_bool_matrix(n_atoms, n_atoms);
  copy_2d_bool_matrix(fragmap_settings->probe->auto_mappings[fragmap_settings->probe->n_mappings], m, n_atoms, n_atoms);
  fragmap_settings->probe->n_mappings++;
}


static void remove_redundant_orientations(ATOM_COORDS **orientations, int *n_orientations, MOLECULE *molecule, BOOLEAN ***mappings, int n_mappings, double rms_cutoff) {
  int i,j;

  int selected_indices[molecule->selection->natoms];
  get_selected_atom_indices(molecule, selected_indices);

  int remove_ids[*n_orientations];
  int n_remove;
  double rms;


  n_remove = 0;

  time_t start;
  time(&start);

  // find which orientations to remove
  for (i=0; i<((*n_orientations)-1); i++) {
    if (is_in_iarray(i, remove_ids, n_remove)) {
      continue;
    }
    for (j=(i+1); j<(*n_orientations); j++) {
      if (is_in_iarray(j, remove_ids, n_remove)) {
        continue;
      }
      rms = rms_coords_atleat(orientations[i], orientations[j], selected_indices, selected_indices, molecule->selection->natoms, mappings, n_mappings, rms_cutoff);
      if (rms < rms_cutoff) {
        //printf("remove %d, rms with %d: %.2f\n", j+1, i+1, rms);
        remove_ids[n_remove] = j;
        n_remove++;
      }
    }
  }
//  write_probe_orientations(molecule, orientations, *n_orientations, "mol");
  // free redundant orientations
  int n_freed = 0;
  for (i=0; i<(*n_orientations); i++) {
    if (is_in_iarray(i, remove_ids, n_remove)) {
      free(orientations[i]);
      orientations[i] = NULL;
      n_freed++;
    } else {
      // move the orientation
      orientations[i-n_freed] = orientations[i];
    }
  }
  *n_orientations -= n_remove;
}

double score_probe_using_fields(MOLECULE_PROBE *probe, SETTINGS *settings) {
  int i, j;
  double total_score = 0.0;
  double score;
  MOLECULE *molecule;
  ATOM *atom;
  MAP *field = NULL;
  molecule = probe->molecule;
  for (i=0, atom=molecule->atom; i<molecule->natoms; i++, atom++) {
    if ((atom->element->id == HYDROGEN) || (atom->flags & SKIP_ATOM)) {
      continue;
    }

    field = probe->atom_maps[i];
    if (field == NULL) {
      error_fn("%s: field for atom type %s not found\n", __func__, atom->type->name);
    }

    score = trilinear_interpolate(field, atom->position, settings->sfunc->bad_score);
    if (score >= settings->sfunc->bad_score) {
      return(score);
    } else {
      total_score += score;
    }
  }
  return(total_score);
}


void calculate_atom_fields_for_probe(MOLECULE_PROBE *probe, SETTINGS *settings, SYSTEM *system) {
  int i, j;
  char atom_field_file[MAX_LINE_LEN];
  MOLECULE *molecule;
  ATOM *atom;
  double mol_rad;
  double dist, buffer_padding;
  int n_types;

  ATOM_TYPE_LIST *type_list;
  MAP *first_atom_field, *field;
  int first_atom_field_used;
  FIELD_PROBE *atom_probe;
  PLI_SFUNC *sfunc;
  GRID *grid;


  molecule = probe->molecule;
  //
  sfunc = settings->sfunc;
  // map pointers for every atom
  probe->atom_maps = calloc(molecule->natoms, sizeof(MAP*));
  for (i=0; i<molecule->natoms; i++) {
    probe->atom_maps[i] = NULL;
  }

  // add buffer to grid
  mol_rad = molecule_radius(molecule);

  if ((ceil(mol_rad) - mol_rad) > fragmap_settings->grid_spacing) {
    buffer_padding = ceil(mol_rad);
  } else{
    buffer_padding = ceil(mol_rad) + fragmap_settings->grid_spacing;
  }

  (params_get_parameter("grid_padding"))->value.d += buffer_padding;

  // allocate the first atom atom field now so that we can use grid bounds for other atoms

  grid = system_grid(system,(params_get_parameter("grid_spacing"))->value.d,(params_get_parameter("grid_padding"))->value.d);
  first_atom_field = new_map("atom map","double",grid);

  mask_molecule(system->protein,first_atom_field->mask,first_atom_field->grid,CONTACT_MASK);
  first_atom_field_used = 0;

  for (i=0; i<molecule->natoms; i++) {
    atom = molecule->atom+i;
    if ((atom->element->id == HYDROGEN) || (atom->flags & SKIP_ATOM)) {
      continue;
    }
    // already calculated field for the atom type
    if (probe->atom_maps[i] != NULL) {
      continue;
    }
    atom_probe = create_atom_type_probe(atom->type, settings, NULL);
    if (!first_atom_field_used) {
      field = first_atom_field;
      first_atom_field_used = 1;
    } else {

      field = new_map("atom map","double",copy_grid(NULL,grid));

      copy_mask_from_field(field->mask, first_atom_field);
    }
    field = calculate_field(system, atom_probe, field, sfunc);
    probe->atom_maps[i] = field;
    // check if the other atoms have the same atom type
    for (j=(i+1); j<molecule->natoms; j++) {
      if ((molecule->atom+j)->type == atom->type) {
        probe->atom_maps[j] = field;
      }
    }
    if (fragmap_settings->save_intermediate_grids) {
      (atom->type->id != -1) ? sprintf(atom_field_file, "atom_type_%d.grd", atom->type->id): sprintf(atom_field_file, "element_type_%d.grd", atom->type->element->id);
      write_insight_map(atom_field_file, field, 0);
    }
    // free probe now
  }

  (params_get_parameter("grid_padding"))->value.d -= buffer_padding;
}


void copy_mask_from_field(unsigned char* mask_to, MAP *field_from) {
  int nchars;

  nchars = (field_from->grid->npoints[0]*field_from->grid->npoints[1]*field_from->grid->npoints[2])/8 + 1;

  if (mask_to == NULL) {
    // alloc
    mask_to = (unsigned char*) calloc(nchars,sizeof(unsigned char));
    if (mask_to == NULL) {
      error_fn("%s: can not allocate mask", __func__);
    }
  }

  memcpy(mask_to, field_from->mask, nchars*sizeof(unsigned char));

}


static void free_molecule_probe(MOLECULE_PROBE *probe, SETTINGS *settings) {
  if (probe) {
    MOLECULE *molecule = probe->molecule;
    if (probe->molecule) {
      unprep_molecule(molecule, settings);
      free(molecule);
    }
    // TODO: free other structs
    free(probe);
  }
}


void copy_field (MAP *source, MAP **targetp) {
  int nx, ny, nz, ix, iy, iz;

  if (*targetp == NULL) {
    
    *targetp = new_map(source->name,source->type->name,copy_grid(NULL,source->grid));

    if (*targetp == NULL) {
      error_fn("%s: out of memory allocating field", __func__);
    }
    (*targetp)->grid = malloc(sizeof(GRID));
    if ((*targetp)->grid == NULL) {
      error_fn("%s: out of memory allocating grid", __func__);
    }
    memcpy((*targetp)->grid, source->grid, sizeof(GRID));
    (*targetp)->type = source->type;
    alloc_map_matrix(*targetp);
  } else {
    if (source->type != (*targetp)->type) {
      error_fn("%s: can not copy from a different type of map", __func__);
    }
    memcpy((*targetp)->grid, source->grid, sizeof(GRID));
  }

  nx = source->grid->npoints[0];
  ny = source->grid->npoints[1];
  nz = source->grid->npoints[2];

  double ***target_matrix = (double***) (*targetp)->matrix;
  double ***source_matrix = (double***) source->matrix;

  copy_3d_double_matrix(source_matrix, target_matrix, source->grid->npoints);
  int nchars = (nx*ny*nz)/8 + 1;
  memcpy((*targetp)->mask, source->mask, nchars*sizeof(unsigned char));
}

void get_min_score_point(MAP* field, int *min_score_point) {
  int ix, iy, iz, nx, ny, nz, index;
  GRID *grid;
  unsigned char *mask;
  double ***matrix;
  double score, min_score;

  min_score = 100000.0;
  matrix = (double***) field->matrix;
  mask = field->mask;
  grid = field->grid;
  nx = grid->npoints[0];
  ny = grid->npoints[1];
  nz = grid->npoints[2];

  for (ix=0;ix<nx;ix++) {
    for (iy=0;iy<ny;iy++) {
      for (iz=0;iz<nz;iz++) {
        index = ny*nz*ix + nz*iy + iz;
        if (!(is_mask_bit_set(mask,index))) {
          score = matrix[ix][iy][iz];
          if (score<min_score) {
            min_score = score;
            min_score_point[0] = ix;
            min_score_point[1] = iy;
            min_score_point[2] = iz;
          }
        }
      }
    }
  }
}


GRID* copy_grid_extend(GRID *to_grid, GRID *from_grid, double extend_by) {
  int i, j, n_extend;

  if (to_grid == NULL) {
    to_grid = (GRID*) malloc(sizeof(GRID));
  }
  if (to_grid == NULL) {
    error_fn("%s: out of memory allocating grid", __func__);
  }

  n_extend = (int) ceil(extend_by/from_grid->spacing);

  for (i=0; i<3; i++) {
    to_grid->limit[i][0] = from_grid->limit[i][0] - n_extend;
    to_grid->limit[i][1] = from_grid->limit[i][1] + n_extend;

    to_grid->flimit[i][0] = from_grid->flimit[i][0] - ((double)n_extend)*from_grid->spacing;
    to_grid->flimit[i][1] = from_grid->flimit[i][1] + ((double)n_extend)*from_grid->spacing;

    to_grid->npoints[i] = from_grid->npoints[i] + n_extend * 2;
  }

  to_grid->padding = from_grid->padding;
  to_grid->spacing = from_grid->spacing;

  return to_grid;
}

void write_probe_orientations (MOLECULE *molecule, ATOM_COORDS **orientations, int n_ori, char *filename_prefix) {
  int i;
  char filename[MAX_LINE_LEN];

  sprintf(filename, "%s.sdf", filename_prefix);
  PLI_FILE *sdf_file = open_file(filename, "w");
  ATOM_COORDS *coords, *initial_coords;
  ATOMLIST *atom_list;
  char molecule_name[50];
  PLI_FILE *file;

  atom_list = molecule->selection;
  initial_coords = molecule2coords(molecule);

  for (i=0; i<n_ori; i++) {
    coords = orientations[i];
    coords2molecule(coords, molecule);
    sprintf(molecule_name, "molecule %d", i+1);
    //sprintf(filename, "%s_%d.pdb", filename_prefix, i+1);
    //file = open_file(filename, "w");
    //write_pdb_atom_list(file, atom_list, OUTPUT_AXES);
    //write_pdb_atom_list(file, atom_list, 0);
    write_mdl_atom_list(sdf_file, atom_list, molecule_name, 0, 0);
    //close_file(file);
  }

  close_file(sdf_file);
  coords2molecule(initial_coords, molecule);
  free(initial_coords);
}

void print_vector(double *v) {
  printf("vector: %.3f, %.3f, %.3f\n", v[0], v[1], v[2]);
}

// add a new scoring function!!
void run_grid_score(SETTINGS *settings) {
  int i, j;
  SYSTEM *system;
  SYSTEM_LIST *syslist;
  MAP *single_map;
  MAP **type_maps;
  MAP *map;

  ATOM_TYPING_SCHEME *scheme;
  MOLECULE_LIST *ligand_list;
  MOLECULE *ligand;
  double score;
  ATOM *atom;
  int gp[3];
  double ***matrix;
  char map_file[MAX_LINE_LEN];

  scheme = get_atom_typing_scheme();
  syslist = settings2syslist(settings);

  single_map = NULL;
  type_maps = NULL;


  if (strcmp(params_get_parameter("gridin")->value.s, "undefined")) {
    single_map = read_insight_map(params_get_parameter("gridin")->value.s);
  } else {
    type_maps = (MAP**) calloc(scheme->n_atom_types, sizeof(MAP*));
    for (i=0; i<scheme->n_atom_types; i++) {
      type_maps[i] = NULL;
    }
  }

  ligand_list = settings->sysmols->ligand_list;

  for (i=0; i<ligand_list->n_molecules; i++) {
    ligand = ligand_list->molecules[i];
    prep_molecule(ligand, settings);
    score = 0.0;
    for (j=0; j<ligand->selection->natoms; j++) {
      atom = ligand->selection->atom[j];
      if ((atom->flags & SKIP_ATOM) || (atom->element->id == HYDROGEN)) {
        continue;
      }
      if (single_map) {
        //pos2grid_round(atom->position, gp, single_map->grid);
        //score += ((double***) single_map->matrix)[gp[0]][gp[1]][gp[2]];
        score += trilinear_interpolate(single_map, atom->position, 1000.0);
      } else {
        map = type_maps[atom->type->i];
        if (map == NULL) {
          (atom->type->id != -1) ? sprintf(map_file, "atom_type_%d.grd", atom->type->id): sprintf(map_file, "element_type_%d.grd", atom->element->id);
          map = read_insight_map(map_file);
        }
        //pos2grid_round(atom->position, gp, map->grid);
        score += trilinear_interpolate(map, atom->position, settings->sfunc->bad_score);
      }
    }
    printf("{\"system\":{\"name\":\"%s\",\"natoms\":%d,\"grid_score\":%.3f}}\n", ligand->name, ligand->selection->natoms, score);
  }
}

