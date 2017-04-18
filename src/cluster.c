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

enum {UNCLUSTERED = -1, CLUSTER_MEMBER = 0, CLUSTER_SEED = 1};

static void cluster_poses(POSE_LIST*, double, int[][2]);
static void update_div_cluster_set(CLUSTER_SET*);
static void set_pose_values(PROBE_POSE*, double, double*, int, MOLECULE*, DOCKED_WATERS*);
static void update_pose_list(POSE_LIST*, int);
static void sort_pose_list(POSE_LIST *pose_list);
static void merge_diverse_clusters (CLUSTER_SET *div_clusters, MOLECULE_PROBE *probe);
static void merge_clusters (POSE_LIST *cluster1, POSE_LIST *cluster2);
void sort_cluster_set(CLUSTER_SET *set);
static double mean_cluster_score(POSE_LIST *list);
double pose_sqr_dist_center(POSE_LIST *cluster, double *center, int min);
double pose_sqr_dist_rms(POSE_LIST *cluster, PROBE_POSE *pose, MOLECULE_PROBE *probe, int min);
double mol_pose_rms(MOLECULE *mol, double **positions, BOOLEAN ***mappings, int n_mappings);
void test_div_cluster_score(CLUSTER_SET *div_cluster_set, SYSTEM *system, ATOM_COORDS **probe_orientations, MOLECULE_PROBE *probe);
static CLUSTER_SETTINGS *cluster_settings = NULL;

CLUSTER_SETTINGS* get_cluster_settings() {
  if (cluster_settings == NULL) {
    init_cluster_settings();
  }
  return cluster_settings;
}

void init_cluster_settings() {

  if (cluster_settings) {
    return;
  }

  cluster_settings = malloc(sizeof(struct ClusterSettings));
  if (cluster_settings == NULL) {
    error_fn("%s: out of memory allocating cluster settings", __func__);
  }

  cluster_settings->save_poses = (params_get_parameter("save_poses"))->value.i;
  cluster_settings->n_saved_poses = (params_get_parameter("n_saved_poses"))->value.i;
  cluster_settings->n_poses_per_orientation = (params_get_parameter("n_poses_per_orientation"))->value.i;
  cluster_settings->n_poses_per_gridpoint = (params_get_parameter("n_poses_per_gridpoint"))->value.i;
  cluster_settings->saved_sdf_file = params_get_parameter("saved_sdf_file")->value.s;
  cluster_settings->saved_data_file = params_get_parameter("saved_data_file")->value.s;

  if (cluster_settings->n_poses_per_gridpoint > 0 && cluster_settings->n_poses_per_orientation > 0) {
    error_fn("%s: only one of the options [n_orientations_per_gridpoint,\
             n_positions_per_orientation] can be > 0\n", __func__);
  }


  if (cluster_settings->n_poses_per_gridpoint > 0) {
    cluster_settings->save_poses = 1;
    if (cluster_settings->n_saved_poses == 0) {
      warning_fn("setting \"n_saved_poses\" to default value (%s)...\n", (params_get_parameter("n_saved_poses"))->default_value);
      sscanf((params_get_parameter("n_saved_poses"))->default_value, "%d", &(cluster_settings->n_saved_poses));
    }
  }


  cluster_settings->save_diverse_clusters = (params_get_parameter("save_diverse_clusters"))->value.i;
  cluster_settings->n_diverse_poses = (params_get_parameter("n_diverse_poses"))->value.i;
  cluster_settings->diverse_cluster_size = (params_get_parameter("diverse_cluster_size"))->value.i;
  cluster_settings->diverse_data_file = params_get_parameter("diverse_data_file")->value.s;
  cluster_settings->diverse_sdf_file_prefix = params_get_parameter("diverse_sdf_file_prefix")->value.s;
  cluster_settings->save_peak_clusters = (params_get_parameter("save_peak_clusters"))->value.i;

  cluster_settings->save_cutoff_atom = params_get_parameter("save_cutoff")->value.d;
  cluster_settings->cluster_cutoff = (params_get_parameter("cluster_cutoff"))->value.d;
  cluster_settings->sqr_cluster_cutoff = sqr((params_get_parameter("cluster_cutoff"))->value.d);
  cluster_settings->dist_type = (strcmp((params_get_parameter("cluster_dist_type"))->value.s, "rmsd") == 0) ? RMSD_DIST: CENTER_DIST;
  cluster_settings->min_dist = params_get_parameter("cluster_min_dist")->value.i;
  cluster_settings->merge_clusters = (params_get_parameter("merge_clusters"))->value.i;
  if (cluster_settings->n_poses_per_gridpoint ||
      cluster_settings->n_poses_per_orientation ||
      cluster_settings->save_diverse_clusters ||
      cluster_settings->save_peak_clusters ||
      cluster_settings->save_poses) {
    cluster_settings->save = 1;
  } else {
    cluster_settings->save = 0;
  }

  cluster_settings->n_peaks = (params_get_parameter("n_peaks"))->value.i;
  cluster_settings->peak_cluster_size = (params_get_parameter("peak_cluster_size"))->value.i;
}


void merge_diverse_clusters (CLUSTER_SET *div_clusters, MOLECULE_PROBE *probe) {
  int i, j, k;
  POSE_LIST *cluster1, *cluster2;
  double dist;
  PROBE_POSE *pose;

  for (i=0; i<(div_clusters->n_alloc_clusters-1); i++) {
    cluster1 = div_clusters->pose_lists[i];
    if ((cluster1 == NULL) || (cluster1->n_poses == 0)) {
      continue;
    }
    for (j=(i+1); j<div_clusters->n_alloc_clusters; j++) {
      cluster2 = div_clusters->pose_lists[j];
      if ((cluster2 == NULL) || (cluster2->n_poses == 0)) {
        continue;
      }
      for (k=0; k<cluster2->n_poses; k++) {
        pose = cluster2->poses[k];
        if (cluster_settings->dist_type == RMSD_DIST) {
          dist = pose_sqr_dist_rms(cluster1, pose, probe, 1);
        }
        else if (cluster_settings->dist_type == CENTER_DIST) {
          dist = pose_sqr_dist_center(cluster1, pose->center, 1);
        } else {
          error_fn("%s: Distance type %d not implemented", __func__, cluster_settings->dist_type);
        }
        if (dist < cluster_settings->cluster_cutoff) {
          merge_clusters(cluster1, cluster2);
          break;
        }
      }
    }
  }
}


void merge_clusters (POSE_LIST *cluster1, POSE_LIST *cluster2) {
  PROBE_POSE *pose1, *pose2;

  int i = 0;
  int cluster1_loc = 0;
  int cluster2_loc = 0;


  POSE_LIST *merged_clster = alloc_pose_list(cluster1->n_poses + cluster2->n_poses);

  while ((cluster1_loc < cluster1->n_poses) || (cluster2_loc < cluster2->n_poses)) {
    if (cluster1_loc == cluster1->n_poses) {
      merged_clster->poses[i] = cluster2->poses[cluster2_loc];
      cluster2_loc += 1;
      i+=1;
      continue;
    }
    else if (cluster2_loc == cluster2->n_poses) {
      merged_clster->poses[i] = cluster1->poses[cluster1_loc];
      cluster1_loc += 1;
      i+=1;
      continue;
    }
    else if (cluster1->poses[cluster1_loc]->score <= cluster2->poses[cluster2_loc]->score) {
      merged_clster->poses[i] = cluster1->poses[cluster1_loc];
      cluster1_loc += 1;
      i+=1;
      continue;
    }
    else {
      merged_clster->poses[i] = cluster2->poses[cluster2_loc];
      cluster2_loc += 1;
      i+=1;
      continue;
    }
  }

  // copy top poses to cluster 1
  for (i=0; i<(cluster1->n_poses + cluster2->n_poses); i++) {
    if (i==cluster1->max_n) {
      break;
    }
    cluster1->poses[i] = merged_clster->poses[i];
  }
  cluster1->n_poses = i;
  cluster1->best_score_pose = cluster1->poses[0];
  cluster1->worst_score_pose = cluster1->poses[i-1];

  cluster2->n_poses = 0;

  free_pose_list(merged_clster);
  free_pose_list(cluster2);
  // cluster2 = NULL;
}

void test_div_cluster_score(CLUSTER_SET *div_cluster_set, SYSTEM *system, ATOM_COORDS **probe_orientations, MOLECULE_PROBE *probe) {
  int i, j, k, l;
  PROBE_POSE *pose;
  POSE_LIST *cluster;
  ATOMLIST *selection;
  ATOM *atom;
  char coord[10];
  MOLECULE *mol;
  ATOM_COORDS *coords;
  double old_pose_recalc, new_pose_recalc;

  selection = system->selection;
  system->selection = probe->molecule->selection;

  char filename_before[MAX_LINE_LEN];
  char filename_after[MAX_LINE_LEN];
  PLI_FILE *file_before;
  PLI_FILE *file_after;


  for (i=0; i<div_cluster_set->n_alloc_clusters; i++) {
    cluster = div_cluster_set->pose_lists[i];
    if (cluster->n_poses == 0) {
      continue;
    }
    for (j=0; j<cluster->n_poses; j++) {
      //printf("i: %d, j: %d\n", i, j);
      pose = cluster->poses[j];
      coords2molecule(probe_orientations[pose->orientation_index], probe->molecule);
      null_vector(probe->molecule->geometric_center);
      move_molecule_centre(probe->molecule, pose->center);
      copy_vector(pose->center, probe->molecule->geometric_center);

      reset_system(system, 0);
      score_system(system, system->settings->sfunc);
      old_pose_recalc = system->score;

      sprintf(filename_before, "probe_i_%d_j_%d_before.sdf", i, j);
      file_before = open_file(filename_before, "w");
      write_mdl_atom_list(file_before, probe->molecule->selection, "mol", 0, 0);
      close_file(file_before);
      for (k=0; k<probe->molecule->natoms; k++) {
        atom = probe->molecule->atom + k;
        for (l=0; l<3; l++) {
          //printf("coords: %.10f", atom->position[l]);
          sprintf(coord, "%.4f", atom->position[l]);
          //printf("\t%s", coord);
          sscanf(coord, "%lf", &(atom->position[l]));
          //printf("\t%.10f\n", atom->position[l]);
        }
      }
      //shake_molecule(probe->molecule);
      sprintf(filename_after, "probe_i_%d_j_%d_after.sdf", i, j);
      file_after = open_file(filename_after, "w");
      write_mdl_atom_list(file_after, probe->molecule->selection, "mol", 0, 0);
      close_file(file_after);


      reset_system(system, 0);
      score_system(system, system->settings->sfunc);
      new_pose_recalc = system->score;

      printf("pose_score: %.3f, old_pose_recalc: %.3f; new_pose_recalc %.3f, diff1: %.3f, diff2: %.3f\n",
             pose->score, old_pose_recalc, new_pose_recalc, pose->score - old_pose_recalc, pose->score - new_pose_recalc);
      printf("\n\n\n\n");

    }
  }
  system->selection = selection;
}


void sort_cluster_set(CLUSTER_SET *set) {
  for (int i=0; i<(set->n_alloc_clusters-1); i++) {
    if (set->pose_lists[i]->n_poses == 0) {
      continue;
    }
    int position = i;

    for (int j=i+1; j<set->n_alloc_clusters; j++) {
      if (set->pose_lists[j]->n_poses == 0) {
        continue;
      }
      if (set->pose_lists[position]->best_score_pose->score > set->pose_lists[j]->best_score_pose->score) {
        position = j;
      }
    }
    if (position != i) {
      POSE_LIST *temp_cluster = set->pose_lists[i];
      set->pose_lists[i] = set->pose_lists[position];
      set->pose_lists[position] = temp_cluster;
    }
  }
}

void write_cluster_set(CLUSTER_SET *div_cluster_set, char *prefix, MOLECULE_PROBE *probe, MOLECULE *water_template) {
  int i, j, k, cluster_id, molecule_id, i_water;
  PLI_FILE *probe_cluster_file, *data_file;
  POSE_LIST *cluster;
  ATOM_COORDS *coords, *initial_coords;
  PROBE_POSE *pose;
  ATOM *atom;
  char mol_name[MAX_LINE_LEN];
  MOLECULE_PROPERTY mp_score, mp_cluster_id, mp_cluster_attr;
  char cluster_attr[20];

  char cluster_file_name[MAX_LINE_LEN];
  char data_file_name[MAX_LINE_LEN];
  int write_hydrogens;

  char water_file_name[MAX_LINE_LEN];
  PLI_FILE *water_cluster_file;

  initial_coords = molecule2coords(probe->molecule);

  sort_cluster_set(div_cluster_set);

  if (cluster_settings->merge_clusters) {
    merge_diverse_clusters(div_cluster_set, probe);
  }


  // write the poses now
  sprintf(data_file_name, "%ss.dat", prefix);
  data_file = open_file(data_file_name, "w");

  write_line(data_file, "%14s%14s%14s%14s%14s\n", "cluster_ID", "n_molecules", "min_score", "max_score", "mean_score");

  cluster_id = 1;

  for (i=0; i<div_cluster_set->n_alloc_clusters; i++) {
    molecule_id = 1;
    cluster = div_cluster_set->pose_lists[i];

    if ((cluster == NULL) || (cluster->n_poses == 0)) {
      continue;
    }
    sort_pose_list(cluster);
    sprintf(cluster_file_name, "%s_%d.sdf", prefix, cluster_id);
    probe_cluster_file = open_file(cluster_file_name, "w");

    if (water_template) {
      water_cluster_file = NULL;
      sprintf(water_file_name, "%s_%d_waters.sdf", prefix, cluster_id);
    }

    for (j=0; j<cluster->n_poses; j++) {
      pose = cluster->poses[j];

      // orient molecule
      if (probe->orientations != NULL) {
        coords = probe->orientations[pose->orientation_index];
        coords2molecule(coords, probe->molecule);
        null_vector(probe->molecule->geometric_center);
        move_molecule_centre(probe->molecule, pose->center);
        copy_vector(pose->center, probe->molecule->geometric_center);
        write_hydrogens = 1;
      } else {
        for (k=0; k<probe->molecule->selection->natoms; k++) {
          copy_vector(pose->position[k], probe->molecule->selection->atom[k]->position);
        }
        write_hydrogens = 0;
      }

      sprintf(mol_name, "cluster_%d_molecule_%d", cluster_id, molecule_id);

      strcpy(mp_score.name, "Score");
      mp_score.type='d';
      mp_score.value.d = pose->score;

      strcpy(mp_cluster_id.name, "Cluster ID");
      mp_cluster_id.type = 'i';
      mp_cluster_id.value.i = cluster_id;

      strcpy(mp_cluster_attr.name, "Cluster attribute");
      mp_cluster_attr.type = 's';
      j == 0 ? strcpy(cluster_attr, "Seed"): strcpy(cluster_attr, "Member");
      strcpy(mp_cluster_attr.value.s, cluster_attr);

      write_mdl_atom_list(probe_cluster_file, probe->molecule->selection, mol_name, write_hydrogens, 3, &mp_cluster_id, &mp_cluster_attr, &mp_score);


      // save bound waters
      if (water_template) {
        // open file for writing if not already opened
        if (pose->n_added_waters) {
          mp_score.value.d = pose->water_score;
          if (! water_cluster_file) {
            water_cluster_file = open_file(water_file_name, "w");
          }

          for (i_water=0; i_water<pose->n_added_waters; i_water++) {
            copy_vector(pose->water_positions[i_water], water_template->atom->position);
            write_mdl_atom_list(water_cluster_file, water_template->selection, mol_name, 0, 3, &mp_cluster_id, &mp_cluster_attr, &mp_score);
          }
        }
      }
      molecule_id++;
    }
    close_file(probe_cluster_file);

    if (water_template && water_cluster_file) {
      close_file(water_cluster_file);
    }

    write_line(data_file, "%14d%14d%14.3f%14.3f%14.3f\n",
                            cluster_id,
                            cluster->n_poses,
                            cluster->best_score_pose->score,
                            cluster->worst_score_pose->score,
                            mean_cluster_score(cluster));
    cluster_id++;
  }

  coords2molecule(initial_coords, probe->molecule);
  free(initial_coords);
  close_file(data_file);
}


double mean_cluster_score(POSE_LIST *list) {
  if (list->n_poses == 0) {
    return 0.0;
  }
  double sum_score = list->poses[0]->score;
  for (int i =1; i<list->n_poses; i++) {
    sum_score += list->poses[i]->score;
  }
  return sum_score/(double) list->n_poses;
}

POSE_LIST* alloc_pose_list(int n_poses) {

  if (n_poses < 1) {
    return NULL;
  }

  POSE_LIST *p_list;

  p_list = (POSE_LIST *) malloc(sizeof(POSE_LIST));

  if (p_list==NULL) {
    error_fn("%s: error allocating probe_position matrix\n", __func__);
  }

  p_list->max_n = n_poses;
  p_list->n_poses = 0;
  p_list->worst_score_pose = NULL;
  p_list->best_score_pose = NULL;
  p_list->poses = calloc(n_poses, sizeof(PROBE_POSE *));

  if (p_list->poses == NULL) {
      error_fn("%s: error allocating poses", __func__);
  }

  return p_list;
}

void free_pose_list(POSE_LIST *pose_list) {
  int i;
  if (pose_list) {
    if (pose_list->poses) {
      free(pose_list->poses);
    }
    free(pose_list);
  }
}

PROBE_POSE* alloc_poses(int n_poses, int n_atoms, int n_waters) {
  int i;
  PROBE_POSE *poses;

  poses = calloc(n_poses, sizeof(PROBE_POSE));

  for (i=0; i<n_poses; i++) {
    (poses + i)->n_added_waters = 0;

    if (n_atoms) {
      (poses + i)->position = alloc_2d_fmatrix(n_atoms, 4);
    } else {
      (poses + i)->position = NULL;
    }

    if (n_waters) {
      (poses + i)->water_positions = alloc_2d_fmatrix(n_waters, 4);
    } else {
      (poses + i)->water_positions = NULL;
    }

  }

  return(poses);
}

void free_poses(PROBE_POSE *allocated_poses) {
  if (allocated_poses) {
    if (allocated_poses->position) {
      free_2d_fmatrix(allocated_poses->position);
    }
    free(allocated_poses);
  }
}



CLUSTER_SET* alloc_div_cluster_set(int n_max_poses, int n_atoms, int n_waters){

  if (n_max_poses < 1) {
    return NULL;
  }

  CLUSTER_SET *cluster_set;
  cluster_set = malloc(sizeof(CLUSTER_SET));
  if (cluster_set == NULL) {
    error_fn("%s: error allocating clusters");
  }
  cluster_set->pose_lists = NULL;
  cluster_set->n_alloc_clusters = 0;
  cluster_set->n_max_poses = n_max_poses;
  cluster_set->n_poses = 0;
  cluster_set->worst_score_cluster = NULL;
  cluster_set->allocated_poses = alloc_poses(n_max_poses, n_atoms, n_waters);
  return cluster_set;
}

void free_div_cluster_set(CLUSTER_SET *div_cluster_set) {
  int i;
  POSE_LIST *cluster;

  if (div_cluster_set) {
    for (i=0; i<div_cluster_set->n_alloc_clusters; i++) {
      cluster = div_cluster_set->pose_lists[i];
      if (cluster) {
        free_pose_list(cluster);
      }
    }
    if (div_cluster_set->allocated_poses) {
      free_poses(div_cluster_set->allocated_poses);
    }
    free(div_cluster_set);
  }
}

void sort_pose_list(POSE_LIST *pose_list) {
  int i, j, position;
  PROBE_POSE *temp_pose;
  for (i=0; i<pose_list->n_poses-1; i++) {
    position = i;
    for (j=i+1; j<pose_list->n_poses; j++) {
      if (pose_list->poses[position]->score > pose_list->poses[j]->score) {
        position = j;
      }
    }
    if (position != i) {
      temp_pose = pose_list->poses[i];
      pose_list->poses[i] = pose_list->poses[position];
      pose_list->poses[position] = temp_pose;
    }
  }
}


void save_diverse_pose(CLUSTER_SET *div_cluster_set, MOLECULE_PROBE *probe, DOCKED_WATERS *docked_waters) {

  int i, n_max_poses_per_cluster;
  double sqr_dist, cur_sqr_dist;

  POSE_LIST *acceptor_cluster, *cluster, *empty_cluster, *donor_cluster;
  PROBE_POSE *pose;

  n_max_poses_per_cluster = cluster_settings->diverse_cluster_size;


  double score = probe->score;
  double *center = probe->molecule->geometric_center;
  MOLECULE *molecule = probe->molecule;
  int ori_ind = probe->orientation_id;

  // return if the div_cluster is full and score is worse than the overall worst score
  if ((div_cluster_set->n_poses == div_cluster_set->n_max_poses)
      && (score > div_cluster_set->worst_score_cluster->worst_score_pose->score)) {
    return;
  }

  // find the cluster to add the current pose to
  sqr_dist = cluster_settings->sqr_cluster_cutoff + 10.0;
  empty_cluster = NULL;
  acceptor_cluster = NULL;   // the cluster where the current pose goes to
  donor_cluster = NULL;      // the cluster which loses a pose to make space for the current pose

  // find acceptor cluster and empty cluster

  for (i=0; i<div_cluster_set->n_alloc_clusters; i++) {
    cluster = div_cluster_set->pose_lists[i];

    if (cluster->n_poses == 0) {
      empty_cluster = cluster;
      continue;
    }

    if (cluster_settings->dist_type == RMSD_DIST) {
      // check the center dist first
      cur_sqr_dist = pose_sqr_dist_center(cluster, probe->molecule->geometric_center, 0);
      cur_sqr_dist = 0.0;
      if (cur_sqr_dist < cluster_settings->sqr_cluster_cutoff) {
        cur_sqr_dist = pose_sqr_dist_rms(cluster, NULL, probe, cluster_settings->min_dist);
      }
    } else if (cluster_settings->dist_type == CENTER_DIST) {
      cur_sqr_dist = pose_sqr_dist_center(cluster, probe->molecule->geometric_center, cluster_settings->min_dist);
    } else {
      error_fn("%s: Distance type %d not implemented", __func__, cluster_settings->dist_type);
    }

    if (cur_sqr_dist < sqr_dist) {
      sqr_dist = cur_sqr_dist;
      acceptor_cluster = cluster;
    }
  }

  if (sqr_dist > cluster_settings->sqr_cluster_cutoff) {
    acceptor_cluster = NULL;
  }

  if (acceptor_cluster != NULL) {
    if (acceptor_cluster->n_poses == acceptor_cluster->max_n && score>acceptor_cluster->worst_score_pose->score) {
      return;
    }
  }

  // if acceptor cluster is not found create a new cluster
  if (acceptor_cluster == NULL) {
    if (empty_cluster != NULL) {
      acceptor_cluster = empty_cluster;
    } else {
      div_cluster_set->pose_lists = realloc(div_cluster_set->pose_lists, (div_cluster_set->n_alloc_clusters+1)*sizeof(POSE_LIST*));
      div_cluster_set->n_alloc_clusters++;
      div_cluster_set->pose_lists[div_cluster_set->n_alloc_clusters-1] = alloc_pose_list(n_max_poses_per_cluster);
      update_div_cluster_set(div_cluster_set);
      acceptor_cluster = div_cluster_set->pose_lists[div_cluster_set->n_alloc_clusters-1];
    }
  }

  // we have an acceptor cluster now

  // the acceptor cluster is already full then replace worst pose of the cluster
  if (acceptor_cluster->n_poses == acceptor_cluster->max_n) {
    pose = acceptor_cluster->worst_score_pose;
    donor_cluster = acceptor_cluster;

  // acceptor cluster is not full
  } else {
    // div_cluster_set is full then replace the global worst pose
    if (div_cluster_set->n_poses == div_cluster_set->n_max_poses) {
      pose = div_cluster_set->worst_score_cluster->worst_score_pose;
      if (div_cluster_set->worst_score_cluster == acceptor_cluster) {
        donor_cluster = acceptor_cluster;
      } else {
        acceptor_cluster->poses[acceptor_cluster->n_poses] = pose;
        donor_cluster = div_cluster_set->worst_score_cluster;
      }

    // div_cluster_set is not full yet create a new pose
    } else {
      pose = div_cluster_set->allocated_poses + div_cluster_set->n_poses;
      acceptor_cluster->poses[acceptor_cluster->n_poses] = pose;
      donor_cluster = NULL;
    }
  }


  // now we have pose, acceptor cluster and the donor cluster
  set_pose_values(pose, probe->score, probe->molecule->geometric_center, probe->orientation_id, probe->molecule, docked_waters);

  if (acceptor_cluster != donor_cluster) {
    acceptor_cluster->n_poses++;
    update_pose_list(acceptor_cluster, 0);
    if (donor_cluster != NULL) {
      update_pose_list(donor_cluster, 1);
    }
  } else {
    update_pose_list(acceptor_cluster, 0);
  }

  if (donor_cluster == NULL) {
    div_cluster_set->n_poses++;
  }
  update_div_cluster_set(div_cluster_set);
}



void update_pose_list(POSE_LIST *pose_list, int remove_worst_pose) {
  int i;
  double worst_score, best_score;
  worst_score = -1E10;
  best_score = 1E10;


  if (remove_worst_pose) {
    for (i=0; i<pose_list->n_poses; i++) {
       if (pose_list->poses[i] == pose_list->worst_score_pose) {
         pose_list->poses[i] = pose_list->poses[pose_list->n_poses-1];
         break;
       }
    }
    pose_list->n_poses--;
  }

  pose_list->worst_score_pose = NULL;
  pose_list->best_score_pose = NULL;

  // TODO: this can be done without all comparisons for best score and for worst score when worst score pose is not removed
  for (i=0; i<pose_list->n_poses; i++) {
    if (pose_list->poses[i]->score > worst_score) {
      pose_list->worst_score_pose = pose_list->poses[i];
      worst_score = pose_list->poses[i]->score;
    }
    if (pose_list->poses[i]->score < best_score) {
      pose_list->best_score_pose = pose_list->poses[i];
      best_score = pose_list->poses[i]->score;
    }
  }
}


void update_div_cluster_set(CLUSTER_SET *div_cluster_set) {
  int i;
  double worst_score;
  POSE_LIST *cluster;

  worst_score = -1E10;

  for (i=0; i<div_cluster_set->n_alloc_clusters; i++) {
    cluster = div_cluster_set->pose_lists[i];
    if (cluster->n_poses == 0) {
      continue;
    }

    if (cluster->worst_score_pose->score > worst_score) {
      div_cluster_set->worst_score_cluster = cluster;
      worst_score = cluster->worst_score_pose->score;
    }
  }
}


void cluster_poses(POSE_LIST *pose_list, double cutoff, int cluster_output[][2]) {
  int i, j, cluster_id;
  PROBE_POSE *pose1, *pose2;

  cluster_id = 1;

  for (i=0; i<pose_list->n_poses; i++) {
    cluster_output[i][1] = UNCLUSTERED;
  }

  for (i=0; i<pose_list->n_poses; i++) {
    pose1 = pose_list->poses[i];
    if (cluster_output[i][1] == UNCLUSTERED) {
      cluster_output[i][1] = CLUSTER_SEED;
      cluster_output[i][0] = cluster_id;
      if (i < (pose_list->n_poses-1)) {
        for (j=i+1; j<pose_list->n_poses; j++) {
          pose2 = pose_list->poses[j];
          if (distance(pose1->center, pose2->center) < cutoff) {
            cluster_output[j][1] = CLUSTER_MEMBER;
            cluster_output[j][0] = cluster_id;
          }
        }
      }
      cluster_id++;
    }
  }
}


double pose_sqr_dist_center(POSE_LIST *cluster, double *center, int min) {
  double min_sqr_dist = sqr_distance(cluster->best_score_pose->center, center);

  if (! min) {
    return min_sqr_dist;
  }

  for (int i=0; i<cluster->n_poses; i++) {
    if (cluster->poses[i] == cluster->best_score_pose) {
      continue;
    }
    double sqr_dist = sqr_distance(cluster->poses[i]->center, center);
    if (sqr_dist < min_sqr_dist) {
      min_sqr_dist = sqr_dist;
    }
  }

  return min_sqr_dist;
}

double pose_sqr_dist_rms(POSE_LIST *cluster, PROBE_POSE *pose, MOLECULE_PROBE *probe, int min) {
  // use orientation index, center, positions from pose if provided, otherwise use molecule for these
  int pose_ori_id = (pose == NULL) ? probe->orientation_id: pose->orientation_index;
  double *pose_center = (pose == NULL) ? probe->molecule->geometric_center: pose->center;

  if (pose_ori_id != -1) {
    double min_rmsd = rms_coords(probe->orientations[pose_ori_id],
                                 probe->orientations[cluster->best_score_pose->orientation_index],
                                 probe->selected_ids,
                                 probe->selected_ids,
                                 probe->molecule->selection->natoms,
                                 probe->auto_mappings,
                                 probe->n_mappings,
                                 pose_center,
                                 cluster->best_score_pose->center);
    if (!min) {
      return sqr(min_rmsd);
    }

    for (int i=0; i<cluster->n_poses; i++) {
      if (cluster->poses[i] == cluster->best_score_pose) {
        continue;
      }
      double rmsd = rms_coords(probe->orientations[pose_ori_id],
                               probe->orientations[cluster->poses[i]->orientation_index],
                               probe->selected_ids,
                               probe->selected_ids,
                               probe->molecule->selection->natoms,
                               probe->auto_mappings,
                               probe->n_mappings,
                               pose_center,
                               cluster->poses[i]->center);

      if (rmsd < min_rmsd) {
        min_rmsd = rmsd;
      }
    }
    return sqr(min_rmsd);
  } else {
    // rmsd with molecule atom positions
    double min_rmsd = mol_pose_rms(probe->molecule, cluster->best_score_pose->position, probe->auto_mappings, probe->n_mappings);
    if (!min) {
      return sqr(min_rmsd);
    }

    for (int i=0; i<cluster->n_poses; i++) {
      if (cluster->poses[i] == cluster->best_score_pose) {
        continue;
      }
      double rmsd = mol_pose_rms(probe->molecule, cluster->poses[i]->position, probe->auto_mappings, probe->n_mappings);
      if (rmsd < min_rmsd) {
        min_rmsd = rmsd;
      }
    }
    return sqr(min_rmsd);
  }
}


double mol_pose_rms(MOLECULE *mol, double **positions, BOOLEAN ***mappings, int n_mappings) {
  double min_rms;
  double rms;

  min_rms = 1E100;
  int i,j,k;

  for (i=0; i<n_mappings; i++) {
    rms = 0.0;
    for (j=0; j<mol->selection->natoms; j++) {
      for (k=0; k<mol->selection->natoms; k++) {
        if (mappings[i][j][k] == 1) {
          rms += sqr_distance(mol->selection->atom[j]->position, positions[k]);
        }
      }
    }
    rms = sqrt(rms/(double) mol->selection->natoms);
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



void write_pose_list(POSE_LIST *saved_pose_list, MOLECULE *molecule, ATOM_COORDS **probe_orientations) {
  int cluster_output[saved_pose_list->n_poses][2];

  ATOMLIST *selection;
  int i,j;
  double new_score;
  ATOM *atom;
  PLI_FILE *sdf_file = NULL;
  PLI_FILE *score_file;
  PROBE_POSE *pose;
  ATOMLIST *atom_list;
  ATOM_COORDS *coords, *initial_coords;
  double cluster_cutoff;
  char *sdf_file_name, *data_file_name;
  int write_hydrogens;

  sdf_file_name = cluster_settings->saved_sdf_file;
  data_file_name = cluster_settings->saved_data_file;
  cluster_cutoff = cluster_settings->cluster_cutoff;


  char mol_name[MAX_LINE_LEN];
  MOLECULE_PROPERTY mp_score, mp_cluster_attr, mp_cluster_id;

  score_file = open_file(data_file_name, "w");
  sdf_file = open_file(sdf_file_name, "w");

  initial_coords = molecule2coords(molecule);
  sort_pose_list(saved_pose_list);
  cluster_poses(saved_pose_list, cluster_cutoff, cluster_output);


  /*if (rescore) {
    reset_system(system,0);
    selection = system->selection;
    add_molecule_to_list(system->molecule_list,molecule);
    system->selection = molecule->selection;
  }
  */

  for (i=0; i<saved_pose_list->n_poses; i++) {
    pose = saved_pose_list->poses[i];

    // orient molecule
    if (probe_orientations != NULL) {
      coords = probe_orientations[pose->orientation_index];
      coords2molecule(coords, molecule);
      null_vector(molecule->geometric_center);
      move_molecule_centre(molecule, pose->center);
      copy_vector(pose->center, molecule->geometric_center);
      write_hydrogens = 1;
    } else {
      for (j=0; j<molecule->selection->natoms; j++) {
        copy_vector(pose->position[j], molecule->selection->atom[j]->position);
      }
      write_hydrogens = 0;
    }

    write_line(score_file, "score %.2f \t Cluster_ID: %d \t Seed/Member: %d\n",
               pose->score, cluster_output[i][0], cluster_output[i][1]);

    // write sdf file
    sprintf(mol_name, "%s_%d", molecule->name, i+1);

    strcpy(mp_score.name, "Score"), mp_score.type='d', mp_score.value.d = pose->score;
    strcpy(mp_cluster_id.name, "Cluster ID"), mp_cluster_id.type = 'i', mp_cluster_id.value.i = cluster_output[i][0];
    strcpy(mp_cluster_attr.name, "Cluster attribute"), mp_cluster_attr.type = 's';
    (cluster_output[i][1] == CLUSTER_SEED) ? strcpy(mp_cluster_attr.value.s, "Seed"): strcpy(mp_cluster_attr.value.s, "Member");

    write_mdl_atom_list(sdf_file, molecule->selection, mol_name, write_hydrogens, 3, &mp_cluster_id, &mp_cluster_attr, &mp_score);
  }
  close_file(sdf_file);
  close_file(score_file);
  coords2molecule(initial_coords, molecule);
  free(initial_coords);
/*
  if (rescore) {
    // reset system
    system->selection = selection;
    remove_molecule_from_list(system->molecule_list, probe->molecule);
    reset_system(system,0);
*/
}


void save_pose(POSE_LIST *pose_list, PROBE_POSE *allocated_poses, MOLECULE_PROBE *probe) {
  int i, j, n;
  PROBE_POSE *pose_alloc_p = NULL;
  ATOM *atom;
  double worst_score;

  // find where to allocate the current pose
  if (pose_list->n_poses < pose_list->max_n) {
    pose_alloc_p = allocated_poses + pose_list->n_poses;
    pose_list->poses[pose_list->n_poses] = pose_alloc_p;
    pose_list->n_poses++;

  } else if (probe->score < pose_list->worst_score_pose->score) {
    pose_alloc_p = pose_list->worst_score_pose;
  }

  if (pose_alloc_p != NULL) {
    set_pose_values(pose_alloc_p, probe->score, probe->molecule->geometric_center, probe->orientation_id, probe->molecule, NULL);
    update_pose_list(pose_list, 0);
  }
}

void set_pose_values(PROBE_POSE *pose, double score, double *center, int orientation_id, MOLECULE *molecule, DOCKED_WATERS *docked_waters) {
  ATOM *atom;

  copy_vector(center, pose->center);
  pose->orientation_index = orientation_id;
  pose->score = score;

  if(orientation_id == -1) {
    for (int i=0; i<molecule->selection->natoms; i++) {
      copy_vector(molecule->selection->atom[i]->position, pose->position[i]);
    }
  }
  if (docked_waters) {
    pose->n_added_waters = docked_waters->n_waters;
    pose->water_score = docked_waters->score;
    for (int i=0; i<docked_waters->n_waters; i++) {
      copy_vector((docked_waters->mols+i)->atom->position, pose->water_positions[i]);
    }
  }
}


void add_pose_to_poselist(POSELIST_LITE *list, int orientation_index, double score) {

  if (list->n_alloc_poses == 0) {

    list->n_alloc_poses = 10;

    list->poses = (PROBE_POSE_LITE*) calloc(list->n_alloc_poses,sizeof(PROBE_POSE_LITE));

  } else if (list->n_poses == list->n_alloc_poses) {

    list->n_alloc_poses *= 2;

    list->poses = (PROBE_POSE_LITE*) realloc(list->poses,(list->n_alloc_poses)*sizeof(PROBE_POSE_LITE));

  }

  if (list->poses == NULL) {

    error_fn("%s: out of memory allocating poses", __func__);
  }

  (list->poses + list->n_poses)->orientation_index = orientation_index;
  (list->poses + list->n_poses)->score = score;

  list->n_poses++;
}


CLUSTER_SET* get_peak_clusters(MAP *poselist_map, MAP *score_map, MOLECULE_PROBE *probe, ATOM_COORDS **probe_orientations, CLUSTER_SETTINGS *cluster_settings, PROBE_POSE *cluster_poses, BOOLEAN ***automorphs, int n_automorphs) {
  printf("getting peak clusters\n");
  int i, j, k, min_index, pose_index, index, nx, ny, nz, nchar, gp[3], min_gp[3], n_gp_window, ix, iy, iz;
  GRID *grid;
  double pos[4], min_pos[4], current_pose_center[4];
  POSELIST_LITE *poselist;
  POSELIST_LITE ***poselist_matrix;
  int n_clusters;
  int cluster_size;
  int n_poses_saved;
  CLUSTER_SET *cluster_set;
  int selected_indices[probe->molecule->selection->natoms];
  double min_score;
  double peak_cluster_diversity = params_get_parameter("peak_cluster_diversity")->value.d;

  PROBE_POSE *cluster_ref_pose;

  MOLECULE *mol;
  ATOM_COORDS *coords;


  POSE_LIST *cluster;
  PROBE_POSE *allocated_pose;
  unsigned char *mask_backup, *mask;
  int cluster_window[3][2];
  double rms;
  double grid_point_dist;
  PROBE_POSE_LITE* current_pose;


  n_clusters = cluster_settings->n_peaks;
  cluster_size = cluster_settings->peak_cluster_size;
  n_poses_saved = 0;

  get_selected_atom_indices(probe->molecule, selected_indices);

  cluster_set = alloc_div_cluster_set(n_clusters * cluster_size, 0, 0);
  cluster_set->pose_lists = (POSE_LIST**) calloc(n_clusters, sizeof(POSE_LIST*));
  cluster_set->n_alloc_clusters = n_clusters;
  for (i=0; i<cluster_set->n_alloc_clusters; i++) {

    // add a function
    cluster_set->pose_lists[i] = (POSE_LIST *) malloc(sizeof(POSE_LIST));
    cluster_set->pose_lists[i]->max_n = cluster_size;
    cluster_set->pose_lists[i]->n_poses = 0;
    cluster_set->pose_lists[i]->worst_score_pose = NULL;
    cluster_set->pose_lists[i]->best_score_pose = NULL;
    cluster_set->pose_lists[i]->poses = (PROBE_POSE**) calloc(cluster_size, sizeof(PROBE_POSE*));
    for (int j=0; j<cluster_size; j++) {
      cluster_set->pose_lists[i]->poses[j] = cluster_set->allocated_poses + (i*cluster_size) + j;
    }
  }


  grid = score_map->grid;
  mask = score_map->mask;
  nx = grid->npoints[0];
  ny = grid->npoints[1];
  nz = grid->npoints[2];

  poselist_matrix = (POSELIST_LITE ***) poselist_map->matrix;

  // copy mask
  nchar = (nx*ny*nz)/8 + 1;
  mask_backup = (unsigned char*) calloc(nchar, sizeof(unsigned char));
  memcpy(mask_backup, score_map->mask, nchar * sizeof(unsigned char));

  n_gp_window = (int) floor(cluster_settings->cluster_cutoff/grid->spacing);


  // for each cluster
  for (i=0; i<n_clusters; i++) {
    cluster = cluster_set->pose_lists[i];
    cluster->n_poses = 0;

    min_index = map_extreme_index(score_map, 1, "min", min_gp, min_gp+1, min_gp+2);

    // break the loop if no point in the grid is found
    if (min_index == -1) {
      break;
    }

    // get position of the minimum point
    for (j=0; j<3; j++) {
      min_pos[j] = grid->flimit[j][0] + ((double) min_gp[j])*(grid->spacing);
    }

    // get grid bounds around the minimum point for clustering
    for (j=0; j<3; j++) {
      cluster_window[j][0] = ((min_gp[j]-n_gp_window) >= 0) ? min_gp[j]-n_gp_window : 0;
      cluster_window[j][1] = ((min_gp[j]+n_gp_window) < (grid->npoints[j]-1)) ? min_gp[j]+n_gp_window : grid->npoints[j]-1;
    }

    // for each pose in the cluster
    for (j=0; j<cluster_size; j++) {
      current_pose = NULL;
      min_score = 10000.0;
      for (ix=cluster_window[0][0]; ix<=cluster_window[0][1]; ix++) {
        pos[0] = grid->flimit[0][0] + ((double) ix)*(grid->spacing);
        for (iy=cluster_window[1][0]; iy<=cluster_window[1][1]; iy++) {
          pos[1] = grid->flimit[1][0] + ((double) iy)*(grid->spacing);
          for (iz=cluster_window[2][0]; iz<=cluster_window[2][1]; iz++) {
            index = ny*nz*ix + nz*iy + iz;

            if (is_mask_bit_set(mask,index)) {
              continue;
            }

            pos[2] = grid->flimit[2][0] + ((double) iz)*(grid->spacing);

            grid_point_dist = distance(pos, min_pos);
            if (grid_point_dist > cluster_settings->cluster_cutoff) {
              continue;
            }

            poselist = &(poselist_matrix[ix][iy][iz]);

            // for each pose in the poselist
            for (k=0; k<poselist->n_poses; k++) {
              if ((poselist->poses + k)->score > min_score) {
                continue;
              }
              if ((peak_cluster_diversity > 0.1) && (j>0)) {
                // for the second pose in the cluster onwards, check rmsd
                for (int l=0; l<cluster->n_poses; l++) {
                  rms = rms_coords(probe_orientations[cluster->poses[l]->orientation_index],
                                   probe_orientations[(poselist->poses + k)->orientation_index],
                                   selected_indices, selected_indices, probe->molecule->selection->natoms,
                                   automorphs, n_automorphs, cluster->poses[l]->center, pos);
                  if (rms < peak_cluster_diversity){
                    break;
                  }
                }
                if (rms < peak_cluster_diversity) {
                  continue;
                }
              }
              current_pose = poselist->poses + k;
              current_pose_center[0] = pos[0]; current_pose_center[1] = pos[1]; current_pose_center[2] = pos[2];
              min_score = current_pose->score;
            } // for each pose in the pose_list
          } // iz
        } // iy
      } // ix
      if (current_pose == NULL) {
        break;
      }
      //
      allocated_pose = cluster_poses+n_poses_saved;
      if (j == 0) {
        cluster->best_score_pose = allocated_pose;
      }
      set_pose_values(allocated_pose, current_pose->score, current_pose_center, current_pose->orientation_index, NULL, NULL);
      cluster->poses[j] = allocated_pose;
      cluster->n_poses++;
      cluster->worst_score_pose = allocated_pose;
      n_poses_saved++;
      // this pose should not be found in the next iteration as min pose
      current_pose->score = 9.99E10;
    } // for each pose in the cluster

    // mask the gp from the current min point
    for (ix=cluster_window[0][0]; ix<=cluster_window[0][1]; ix++) {
      for (iy=cluster_window[1][0]; iy<=cluster_window[1][1]; iy++) {
        for (iz=cluster_window[2][0]; iz<=cluster_window[2][1]; iz++) {
          index = ny*nz*ix + nz*iy + iz;
          set_mask_bit(mask, index);
        }
      }
    }
  }
  memcpy(mask, mask_backup, nchar);
  free(mask_backup);
  return(cluster_set);
}


void run_cluster(SETTINGS *settings) {
  MOLECULE *ligand;
  MOLECULE_LIST *ligand_list;
  CLUSTER_SET *cluster_set;
  MOLECULE_PROBE *probe;
  int n_atoms;
  ULLMAN_MAPPINGS *mappings;

  ligand_list = settings->sysmols->ligand_list;
  ligand = ligand_list->molecules[0];
  prep_molecule(ligand, settings);
  n_atoms = ligand->selection->natoms;
  cluster_set = alloc_div_cluster_set(ligand_list->n_molecules, n_atoms, 0);
  probe = alloc_molecule_probe();

  // get scores
  float scores[ligand_list->n_molecules];
  char *scores_file = params_get_parameter("scores_file")->value.s;
  PLI_FILE *f = open_file(scores_file, "r");
  if (!f) {
    error_fn("%s: File %s not found", __func__, scores_file);
  }
  char line[MAX_LINE_LEN];
  char tag[] = "\"pliff_score\":";

  for (int i=0; i<ligand_list->n_molecules; i++) {
    if ((f == NULL) || (read_line(line, MAX_LINE_LEN, f) == NULL)) {
      error_fn("%s: Not enough score lines in the file %s", __func__, scores_file);
    }

    char *ptr = strstr(line, tag);
    if (ptr == NULL) {
      error_fn("%s: tag (%s) not found in the line %s", __func__, tag, line);
    }
    char *token = strtok(ptr + strlen(tag), ",");
    if (token == NULL) {
      error_fn("%s: score can not be read from the line %s", __func__, line);
    }
    int success = sscanf(token, "%f", &(scores[i]));
    if (! success) {
      error_fn("%s: score can not be read from the token %s", __func__, token);
    }
  }


  for (int i=0; i<ligand_list->n_molecules; i++) {
    ligand = ligand_list->molecules[i];
    prep_molecule(ligand, settings);
    add_molecule_center(ligand, 0);
    if (ligand->selection->natoms != n_atoms) {
      error_fn("Ligand %d in the file has different number of atoms", i+1);
    }
    probe->molecule = ligand;
    probe->score = scores[i];
    if (i==0) {
      mappings = find_ullman_mappings(ligand, ligand, MATCH_ALL, MATCH_ATOM_TYPE);
      probe->auto_mappings = mappings->mappings;
      probe->n_mappings = mappings->n_mappings;
    }

    save_diverse_pose(cluster_set, probe, NULL);
    unprep_molecule(ligand, settings);
  }
  prep_molecule(probe->molecule, settings);
  write_cluster_set(cluster_set, "diverse_cluster", probe, NULL);
}
