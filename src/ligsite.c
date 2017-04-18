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



static int ls_vectors[13][3] = { { 1, 0, 0 },
				 { 0, 1, 0 },
				 { 0, 0, 1 },
				 { 1, 1, 1 },
				 { 1, 1,-1 },
				 { 1,-1, 1 },
				 { 1,-1,-1 },
				 { 0, 1, 1 },
				 { 0, 1,-1 },
				 { 1, 0, 1 },
				 { 1, 0,-1 },
				 { 1, 1, 0 },
				 { 1,-1, 0 } };

static double ls_distances[13] = { 1.0, 1.0, 1.0, 1.732, 1.732, 1.732, 1.732, 1.414, 1.414, 1.414, 1.414, 1.414, 1.414 };

static int ls_max_steps[13];
static int ls_n_vectors = 13;

static int ls_enclosed = 1;
static int ls_smoothe = 0;
static int ls_cap_dist = 0;
static double ls_max_dist = 8.0;

static double ls_grid_spacing = 0.5;
static double ls_grid_padding = 2.0;



static void ligsite_init(void);
static double ligsite_score(int,int,int,unsigned char*,GRID*);
static int ligsite_bump(int,int,int,int,unsigned char*,GRID*,int);



MAP* ligsite_map(SYSTEM *system) {

  int i,ix,iy,iz,nx,ny,nz,index;
  BYTE *volume_mask,*clash_mask;
  double ***matrix;
  GRID *grid;
  MAP *map;

  ligsite_init();

  reset_system(system,0);

  grid = system_grid(system,(params_get_parameter("grid_spacing"))->value.d,(params_get_parameter("grid_padding"))->value.d);

  map = new_map("ligsite","double",grid);

  volume_mask = mask_system_volume(system,NULL,grid,PROBE_MASK,6.0);
  clash_mask = mask_system_atoms(system,NULL,grid,PROBE_MASK,1.0);

  matrix = (double***) map->matrix;

  nx = grid->npoints[0];
  ny = grid->npoints[1];
  nz = grid->npoints[2];

  for (ix=0;ix<nx;ix++) {

    for (iy=0;iy<ny;iy++) {

      for (iz=0;iz<nz;iz++) {

	index = ny*nz*ix + nz*iy + iz;

	if ((!(is_mask_bit_set(clash_mask,index))) && (is_mask_bit_set(volume_mask,index))) {

	  matrix[ix][iy][iz] = ligsite_score(ix,iy,iz,clash_mask,grid);

	} else {

	  matrix[ix][iy][iz] = 0.0;
	}
      }
    }
  }

  free(volume_mask);
  free(clash_mask);

  for (i=0;i<ls_smoothe;i++) {

    smoothe_map(map);
  }

  return(map);
}



double ligsite_score(int ix,int iy,int iz,unsigned char* mask,GRID *grid) {

  int i;
  double score;

  score = 0.0;

  for (i=0;i<ls_n_vectors;i++) {

    if (ls_enclosed) {

      if ((ligsite_bump(ix,iy,iz,i,mask,grid,1.0)) && (ligsite_bump(ix,iy,iz,i,mask,grid,-1.0))) {

	score += 1.0;
      }

    } else {
      
      if (ligsite_bump(ix,iy,iz,i,mask,grid,1.0)) {

	score += 0.5;
      }

      if (ligsite_bump(ix,iy,iz,i,mask,grid,-1.0)) {

	score += 0.5;
      }
    }
  }

  return(score);
}



static void ligsite_init(void) {

  int i;
  double grid_spacing;

  ls_enclosed = (params_get_parameter("ligsite_enclosed"))->value.i;
  ls_smoothe = (params_get_parameter("ligsite_smoothe"))->value.i;
  ls_cap_dist = (params_get_parameter("ligsite_cap_dist"))->value.i;
  ls_max_dist = (params_get_parameter("ligsite_max_dist"))->value.d;

  grid_spacing = (params_get_parameter("grid_spacing"))->value.d;

  if (ls_cap_dist) {

    for (i=0;i<ls_n_vectors;i++) {

      ls_max_steps[i] = round(ls_max_dist/(grid_spacing*ls_distances[i]));
    }
  }
}


static int ligsite_bump(int ix,int iy,int iz,int iv,unsigned char* mask,GRID *grid,int sign) {

  int steps,index,*lsv;

  lsv = ls_vectors[iv];

  steps = 0;

  do {

    steps += 1;

    ix += sign*lsv[0];
    iy += sign*lsv[1];
    iz += sign*lsv[2];

    index = map_grid2index(ix,iy,iz,grid);

    if (index == -1) {

      return(0);
    }

    if (ls_cap_dist) {

      if (steps > ls_max_steps[iv]) {

	return(0);
      }
    }
    
  } while  (!is_mask_bit_set(mask,index));

  return(1);
}
