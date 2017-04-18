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



#define DS_N_VECTORS 13
#define DS_N_SPECTRA 3
#define DS_N_SPECTRUM_BINS 13
#define DS_MAX_DIST 6.0
#define DS_SPECTRUM_WIDTH 0.50



static double ds_spectra_counts[3][13] = { {  0.393,  2.405,  5.102,  5.452,  3.803,  2.250,  1.356,  1.142,  1.059,  0.699,  0.430,  0.413,  1.497  },
					   {  0.261,  1.505,  3.511,  4.889,  4.960,  3.650,  1.830,  0.969,  0.800,  0.528,  0.393,  0.494,  2.211  },
					   {  1.500,  4.125,  4.539,  3.433,  2.991,  2.536,  1.622,  1.072,  0.939,  0.656,  0.404,  0.336,  1.847  } };

static double ds_spectra_sigmas[3][13] = { {  0.32,  0.70,  0.69,  0.69,  0.73,  0.48,  0.42,  0.45,  0.48,  0.39,  0.26,  0.36,  0.51  },
					   {  0.28,  0.59,  0.45,  0.58,  0.56,  0.40,  0.39,  0.39,  0.36,  0.25,  0.23,  0.30,  0.61  },
					   {  0.44,  0.72,  0.57,  0.59,  0.57,  0.37,  0.43,  0.51,  0.43,  0.30,  0.25,  0.30,  0.59  } };

static int ds_vectors[13][3] = { { 1, 0, 0 },  { 0, 1, 0 },  { 0, 0, 1 },  { 1, 1, 1 },  { 1, 1,-1 },
				 { 1,-1, 1 },  { 1,-1,-1 },  { 0, 1, 1 },  { 0, 1,-1 },  { 1, 0, 1 },
				 { 1, 0,-1 },  { 1, 1, 0 },  { 1,-1, 0 } };

static double ds_distances[13] = { 1.0, 1.0, 1.0, 1.732, 1.732, 1.732, 1.732, 1.414, 1.414, 1.414, 1.414, 1.414, 1.414 };

static int ds_max_steps[13];

static double ds_gridspacing;
static int ds_smoothe;



static void digsite_init(void);
static double digsite_score(int,int,int,unsigned char*,GRID*);
static void digsite_update_spectrum(double*,int);
static int digsite_bump(int,int,int,int,unsigned char*,GRID*,int);
static double digsite_zscore(double*,double*,double*);



MAP* digsite_map(SYSTEM *system) {

  int i,ix,iy,iz,nx,ny,nz,index;
  double ***matrix;
  GRID *grid;
  MAP *map;
  BYTE *volume_mask,*clash_mask;

  digsite_init();

  reset_system(system,0);

  grid = system_grid(system,(params_get_parameter("grid_spacing"))->value.d,(params_get_parameter("grid_padding"))->value.d);

  map = new_map("digsite","double",grid);

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

	  matrix[ix][iy][iz] = digsite_score(ix,iy,iz,clash_mask,grid);

	} else {

	  matrix[ix][iy][iz] = 0.0;
	}
      }
    }
  }

  free(volume_mask);
  free(clash_mask);

  inflate_map(map);
  average_map(map);

  return(map);
}



static void digsite_init(void) {

  int i;
  
  ds_smoothe = (params_get_parameter("digsite_smoothe"))->value.i;
  ds_gridspacing = (params_get_parameter("grid_spacing"))->value.d;

  for (i=0;i<DS_N_VECTORS;i++) {
    
    ds_max_steps[i] = round(DS_MAX_DIST/(ds_gridspacing*ds_distances[i]));
  }
}



static double digsite_score(int ix,int iy,int iz,unsigned char* mask,GRID *grid) {

  int i,bin1,bin2;
  double score,zscore,zmin,counts[DS_N_SPECTRUM_BINS],*ds_counts,*ds_sigmas;

  for (i=0;i<DS_N_SPECTRUM_BINS;i++) {

    counts[i] = 0.0;
  }

  for (i=0;i<DS_N_VECTORS;i++) {
      
    bin1 = digsite_bump(ix,iy,iz,i,mask,grid,1.0);
    bin2 = digsite_bump(ix,iy,iz,i,mask,grid,-1.0);

    digsite_update_spectrum(counts,bin1);
    digsite_update_spectrum(counts,bin2);
  }

  zmin = 1.0E10;

  for (i=0;i<DS_N_SPECTRA;i++) {

    ds_counts = ds_spectra_counts[i];
    ds_sigmas = ds_spectra_sigmas[i];

    zscore = digsite_zscore(counts,ds_counts,ds_sigmas);

    if (zscore < zmin) {

      zmin = zscore;
    }
  }

  score = 2.0 - zmin;

  if (score < 0.0) {

    score = 0.0;
  }

  return(score);
}



static double digsite_zscore(double *counts,double *ds_counts,double *ds_sigmas) {

  int i;
  double sumz,z;

  sumz = 0.0;

  for (i=0;i<DS_N_SPECTRUM_BINS;i++) {

    sumz += fabs(counts[i]-ds_counts[i])/ds_sigmas[i];
  }

  z = sumz/DS_N_SPECTRUM_BINS;

  return(z);
}



static void digsite_update_spectrum(double *counts,int bin) {

  if (bin <= 0) {

    counts[0] += 0.75;
    counts[1] += 0.25;

  } else if (bin < DS_N_SPECTRUM_BINS-2) {

    counts[bin-1] += 0.25;
    counts[bin] += 0.50;
    counts[bin+1] += 0.25;

  } else if (bin == DS_N_SPECTRUM_BINS-2) {
    
    counts[bin-1] += 0.25;
    counts[bin] += 0.75;

  } else if (bin >= DS_N_SPECTRUM_BINS-1) {

    counts[DS_N_SPECTRUM_BINS-1] += 1.0;
  }
}



static int digsite_bump(int ix,int iy,int iz,int iv,unsigned char* mask,GRID *grid,int sign) {

  int steps,index,*lsv,bin;
  double dist;

  lsv = ds_vectors[iv];

  steps = 0;

  do {

    steps += 1;

    ix += sign*lsv[0];
    iy += sign*lsv[1];
    iz += sign*lsv[2];

    index = map_grid2index(ix,iy,iz,grid);

    if (index == -1) {

      return(12);
    }

    if (steps > ds_max_steps[iv]) {

      return(12);
    }
    
  } while  (!is_mask_bit_set(mask,index));

  dist = ((double) steps)*ds_gridspacing*ds_distances[iv];

  bin = floor(dist/DS_SPECTRUM_WIDTH);

  return(bin);
}
