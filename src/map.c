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



static double*** alloc_3d_double_matrix(int,int,int);
static ATOMLIST*** alloc_3d_atomlist_matrix(int,int,int);
static unsigned char* alloc_3d_mask(GRID*);
static int map_max_index(MAP*,int,int*,int*,int*);
static void init_map_island(MAP_ISLAND*,int);
static void grow_map_island(MAP_ISLAND*,MAP*,double,int);
static void add_point_to_map_island(MAP_ISLAND*,int);



MAP* new_map(char *name) {

  MAP *map;

  map = (MAP*) malloc(sizeof(MAP));

  if (map == NULL) {

    error_fn("new_map: out of memory allocating map");
  }

  map->grid = NULL;
  map->matrix = NULL;
  map->mask = NULL;

  strcpy(map->name,name);

  return(map);
}



MAP* new_system_map(char *name,SYSTEM *system,GRID *grid,enum MAP_TYPE type) {

  double spacing,padding;
  MAP *map;

  map = new_map(name);

  if (map == NULL) {

    error_fn("new_system_map: out of memory allocating map");
  }

  if (grid) {

    map->grid = (GRID*) malloc(sizeof(GRID));

    if (map->grid == NULL) {

      error_fn("new_system_map: out of memory allocating grid");
    }

    memcpy(map->grid,grid,sizeof(GRID));

  } else {

    spacing = (params_get_parameter("grid_spacing"))->value.d;
    padding = (params_get_parameter("grid_padding"))->value.d;

    map->grid = atomlist2grid(system->selection,spacing,padding,NULL);
  }

  map->type = type;

  alloc_map_matrix(map);

  return(map);
}



void free_map(MAP *map) {

  if (map) {

    if (map->matrix != NULL) {

      free(**(map->matrix));
      free(*(map->matrix));
      free(map->matrix);
    }

    if (map->mask != NULL) {

      free(map->mask);
    }
  }
}



void alloc_map_matrix(MAP *map) {

  map->matrix = (void***) alloc_3d_matrix(map->grid,map->type);

  if (map->matrix == NULL) {

    error_fn("alloc_map_matrix: out of memory allocating map matrix");
  }

  map->mask = (unsigned char*) alloc_3d_mask(map->grid);

  init_3d_matrix(map);
  init_3d_mask(map);
}



GRID* new_grid(void) {

  GRID *grid;

  grid = (GRID*) malloc(sizeof(GRID));

  if (grid == NULL) {

    error_fn("new_grid: out of memory allocating grid");
  }

  return(grid);
}



void init_grid(GRID *grid,double spacing,double padding) {

  int i;

  grid->spacing = spacing;
  grid->padding = padding;

  for (i=0;i<3;i++) {

    grid->flimit[i][0] =  1.0E30;
    grid->flimit[i][1] = -1.0E30;

    grid->limit[i][0] =  99999;
    grid->limit[i][1] = -99999;
  }
}



void write_grid(FILE *file,GRID *grid) {

  int i;

  for (i=0;i<3;i++) {

    fprintf(file,"%8.3lf %8.3lf   %5d\n",grid->flimit[i][0],grid->flimit[i][1],grid->npoints[i]);
  }
}



void update_gridlimits(GRID *grid,double *v) {

  int i,p1,p2;

  for (i=0;i<3;i++) {

    p1 = floor((v[i] - grid->padding)/grid->spacing);
    p2 = floor((v[i] + grid->padding)/grid->spacing);

    if (p1 < grid->limit[i][0]) {

      grid->limit[i][0] = p1;
      grid->flimit[i][0] = ((double) p1)*(grid->spacing);
    }

    if (p2 > grid->limit[i][1]) {

      grid->limit[i][1] = p2;
      grid->flimit[i][1] = ((double) p2)*(grid->spacing);
    }
  }
}



int pos2grid(double *pos,int *ipos,GRID *grid) {

  int i,*iv;
  double *v;
  
  for (i=0,v=pos,iv=ipos;i<3;i++,v++,iv++) {

    *iv = floor((*v - grid->flimit[i][0])/grid->spacing);

    if ((*iv < 0) || (*iv >= grid->npoints[i])) {

      return(1);
    }
  }

  return(0);
}



int pos2grid_round(double *pos,int *ipos,GRID *grid) {

  int i,*iv;  double *v;

  for (i=0,v=pos,iv=ipos;i<3;i++,v++,iv++) {

    *iv = round((*v - grid->flimit[i][0])/grid->spacing);

    if ((*iv < 0) || (*iv >= grid->npoints[i])) {

      return(1);
    }
  }

  return(0);
}



int map_index2grid(int index,int *ix,int *iy,int *iz,GRID *grid) {

  int nx,ny,nz;

  nx = grid->npoints[0];
  ny = grid->npoints[1];
  nz = grid->npoints[2];

  if ((index < 0) || (index >= (nx*ny*nz))) {

    return(1);
  }

  *ix = (index/(ny*nz)) % nx;
  *iy = (index/nz) % ny;
  *iz = index % nz;

  return(0);
}



int map_grid2index(int ix,int iy,int iz,GRID *grid) {

  int nx,ny,nz;

  nx = grid->npoints[0];
  ny = grid->npoints[1];
  nz = grid->npoints[2];

  if ((ix < 0) || (ix >= nx) || (iy < 0) || (iy >= ny) || (iz < 0) || (iz >= nz)) {

    return(-1);
  }

  return(ny*nz*ix + nz*iy + iz);
}



void*** alloc_3d_matrix(GRID *grid,enum MAP_TYPE type) {

  int nx,ny,nz;
  void ***matrix;

  nx = grid->npoints[0];
  ny = grid->npoints[1];
  nz = grid->npoints[2];

  if (type == ATOMLIST_MAP) {

    matrix = (void***) alloc_3d_atomlist_matrix(nx,ny,nz);

  } else if (type == DOUBLE_MAP) {

    matrix = (void***) alloc_3d_double_matrix(nx,ny,nz);    

  } else {

    error_fn("alloc_3d_matrix: unkown map type");
  }

  return(matrix);
}



void init_3d_matrix(MAP *map) {

  int i,j,k,nx,ny,nz;
  double ***dm;
  ATOMLIST ***am;

  nx = map->grid->npoints[0];
  ny = map->grid->npoints[1];
  nz = map->grid->npoints[2];

  if (map->type == ATOMLIST_MAP) {

    am = (ATOMLIST***) map->matrix;

    for (i=0;i<nx;i++) {

      for (j=0;j<ny;j++) {

	for (k=0;k<nz;k++) {

	  init_atomlist(&(am[i][j][k]));
	}
      }
    }

  } else if (map->type == DOUBLE_MAP) {

    dm = (double***) map->matrix;

    for (i=0;i<nx;i++) {

      for (j=0;j<ny;j++) {

	for (k=0;k<nz;k++) {

	  dm[i][j][k] = 0.0;
	}
      }
    }

  } else {

    error_fn("init_3d_matrix: unknown map type");
  }
}



static unsigned char* alloc_3d_mask(GRID *grid) {

  int nx,ny,nz,nchars;
  unsigned char *mask;

  nx = grid->npoints[0];
  ny = grid->npoints[1];
  nz = grid->npoints[2];

  nchars = (nx*ny*nz)/8 + 1;

  mask = (unsigned char*) calloc(nchars,sizeof(unsigned char));

  if (mask == NULL) {

    error_fn("alloc_3d_mask: out of memory allocating mask");
  }

  return(mask);
}



void init_3d_mask(MAP *map) {

  int nx,ny,nz,nchars;
  GRID *grid;

  grid = map->grid;

  nx = grid->npoints[0];
  ny = grid->npoints[1];
  nz = grid->npoints[2];

  nchars = (nx*ny*nz)/8 + 1;

  memset(map->mask,0,nchars);
}



static double*** alloc_3d_double_matrix(int nx,int ny,int nz) {

  int i,j,address = 0;
  double ***matrix;

  matrix = (double***) calloc(nx,sizeof(double**));

  if (matrix == NULL)
    return(NULL);

  matrix[0] = (double**) calloc(nx*ny,sizeof(double*));

  if (matrix[0] == NULL)
    return(NULL);
 
  for (i=1;i<nx;i++)
    matrix[i] = matrix[i-1] + ny;

  matrix[0][0] = (double*) calloc(nx*ny*nz,sizeof(double));

  if (matrix[0][0] == NULL)
    return(NULL);

  for (i=0;i<nx;i++) {

    for (j=0;j<ny;j++) {

      matrix[i][j] = matrix[0][0] + address;
      address += nz;
    }
  }

  return(matrix);
}



static ATOMLIST*** alloc_3d_atomlist_matrix(int nx,int ny,int nz) {

  int i,j,address = 0;
  ATOMLIST ***matrix;

  matrix = (ATOMLIST***) calloc(nx,sizeof(ATOMLIST**));

  if (matrix == NULL)
    return(NULL);

  matrix[0] = (ATOMLIST**) calloc(nx*ny,sizeof(ATOMLIST*));

  if (matrix[0] == NULL)
    return(NULL);
 
  for (i=1;i<nx;i++)
    matrix[i] = matrix[i-1] + ny;

  matrix[0][0] = (ATOMLIST*) calloc(nx*ny*nz,sizeof(ATOMLIST));

  if (matrix[0][0] == NULL)
    return(NULL);

  for (i=0;i<nx;i++) {

    for (j=0;j<ny;j++) {

      matrix[i][j] = matrix[0][0] + address;
      address += nz;
    }
  }

  return(matrix);
}



void mask_map_molecule(MAP *map,MOLECULE *molecule,int type) {

  int i,j,imax,iv_min[3],iv_max[3],iv[3],ix,iy,iz,ny,nz,index;
  double *v,dmax,sqr_dmax,dx,dy,dz;
  ATOM *atom;
  unsigned char *mask;
  GRID *grid;

  grid = map->grid;
  mask = map->mask;

  ny = grid->npoints[1];
  nz = grid->npoints[2];

  for (i=0,atom=molecule->atom;i<molecule->natoms;i++,atom++) {

    if ((!(atom->flags & SKIP_ATOM)) && (atom->element->id != HYDROGEN)) {

      v = atom->position;

      if (type == VDW_MASK) {

	dmax = (atom->element->flags & METAL_ELEMENT) ? 2.0 : atom->vdw_radius_H2O;

      } else {

	if (atom->type->flags & ANY_HBOND_ATOM_TYPE) {

	  dmax = (type == CONTACT_MASK) ? 2.0 : 1.2;

	} else if (atom->element->flags & METAL_ELEMENT) {

	  dmax = (type == CONTACT_MASK) ? 1.5 : 0.8;

	} else {

	  dmax = (type == CONTACT_MASK) ? atom->vdw_radius + 0.8 : atom->vdw_radius - 0.2;
	}
      }

      if (point_in_grid(v,grid,dmax)) {

	point_to_gridpoint(v,iv,grid);

	imax = ceil(dmax/grid->spacing);

	for (j=0;j<3;j++) {

	  iv_min[j] = (iv[j] - imax < 0) ? 0 : iv[j] - imax;
	  iv_max[j] = (iv[j] + imax >= grid->npoints[j]) ? grid->npoints[j] - 1 : iv[j] + imax;
	}

	sqr_dmax = dmax*dmax;

	for (ix=iv_min[0];ix<=iv_max[0];ix++) {

	  dx = grid->flimit[0][0] + ((double) ix)*(grid->spacing) - v[0];

	  for (iy=iv_min[1];iy<=iv_max[1];iy++) {

	    dy = grid->flimit[1][0] + ((double) iy)*(grid->spacing) - v[1];

	    for (iz=iv_min[2];iz<=iv_max[2];iz++) {

	      index = ny*nz*ix + nz*iy + iz;

	      dz = grid->flimit[2][0] + ((double) iz)*(grid->spacing) - v[2];
	      
	      if (!is_mask_bit_set(mask,index)) {

		if (((dx*dx)+(dy*dy)+(dz*dz)) < sqr_dmax) {

		  set_mask_bit(mask,index);
		}
	      }
	    }
	  }  
	}
      }
    }
  }
}



int point_in_grid(double *v,GRID *grid,double border) {

  int i;

  for (i=0;i<3;i++) {

    if (((v[i]+border) < grid->flimit[i][0]) || ((v[i]-border) > grid->flimit[i][1])) {

      return(0);
    }
  }
  return(1);
}


void point_to_gridpoint(double *v,int *iv,GRID *grid) {

  int i;

  for (i=0;i<3;i++) {

    iv[i] = round((v[i]-grid->flimit[i][0])/(grid->spacing));
  }
}



void gridpoint_to_point(int *iv,double *v,GRID *grid) {

  int i;
  
  for (i=0;i<3;i++) {

    v[i] = grid->flimit[i][0] + ((double) iv[i])*grid->spacing;
  }

  v[3] = 1.0;
}



void set_mask_bit(unsigned char *mask,int index) {

  int minor,major;
  unsigned char *byte;

  major = index / 8;
  minor = index % 8;

  byte = mask + major;

  *byte |= 01<<minor;
}



int is_mask_bit_set(unsigned char *mask,int index) {

  int minor,major;
  unsigned char *byte;

  major = index / 8;

  byte = mask + major;

  if (!(*byte)) {

    return(0);
  }

  minor = index % 8;

  if (*byte & 01<<minor) {

    return(1);
  }

  return(0);
}



int point_in_island(int index,MAP_ISLAND *island) {

  int i,*point;

  for (i=0,point=island->points;i<island->n_points;i++,point++) {

    if (index == *point) {

      return(1);
    }
  }

  return(0);
}



MAP_ISLAND* map2islands(MAP *map,double level,int *n_islands) {

  int i,nx,ny,nz,index,ix,iy,iz,n_alloc_islands;
  double ***dm;
  GRID *grid;
  MAP_ISLAND *islands,*island;

  grid = map->grid;

  nx = grid->npoints[0];
  ny = grid->npoints[1];
  nz = grid->npoints[2];

  init_3d_mask(map);

  if (map->type == DOUBLE_MAP) {

    dm = (double***) map->matrix;

    n_alloc_islands = 100;

    islands = (MAP_ISLAND*) calloc(n_alloc_islands,sizeof(MAP_ISLAND));

    if (islands == NULL) {

      error_fn("map2islands: out of memory allocating islands");
    }

    *n_islands = 0;
    
    while (((index = map_max_index(map,1,&ix,&iy,&iz)) != -1) && (dm[ix][iy][iz] >= level)) {
      
      if (*n_islands == n_alloc_islands) {

	n_alloc_islands *= 2;

	islands = (MAP_ISLAND*) realloc(islands,n_alloc_islands*sizeof(MAP_ISLAND));

	if (islands == NULL) {

	  error_fn("map2islands: out of memory re-allocating islands");
	}
      }

      island = islands + (*n_islands);
      
      init_map_island(island,1000);
      
      grow_map_island(island,map,level,index);
      
      (*n_islands)++;
    }

    if (*n_islands > 0) {

      islands = (MAP_ISLAND*) realloc(islands,(*n_islands)*sizeof(MAP_ISLAND));

      for (i=0,island=islands;i<*n_islands;i++,island++) {

	if (island->n_points) {

	  island->n_alloc_points = island->n_points;

	  island->points = (int*) realloc(island->points,(island->n_points)*sizeof(int));
	}

	island->id = i+1;

	calc_map_island_properties(island,map);
      }
    }
  }

  return(islands);
}



MAP_ISLAND* index2island(int index,MAP_ISLAND *islands,int n_islands) {

  int i;
  MAP_ISLAND *island;

  for (i=0,island=islands;i<n_islands;i++,island++) {

    if (point_in_island(index,island)) {

      return(island);
    }
  }

  return(NULL);
}



double island_mask_overlap(MAP_ISLAND *island,MAP *map,double level) {

  int i,n_total,n_overlap,*point,index,ix,iy,iz;
  double ***dm;
  unsigned char *mask;
  GRID *grid;

  if (map->type != DOUBLE_MAP) {

    error_fn("island_mask_overlap: map type '%d' unsupported",(int) map->type);
  }

  if (island->n_points <= 0) {

    return(0.0);
  }

  grid = map->grid;
  mask = map->mask;

  dm = (double***) map->matrix;

  n_total = n_overlap = 0;

  for (i=0,point=island->points;i<island->n_points;i++,point++) {

    index = *point;

    if (!map_index2grid(index,&ix,&iy,&iz,grid)) {

      if (dm[ix][iy][iz] > level) {

	if (is_mask_bit_set(mask,index)) {

	  n_overlap++;
	}

	n_total++;
      }
    }
  }

  if (n_total > 0) {
    
    return(((double) n_overlap)/((double) n_total));
  }

  return(0.0);
}



double island_map_overlap(MAP_ISLAND *island,MAP *map) {

  int i,*point,index,ix,iy,iz;
  double ***dm,sumY;
  unsigned char *mask;
  GRID *grid;

  if (map->type != DOUBLE_MAP) {

    error_fn("island_mask_overlap: map type '%d' unsupported",(int) map->type);
  }

  if (island->n_points <= 0) {

    return(0.0);
  }

  grid = map->grid;
  mask = map->mask;

  dm = (double***) map->matrix;

  sumY = 0.0;

  for (i=0,point=island->points;i<island->n_points;i++,point++) {

    index = *point;

    if (!map_index2grid(index,&ix,&iy,&iz,grid)) {

      if (is_mask_bit_set(mask,index)) {

	sumY += dm[ix][iy][iz];
      }
    }
  }

  return(sumY*pow(grid->spacing,3));
}



static void grow_map_island(MAP_ISLAND *island,MAP *map,double level,int index) {

  int ix,iy,iz;
  double ***dm;
  GRID *grid;
  unsigned char *mask;

  if (index == -1) {

    return;
  }

  mask = map->mask;

  if (is_mask_bit_set(mask,index)) {

    return;
  }

  grid = map->grid;

  if (map_index2grid(index,&ix,&iy,&iz,grid)) {

    return;
  }

  if (map->type == DOUBLE_MAP) {

    dm = (double***) map->matrix;

    if (dm[ix][iy][iz] < level) {

      return;
    }

    set_mask_bit(mask,index);

    add_point_to_map_island(island,index);

    grow_map_island(island,map,level,map_grid2index(ix-1,iy,iz,grid));
    grow_map_island(island,map,level,map_grid2index(ix+1,iy,iz,grid));
    grow_map_island(island,map,level,map_grid2index(ix,iy-1,iz,grid));
    grow_map_island(island,map,level,map_grid2index(ix,iy+1,iz,grid));
    grow_map_island(island,map,level,map_grid2index(ix,iy,iz-1,grid));
    grow_map_island(island,map,level,map_grid2index(ix,iy,iz+1,grid));

    grow_map_island(island,map,level,map_grid2index(ix-1,iy-1,iz,grid));
    grow_map_island(island,map,level,map_grid2index(ix-1,iy+1,iz,grid));
    grow_map_island(island,map,level,map_grid2index(ix+1,iy-1,iz,grid));
    grow_map_island(island,map,level,map_grid2index(ix+1,iy+1,iz,grid));

    grow_map_island(island,map,level,map_grid2index(ix-1,iy,iz-1,grid));
    grow_map_island(island,map,level,map_grid2index(ix-1,iy,iz+1,grid));
    grow_map_island(island,map,level,map_grid2index(ix+1,iy,iz-1,grid));
    grow_map_island(island,map,level,map_grid2index(ix+1,iy,iz+1,grid));

    grow_map_island(island,map,level,map_grid2index(ix,iy-1,iz-1,grid));
    grow_map_island(island,map,level,map_grid2index(ix,iy-1,iz+1,grid));
    grow_map_island(island,map,level,map_grid2index(ix,iy+1,iz-1,grid));
    grow_map_island(island,map,level,map_grid2index(ix,iy+1,iz+1,grid));

    grow_map_island(island,map,level,map_grid2index(ix-1,iy-1,iz-1,grid));
    grow_map_island(island,map,level,map_grid2index(ix+1,iy+1,iz+1,grid));

    grow_map_island(island,map,level,map_grid2index(ix-1,iy-1,iz+1,grid));
    grow_map_island(island,map,level,map_grid2index(ix-1,iy+1,iz-1,grid));
    grow_map_island(island,map,level,map_grid2index(ix+1,iy-1,iz-1,grid));

    grow_map_island(island,map,level,map_grid2index(ix+1,iy+1,iz-1,grid));
    grow_map_island(island,map,level,map_grid2index(ix+1,iy-1,iz+1,grid));
    grow_map_island(island,map,level,map_grid2index(ix-1,iy+1,iz+1,grid));

  } else {

    error_fn("grow_map_island: unsupported map type '%d'",map->type);
  }
}



void calc_map_island_properties(MAP_ISLAND *island,MAP *map) {

  int i,*index,ix,iy,iz;
  double ***dm,Y,maxY,sumY;
  GRID *grid;

  if (island->n_points > 0) {

    grid = map->grid;

    if (map->type == DOUBLE_MAP) {

      dm = (double***) map->matrix;

      sumY = 0.0;
      maxY = -999.9E10;

      for (i=0,index=island->points;i<island->n_points;i++,index++) {

	if (!map_index2grid(*index,&ix,&iy,&iz,grid)) {

	  Y = dm[ix][iy][iz];

	  if (Y > maxY) {

	    maxY = Y;
	  }

	  sumY += Y;
	}
      }

      island->maxY = maxY;
      island->sumY = sumY;
      island->avgY = sumY/((double) island->n_points);
      island->integral = sumY*pow(grid->spacing,3);
      island->volume = ((double) island->n_points)*pow(grid->spacing,3);
    }
  }
}



static void init_map_island(MAP_ISLAND *island,int n_alloc_points) {

  island->id = 0;

  island->n_points = 0;
  island->n_alloc_points = n_alloc_points;
  island->points = NULL;

  island->maxY = -999.9E10;
  island->sumY = 0.0;
  island->avgY = 0.0;
  island->integral = 0.0;
  island->volume = 0.0;

  if (n_alloc_points) {

    island->points = (int*) calloc(n_alloc_points,sizeof(int));

    if (island->points == NULL) {

      error_fn("init_map_island: out of memory allocating points");
    }
  }
}



static void add_point_to_map_island(MAP_ISLAND *island,int index) {

  if (island->n_points == island->n_alloc_points) {

    island->n_alloc_points *= 2;

    island->points = (int*) realloc(island->points,(island->n_alloc_points)*sizeof(int));

    if (island->points == NULL) {

      error_fn("add_point_to_map_island: out of memory re-allocating points");
    }
  }

  island->points[island->n_points] = index;

  island->n_points++;
}


static int map_max_index(MAP *map,int use_mask,int *mix,int *miy,int *miz) {

  int nx,ny,nz,ix,iy,iz,index,mindex;
  double ***dm,value,max_value;
  GRID *grid;
  unsigned char *mask;

  grid = map->grid;
  mask = map->mask;

  nx = grid->npoints[0];
  ny = grid->npoints[1];
  nz = grid->npoints[2];

  max_value = -9.99E10;

  mindex = -1;

  if (map->type == DOUBLE_MAP) {

    dm = (double***) map->matrix;

    for (ix=0;ix<nx;ix++) {

      for (iy=0;iy<ny;iy++) {

	for (iz=0;iz<nz;iz++) {

	  index = ny*nz*ix + nz*iy + iz;

	  value = dm[ix][iy][iz];

	  if ((value > max_value) && ((!use_mask) || (!is_mask_bit_set(mask,index)))) {

	    mindex = index;

	    *mix = ix;
	    *miy = iy;
	    *miz = iz;

	    max_value = value;
	  }
	}
      }
    }

  } else {

    error_fn("map_max_index: map type '%d' unsuitable",(int) map->type);
  }

  return(mindex);
}


