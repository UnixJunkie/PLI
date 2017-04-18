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



static MAP_TYPE map_types[] = {
  { "int",      sizeof(int),           v_init_int,      NULL },
  { "double",   sizeof(double),        v_init_double,   NULL },
  { "atomlist", sizeof(ATOMLIST),      v_init_atomlist, NULL },
  { "alist",    sizeof(LIST),          v_init_alist,    v_free_list_items },
  { "poselist", sizeof(POSELIST_LITE), NULL,            NULL },
  { "null",     0,                     NULL,            NULL },
  { "last",     0,                     NULL,            NULL }
};



static MAP_TYPE* get_map_type(char*);
static void*** alloc_3d_void_matrix(int,int,int,int);
static LIST*** alloc_3d_list_matrix(int,int,int);
static ATOMLIST*** alloc_3d_atomlist_matrix(int,int,int);
static void init_poselist_lite(POSELIST_LITE*);



MAP* new_map(char *name,char *type,GRID *grid) {

  MAP *map;

  map = (MAP*) malloc(sizeof(MAP));

  if (map == NULL) {

    error_fn("%s: out of memory allocating map",__func__);
  }

  map->type = get_map_type(type);

  map->grid = grid;

  if (grid) {

    alloc_map_matrix(map);

  } else {

    map->matrix = NULL;
    map->mask = NULL;
  }

  strcpy(map->name,name);

  return(map);
}



static MAP_TYPE* get_map_type(char *name) {

  MAP_TYPE *type;

  type = map_types;

  while (strcmp(type->name,"last")) {

    if (!strcmp(type->name,name)) {

      return(type);
    }

    type++;
  }

  error_fn("%s: unknown map type '%s'",__func__,name);
}



GRID* new_grid(double spacing,double padding) {

  GRID *grid;

  grid = (GRID*) malloc(sizeof(GRID));
  
  if (grid == NULL) {
    
    error_fn("%s: out of memory allocating grid",__func__);
  }
  
  init_grid(grid,spacing,padding);
  
  return(grid);
}



GRID* copy_grid(GRID *to_grid,GRID *from_grid) {

  to_grid = (to_grid) ? to_grid : new_grid(from_grid->spacing,from_grid->padding);

  memcpy(to_grid,from_grid,sizeof(GRID));

  return(to_grid);
}



MAP* copy_map(MAP *to_map,MAP *from_map) {

  int n_items,n_bytes;
  GRID *grid;

  grid = from_map->grid;

  to_map = (to_map) ? to_map : new_map(from_map->name,from_map->type->name,copy_grid(NULL,grid));

  n_items = grid_points(grid);

  if ((n_items) && (to_map->type->item_size)) {

    memcpy(to_map->matrix[0][0],from_map->matrix[0][0],n_items*(to_map->type->item_size));
  }

  copy_mask(to_map->mask,from_map->mask,grid);
}



int grid_points(GRID *grid) {

  int nx,ny,nz;

  nx = grid->npoints[0];
  ny = grid->npoints[1];
  nz = grid->npoints[2];

  return(nx*ny*nz);
}



void increment_map(MAP *map,unsigned char *mask,double value) {

  int nx,ny,nz,ix,iy,iz,index;
  GRID *grid;
  double ***matrix;

  if (strcmp(map->type->name,"double")) {

    error_fn("%s: map '%s' not of type double",__func__,map->name);
  }
  
  grid = map->grid;
  
  matrix = (double***) map->matrix;
  
  nx = grid->npoints[0];
  ny = grid->npoints[1];
  nz = grid->npoints[2];
  
  for (ix=0;ix<nx;ix++) {
    
    for (iy=0;iy<ny;iy++) {
      
      for (iz=0;iz<nz;iz++) {
	
	if (mask == NULL) {
	  
	  matrix[ix][iy][iz] += value;
	  
	} else {

	  index = ny*nz*ix + nz*iy + iz;
	  
	  if (is_mask_bit_set(mask,index)) {
	    
	    matrix[ix][iy][iz] += value;
	  }
	}
      }
    }
  }
}



void free_map(MAP *map) {

  int i,n_items,item_size;
  void *items,*item,(*free_item)(void*);

  if (map) {

    if (map->matrix != NULL) {

      free_item = map->type->free;

      if (free_item) {

	item_size = map->type->item_size;

	n_items = grid_points(map->grid);

	items = map->matrix[0][0];

	for (i=0;i<n_items;i++) {

	  item = (void*) (((char*) items) + (i*item_size)); 
	
	  free_item(item);
	}
      }

      free(**(map->matrix));
      free(*(map->matrix));
      free(map->matrix);
    }

    if (map->mask != NULL) {

      free(map->mask);
    }

    free(map);
  }
}



void alloc_map_matrix(MAP *map) {

  int item_size;
  GRID *grid;

  grid = map->grid;
  item_size = map->type->item_size;

  map->matrix = (void***) alloc_3d_matrix(grid,item_size);

  if ((item_size != 0) && (map->matrix == NULL)) {
      
    error_fn("%s: out of memory allocating map matrix for %s",__func__,map->name);
  }
    
  init_3d_matrix(map);

  map->mask = alloc_mask(grid);
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


void move_grid_center(GRID *grid, double new_center[4]) {
  int i;

  double old_center[4];
  for (i=0; i<3; i++) {
    old_center[i] = grid->flimit[i][0] + (grid->flimit[i][1] - grid->flimit[i][0]) / 2.0;
  }

  for (i=0; i<3; i++) {
    grid->flimit[i][0] += (new_center[i] - old_center[i]);
    grid->flimit[i][1] += (new_center[i] - old_center[i]);
    grid->limit[i][0] += (int) round((new_center[i] - old_center[i])/grid->spacing);
    grid->limit[i][1] += (int) round((new_center[i] - old_center[i])/grid->spacing);
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



double trilinear_interpolate(MAP *field, double *p, double max_score) {

  double v000, v100, v001, v101, v010, v110, v011, v111;
  double v00, v10, v01, v11, v0, v1, v, x0, x1, y0, y1, z0, z1, xd, yd, zd;
  double xd0, xd1, yd0, yd1, zd0, zd1;
  int ix, iy, iz, nx, ny, nz, ix1, iy1, iz1;

  GRID *grid;
  double ***matrix;
  matrix = (double***) field->matrix;
  grid = field->grid;

  nx = grid->npoints[0];
  ny = grid->npoints[1];
  nz = grid->npoints[2];

  if (! point_in_grid(p, grid, 0.0)) {
    error_fn("%s: [%.2f %.2f %.2f] point lies outside the grid\n", __func__, p[0], p[1], p[2]);
  }

  xd0 = fmod((p[0] - grid->flimit[0][0]), grid->spacing);
  //xd1 = grid->spacing - xd0;
  ix = floor((p[0] - grid->flimit[0][0]) / grid->spacing);
  xd = xd0/grid->spacing;

  yd0 = fmod((p[1] - grid->flimit[1][0]), grid->spacing);
  //yd1 = grid->spacing - yd0;
  iy = floor((p[1] - grid->flimit[1][0]) / grid->spacing);
  yd = yd0/grid->spacing;

  zd0 = fmod((p[2] - grid->flimit[2][0]), grid->spacing);
  //zd1 = grid->spacing - zd0;
  iz = floor((p[2] - grid->flimit[2][0]) / grid->spacing);
  zd = zd0/grid->spacing;

  ix1 = (ix==(nx-1))? ix: ix+1;
  iy1 = (iy==(ny-1))? iy: iy+1;
  iz1 = (iz==(nz-1))? iz: iz+1;

  v000 = matrix[ix][iy][iz];
  v100 = matrix[ix1][iy][iz];

  v010 = matrix[ix][iy1][iz];
  v110 = matrix[ix1][iy1][iz];

  v001 = matrix[ix][iy][iz1];
  v101 = matrix[ix1][iy][iz1];

  v011 = matrix[ix][iy1][iz1];
  v111 = matrix[ix1][iy1][iz1];

  if (v000 >= max_score || v100 >= max_score || v010 >= max_score || v110 >= max_score || \
      v001 >= max_score || v101 >= max_score || v011 >= max_score || v111 >= max_score) {
    return (max_score);
  }

  v00 = linear_interpolate(v000, v100, xd);
  v10 = linear_interpolate(v010, v110, xd);
  v01 = linear_interpolate(v001, v101, xd);
  v11 = linear_interpolate(v011, v111, xd);

  v0 = linear_interpolate(v00, v10, yd);
  v1 = linear_interpolate(v01, v11, yd);

  v = linear_interpolate(v0, v1, zd);

  return (v);
}


// TODO: one routine can now deal with any map type 
// should retire all map type specific routines for
// allocating memory etc.

void*** alloc_3d_matrix(GRID *grid,int item_size) {

  int nx,ny,nz;
  void ***matrix;

  nx = grid->npoints[0];
  ny = grid->npoints[1];
  nz = grid->npoints[2];

  matrix = (void***) alloc_3d_void_matrix(nx,ny,nz,item_size);

  return(matrix);
}



void init_3d_matrix(MAP *map) {

  int n_items,item_size,i;
  void (*init)(void*),*items,*item;

  item_size = map->type->item_size;

  if (item_size > 0) {

    init = map->type->init;

    if (init) {

      n_items = grid_points(map->grid);

      items = map->matrix[0][0];

      for (i=0;i<n_items;i++) {

	item = (void*) (((char*) items) + (i*item_size)); 
	
	init(item);
      }
    }
  }
}



double*** alloc_3d_double_matrix(int nx,int ny,int nz) {

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



static void*** alloc_3d_void_matrix(int nx,int ny,int nz,int item_size) {

  int i,j,address = 0;
  void ***matrix;

  if (item_size == 0) {

    return(NULL);
  }

  matrix = (void***) calloc(nx,sizeof(void**));

  if (matrix == NULL)
    return(NULL);

  matrix[0] = (void**) calloc(nx*ny,sizeof(void*));

  if (matrix[0] == NULL)
    return(NULL);

  for (i=1;i<nx;i++) {

    matrix[i] = matrix[i-1] + ny;
  }


  matrix[0][0] = (void*) calloc(nx*ny*nz,item_size);

  if (matrix[0][0] == NULL)
    return(NULL);

  for (i=0;i<nx;i++) {

    for (j=0;j<ny;j++) {

      matrix[i][j] = (void*) (((char*) matrix[0][0]) + (address*item_size));

      address += nz;
    }
  }

  return(matrix);
}



static LIST*** alloc_3d_list_matrix(int nx,int ny,int nz) {

  int i,j,address = 0;
  LIST ***matrix;

  matrix = (LIST***) calloc(nx,sizeof(LIST**));

  if (matrix == NULL)
    return(NULL);

  matrix[0] = (LIST**) calloc(nx*ny,sizeof(LIST*));

  if (matrix[0] == NULL)
    return(NULL);
 
  for (i=1;i<nx;i++)
    matrix[i] = matrix[i-1] + ny;

  matrix[0][0] = (LIST*) calloc(nx*ny*nz,sizeof(LIST));

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


POSELIST_LITE*** alloc_3d_poselist_matrix(int nx,int ny,int nz) {

  int i,j,address = 0;
  POSELIST_LITE ***matrix;

  matrix = (POSELIST_LITE***) calloc(nx,sizeof(POSELIST_LITE**));

  if (matrix == NULL)
    return(NULL);

  matrix[0] = (POSELIST_LITE**) calloc(nx*ny,sizeof(POSELIST_LITE*));

  if (matrix[0] == NULL)
    return(NULL);

  for (i=1;i<nx;i++)
    matrix[i] = matrix[i-1] + ny;

  matrix[0][0] = (POSELIST_LITE*) calloc(nx*ny*nz,sizeof(POSELIST_LITE));

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


void init_3d_poselist_matrix(MAP *map) {

  int i,j,k,nx,ny,nz;
  POSELIST_LITE ***pm;

  nx = map->grid->npoints[0];
  ny = map->grid->npoints[1];
  nz = map->grid->npoints[2];

  if (strcmp(map->type->name,"poselist")) {

    error_fn("%s: map '%s' not of type poselist",__func__,map->name);
  }

  pm = (POSELIST_LITE***) map->matrix;

  for (i=0;i<nx;i++) {

    for (j=0;j<ny;j++) {

      for (k=0;k<nz;k++) {

        init_poselist_lite(&(pm[i][j][k]));
      }
    }
  }
}

static void init_poselist_lite(POSELIST_LITE* l) {
  l->n_alloc_poses = 0;
  l->n_poses = 0;
  l->poses = NULL;
}


void mask_distant_points(MAP *map, SYSTEM *system, MOLECULE *molecule, ATOMLIST *atomlist, double max_dist, double angle) {
  /* masks all grid points that do not have any atom in the system,
   * molecule, and atomlist (whichever is not null)
   * within max_dist */



  int i,j,ix,iy,iz,nx,ny,nz,index;
  int mask_system, mask_molecule, mask_atomlist;
  double *v,dmax,sqr_dmax,dx,dy,dz;
  ATOM *atom;
  unsigned char *mask;
  GRID *grid;


  mask_system = (system) ? 1: 0;
  mask_molecule = (molecule) ? 1: 0;
  mask_atomlist = (atomlist) ? 1: 0;

  if ((mask_atomlist + mask_molecule + mask_system) != 1) {
    error_fn("%s: (only) one of the system, molecule, atomlist should be non-null");
  }

  grid = map->grid;
  mask = map->mask;

  nx = grid->npoints[0];
  ny = grid->npoints[1];
  nz = grid->npoints[2];

  int irange[3][2];
  double pos[3];

  for (ix=0; ix<nx; ix++) {
    pos[0] = grid->flimit[0][0] + ix*grid->spacing;
    for (iy=0; iy<ny; iy++) {
      pos[1] = grid->flimit[1][0] + iy*grid->spacing;
      for (iz=0; iz<nz; iz++) {
        index = ny*nz*ix + nz*iy + iz;
        if (! is_mask_bit_set(mask, index)) {
          pos[2] = grid->flimit[2][0] + iz*grid->spacing;
          if (system) {
            if (!position_within_dist_system(pos, system, max_dist, angle)) {
              set_mask_bit(mask, index);
            }
          }
          else if (molecule) {
            if (!position_within_dist_molecule(pos, molecule, max_dist, angle)) {
              set_mask_bit(mask, index);
            }
          }
          else if (atomlist) {
            if (! position_within_dist_atomlist(pos, atomlist, max_dist, angle)) {
              set_mask_bit(mask, index);
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


void integrate_field(MAP *field, double range, double ignore_value) {
  int range_gp, index, irange[3][2], i, j, k;
  GRID *grid = field->grid;
  unsigned char *mask = field->mask;
  double ***score_matrix = (double***) field->matrix;
  
  range_gp = (int) ceil(range/grid->spacing);

  int nx, ny, nz, ix, iy, iz;

  nx = grid->npoints[0];
  ny = grid->npoints[1];
  nz = grid->npoints[2];

  double ***integrated_field = alloc_3d_double_matrix(nx, ny, nz);


  for (ix=0; ix<nx; ix++) {
    irange[0][0] = ((ix-range_gp) < 0) ? 0 : ix-range_gp;
    irange[0][1] = ((ix+range_gp) > (nx-1)) ? nx-1 : ix+range_gp;
    for (iy=0; iy<ny; iy++) {
      irange[1][0] = ((iy-range_gp) < 0) ? 0 : iy-range_gp;
      irange[1][1] = ((iy+range_gp) > (ny-1)) ? ny-1 : iy+range_gp;
      for (iz=0; iz<nz; iz++) {
	index = ix*ny*nz + iy*nz + iz;
	if (!(is_mask_bit_set(mask, index))) {
	  integrated_field[ix][iy][iz] = 0.0;
	  irange[2][0] = ((iz-range_gp) < 0) ? 0 : iz-range_gp;
	  irange[2][1] = ((iz+range_gp) > (nz-1)) ? nz-1 : iz+range_gp;
	  for (i=irange[0][0]; i<=irange[0][1]; i++) {
	    for (j=irange[1][0]; j<=irange[1][1]; j++) {
	      for (k=irange[2][0]; k<=irange[2][1]; k++) {
		index = i*ny*nz + j*nz + k;
		if ((! is_mask_bit_set(mask, index)) &&
		    (fabs(score_matrix[i][j][k] - ignore_value) > 1e-10)) {
		  integrated_field[ix][iy][iz] += score_matrix[i][j][k];
		}
	      }
	    }
	  }
	}
      }
    }
  }

  for (ix=0; ix<nx; ix++) {
    for (iy=0; iy<ny; iy++) {
      for (iz=0; iz<nz; iz++) {
	index = ix*ny*nz + iy*nz + iz;
	if (!(is_mask_bit_set(mask, index))) {
	  score_matrix[ix][iy][iz] = integrated_field[ix][iy][iz];
	}
      }
    }
  }
  free_3d_matrix((void***)integrated_field);
}



MAP* copy_fmap(MAP *map) {

  int ix,iy,iz,nx,ny,nz,npts;
  double ***matrix,***cmatrix;
  MAP *cmap;
  GRID *grid;

  if (strcmp(map->type->name,"double")) {

    error_fn("%s: map '%s' not of type double",__func__,map->name);
  }

  grid = map->grid;

  cmap = new_map(map->name,map->type->name,NULL);

  cmap->grid = (GRID*) malloc(sizeof(GRID));

  if (cmap->grid == NULL) {

    error_fn("%s: out of memory allocating grid",__func__);
  }

  memcpy(cmap->grid,grid,sizeof(GRID));

  cmap->type = map->type;

  alloc_map_matrix(cmap);

  nx = grid->npoints[0];
  ny = grid->npoints[1];
  nz = grid->npoints[2];

  npts = (nx*ny*nz);

  matrix = (double***) map->matrix;
  cmatrix = (double***) cmap->matrix;

  for (ix=0;ix<nx;ix++) {

    for (iy=0;iy<ny;iy++) {

      for (iz=0;iz<nz;iz++) {

	cmatrix[ix][iy][iz] = matrix[ix][iy][iz];
      }
    }
  }

  memcpy(cmap->mask,map->mask,((npts/8)+1)*sizeof(unsigned char));

  return(cmap);
}



void inflate_map(MAP *map) {

  int ix,iy,iz,nx,ny,nz,jx,jy,jz,jx1,jx2,jy1,jy2,jz1,jz2,index;
  double max_value,value,***matrix,***cmatrix;
  unsigned char *mask;
  GRID *grid;
  MAP *cmap;
  
  cmap = copy_fmap(map);

  grid = map->grid;

  matrix = (double***) map->matrix;
  cmatrix = (double***) cmap->matrix;

  mask = map->mask;

  nx = grid->npoints[0];
  ny = grid->npoints[1];
  nz = grid->npoints[2];

  for (ix=0;ix<nx;ix++) {

    jx1 = (ix == 0) ? 0 : -1;
    jx2 = (ix+1 == nx) ? 0 : 1;

    for (iy=0;iy<ny;iy++) {

      jy1 = (iy == 0) ? 0 : -1;
      jy2 = (iy+1 == ny) ? 0 : 1;

      for (iz=0;iz<nz;iz++) {

	index = ny*nz*ix + nz*iy + iz;

	//if (!is_mask_bit_set(mask,index)) {

	  jz1 = (iz == 0) ? 0 : -1;
	  jz2 = (iz+1 == nz) ? 0 : 1;

	  max_value = 0.0;

	  for (jx=jx1;jx<=jx2;jx++) {

	    for (jy=jy1;jy<=jy2;jy++) {
	      
	      for (jz=jz1;jz<=jz2;jz++) {

		index = ny*nz*(ix+jx) + nz*(iy+jy) + (iz+jz);

		//if (!is_mask_bit_set(mask,index)) {
		  
		  value = cmatrix[ix+jx][iy+jy][iz+jz];

		  if (value > max_value) {

		    max_value = value;
		  }
		  //}
	      }
	    }
	  }

	  matrix[ix][iy][iz] = max_value;
	  //}
      }
    }
  }

  free_map(cmap);
}



void average_map(MAP *map) {

  int ix,iy,iz,nx,ny,nz,jx,jy,jz,jx1,jx2,jy1,jy2,jz1,jz2,index;
  double value,weight,***matrix,***cmatrix;
  unsigned char *mask;
  GRID *grid;
  MAP *cmap;
  
  cmap = copy_fmap(map);

  grid = map->grid;

  matrix = (double***) map->matrix;
  cmatrix = (double***) cmap->matrix;

  mask = map->mask;

  nx = grid->npoints[0];
  ny = grid->npoints[1];
  nz = grid->npoints[2];

  for (ix=0;ix<nx;ix++) {

    jx1 = (ix == 0) ? 0 : -1;
    jx2 = (ix+1 == nx) ? 0 : 1;

    for (iy=0;iy<ny;iy++) {

      jy1 = (iy == 0) ? 0 : -1;
      jy2 = (iy+1 == ny) ? 0 : 1;

      for (iz=0;iz<nz;iz++) {

	index = ny*nz*ix + nz*iy + iz;

	//if (!is_mask_bit_set(mask,index)) {

	  jz1 = (iz == 0) ? 0 : -1;
	  jz2 = (iz+1 == nz) ? 0 : 1;

	  weight = 0.0;
	  value = 0.0;

	  for (jx=jx1;jx<=jx2;jx++) {

	    for (jy=jy1;jy<=jy2;jy++) {
	      
	      for (jz=jz1;jz<=jz2;jz++) {

		index = ny*nz*(ix+jx) + nz*(iy+jy) + (iz+jz);

		//if (!is_mask_bit_set(mask,index)) {
		  		  
		  value += cmatrix[ix+jx][iy+jy][iz+jz];

		  weight += 1.0;
		  //}
	      }
	    }
	  }

	  matrix[ix][iy][iz] = (value/weight);
	  //}
      }
    }
  }

  free_map(cmap);
}



void smoothe_map(MAP *map) {

  int ix,iy,iz,nx,ny,nz,jx,jy,jz,jx1,jx2,jy1,jy2,jz1,jz2,index;
  double value,weight,wt,***matrix,***cmatrix;
  unsigned char *mask;
  GRID *grid;
  MAP *cmap;
  
  cmap = copy_fmap(map);

  grid = map->grid;

  matrix = (double***) map->matrix;
  cmatrix = (double***) cmap->matrix;

  mask = map->mask;

  nx = grid->npoints[0];
  ny = grid->npoints[1];
  nz = grid->npoints[2];

  for (ix=0;ix<nx;ix++) {

    jx1 = (ix == 0) ? 0 : -1;
    jx2 = (ix+1 == nx) ? 0 : 1;

    for (iy=0;iy<ny;iy++) {

      jy1 = (iy == 0) ? 0 : -1;
      jy2 = (iy+1 == ny) ? 0 : 1;

      for (iz=0;iz<nz;iz++) {

	index = ny*nz*ix + nz*iy + iz;

	if (!is_mask_bit_set(mask,index)) {

	  jz1 = (iz == 0) ? 0 : -1;
	  jz2 = (iz+1 == nz) ? 0 : 1;

	  weight = 0.0;
	  value = 0.0;

	  for (jx=jx1;jx<=jx2;jx++) {

	    for (jy=jy1;jy<=jy2;jy++) {
	      
	      for (jz=jz1;jz<=jz2;jz++) {

		index = ny*nz*(ix+jx) + nz*(iy+jy) + (iz+jz);

		if (!is_mask_bit_set(mask,index)) {
		  
		  wt = abs(jx) + abs(jy) + abs(jz);
		  
		  value += wt*cmatrix[ix+jx][iy+jy][iz+jz];

		  weight += wt;
		}
	      }
	    }
	  }

	  matrix[ix][iy][iz] = 0.5*cmatrix[ix][iy][iz] + 0.5*(value/weight);
	}
      }
    }
  }

  free_map(cmap);
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

  init_mask(map->mask,grid);

  if (strcmp(map->type->name,"double")) {

    error_fn("%s: map '%s' not of type double",__func__,map->name);
  }

  dm = (double***) map->matrix;
  
  n_alloc_islands = 100;
  
  islands = (MAP_ISLAND*) calloc(n_alloc_islands,sizeof(MAP_ISLAND));
  
  if (islands == NULL) {
    
    error_fn("map2islands: out of memory allocating islands");
  }
  
  *n_islands = 0;
  
  while (((index = map_extreme_index(map,1,"max",&ix,&iy,&iz)) != -1) && (dm[ix][iy][iz] >= level)) {
    
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

  return(islands);
}



void free_map_islands(MAP_ISLAND *islands,int n_islands) {

  int i;
  MAP_ISLAND *island;

  if (n_islands == 0) {

    return;
  }

  for (i=0,island=islands;i<n_islands;i++,island++) {

    if (island->n_alloc_points) {

      free(island->points);
    }

    if (island->site_atoms) {

      free_atomlist(island->site_atoms);
    }
  }

  free(islands);
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

  if (strcmp(map->type->name,"double")) {

    error_fn("%s: map '%s' not of type double",__func__,map->name);
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

  if (strcmp(map->type->name,"double")) {

    error_fn("%s: map '%s' not of type double",__func__,map->name);
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



void grow_map_island(MAP_ISLAND *island,MAP *map,double level,int index) {

  int ix,iy,iz;
  double ***dm;
  GRID *grid;
  unsigned char *mask;

  if (strcmp(map->type->name,"double")) {

    error_fn("%s: map '%s' not of type double",__func__,map->name);
  }

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

  dm = (double***) map->matrix;
  
  if (dm[ix][iy][iz] < level) {
    
    island->ix = ix;
    island->iy = iy;
    island->iz = iz;
    
    return;
  }
  
  set_mask_bit(mask,index);
  
  add_point_to_map_island(island,index,0);
  
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
}



void calc_map_island_properties(MAP_ISLAND *island,MAP *map) {

  int i,*index,ix,iy,iz,frag_n_points;
  double ***dm,Y,maxY,sumY,frag_sumY,gridvol,frag_volume;
  GRID *grid;

  if (strcmp(map->type->name,"double")) {

    error_fn("%s: map '%s' not of type double",__func__,map->name);
  }

  if (island->n_points > 0) {

    grid = map->grid;
    
    gridvol = pow(grid->spacing,3.0);
    
    frag_volume = (params_get_parameter("pocket_frag_volume"))->value.d;
    
    frag_n_points = floor(frag_volume/gridvol);
    
    dm = (double***) map->matrix;
    
    sumY = frag_sumY = 0.0;
    maxY = -999.9E10;
    
    for (i=0,index=island->points;i<island->n_points;i++,index++) {
      
      if (!map_index2grid(*index,&ix,&iy,&iz,grid)) {
	
	Y = dm[ix][iy][iz];
	
	if (Y > maxY) {
	  
	  maxY = Y;
	}
	
	sumY += Y;
	
	if (i < frag_n_points) {
	  
	  frag_sumY += Y;
	}
      }
    }
    
    island->maxY = maxY;
    island->sumY = sumY;
    island->avgY = sumY/((double) island->n_points);
    island->integral = sumY*pow(grid->spacing,3);
    island->n_frag_points = frag_n_points;
    island->frag_integral = frag_sumY*pow(grid->spacing,3);
    island->volume = ((double) island->n_points)*pow(grid->spacing,3);
  }
}



void init_map_island(MAP_ISLAND *island,int n_alloc_points) {

  island->id = 0;

  island->n_points = 0;
  island->n_edge_points = 0;
  island->n_alloc_points = n_alloc_points;
  island->points = NULL;


  island->maxY = -999.9E10;
  island->sumY = 0.0;
  island->avgY = 0.0;
  island->integral = 0.0;
  island->n_frag_points = 0;
  island->frag_integral = 0.0;
  island->volume = 0.0;

  island->ix = -1;
  island->iy = -1;
  island->iz = -1;

  if (n_alloc_points) {

    island->points = (int*) calloc(n_alloc_points,sizeof(int));

    if (island->points == NULL) {

      error_fn("%s: out of memory allocating points",__func__);
    }
  }

  island->site_atoms = NULL;
}



void add_point_to_map_island(MAP_ISLAND *island,int index,int edge) {

  int n_points;

  n_points = island->n_points + island->n_edge_points;

  if (n_points == island->n_alloc_points) {

    island->n_alloc_points *= 2;

    island->points = (int*) realloc(island->points,(island->n_alloc_points)*sizeof(int));

    if (island->points == NULL) {

      error_fn("add_point_to_map_island: out of memory re-allocating points");
    }
  }

  // WARNING: this routine does not allow for the possibility of adding normal
  // points when edge points already exist - i.e. it assumes that either
  //   (1) only normal points are used
  //   (2) when a mixture of normal and edge points is used, points are always added as edge points first

  island->points[n_points] = index;

  if (edge) {

    island->n_edge_points++;

  } else {

    island->n_points++;
  }
}



void remove_small_islands(MAP *map,double level,double volume) {

  int i,j,*idx,ix,iy,iz,nx,ny,nz,n_islands;
  double ***matrix,***cmatrix;
  GRID *grid;
  MAP *cmap;
  MAP_ISLAND *islands,*island;
  
  cmap = copy_fmap(map);

  grid = map->grid;

  matrix = (double***) map->matrix;
  cmatrix = (double***) cmap->matrix;

  nx = grid->npoints[0];
  ny = grid->npoints[1];
  nz = grid->npoints[2];
   
  for (ix=0;ix<nx;ix++) {
      
    for (iy=0;iy<ny;iy++) {
	
      for (iz=0;iz<nz;iz++) {
	
	matrix[ix][iy][iz] = 0.0;
      }
    }
  }
 
  islands = map2islands(cmap,level,&n_islands);

  for (i=0,island=islands;i<n_islands;i++,island++) {

    calc_map_island_properties(island,map);
      
    if (island->volume > volume) {

      for (j=0,idx=island->points;j<island->n_points;j++,idx++) {

	if (!map_index2grid(*idx,&ix,&iy,&iz,grid)) {
	  
	  matrix[ix][iy][iz] = cmatrix[ix][iy][iz];
	}
      }
    }
  }
  
  free_map_islands(islands,n_islands);

  free_map(cmap);
}


int map_extreme_index(MAP *map,int use_mask, char *extreme_type, int *mix, int *miy, int *miz) {

  int nx,ny,nz,ix,iy,iz,index,mindex, find_max;
  double ***dm,value,max_value,min_value;
  GRID *grid;
  unsigned char *mask;

  if (strcmp(map->type->name,"double")) {

    error_fn("%s: map '%s' not of type double",__func__,map->name);
  }

  grid = map->grid;
  mask = map->mask;

  nx = grid->npoints[0];
  ny = grid->npoints[1];
  nz = grid->npoints[2];

  max_value = -9.99E10;
  min_value = 9.99E10;

  mindex = -1;

  find_max = (strcmp(extreme_type, "max") == 0) ? 1:0;

  dm = (double***) map->matrix;
  
  for (ix=0;ix<nx;ix++) {
    
    for (iy=0;iy<ny;iy++) {
      
      for (iz=0;iz<nz;iz++) {
	
	index = ny*nz*ix + nz*iy + iz;
	
	value = dm[ix][iy][iz];
	
	if ((use_mask) && (is_mask_bit_set(mask,index))) {
	  continue;
	}
	if (find_max) {
	  if (value > max_value) {
	    mindex = index;
	    *mix = ix; *miy = iy; *miz = iz;
	    max_value = value;
	  }
	} else if (value < min_value){
	  mindex = index;
	  *mix = ix; *miy = iy; *miz = iz;
	  min_value = value;
	}
      }
    }
  }

  return(mindex);
}

void print_grid(GRID *grid) {
  int i;
  printf("spacing: %.2f\n", grid->spacing);
  printf("padding: %.2f\n", grid->padding);
  for (i=0; i<3; i++) {
    printf ("%d: grid_n_points=%d | grid_limit=%d, %d | grid_flimit = %.2f, %.2f\n",
            i, grid->npoints[i], grid->limit[i][0], grid->limit[i][1],
            grid->flimit[i][0], grid->flimit[i][1]);
  }
}

void mask_map_cutoff(MAP *map, double cutoff, char *side) {

  if (strcmp(map->type->name,"double")) {

    error_fn("%s: map '%s' is not of type double",__func__,map->name);
  }

  int nx,ny,nz,ix,iy,iz,index;
  int mask_above;
  GRID *grid;
  unsigned char *mask;
  double ***matrix;

  grid = map->grid;
  mask = map->mask;

  nx = grid->npoints[0];
  ny = grid->npoints[1];
  nz = grid->npoints[2];

  mask_above = (! strcmp(side, "above")) ? 1: 0;



  matrix = (double***) map->matrix;

  for (ix=0;ix<nx;ix++) {

    for (iy=0;iy<ny;iy++) {

      for (iz=0;iz<nz;iz++) {
        index = ny*nz*ix + nz*iy + iz;
        if (! is_mask_bit_set(mask, index)) {
          if (mask_above) {
            if (matrix[ix][iy][iz] > cutoff) {
              set_mask_bit(mask, index);
            }
          } else if (matrix[ix][iy][iz] < cutoff) {
            set_mask_bit(mask, index);
          }
        }
      }
    }
  }
}



void v_init_double(void *value) {

  *((double*) value) = 0.0;
}



void v_init_int(void *value) {

  *((int*) value) = 0;
}



void v_init_atomlist(void *list) {

  init_atomlist((ATOMLIST*) list);
}



void v_init_alist(void *vlist) {

  LIST *list;

  list = (LIST*) vlist;

  init_list(list,"",sizeof(ATOM*));
}



void v_free_list_items(void *vlist) {

  LIST *list;
  void *items;

  list = (LIST*) vlist;

  items = list->items;

  if (items) {
    
    free(items);
  }
}
