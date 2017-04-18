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



static BYTE* mask_sphere(double*,BYTE*,GRID*,double);
static double mask_dmax(ATOM*,int,double);
static int is_byte_bit_set(BYTE*,int);
static void set_byte_bit(BYTE*,int);
static void unset_byte_bit(BYTE*,int);



int mask_bytes(GRID *grid) {

  int nx,ny,nz;

  nx = grid->npoints[0];
  ny = grid->npoints[1];
  nz = grid->npoints[2];

  return((nx*ny*nz)/8 + 1);
}



BYTE* mask_system_atoms(SYSTEM *system,BYTE *mask,GRID *grid,int type,double dmax) {

  int i;
  MOLECULE **molecule;

  mask = (mask) ? mask : alloc_mask(grid);

  for (i=0,molecule=system->molecule_list->molecules;i<system->molecule_list->n_molecules;i++,molecule++) {

    mask_atoms((*molecule)->nonh_atoms,mask,grid,type,dmax);
  }
  
  return(mask);
}



BYTE* mask_system_volume(SYSTEM *system,BYTE *mask,GRID *grid,int type,double dmax) {

  int i;
  BYTE *site_mask;
  MOLECULE **molecule;
  MAP *site_map;

  mask = (mask) ? mask : alloc_mask(grid);

  for (i=0,molecule=system->molecule_list->molecules;i<system->molecule_list->n_molecules;i++,molecule++) {

    mask_atoms((*molecule)->active_atoms,mask,grid,type,dmax);
  }
  
  if (system->site_ligand) {

    site_mask = mask_atoms(system->site_ligand->nonh_atoms,NULL,grid,FIXED_MASK,(params_get_parameter("site_radius"))->value.d);

    and_mask(mask,site_mask,mask,grid);

    free(site_mask);

    site_map = mask_map(mask,NULL,grid);

    //write_insight_map("site_map.grd",site_map,0);

    free_map(site_map);
  }

  return(mask);
}



// TODO: this routine should be retired and replaced with the mask_atoms() one

BYTE* mask_molecule(MOLECULE *molecule,BYTE *mask,GRID *grid,int type) {

  int i;
  double dmax;
  ATOM *atom;

  mask = (mask) ? mask : alloc_mask(grid);

  dmax = (params_get_parameter("custom_mask_dmax"))->value.d;

  for (i=0,atom=molecule->atom;i<molecule->natoms;i++,atom++) {

    mask_atom(atom,mask,grid,type,dmax);
  }

  return(mask);
}



BYTE* mask_atoms(LIST *list,BYTE *mask,GRID *grid,int type,double dmax) {

  int i;
  ATOM **atom;

  mask = (mask) ? mask : alloc_mask(grid);

  for (i=0,atom=(ATOM**) list->items;i<list->n_items;i++,atom++) {

    mask_atom(*atom,mask,grid,type,dmax);
  }

  return(mask);
}



BYTE *mask_atom(ATOM *atom,BYTE *mask,GRID *grid,int type,double dmax) {

  mask = (mask) ? mask : alloc_mask(grid);

  if ((!(atom->flags & SKIP_ATOM)) && (atom->element->id != HYDROGEN)) {

    dmax = mask_dmax(atom,type,dmax);

    mask_sphere(atom->position,mask,grid,dmax);
  }

  return(mask);
}



static BYTE* mask_sphere(double *v,BYTE *mask,GRID *grid,double dmax) {

  int i,imax,iv[3],iv_min[3],iv_max[3],ny,nz,ix,iy,iz,index;
  double sqr_dmax,dx,dy,dz;

  mask = (mask) ? mask : alloc_mask(grid);

  if (point_in_grid(v,grid,dmax)) {

    point_to_gridpoint(v,iv,grid);
    
    imax = ceil(dmax/grid->spacing);
    
    for (i=0;i<3;i++) {
      
      iv_min[i] = (iv[i] - imax < 0) ? 0 : iv[i] - imax;
      iv_max[i] = (iv[i] + imax >= grid->npoints[i]) ? grid->npoints[i] - 1 : iv[i] + imax;
    }
  
    ny = grid->npoints[1];
    nz = grid->npoints[2];
  
    sqr_dmax = dmax*dmax;
    
    for (ix=iv_min[0];ix<=iv_max[0];ix++) {
      
      dx = grid->flimit[0][0] + ((double) ix)*(grid->spacing) - v[0];
      
      for (iy=iv_min[1];iy<=iv_max[1];iy++) {
	
	dy = grid->flimit[1][0] + ((double) iy)*(grid->spacing) - v[1];
	
	for (iz=iv_min[2];iz<=iv_max[2];iz++) {
	  
	  dz = grid->flimit[2][0] + ((double) iz)*(grid->spacing) - v[2];
	  
	  index = ny*nz*ix + nz*iy + iz;

	  if ((!is_mask_bit_set(mask,index)) && (((dx*dx)+(dy*dy)+(dz*dz)) < sqr_dmax)) {

	    set_mask_bit(mask,index);
	  }
	}
      }  
    }
  }

  return(mask);
}



BYTE* alloc_mask(GRID *grid) {

  int n_bytes;
  BYTE *mask;

  n_bytes = mask_bytes(grid);

  mask = (BYTE*) calloc(n_bytes,sizeof(BYTE));

  if (mask == NULL) {

    error_fn("%s: out of memory allocating mask",__func__);
  }

  memset(mask,0,n_bytes);

  return(mask);
}



void init_mask(BYTE *mask,GRID *grid) {

  int n_bytes;

  n_bytes = mask_bytes(grid);

  memset(mask,0,n_bytes);
}



BYTE* copy_mask(BYTE *to_mask,BYTE *from_mask,GRID *grid) {

  int n_bytes;

  to_mask = (to_mask) ? to_mask : alloc_mask(grid);

  n_bytes = mask_bytes(grid);

  memcpy(to_mask,from_mask,n_bytes);
}



BYTE* invert_mask(BYTE *mask1,BYTE *mask,GRID *grid) {

  int i,n_bytes,bit;
  BYTE *byte1,*byte;

  mask = (mask) ? mask : alloc_mask(grid);  

  n_bytes = mask_bytes(grid);

  for (i=0,byte1=mask1,byte=mask;i<n_bytes;i++,byte1++,byte++) {

    for (bit=0;bit<8;bit++) {

      if (is_byte_bit_set(byte1,bit)) {

	unset_byte_bit(byte,bit);

      } else {

	set_byte_bit(byte,bit);
      }
    }
  }

  return(mask);
}


BYTE* and_mask(BYTE *mask1,BYTE *mask2,BYTE *mask,GRID *grid) {

  int i,n_bytes,bit;
  BYTE *byte1,*byte2,*byte;

  mask = (mask) ? mask : alloc_mask(grid);

  n_bytes = mask_bytes(grid);

  for (i=0,byte1=mask1,byte2=mask2,byte=mask;i<n_bytes;i++,byte1++,byte2++,byte++) {

    for (bit=0;bit<8;bit++) {

      if ((is_byte_bit_set(byte1,bit)) && (is_byte_bit_set(byte2,bit))) {

	set_byte_bit(byte,bit);

      } else {

	unset_byte_bit(byte,bit);
      }
    }
  }

  return(mask);
}



BYTE* or_mask(BYTE *mask1,BYTE *mask2,BYTE *mask,GRID *grid) {

  int i,n_bytes,bit;
  BYTE *byte1,*byte2,*byte;

  mask = (mask) ? mask : alloc_mask(grid);

  n_bytes = mask_bytes(grid);

  for (i=0,byte1=mask1,byte2=mask2,byte=mask;i<n_bytes;i++,byte1++,byte2++,byte++) {

    for (bit=0;bit<8;bit++) {

      if ((is_byte_bit_set(byte1,bit)) || (is_byte_bit_set(byte2,bit))) {

	set_byte_bit(byte,bit);

      } else {

	unset_byte_bit(byte,bit);
      }
    }
  }

  return(mask);
}


void set_mask_bit(BYTE *mask,int index) {

  int minor,major;

  major = index / 8;
  minor = index % 8;

  set_byte_bit(mask+major,minor);
}


void unset_mask_bit(BYTE *mask,int index) {

  int minor,major;

  major = index / 8;
  minor = index % 8;

  unset_byte_bit(mask+major,minor);
}


int is_mask_bit_set(BYTE *mask,int index) {

  int minor,major;
  BYTE *byte,value;

  major = index / 8;

  byte = mask + major;

  value = *byte;

  if (!value) {

    return(0);
  }

  minor = index % 8;

  if (value & 01<<minor) {

    return(1);
  }

  return(0);
}



MAP* mask_map(BYTE *mask,MAP *map,GRID *grid) {

  int nx,ny,nz,ix,iy,iz;
  double ***matrix;

  map = (map) ? map : new_map("mask","double",NULL);

  map->grid = copy_grid(map->grid,grid);

  alloc_map_matrix(map);
    
  matrix = (double***) map->matrix;

  nx = grid->npoints[0];
  ny = grid->npoints[1];
  nz = grid->npoints[2];

  for (ix=0;ix<nx;ix++) {

    for (iy=0;iy<ny;iy++) {

      for (iz=0;iz<nz;iz++) {

	if (is_mask_bit_set(mask,ny*nz*ix + nz*iy + iz)) {

	  matrix[ix][iy][iz] = 1.0;
	}
      }
    }
  }

  return(map);
}



// TODO: tidy this routine up - do we really need all these types?

static double mask_dmax(ATOM *atom,int type,double dmax) {

  if (type == FIXED_MASK) {
    
    return(dmax);
    
  } else if (type == PROBE_MASK) {
    
    return(atom->vdw_radius + dmax);
    
  } else if (type == VDW_MASK) {
    
    return((atom->element->flags & METAL_ELEMENT) ? 2.0 : atom->vdw_radius_H2O);
    
  } else {
    
    if (atom->type->flags & ANY_HBOND_ATOM_TYPE) {
      
      return((type == CONTACT_MASK) ? 2.0 : 1.2);
      
    } else if (atom->element->flags & METAL_ELEMENT) {
      
      return((type == CONTACT_MASK) ? 1.5 : 0.8);
      
    } else {
      
      return((type == CONTACT_MASK) ? atom->vdw_radius + 0.8 : atom->vdw_radius - 0.2);
    }
  }
}



static int is_byte_bit_set(BYTE *byte,int bit) {

  BYTE value;

  value = *byte;

  if (!value) {

    return(0);
  }

  if (value & 01<<bit) {

    return(1);
  }

  return(0);
}



static void set_byte_bit(BYTE *byte,int bit) {

  *byte |= 01<<bit;
}



static void unset_byte_bit(BYTE *byte,int bit) {

  *byte &= ~(1 << bit);
}
