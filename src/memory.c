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



#define MEM_BLOCK_SIZE 1000

typedef struct MemBlock MEM_BLOCK;

struct MemBlock {
  MEM_BLOCK *block;
  void *address;
  int position;
  int size;
};



static LIST *mem_blocks = NULL;
static LIST *mem_free = NULL;
static LIST *mem_claimed = NULL;



static MEM_BLOCK* mem_find_free(int);
static MEM_BLOCK* mem_find_claimed(void*);
static MEM_BLOCK* mem_add_block(int);
static void mem_add_free(MEM_BLOCK*);
static void mem_claim_block(MEM_BLOCK*,int);
static void mem_release_block(MEM_BLOCK*);



void* mem_claim(int size) {

  MEM_BLOCK *block;

  if (!mem_blocks) {

    mem_blocks = new_list("memory blocks",sizeof(MEM_BLOCK),0);
  }

  block = mem_find_free(size);

  if (!block) {

    block = (size > MEM_BLOCK_SIZE) ? mem_add_block(size) : mem_add_block(MEM_BLOCK_SIZE);
  }

  mem_claim_block(block,size);

  return(block->address);
}



void mem_release(void *address) {

  MEM_BLOCK *block;

  block = mem_find_claimed(address);

  if (!block) {

    error_fn("%s: could not find claimed memory",__func__);
  }

  mem_release_block(block);
}



static MEM_BLOCK* mem_find_free(int size) {

  int i;
  MEM_BLOCK *block;

  for (i=0,block=(MEM_BLOCK*)mem_free->items;i<mem_free->n_items;i++,block++) {

    if (block->size >= size) {

      return(block);
    }
  }

  return(NULL);
}



static MEM_BLOCK* mem_find_claimed(void *address) {

  int i;
  MEM_BLOCK *block;

  if (mem_claimed) {

    for (i=0,block=(MEM_BLOCK*)mem_claimed->items;i<mem_claimed->n_items;i++,block++) {

      if (block->address == address) {

	return(block);
      }
    }
  }

  return(NULL);
}



static MEM_BLOCK* mem_add_block(int size) {

  MEM_BLOCK *block;

  block = (MEM_BLOCK*) add_list_item(mem_blocks);

  block->address = (void*) malloc(size);

  if (!block->address) {

    error_fn("%s: out of memory allocating block of size %d",__func__,size);
  }

  block->block = block;
  block->size = size;

  mem_add_free(block);

  return(block);
}



static void mem_add_free(MEM_BLOCK *block) {

  MEM_BLOCK *item;

  if (!mem_free) {

    mem_free = new_list("memory free",sizeof(MEM_BLOCK),0);
  }

  item = (MEM_BLOCK*) add_list_item(mem_free);

  item->block = block;
  item->address = block->address;
  item->position = 0;
  item->size = block->size;
}



static void mem_claim_block(MEM_BLOCK *block,int size) {

  MEM_BLOCK *item;

  if (!mem_claimed) {

    mem_claimed = new_list("memory claimed",sizeof(MEM_BLOCK),0);
  }

  item = add_list_item(mem_claimed);

  item->block = block->block;
  item->address = block->address;
  item->position = block->position;
  item->size = size;
  
  if (block->size > size) {

    block->position += size;
    block->address = (void*) (((char*) block->address) + size);
    block->size -= size;

  } else {
    
    remove_list_item(mem_free,block);
  }
}



static void mem_release_block(MEM_BLOCK *block) {

}



int** alloc_2d_imatrix(int nrows, int ncolumns) {

  int **matrix;
  int i;

  matrix = (int**) calloc(nrows,sizeof(int*));

  if (matrix == NULL) {

    error_fn("alloc_2d_imatrix: out of memory allocating matrix");
  }

  matrix[0] = (int*) calloc(nrows*ncolumns,sizeof(int));

  if (matrix[0] == NULL) {

    error_fn("alloc_2d_imatrix: out of memory allocating matrix content");
  }
 
  for (i=1;i<nrows;i++) {

    matrix[i] = matrix[i-1] + ncolumns;
  }

  return matrix;
}


void free_3d_int_matrix (int ***matrix) {

  if (matrix) {

    free(**matrix);
    free(*matrix);
    free(matrix);
  }
}


double** alloc_2d_fmatrix(int nrows, int ncolumns) {

  double **matrix;
  int i;

  matrix = (double**) calloc(nrows,sizeof(double*));

  if (matrix == NULL) {

    error_fn("alloc_2d_fmatrix: out of memory allocating matrix");
  }

  matrix[0] = (double*) calloc(nrows*ncolumns,sizeof(double));

  if (matrix[0] == NULL) {

    error_fn("alloc_2d_fmatrix: out of memory allocating matrix content");
  }
 
  for (i=1;i<nrows;i++) {

    matrix[i] = matrix[i-1] + ncolumns;
  }

  return matrix;
}

BOOLEAN **alloc_2d_bool_matrix (int n_rows, int n_cols)
{
  BOOLEAN **matrix;
  int i;

  matrix = (BOOLEAN **) calloc (n_rows, sizeof (BOOLEAN *));
  if (matrix == NULL) {
    error_fn("%s: out of memory allocating matrix", __func__);
  }
  matrix[0] = (BOOLEAN *) calloc(n_rows*n_cols, sizeof (BOOLEAN));
  if (matrix[0] == NULL)
    error_fn ("%s: out of memory allocating matrix contents", __func__);

  for (i=1; i<n_rows; i++)
    matrix[i] = matrix[i-1] + n_cols;

  return matrix;
}


void free_2d_imatrix (int **matrix) {

  if (matrix) {

    free(*matrix);
    free(matrix);
  }
}



void free_2d_fmatrix(double **matrix) {

  if (matrix) {

    free(*matrix);
    free(matrix);
  }
}


void free_3d_matrix(void ***matrix) {
  if (matrix) {
    free(**matrix);
    free(*matrix);
    free(matrix);
  }
}


void free_2d_bool_matrix (BOOLEAN **m)
{
  free (*m);
  free (m);
}

void copy_2d_bool_matrix (BOOLEAN **to_matrix, BOOLEAN **from_matrix, int nrows, int ncols)
{
  memcpy (*to_matrix, *from_matrix, sizeof (BOOLEAN) * nrows * ncols);
}

BOOLEAN* alloc_bool_array(int n) {
  BOOLEAN *arr;

  arr = (BOOLEAN*) calloc (n, sizeof (BOOLEAN));
  if (arr == NULL) {
    error_fn ("%s: out of memory allocating array", __func__);
  }

  return arr;
}

void copy_bool_array(BOOLEAN *to_array, BOOLEAN *from_array, int n) {
  memcpy(to_array, from_array, sizeof(BOOLEAN));
}


int*** alloc_3d_int_matrix(int nx,int ny,int nz) {

  int i,j,address = 0;
  int ***matrix;

  matrix = (int***) calloc(nx,sizeof(int**));

  if (matrix == NULL)
    return(NULL);

  matrix[0] = (int**) calloc(nx*ny, sizeof(int*));

  if (matrix[0] == NULL)
    return(NULL);

  for (i=1;i<nx;i++)
    matrix[i] = matrix[i-1] + ny;

  matrix[0][0] = (int*) calloc(nx*ny*nz,sizeof(int));

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


void copy_3d_double_matrix (double ***source, double ***target, int *dimension) {
  int ix, iy, iz;
  for (ix=0; ix<dimension[0]; ix++) {
    for (iy=0; iy<dimension[1]; iy++) {
      for (iz=0; iz<dimension[2]; iz++) {
	target[ix][iy][iz] = source[ix][iy][iz];
      }
    }
  }
}

// TODO: use memcpy
void copy_3d_int_matrix (int ***source, int ***target, int *dimension) {
  int ix, iy, iz;
  for (ix=0; ix<dimension[0]; ix++) {
    for (iy=0; iy<dimension[1]; iy++) {
      for (iz=0; iz<dimension[2]; iz++) {
	target[ix][iy][iz] = source[ix][iy][iz];
      }
    }
  }
}

void copy_2d_int_matrix(int **source, int **target, int nrows, int ncols) {
  int i, j;
  for (i=0; i<nrows; i++) {
    for (j=0; j<ncols; j++) {
      target[i][j] = source[i][j];
    }
  }
}




char** alloc_2d_cmatrix(int nrows, int ncolumns) {

  char **matrix;
  int i;

  matrix = (char**) calloc(nrows,sizeof(int*));

  if (matrix == NULL) {

    error_fn("%s: out of memory allocating matrix",__func__);
  }

  matrix[0] = (char*) calloc(nrows*ncolumns,sizeof(int));

  if (matrix[0] == NULL) {

    error_fn("%s: out of memory allocating matrix content",__func__);
  }
 
  for (i=1;i<nrows;i++) {

    matrix[i] = matrix[i-1] + ncolumns;
  }

  return matrix;
}



void free_2d_cmatrix(char **matrix) {

  if (matrix) {

    free(*matrix);
    free(matrix);
  }
}

void free_ring (RING *ring) {
  if (ring == NULL) {
    return;
  }

  ATOM * atom;

  for (int i=0; i<ring->size; i++) {
    atom = ring->atom[i];
    atom->ring = NULL;
  }
  free(ring);
}
