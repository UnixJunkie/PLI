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
