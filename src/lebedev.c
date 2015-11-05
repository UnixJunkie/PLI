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



double** read_lebedev_sphere_points(int n_points) {

  int i,flag,order,n_words;
  char *pli_dir,filename[MAX_LINE_LEN],line[MAX_LINE_LEN],word[MAX_LINE_LEN];
  double **points;
  PLI_FILE *file;

  pli_dir = get_pli_dir();

  sprintf(filename,"%s/params/lebedev.pli",pli_dir);

  file = open_file(filename,"r");

  if (file == NULL) {

    error_fn("read_lebedev_sphere_points: couldn't open file '%s'",filename);
  }

  points = NULL;

  while (!end_of_file(file)) {

    if (read_line(line,MAX_LINE_LEN,file) == NULL)
      break;

    flag = sscanf(line,"%s",word);

    if (flag != EOF) {

      if (!strcmp(word,"Order")) {

	sscanf(line,"%*s %*s %d",&order);

	if (order == n_points) {

	  points = (double**) alloc_2d_fmatrix(n_points,5);

	  for (i=0;i<n_points;i++) {

	    if (read_line(line,MAX_LINE_LEN,file) == NULL) {

	      error_fn("read_lebedev_sphere_points: not enough points");
	    }

	    n_words = sscanf(line,"%lf %lf %lf %lf %lf",&(points[i][0]),&(points[i][1]),&(points[i][2]),&(points[i][3]),&(points[i][4]));

	    if (n_words != 5) {

	      error_fn("read_lebedev_sphere_points: corrupted line:\n%s",line);
	    }
	  }
	}
      }
    }
  }

  close_file(file);

  if (points == NULL) {

    error_fn("read_lebedev_sphere_points: no sphere points for order=%d",n_points);
  }

  return(points);
}
