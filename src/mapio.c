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




void write_insight_map(char *filename,MAP *map,int save_mask) { 
  
  GRID *grid;
  unsigned char *mask;
  PLI_FILE *file;
  int nx,ny,nz,ix,iy,iz,index;
  int ***imatrix;
  double ***fmatrix;

  if ((file = open_file(filename,"w")) == NULL) {

    error_fn("write_insight_map: failed to open file %s",filename);
  }

  grid = map->grid;

  write_line(file,"%s\n",map->name);

  write_line(file,"(1F10.3)\n");

  write_line(file," %9.3lf %9.3lf %9.3lf %9.3lf %9.3lf %9.3lf\n",
	     grid->spacing*(grid->limit[0][1]-grid->limit[0][0]),
	     grid->spacing*(grid->limit[1][1]-grid->limit[1][0]),
	     grid->spacing*(grid->limit[2][1]-grid->limit[2][0]),
	     90.0,90.0,90.0);

  write_line(file," %9d %9d %9d\n",
	     grid->npoints[0]-1,
	     grid->npoints[1]-1,
	     grid->npoints[2]-1);

  write_line(file," %9d %9d %9d %9d %9d %9d %9d\n",
	     1,
	     grid->limit[0][0],grid->limit[0][1],
	     grid->limit[1][0],grid->limit[1][1],
	     grid->limit[2][0],grid->limit[2][1]);
    
  nx = grid->npoints[0];
  ny = grid->npoints[1];
  nz = grid->npoints[2];

  if (map->type == INTEGER_MAP) {

    imatrix = (int***) map->matrix;

  } else if (map->type == DOUBLE_MAP) {

    fmatrix = (double***) map->matrix;

  } else {

    error_fn("write_insight_map: cannot save this map type");
  }

  mask = map->mask;

  for (iz=0;iz<nz;iz++) {

    for (iy=0;iy<ny;iy++) {

      for (ix=0;ix<nx;ix++) {

	if (save_mask) {

	  index = ny*nz*ix + nz*iy + iz;

	  if (is_mask_bit_set(mask,index)) {

	    write_line(file,"%10.3lf\n",1.0);

	  } else {

	    write_line(file,"%10.3lf\n",0.0);
	  }

	} else {

	  if (map->type == DOUBLE_MAP) {

	    write_line(file,"%10.3lf\n",fmatrix[ix][iy][iz]);

	  } else {

	    write_line(file,"%10.3lf\n",(double) imatrix[ix][iy][iz]);
	  }
	}
      }
    }
  }

  close_file(file);
}



MAP* read_insight_map(char *filename) { 
 
  PLI_FILE *file;
  char line[MAX_LINE_LEN],title[MAX_LINE_LEN];
  double box_size[3],box_angle[3],***dm;
  int i,ix,iy,iz;
  int fastest_value;
  GRID *grid;
  MAP *map;
  
  if ((file = open_file(filename,"r")) == NULL) {

    error_fn("read_insight_map: failed to open file %s",filename);
  }

  if (read_line(line,MAX_LINE_LEN,file) == NULL) {

    error_fn("read_insight_map: failed to read title from %s",filename);
  }

  sscanf(line,"%[^\n]",title);

  if (read_line(line,MAX_LINE_LEN,file) == NULL) {

    error_fn("read_insight_map: failed to read format from %s",filename);
  }

  map = new_map(title);

  map->type = DOUBLE_MAP;

  map->grid = new_grid();

  grid = map->grid;

  if (read_line(line,MAX_LINE_LEN,file) == NULL) {

    error_fn("read_insight_map: failed to read box size and angles from %s",filename);
  }

  if (sscanf(line," %lf %lf %lf %lf %lf %lf",
	     &box_size[0],&box_size[1],&box_size[2],&box_angle[0],&box_angle[1],&box_angle[2]) != 6) {

    error_fn("read_insight_map: failed to read box size and angles from %s",filename);    
  }
  
  if (read_line(line,MAX_LINE_LEN,file) == NULL) {

    error_fn("read_insight_map: failed to read step numbers from %s",filename);
  }

  if (sscanf(line,"%d %d %d",grid->npoints,grid->npoints+1,grid->npoints+2) != 3) {

    error_fn("read_insight_map: failed to read step numbers from %s",filename);
  }

  if (read_line(line,MAX_LINE_LEN,file) == NULL) {

    error_fn("read_insight_map: failed to read grid limits from %s",filename);
  }

  if (sscanf(line,"%d %d %d %d %d %d %d",
	     &fastest_value,
	     &grid->limit[0][0],&grid->limit[0][1],
	     &grid->limit[1][0],&grid->limit[1][1],
	     &grid->limit[2][0],&grid->limit[2][1]) != 7) {

    error_fn("read_insight_map: failed to read grid limits from %s",filename);
  }

  if (fastest_value != 1) {

    error_fn("read_insight_map: failed to readi grid limits from %s",filename);
  }

  grid->spacing = 0.0;

  for (i=0;i<3;i++) {

    grid->spacing += (box_size[i]/((double) grid->npoints[i]));

    grid->npoints[i]++;
  }

  grid->spacing /= 3.0;

  for (i=0;i<3;i++) {

    grid->flimit[i][0] = grid->spacing*((double) grid->limit[i][0]);
    grid->flimit[i][1] = grid->spacing*((double) grid->limit[i][1]);
  }

  alloc_map_matrix(map);
  
  dm = (double***) map->matrix;

  for (iz=0;iz<grid->npoints[2];iz++) {

    for (iy=0;iy<grid->npoints[1];iy++) {

      for (ix=0;ix<grid->npoints[0];ix++) {

	if (read_line(line,MAX_LINE_LEN,file) == NULL) {

	  error_fn("read_insight_map: failed to read map value from %s",filename);
	}

	if (sscanf(line,"%lf",&(dm[ix][iy][iz])) != 1) {

	  error_fn("read_insight_map: failed to read map value from %s",filename);
	}
      }
    }
  }

  close_file(file);

  return(map);
}

