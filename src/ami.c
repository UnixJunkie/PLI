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



MAP* ami_site_map(SYSTEM *system) {

  int i,ix,iy,iz,nx,ny,nz,index;
  BYTE *volume_mask,*clash_mask;
  double ***matrix,pos[4],area,value,dr,probe_area,radius;
  FIELD_PROBE *probe;
  ATOM *atom1,*atom2;
  CONTACTLIST *contactlist;
  CONTACT *contact;
  HISTOGRAM *ami_hist;
  GRID *grid;
  MAP *map;

  reset_system(system,0);

  calc_system_simb(system,1);

  grid = system_grid(system,(params_get_parameter("grid_spacing"))->value.d,(params_get_parameter("grid_padding"))->value.d);

  map = new_map("amisite","double",grid);

  probe = get_field_probe("H2O");

  probe->type = ATOM_PROBE;

  radius = 2.0*(probe->atom->vdw_radius_H2O) - (probe->atom->vdw_radius);

  volume_mask = mask_system_volume(system,NULL,grid,PROBE_MASK,radius);
  clash_mask = mask_system_atoms(system,NULL,grid,CONTACT_MASK,0.0);

  matrix = (double***) map->matrix;

  init_field_system(system,probe);

  atom1 = probe->atom;

  dr = atom1->vdw_radius_H2O - atom1->vdw_radius;

  atom1->vdw_radius = (params_get_parameter("amisite_probe_radius"))->value.d;
  atom1->vdw_radius_H2O = atom1->vdw_radius + dr;

  probe_area = 4.0*PI*sqr(atom1->vdw_radius_H2O);

  pos[3] = 1.0;

  ami_hist = get_histogram(system->settings->force_field->AMI_HISTOGRAM_ID);

  nx = grid->npoints[0];
  ny = grid->npoints[1];
  nz = grid->npoints[2];

  for (ix=0;ix<nx;ix++) {

    pos[0] = grid->flimit[0][0] + ((double) ix)*(grid->spacing);

    for (iy=0;iy<ny;iy++) {

      pos[1] = grid->flimit[1][0] + ((double) iy)*(grid->spacing);

      for (iz=0;iz<nz;iz++) {

	index = ny*nz*ix + nz*iy + iz;

	if ((!(is_mask_bit_set(clash_mask,index))) && (is_mask_bit_set(volume_mask,index))) {

	  value = 0.0;
	    
	  pos[2] = grid->flimit[2][0] + ((double) iz)*(grid->spacing);
	    
	  move_probe_to_new_position(probe,pos);
	    
	  set_contacts_atom(atom1,system,0);
	    
	  contactlist = atom1->contactlist;
	    
	  if (contactlist != NULL) {
	      
	    for (i=0,contact=contactlist->contacts;i<contactlist->ncontacts;i++,contact++) {

	      atom2 = contact->atom2;
		
	      area = contact->area;
		
	      value += area*histogram_X2Y(ami_hist,atom2->tdf);
	    }
	  }

	  value += (atom1->exposed_area)*histogram_X2Y(ami_hist,2.8);

	  matrix[ix][iy][iz] = value/probe_area;

	} else {
	  
	  matrix[ix][iy][iz] = 0.0;
	}
      }
    }
  }

  finish_field_system(system,probe);

  free(volume_mask);
  free(clash_mask);

  return(map);
}
