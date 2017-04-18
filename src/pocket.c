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



static void pocket_set_oflags(void);
static void pocket_mask_protein(SYSTEM*,MAP*);
static MAP_ISLAND* pocket_map2sites(MAP*,int*);
static void pocket_calc_site_atoms(MAP_ISLAND*,SYSTEM*,GRID*);
static void pocket_grow_site(MAP*,MAP_ISLAND*,int);
static void pocket_site2map(MAP_ISLAND*,MAP*);
static void pocket_mask_site(MAP*,MAP_ISLAND*);
static int pocket_point_in_list(int,int*,int);
static void pocket_output(MAP*,MAP_ISLAND*,int);
static void pocket_output_map(MAP*);
static void pocket_output_sites(MAP*,MAP_ISLAND*,int);
static void pocket_output_site(MAP*,MAP_ISLAND*);
static void pocket_output_site_atoms(MAP_ISLAND*);
static void pocket_output_site_residues(MAP_ISLAND*);
static void pocket_output_site_points(MAP*,MAP_ISLAND*);
static void pocket_output_site_map(MAP*,MAP_ISLAND*);



static POCKET_MODE pocket_modes[] = { { "ligsite",      ligsite_map },
				      { "amisite",      ami_site_map },
				      { "digsite",      digsite_map },
				      { "last",         NULL } };



static char *pocket_oflags;
static int pocket_calc_sites;
static double pocket_min_integral;



void pocket_run(SETTINGS *settings) {

  int i,n_sites;
  char *output,filename[MAX_LINE_LEN],name[MAX_LINE_LEN];
  SYSTEM *system;
  POCKET_MODE *mode;
  MAP *map,*sites_map;
  MAP_ISLAND *sites,*site;
  
  pocket_set_oflags();

  system = settings2system(settings);

  prep_system(system,ANY_MOLECULE);

  mode = pocket_mode(NULL);

  sprintf(name,"%s_min_integral",mode->name);

  pocket_min_integral = (params_get_parameter(name))->value.d;

  // generate map:

  map = mode->map(system);

  // generate sites:

  if (pocket_calc_sites) {

    pocket_mask_protein(system,map);

    sites = pocket_map2sites(map,&n_sites);

    for (i=0,site=sites;i<n_sites;i++,site++) {

      pocket_calc_site_atoms(site,system,map->grid);
    }
  }

  // output requested data:

  pocket_output(map,sites,n_sites);

  // admin:

  free_map(map);

  unprep_system(system,ANY_MOLECULE);

  free_system(system);
}



double pocket_score(SYSTEM *system) {

  int i;
  double score;
  ATOM **atomp,*atom;
  ATOMLIST *list;

  list = system->selection;

  score = 0.0;

  for (i=0,atomp=list->atom;i<list->natoms;i++,atomp++) {

    atom = *atomp;

    atom->score = atom_pocket_score(atom);

    score += atom->score;
  }

  system->score = score;

  return(score);
}



POCKET_MODE* pocket_mode(char *mode_name) {

  char *name;
  POCKET_MODE *mode;

  name = (mode_name) ? mode_name : (params_get_parameter("pocket_mode"))->value.s;

  mode = pocket_modes;

  while (strcmp(mode->name,"last")) {

    if (!strcmp(mode->name,name)) {

      return(mode);
    }

    mode++;
  }

  error_fn("%s: unknown pocket mode '%s'",__func__,name);
}



void set_atom_depth_maps(SYSTEM *system) {

  int i;
  ATOM **atomp,*atom;
  MAP *depth_map;
  ATOMLIST *selection;

  depth_map = (system->protein) ? system->protein->depth_map : NULL;

  selection = system->selection;

  for (i=0,atomp=selection->atom;i<selection->natoms;i++,atomp++) {

    atom = *atomp;

    if (atom->molecule->flags & LIGAND_MOLECULE) {

      atom->depth_map = depth_map;
 
    } else {

      atom->depth_map = NULL;
    }
  }
}



double atom_pocket_score(ATOM *atom) {

  double *pos,score;
  MAP *map;

  map = atom->depth_map;

  if (map == NULL) {

    return(0.0);
  }

  pos = atom->position;

  score = trilinear_interpolate(map,pos,1.0E30);

  return(score);
}



double atom_depth_weight(ATOM *atom) {

  double *pos,weight;
  MAP *map;

  map = atom->depth_map;

  if (map == NULL) {

    return(1.0);
  }

  pos = atom->position;

  weight = trilinear_interpolate(map,pos,1.0E30);

  return(weight);
}



static void pocket_set_oflags(void) {

  pocket_oflags = (params_get_parameter("output"))->value.s;

  if ((word_in_text("sites",pocket_oflags,',')) || (word_in_text("atoms",pocket_oflags,',')) || (word_in_text("residues",pocket_oflags,',')) ||
      (word_in_text("sitepoints",pocket_oflags,',')) || (word_in_text("sitemaps",pocket_oflags,','))) {
    
    pocket_calc_sites = 1;

  } else {

    pocket_calc_sites = 0;
  }
}



static void pocket_output(MAP *map,MAP_ISLAND *sites,int n_sites) {

  if (word_in_text("map",pocket_oflags,',')) {

    pocket_output_map(map);
  }
  
  if (pocket_calc_sites) {
    
    pocket_output_sites(map,sites,n_sites);
  }
}



static void pocket_output_map(MAP *map) {

  char filename[MAX_LINE_LEN];

  sprintf(filename,"%s.grd.gz",map->name);
  
  write_insight_map(filename,map,0);
}



static void pocket_output_sites(MAP *map,MAP_ISLAND *sites,int n_sites) {

  int i;
  MAP_ISLAND *site;

  for (i=0,site=sites;i<n_sites;i++,site++) {
    
    pocket_output_site(map,site);
  }
}



static void pocket_output_site(MAP *map,MAP_ISLAND *site) {

  char *oformat,text[MAX_LINE_LEN];

  oformat = (params_get_parameter("oformat"))->value.s;

  if (!strcmp(oformat,"json")) {

    set_siteio(NULL,"basic");
    set_residueio(NULL,"basic");
    set_atomio(NULL,"id");

  } else {

    set_siteio(NULL,"header,basic");
    set_residueio(NULL,"header,basic");
    set_atomio(NULL,"header,id");
  }

  site2text(site,text);

  write_line(PLI_STDOUT,"%s\n",text);
  
  if (word_in_text("residues",pocket_oflags,',')) {
    
    pocket_output_site_residues(site);
  }
  
  if (word_in_text("atoms",pocket_oflags,',')) {

    pocket_output_site_atoms(site);
  }
    
  if (word_in_text("sitepoints",pocket_oflags,',')) {

    pocket_output_site_points(map,site);
  }
    
  if (word_in_text("sitemaps",pocket_oflags,',')) {

    pocket_output_site_map(map,site);
  }

  set_siteio(NULL,NULL);
  set_residueio(NULL,NULL);
  set_atomio(NULL,NULL);
}



static void pocket_output_site_residues(MAP_ISLAND *site) {

  int i,n_residues;
  RESIDUE *residues,*residue;

  residues = atomlist2residues(site->site_atoms,&n_residues);
  
  if (residues) {
    
    for (i=0,residue=residues;i<n_residues;i++,residue++) {
      
      write_residue(PLI_STDOUT,residue);
    }
  }
  
  free(residues);
}



static void pocket_output_site_atoms(MAP_ISLAND *site) {

  int i;
  ATOMLIST *list;
  ATOM **atomp;

  list = site->site_atoms;
  
  if (list) {
    
    for (i=0,atomp=list->atom;i<list->natoms;i++,atomp++) {
      
      write_atom(PLI_STDOUT,*atomp);
    }
  }
}



static void pocket_output_site_points(MAP *map,MAP_ISLAND *site) {

  int i,ix,iy,iz,index,*point;
  double *pos;
  char filename[MAX_LINE_LEN];
  PLI_FILE *file;
  ATOM *atom;
  GRID *grid;

  sprintf(filename,"%s_%d.pdb",map->name,site->id);

  file = open_file(filename,"w");

  atom = (get_field_probe("H2O"))->atom;

  atom->id = 1;
  atom->subid = 1;
  atom->occupancy = 1.0;
  atom->bfactor = 20.0;

  strcpy(atom->name,"O");
  strcpy(atom->subname,"HOH");

  grid = map->grid;
  pos = atom->position;

  for (i=0,point=site->points;i<site->n_points;i++,point++) {
      
    index = *point;
    
    if (!map_index2grid(index,&ix,&iy,&iz,grid)) {

      pos[0] = grid->flimit[0][0] + ((double) ix)*(grid->spacing);
      pos[1] = grid->flimit[1][0] + ((double) iy)*(grid->spacing);
      pos[2] = grid->flimit[2][0] + ((double) iz)*(grid->spacing);

      write_pdb_atom(file,atom);

      atom->id++;
      atom->subid++;
    }
  }

  close_file(file);
}



static void pocket_output_site_map(MAP *map,MAP_ISLAND *site) {

  char filename[MAX_LINE_LEN],name[MAX_LINE_LEN];

  pocket_site2map(site,map);
  
  strcpy(name,map->name);

  sprintf(map->name,"%s_%d",name,site->id);

  sprintf(filename,"%s_%d.grd.gz",name,site->id);
  
  write_insight_map(filename,map,0);

  strcpy(map->name,name);
}



static void pocket_mask_protein(SYSTEM *system,MAP *map) {

  int nx,ny,nz,ix,iy,iz,index;
  unsigned char *mask;
  double ***matrix;
  GRID *grid;

  init_mask(map->mask,map->grid);

  mask_molecule(system->protein,map->mask,map->grid,VDW_MASK);

  grid = map->grid;
  mask = map->mask;
  matrix = (double***) map->matrix;

  nx = grid->npoints[0];
  ny = grid->npoints[1];
  nz = grid->npoints[2];

  for (ix=0;ix<nx;ix++) {

    for (iy=0;iy<ny;iy++) {

      for (iz=0;iz<nz;iz++) {

	index = ny*nz*ix + nz*iy + iz;

	if (is_mask_bit_set(mask,index)) {

	  matrix[ix][iy][iz] = 0.0;
	}
      }
    }
  }
}



static MAP_ISLAND* pocket_map2sites(MAP *map,int *n_sites) {

  int i,nx,ny,nz,index,ix,iy,iz,n_alloc_sites,finished,n_points,max_n_points;
  unsigned char *mask;
  double ***dm,max_level,max_volume;
  double max_integral,gridvol;
  GRID *grid;
  MAP_ISLAND *sites,*site;

  if (strcmp(map->type->name,"double")) {

    error_fn("%s: map '%s' is not of type double",__func__,map->name);
  }

  grid = map->grid;

  gridvol = pow(grid->spacing,3.0);

  max_volume = (params_get_parameter("pocket_drug_volume"))->value.d;

  max_n_points = floor(max_volume/gridvol);

  nx = grid->npoints[0];
  ny = grid->npoints[1];
  nz = grid->npoints[2];

  mask = map->mask;

  init_mask(mask,grid);
  
  dm = (double***) map->matrix;
  
  n_alloc_sites = 100;
  
  sites = (MAP_ISLAND*) calloc(n_alloc_sites,sizeof(MAP_ISLAND));
  
  if (sites == NULL) {
    
    error_fn("%s: out of memory allocating sites",__func__);
  }
  
  *n_sites = 0;
  
  finished = 0;
  
  while (((index = map_extreme_index(map,1,"max",&ix,&iy,&iz)) != -1) && (!finished)) {
    
    max_level = dm[ix][iy][iz];
    
    max_integral = ((double) max_n_points)*max_level*gridvol;
    
    if (max_integral > pocket_min_integral) {
      
      if (*n_sites == n_alloc_sites) {
	
	n_alloc_sites *= 2;
	sites = (MAP_ISLAND*) realloc(sites,n_alloc_sites*sizeof(MAP_ISLAND));
	
	if (sites == NULL) {
	  
	  error_fn("%s: out of memory re-allocating sites",__func__);
	}
      }
      
      site = sites + (*n_sites);
      
      init_map_island(site,1000);
      
      add_point_to_map_island(site,index,1);
      
      pocket_grow_site(map,site,max_n_points);
      
      calc_map_island_properties(site,map);
      
      if (site->integral > pocket_min_integral) {
	
	pocket_mask_site(map,site);
	
	site->id = *n_sites;
	
	(*n_sites)++;
      }
      
    } else {
      
      finished = 1;
    }
  }
  
  if (*n_sites > 0) {

    sites = (MAP_ISLAND*) realloc(sites,(*n_sites)*sizeof(MAP_ISLAND));
    
    for (i=0,site=sites;i<*n_sites;i++,site++) {
      
      if (site->n_points) {
	
	site->n_alloc_points = site->n_points;
	
	site->points = (int*) realloc(site->points,(site->n_points)*sizeof(int));
      }
    }
  }

  return(sites);
}



static void pocket_calc_site_atoms(MAP_ISLAND *site,SYSTEM *system,GRID *grid) {

  int i,j,ix,iy,iz,*point;
  double pos[4],dr;
  ATOM *atom1,*atom2;
  ATOMLIST *list;
  CONTACT *contact;
  CONTACTLIST *contactlist;
  FIELD_PROBE *probe;

  probe = get_field_probe("H2O");

  probe->type = ATOM_PROBE;

  init_field_system(system,probe);

  pos[3] = 1.0;

  atom1 = probe->atom;

  dr = atom1->vdw_radius_H2O - atom1->vdw_radius;

  atom1->vdw_radius = (params_get_parameter("pocket_probe_radius"))->value.d;
  atom1->vdw_radius_H2O = atom1->vdw_radius + dr;

  list = alloc_atomlist();

  init_atomlist(list);

  for (i=0,point=site->points;i<site->n_points;i++,point++) {

    if (!map_index2grid(*point,&ix,&iy,&iz,grid)) {
    
      pos[0] = grid->flimit[0][0] + ((double) ix)*(grid->spacing);
      pos[1] = grid->flimit[1][0] + ((double) iy)*(grid->spacing);
      pos[2] = grid->flimit[2][0] + ((double) iz)*(grid->spacing);

      move_probe_to_new_position(probe,pos);
	    
      set_contacts_atom(atom1,system,0);
	    
      contactlist = atom1->contactlist;
	    
      if (contactlist != NULL) {

	for (j=0,contact=contactlist->contacts;j<contactlist->ncontacts;j++,contact++) {
		
	  atom2 = contact->atom2;
	  
	  if (!atom_in_list(list,atom2)) {

	    add_atom_to_list(list,atom2,0);
	  }
	}
      }
    }    
  }
  
  site->site_atoms = list;

  finish_field_system(system,probe);
}



static void pocket_grow_site(MAP *map,MAP_ISLAND *site,int max_points) {

  int i,n_points,n_edge_points,n_total_points,ix,iy,iz,ixc,iyc,izc,*point,index,maxi;
  unsigned char* mask;
  double ***matrix,f,maxf;
  GRID *grid;

  grid = map->grid;
  matrix = (double***) map->matrix;
  mask = (unsigned char*) map->mask;

  n_points = site->n_points;
  n_edge_points = site->n_edge_points;

  if ((n_points >= max_points) || (n_edge_points == 0)) {

    return;
  }

  // find edge point with maximum value:

  n_total_points = n_points + n_edge_points;

  maxf = -1.0E10;
  maxi = -1;

  for (i=n_points,point=site->points+n_points;i<n_total_points;i++,point++) {

    index = *point;

    if (!map_index2grid(index,&ix,&iy,&iz,grid)) {
    
      f = matrix[ix][iy][iz];

      if (f > maxf) {

	maxf = f;

	maxi = i;
      }
    }    
  }

  if (maxi == -1) {

    error_fn("%s: unexpected error occurred (1)",__func__);
  }

  index = site->points[maxi];

  if (map_index2grid(index,&ixc,&iyc,&izc,grid)) {

    error_fn("%s: unexpected error occurred (2)",__func__);
  }

  // add best edge point to the site:

  site->points[maxi] = site->points[n_points];
  site->points[n_points] = index;

  set_mask_bit(mask,index);

  site->n_points++;
  site->n_edge_points--;

  // add edge points around added point:

  for (ix=ixc-1;ix<=ixc+1;ix++) {

    for (iy=iyc-1;iy<=iyc+1;iy++) {

      for (iz=izc-1;iz<=izc+1;iz++) {

	index = map_grid2index(ix,iy,iz,grid);

	if (index != -1) {

	  if (!is_mask_bit_set(mask,index)) {

	    // TODO: don't need to check all points:

	    if (!pocket_point_in_list(index,site->points,n_total_points)) {

	      add_point_to_map_island(site,index,1);
	    }
	  }
	}
      }
    }
  }

  // recurse:

  pocket_grow_site(map,site,max_points);
}



static void pocket_site2map(MAP_ISLAND *site,MAP *map) {

  int i,index,*point,nx,ny,nz,ix,iy,iz;
  double ***matrix;
  GRID *grid;

  grid = map->grid;

  nx = grid->npoints[0];
  ny = grid->npoints[1];
  nz = grid->npoints[2];
   
  matrix = (double***) map->matrix;

  for (ix=0;ix<nx;ix++) {
      
    for (iy=0;iy<ny;iy++) {
	
      for (iz=0;iz<nz;iz++) {
	
	matrix[ix][iy][iz] = 0.0;
      }
    }
  }

  for (i=0,point=site->points;i<site->n_points;i++,point++) {
      
    index = *point;
    
    if (!map_index2grid(index,&ix,&iy,&iz,grid)) {
      
      matrix[ix][iy][iz] = 1.0;
    }
  }
}



static void pocket_mask_site(MAP *map,MAP_ISLAND *site) {

  int n_points,i,index,*point,n,ixc,iyc,izc,ix,iy,iz;
  unsigned char *mask;
  GRID *grid;

  n_points = ((params_get_parameter("pocket_allow_overlap"))->value.i) ? site->n_frag_points : site->n_points;

  grid = map->grid;
  mask = map->mask;

  n = 3;

  for (i=0,point=site->points;i<n_points;i++,point++) {
      
    if (!map_index2grid(*point,&ixc,&iyc,&izc,grid)) {
      
      for (ix=ixc-n;ix<=ixc+n;ix++) {

	for (iy=iyc-n;iy<=iyc+n;iy++) {

	  for (iz=izc-n;iz<=izc+n;iz++) {

	    index = map_grid2index(ix,iy,iz,grid);

	    if (index) {

	      set_mask_bit(mask,index);
	    }
	  }
	}
      }
    }
  }
}



static int pocket_point_in_list(int index,int *points,int n_points) {

  int i,*point;
  
  for (i=0,point=points;i<n_points;i++,point++) {

    if (*point == index) {

      return(1);
    }
  }

  return(0);
}
