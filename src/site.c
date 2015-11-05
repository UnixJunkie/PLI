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



static struct SiteSettings {
  enum { SITE_FROM_FILE,SITE_FROM_LIGAND } def_type;
  double max_dist;
  char ligand_file[MAX_LINE_LEN];
  char site_file[MAX_LINE_LEN];
  MOLECULE *ligand;
  int n_atom_ids;
  int *atom_ids;
} *site_settings = NULL;



static void read_site_file(char*);



void init_site_settings(void) {

  char *def_type;

  if (site_settings == NULL) {

    site_settings = (struct SiteSettings*) malloc(sizeof(struct SiteSettings));
  }

  if (site_settings == NULL) {

    error_fn("init_site_settings: out of memory allocating settings");
  }

  // site definition:

  def_type = (params_get_parameter("site_def_type"))->value.s;

  if (!strcmp(def_type,"ligand")) {

    site_settings->def_type = SITE_FROM_LIGAND;

  } else if (!strcmp(def_type,"file")) {

    site_settings->def_type = SITE_FROM_FILE;

  } else {

    error_fn("init_site_settings: no such site_def_type '%s'",def_type);
  }

  // maximum distance for atoms to be part of site:

  site_settings->max_dist = (params_get_parameter("site_max_dist"))->value.d;

  // ligand file:

  strcpy(site_settings->ligand_file,(params_get_parameter("site_ligand_file"))->value.s);

  // site file (containst site atoms):

  strcpy(site_settings->site_file,(params_get_parameter("site_file"))->value.s);

  if (strcmp(site_settings->site_file,"undefined")) {

    site_settings->def_type = SITE_FROM_FILE;
  }

  site_settings->ligand = NULL;

  site_settings->n_atom_ids = 0;
  site_settings->atom_ids = NULL;
}



void init_site(SETTINGS *settings) {

  if (!strcmp(settings->protein_file,"undefined")) {

    return;
  }

  if ((strcmp(settings->selection,"site")) && (strcmp(settings->selection,"complex")) && (strcmp(settings->resolve_side_chains,"site"))) {

    return;
  }

  if (site_settings->def_type == SITE_FROM_LIGAND) {

    if (!strcmp(site_settings->ligand_file,"undefined")) {

      if (!strcmp(settings->ligand_file,"undefined")) {
	error_fn("init_site: site ligand file undefined");
      }
    
      warning_fn("init_site: using ligand file to define site");

      strcpy(site_settings->ligand_file,settings->ligand_file);
    }

    site_settings->ligand = read_molecule(site_settings->ligand_file);

    if (site_settings->ligand == NULL) {

      error_fn("init_site: site ligand undefined");
    }

  } else if (site_settings->def_type == SITE_FROM_FILE) {

    if (!strcmp(site_settings->site_file,"undefined")) {

      error_fn("init_site: site file undefined");      
    }

    read_site_file(site_settings->site_file);

  } else {

    error_fn("init_site: option not implemented");
  }
}



ATOMLIST* molecule2site(MOLECULE *molecule) {

  ATOMLIST *site_atoms;
  MOLECULE *ligand;

  site_atoms = NULL;

  if (site_settings->def_type == SITE_FROM_LIGAND) {

    ligand = NULL;

    if (site_settings->ligand) {

      ligand = site_settings->ligand;

    } else if ((molecule->system) && (molecule->system->ligand)) {
      
      ligand = molecule->system->ligand;
    }

    if (ligand == NULL) {

      error_fn("molecule2site: site ligand undefined");
    }

    site_atoms = molecule_sphere_to_atom_list(molecule,ligand,(params_get_parameter("site_max_dist"))->value.d);

  } else if (site_settings->def_type == SITE_FROM_FILE) {

    if (molecule->flags & PROTEIN_MOLECULE) {

      site_atoms = molecule_atom_ids_to_atom_list(molecule,site_settings->atom_ids,site_settings->n_atom_ids);
    }
  }

  if (site_atoms == NULL) {

    error_fn("molecule2site: no site atoms found");
  }

  return(site_atoms);
}



static void read_site_file(char *filename) {

  int atom_id,n_atom_ids,n_alloc_atom_ids;
  char line[MAX_LINE_LEN];
  PLI_FILE *site_file;

  site_file = open_file(filename,"r");

  if (site_file == NULL) {

    error_fn("read_site_file: error opening site file \'%s\'",filename);
  }

  if (site_settings->atom_ids != NULL) {

    free(site_settings->atom_ids);
  }

  n_atom_ids = n_alloc_atom_ids = 0;

  while (!end_of_file(site_file)) {

    if (read_line(line,MAX_LINE_LEN,site_file) == NULL)
      break;

    if (sscanf(line,"%d",&atom_id) == 1) {

      if (n_alloc_atom_ids == 0) {

	n_alloc_atom_ids = 1000;

	site_settings->atom_ids = (int*) calloc(n_alloc_atom_ids,sizeof(int));

      } else if (site_settings->n_atom_ids == n_alloc_atom_ids) {

	n_alloc_atom_ids += 1000;

	site_settings->atom_ids = (int*) realloc(site_settings->atom_ids,n_alloc_atom_ids*sizeof(int));
      }

      if (site_settings->atom_ids == NULL) {

	error_fn("read_site_file: out of memory (re)allocating site_atoms");
      }

      site_settings->atom_ids[n_atom_ids++] = atom_id;
    }
  }
  
  close_file(site_file);

  site_settings->n_atom_ids = n_atom_ids;
}

