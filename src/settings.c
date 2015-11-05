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



static void init_settings(SETTINGS*);
static unsigned int get_output_flags(char*);



SETTINGS* get_settings(int nargs,char **args) {

  SETTINGS *settings;

  // set all parameters from command line and input file(s):

  params_init(nargs,args);

  // io setup:

  init_io();

  // allocate memory for settings:

  settings = (SETTINGS*) malloc(sizeof(SETTINGS));

  if (settings == NULL) {

    error_fn("get_settings: out of memory allocating settings");
  }

  // initialise "global" parameters and settings:

  init_settings(settings);

  // initialise parameters in other files:

  init_histogram_settings();
  init_site_settings();
  init_field_settings();
  init_contact_settings();
  init_ff_settings();
  init_minimise_settings();
  init_pliff_settings();

  // other settings:

  init_site(settings);
  init_ff(settings);
  init_score(settings);
  init_field_probes(settings);

  // load system molecules:

  settings->sysmols = load_sysmols(settings);

  return(settings);
}


// init_settings()
//
// the settings structure predates the pli_params (params.c)
// now, many of these could be replaced with simple params_get_parameter calls
// wherever they are needed; will get done at some point

static void init_settings(SETTINGS *settings) {

  int i;
  double max_cov_radius,max_vdw_radius;
  ELEMENT *element;
  UNITED_ATOM *uatom;
  ATOM_TYPING_SCHEME *scheme;

  settings->pli_dir = get_pli_dir();

  sprintf(settings->settings_file,"%s/params/settings.pli",settings->pli_dir);
  sprintf(settings->types_file,"%s/params/types.pli",settings->pli_dir);

  // general parameters:

  settings->mode = get_mode((params_get_parameter("mode"))->value.s);

  settings->sfunc = get_sfunc((params_get_parameter("sfunc"))->value.s);

  strcpy(settings->jobname,(params_get_parameter("jobname"))->value.s);

  strcpy(settings->protein_file,(params_get_parameter("protein"))->value.s);
  strcpy(settings->ligand_file,(params_get_parameter("ligand"))->value.s);
  strcpy(settings->symmetry_file,(params_get_parameter("symmetry"))->value.s);
  strcpy(settings->selection,(params_get_parameter("selection"))->value.s);

  // output control:

  settings->verbose = (params_get_parameter("verbose"))->value.i;

  settings->oflags = get_output_flags((params_get_parameter("output"))->value.s);
  settings->oformat = get_output_format((params_get_parameter("oformat"))->value.s);
  settings->astyle = get_atom_style((params_get_parameter("astyle"))->value.s);

  strcpy(settings->oatoms,(params_get_parameter("oatoms"))->value.s);

  settings->report_alt_types = (params_get_parameter("report_alt_types"))->value.i;

  // atom typing:

  settings->atom_typing_scheme = alloc_atom_typing_scheme();

  settings->protein_use_pdb_dict = (params_get_parameter("pdict"))->value.i;
  settings->ligand_use_pdb_dict = (params_get_parameter("ldict"))->value.i;

  settings->protein_use_hydrogens = (params_get_parameter("phs"))->value.i;
  settings->ligand_use_hydrogens = (params_get_parameter("lhs"))->value.i;

  strcpy(settings->resolve_side_chains,(params_get_parameter("resolve_side_chains"))->value.s);
  strcpy(settings->alt_type_assignment,(params_get_parameter("alt_type_assignment"))->value.s);

  settings->resolve_metals = (params_get_parameter("resolve_metals"))->value.i;

  // test parameters:

  settings->id1 = (params_get_parameter("id1"))->value.i;
  settings->id2 = (params_get_parameter("id2"))->value.i;

  strcpy(settings->hist_name,(params_get_parameter("hist_name"))->value.s);

  // other parameters:

  settings->minimise = (params_get_parameter("minimise"))->value.i;
  strcpy(settings->keep_waters,(params_get_parameter("keep_waters"))->value.s);

  // a few fixed force field parameters:

  settings->force_field = NULL;

  settings->water_vdw_radius = 1.4;
  settings->max_covalent_dist = 2.5;
  settings->max_contact_dist = 6.5;

  // atom typing scheme:

  read_atom_typing_scheme(settings->types_file,settings->atom_typing_scheme);
  scheme = settings->atom_typing_scheme;
  // set water radius:

  uatom = get_united_atom("H2O",scheme);

  if (uatom == NULL) {

    error_fn("get_settings: water united atom undefined");
  }

  settings->water_vdw_radius = uatom->vdw_radius;

  // set maximum covalent distance:

  max_cov_radius = -1.0;

  for (i=0,element=scheme->elements;i<scheme->n_elements;i++,element++) {

    if (element->cov_radius > max_cov_radius) {

      max_cov_radius = element->cov_radius;
    }
  }

  settings->max_covalent_dist = 2.0*max_cov_radius + COVALENT_TOLERANCE;

  // set maximum contact distance:

  max_vdw_radius = -1.0;

  for (i=0,uatom=scheme->united_atoms;i<scheme->n_united_atoms;i++,uatom++) {

    if (uatom->vdw_radius > max_vdw_radius) {

      max_vdw_radius = uatom->vdw_radius;
    }
  }

  settings->max_contact_dist = 2.0*(max_vdw_radius + settings->water_vdw_radius);

  //settings->max_contact_dist += 3.0;

  // set selection:

  if (!strcmp(settings->selection,"undefined")) {

    strcpy(settings->selection,settings->mode->selection);
  }

  if ((settings->mode->selection_deviation_warning) && (strcmp(settings->selection,settings->mode->selection))) {

    warning_fn("get_settings: running with selection '%s' instead of default selection '%s'",settings->selection,settings->mode->selection);
  }
}



static unsigned int get_output_flags(char *flags) {

  unsigned int oflags;
  char *flag;

  oflags = OUTPUT_NONE;

  flag = strtok(flags,",");

  while (flag) {

    if (!strcmp(flag,"system")) {

      oflags |= OUTPUT_SYSTEM;
      
    } else if (!strcmp(flag,"atoms")) {
      
      oflags |= OUTPUT_ATOMS;
      
    } else if (!strcmp(flag,"contacts")) {
      
      oflags |= OUTPUT_CONTACTS;
      
    } else if (!strcmp(flag,"hbonds")) {
      
      oflags |= OUTPUT_HBONDS;
      
    } else if (!strcmp(flag,"clashes")) {
      
      oflags |= OUTPUT_CLASHES;
      
    } else if (!strcmp(flag,"scores")) {
      
      oflags |= OUTPUT_SCORES;
      
    } else if (!strcmp(flag,"pdb")) {
            
      oflags |= OUTPUT_PDB;
      
    } else if (!strcmp(flag,"mol")) {
            
      oflags |= OUTPUT_MOL;
      
    } else if (!strcmp(flag,"axes")) {
      
      oflags |= OUTPUT_AXES;
      
    } else if (!strcmp(flag,"vpts")) {
      
      oflags |= OUTPUT_VPTS;
   
    } else if (!strcmp(flag,"geometries")) {
      
      oflags |= OUTPUT_GEOMETRIES;
      
    } else if (!strcmp(flag,"covalent")) {

      oflags |= OUTPUT_COVALENT;

    } else if (!strcmp(flag,"intra")) {

      oflags |= OUTPUT_INTRAMOLECULAR;
      
    } else if (!strcmp(flag,"sas_stats")) {

      oflags |= OUTPUT_SAS_STATS;
      
    } else if (!strcmp(flag,"ligand")) {

      oflags |= OUTPUT_LIGAND;

    } else if (!strcmp(flag,"simb")) {

      oflags |= OUTPUT_SIMB;

    } else if (!strcmp(flag,"stypes")) {

      oflags |= OUTPUT_STYPES;  
 
    } else if (!strcmp(flag,"site_stats")) {

      oflags |= OUTPUT_SITE_STATS;  
 
    } else if (strcmp(flag,"none")) {

      error_fn("set_output_flags: unknown output option '%s'",flag);
    }

    flag = strtok(NULL,",");
  }

  return(oflags);
}


