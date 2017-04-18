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



static double (*last_score_fp)(struct System*) = NULL;






static void prep_score_molecule(MOLECULE*);
static void unprep_score_molecule(MOLECULE*);
static void check_score_function(SYSTEM*,double (*score_fp)(struct System*));
static void init_system_ref_scores(SYSTEM*);



void run_score(SETTINGS *settings) {

  int i;
  SYSTEM *system;
  SYSTEM_LIST *syslist;

  syslist = settings2syslist(settings);

  // loop over systems:

  for (i=0,system=syslist->systems;i<syslist->n_systems;i++,system++) {

    score_system(system,settings->sfunc);

    output_system(PLI_STDOUT,system);

    unprep_system(system,LIGAND_MOLECULE);
  }

  // free up stuff:

  for (i=0,system=syslist->systems;i<syslist->n_systems;i++,system++) {

    unprep_system(system,ANY_MOLECULE);

    free_system(system);
  }

  free(syslist);
}



void prep_score_system(SYSTEM *system,PLI_SFUNC *sfunc) {

  int i;
  MOLECULE **moleculep,*molecule;
  MOLECULE_LIST *list;

  if (system->score_prepped) {

    return;
  }

  prep_system(system,ANY_MOLECULE);

  list = system->molecule_list;

  for (i=0,moleculep=list->molecules;i<list->n_molecules;i++,moleculep++) {

    molecule = *moleculep;
    
    prep_score_molecule(molecule);
  }

  set_atom_depth_maps(system);

  calc_system_ref_scores(system,sfunc);

  system->score_prepped = TRUE;
}



void unprep_score_system(SYSTEM *system,FLAG flags,PLI_SFUNC *sfunc) {

  int i;
  MOLECULE **moleculep,*molecule;
  MOLECULE_LIST *list;

  list = system->molecule_list;

  for (i=0,moleculep=list->molecules;i<list->n_molecules;i++,moleculep++) {

    molecule = *moleculep;

    if (molecule->flags & flags) {

      unprep_score_molecule(molecule);
    }
  }

  system->score_prepped = FALSE;
}



static void prep_score_molecule(MOLECULE *molecule) {

  char *depth_weighting;
  POCKET_MODE *mode;

  if (molecule->flags & PROTEIN_MOLECULE) {

    depth_weighting = (params_get_parameter("depth_weighting"))->value.s;

    if (strcmp(depth_weighting,"none")) {

      if (molecule->depth_map == NULL) {

        mode = pocket_mode(depth_weighting);

        molecule->depth_map = mode->map(molecule->molsystem);
      }
    }
  }

  if (molecule->use_grids) {

    if (molecule->conmap == NULL) {
      
      //molecule->conmap = molecule2contactmap(molecule);
    }
  }
}



static void unprep_score_molecule(MOLECULE *molecule) {

  if (molecule->flags & PROTEIN_MOLECULE) {

    if (molecule->depth_map != NULL) {

      free_map(molecule->depth_map);

      molecule->depth_map = NULL;
    }
  }

  if (molecule->conmap) {

    free_map(molecule->conmap);

    molecule->conmap = NULL;
  }
}



void score_system(SYSTEM *system,PLI_SFUNC *sfunc) {

  // TODO: currently the setting of contacts (set_contacts_system) is 
  // linked to the pliff scoring function specifically
  //
  // it may be a good idea to have a generic function that deals with this
  // so other scoring functions could also (choose to) use it
  //
  // this could then also deal with the setting and unsetting of atom->status flags,
  // which is currently also exclusive to the pliff scoring function

  prep_score_system(system,sfunc);

  check_score_function(system,sfunc->score_system);

  sfunc->score_system(system);

  score_system_constraints(system);

  system->score += system->constraint_score;
}



void minimise_system(SYSTEM *system,PLI_SFUNC *sfunc) {

  if (sfunc->minimise_system) {

    sfunc->minimise_system(system);

  } else {

    minimise(system,sfunc,-1);
  }
}



void set_dof_score_gradient(DOF *variable,SYSTEM *system,PLI_SFUNC *sfunc) {

  if (sfunc->score_gradient) {

    variable->fd = sfunc->score_gradient(variable,system);

  } else {

    variable->fd = numerical_score_gradient(variable,system,sfunc->score_system);
  }

  variable->fd += constraint_score_gradient(variable,system);
}



// numerical score gradient
// this could be faster if the central point's score was re-used

double numerical_score_gradient(DOF *variable,SYSTEM *system,double (*score_fp)(struct System*)) {

  double score,score1,score2,fd1,fd2;

  check_score_function(system,score_fp);

  store_positions(variable->atomlist,variable->pos,variable->u,variable->v,variable->w);

  apply_dof_shift(variable,-variable->fd_shift);

  score1 = score_fp(system);

  restore_positions(variable->atomlist,variable->pos,variable->u,variable->v,variable->w);

  apply_dof_shift(variable,variable->fd_shift);

  score2 = score_fp(system);

  restore_positions(variable->atomlist,variable->pos,variable->u,variable->v,variable->w);

  score = score_fp(system);

  fd1 = (score - score1)/0.00001;
  fd2 = (score2 - score)/0.00001;

  if ((fd1 < 0.0) && (fd2 > 0.0)) {

    // minimum:

    return(0.0);

  } else if ((fd1 > 0.0) && (fd2 < 0.0)) {

    // maximum:

    return(0.0);

  } else if ((fd1 > 0.0) && (fd2 > 0.0)) {

    // uphill:

    return(fd1);

  } else {

    // downhill:

    return(fd2);
  }
}



static void check_score_function(SYSTEM *system,double (*score_fp)(struct System*)) {

  if (score_fp != last_score_fp) {

    reset_system(system,0);
  }

  last_score_fp = score_fp;
}



void calc_system_ref_scores(SYSTEM *system,PLI_SFUNC *sfunc) {

  int i;
  ATOM **atomp;
  MOLECULE **sysmolp;
  SYSTEM *molsystem;

  init_system_scores(system);

  init_system_ref_scores(system);

  reset_system(system,0);

  int use_ref_score = params_get_parameter("use_ref_score")->value.i;
  for (i=0,sysmolp=system->molecule_list->molecules;i<system->molecule_list->n_molecules;i++,sysmolp++) {

    molsystem = (*sysmolp)->molsystem;

    init_system_ref_scores(molsystem);

    // TODO: refrence scores are turned off by default because docking or mapping poses will each get
    // a different reference score when rescored, and then their (re)scores will not match to the
    // scores from the mapping.
    if (use_ref_score) {
      score_system(molsystem,sfunc);

      add_system_scores(system,molsystem);
    }
  }

  store_system_ref_scores(system);

  for (i=0,atomp=system->selection->atom;i<system->selection->natoms;i++,atomp++) {

    store_atom_ref_scores(*atomp);
  }

  reset_system(system,0);
}



void write_system_scores(PLI_FILE *file,SYSTEM *system) {

  SETTINGS *settings;
  PLI_SFUNC *sfunc;

  settings = get_settings();

  sfunc = settings->sfunc;

  if (sfunc->write_system_scores) {

    sfunc->write_system_scores(file,system);

  } else {

    if (settings->oformat == JSON) {

      write_line(file,"\"score\":%lf",system->score);

    } else {

      write_line(file," %10.4lf",system->score);
    }
  }
}



void write_atom_scores(PLI_FILE *file,ATOM *atom) {

  SETTINGS *settings;
  PLI_SFUNC *sfunc;

  settings = get_settings();

  sfunc = settings->sfunc;

  if (sfunc->write_atom_scores) {

    sfunc->write_atom_scores(file,atom);

  } else {

    if (settings->oformat == JSON) {

      write_line(file,"\"score\":%lf",atom->score);

    } else {

      write_line(file," %10.4lf",atom->score);
    }
  }
}



void write_contact_scores(PLI_FILE *file,CONTACT *contact) {

  SETTINGS *settings;
  PLI_SFUNC *sfunc;

  settings = get_settings();

  sfunc = settings->sfunc;

  if (sfunc->write_contact_scores) {

    sfunc->write_contact_scores(file,contact);

  } else {

    if (settings->oformat == JSON) {

      write_line(file,"\"score\":%lf",contact->score);

    } else {

      write_line(file," %10.4lf",contact->score);
    }
  }
}



void init_system_scores(SYSTEM *system) {

  int i;

  system->score = 0.0;
  system->constraint_score = 0.0;

  for (i=0;i<MAX_SYSTEM_SCORES;i++) {

    system->scores[i] = 0.0;
  }
}



void init_atom_scores(ATOM *atom) {

  int i;

  atom->score = 0.0;
  atom->constraint_score = 0.0;

  for (i=0;i<MAX_ATOM_SCORES;i++) {

    atom->scores[i] = 0.0;
  }
}



void add_system_scores(SYSTEM *system1,SYSTEM *system2) {

  int i;

  system1->score += system2->score;

  for (i=0;i<MAX_SYSTEM_SCORES;i++) {

    system1->scores[i] += system2->scores[i];
  }
}



void store_system_ref_scores(SYSTEM *system) {

  int i;

  system->ref_score = system->score;

  for (i=0;i<MAX_SYSTEM_SCORES;i++) {

    system->ref_scores[i] = system->scores[i];
  }
}



void store_atom_ref_scores(ATOM *atom) {

  int i;

  atom->ref_score = atom->score;

  for (i=0;i<MAX_ATOM_SCORES;i++) {

    if (atom->flags & SKIP_ATOM) {

      atom->ref_scores[i] = 0.0;

    } else {

      atom->ref_scores[i] = atom->scores[i];
    }
  }
}



void mirror_contact_scores(CONTACT *contact1,CONTACT *contact2) {

  int i;
  double *scores1,*scores2;

  scores1 = contact1->scores;
  scores2 = contact2->scores;

  for (i=0;i<MAX_CONTACT_SCORES;i++) {

    scores2[i] = scores1[i];
  }

  contact2->score = contact1->score;
  contact2->contact_score = contact1->contact_score;
  contact2->geometry_score = contact1->geometry_score;
}



static void init_system_ref_scores(SYSTEM *system) {

  int i;

  system->ref_score = 0.0;

  for (i=0;i<MAX_SYSTEM_SCORES;i++) {

    system->ref_scores[i] = 0.0;
  }
}
