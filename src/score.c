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


static void (*score_system_fp)(struct System*);
static void (*minimise_system_fp)(struct System*);
static void (*write_system_scores_fp)(PLI_FILE*,struct System*,enum OUTPUT_FORMAT,unsigned long int);
static void (*write_atom_scores_fp)(PLI_FILE*,struct Atom*,enum OUTPUT_FORMAT,unsigned long int);
static void (*write_contact_scores_fp)(PLI_FILE*,struct Contact*,enum OUTPUT_FORMAT,unsigned long int);

static void init_system_ref_scores(SYSTEM*);



void run_score(SETTINGS *settings) {

  int i;
  SYSTEM *system;
  SYSTEM_LIST *syslist;

  syslist = settings2syslist(settings);

  // loop over systems:

  for (i=0,system=syslist->systems;i<syslist->n_systems;i++,system++) {

    prep_score(system,settings->sfunc);
    
    score_system(system,settings->sfunc);

    unprep_score(system);
  }
}



void prep_score(SYSTEM *system,PLI_SFUNC *sfunc) {

  prep_system(system,ANY_MOLECULE);

  calc_system_ref_scores(system,sfunc);
}



void unprep_score(SYSTEM *system) {

  output_system(PLI_STDOUT,system);

  unprep_system(system,LIGAND_MOLECULE);
}



void score_system(SYSTEM *system,PLI_SFUNC *sfunc) {

  sfunc->score_system(system);

  score_system_constraints(system);
}



void minimise_system(SYSTEM *system,PLI_SFUNC *sfunc) {

  if (sfunc->minimise_system) {

    sfunc->minimise_system(system);

  } else {

    minimise(system,sfunc,-1);
  }
}



void calc_system_ref_scores(SYSTEM *system,PLI_SFUNC *sfunc) {

  int i;
  ATOM **atomp;
  MOLECULE **sysmolp;
  SYSTEM *molsystem;

  init_system_scores(system);

  init_system_ref_scores(system);

  reset_system(system,0);

  for (i=0,sysmolp=system->molecule_list->molecules;i<system->molecule_list->n_molecules;i++,sysmolp++) {

    molsystem = (*sysmolp)->molsystem;

    init_system_ref_scores(molsystem);

    score_system(molsystem,sfunc);

    add_system_scores(system,molsystem);
  }

  store_system_ref_scores(system);

  for (i=0,atomp=system->selection->atom;i<system->selection->natoms;i++,atomp++) {

    store_atom_ref_scores(*atomp);
  }

  reset_system(system,0);
}



void write_system_scores(PLI_FILE *file,SYSTEM *system,enum OUTPUT_FORMAT oformat,unsigned long int oflags) {

  if (write_system_scores_fp) {

    write_system_scores_fp(file,system,oformat,oflags);

  } else {

    write_line(file,(oformat == JSON) ? "\"score\":%lf" : " %10.4lf",system->score);
  }
}



void write_atom_scores(PLI_FILE *file,ATOM *atom,enum OUTPUT_FORMAT oformat,unsigned long int oflags) {

  if (write_atom_scores_fp) {

    write_atom_scores_fp(file,atom,oformat,oflags);

  } else {
    
    write_line(file,(oformat == JSON) ? "\"score\":%lf" : " %10.4lf",atom->score);
  }
}



void write_contact_scores(PLI_FILE *file,CONTACT *contact,enum OUTPUT_FORMAT oformat,unsigned long int oflags) {

  if (write_contact_scores_fp) {

    write_contact_scores_fp(file,contact,oformat,oflags);

  } else {
    
    write_line(file,(oformat == JSON) ? "\"score\":%lf" : " %10.4lf",contact->score);
  }
}



void init_score(SETTINGS *settings) {

  score_system_fp = pliff_score_system;
  minimise_system_fp = pliff_minimise_system;
  write_system_scores_fp = pliff_write_system_scores;
  write_atom_scores_fp = pliff_write_atom_scores;
  write_contact_scores_fp = pliff_write_contact_scores;
}



void init_system_scores(SYSTEM *system) {

  int i;

  system->score = 0.0;
  system->constraint_score = 0.0;

  for (i=0;i<MAX_SYSTEM_SCORES;i++) {

    system->scores[i] = 0.0;
  }
}



static void init_system_ref_scores(SYSTEM *system) {

  int i;

  system->ref_score = 0.0;

  for (i=0;i<MAX_SYSTEM_SCORES;i++) {

    system->ref_scores[i] = 0.0;
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
