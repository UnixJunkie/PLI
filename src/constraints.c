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



double score_system_constraints(SYSTEM *system) {

  int i;
  double score,*tether_pos;
  ATOM **atomp,*atom;
  ATOMLIST *selection;

  selection = system->selection;

  if (selection == NULL) {

    return(0.0);
  }

  score = 0.0;

  // atom position tethers:

  for (i=0,atomp=selection->atom;i<selection->natoms;i++,atomp++) {
  
    atom = *atomp;

    tether_pos = atom->tether_position;

    if (tether_pos) {

      atom->constraint_score = (atom->tether_k)*sqr_distance(atom->position,tether_pos);

      score += atom->constraint_score;

    } else {

      atom->constraint_score = 0.0;
    }
  }

  system->constraint_score = score;

  return(score);
}



double constraint_score_gradient(DOF *variable,SYSTEM *system) {

  double fd;

  fd = numerical_score_gradient(variable,system,score_system_constraints);

  return(fd);
}



void prep_min_atom_tethers(SYSTEM *system) {

  int i;
  double tether_k;
  ATOM *atom;
  MOLECULE *ligand;

  tether_k = (params_get_parameter("min_tether_k"))->value.d;

  // TODO: probably need to change this for unrestrained 'docking'

  ligand = system->ligand;
  
  if (ligand) {
    
    for (i=0,atom=ligand->atom;i<ligand->natoms;i++,atom++) {
      
      atom->tether_position = atom->original_position;
      
      atom->tether_k = tether_k;
    }
  }
}



void unprep_min_atom_tethers(SYSTEM *system) {

  int i;
  ATOM *atom;
  MOLECULE *ligand;

  // TODO: probably need to change this for unrestrained 'docking'

  ligand = system->ligand;
  
  if (ligand) {
    
    for (i=0,atom=ligand->atom;i<ligand->natoms;i++,atom++) {
      
      atom->tether_position = NULL;
    }
  }
}



void prep_template_atom_tethers(SYSTEM *system) {

  int i;
  double tether_k;
  ATOM *atom1,*atom2;
  MATCH *match;

  match = system->template_match;

  if (match) {

    tether_k = (params_get_parameter("min_tether_k"))->value.d;

    for (i=0;i<match->n_atoms;i++) {

      atom1 = match->atoms1[i];
      atom2 = match->atoms2[i];

      atom2->tether_position = atom1->position;

      atom2->tether_k = tether_k;
    }
  }
}



void unprep_template_atom_tethers(SYSTEM *system) {

  int i;
  ATOM *atom;
  MATCH *match;

  match = system->template_match;

  if (match) {

    for (i=0;i<match->n_atoms;i++) {

      atom = match->atoms2[i];

      atom->tether_position = NULL;
    }
  }
}
