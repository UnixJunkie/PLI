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



void score_system_constraints(SYSTEM *system) {

  int i;
  double *tether_pos;
  ATOM **atomp,*atom;
  ATOMLIST *selection;

  selection = system->selection;

  if (selection == NULL) {

    return;
  }

  system->constraint_score = 0.0;

  // atom position tethers:

  for (i=0,atomp=selection->atom;i<selection->natoms;i++,atomp++) {
  
    atom = *atomp;

    tether_pos = atom->tether_position;

    if (tether_pos) {

      atom->constraint_score = (atom->tether_k)*sqr_distance(atom->position,tether_pos);

      system->constraint_score += atom->constraint_score;

    } else {

      atom->constraint_score = 0.0;
    }
  }

  system->score += system->constraint_score;
}



void setup_list_atom_tethers(ATOMLIST *list,double tether_k) {

  int i;
  ATOM **atomp,*atom;

  for (i=0,atomp=list->atom;i<list->natoms;i++,atomp++) {
  
    atom = *atomp;

    atom->tether_position = atom->original_position;

    atom->tether_k = tether_k;
  }
}
