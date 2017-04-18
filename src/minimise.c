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



#define GOLDEN_SECTION_TOLERANCE 1.0E-4



static struct MinimiseSettings {
  unsigned int space;
  enum MINIMISATION_ALGORITHM algorithm;
  int max_iter;
  int tether_atoms;
  double tether_k;
} *minimise_settings = NULL;

typedef struct MinBracket {
  double ax;
  double bx;
  double cx;
  int error;
} MIN_BRACKET;



static int take_minimise_step(DOF_LIST*,SYSTEM*,enum MINIMISATION_ALGORITHM,PLI_SFUNC*);
static int take_minimise_step_sd(DOF_LIST*,SYSTEM*,PLI_SFUNC*);
static int take_minimise_step_cg(DOF_LIST*,SYSTEM*,PLI_SFUNC*);
static int check_gradients(DOF_LIST*);
static int minimise_along_gradient(DOF_LIST*,SYSTEM*,PLI_SFUNC*);
static MIN_BRACKET* bracket_minimum(DOF_LIST*,SYSTEM*,PLI_SFUNC*);
static double golden_section_search(DOF_LIST*,SYSTEM*,MIN_BRACKET*,PLI_SFUNC*);
static void init_bracket(MIN_BRACKET*,DOF_LIST*);
static double score_line_shift(DOF_LIST*,SYSTEM*,double,PLI_SFUNC*);
static void set_dof_list_gradients(DOF_LIST*,SYSTEM*,PLI_SFUNC*);
static int set_dof_list_conjugate_gradients(DOF_LIST*,SYSTEM*);



void minimise(SYSTEM *system,PLI_SFUNC *sfunc,int max_iter) {

  int i,j,flag,maxi;
  double score;
  DOF *variable;
  DOF_LIST *list;
  ATOMLIST *selection;
  PLI_FILE *file;
  ATOM **atomp,*atom;
  ATOM_COORDS *coords,*coord;

  score_system(system,sfunc);

  score = system->score;

  maxi = (max_iter > 0) ? max_iter : minimise_settings->max_iter;

  selection = system->selection;

  coords = get_list_coords(selection);

  list = setup_dof_list(system);

  if ((list) && (list->n_variables)) {

   flag = 0;
    i = 0;

    do {

      flag = take_minimise_step(list,system,minimise_settings->algorithm,sfunc);

      i++;

    } while ((!flag) && (i <= maxi));
  }

  free_dof_list(list);

  system->min_ds = system->score - score;

  system->min_dx = 0.0;

  for (i=0,atomp=selection->atom,coord=coords;i<selection->natoms;i++,atomp++,coord++) {

    atom = *atomp;

    system->min_dx += sqr_distance(atom->position,coord->position);
  }

  system->min_dx = sqrt(system->min_dx/selection->natoms);

  free(coords);
}



void init_minimise_settings(void) {

  char *algorithm;

  if (minimise_settings == NULL) {

    minimise_settings = (struct MinimiseSettings*) malloc(sizeof(struct MinimiseSettings));
  }

  if (minimise_settings == NULL) {

    error_fn("init_minimise_settings: out of memory allocating settings");
  }

  // minimisation algorithm:

  algorithm = (params_get_parameter("min_algorithm"))->value.s;

  if (!strcmp(algorithm,"sd")) {

    minimise_settings->algorithm = STEEPEST_DESCENT;

  } else if (!strcmp(algorithm,"cg")) {

    minimise_settings->algorithm = CONJUGATE_GRADIENT;

  } else {

    error_fn("init_minimise_settings: no such algorithm '%s'",algorithm);
  }

  // degrees of freedom:

  minimise_settings->space = dofs_get_space((params_get_parameter("min_space"))->value.s);

  // other parameters:

  minimise_settings->max_iter = (params_get_parameter("min_max_iter"))->value.i;
  minimise_settings->tether_atoms = (params_get_parameter("min_tether_atoms"))->value.i;
  minimise_settings->tether_k = (params_get_parameter("min_tether_k"))->value.d;
}



static int take_minimise_step(DOF_LIST *list,SYSTEM *system,enum MINIMISATION_ALGORITHM algorithm,PLI_SFUNC *sfunc) {

  int flag;

  flag = 0;

  if (algorithm == STEEPEST_DESCENT) {

    flag = take_minimise_step_sd(list,system,sfunc);

  } else if (algorithm == CONJUGATE_GRADIENT) {

    flag = take_minimise_step_cg(list,system,sfunc);
  }

  return(flag);
}



static int take_minimise_step_sd(DOF_LIST *list,SYSTEM *system,PLI_SFUNC *sfunc) {

  set_dof_list_gradients(list,system,sfunc);

  if (!check_gradients(list)) {

    return(1);
  }

  if (minimise_along_gradient(list,system,sfunc)) {

    return(2);
  }

  return(0);
}



static int take_minimise_step_cg(DOF_LIST *list,SYSTEM *system,PLI_SFUNC *sfunc) {

  set_dof_list_gradients(list,system,sfunc);

  if (!check_gradients(list)) {

    return(1);
  }

  if (set_dof_list_conjugate_gradients(list,system)) {

    return(2);
  }

  if (minimise_along_gradient(list,system,sfunc)) {

    return(3);
  }

  return(0);
}



static int check_gradients(DOF_LIST *list) {

  int i;
  DOF *variable;

  for (i=0,variable=list->variables;i<list->n_variables;i++,variable++) {

    if (fabs(variable->fd) > 1.0E-10) {

      return(1);
    }
  }

  return(0);
}



static int minimise_along_gradient(DOF_LIST *list,SYSTEM *system,PLI_SFUNC *sfunc) {

  int i;
  double shift,f;
  MIN_BRACKET *bracket;
  DOF *variable;

  bracket = bracket_minimum(list,system,sfunc);

  if (!bracket->error) {

    shift = golden_section_search(list,system,bracket,sfunc);

    f = score_line_shift(list,system,shift,sfunc);

    set_dof_list_shifts(list,shift);

    apply_dof_list_shifts(list);

    score_system(system,sfunc);
    
    free(bracket);

    return(0);

  } else {

    free(bracket);

    return(1);
  }
}



static MIN_BRACKET* bracket_minimum(DOF_LIST *list,SYSTEM *system,PLI_SFUNC *sfunc) {

  int n;
  double fa,fb,fc,ax,bx,cx;
  MIN_BRACKET *bracket;

  bracket = (MIN_BRACKET*) malloc(sizeof(MIN_BRACKET));

  if (bracket == NULL) {

    error_fn("bracket_minimum: out of memory allocating bracket");
  }

  init_bracket(bracket,list);

  ax = bracket->ax;
  bx = bracket->bx;
  cx = bracket->cx;

  fa = system->score;
  fb = score_line_shift(list,system,bx,sfunc);

  if (fb > fa) {

    bracket->error = 1;

    return(bracket);
  }

  fc = score_line_shift(list,system,cx,sfunc);

  n = 0;

  while (fb >= fc) {

    ax = bx;
    bx = cx;
    cx = bx + (1.0+GOLDEN_RATIO)*(bx-ax);

    fa = fb;
    fb = fc;
    fc = score_line_shift(list,system,cx,sfunc);

    if (n > 10000) {

      bracket->error = 2;

      return(bracket);
    }

    n++;
  };

  bracket->ax = ax;
  bracket->bx = bx;
  bracket->cx = cx;

 return(bracket);
}



static double golden_section_search(DOF_LIST *list,SYSTEM *system,MIN_BRACKET *bracket,PLI_SFUNC *sfunc) {

  int i,j,n;
  double ax,bx,cx,x0,x1,x2,x3,f1,f2,C;
  DOF *variable;
  ATOM *atom;

  C = 1.0 - GOLDEN_RATIO;

  ax = bracket->ax;
  bx = bracket->bx;
  cx = bracket->cx;

  x0 = bracket->ax;
  x3 = bracket->cx;

  if (fabs(cx-bx) > fabs(bx-ax)) {

    x1 = bx;
    x2 = bx + C*(cx-bx);

  } else {

    x2 = bx;
    x1 = bx - C*(bx-ax);
  }

  f1 = score_line_shift(list,system,x1,sfunc);
  f2 = score_line_shift(list,system,x2,sfunc);

  n = 0;

  while ((fabs(x3-x0) > GOLDEN_SECTION_TOLERANCE*(fabs(x1)+fabs(x2))) && (fabs(x3-x0) > 1.0E-10)) {

    if (f2 < f1) {

      x0 = x1;
      x1 = x2;
      x2 = GOLDEN_RATIO*x1 + C*x3;
      
      f1 = f2;
      f2 = score_line_shift(list,system,x2,sfunc);

    } else {

      x3 = x2;
      x2 = x1;
      x1 = GOLDEN_RATIO*x2 + C*x0;

      f2 = f1;
      f1 = score_line_shift(list,system,x1,sfunc);
    }

    n++;
  }

  if (f1 < f2) {

    return(x1);
  }

  return(x2);
}



static void init_bracket(MIN_BRACKET *bracket,DOF_LIST *list) {

  int i;
  double bx2;
  DOF *variable;

  bracket->ax = 0.0;

  bx2 = 0.0;

  for (i=0,variable=list->variables;i<list->n_variables;i++,variable++) {

    bx2 += sqr(variable->fd_shift);
  }

  bracket->bx = sqrt(bx2);
  bracket->cx = (1.0+GOLDEN_RATIO)*(bracket->bx);

  bracket->error = 0;
}



static double score_line_shift(DOF_LIST *list,SYSTEM *system,double shift,PLI_SFUNC *sfunc) {

  int i;
  DOF *variable;
  double score1,score2;

  set_dof_list_shifts(list,shift);

  score1 = system->score;

  store_positions(list->atomlist,list->pos,list->u,list->v,list->w);

  apply_dof_list_shifts(list);

  score_system(system,sfunc);

  score2 = system->score;

  restore_positions(list->atomlist,list->pos,list->u,list->v,list->w);

  system->score = score1;

  return(score2);
}



static void set_dof_list_gradients(DOF_LIST *list,SYSTEM *system,PLI_SFUNC *sfunc) {

  int i;
  DOF *variable;

  for (i=0,variable=list->variables;i<list->n_variables;i++,variable++) {

    set_dof_score_gradient(variable,system,sfunc);
  }
}



static int set_dof_list_conjugate_gradients(DOF_LIST *list,SYSTEM *system) {

  int i;
  double gg,dgg,gamma;
  DOF *variable;

  // initialise if needed:

  if (!list->cg_initialised) {

    for (i=0,variable=list->variables;i<list->n_variables;i++,variable++) {

      variable->fd_prev = variable->fd;
      variable->cg = variable->fd;
    }

    list->cg_initialised = 1;
  }

  // calculate Polak-Ribiere gamma for conjugate gradient:

  gg = 0.0;
  dgg = 0.0;

  for (i=0,variable=list->variables;i<list->n_variables;i++,variable++) {

    gg += sqr(variable->fd_prev);
    dgg += (variable->fd + variable->fd_prev)*(variable->fd);
  }
  if (gg < 1.0E-20) {

    return(1);
  }
  
  gamma = dgg/gg;

  // calculate conjugate gradients:

  for (i=0,variable=list->variables;i<list->n_variables;i++,variable++) {

    variable->cg = variable->fd + gamma*(variable->cg);
  }
  
  return(0);
}
