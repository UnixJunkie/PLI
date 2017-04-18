// Copyright  2005 CCDC Software Ltd and 2017 Astex Therapeutics Ltd.
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

/*******************************************************************************
*  This code was adapted directly from ullman.c source file from GOLD program
*  and is © CCDC Software Ltd 2005.
*  CCDC Software Ltd's address is 12 Union Road, Cambridge CB2 1EZ, UK.
******************************************************************************/



#include "pli.h"

#define LIMIT_ISOMORPHISMS
#define MAX_ISOMORPHISMS 10000

static int TIME_LIMIT = 1;
static int MAX_TIME = 180; // seconds

static BOOLEAN **ull_mat, **a_mat, **b_mat;
static int palpha, pbeta, isomorphisms;
/* static int ull_flags; */
static void (*app_func) (BOOLEAN **);

// Limit the time spent in the ullman algorithm.
static time_t start_time;
static int depthfirst_counter = 0;
static int times_up = 0;

static ULLMAN_MATCH_LIST *ullman_list = NULL;
static ULLMAN_MAPPINGS *ullman_mappings = NULL;
static enum MATCH_PROTOCOL ullman_match_protocol = MATCH_ALL;
static enum MATCH_PROTOCOL ullman_match_quality = MATCH_ATOM_TYPE;

static MOLECULE *ullman_query = NULL;
static MOLECULE *ullman_target= NULL;

static void init_ullman_settings(MOLECULE*, MOLECULE*, enum MATCH_PROTOCOL, enum MATCH_QUALITY);
static void init_ullman_matches(MOLECULE*);
static void init_ullman_mappings(MOLECULE*, MOLECULE*);
static void reset_ullman_settings(void);

static void finish_ullman_matches(void);
static void record_ullman_match(BOOLEAN**);
static void record_ullman_mapping(BOOLEAN**);
static BOOLEAN **create_ullman_matrix (MOLECULE *, MOLECULE *);
static int ullman_atom_match(ATOM*,ATOM*);
static void depthfirst (int, BOOLEAN **, BOOLEAN *);
static void refine (BOOLEAN **, int *);
static int row_and_column (BOOLEAN **m, int x, int j);
static int corresponds (BOOLEAN **, int, int);
static int all_false (BOOLEAN **, int);
static void dataout (BOOLEAN **);



ULLMAN_MATCH_LIST* find_ullman_matches(MOLECULE *query,MOLECULE *target,enum MATCH_PROTOCOL protocol,enum MATCH_QUALITY quality) {

  init_ullman_settings(query, target, protocol, quality);

  init_ullman_matches(query);

  do_ullman(query,target,record_ullman_match);

  finish_ullman_matches();

  reset_ullman_settings();

  return(ullman_list);
}


ULLMAN_MAPPINGS* find_ullman_mappings(MOLECULE *query,MOLECULE *target,enum MATCH_PROTOCOL protocol,enum MATCH_QUALITY quality) {

  init_ullman_settings(query,target,protocol,quality);

  init_ullman_mappings(query, target);

  do_ullman(query,target,record_ullman_mapping);

  //finish_ullman_mappings();

  reset_ullman_settings();

  return ullman_mappings;
}

/*
do_ullman

  This is the entry point to the Ullman algorithm.  Mola is the query
  molecule and molb the target. Func is called when a match is found.
Match_type allow the seeting of some options. */
void do_ullman (MOLECULE *mola, MOLECULE *molb, void (*func) (BOOLEAN **)) {
  BOOLEAN *k, *p;
  int  i;

  // Initialise static early exit vars if the Ullman is taking too long.
  depthfirst_counter = 0;
  time(&start_time);
  times_up = 0;

  /* Initilialisation */
  if (mola->selection == NULL || molb->selection == NULL) {
    error_fn("%s: at least one of the molecules has not got selection list allocated", __func__);
  }
  if ((!mola->selection->natoms) || (!molb->selection->natoms)) {
    error_fn("%s: at least one of the molecule has no atoms selected", __func__);
  }
  if (mola->selection->natoms > molb->selection->natoms) {
    error_fn ("%s: query size larger than target\n", __func__);
  }

  a_mat = create_adjacency_matrix(mola);
  b_mat = create_adjacency_matrix(molb);
  /* don't count dummy atoms as neighbours */
  ull_mat = create_ullman_matrix(mola, molb);
  palpha = mola->selection->natoms;
  pbeta = molb->selection->natoms;

  app_func = func;
  k = alloc_bool_array(pbeta);
  isomorphisms = 0;

  for (i=0, p=k; i<pbeta; i++, p++) {
    *p = 0;
  }
  /* do the business */

  depthfirst (0, ull_mat, k);
  //warning_fn ("%d Isomorphisms found\n", isomorphisms);

  free_2d_bool_matrix(a_mat);
  free_2d_bool_matrix(b_mat);
  free_2d_bool_matrix(ull_mat);
  free(k);
}


static void init_ullman_settings(MOLECULE *query,MOLECULE *target,enum MATCH_PROTOCOL protocol,enum MATCH_QUALITY quality) {
  ullman_match_protocol = protocol;
  ullman_match_quality = quality;
  ullman_query = query;
  ullman_target = target;
}


static void init_ullman_matches(MOLECULE *query) {

  ATOMLIST *selection;

  selection = query->selection;

  if ((!selection) || (selection->natoms == 0)) {

    error_fn("%s: no query atoms",__func__);
  }

  if (selection->natoms > MAX_MATCH_ATOMS) {

    error_fn("%s: exceeded maximum number of query atoms (%d)",__func__,MAX_MATCH_ATOMS);
  }

  ullman_list = (ULLMAN_MATCH_LIST*) malloc(sizeof(ULLMAN_MATCH_LIST));

  if (ullman_list == NULL) {

    error_fn("%s: out of memory allocating ullman list",__func__);
  }

  ullman_list->n_matches = 0;
  ullman_list->n_atoms = selection->natoms;
  ullman_list->n_alloc_matches = 100;

  ullman_list->matches = (MATCH*) calloc(ullman_list->n_alloc_matches,sizeof(MATCH));

  if (ullman_list->matches == NULL) {

    error_fn("%s: out of memory allocating ullman matches",__func__);
  }
}


static void reset_ullman_settings(void) {
  ullman_match_protocol = MATCH_ALL;
  ullman_match_quality = MATCH_ATOM_TYPE;
}

static void finish_ullman_matches(void) {
  if (ullman_list->n_matches == 0) {

    free(ullman_list->matches);

    free(ullman_list);

    ullman_list = NULL;
  }
}


static void record_ullman_match(BOOLEAN **m) {

  int i,j,atom_id;
  MATCH *match;

  if (ullman_list->n_matches == ullman_list->n_alloc_matches) {
    
    ullman_list->n_alloc_matches += 100;

    ullman_list->matches = (MATCH*) realloc(ullman_list->matches,(ullman_list->n_alloc_matches)*sizeof(MATCH));

    if (ullman_list->matches == NULL) {

      error_fn("%s: out of memory reallocating ullman matches",__func__);
    }
  }

  match = ullman_list->matches + ullman_list->n_matches;

  match->n_atoms = ullman_list->n_atoms;

  match->score = 0;
  match->flags = 00;

  atom_id = 0;

  for (i=0; i<palpha; i++) {

    for (j=0; j<pbeta; j++) {

      if (m[i][j]) {

	if (atom_id == ullman_list->n_atoms) {

	  error_fn("%s: more matching atoms than expected",__func__);
	}

	match->atoms1[atom_id] = *(ullman_query->selection->atom + i);
	match->atoms2[atom_id] = *(ullman_target->selection->atom + j);
	//match[atom_id][0] = i;
	//match[atom_id][1] = j;

	atom_id++;
      }
    }
  }

  ullman_list->n_matches++;
}

static void init_ullman_mappings(MOLECULE *query, MOLECULE *target) {
  if ((!query->selection) || (query->selection->natoms == 0) || (!target->selection) || (target->selection->natoms == 0)) {
    error_fn("One of the molecules has no atoms selected");
  }

  ullman_mappings = (ULLMAN_MAPPINGS*) malloc(sizeof(ULLMAN_MAPPINGS));
  ullman_mappings->n_atoms1 = query->selection->natoms;
  ullman_mappings->n_atoms2 = target->selection->natoms;
  ullman_mappings->n_mappings = 0;
  ullman_mappings->n_allocated = 0;
  ullman_mappings->mappings = NULL;
}

static void record_ullman_mapping(BOOLEAN **m) {
  MOLECULE *molecule;
  BOOLEAN *automorph_matrix;
  int block_size;

  block_size = 2;

  if (ullman_mappings->n_mappings == ullman_mappings->n_allocated) {
    ullman_mappings->mappings = (BOOLEAN***) realloc(ullman_mappings->mappings,
                                                     (ullman_mappings->n_mappings + block_size) * sizeof(BOOLEAN **));
    ullman_mappings->n_allocated += block_size;
    for (int i=0; i<block_size; i++) {
      ullman_mappings->mappings[ullman_mappings->n_mappings + i] = alloc_2d_bool_matrix(ullman_mappings->n_atoms1, ullman_mappings->n_atoms2);
    }
  }

  copy_2d_bool_matrix(ullman_mappings->mappings[ullman_mappings->n_mappings], m, ullman_mappings->n_atoms1, ullman_mappings->n_atoms2);
  ullman_mappings->n_mappings++;
}


static BOOLEAN **create_ullman_matrix (MOLECULE *mola, MOLECULE *molb)
{
  BOOLEAN **matrix;

  ATOMLIST *lista, *listb;
  ATOM *atom1, *atom2;

  int i, j;

  lista = mola->selection;
  listb = molb->selection;

  matrix = alloc_2d_bool_matrix (lista->natoms, listb->natoms);

  /* Set matrix[i][j] if atom i from the query and atom j from the
  target have the same atom type and atom j has the same of more
  neighbours as atom i. */

  for (i=0; i<lista->natoms; i++) {

    atom1 = lista->atom[i];

    for (j=0; j<listb->natoms; j++) {

      atom2 = listb->atom[j];

      matrix[i][j] = ullman_atom_match(atom1,atom2);
    }
  }
  return (matrix);
}



static int ullman_atom_match(ATOM *atom1,ATOM *atom2) {

  if (ullman_match_quality == MATCH_ATOM_TYPE) {

      if ((atom1->connections->natoms <= atom2->connections->natoms) &&
          (atom1->type == atom2->type)) {

	return(1);
      }
  } else if (ullman_match_quality == MATCH_NODE) {

  } else if (ullman_match_quality == MATCH_ELEMENT) {

    if (atom1->element == atom2->element) {

      return(1);
    }

  } else if (ullman_match_quality == MATCH_ANY) {

    return(1);
  }

  return(0);
}


/*
depthfirst(), refine(), row_and_column(),
corresponds(), all_false(), dataout()

  are C versions of the original fortran/pascal routines.  If you
  want to know more about these then have a look at some of Pete
Willetts publications.  */

static void depthfirst (int d, BOOLEAN **m, BOOLEAN *f)
{
  int node, j, finish;
  BOOLEAN **temp1;

  finish = 0;
  temp1 = alloc_2d_bool_matrix (palpha, pbeta);
  copy_2d_bool_matrix(temp1, m, palpha, pbeta);
  node = 0;

  for (node =0; node<pbeta; node++) {

    if (m[d][node] && !f[node]) {
      /* m at this depth and node has not been visited */

      for (j=0; j<pbeta; j++)
        m[d][j] = 0;

      m[d][node] = 1;
      f[node] = 1;

      refine (m, &finish);
#ifdef LIMIT_ISOMORPHISMS

      ++depthfirst_counter;
      if( TIME_LIMIT &&
          depthfirst_counter % 1000 == 0 &&
          time(NULL) - start_time >= MAX_TIME )
      {
        warning_fn("Exiting Ullman subgraph isomorphism algorithm. "
            "Taking too long (%d secs so far).", time(NULL) - start_time);
        times_up = 1;
      }

      if (!finish && isomorphisms < MAX_ISOMORPHISMS && !times_up) {
#else
        if (!finish) {
#endif
          if (d == (palpha-1)) {
            isomorphisms++;
            //dataout (m);

	    if (app_func) {

	      (*app_func) (m);
	    }
          }
          else {
            depthfirst (d+1, m, f);
          }
        }

        f[node] = 0;

        copy_2d_bool_matrix(m, temp1, palpha, pbeta);
      }

    }

    free_2d_bool_matrix(temp1);
  }




static void refine (BOOLEAN **m, int *finish)
{
  int i, j, no_change;

  *finish = 0;
  no_change = 0;

  while (!no_change) {
    no_change = 1;

    for (i=0; i<palpha; i++) {
      for (j=0; j<pbeta; j++) {

        if (m[i][j] && !corresponds (m, i, j)) {

          m[i][j] = 0;
          no_change = 0;
          if (all_false (m, i)) {
            *finish = 1;
            return;
          }
        }
      }
    }
  }
}


static int row_and_column (BOOLEAN **m, int x, int j)
{
  int y;

  for (y=0; y<pbeta; y++)
    if (m[x][y] && b_mat[y][j])
      return 1;

  return 0;
}


static int corresponds (BOOLEAN **m, int i, int j)
{
  int x;

  for (x=0; x<palpha; x++)
    if (a_mat[i][x] && !row_and_column (m, x, j))
      return 0;

  return 1;
}


static int all_false (BOOLEAN **m, int row)
{
  int col;

  for (col=0; col<pbeta; col++)
    if (m[row][col])
      return 0;

  return 1;
}

static void dataout (BOOLEAN **m)
{
  int i, j;

  for (i=0; i<palpha; i++) {
    printf ("atom no %2d --> ", i+1);
    for (j=0; j<pbeta; j++) {
      if (m[i][j])
        printf ("%2d ", j+1);
    }
    printf ("\n");
  }
}
