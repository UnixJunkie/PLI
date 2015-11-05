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



static ATOM *atom_tree[MAX_RING_SIZE+1];



static void grow_ring_tree(ATOM*,ATOM*,int,int);
static int aromatic_ring(RING*);
static int flat_ring(RING*);
static RING_LIST* alloc_ring_list(void);
static void init_ring_list(RING_LIST*);
static void add_ring_to_list(RING_LIST*,RING*);
static int ring_in_list(RING_LIST*,RING*);
static int same_ring(RING*,RING*);



RING_LIST* atomlist2ringlist(ATOMLIST *atom_list,int max_ring_size) {

  int i;
  ATOM **atomp,*atom;
  RING *ring;
  RING_LIST *ring_list;

  ring_list = alloc_ring_list();

  for (i=0,atomp=atom_list->atom;i<atom_list->natoms;i++,atomp++) {

    atom = *atomp;

    find_ring(atom,max_ring_size);

    ring = atom->ring;

    if ((ring->size > 0) && (!ring_in_list(ring_list,ring))) {

      add_ring_to_list(ring_list,ring);
    }
  }

  return(ring_list);
}



void find_ring(ATOM *atom,int max_ring_size) {
  
  int i;
  ATOM **batom,**ratom;
  MAP *covalent_map;

  if (max_ring_size > MAX_RING_SIZE) {

    error_fn("find_ring: maximum possible ring size exceeded");
  }

  if (atom->ring != NULL) {

    return;
  }

  if (atom->ring == NULL) {

    atom->ring = (RING*) malloc(sizeof(RING));
  }

  if (atom->ring == NULL) {

    error_fn("find_rings: out of memory allocating ring");
  }

  atom->ring->size = 0;
  atom->ring->aromatic = 0;

  if (atom->connections != NULL) {

    if (atom->connections->natoms > 1) {

      atom_tree[0] = atom;

      for (i=0,batom=atom->connections->atom;i<atom->connections->natoms;i++,batom++) {

	if (atom->element->id != HYDROGEN) {

	  grow_ring_tree(*batom,atom,1,max_ring_size);
	}
      }
    }
  }

  if ((atom->ring->size == 5) || (atom->ring->size == 6)) {

    if (aromatic_ring(atom->ring)) {

      atom->ring->aromatic = 1;

      for (i=0,ratom=atom->ring->atom;i<atom->ring->size;i++,ratom++) {

	(*ratom)->flags |= AROMATIC_ATOM;
      }
    }
  }
}



void calc_ring_center(RING *ring) {

  int i,j;
  double *v;
  ATOM **atomp,*atom;

  if (ring->size == 0) {

    error_fn("calc_ring_center: no atoms in ring");
  }

  null_vector(ring->center);

  for (i=0,atomp=ring->atom;i<ring->size;i++,atomp++) {

    atom = *atomp;

    v = atom->position;

    for (j=0;j<3;j++) {

      ring->center[j] += v[j];
    }
  }

  for (i=0;i<3;i++) {

    ring->center[i] /= ring->size;
  }
}



void calc_ring_normal(RING *ring) {

  int i;
  double *center,*normal,*pos1,*pos2,v1[4],v2[4],cp[4];
  ATOM *atom1,*atom2;

  center = ring->center;
  normal = ring->normal;

  null_vector(ring->normal);

  for (i=0;i<ring->size;i++) {

    atom1 = ring->atom[i];
    atom2 = ring->atom[(i+1) % ring->size];

    pos1 = atom1->position;
    pos2 = atom2->position;

    calc_vector(center,pos1,v1);
    calc_vector(center,pos2,v2);

    calc_crossproduct(v1,v2,cp);

    sum_vector(normal,cp,normal);
  }

  scale_vector(normal,1.0);
}



static void grow_ring_tree(ATOM *atom,ATOM *prev_atom,int tree_pos,int max_ring_size) {

  int i,size;
  ATOM **batom;

  if (tree_pos > max_ring_size) {

    return;
  }

  if (tree_pos > 1) {

    for (i=1;i<tree_pos;i++) {

      if (atom_tree[i] == atom) {

	return;
      }
    }
  }

  if (atom == atom_tree[0]) {

    if ((atom_tree[0]->ring->size == 0) || (tree_pos < atom_tree[0]->ring->size)) {

      atom_tree[0]->ring->size = tree_pos;

      for (i=0;i<tree_pos;i++) {

	atom_tree[0]->ring->atom[i] = atom_tree[i];
      }
    }

    return;
  }

  atom_tree[tree_pos] = atom;

  if (atom->connections == NULL ) {

    return;
  }

  if (atom->connections->natoms == 0) {

    return;
  }

  for (i=0,batom=atom->connections->atom;i<atom->connections->natoms;i++,batom++) {

    if ((*batom)->element->id != HYDROGEN) {

      if (*batom != prev_atom) {

	grow_ring_tree(*batom,atom,tree_pos+1,max_ring_size);
      }
    }
  }
}



static int aromatic_ring(RING *ring) {

  int i,el;
  int n_carbon,n_nitrogen,n_oxygen,n_sulfur;
  ATOM_NODE *node;
  ATOM **atomp,*atom;

  n_carbon = n_nitrogen = n_oxygen = n_sulfur = 0;

  for (i=0,atomp=ring->atom;i<ring->size;i++,atomp++) {

    atom = *atomp;

    node = atom->node;

    el = atom->element->id;

    if (el == CARBON) {

      if ((!node) || (node->hybridisation != 2)) {

	return(0);
      }

      n_carbon++;

    } else if (el == NITROGEN) {

      n_nitrogen++;

    } else if (el == OXYGEN) {

      n_oxygen++;

    } else if (el == SULFUR) {

      n_sulfur++;
    }
  }

  if (n_carbon + n_nitrogen + n_oxygen + n_sulfur != ring->size) {

    return(0);
  }

  if (ring->size == 5) {

    if (n_oxygen + n_sulfur > 1) {

      return(0);
    }

    return(1);

  } else if (ring->size == 6) {

    if (n_oxygen + n_sulfur > 0) {

      return(0);
    }
    
    return(1);
  }

  return(0);
}



static int flat_ring(RING *ring) {

  int i;
  double torsion;
  ATOM **atom1,**atom2,**atom3,**atom4;

  for (i=0;i<ring->size;i++) {

    atom1 = ring->atom + i;
    atom2 = (i+1<ring->size) ? ring->atom + i + 1 : ring->atom + i + 1 - ring->size;
    atom3 = (i+2<ring->size) ? ring->atom + i + 2 : ring->atom + i + 2 - ring->size;
    atom4 = (i+3<ring->size) ? ring->atom + i + 3 : ring->atom + i + 3 - ring->size;

    torsion = (180/PI)*torsion_angle((*atom1)->position,
				     (*atom2)->position,
				     (*atom3)->position,
				     (*atom4)->position);

    if (fabs(torsion) > 6.0)
      return(0);
  }

  return(1);
}



static RING_LIST* alloc_ring_list(void) {

  RING_LIST *list;

  list = (RING_LIST*) malloc(sizeof(RING_LIST));

  if (list == NULL) {

    error_fn("alloc_ring_list: out of memory allocating ring list");
  }

  init_ring_list(list);

  return(list);
}



static void init_ring_list(RING_LIST *list) {

  list->n_rings = 0;
  list->rings = NULL;
  list->n_alloc_rings = 0;
}



static int ring_in_list(RING_LIST *list,RING *ring) {

  int i;
  RING **ringp;

  for (i=0,ringp=list->rings;i<list->n_rings;i++,ringp++) {

    if (same_ring(ring,*ringp)) {

      return(1);
    }
  }

  return(0);
}



static void add_ring_to_list(RING_LIST *list,RING *ring) {

  if (list->n_alloc_rings == 0) {

    list->n_alloc_rings = 10;

    list->rings = (RING**) calloc(list->n_alloc_rings,sizeof(RING*));
    
  } else if (list->n_rings == list->n_alloc_rings) {

    list->n_alloc_rings *= 2;

    list->rings = (RING**) realloc(list->rings,(list->n_alloc_rings)*sizeof(RING*));
  }

  if (list->rings == NULL) {

    error_fn("add_ring_to_list: out of memory allocating rings");
  }

  list->rings[list->n_rings] = ring;

  list->n_rings++;
}



void free_ring_list(RING_LIST *ring_list) {

  if (ring_list) {

    if (ring_list->rings) {

      free(ring_list->rings);
    }

    free(ring_list);
  }
}



static int same_ring(RING *ring1,RING *ring2) {

  int i,j;
  ATOM **atom1p,**atom2p,*atom1,*atom2;

  if (ring1->size != ring2->size) {

    return(0);
  }

  for (i=0,atom1p=ring1->atom;i<ring1->size;i++,atom1p++) {

    atom1 = *atom1p;
    atom2 = NULL;

    for (j=0,atom2p=ring2->atom;(j<ring2->size)&&(atom1!=atom2);j++,atom2p++) {

      atom2 = *atom2p;
    }

    if (atom1 != atom2) {

      return(0);
    }
  }

  return(1);
}
