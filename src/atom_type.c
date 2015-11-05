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



static ATOM_TYPING_SCHEME *SCHEME = NULL;



static void check_atom_type(ATOM*,ATOM_TYPING_SCHEME*);
static void check_atom_type_vs_hydrogens(ATOM*,ATOM_TYPING_SCHEME*);
static ATOM_TYPE* tautomer_atom_type(ATOM*,ATOM_TYPING_SCHEME*);
static ATOM_TYPE* derive_atom_type(ATOM*,ATOM_TYPING_SCHEME*);
static int get_atom_set(ATOM*);
static int charged_aromatic_n(ATOM*);
static ATOM_TYPE* ambiguous_aromatic_n(ATOM*,ATOM_TYPING_SCHEME*);
static int ambiguous_pyridone_n(ATOM*);
static int ambiguous_pyridone_o(ATOM*);
static int ambiguous_pyridine_n(ATOM*);
static int pyridine_ring(RING*);
static int aromatic_nc_6_ring(RING*);
static int pyridone_atom_pair(ATOM*,ATOM*);
static int amino_pyridine_atom_pair(ATOM*,ATOM*);
static int atom_type_match(ATOM*,ATOM_TYPE_DEF*);

static void read_elements(PLI_FILE*,ATOM_TYPING_SCHEME*);
static void read_atom_nodes(PLI_FILE*,ATOM_TYPING_SCHEME*);
static void read_united_atoms(PLI_FILE*,ATOM_TYPING_SCHEME*);
static void read_atom_types(PLI_FILE*,ATOM_TYPING_SCHEME*);
static void read_ambiguous_types(PLI_FILE*,ATOM_TYPING_SCHEME*);
static void read_residue_atoms(PLI_FILE*,ATOM_TYPING_SCHEME*);
static void read_atom_type_def(PLI_FILE*,ATOM_TYPING_SCHEME*);
static void read_atom_node_def(PLI_FILE*,ATOM_TYPE_DEF*,ATOM_TYPING_SCHEME*);

static void finalise_scheme(ATOM_TYPING_SCHEME*);

static void init_scheme(ATOM_TYPING_SCHEME*);
static void init_element(ELEMENT*);
static void init_atom_node(ATOM_NODE*);
static void init_united_atom(UNITED_ATOM*);
static void init_atom_type(ATOM_TYPE*);
static void init_residue_atom(RESIDUE_ATOM*);
static void init_atom_type_def(ATOM_TYPE_DEF*);

static void alloc_elements(ATOM_TYPING_SCHEME*);
static void alloc_atom_nodes(ATOM_TYPING_SCHEME*);
static void alloc_united_atoms(ATOM_TYPING_SCHEME*);
static void alloc_atom_types(ATOM_TYPING_SCHEME*);
static void alloc_residue_atoms(ATOM_TYPING_SCHEME*);
static void alloc_atom_type_defs(ATOM_TYPING_SCHEME*);
static void alloc_atom_node_defs(ATOM_TYPE_DEF*);







void run_type(SETTINGS *settings) {

  int i;
  ATOM **atom;
  ATOMLIST *selection;
  SYSTEM *system;
  PLI_FILE *file;
  ATOM_TYPING_SCHEME *scheme;

  system = settings2system(settings);

  prep_system(system,ANY_MOLECULE);

  output_system(PLI_STDOUT,system);

  unprep_system(system,ANY_MOLECULE);

  free_system(system);
}



void read_atom_typing_scheme(char *filename,ATOM_TYPING_SCHEME *scheme) {

  int flag;
  char line[MAX_LINE_LEN],word[MAX_LINE_LEN];
  PLI_FILE *file;

  file = open_file(filename,"r");

  if (file == NULL) {

    error_fn("couldn't open param file '%s'",filename);
  }

  while (!end_of_file(file)) {

    if (read_line(line,MAX_LINE_LEN,file) == NULL)
      break;

    flag = sscanf(line,"%s",word);

    if (flag != EOF) {

      if (!strcmp(word,"elements")) {

	read_elements(file,scheme);

      } else if (!strcmp(word,"atom_nodes")) {

	read_atom_nodes(file,scheme);

      } else if (!strcmp(word,"united_atoms")) {

	read_united_atoms(file,scheme);

      } else if (!strcmp(word,"atom_geometries")) {

	read_atom_geometries(file,scheme);

      } else if (!strcmp(word,"atom_types")) {

	read_atom_types(file,scheme);

      } else if (!strcmp(word,"ambiguous_types")) {

	read_ambiguous_types(file,scheme);

      } else if (!strcmp(word,"residue_atoms")) {

	read_residue_atoms(file,scheme);

      } else if (!strcmp(word,"atom_type_def")) {

	read_atom_type_def(file,scheme);
      }
    }
  }

  close_file(file);

  finalise_scheme(scheme);

  SCHEME = scheme;
}



static void read_elements(PLI_FILE *paramfile,ATOM_TYPING_SCHEME *scheme) {

  int flag,n_read;
  char line[MAX_LINE_LEN],word[MAX_LINE_LEN];
  ELEMENT *element;

  strcpy(word,"");

  while ((!end_of_file(paramfile)) && (strcmp(word,"end_elements"))) {
  
    if (read_line(line,MAX_LINE_LEN,paramfile) == NULL)
      break;

    flag = sscanf(line,"%s",word);

    if (flag != EOF) {

      if (strcmp(word,"end_elements")) {

	alloc_elements(scheme);

	element = scheme->elements + scheme->n_elements;

	init_element(element);

	n_read = sscanf(line,"%d \"%[^\"]\" %lf %lf %d",&(element->id),element->name,&(element->cov_radius),&(element->vdw_radius),&(element->flags));

	if (n_read == 5) {

	  scheme->n_elements++;

	} else {

	  warning_fn("read_elements: skipping corrupted element line:\n%s",line);
	}
      }

    } else {

      strcpy(word,"");
    }
  }
}



static void read_atom_nodes(PLI_FILE *paramfile,ATOM_TYPING_SCHEME *scheme) {

  int flag,n_read;
  char line[MAX_LINE_LEN],word[MAX_LINE_LEN];
  ATOM_NODE *atom_node;

  strcpy(word,"");

  while ((!end_of_file(paramfile)) && (strcmp(word,"end_atom_nodes"))) {
  
    if (read_line(line,MAX_LINE_LEN,paramfile) == NULL)
      break;

    flag = sscanf(line,"%s",word);

    if (flag != EOF) {

      if (strcmp(word,"end_atom_nodes")) {

	alloc_atom_nodes(scheme);

	atom_node = scheme->atom_nodes + scheme->n_atom_nodes;

	init_atom_node(atom_node);

	n_read = sscanf(line,"%*d \"%[^\"]\" \"%[^\"]\" %d %d %d %d %d",
			atom_node->name,atom_node->element_name,
			&(atom_node->n_single),&(atom_node->n_double),&(atom_node->n_triple),
			&(atom_node->hybridisation),&(atom_node->formal_charge));

	if (n_read == 7) {

	  scheme->n_atom_nodes++;

	} else {

	  warning_fn("read_atom_nodes: skipping corrupted atom node line:\n%s",line);
	}
      }

    } else {

      strcpy(word,"");
    }
  }
}



static void read_united_atoms(PLI_FILE *paramfile,ATOM_TYPING_SCHEME *scheme) {

  int flag,n_read;
  char line[MAX_LINE_LEN],word[MAX_LINE_LEN];
  UNITED_ATOM *united_atom;

  strcpy(word,"");

  while ((!end_of_file(paramfile)) && (strcmp(word,"end_united_atoms"))) {
  
    if (read_line(line,MAX_LINE_LEN,paramfile) == NULL)
      break;

    flag = sscanf(line,"%s",word);

    if (flag != EOF) {

      if (strcmp(word,"end_united_atoms")) {

	alloc_united_atoms(scheme);

	united_atom = scheme->united_atoms + scheme->n_united_atoms;

	init_united_atom(united_atom);

	n_read = sscanf(line,"%d \"%[^\"]\" \"%[^\"]\" %d %d %d %lf %lf %lf",
			&(united_atom->id),united_atom->name,united_atom->atom_node_name,
			&(united_atom->hybridisation),&(united_atom->n_hydrogens),
			&(united_atom->formal_charge),&(united_atom->vdw_radius),
			&(united_atom->ideal_hbond_alpha),&(united_atom->ideal_hbond_beta));

	if (n_read == 9) {

	  scheme->n_united_atoms++;

	} else {

	  warning_fn("read_united_atoms: skipping corrupted united atom line:\n%s",line);
	}
      }

    } else {

      strcpy(word,"");
    }
  }
}



static void read_atom_types(PLI_FILE *paramfile,ATOM_TYPING_SCHEME *scheme) {

  int flag,n_read;
  char line[MAX_LINE_LEN],word[MAX_LINE_LEN];
  ATOM_TYPE *atom_type;

  strcpy(word,"");

  while ((!end_of_file(paramfile)) && (strcmp(word,"end_atom_types"))) {
  
    if (read_line(line,MAX_LINE_LEN,paramfile) == NULL)
      break;

    flag = sscanf(line,"%s",word);

    if (flag != EOF) {

      if (strcmp(word,"end_atom_types")) {

	alloc_atom_types(scheme);

	atom_type = scheme->atom_types + scheme->n_atom_types;

	init_atom_type(atom_type);

	n_read = sscanf(line,"%d \"%[^\"]\" \"%[^\"]\" \"%[^\"]\" %d %*d %*d %*d %d %d %s",
			&(atom_type->id),atom_type->name,
			atom_type->united_atom_name,atom_type->atom_geometry_name,
			&(atom_type->flags),
			&(atom_type->n_hs),&(atom_type->n_lps),atom_type->hbond_geometry_name);

	if (n_read == 8) {

	  scheme->n_atom_types++;

	} else {

	  warning_fn("read_atom_types: skipping corrupted atom type line:\n%s",line);
	}
      }

    } else {

      strcpy(word,"");
    }
  }
}



static void read_ambiguous_types(PLI_FILE *paramfile,ATOM_TYPING_SCHEME *scheme) {

  int flag,n_read,status,assigned_status;
  double probability;
  char line[MAX_LINE_LEN],word[MAX_LINE_LEN];
  char type_str[MAX_LINE_LEN],alt_type_str[MAX_LINE_LEN];
  ATOM_TYPE *type,*alt_type_type;
  ALT_ATOM_TYPE *alt_type;

  strcpy(word,"");

  while ((!end_of_file(paramfile)) && (strcmp(word,"end_ambiguous_types"))) {
  
    if (read_line(line,MAX_LINE_LEN,paramfile) == NULL)
      break;

    flag = sscanf(line,"%s",word);

    if (flag != EOF) {

      if (strcmp(word,"end_ambiguous_types")) {

	n_read = sscanf(line," \"%[^\"]\" \"%[^\"]\" %d %lf %d",type_str,alt_type_str,&status,&probability,&assigned_status);

	if (n_read == 5) {

	  type = get_atom_type(type_str,scheme);

	  if (type) {

	    alt_type_type = get_atom_type(alt_type_str,scheme);

	    if (alt_type_type) {

	      if (type->n_alt_types < MAX_ALT_TYPES) {

		alt_type = type->alt_types + type->n_alt_types;

		alt_type->type = alt_type_type;
		alt_type->group_status = status;
		alt_type->assigned_group_status = assigned_status;
		alt_type->probability = probability;

		type->n_alt_types++;

	      } else {

		warning_fn("read_ambiguous_types: more than %d ambiguous types",MAX_ALT_TYPES);
	      }

	    } else {

	      warning_fn("read_ambiguous_types: no such type '%s'",alt_type_str);
	    }

	  } else {

	    warning_fn("read_ambiguous_types: no such type '%s'",type_str);
	  }

	} else {

	  warning_fn("read_ambiguous_types: skipping corrupted ambiguous type line:\n%s",line);
	}
      }

    } else {

      strcpy(word,"");
    }
  }
}



static void read_residue_atoms(PLI_FILE *paramfile,ATOM_TYPING_SCHEME *scheme) {

  int flag,n_read;
  char line[MAX_LINE_LEN],word[MAX_LINE_LEN];
  RESIDUE_ATOM *residue_atom;

  strcpy(word,"");

  while ((!end_of_file(paramfile)) && (strcmp(word,"end_residue_atoms"))) {
  
    if (read_line(line,MAX_LINE_LEN,paramfile) == NULL)
      break;

    flag = sscanf(line,"%s",word);

    if (flag != EOF) {

      if (strcmp(word,"end_residue_atoms")) {

	alloc_residue_atoms(scheme);

	residue_atom = scheme->residue_atoms + scheme->n_residue_atoms;

	init_residue_atom(residue_atom);

	n_read = sscanf(line,"%d \"%[^\"]\" \"%[^\"]\" \"%[^\"]\"",
			&(residue_atom->id),residue_atom->subname,
			residue_atom->name,residue_atom->atom_type_name);

	if (n_read == 4) {

	  scheme->n_residue_atoms++;

	} else {

	  warning_fn("read_residue_atoms: skipping corrupted residue atom line:\n%s",line);
	}
      }

    } else {

      strcpy(word,"");
    }
  }
}



static void read_atom_type_def(PLI_FILE *paramfile,ATOM_TYPING_SCHEME *scheme) {

  int flag;
  char line[MAX_LINE_LEN],word[MAX_LINE_LEN],atom_type_name[MAX_LINE_LEN];
  ATOM_TYPE *type;
  ATOM_TYPE_DEF *def;

  type = NULL;

  strcpy(word,"");

  while ((!end_of_file(paramfile)) && (strcmp(word,"end_atom_type_def"))) {
  
    if (read_line(line,MAX_LINE_LEN,paramfile) == NULL)
      break;

    flag = sscanf(line,"%s",word);

    if (flag != EOF) {

      if (!strcmp(word,"name")) {

	sscanf(line,"%*s \"%[^\"]\"",atom_type_name);

	type = get_atom_type(atom_type_name,scheme);

	if (type == NULL) {

	  error_fn("read_atom_type_def: undefined atom type '%s'",atom_type_name);
	}

      } else if (!strcmp(word,"atom_node_def")) {

	if (type == NULL) {

	  error_fn("read_atom_type_def: atom type needs to be specified before atom nodes");
	}

	alloc_atom_type_defs(scheme);

	def = scheme->atom_type_defs + scheme->n_atom_type_defs;

	init_atom_type_def(def);

	def->type = type;

	alloc_atom_node_defs(def);

	read_atom_node_def(paramfile,def->defs + def->n_defs,scheme);

	def->n_defs++;

	scheme->n_atom_type_defs++;
      }

    } else {

      strcpy(word,"");
    }
  }
}



static void read_atom_node_def(PLI_FILE *paramfile,ATOM_TYPE_DEF *def,ATOM_TYPING_SCHEME *scheme) {

  int flag;
  char line[MAX_LINE_LEN],word[MAX_LINE_LEN];

  init_atom_type_def(def);

  strcpy(word,"");

  while ((!end_of_file(paramfile)) && (strcmp(word,"end_atom_node_def"))) {

    if (read_line(line,MAX_LINE_LEN,paramfile) == NULL)
      break;

    flag = sscanf(line,"%s",word);

    if (flag != EOF) {

      if (!strcmp(word,"name")) {

	sscanf(line,"%*s \"%[^\"]\"",def->node_name);

	def->node = get_atom_node(def->node_name,scheme);

	if (def->node == NULL) {

	  error_fn("read_atom_node_def: atom node '%s' undefined",def->node_name);
	}

      } else if (!strcmp(word,"element")) {

	sscanf(line,"%*s \"%[^\"]\"",def->element_name);

	def->element = get_element(def->element_name,scheme);

	if (def->element == NULL) {

	  error_fn("read_atom_node_def: element '%s' undefined",def->element_name);
	}
	
      } else if (!strcmp(word,"hybridisation")) {

	sscanf(line,"%*s %d",&(def->hybridisation));

      } else if (!strcmp(word,"aromatic")) {

	sscanf(line,"%*s %d",&(def->aromatic));

      } else if (!strcmp(word,"planar")) {

	sscanf(line,"%*s %d",&(def->planar));

      } else if (!strcmp(word,"atom_node_def")) {

	alloc_atom_node_defs(def);

	read_atom_node_def(paramfile,def->defs + def->n_defs,scheme);

	def->n_defs++;
      }

    } else {

      strcpy(word,"");
    }
  }
}



static void finalise_scheme(ATOM_TYPING_SCHEME *scheme) {

  int i;
  ELEMENT *element;
  ATOM_NODE *node;
  UNITED_ATOM *uatom;
  ATOM_TYPE *type;
  ATOM_TYPE_DEF *def;

  if (scheme == NULL) {

    error_fn("init_atom_typing_scheme: atom typing scheme undefined");
  }

  for (i=0,node=scheme->atom_nodes;i<scheme->n_atom_nodes;i++,node++) {

    node->element = get_element(node->element_name,scheme);

    if (node->element == NULL) {

      error_fn("init_atom_typing_scheme: element '%s' undefined",node->element_name);
    }      
  }

  for (i=0,uatom=scheme->united_atoms;i<scheme->n_united_atoms;i++,uatom++) {

    uatom->atom_node = get_atom_node(uatom->atom_node_name,scheme);

    if (uatom->atom_node == NULL) {

      error_fn("init_atom_typing_scheme: atom_node '%s' undefined",uatom->atom_node_name);
    }      
  }
  
  for (i=0,type=scheme->atom_types;i<scheme->n_atom_types;i++,type++) {

    type->united_atom = get_united_atom(type->united_atom_name,scheme);

    if (type->united_atom == NULL) {

      error_fn("init_atom_typing_scheme: uatom '%s' undefined",type->united_atom_name);
    }

    node = type->united_atom->atom_node;

    if (node) {

      element = node->element;

      if ((element) && (element->flags & METAL_ELEMENT)) {

	type->flags |= METAL_ATOM_TYPE;
      }
    }

    if ((type->flags & HBOND_ACCEPTOR_ATOM_TYPE) || (!strcmp(type->name,"-SH")) || (!strcmp(type->name,">S"))) {

      type->flags |= METAL_ACCEPTOR_ATOM_TYPE;
    }

    type->geometry = get_atom_geometry(type->atom_geometry_name,scheme);

    type->hbond_geometry = get_hbond_geometry(type->hbond_geometry_name);
  }

  for (i=0,element=scheme->elements;i<scheme->n_elements;i++,element++) {

    alloc_atom_types(scheme);

    type = scheme->atom_types + scheme->n_atom_types;

    init_atom_type(type);

    strcpy(type->name,element->name);

    type->element = element;

    if (element->flags & METAL_ELEMENT) {

      type->hbond_geometry = get_hbond_geometry("POINT");

      type->flags |= METAL_ATOM_TYPE;
    }

    type->geometry = get_atom_geometry("X",scheme);

    scheme->n_atom_types++;

    alloc_atom_type_defs(scheme);

    def = scheme->atom_type_defs + scheme->n_atom_type_defs;

    init_atom_type_def(def);

    def->type = type;

    alloc_atom_node_defs(def);

    init_atom_type_def(def->defs);

    strcpy(def->defs->element_name,element->name);

    def->defs->element = element;

    scheme->n_atom_type_defs++;
  }

  for (i=0,type=scheme->atom_types;i<scheme->n_atom_types;i++,type++) {

    type->i = i;
  }
}



ATOM_TYPING_SCHEME* alloc_atom_typing_scheme(void) {

  ATOM_TYPING_SCHEME *scheme;

  scheme = (ATOM_TYPING_SCHEME*) malloc(sizeof(ATOM_TYPING_SCHEME));

  if (scheme == NULL) {

    error_fn("alloc_atom_typing_scheme: out of memory allocating scheme");
  }

  init_scheme(scheme);

  return(scheme);
}



static void init_scheme(ATOM_TYPING_SCHEME *scheme) {

  scheme->n_alloc_elements = 0;
  scheme->n_elements = 0;
  scheme->elements = NULL;
  scheme->n_alloc_atom_nodes = 0;
  scheme->n_atom_nodes = 0;
  scheme->atom_nodes = NULL;
  scheme->n_alloc_united_atoms = 0;
  scheme->n_united_atoms = 0;
  scheme->united_atoms = NULL;
  scheme->n_alloc_atom_geometries = 0;
  scheme->n_atom_geometries = 0;
  scheme->atom_geometries = NULL;
  scheme->n_alloc_atom_types = 0;
  scheme->n_atom_types = 0;
  scheme->atom_types = NULL;
  scheme->n_alloc_residue_atoms = 0;
  scheme->n_residue_atoms = 0;
  scheme->residue_atoms = NULL;
  scheme->n_alloc_atom_type_defs = 0;
  scheme->n_atom_type_defs = 0;
  scheme->atom_type_defs = NULL;
}



static void init_element(ELEMENT *element) {

  strcpy(element->name,"");

  element->id = 0;
  element->cov_radius = 0.70;
  element->vdw_radius = 1.60;
}



static void init_atom_node(ATOM_NODE *atom_node) {

  strcpy(atom_node->name,"");
  strcpy(atom_node->element_name,"");

  atom_node->element = NULL;

  atom_node->n_single = 0;
  atom_node->n_double = 0;
  atom_node->n_triple = 0;

  atom_node->hybridisation = -1;
  atom_node->formal_charge = 0;
}



static void init_united_atom(UNITED_ATOM *uatom) {

  strcpy(uatom->name,"");
  strcpy(uatom->atom_node_name,"");

  uatom->atom_node = NULL;

  uatom->id = 0;
  uatom->hybridisation = 0;
  uatom->n_hydrogens = 0;
  uatom->formal_charge = 0;
  uatom->vdw_radius = 1.60;
}



static void init_atom_type(ATOM_TYPE *type) {

  type->i = -1;
  type->id = -1;

  strcpy(type->name,"");
  strcpy(type->united_atom_name,"");
  strcpy(type->atom_geometry_name,"");

  type->united_atom = NULL;
  type->element = NULL;
  type->geometry = NULL;
  type->hbond_geometry = NULL;

  type->n_hs = 0;
  type->n_lps = 0;

  type->n_alt_types = 0;

  type->flags = 0;
}



static void init_residue_atom(RESIDUE_ATOM *resatom) {

  strcpy(resatom->subname,"");
  strcpy(resatom->name,"");
  strcpy(resatom->atom_type_name,"");

  resatom->id = 0;
}



static void init_atom_type_def(ATOM_TYPE_DEF *def) {

  strcpy(def->node_name,"");
  strcpy(def->element_name,"");

  def->type = NULL;
  def->node = NULL;
  def->element = NULL;

  def->hybridisation = -1;
  def->aromatic = -1;
  def->planar = -1;

  def->n_alloc_defs = 0;
  def->n_defs = 0;
  def->defs = NULL;
}



ATOM_TYPING_SCHEME* get_atom_typing_scheme(void) {

  if (SCHEME == NULL) {

    error_fn("get_atom_typing_scheme: SCHEME undefined");
  }

  return(SCHEME);
}



ELEMENT* get_element(char *name,ATOM_TYPING_SCHEME *scheme) {

  int i;
  ELEMENT *element;

  if (scheme == NULL) {

    scheme = get_atom_typing_scheme();
  }

  for (i=0,element=scheme->elements;i<scheme->n_elements;i++,element++) {

    if (!strcmp(name,element->name)) {

      return(element);
    }
  }

  return(NULL);
}



ELEMENT* get_element_by_id(int id,ATOM_TYPING_SCHEME *scheme) {

  int i;
  ELEMENT *element;

  if (scheme == NULL) {

    scheme = get_atom_typing_scheme();
  }

  for (i=0,element=scheme->elements;i<scheme->n_elements;i++,element++) {

    if (element->id == id) {

      return(element);
    }
  }

  return(NULL);
}



static void alloc_elements(ATOM_TYPING_SCHEME *scheme) {

  if (scheme->n_alloc_elements == 0) {

    scheme->n_alloc_elements = 200;

    scheme->elements = (ELEMENT*) calloc(scheme->n_alloc_elements,sizeof(ELEMENT));

    if (scheme->elements == NULL) {

      error_fn("alloc_elements: out of memory allocating elements");
    }

  } else if (scheme->n_elements == scheme->n_alloc_elements) {

    scheme->n_alloc_elements += 200;

    scheme->elements = (ELEMENT*) realloc(scheme->elements,scheme->n_alloc_elements*sizeof(ELEMENT));

    if (scheme->elements == NULL) {

      error_fn("alloc_elements: out of memory re-allocating elements");
    }
  }
}



static void alloc_atom_nodes(ATOM_TYPING_SCHEME *scheme) {

  if (scheme->n_alloc_atom_nodes == 0) {

    scheme->n_alloc_atom_nodes = 200;

    scheme->atom_nodes = (ATOM_NODE*) calloc(scheme->n_alloc_atom_nodes,sizeof(ATOM_NODE));

    if (scheme->atom_nodes == NULL) {

      error_fn("alloc_atom_nodes: out of memory allocating atom_nodes");
    }

  } else if (scheme->n_atom_nodes == scheme->n_alloc_atom_nodes) {

    scheme->n_alloc_atom_nodes += 200;

    scheme->atom_nodes = (ATOM_NODE*) realloc(scheme->atom_nodes,scheme->n_alloc_atom_nodes*sizeof(ATOM_NODE));

    if (scheme->atom_nodes == NULL) {

      error_fn("alloc_atom_nodes: out of memory re-allocating atom_nodes");
    }
  }
}



static void alloc_united_atoms(ATOM_TYPING_SCHEME *scheme) {

  if (scheme->n_alloc_united_atoms == 0) {

    scheme->n_alloc_united_atoms = 200;

    scheme->united_atoms = (UNITED_ATOM*) calloc(scheme->n_alloc_united_atoms,sizeof(UNITED_ATOM));

    if (scheme->united_atoms == NULL) {

      error_fn("alloc_united_atoms: out of memory allocating united_atoms");
    }

  } else if (scheme->n_united_atoms == scheme->n_alloc_united_atoms) {

    scheme->n_alloc_united_atoms += 200;

    scheme->united_atoms = (UNITED_ATOM*) realloc(scheme->united_atoms,scheme->n_alloc_united_atoms*sizeof(UNITED_ATOM));

    if (scheme->united_atoms == NULL) {

      error_fn("alloc_united_atoms: out of memory re-allocating united_atoms");
    }
  }
}



static void alloc_atom_types(ATOM_TYPING_SCHEME *scheme) {

  if (scheme->n_alloc_atom_types == 0) {

    scheme->n_alloc_atom_types = 200;

    scheme->atom_types = (ATOM_TYPE*) calloc(scheme->n_alloc_atom_types,sizeof(ATOM_TYPE));

    if (scheme->atom_types == NULL) {

      error_fn("alloc_atom_types: out of memory allocating atom_types");
    }

  } else if (scheme->n_atom_types == scheme->n_alloc_atom_types) {

    scheme->n_alloc_atom_types += 200;

    scheme->atom_types = (ATOM_TYPE*) realloc(scheme->atom_types,scheme->n_alloc_atom_types*sizeof(ATOM_TYPE));

    if (scheme->atom_types == NULL) {

      error_fn("alloc_atom_types: out of memory re-allocating atom_types");
    }
  }
}



static void alloc_residue_atoms(ATOM_TYPING_SCHEME *scheme) {

  if (scheme->n_alloc_residue_atoms == 0) {

    scheme->n_alloc_residue_atoms = 200;

    scheme->residue_atoms = (RESIDUE_ATOM*) calloc(scheme->n_alloc_residue_atoms,sizeof(RESIDUE_ATOM));

    if (scheme->residue_atoms == NULL) {

      error_fn("alloc_residue_atoms: out of memory allocating residue_atoms");
    }

  } else if (scheme->n_residue_atoms == scheme->n_alloc_residue_atoms) {

    scheme->n_alloc_residue_atoms += 200;

    scheme->residue_atoms = (RESIDUE_ATOM*) realloc(scheme->residue_atoms,scheme->n_alloc_residue_atoms*sizeof(RESIDUE_ATOM));

    if (scheme->residue_atoms == NULL) {

      error_fn("alloc_residue_atoms: out of memory re-allocating residue_atoms");
    }
  }
}



static void alloc_atom_type_defs(ATOM_TYPING_SCHEME *scheme) {

  if (scheme->n_alloc_atom_type_defs == 0) {

    scheme->n_alloc_atom_type_defs = 200;

    scheme->atom_type_defs = (ATOM_TYPE_DEF*) calloc(scheme->n_alloc_atom_type_defs,sizeof(ATOM_TYPE_DEF));

    if (scheme->atom_type_defs == NULL) {

      error_fn("alloc_atom_type_defs: out of memory allocating atom_type_defs");
    }

  } else if (scheme->n_atom_type_defs == scheme->n_alloc_atom_type_defs) {

    scheme->n_alloc_atom_type_defs += 200;

    scheme->atom_type_defs = (ATOM_TYPE_DEF*) realloc(scheme->atom_type_defs,scheme->n_alloc_atom_type_defs*sizeof(ATOM_TYPE_DEF));

    if (scheme->atom_type_defs == NULL) {

      error_fn("alloc_atom_type_defs: out of memory re-allocating atom_type_defs");
    }
  }
}



static void alloc_atom_node_defs(ATOM_TYPE_DEF* def) {

  if (def->n_alloc_defs == 0) {

    def->n_alloc_defs = 5;

    def->defs = (ATOM_TYPE_DEF*) calloc(def->n_alloc_defs,sizeof(ATOM_TYPE_DEF));

    if (def->defs == NULL) {

      error_fn("alloc_atom_node_defs: out of memory allocating defs");
    }

  } else if (def->n_defs == def->n_alloc_defs) {

    def->n_alloc_defs += 5;

    def->defs = (ATOM_TYPE_DEF*) realloc(def->defs,(def->n_alloc_defs)*sizeof(ATOM_TYPE_DEF));

    if (def->defs == NULL) {

      error_fn("alloc_atom_node_defs: out of memory re-allocating defs");
    }
  }
}



ATOM_NODE* get_atom_node(char *name,ATOM_TYPING_SCHEME *scheme) {

  int i;
  ATOM_NODE *node;

  for (i=0,node=scheme->atom_nodes;i<scheme->n_atom_nodes;i++,node++) {

    if (!strcmp(node->name,name)) {

      return(node);
    }
  }

  return(NULL);
}



UNITED_ATOM* get_united_atom(char *name,ATOM_TYPING_SCHEME *scheme) {

  int i;
  UNITED_ATOM *uatom;

  for (i=0,uatom=scheme->united_atoms;i<scheme->n_united_atoms;i++,uatom++) {

    if (!strcmp(uatom->name,name)) {

      return(uatom);
    }
  }

  return(NULL);
}



void set_molecule_atom_types(MOLECULE *molecule,ATOM_TYPING_SCHEME *scheme) {

  int i;
  ATOM *atom;

  for (i=0,atom=molecule->atom;i<molecule->natoms;i++,atom++) {

    atom_type_atom(atom,scheme);
  }
}



void atom_type_atom(ATOM *atom,ATOM_TYPING_SCHEME *scheme) {

  int n_bonds_atom,n_bonds_node;
  UNITED_ATOM *uatom;
  ATOM_NODE *node;

  if ((atom->flags & ATOM_TYPE_ATTEMPTED) || (atom->type != NULL)) {

    return;
  }

  atom->flags |= ATOM_TYPE_ATTEMPTED;

  atom->set = get_atom_set(atom);

  atom->type = resatoms2type(atom,scheme);

  if (atom->type) {

    if (atom->node) {

      uatom = atom->type->united_atom;

      if (uatom) {

	node = uatom->atom_node;

	n_bonds_atom = atom->node->n_single + atom->node->n_double + atom->node->n_triple;
	n_bonds_node = node->n_single + node->n_double + node->n_triple;

	if (n_bonds_atom == n_bonds_node) {

	  check_atom_type(atom,scheme);

	  if (atom->type->flags & AROMATIC_ATOM_TYPE) {

	    atom->flags |= AROMATIC_ATOM;
	  }

	  return;
	}
      }
    }

    atom->type = NULL;
  }
  
  atom->type = tautomer_atom_type(atom,scheme);

  if (atom->type != NULL) {

    check_atom_type(atom,scheme);

    if (atom->type->flags & AROMATIC_ATOM_TYPE) {

      atom->flags |= AROMATIC_ATOM;
    }

    return;
  }

  atom->type = derive_atom_type(atom,scheme);
 
  check_atom_type(atom,scheme);

  if ((atom->type != NULL) && (atom->type->flags & AROMATIC_ATOM_TYPE)) {

    atom->flags |= AROMATIC_ATOM;
  }
}



void atom_type_heme_nitrogens(SYSTEM *system) {

  int i,j,k,l,n;
  ATOM *atom,*atom1,*atom2;
  ATOM_TYPE *type,*heme_nitrogen;
  ATOMLIST *nlist;
  CONTACTLIST *contactlist;
  CONTACT *contact1,*contact2;
  MOLECULE **moleculep,*molecule;
  MOLECULE_LIST *list;

  list = system->molecule_list;

  if (list == NULL) {

    return;
  }

  nlist = alloc_atomlist();

  init_atomlist(nlist);

  heme_nitrogen = get_atom_type("Aromatic =N-",system->settings->atom_typing_scheme);

  for (i=0,moleculep=list->molecules;i<list->n_molecules;i++,moleculep++) {

    molecule = *moleculep;

    for (j=0,atom=molecule->atom;j<molecule->natoms;j++,atom++) {

      type = atom->type;

      if ((type) && (type->flags & METAL_ATOM_TYPE)) {

	set_contacts_atom(atom,system,1);

	contactlist = atom->contactlist;

	if ((contactlist) && (contactlist->ncontacts > 0)) {

	  nlist->natoms = 0;

	  for (k=0,contact1=contactlist->contacts;k<contactlist->ncontacts;k++,contact1++) {

	    atom1 = contact1->atom2;

	    if ((atom1->element->id == NITROGEN) && (hbond_geometry_score(contact1) > 0.2)) {

	      n = 1;

	      for (l=0,contact2=contactlist->contacts;l<contactlist->ncontacts;l++,contact2++) {

		if (contact1 != contact2) {

		  atom2 = contact2->atom2;

		  if ((atom2->element->id == NITROGEN) && (hbond_geometry_score(contact2) > 0.2)) {

		    if (same_residue_atoms(atom1,atom2)) {

		      n++;
		    }
		  }
		}
	      }

	      if (n == 4) {

		add_atom_to_list(nlist,atom1,0);
	      }
	    }
	  }

	  if (nlist->natoms == 4) {

	    for (k=0;k<4;k++) {

	      (*(nlist->atom+k))->type = heme_nitrogen;
	    }
	  }
	}
      }
    }
  }

  free_atomlist(nlist);
}



static int get_atom_set(ATOM *atom) {

  MOLECULE *molecule;

  molecule = atom->molecule;

  if (molecule->flags & LIGAND_MOLECULE) {

    return(LIGAND_ATOM_SET);

  } else if ((molecule->flags & PROTEIN_MOLECULE) && (atom->flags & AMINO_ACID_ATOM)) {

    return(PROTEIN_ATOM_SET);

  } else if (atom->flags & WATER_OXYGEN) {

    return(WATER_ATOM_SET);
  }

  return(OTHER_ATOM_SET);
}



static void check_atom_type(ATOM *atom,ATOM_TYPING_SCHEME *scheme) {

  ATOM_TYPE *type;
  UNITED_ATOM *uatom;
  ATOM_NODE *node;
  ELEMENT *element;
  MOLECULE *molecule;

  molecule = atom->molecule;

  if ((molecule) && (molecule->use_hydrogens)) {

    check_atom_type_vs_hydrogens(atom,scheme);
  }

  type = atom->type;

  if (type == NULL) {

    return;
  }

  uatom = type->united_atom;

  if (uatom) {

    node = uatom->atom_node;

    element = (node) ? node->element : NULL;

  } else {

    node = NULL;

    element = atom->type->element;
  }

  if (element == NULL) {

    error_fn("check_atom_type: element undefined for atom type %d",type->id);
  }

  if (atom->element == NULL) {

    warning_fn("check_atom_type: setting element from atom type for atom %d",atom->id);

    atom->element = element;

  } else if (atom->element != element) {

    warning_fn("check_atom_type: atom element does not match atom type element for atom %d (using atom type element)",atom->id);

    atom->element = element;
  }

  if (node) {

    atom->node = node;

    //printf("%-10s",atom->node->name);
    //write_atom(PLI_STDOUT,atom,TYPE_ASTYLE,0,0);
  }

  if ((atom->type != NULL) && (atom->flags & AMINO_ACID_ATOM) && (!(atom->type->flags & AMINO_ACID_ATOM_TYPE))) {

    atom->error_flags |= ATOM_TYPING_ERROR;
  }
}



static void check_atom_type_vs_hydrogens(ATOM *atom,ATOM_TYPING_SCHEME *scheme) {

  int i,j,use_type,taro,aaro,n_tbonds,n_abonds;
  double max_p;
  ATOM_NODE *node,*anode;
  UNITED_ATOM *uatom;
  ATOM_TYPE *type,*best_type;
  ALT_ATOM_TYPE *alt_type,*best_alt_type;
  RESIDUE_ATOM *resatom;

  type = atom->type;

  if ((type == NULL) || (type->n_alt_types == 0)) {

    return;
  }

  max_p = 0.0;
  best_alt_type = NULL;

  for (i=0,alt_type=type->alt_types;i<type->n_alt_types;i++,alt_type++) {

    if (alt_type->group_status == atom->group_status) {

      uatom = alt_type->type->united_atom;

      if ((uatom) && (uatom->n_hydrogens == atom->n_hydrogens)) {

	if (alt_type->probability > max_p) {
	  
	  best_alt_type = alt_type;
	}
      }
    }
  }

  if (best_alt_type) {

    update_group_statuses(atom,best_alt_type,scheme);

    atom->type = best_alt_type->type;
    atom->type_probability = best_alt_type->probability;

  } else {

    best_type = NULL;

    for (i=0,type=scheme->atom_types;i<scheme->n_atom_types;i++,type++) {

      // check against atom nodes, aromaticity, etc.
      // AND against n_hydrogens
      // dont use residue atom types

      uatom = type->united_atom;

      if (uatom) {

	node = uatom->atom_node;
	anode = atom->node;

	if ((node) && (anode)) {

	  taro = (type->flags & AROMATIC_ATOM_TYPE) ? 1 : 0;
	  aaro = (atom->type->flags & AROMATIC_ATOM_TYPE) ? 1 : 0;
	  
	  n_tbonds = node->n_single + node->n_double + node->n_triple;
	  n_abonds = anode->n_single + anode->n_double + anode->n_triple;
	  
	  if ((n_tbonds == n_abonds) && (node->element == anode->element) &&
	      (uatom->n_hydrogens == atom->n_hydrogens) && 
	      (taro == aaro)) {
	    
	    use_type = 1;
	    
	    for (j=0,resatom=scheme->residue_atoms;j<scheme->n_residue_atoms;j++,resatom++) {
	      
	      if (!strcmp(resatom->atom_type_name,type->name)) {
		
		use_type = 0;
	      }
	    }

	    if ((use_type) && (!best_type)) {

	      best_type = type;
	    }
	  }
	}
      }
    }

    if (best_type) {

      atom->type = best_type;
    
    } else {

      best_alt_type = NULL;
      max_p = 0.0;
      type = atom->type;

      for (i=0,alt_type=type->alt_types;i<type->n_alt_types;i++,alt_type++) {

	if (alt_type->group_status == atom->group_status) {

	  if (alt_type->probability > max_p) {
	  
	    best_alt_type = alt_type;
	  }
	}
      }

      if (best_alt_type) {

	atom->type = best_alt_type->type;
        atom->type_probability = best_alt_type->probability;

      } else {

	warning_fn("atom type could not be resolved from connected hydrogens for atom %d (%s,%s,%d,%s)\n",atom->id,atom->name,atom->subname,atom->subid,atom->type->name);
      }
    }
  }
}


static ATOM_TYPE* tautomer_atom_type(ATOM *atom,ATOM_TYPING_SCHEME *scheme) {

  ATOM_TYPE *type;

  if (atom->node == NULL) {

    return(NULL);
  }

  if (charged_aromatic_n(atom)) {

    return(get_atom_type("Aromatic =NH-",scheme));
  }

  if (type = ambiguous_aromatic_n(atom,scheme)) {

    return(type);
  }

  if (ambiguous_pyridone_n(atom)) {

    return(get_atom_type("Pyridone N?",scheme));
  }

  if (ambiguous_pyridone_o(atom)) {

    return(get_atom_type("Pyridone O?",scheme));
  }

  if (ambiguous_pyridine_n(atom)) {

    return(get_atom_type("Pyridine N?",scheme));
  }

  return(NULL);
}



static ATOM_TYPE* derive_atom_type(ATOM *atom,ATOM_TYPING_SCHEME *scheme) {

  int i,j;
  ATOM_TYPE_DEF *def;

  for (i=0,def=scheme->atom_type_defs;i<scheme->n_atom_type_defs;i++,def++) {

    if (atom_type_match(atom,def->defs)) {

      return(def->type);
    }
  }

  return(NULL);
}



static int charged_aromatic_n(ATOM *atom) {

  RING *ring;

  if (atom->node == NULL) {

    return(0);
  }

  if (strcmp(atom->node->name,"=N-")) {

    return(0);
  }

  find_ring(atom,6);

  ring = atom->ring;

  if (!pyridine_ring(ring)) {

    return(0);
  }

  if (amino_pyridine_atom_pair(atom,*(ring->atom+3))) {

    return(1);
  }

  return(0);
}



static ATOM_TYPE* ambiguous_aromatic_n(ATOM *atom,ATOM_TYPING_SCHEME *scheme) {

  int i,n_carbon,n_nitrogen,posN[5],diff;
  ATOM **ratom;
  ATOM_NODE *rnode;
  RING *ring;

  if (atom->node == NULL) {

    return(NULL);
  }

  if ((strcmp(atom->node->name,"=N-")) && (strcmp(atom->node->name,">N"))) {

    return(NULL);
  }

  find_ring(atom,5);

  ring = atom->ring;

  if ((!ring) || (ring->aromatic == 0) || (ring->size != 5)) {

    return(NULL);
  }

  n_carbon = n_nitrogen = 0;

  for (i=0,ratom=ring->atom;i<ring->size;i++,ratom++) {

    if ((rnode = (*ratom)->node) == NULL) {

      return(NULL);
    }

    if ((!strcmp(rnode->name,"=C-")) || (!strcmp(rnode->name,"=C<"))) {

      n_carbon++;

    } else if ((!strcmp(rnode->name,"=N-")) || (!strcmp(rnode->name,">N"))) {

      posN[n_nitrogen] = i;

      n_nitrogen++;
    }
  }

  if ((n_carbon == 3) && (n_nitrogen == 2)) {

    diff = posN[1] - posN[0];

    if ((diff == 1) || (diff == 4)) {

      return(get_atom_type("Pyrazole N?",scheme));

    } else {

      return(get_atom_type("Imidazole N?",scheme));
    }

  } else if ((n_carbon == 2) && (n_nitrogen == 3)) {

    return(get_atom_type("Triazole N?",scheme));
  }

  return(NULL);
}



static int ambiguous_pyridone_n(ATOM *atom) {

  ATOM *oxygen,*nitrogen;

  if (atom->node == NULL) {

    return(0);
  }

  if ((strcmp(atom->node->name,"=N-")) && (strcmp(atom->node->name,">N"))) {

    return(0);
  }

  oxygen = pyridone_o_from_n(atom);

  if (oxygen) {

    nitrogen = pyridone_n_from_o(oxygen);

    if (nitrogen == atom) {

      return(1);
    }
  }

  return(0);
}




static int ambiguous_pyridone_o(ATOM *atom) {

  ATOM *nitrogen,*oxygen;

  if (atom->node == NULL) {

    return(0);
  }

  if ((strcmp(atom->node->name,"=O")) && (strcmp(atom->node->name,"-O"))) {

    return(0);
  }

  nitrogen = pyridone_n_from_o(atom);

  if (nitrogen) {

    oxygen = pyridone_o_from_n(nitrogen);

    if (oxygen == atom) {

      return(1);
    }
  }

  return(0);
}



ATOM* pyridone_o_from_n(ATOM *nitrogen) {

  int i,n_pairs;
  ATOM *carbon,**atomp,*atom;
  RING *ring;

  if (!nitrogen->ring) {

    find_ring(nitrogen,6);
  }

  ring = nitrogen->ring;

  if (!aromatic_nc_6_ring(ring)) {

    return(0);
  }

  n_pairs = 0;

  carbon = NULL;

  if (pyridone_atom_pair(nitrogen,*(ring->atom+1))) {

    n_pairs++;

    carbon = *(ring->atom+1);
  }

  if (pyridone_atom_pair(nitrogen,*(ring->atom+3))) {

    n_pairs++;

    carbon = *(ring->atom+3);
  }

  if (pyridone_atom_pair(nitrogen,*(ring->atom+5))) {

    n_pairs++;

    carbon = *(ring->atom+5);
  }

  if ((n_pairs != 1) || (!carbon)) {

    return(NULL);
  }

  for (i=0,atomp=carbon->connections->atom;i<carbon->connections->natoms;i++,atomp++) {

    atom = *atomp;

    if (atom->node == NULL) {

      return(NULL);
    }

    if (((!strcmp(atom->node->name,"=O")) && (!strcmp(nitrogen->node->name,">N"))) ||
	((!strcmp(atom->node->name,"-O")) && (!strcmp(nitrogen->node->name,"=N-")))) {

      return(atom);
    }
  }
 
  return(NULL);
}



ATOM* pyridone_n_from_o(ATOM *oxygen) {

  int i,n_pairs;
  ATOM *carbon,*nitrogen,*oxygen1;
  RING *ring;

  if (!(oxygen->connections) || (oxygen->connections->natoms != 1)) {

    return(NULL);
  }

  carbon = *(oxygen->connections->atom);

  if (!carbon) {

    return(NULL);
  }

  if (!carbon->ring) {

    find_ring(carbon,6);
  }

  ring = carbon->ring;

  if (!aromatic_nc_6_ring(ring)) {

    return(0);
  }

  n_pairs = 0;

  nitrogen = NULL;

  if (pyridone_atom_pair(*(ring->atom+1),carbon)) {

    n_pairs++;

    nitrogen = *(ring->atom+1);
  }

  if (pyridone_atom_pair(*(ring->atom+3),carbon)) {

    n_pairs++;

    nitrogen = *(ring->atom+3);
  }

  if (pyridone_atom_pair(*(ring->atom+5),carbon)) {

    n_pairs++;

    nitrogen = *(ring->atom+5);
  }

  if ((n_pairs == 1) && (nitrogen)) {

    oxygen1 = pyridone_o_from_n(nitrogen);

    if (oxygen1 == oxygen) {

      return(nitrogen);
    }
  }

  return(NULL);
}



ATOM *carboxyl_o_from_o(ATOM *oxygen1) {

  int i;
  ATOM *carbon,**atomp,*atom,*oxygen2;

  if (!(oxygen1->connections) || (oxygen1->connections->natoms != 1)) {

    return(NULL);
  }

  carbon = *(oxygen1->connections->atom);

  if ((!carbon) || (!carbon->connections)) {

    return(NULL);
  }

  oxygen2 = NULL;

  for (i=0,atomp=carbon->connections->atom;i<carbon->connections->natoms;i++,atomp++) {

    atom = *atomp;

    if ((atom != oxygen1) && (atom->element) && (atom->element->id == OXYGEN) &&
	(atom->connections) && (atom->connections->natoms == 1)) {

      if (oxygen2) {

	return(NULL);
      }

      oxygen2 = atom;
    }
  }

  return(oxygen2);
}



static int ambiguous_pyridine_n(ATOM *atom) {

  RING *ring;

  if (atom->node == NULL) {

    return(0);
  }

  if ((strcmp(atom->node->name,"=N-")) && (strcmp(atom->node->name,">N"))) {

    return(0);
  }

  find_ring(atom,6);

  ring = atom->ring;

  if (!pyridine_ring(ring)) {

    return(0);
  }

  if (amino_pyridine_atom_pair(atom,*(ring->atom+1))) {

    return(1);
  }

  if (amino_pyridine_atom_pair(atom,*(ring->atom+5))) {

    return(1);
  }

  return(0);
}



static int pyridine_ring(RING *ring) {

  int i,n_carbon,n_nitrogen;
  ATOM **ratom;
  ATOM_NODE *rnode;

  if (ring == NULL) {

    return(0);
  }

  if ((ring->aromatic == 0) || (ring->size != 6)) {

    return(0);
  }

  n_carbon = n_nitrogen = 0;

  for (i=0,ratom=ring->atom;i<ring->size;i++,ratom++) {

    if ((rnode = (*ratom)->node) == NULL) {

      return(0);
    }

    if ((!strcmp(rnode->name,"=C-")) || (!strcmp(rnode->name,"=C<"))) {

      n_carbon++;

    } else if ((!strcmp(rnode->name,"=N-")) || (!strcmp(rnode->name,">N"))) {

      n_nitrogen++;
    }
  }

  if ((n_carbon == 5) && (n_nitrogen == 1)) {

    return(1);
  }

  return(0);
}



static int aromatic_nc_6_ring(RING *ring) {

  int i,n_carbon,n_nitrogen;
  ATOM **ratom;
  ATOM_NODE *rnode;

  if (ring == NULL) {

    return(0);
  }

  if ((ring->aromatic == 0) || (ring->size != 6)) {

    return(0);
  }

  n_carbon = n_nitrogen = 0;

  for (i=0,ratom=ring->atom;i<ring->size;i++,ratom++) {

    if ((rnode = (*ratom)->node) == NULL) {

      return(0);
    }

    if ((!strcmp(rnode->name,"=C-")) || (!strcmp(rnode->name,"=C<"))) {

      n_carbon++;

    } else if ((!strcmp(rnode->name,"=N-")) || (!strcmp(rnode->name,">N"))) {

      n_nitrogen++;
    }
  }

  if (n_carbon + n_nitrogen == 6) {

    return(1);
  }

  return(0);
}




static int pyridone_atom_pair(ATOM *n,ATOM *c) {

  int i;
  ATOM **atom;

  if ((n->node == NULL) || (c->node == NULL)) {

    return(0);
  }

  if ((n->element->id != NITROGEN) || (c->element->id != CARBON)) {

    return(0);
  }

  if (c->connections->natoms != 3) {

    return(0);
  }

  for (i=0,atom=c->connections->atom;i<c->connections->natoms;i++,atom++) {

    if ((*atom)->node == NULL) {

      return(0);
    }

    if (((!strcmp((*atom)->node->name,"=O")) && (!strcmp(n->node->name,">N"))) ||
	((!strcmp((*atom)->node->name,"-O")) && (!strcmp(n->node->name,"=N-")))) {

      return(1);
    }
  }

  return(0);
}



static int amino_pyridine_atom_pair(ATOM *n,ATOM *c) {

  int i;
  ATOM **atom1,*atom2;
  ATOMLIST *conns;

  if ((n->node == NULL) || (c->node == NULL)) {

    return(0);
  }

  if ((strcmp(n->node->name,"=N-")) || (strcmp(c->node->name,"=C<"))) {

    return(0);
  }

  if (c->connections->natoms != 3) {

    return(0);
  }

  for (i=0,atom1=c->connections->atom;i<c->connections->natoms;i++,atom1++) {

    if ((*atom1)->node == NULL) {

      return(0);
    }

    if (!strcmp((*atom1)->node->name,"-N")) {

      return(1);
    }

    if (!strcmp((*atom1)->node->name,">N")) {

      conns = (*atom1)->connections;

      if (conns->natoms == 2) {

	atom2 = (conns->atom[0] == c) ? conns->atom[1] : conns->atom[0];

	if (atom2->node == NULL) {

	  return(0);
	}

	if ((atom2->node->element->id == CARBON) && (atom2->node->hybridisation == 3)) {

	  return(1);
	}
      }
    }
  }

  return(0);
}



static int atom_type_match(ATOM *atom,ATOM_TYPE_DEF *def) {

  int i,j,*atom_matched,def_matched;
  ATOM **batom;
  ATOM_TYPE_DEF *bdef;

  if (strcmp(def->element_name,"")) {

    if (strcmp(def->element_name,atom->element->name)) {

      return(0);
    }
  }

  if (strcmp(def->node_name,"")) {

    if ((atom->node == NULL) || (strcmp(def->node_name,atom->node->name))) {

      return(0);
    }
  }

  if (def->hybridisation != -1) {

    if ((atom->node == NULL) || (def->hybridisation != atom->node->hybridisation)) {

      return(0);
    }
  }

  if (def->planar != -1) {

    if (def->planar != planar_atom(atom)) {

      return(0);
    }
  }

  if (def->aromatic != -1) {
    
    find_ring(atom,6);

    if (((def->aromatic == 1) && (!(atom->flags & AROMATIC_ATOM))) ||
	((def->aromatic == 0) && (atom->flags & AROMATIC_ATOM))) {

      return(0);
    }
  }

  if (def->n_defs > 0) {

    if (atom->connections == NULL) {

      return(0);
    }

    if (atom->connections->natoms == 0) {

      return(0);
    }

    atom_matched = (int*) calloc(atom->connections->natoms,sizeof(int));

    for (i=0;i<atom->connections->natoms;i++) {

      atom_matched[i] = 0;
    }

    for (i=0,bdef=def->defs;i<def->n_defs;i++,bdef++) {

      def_matched = 0;

      for (j=0,batom=atom->connections->atom;j<atom->connections->natoms;j++,batom++) {

	if (!atom_matched[j]) {

	  if (atom_type_match(*batom,bdef)) {

	    atom_matched[j] = 1;

	    def_matched = 1;

	    break;
	  }
	}
      }

      if (!def_matched) {

	return(0);
      }
    }

    free(atom_matched);
  }

  return(1);
}



void set_system_atom_nodes(SYSTEM *system) {

  int i,j;
  ATOM *atom;
  MOLECULE **moleculep,*molecule;
  MOLECULE_LIST *list;

  list = system->molecule_list;

  if (list == NULL) {

    return;
  }

  for (i=0,moleculep=list->molecules;i<list->n_molecules;i++,moleculep++) {

    molecule = *moleculep;

    for (j=0,atom=molecule->atom;j<molecule->natoms;j++,atom++) {

      set_atom_node(atom,system->settings->atom_typing_scheme);
    }
  }
}



void set_molecule_atom_nodes(MOLECULE *molecule,ATOM_TYPING_SCHEME *scheme) {

  int i;
  ATOM *atom;

  for (i=0,atom=molecule->atom;i<molecule->natoms;i++,atom++) {

    set_atom_node(atom,scheme);
  }
}



ATOM_NODE* set_atom_node(ATOM *atom,ATOM_TYPING_SCHEME *scheme) {

  int i,nbonds[4];
  BOND *bond;
  ELEMENT *element;
  ATOM **batom;
  ATOM_NODE *node;
  MOLECULE *molecule;

  if ((atom->node != NULL) || (atom->flags & NODE_SET_ATTEMPTED)) {

    return(atom->node);
  }

  atom->flags |= NODE_SET_ATTEMPTED;

  if (atom->element->id == HYDROGEN) {

    return(atom->node);
  }

  molecule = atom->molecule;

  nbonds[1] = nbonds[2] = nbonds[3] = 0;

  // find bonds:

  if (atom->connections != NULL) {

    for (i=0,batom=atom->connections->atom;i<atom->connections->natoms;i++,batom++) {

      if ((*batom)->element->id != HYDROGEN) {

	bond = get_bond(molecule,atom,*batom);

	if (bond == NULL) {

	  error_fn("atom_type_match: bond undefined for atoms %s - %s",atom->name,(*batom)->name);
	}
	
	nbonds[bond->type]++;  
      }
    }
  }

  for (i=0,node=scheme->atom_nodes;i<scheme->n_atom_nodes;i++,node++) {

    element = get_element(node->element_name,scheme);

    if ((element == atom->element) &&
	(nbonds[1] == node->n_single) && (nbonds[2] == node->n_double) && (nbonds[3] == node->n_triple)) {

      atom->node = node;

      return(atom->node);
    }
  }

  warning_fn("atom_type_match: node undefined for atom %s (%d) in residue %s (%s)",atom->name,atom->id,atom->subname,atom->molecule->name);


  return(NULL);
}



ATOM_TYPE* resatoms2type(ATOM *atom,ATOM_TYPING_SCHEME *scheme) {

  int i;
  ATOM_TYPE *type;
  RESIDUE_ATOM *resatom;

  for (i=0,resatom=scheme->residue_atoms;i<scheme->n_residue_atoms;i++,resatom++) {

    if ((!strcmp(resatom->subname,atom->subname)) && (!strcmp(resatom->name,atom->name))) {

      type = get_atom_type(resatom->atom_type_name,scheme);

      if (type == NULL) {

	error_fn("resatoms2type: unknown atom type '%s'",resatom->atom_type_name);
      }

      return(type);
    }
  }

  return(NULL);
}



void set_atom_type(ATOM *atom,int id,ATOM_TYPING_SCHEME *scheme) {

  ATOM_TYPE *type;
  UNITED_ATOM *uatom;
  ATOM_NODE *node;

  type = get_atom_type_by_id(id,scheme);

  if (type->element != NULL) {

    atom->element = type->element;

  } else {

    uatom = type->united_atom;
    node = uatom->atom_node;

    atom->element = node->element;
  }

  atom->type = type;
}



ATOM_TYPE* get_atom_type(char *name,ATOM_TYPING_SCHEME *scheme) {

  int i;
  ATOM_TYPE *type;

  for (i=0,type=scheme->atom_types;i<scheme->n_atom_types;i++,type++) {

    if (!strcmp(type->name,name)) {

      return(type);
    }
  }

  return(NULL);
}



ATOM_TYPE* get_atom_type_by_id(int id,ATOM_TYPING_SCHEME *scheme) {

  int i;
  ATOM_TYPE *type;

  for (i=0,type=scheme->atom_types;i<scheme->n_atom_types;i++,type++) {

    if (type->id == id) {

      return(type);
    }
  }

  return(NULL);
}



ALT_ATOM_TYPE *get_alt_type(ATOM_TYPE *type,char *name,int group_status) {

  int i;
  ALT_ATOM_TYPE *alt_type;

  if ((type) && (type->n_alt_types)) {

    for (i=0,alt_type=type->alt_types;i<type->n_alt_types;i++,alt_type++) {

      if ((alt_type->group_status == group_status) && (!strcmp(alt_type->type->name,name))) {

	return(alt_type);
      }
    }
  }

  return(NULL);
}



double get_atom_type_vdw_radius(ATOM_TYPE *type) {

  if (type->id == -1) {

    return(type->element->vdw_radius);
  }

  return(type->united_atom->vdw_radius);
}



void change_atom_type(ATOM *atom,ATOM_TYPE *type) {

  atom->type = type;

  if (type->id == -1) {

    atom->element = type->element;

  } else {

    atom->element = type->united_atom->atom_node->element;
  }

  set_atom_vdw_radii(atom,1.40);

  atom->status = ATOM_CHANGED;
  atom->cstatus = ATOM_NOTHING_CALCULATED;
}
