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

#define PLIFF_SYSTEM_SCORE 0
#define PLIFF_SYSTEM_NB_SCORE 1
#define PLIFF_SYSTEM_CONTACT_SCORE 2
#define PLIFF_SYSTEM_GEOMETRY_SCORE 3
#define PLIFF_SYSTEM_TRIPLET_SCORE 4
#define PLIFF_SYSTEM_RB_SCORE 5
#define PLIFF_SYSTEM_TC_SCORE 6
#define PLIFF_SYSTEM_LIGAND_SCORE 7

#define PLIFF_ATOM_SCORE 0
#define PLIFF_ATOM_NB_SCORE 1
#define PLIFF_ATOM_CONTACT_SCORE 2
#define PLIFF_ATOM_GEOMETRY_SCORE 3
#define PLIFF_ATOM_TRIPLET_SCORE 4
#define PLIFF_ATOM_Z_COORD_SCORE 5

#define PLIFF_CONTACT_SCORE 0
#define PLIFF_R_SCORE 1
#define PLIFF_ALPHA1_SCORE 2
#define PLIFF_BETA1_SCORE 3
#define PLIFF_ALPHA2_SCORE 4
#define PLIFF_BETA2_SCORE 5
#define PLIFF_GEOMETRY_SCORE 6
#define PLIFF_HBOND_GEOMETRY_SCORE 7
#define PLIFF_CONTACT_SCALE 8



static double pliff_zero_coeff = 0.0;
static double pliff_contact_coeff = 0.0034;
static double pliff_geometry_coeff = 0.0034;

static int pliff_recalc_contacts = 1;
static int pliff_use_areas = 1;
static int pliff_ignore_long_contacts = 0;



static void pliff_score_contacts(ATOMLIST*,ATOM_TYPING_SCHEME*,FORCE_FIELD*);
static void pliff_scale_contacts(ATOMLIST*,ATOM_TYPING_SCHEME*,FORCE_FIELD*);
static void pliff_score_atoms(ATOMLIST*);
static void pliff_scale_atom_contacts(ATOM*,ATOM_TYPING_SCHEME*,FORCE_FIELD*);
static int pliff_get_da_hbonds(CONTACTLIST*,int,int,CONTACT**,int*);
static void pliff_permute_da_hbonds(CONTACT**,int,int,int,int*,int,int,double*,double*,double*);
static double pliff_score_da_permutation(CONTACT**,int,int*,double*);
static double pliff_hbonds_triplet_score(CONTACT**,int,int*);
static double pliff_triplet_score(double);
static void pliff_contactlist2hbonds(CONTACTLIST*,CONTACT**,int*,double,double);
static double pliff_score_vs_histogram(int,double,int);
static ATOM_TYPE* pliff_qscore_get_atom_type(int,int,ATOM_TYPING_SCHEME*);
static double pliff_init_contact_scale(CONTACT*,NONBONDED_FF*,HISTOGRAM*,HISTOGRAM*);



void init_pliff_settings(void) {

  pliff_zero_coeff = (params_get_parameter("pliff_zero_coeff"))->value.d;
  pliff_contact_coeff = (params_get_parameter("pliff_contact_coeff"))->value.d;
  pliff_geometry_coeff = (params_get_parameter("pliff_geometry_coeff"))->value.d;

  pliff_ignore_long_contacts = (params_get_parameter("pliff_ignore_long_contacts"))->value.i;
}



void run_pliff_qscore(SETTINGS *settings) {

  int type_id1,type_id2,element_id1,element_id2,n_read,key;
  double area,distance,alpha1,beta1,alpha2,beta2;
  double vdw_radius1,vdw_radius2;
  double contact_score,geometry_score,score;
  char line[MAX_LINE_LEN],word[MAX_LINE_LEN];
  ATOM_TYPING_SCHEME *scheme;
  FORCE_FIELD *ff;
  ATOM_TYPE *type1,*type2;
  NONBONDED_FF *nonbonded_ff;

  scheme = settings->atom_typing_scheme;
  ff = settings->force_field;

  strcpy(word,"");

  while (strcmp(word,"END")) {
  
    if (fgets(line,MAX_LINE_LEN,stdin) == NULL) {

      break;
    }

    if (sscanf(line,"%s",word) == EOF) {

      break;
    }

    if (!strcmp(word,"CONTACT")) {

      n_read = sscanf(line,"%*s %d %d %d %d %lf %lf %lf %lf %lf %lf %d",
		      &type_id1,&element_id1,&type_id2,&element_id2,
		      &area,&distance,&alpha1,&beta1,&alpha2,&beta2,&key);

      if (n_read == 11) {

	type1 = pliff_qscore_get_atom_type(type_id1,element_id1,scheme);
	type2 = pliff_qscore_get_atom_type(type_id2,element_id2,scheme);

	if ((type1 != NULL) && (type2 != NULL)) {

	  vdw_radius1 = get_atom_type_vdw_radius(type1);
	  vdw_radius2 = get_atom_type_vdw_radius(type2);

	  nonbonded_ff = get_nonbonded_ff(ff->nonbonded,type1,type2);

	  contact_score = nonbonded_ff->lnP;

	  geometry_score = pliff_score_vs_histogram(nonbonded_ff->R_histogram_id,distance-vdw_radius1-vdw_radius2,0);

	  if ((alpha1 >= 0.0) && (nonbonded_ff->ALPHA1_histogram_id != -1)) {

	    geometry_score += pliff_score_vs_histogram(nonbonded_ff->ALPHA1_histogram_id,alpha1,0);
	  }

	  if ((beta1 >= 0.0) && (nonbonded_ff->BETA1_histogram_id != -1)) {

	    geometry_score += pliff_score_vs_histogram(nonbonded_ff->BETA1_histogram_id,beta1,0);
	  }

	  if ((alpha2 >= 0.0) && (nonbonded_ff->ALPHA2_histogram_id != -1)) {

	    geometry_score += pliff_score_vs_histogram(nonbonded_ff->ALPHA2_histogram_id,alpha2,0);
	  }

	  if ((beta2 >= 0.0) && (nonbonded_ff->BETA2_histogram_id != -1)) {

	    geometry_score += pliff_score_vs_histogram(nonbonded_ff->BETA2_histogram_id,beta2,0);
	  }

	  score = area*(contact_score + geometry_score);

	  printf("SCORE %10.4lf %10.4lf %10.4lf %10d\n",contact_score,geometry_score,score,key);

	} else {

	  warning_fn("qscore: atom type(s) undefined");
	}

      } else {

	warning_fn("qscore: line does not have correct number of fields");
      }
    }
  }
}



void pliff_minimise_system(SYSTEM *system) {

  PLI_SFUNC *sfunc;

  sfunc = get_sfunc("pliff");

  minimise(system,sfunc,-1);
}



void pliff_score_system(SYSTEM *system) {

  int i,j,id1,id2,n_contacts;
  unsigned int cflags;
  double *scores,*cscores,p;
  ATOM **atomp,*atom;
  ATOMLIST *selection;
  CONTACT *contact;
  CONTACTLIST *contactlist;
  SETTINGS *settings;
  ATOM_TYPING_SCHEME *scheme;
  FORCE_FIELD *ff;

  settings = system->settings;
  scheme = settings->atom_typing_scheme;
  ff = settings->force_field;
  selection = system->selection;

  if (selection == NULL) {

    return;
  }

  // calculate (new) contacts:

  if (pliff_recalc_contacts) {

    set_contacts_system(system,0);
  }

  // calculate (new) geometries and scores:

  pliff_score_contacts(selection,scheme,ff);

  pliff_scale_contacts(selection,scheme,ff);

  pliff_score_atoms(selection);

  // add up the scores:

  init_system_scores(system);

  scores = system->scores;

  p = 1.0;

  n_contacts = 0;

  for (i=0,atomp=selection->atom;i<selection->natoms;i++,atomp++) {
  
    atom = *atomp;

    if (!(atom->flags & SKIP_ATOM)) {

      contactlist = atom->contactlist;

      if (contactlist != NULL) {

	id1 = atom->unique_id;
	
	for (j=0,contact=contactlist->contacts;j<contactlist->ncontacts;j++,contact++) {

	  cflags = contact->flags;
	  
	  if ((!(cflags & COVALENT_CONTACT)) && (!(cflags & INTRAMOLECULAR_CONTACT)) && (!(cflags & SECONDARY_CONTACT))) {
	    
	    id2 = contact->atom2->unique_id;
	    
	    if ((id2 > id1) || (!atom_in_list(selection,contact->atom2))) {

	      cscores = contact->scores;
	      scores[PLIFF_SYSTEM_CONTACT_SCORE] += (cscores[PLIFF_CONTACT_SCALE]*cscores[PLIFF_CONTACT_SCORE]);
	      scores[PLIFF_SYSTEM_GEOMETRY_SCORE] += (cscores[PLIFF_CONTACT_SCALE]*cscores[PLIFF_GEOMETRY_SCORE]);

	      n_contacts++;
	    }
	  }
	}
      }

      scores[PLIFF_SYSTEM_TRIPLET_SCORE] += atom->scores[PLIFF_ATOM_TRIPLET_SCORE];

      p *= atom->type_probability;
    }
  }

  // this is a bit of a hack - molsystems won't have a ligand, so won't get scored here
  // the idea being that a ligand in isolation can adopt a conformation with zero
  // clashes and torsional issues. this will need to be addressed at some point, because:
  // (1) the assumption won't always hold true and (2) it won't allow us to minimise a
  // ligand molecule in isolation

  if (system->ligand) {

    scores[PLIFF_SYSTEM_LIGAND_SCORE] = molecule_internal_energy(system->ligand,scheme,ff);
  }

  //scores[PLIFF_SYSTEM_TC_SCORE] = -0.5925*log(p);

  scores[PLIFF_SYSTEM_GEOMETRY_SCORE] += scores[PLIFF_SYSTEM_TRIPLET_SCORE];

  scores[PLIFF_SYSTEM_NB_SCORE] = scores[PLIFF_SYSTEM_CONTACT_SCORE] + scores[PLIFF_SYSTEM_GEOMETRY_SCORE];

  scores[PLIFF_SYSTEM_SCORE] = scores[PLIFF_SYSTEM_CONTACT_SCORE] + scores[PLIFF_SYSTEM_GEOMETRY_SCORE] + scores[PLIFF_SYSTEM_RB_SCORE] + scores[PLIFF_SYSTEM_LIGAND_SCORE];

  system->score = scores[PLIFF_SYSTEM_SCORE];
}



static void pliff_score_contacts(ATOMLIST *selection,ATOM_TYPING_SCHEME *scheme,FORCE_FIELD *ff) {

  int i,j,id1,id2;
  unsigned int flags1,flags2,status1,status2,cflags;
  ATOM **atomp1,*atom1,*atom2;
  CONTACT *contact,*icontact;
  CONTACTLIST *contactlist;

  for (i=0,atomp1=selection->atom;i<selection->natoms;i++,atomp1++) {
  
    atom1 = *atomp1;

    id1 = atom1->unique_id;
    flags1 = atom1->flags;
    status1 = atom1->status;

    contactlist = atom1->contactlist;

    if (contactlist != NULL) {

      for (j=0,contact=contactlist->contacts;j<contactlist->ncontacts;j++,contact++) {

	cflags = contact->flags;

	atom2 = contact->atom2;
	id2 = atom2->unique_id;
	flags2 = atom2->flags;
	status2 = atom2->status;

	if ((!(cflags & COVALENT_CONTACT)) && (!(cflags & SECONDARY_CONTACT)) && ((!(cflags & INTRAMOLECULAR_CONTACT)) || (hbond_flags_match(flags1,flags2)))) {

	  if ((id2 > id1) || (!atom_in_list(selection,atom2))) {

	    if ((status1 != ATOM_SCORED) || ((atom_in_list(selection,atom2)) && (status2 != ATOM_SCORED))) {

	      calc_contact_geometry(contact);

	      pliff_score_contact(contact,scheme,ff);

	      icontact = find_contact(atom2,atom1);

	      if (icontact) {

		mirror_contact_geometry(contact,icontact);

		mirror_contact_scores(contact,icontact);
	      }
	    }
	  }
	}
      }
    }
  }
}



static void pliff_score_atoms(ATOMLIST *selection) {

  int i,j;
  unsigned int cflags;
  double *scores,*cscores;
  ATOM **atomp,*atom;
  CONTACT *contact;
  CONTACTLIST *contactlist;

  for (i=0,atomp=selection->atom;i<selection->natoms;i++,atomp++) {
  
    atom = *atomp;

    if (atom->status != ATOM_SCORED) {

      scores = atom->scores;

      scores[PLIFF_ATOM_CONTACT_SCORE] = 0.0;
      scores[PLIFF_ATOM_GEOMETRY_SCORE] = 0.0;

      contactlist = atom->contactlist;

      if (contactlist != NULL) {

	for (j=0,contact=contactlist->contacts;j<contactlist->ncontacts;j++,contact++) {

	  cflags = contact->flags;

	  if ((!(cflags & COVALENT_CONTACT)) && (!(cflags & INTRAMOLECULAR_CONTACT)) && (!(cflags & SECONDARY_CONTACT))) {
	    
	    cscores = contact->scores;

	    scores[PLIFF_ATOM_CONTACT_SCORE] += (cscores[PLIFF_CONTACT_SCALE]*cscores[PLIFF_CONTACT_SCORE]);
	    scores[PLIFF_ATOM_GEOMETRY_SCORE] += (cscores[PLIFF_CONTACT_SCALE]*cscores[PLIFF_GEOMETRY_SCORE]);
	  }
	}
      }

      scores[PLIFF_ATOM_GEOMETRY_SCORE] += scores[PLIFF_ATOM_TRIPLET_SCORE];

      scores[PLIFF_ATOM_NB_SCORE] = scores[PLIFF_ATOM_CONTACT_SCORE] + scores[PLIFF_ATOM_GEOMETRY_SCORE];

      scores[PLIFF_ATOM_SCORE] = scores[PLIFF_ATOM_CONTACT_SCORE] + scores[PLIFF_ATOM_GEOMETRY_SCORE];

      atom->score = scores[PLIFF_ATOM_SCORE];

      atom->cstatus |= ATOM_SCORES_CALCULATED;

      atom->status = ATOM_SCORED;
    }
  }
}



static void pliff_scale_contacts(ATOMLIST *selection,ATOM_TYPING_SCHEME *scheme,FORCE_FIELD *ff) {

  int i,j;
  ATOM **atomp,*atom;
  ATOM_TYPE *type1,*type2;
  CONTACTLIST *contactlist;
  CONTACT *contact;
  NONBONDED_FF *nonbonded_ff;
  HISTOGRAM *hist_ALPHA1,*hist_ALPHA2;

  for (i=0,atomp=selection->atom;i<selection->natoms;i++,atomp++) {
  
    atom = *atomp;

    contactlist = atom->contactlist;

    if (contactlist != NULL) {

      type1 = atom->type;

      for (j=0,contact=contactlist->contacts;j<contactlist->ncontacts;j++,contact++) {

	type2 = contact->atom2->type;

	nonbonded_ff = get_nonbonded_ff(ff->nonbonded,type1,type2);

	hist_ALPHA1 = get_histogram(nonbonded_ff->ALPHA1_histogram_id);
	hist_ALPHA2 = get_histogram(nonbonded_ff->ALPHA2_histogram_id);

	contact->scores[PLIFF_CONTACT_SCALE] = pliff_init_contact_scale(contact,nonbonded_ff,hist_ALPHA1,hist_ALPHA2);
      }
    }
  }

  for (i=0,atomp=selection->atom;i<selection->natoms;i++,atomp++) {
  
    atom = *atomp;

    pliff_scale_atom_contacts(atom,scheme,ff);
  }
}



void pliff_write_system_scores(PLI_FILE *file,SYSTEM *system,enum OUTPUT_FORMAT oformat,unsigned long int oflags) {

  double *s,*r;

  s = system->scores;
  r = system->ref_scores;

  write_line(file,(oformat == JSON) ? 
	     "\"pliff_score\":%lf,\"pliff_nb_score\":%lf,\"pliff_cscore\":%lf,\"pliff_gscore\":%lf,\"pliff_tscore\":%lf,\"pliff_iscore\":%lf" : 
	     " %10.4lf %10.4lf %10.4lf %10.4lf %10.4lf %10.4lf",
	     s[PLIFF_SYSTEM_SCORE]-r[PLIFF_SYSTEM_SCORE],
             s[PLIFF_SYSTEM_NB_SCORE]-r[PLIFF_SYSTEM_NB_SCORE],
	     s[PLIFF_SYSTEM_CONTACT_SCORE]-r[PLIFF_SYSTEM_CONTACT_SCORE],
	     s[PLIFF_SYSTEM_GEOMETRY_SCORE]-r[PLIFF_SYSTEM_GEOMETRY_SCORE],
	     s[PLIFF_SYSTEM_TRIPLET_SCORE]-r[PLIFF_SYSTEM_TRIPLET_SCORE],
	     s[PLIFF_SYSTEM_LIGAND_SCORE]
             );
}



void pliff_write_atom_scores(PLI_FILE *file,ATOM *atom,enum OUTPUT_FORMAT oformat,unsigned long int oflags) {

  double *s,*r;

  s = atom->scores;
  r = atom->ref_scores;

  write_line(file,(oformat == JSON) ? 
	     "\"pliff_score\":%lf,\"pliff_cscore\":%lf,\"pliff_gscore\":%lf,\"pliff_tscore\":%lf,\"constraint_score\":%lf" : 
	     " %10.4lf %10.4lf %10.4lf %10.4lf %10.4lf",
	     s[PLIFF_ATOM_SCORE] - r[PLIFF_ATOM_SCORE],
	     s[PLIFF_ATOM_CONTACT_SCORE] - r[PLIFF_ATOM_CONTACT_SCORE],
	     s[PLIFF_ATOM_GEOMETRY_SCORE] - r[PLIFF_ATOM_GEOMETRY_SCORE],
	     s[PLIFF_ATOM_TRIPLET_SCORE] - r[PLIFF_ATOM_TRIPLET_SCORE],
	     atom->constraint_score);
}



void pliff_write_contact_scores(PLI_FILE *file,CONTACT *contact,enum OUTPUT_FORMAT oformat,unsigned long int oflags) {

  double *s;

  s = contact->scores;

  write_line(file,(oformat == JSON) ? 
	     "\"pliff_score\":%lf,\"pliff_cscore\":%lf,\"pliff_gscore\":%lf" : 
	     " %10.4lf %10.4lf %10.4lf",
	     contact->score,
	     s[PLIFF_CONTACT_SCALE]*s[PLIFF_CONTACT_SCORE],
	     s[PLIFF_CONTACT_SCALE]*s[PLIFF_GEOMETRY_SCORE]);
}



void pliff_score_contact(CONTACT *contact,ATOM_TYPING_SCHEME *scheme,FORCE_FIELD *ff) {

  int within_limits1,within_limits2,i,j;
  double D,A,f_r_alpha,f_r_beta,*scores,fg[MAX_CONTACT_SCORES];
  ATOM *atom1,*atom2;
  ATOM_TYPE *type1,*type2;
  ATOM_GEOMETRY *geometry1,*geometry2;
  NONBONDED_FF *nonbonded_ff;
  HISTOGRAM *hist_R,*hist_ALPHA1,*hist_BETA1,*hist_ALPHA2,*hist_BETA2;

  atom1 = contact->atom1;
  atom2 = contact->atom2;

  type1 = atom1->type;
  type2 = atom2->type;

  D = (contact->distance)-(atom1->vdw_radius)-(atom2->vdw_radius);
  A = (pliff_use_areas) ? (contact->area)+(contact->iarea) : 1.0;

  nonbonded_ff = get_nonbonded_ff(ff->nonbonded,type1,type2);

  geometry1 = atom1->geometry;
  geometry2 = atom2->geometry;

  if (geometry1 == NULL) {

    error_fn("pliff_score_contact: atom geometry undefined for atom %s (%d)",atom1->name,atom1->id);
  }

  if (geometry2 == NULL) {

    error_fn("pliff_score_contact: atom geometry undefined for atom %s (%d)",atom2->name,atom2->id);
  }

  within_limits1 = (atom1->flags & SINGLE_ATOM_PROBE) ? 1 : 0;
  within_limits2 = (atom2->flags & SINGLE_ATOM_PROBE) ? 1 : 0;

  within_limits1 = within_limits2 = 0;

  f_r_alpha = r_alpha_correction(D);
  f_r_beta = r_beta_correction(D);

  // set pointers to histograms:

  hist_R      = get_histogram(nonbonded_ff->R_histogram_id);
  hist_ALPHA1 = get_histogram(nonbonded_ff->ALPHA1_histogram_id);
  hist_BETA1  = get_histogram(nonbonded_ff->BETA1_histogram_id);
  hist_ALPHA2 = get_histogram(nonbonded_ff->ALPHA2_histogram_id);
  hist_BETA2  = get_histogram(nonbonded_ff->BETA2_histogram_id);

  if (hist_R->maxY < 0.0) {

    printf("max Y value negative for R (%d,%d)\n",type1->id,type2->id);
  }

  // pointer to scores, and initialise:

  scores = contact->scores;

  for (i=0;i<MAX_CONTACT_SCORES;i++) {

    scores[i] = fg[i] = 0.0;
  }

  // calculate R score:

  if (((contact->distance < nonbonded_ff->Dc) && (nonbonded_ff->use_clash_potential)) || (D < hist_R->startX)) {

    scores[PLIFF_R_SCORE] = lennard_jones_energy(contact->distance,nonbonded_ff->LJ_A,nonbonded_ff->LJ_B);

  } else {

    scores[PLIFF_R_SCORE] = pliff_score_vs_histogram(nonbonded_ff->R_histogram_id,D,0);
  }

  // return here for long contacts if they are to be ignored:

  if ((pliff_ignore_long_contacts) && (contact->distance > nonbonded_ff->Do) && (scores[PLIFF_R_SCORE] < 0.0)) {

    scores[PLIFF_R_SCORE] = 0.0;

    return;
  }

  fg[PLIFF_R_SCORE] = 1.0;

  // calculate ALPHA1 score:

  if ((geometry1->u_axis) && (nonbonded_ff->ALPHA1_histogram_id != -1)) {

    if (hist_ALPHA1->maxY < 0.0) {

      printf("max Y value negative for ALPHA1 (%d,%d)\n",type1->id,type2->id);
    }

    scores[PLIFF_ALPHA1_SCORE] = pliff_score_vs_histogram(nonbonded_ff->ALPHA1_histogram_id,contact->alpha1,within_limits1);

    fg[PLIFF_ALPHA1_SCORE] = f_r_alpha;
  }

  // calculate BETA1 score:

  if ((geometry1->v_axis) && (nonbonded_ff->BETA1_histogram_id != -1)) {

    if (hist_BETA1->maxY < 0.0) {

      printf("max Y value negative for BETA1 (%d,%d)\n",type1->id,type2->id);
    }

    scores[PLIFF_BETA1_SCORE] = pliff_score_vs_histogram(nonbonded_ff->BETA1_histogram_id,contact->beta1,within_limits1);

    fg[PLIFF_BETA1_SCORE] = f_r_beta*alpha_beta_correction(contact->alpha1);
  }

  // calculate ALPHA2 score:

  if ((geometry2->u_axis) && (nonbonded_ff->ALPHA2_histogram_id != -1)) {

    if (hist_ALPHA2->maxY < 0.0) {

      printf("max Y value negative for ALPHA2 (%d,%d)\n",type1->id,type2->id);
    }

    scores[PLIFF_ALPHA2_SCORE] = pliff_score_vs_histogram(nonbonded_ff->ALPHA2_histogram_id,contact->alpha2,within_limits2);

    fg[PLIFF_ALPHA2_SCORE] = f_r_alpha;
  }

  // calculate BETA2 score:

  if ((geometry2->v_axis) && (nonbonded_ff->BETA2_histogram_id != -1)) {

    if (hist_BETA2->maxY < 0.0) {

      printf("max Y value negative for BETA2 (%d,%d)\n",type1->id,type2->id);
    }

    scores[PLIFF_BETA2_SCORE] = pliff_score_vs_histogram(nonbonded_ff->BETA2_histogram_id,contact->beta2,within_limits2);

    fg[PLIFF_BETA2_SCORE] = f_r_beta*alpha_beta_correction(contact->alpha2);
  }

  // combine geometry scores:
  //
  // each geometry parameter has a geometry correction factor fg. for example, at longer distance, the alpha and beta distributions
  // become less pronounced, simply because the effect will be weaker.

  scores[PLIFF_R_SCORE] *= -(pliff_geometry_coeff*A);

  scores[PLIFF_GEOMETRY_SCORE] = scores[PLIFF_R_SCORE];

  for (i=PLIFF_ALPHA1_SCORE;i<=PLIFF_BETA2_SCORE;i++) {
    
    scores[i] *= -(pliff_geometry_coeff*A*fg[i]);

    scores[PLIFF_GEOMETRY_SCORE] += scores[i];
  }

  // calculate contact score:

  scores[PLIFF_CONTACT_SCORE] = -pliff_contact_coeff*A*(nonbonded_ff->lnP);

  // for hydrogen-bonding contacts, also calculate hbond geometry score:

  if (hbond_flags_match(type1->flags,type2->flags)) {

    scores[PLIFF_HBOND_GEOMETRY_SCORE] = hbond_geometry_score(contact);
  }

  // sum terms:

  contact->contact_score = scores[PLIFF_CONTACT_SCORE];
  contact->geometry_score = scores[PLIFF_GEOMETRY_SCORE];

  contact->score = scores[PLIFF_CONTACT_SCORE] + scores[PLIFF_GEOMETRY_SCORE];
}



static void pliff_scale_atom_contacts(ATOM *atom,ATOM_TYPING_SCHEME *scheme,FORCE_FIELD *ff) {

  int i,n_hs,n_lps,n_max,n_hbonds,hbonds_on[MAX_Z];
  unsigned int flags;
  double score,*scores,Q,EtQ,Qhbs[MAX_Z],f;
  ATOM_TYPE *type;
  CONTACTLIST *contactlist;
  CONTACT **hbondp,*hbonds[MAX_Z],*hbond,*ihbond;

  type = atom->type;
  flags = type->flags;

  scores = atom->scores;

  scores[PLIFF_ATOM_TRIPLET_SCORE] = 0.0;

  contactlist = atom->contactlist;

  if (!contactlist) {

    return;
  }

  if ((flags & HBOND_DA_ATOM_TYPE) == HBOND_DA_ATOM_TYPE) {

    n_hs = type->n_hs;
    n_lps = type->n_lps;
      
    if ((n_hs < 1) || (n_lps < 1)) {
      
      error_fn("pliff_scale_atom_contacts: hydrogens or lone pairs unknown for atom %d in molecule %s (do you need to resolve types?)",
	       atom->id,atom->molecule->name);
    }
 
    n_max = pliff_get_da_hbonds(contactlist,n_hs,n_lps,hbonds,&n_hbonds);

    if ((n_max > 0) && (n_max < n_hbonds)) {
            
      for (i=0;i<n_hbonds;i++) {
	
	hbonds_on[i] = 0;

	Qhbs[i] = 0.0;
      }
      
      Q = EtQ = 0.0;
      
      pliff_permute_da_hbonds(hbonds,n_hbonds,n_hs,n_lps,hbonds_on,0,n_max,Qhbs,&Q,&EtQ);
      
      scores[PLIFF_ATOM_TRIPLET_SCORE] = EtQ/Q;

      for (i=0,hbondp=hbonds;i<n_hbonds;i++,hbondp++) {

	hbond = *hbondp;

	f = Qhbs[i]/Q;

	hbond->scores[PLIFF_CONTACT_SCALE] *= f;

	ihbond = find_contact(hbond->atom2,hbond->atom1);

	if (ihbond) {

	  ihbond->scores[PLIFF_CONTACT_SCALE] *= f;
	}
      }

    } else if (n_hbonds > 1) {

      scores[PLIFF_ATOM_TRIPLET_SCORE] = pliff_hbonds_triplet_score(hbonds,n_hbonds,NULL);
    }

  } else if (flags & METAL_ATOM_TYPE) {

    // need to do something with metals here to sort out Acc...Me...Acc angles and coordination numbers
  }
}



static int pliff_get_da_hbonds(CONTACTLIST *list,int n_hs,int n_lps,CONTACT **contacts,int *n_contacts) {

  int i,n_da,n_acc,n_don,n_h_hbonds,n_lp_hbonds,n_hbonds,max_hbonds;
  unsigned int flags;
  CONTACT **contactp,*contact;

  pliff_contactlist2hbonds(list,contacts,n_contacts,0.0,0.0);

  n_da = n_acc = n_don = 0;

  for (i=0,contactp=contacts;i<*n_contacts;i++,contactp++) {

    contact = *contactp;

    flags = contact->atom2->type->flags;

    if (flags & HBOND_ACCEPTOR_ATOM_TYPE) {

      if (flags & HBOND_DONOR_ATOM_TYPE) {

	n_da++;

      } else {

	n_acc++;
      }
	
    } else if ((flags & HBOND_DONOR_ATOM_TYPE) || (flags & CHBOND_DONOR_ATOM_TYPE) || (flags & METAL_ATOM_TYPE)) {

      n_don++;
    }
  }

  n_h_hbonds = (n_acc > n_hs) ? n_hs : n_acc;
  n_lp_hbonds = (n_don > n_lps) ? n_lps : n_don;

  n_hbonds = n_h_hbonds + n_lp_hbonds;

  max_hbonds = n_hs + n_lps;

  if (n_hbonds < max_hbonds) {

    n_hbonds = (n_hbonds+n_da > max_hbonds) ? max_hbonds : n_hbonds + n_da;
  }

  return(n_hbonds);
}



static void pliff_permute_da_hbonds(CONTACT **contacts,int n_contacts,int n_hs,int n_lps,int *hbonds,int cid,int max_hbonds,double *Qhbs,double *Q,double *EtQ) {

  int i,n_h_hbonds,n_lp_hbonds,n_hbonds;
  double E,Etriplet,q;
  unsigned int flags;
  CONTACT **contactp,*contact;

  // count number of hbonds of each type:

  n_h_hbonds = n_lp_hbonds = 0;

  for (i=0;i<cid;i++) {

    if (hbonds[i] == 1) {

      n_h_hbonds++;

    } else if (hbonds[i] == 2) {

      n_lp_hbonds++;
    }
  }

  n_hbonds = n_h_hbonds + n_lp_hbonds;

  // return if current permutation corresponds to n_contacts:

  if (cid == n_contacts) {

    if (n_hbonds == max_hbonds) {
      
      E = pliff_score_da_permutation(contacts,n_contacts,hbonds,&Etriplet);

      q = exp(-E/RT);

      for (i=0;i<n_contacts;i++) {
	
	if (hbonds[i]) {
	  
	  Qhbs[i] += q;
	}
      }
      
      (*Q) += q;
      
      (*EtQ) += (Etriplet*q);
    }

    return;
  }

  contact = contacts[cid];

  flags = contact->atom2->type->flags;

  // set hbond type for current contact to none:

  hbonds[cid] = 0;

  pliff_permute_da_hbonds(contacts,n_contacts,n_hs,n_lps,hbonds,cid+1,max_hbonds,Qhbs,Q,EtQ);

  if ((n_h_hbonds < n_hs) && (flags & HBOND_ACCEPTOR_ATOM_TYPE)) {

    // set hbond type for current contact to donor:

    hbonds[cid] = 1;

    pliff_permute_da_hbonds(contacts,n_contacts,n_hs,n_lps,hbonds,cid+1,max_hbonds,Qhbs,Q,EtQ);
  }

  if ((n_lp_hbonds < n_lps) && ((flags & HBOND_DONOR_ATOM_TYPE) || (flags & CHBOND_DONOR_ATOM_TYPE) || (flags & METAL_ATOM_TYPE))) {

    // set hbond type for current contact to acceptor:

    hbonds[cid] = 2;

    pliff_permute_da_hbonds(contacts,n_contacts,n_hs,n_lps,hbonds,cid+1,max_hbonds,Qhbs,Q,EtQ);
  }
}




// score this da hbonding configuration
// currently non-hbonding contacts are not counted
// should probably apply a penalty to these

static double pliff_score_da_permutation(CONTACT **contacts,int n_contacts,int *hbonds,double *Etriplet) {

  int i;
  double score,triplet_score;
  CONTACT **contactp;

  score = 0.0;

  for (i=0,contactp=contacts;i<n_contacts;i++,contactp++) {

    if (hbonds[i]) {

      score += (*contactp)->score;
    }
  }

  triplet_score = pliff_hbonds_triplet_score(contacts,n_contacts,hbonds);

  *Etriplet = triplet_score;

  return(score + triplet_score);
}



static double pliff_hbonds_triplet_score(CONTACT **contacts,int n_contacts,int *hbonds) {

  int i,j;
  double triplet_score,area1,area2,*pos0,*pos1,*pos2,v1[4],v2[4],angle,f1,f2;
  CONTACT **contactp1,**contactp2,*contact1,*contact2;

  // calculate triplet score:

  triplet_score = 0.0;

  pos0 = (*contacts)->atom1->position;

  for (i=0,contactp1=contacts;i<n_contacts-1;i++,contactp1++) {

    if ((!hbonds) || (hbonds[i])) {

      contact1 = *contactp1;

      area1 = (pliff_use_areas) ? contact1->area + contact1->iarea : 1.0;

      pos1 = contact1->atom2->position;
 
      calc_vector(pos0,pos1,v1);

      f1 = contact1->scores[PLIFF_HBOND_GEOMETRY_SCORE];

      for (j=i+1,contactp2=contacts+i+1;j<n_contacts;j++,contactp2++) {
	
	if ((!hbonds) || (hbonds[j])) {
	  
	  contact2 = *contactp2;

	  area2 = (pliff_use_areas) ? contact2->area + contact2->iarea : 1.0;

	  pos2 = contact2->atom2->position;

	  f2 = contact2->scores[PLIFF_HBOND_GEOMETRY_SCORE];

	  calc_vector(pos0,pos2,v2);

	  angle = vector_angle(v1,v2);

	  triplet_score += 0.5*(area1+area2)*f1*f2*pliff_triplet_score(angle);
	}
      }
    }
  }

  return(pliff_geometry_coeff*triplet_score);
}



static double pliff_triplet_score(double angle) {

  return(sqr((105-angle)/40)-0.92);
}



static void pliff_contactlist2hbonds(CONTACTLIST *list,CONTACT **hbonds,int *n_hbonds,double min_geom_score,double max_score) {

  int i;
  CONTACT *contact;

  *n_hbonds = 0;

  for (i=0,contact=list->contacts;i<list->ncontacts;i++,contact++) {

    if (hbond_flags_match(contact->atom1->type->flags,contact->atom2->type->flags)) { 

      if ((contact->score < max_score) && (contact->scores[PLIFF_HBOND_GEOMETRY_SCORE] > min_geom_score)) {

	if (*n_hbonds == MAX_Z) {

	  error_fn("contactlist2hbonds: too many hydrogen bonds");
	}

	hbonds[*n_hbonds] = contact;

	(*n_hbonds)++;
      }
    }
  }
}


static double pliff_score_vs_histogram(int histogram_id,double X,int within_limits) {

  int i,bin;
  double X1,X2,Y1,Y2,Y;
  HISTOGRAM *histogram;
  HISTOGRAM_POINT *points,*point,*min_point;

  histogram = get_histogram(histogram_id);

  if (histogram == NULL) {

    error_fn("score_vs_histogram: histogram %d undefined",histogram_id);
  }

  if ((within_limits) && ((X < histogram->limits[0]) || (X > histogram->limits[1]))) {

    return(0.0);
  }

  Y = histogram_X2Y(histogram,X);
  
  return(Y);
}



static ATOM_TYPE* pliff_qscore_get_atom_type(int type_id,int element_id,ATOM_TYPING_SCHEME *scheme) {

  int i;
  ELEMENT *element;
  ATOM_TYPE *type;

  if (type_id == -1) {

    element = get_element_by_id(element_id,scheme);

    for (i=0,type=scheme->atom_types;i<scheme->n_atom_types;i++,type++) {

      if ((type->id == -1) && (type->element == element)) {

	return(type);
      }
    }

  } else {

    type = get_atom_type_by_id(type_id,scheme);

    return(type);
  }

  return(NULL);
}



static double pliff_init_contact_scale(CONTACT *contact,NONBONDED_FF *nonbonded_ff,HISTOGRAM *hist_ALPHA1,HISTOGRAM *hist_ALPHA2) {

  int atom_probe1,atom_probe2,within_limits1,within_limits2;
  double alpha1,alpha2,Do,fR;

  atom_probe1 = (int) (contact->atom1->flags & SINGLE_ATOM_PROBE);
  atom_probe2 = (int) (contact->atom2->flags & SINGLE_ATOM_PROBE);

  if ((!atom_probe1) && (!atom_probe2)) {

    return(1.0);
  }

  if ((contact->distance > nonbonded_ff->Do) && (contact->scores[PLIFF_R_SCORE] < 0.0)) {

    return(0.0);
  }

  alpha1 = contact->alpha1;
  alpha2 = contact->alpha2;

  within_limits1 = ((atom_probe1) && (hist_ALPHA1) && ((alpha1 < hist_ALPHA1->limits[0]) || (alpha1 > hist_ALPHA1->limits[1]))) ? 0 : 1;
  within_limits2 = ((atom_probe2) && (hist_ALPHA2) && ((alpha2 < hist_ALPHA2->limits[0]) || (alpha2 > hist_ALPHA2->limits[1]))) ? 0 : 1;

  if ((!within_limits1) || (!within_limits2)) {
   
    Do = nonbonded_ff->Do;

    fR = ramp_function(contact->distance,Do+0.0,Do+0.5,1.0,0.0);

    return(fR);
  }

  return(1.0);
}
