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


#define MIN_SURROGATE_SHIFT_Y 0.10
#define MIN_AMBIGUOUS_PEAK_DIST 0.30
#define MIN_AMBIGUOUS_PEAK_RATIO 0.7
#define WARNING_AMBIGUOUS_PEAK_DIST 0.25
#define WARNING_AMBIGUOUS_PEAK_RATIO 0.6
#define MIN_K1 2.0
#define MIN_K2 1.0
#define MAX_LR_A 5.0
#define MIN_K_RATIO 1.0
#define MAX_K_RATIO 4.0
#define MAX_REFINE_SHIFT_X 0.2
#define MAX_REFINE_SHIFT_Y 1.0
#define REFINE_DX 0.01
#define REFINE_DY 0.05
#define DELTA_D_SLOPE 0.75



static struct FFSettings {
  int use_alpha_potentials;
  int use_beta_potentials;
  double (*g_r_alpha)(double);
  double (*g_r_beta)(double);
  double (*g_alpha_beta)(double);
  enum { BLOCK_R_ALPHA_CORRECTION,NO_R_ALPHA_CORRECTION } r_alpha_correction;
  enum { BLOCK_R_BETA_CORRECTION,NO_R_BETA_CORRECTION } r_beta_correction;
  enum { COSINE_ALPHA_BETA_CORRECTION,BLOCK_ALPHA_BETA_CORRECTION,NO_ALPHA_BETA_CORRECTION } alpha_beta_correction;
  double r_alpha_correction_r1;
  double r_alpha_correction_r2;
  double r_beta_correction_r1;
  double r_beta_correction_r2;
  double alpha_beta_correction_alpha1;
  double alpha_beta_correction_alpha2;
  int max_stored_potentials;
} *ff_settings = NULL;



static FORCE_FIELD *pliff_ff = NULL;
static NONBONDED_FF *last_read_ff = NULL;
static int n_stored_potentials = 0;



static void set_surrogate_R_potential(NONBONDED_FF*);
static void copy_R_potential_params(NONBONDED_FF*,NONBONDED_FF*,double);
static NONBONDED_FF* get_surrogate_nbff(NONBONDED_FF*);
static void reset_nonbonded_ff(NONBONDED_FF*);
static NONBONDED_FF* get_least_used_nonbonded_ff(FORCE_FIELD*);
static HISTOGRAM* process_atom_ff_histogram(HISTOGRAM*,ATOM_FF*);
static void process_atom_ff_histogram_old(HISTOGRAM*,ATOM_FF*);
static void process_nonbonded_ff_histogram(HISTOGRAM*,NONBONDED_FF*);

static int process_R_potential(NONBONDED_FF*);
static void init_R_potential(NONBONDED_FF*);
static void finish_R_potential(NONBONDED_FF*);
static int smoothe_R_potential(NONBONDED_FF*);
static int fit_R_potential(NONBONDED_FF*);

static void apply_R_potential_softness(NONBONDED_FF*);
static void calc_R_potential_fit_quality(HISTOGRAM*,NONBONDED_FF*);
static void normalise_R_histogram(HISTOGRAM*,NONBONDED_FF*);
static double calc_fitted_R_potential(double,NONBONDED_FF*);

static int fit_R_potential_full(HISTOGRAM*,NONBONDED_FF*);
static int get_R_optimum_from_points(HISTOGRAM*,NONBONDED_FF*);
static int refine_R_optimum(HISTOGRAM*,NONBONDED_FF*);
static int fit_short_range_potential(HISTOGRAM*,NONBONDED_FF*);
static int fit_long_range_potential(HISTOGRAM*,NONBONDED_FF*);
static double calc_short_range_k(HISTOGRAM*,double,double);

static void recalc_R_histogram_errors(HISTOGRAM*);
static double calc_2d_correction(int,double (*g)(double));
static void apply_2d_correction(int,double);

static void read_atom_ff(ATOM_FF*);
static int read_nonbonded_ff(NONBONDED_FF*);
static void copy_nonbonded_ff(NONBONDED_FF*,NONBONDED_FF*);
static void init_ami_ff(FORCE_FIELD*);
static void init_atom_ff_list(ATOM_FF*,ATOM_TYPING_SCHEME*,FORCE_FIELD*);
static void init_nonbonded_ff_matrix(NONBONDED_FF**,ATOM_TYPING_SCHEME*,FORCE_FIELD*);
static void init_atom_ff(ATOM_FF*,ATOM_TYPE*,FORCE_FIELD*);
static void init_nonbonded_ff(NONBONDED_FF*,ATOM_TYPE*,ATOM_TYPE*,FORCE_FIELD*);
static FORCE_FIELD* alloc_ff(ATOM_TYPING_SCHEME*);
static NONBONDED_FF** alloc_nonbonded_ff_matrix(ATOM_TYPING_SCHEME*);
static ATOM_FF* alloc_atom_ff_list(ATOM_TYPING_SCHEME*);
static void free_ff(FORCE_FIELD*);
static void free_nonbonded_ff_matrix(NONBONDED_FF**);



FORCE_FIELD* get_ff(void) {

  return(pliff_ff);
}



void init_ff_settings(void) {

  char *corr;

  if (ff_settings == NULL) {

    ff_settings = (struct FFSettings*) malloc(sizeof(struct FFSettings));
  }

  if (ff_settings == NULL) {

    error_fn("init_ff_settings: out of memory allocating settings");
  }

  // use alpha, beta and clash potentials?

  ff_settings->use_alpha_potentials = (params_get_parameter("use_alpha_potentials"))->value.i;
  ff_settings->use_beta_potentials = (params_get_parameter("use_beta_potentials"))->value.i;

  // pointers to correction functions:

  ff_settings->g_r_alpha = r_alpha_correction;
  ff_settings->g_r_beta = r_beta_correction;
  ff_settings->g_alpha_beta = alpha_beta_correction;

  // r-alpha correction:

  corr = (params_get_parameter("r_alpha_correction"))->value.s;

  if (!strcmp(corr,"block")) {

    ff_settings->r_alpha_correction = BLOCK_R_ALPHA_CORRECTION;

  } else if (!strcmp(corr,"none")) {

    ff_settings->r_alpha_correction = NO_R_ALPHA_CORRECTION;
  }

  // r-beta correction:

  corr = (params_get_parameter("r_beta_correction"))->value.s;

  if (!strcmp(corr,"block")) {

    ff_settings->r_beta_correction = BLOCK_R_BETA_CORRECTION;

  } else if (!strcmp(corr,"none")) {

    ff_settings->r_beta_correction = NO_R_BETA_CORRECTION;
  }

  // alpha-beta correction:

  corr = (params_get_parameter("alpha_beta_correction"))->value.s;

  if (!strcmp(corr,"block")) {

    ff_settings->alpha_beta_correction = BLOCK_ALPHA_BETA_CORRECTION;

  } else if (!strcmp(corr,"cosine")) {

    ff_settings->alpha_beta_correction = COSINE_ALPHA_BETA_CORRECTION;

  } else if (!strcmp(corr,"none")) {

    ff_settings->alpha_beta_correction = NO_ALPHA_BETA_CORRECTION;
  }

  // block function parameters:

  ff_settings->r_alpha_correction_r1 = (params_get_parameter("r_alpha_correction_r1"))->value.d;
  ff_settings->r_alpha_correction_r2 = (params_get_parameter("r_alpha_correction_r2"))->value.d;

  ff_settings->r_beta_correction_r1 = (params_get_parameter("r_beta_correction_r1"))->value.d;
  ff_settings->r_beta_correction_r2 = (params_get_parameter("r_beta_correction_r2"))->value.d;

  ff_settings->alpha_beta_correction_alpha1 = (params_get_parameter("alpha_beta_correction_alpha1"))->value.d;
  ff_settings->alpha_beta_correction_alpha2 = (params_get_parameter("alpha_beta_correction_alpha2"))->value.d;

  // limit memory use by restricting number of potentials kept in memory:

  ff_settings->max_stored_potentials = (params_get_parameter("ff_max_stored_potentials"))->value.i;
}



void init_ff(SETTINGS *settings) {

  ATOM_TYPING_SCHEME *scheme;
  FORCE_FIELD *ff;
  HISTOGRAM *histogram;

  scheme = settings->atom_typing_scheme;

  free_ff(settings->force_field);

  settings->force_field = alloc_ff(scheme);

  ff = settings->force_field;

  ff->settings = settings;

  init_atom_ff_list(ff->atom,scheme,ff);

  init_nonbonded_ff_matrix(ff->nonbonded,scheme,ff);

  // set pointer to ff:

  pliff_ff = ff;

  // load ami propensities:

  init_ami_ff(ff);
}



ATOM_FF* get_atom_ff(ATOM_FF *atom_ff_list,ATOM_TYPE *type) {

  ATOM_FF *atom_ff;

  if (type == NULL) {

    error_fn("get_atom_ff: type undefined");
  }

  atom_ff = atom_ff_list + type->i;

  if (atom_ff->ready == 0) {

    read_atom_ff(atom_ff);
  }

  return(atom_ff);
}



NONBONDED_FF* get_nonbonded_ff(NONBONDED_FF **nonbonded_ff_matrix,ATOM_TYPE *type1,ATOM_TYPE *type2) {

  NONBONDED_FF *nonbonded_ff,*paired_ff,*surr_ff;

  if ((type1 == NULL) || (type2 == NULL)) {

    error_fn("get_nonbonded_ff: type(s) undefined");
  }

  nonbonded_ff = nonbonded_ff_matrix[type1->i] + type2->i;
  paired_ff = (type1 == type2) ? (NULL) : nonbonded_ff_matrix[type2->i] + type1->i;

  if (nonbonded_ff->ready == 0) {

    if (read_nonbonded_ff(nonbonded_ff)) {

      surr_ff = get_surrogate_nbff(nonbonded_ff);

      copy_nonbonded_ff(surr_ff,nonbonded_ff);

      nonbonded_ff->surrogate = 3;
    }

    if (paired_ff) {

      copy_nonbonded_ff(nonbonded_ff,paired_ff);
    }
  }

  nonbonded_ff->usage_count++;

  if (paired_ff) {

    paired_ff->usage_count++;
  }

  return(nonbonded_ff);
}



static void read_atom_ff(ATOM_FF *atom_ff) {

  int flag;
  char *pli_dir,name[10],filename[1000],line[MAX_LINE_LEN],word[MAX_LINE_LEN];
  ATOM_TYPE *type;
  HISTOGRAM *histogram,*fhist;
  PLI_FILE *file;

  pli_dir = get_pli_dir();

  type = atom_ff->type;

  (type->id == -1) ? sprintf(name,"etype%d",type->element->id) : sprintf(name,"type%d",type->id);

  sprintf(filename,"%s/ff/nonbonded/%s/nonbonded.pliff",pli_dir,name);

  file = open_file(filename,"r");

  if (file == NULL) {

    error_fn("read_atom_ff: failed to open ff file '%s'",filename);
  }

  while (!end_of_file(file)) {

    if (read_line(line,MAX_LINE_LEN,file) == NULL)
      break;

    flag = sscanf(line,"%s",word);

    if (flag != EOF) {

      if (!strcmp(word,"avg_total_area")) {

	sscanf(line,"%*s %lf %lf",&(atom_ff->avg_total_area),&(atom_ff->std_total_area));

      } else if (!strcmp(word,"avg_inter_area")) {

	sscanf(line,"%*s %lf %lf",&(atom_ff->avg_inter_area),&(atom_ff->std_inter_area));

      } else if (!strcmp(word,"max_polar_area")) {

	sscanf(line,"%*s %lf",&(atom_ff->max_polar_area));

      } else if (!strcmp(word,"histogram")) {

	histogram = read_histogram(file,line);

	if (!strcmp(histogram->name,"Z")) {

	  histogram->max_peaks = 2;

	  atom_ff->Zp_histogram_id = histogram->id;
	}

	fhist = process_atom_ff_histogram(histogram,atom_ff);

	atom_ff->Zf_histogram_id = fhist->id;
      }
    }
  }

  close_file(file);

  atom_ff->ready = 1;
}



static int read_nonbonded_ff(NONBONDED_FF *nonbonded_ff) {

  int i,i1,i2,flag,id,use_ami_pts;
  char *pli_dir,name1[10],name2[10],filename[1000],line[MAX_LINE_LEN],word[MAX_LINE_LEN];
  double f_r_alpha,f_r_beta,f_alpha1_beta,f_alpha2_beta;
  ATOM_TYPE *type1,*type2;
  HISTOGRAM *histogram;
  HISTOGRAM_POINT *point;
  NONBONDED_FF *nb_ff;
  PLI_FILE *file;

  use_ami_pts = (params_get_parameter("ff_use_ami_pts"))->value.i;

  // release some memory if needed:

  if (n_stored_potentials == ff_settings->max_stored_potentials) {

    //nb_ff = get_least_used_nonbonded_ff(nonbonded_ff->force_field);

    if (last_read_ff) {

      reset_nonbonded_ff(last_read_ff);
    }
  }

  pli_dir = get_pli_dir();

  type1 = nonbonded_ff->type1;
  type2 = nonbonded_ff->type2;

  i1 = type1->i;
  i2 = type2->i;

  (type1->id == -1) ? sprintf(name1,"etype%d",type1->element->id) : sprintf(name1,"type%d",type1->id);
  (type2->id == -1) ? sprintf(name2,"etype%d",type2->element->id) : sprintf(name2,"type%d",type2->id);

  (i1 < i2) ? sprintf(filename,"%s/ff/nonbonded/%s/%s.pliff",pli_dir,name1,name2) : sprintf(filename,"%s/ff/nonbonded/%s/%s.pliff",pli_dir,name2,name1);
  file = open_file(filename,"r");

  if (file == NULL) {

    warning_fn("%s: Using surrogate file as %s could not be opened", __func__,filename);

    return(1);
  }

  while (!end_of_file(file)) {

    if (read_line(line,MAX_LINE_LEN,file) == NULL)
      break;

    flag = sscanf(line,"%s",word);

    if (flag != EOF) {

      if ((!use_ami_pts) && (!strcmp(word,"P"))) {

	sscanf(line,"%*s %lf %lf",&(nonbonded_ff->P),&(nonbonded_ff->sP));

      } else if ((use_ami_pts) && (!strcmp(word,"Pu"))) {

	sscanf(line,"%*s %lf %lf",&(nonbonded_ff->P),&(nonbonded_ff->sP));

      } else if ((!use_ami_pts) && (!strcmp(word,"lnP"))) {

	sscanf(line,"%*s %lf %lf",&(nonbonded_ff->lnP),&(nonbonded_ff->slnP));

      } else if ((use_ami_pts) && (!strcmp(word,"lnPu"))) {

	sscanf(line,"%*s %lf %lf",&(nonbonded_ff->lnP),&(nonbonded_ff->slnP));

      } else if (!strcmp(word,"Precision")) {

	sscanf(line,"%*s %s",nonbonded_ff->precision);

      } else if (!strcmp(word,"histogram")) {

	histogram = read_histogram(file,line);

	if (!strcmp(histogram->name,"R")) {

	  nonbonded_ff->R_histogram_id = histogram->id;

	} else if ((!strcmp(histogram->name,"ALPHA1")) && (ff_settings->use_alpha_potentials)) {

	  histogram->max_peaks = 2;

	  if (i1 < i2) {

	    nonbonded_ff->ALPHA1_histogram_id = histogram->id;

	  } else {

	    nonbonded_ff->ALPHA2_histogram_id = histogram->id;
	  }

	  process_nonbonded_ff_histogram(histogram,nonbonded_ff);

	} else if ((!strcmp(histogram->name,"BETA1")) && (ff_settings->use_beta_potentials)) {

	  histogram->max_peaks = 2;

	  if (i1 < i2) {

	    nonbonded_ff->BETA1_histogram_id = histogram->id;

	  } else {

	    nonbonded_ff->BETA2_histogram_id = histogram->id;
	  }

	  process_nonbonded_ff_histogram(histogram,nonbonded_ff);

	} else if ((!strcmp(histogram->name,"ALPHA2")) && (ff_settings->use_alpha_potentials)) {

	  histogram->max_peaks = 2;

	  if (i1 < i2) {

	    nonbonded_ff->ALPHA2_histogram_id = histogram->id;

	  } else {

	    nonbonded_ff->ALPHA1_histogram_id = histogram->id;
	  }

	  process_nonbonded_ff_histogram(histogram,nonbonded_ff);

	} else if ((!strcmp(histogram->name,"BETA2")) && (ff_settings->use_beta_potentials)) {

	  histogram->max_peaks = 2;

	  if (i1 < i2) {

	    nonbonded_ff->BETA2_histogram_id = histogram->id;

	  } else {

	    nonbonded_ff->BETA1_histogram_id = histogram->id;
	  }

	  process_nonbonded_ff_histogram(histogram,nonbonded_ff);
	}
      }
    }
  }

  close_file(file);
  
  if (nonbonded_ff->P < -0.5) {

    error_fn("read_nonbonded_ff: type propensity could not be read from %s",filename);
  }

  process_R_potential(nonbonded_ff);

  f_r_alpha = calc_2d_correction(nonbonded_ff->R_histogram_id,ff_settings->g_r_alpha);
  f_r_beta = calc_2d_correction(nonbonded_ff->R_histogram_id,ff_settings->g_r_beta);

  f_alpha1_beta = calc_2d_correction(nonbonded_ff->ALPHA1_histogram_id,ff_settings->g_alpha_beta);
  f_alpha2_beta = calc_2d_correction(nonbonded_ff->ALPHA2_histogram_id,ff_settings->g_alpha_beta);

  apply_2d_correction(nonbonded_ff->ALPHA1_histogram_id,f_r_alpha);
  apply_2d_correction(nonbonded_ff->ALPHA2_histogram_id,f_r_alpha);

  apply_2d_correction(nonbonded_ff->BETA1_histogram_id,f_r_beta*f_alpha1_beta);
  apply_2d_correction(nonbonded_ff->BETA2_histogram_id,f_r_beta*f_alpha2_beta);

  nonbonded_ff->ready = 1;

  last_read_ff = nonbonded_ff;

  n_stored_potentials++;

  return(0);
}



void copy_nonbonded_ff(NONBONDED_FF *nonbonded_ff1,NONBONDED_FF *nonbonded_ff2) {

  if (nonbonded_ff1 == nonbonded_ff2) {

    return;
  }

  nonbonded_ff2->surrogate = nonbonded_ff1->surrogate;
  nonbonded_ff2->ready = 1;

  nonbonded_ff2->P = nonbonded_ff1->P;
  nonbonded_ff2->sP = nonbonded_ff1->sP;
  nonbonded_ff2->lnP = nonbonded_ff1->lnP;
  nonbonded_ff2->slnP = nonbonded_ff1->slnP;

  nonbonded_ff2->Dvdw = nonbonded_ff1->Dvdw;

  nonbonded_ff2->Do = nonbonded_ff1->Do;
  nonbonded_ff2->Eo = nonbonded_ff1->Eo;

  nonbonded_ff2->Dc = nonbonded_ff1->Dc;

  nonbonded_ff2->k1 = nonbonded_ff1->k1;
  nonbonded_ff2->k2 = nonbonded_ff1->k2;

  nonbonded_ff2->LJ_n = nonbonded_ff1->LJ_n;

  nonbonded_ff2->LR_A = nonbonded_ff1->LR_A;
  nonbonded_ff2->LR_B = nonbonded_ff1->LR_B;
  nonbonded_ff2->LR_D = nonbonded_ff1->LR_D;

  nonbonded_ff2->opt_Z = nonbonded_ff1->opt_Z;
  nonbonded_ff2->opt_sY = nonbonded_ff1->opt_sY;

  nonbonded_ff2->R_histogram_id = nonbonded_ff1->R_histogram_id;
  nonbonded_ff2->ALPHA1_histogram_id = nonbonded_ff1->ALPHA2_histogram_id;
  nonbonded_ff2->BETA1_histogram_id = nonbonded_ff1->BETA2_histogram_id;
  nonbonded_ff2->ALPHA2_histogram_id = nonbonded_ff1->ALPHA1_histogram_id;
  nonbonded_ff2->BETA2_histogram_id = nonbonded_ff1->BETA1_histogram_id;

  nonbonded_ff2->force_field = nonbonded_ff1->force_field;
}



double r_alpha_correction(double r) {

  double r1,r2;

  if (ff_settings->r_alpha_correction == BLOCK_R_ALPHA_CORRECTION) {

    r1 = ff_settings->r_alpha_correction_r1;
    r2 = ff_settings->r_alpha_correction_r2;

    if (r < r1) {

      return(1.0);

    } else if ((r >= r1) && (r < r2)) {

      return(1.0 - (r - r1)/(r2 - r1));

    } else {

      return(0.0);
    }
  }

  return(1.0);
}



double r_beta_correction(double r) {

  double r1,r2;

  if (ff_settings->r_beta_correction == BLOCK_R_BETA_CORRECTION) {

    r1 = ff_settings->r_beta_correction_r1;
    r2 = ff_settings->r_beta_correction_r2;

    if (r < r1) {

      return(1.0);

    } else if ((r >= r1) && (r < r2)) {

      return(1.0 - (r - r1)/(r2 - r1));

    } else {

      return(0.0);
    }
  }

  return(1.0);
}



double alpha_beta_correction(double alpha) {

  double alpha1,alpha2,salpha;

  if (ff_settings->alpha_beta_correction == COSINE_ALPHA_BETA_CORRECTION) {

    return(0.5*(1.0 - cos((PI/90.0)*alpha)));

  } else if (ff_settings->alpha_beta_correction == BLOCK_ALPHA_BETA_CORRECTION) {

    alpha1 = ff_settings->alpha_beta_correction_alpha1;
    alpha2 = ff_settings->alpha_beta_correction_alpha2;

    salpha = (alpha < 90.0) ? alpha : 180.0 - alpha;

    if (salpha < alpha1) {

      return(0.0);

    } else if ((salpha >= alpha1) && (salpha < alpha2)) {

      return((salpha - alpha1)/(alpha2 - alpha1));

    } else {

      return(1.0);
    }
  }

  return(1.0);
}



static void set_surrogate_R_potential(NONBONDED_FF *nbff) {

  int i;
  HISTOGRAM *histogram,*shistogram;
  HISTOGRAM_POINT *point;
  NONBONDED_FF *snbff;

  snbff = get_surrogate_nbff(nbff);

  copy_R_potential_params(nbff,snbff,DELTA_D_SLOPE*(nbff->Dvdw - snbff->Dvdw));

  // sort out histogram:

  if (nbff->surrogate == 1) {

    // histogram was present, but looked dodgy for some reason:

    warning_fn("%s: replacing R potential for %s (%d) - %s (%d) with surrogate potential",__func__,
	       nbff->type1->name,nbff->type1->id,nbff->type2->name,nbff->type2->id);

    histogram = get_histogram(nbff->R_histogram_id);

    for (i=0,point=histogram->points;i<histogram->n_points;i++,point++) {

      point->Y = calc_fitted_R_potential(point->X + nbff->Dvdw,nbff);
    }

  } else if (nbff->surrogate == 2) {

    // no histogram at all, take a copy of the surrogate histogram:

    warning_fn("%s: no R potential found for %s (%d) - %s (%d) using surrogate potential",__func__,
	       nbff->type1->name,nbff->type1->id,nbff->type2->name,nbff->type2->id);

    shistogram = get_histogram(snbff->R_histogram_id);

    histogram = clone_histogram(shistogram,1);

    nbff->R_histogram_id = histogram->id;

  } else {

    error_fn("%s: unknown surrogate type (%d)",__func__,nbff->surrogate);
  }
}



static void copy_R_potential_params(NONBONDED_FF *nbff1,NONBONDED_FF *nbff2,double dD) {

  // copy some params directly:

  strcpy(nbff1->shape,nbff2->shape);

  nbff1->Po = nbff2->Po;
  nbff1->Eo = nbff2->Eo;

  nbff1->k1 = nbff2->k1;
  nbff1->k2 = nbff2->k2;

  nbff1->LJ_n = nbff2->LJ_n;

  nbff1->LR_A = nbff2->LR_A;
  nbff1->LR_B = nbff2->LR_B;
  nbff1->LR_n = nbff2->LR_n;

  // adjust other parameters according to the respective vdw radii:

  nbff1->Do = nbff2->Do + dD;
  nbff1->Dc = (nbff1->Do)/(pow(2,1.0/nbff1->LJ_n));

  nbff1->LR_D = nbff2->LR_D + dD;
}



static NONBONDED_FF* get_surrogate_nbff(NONBONDED_FF *nbff) {

  SETTINGS *settings;
  FORCE_FIELD *ff;
  unsigned int flags1,flags2;
  ATOM_TYPE *type1,*type2;
  NONBONDED_FF *snbff;

  settings = get_settings();
  ff = settings->force_field;

  flags1 = nbff->type1->flags;
  flags2 = nbff->type2->flags;

  if (hbond_flags_match(flags1,flags2,0)) {

    type1 = get_atom_type("Amide >NH",settings->atom_typing_scheme);
    type2 = get_atom_type("Threonine -OH",settings->atom_typing_scheme);

  } else if ((flags1 & HBOND_ACCEPTOR_ATOM_TYPE) && (flags2 & HBOND_ACCEPTOR_ATOM_TYPE)) {

    type1 = type2 = get_atom_type("=O",settings->atom_typing_scheme);

  } else if ((flags1 & HBOND_DONOR_ATOM_TYPE) && (flags2 & HBOND_DONOR_ATOM_TYPE)) {

    type1 = type2 = get_atom_type("=NH-",settings->atom_typing_scheme);

  } else {

    type1 = type2 = get_atom_type("-CH3",settings->atom_typing_scheme);
  }

  snbff = get_nonbonded_ff(ff->nonbonded,type1,type2);

  if ((snbff == NULL) || (snbff->R_histogram_id == -1)) {

    error_fn("%s: surrogate R distribution undefined",__func__);
  }

  return(snbff);
}



static HISTOGRAM* process_atom_ff_histogram(HISTOGRAM *histogram,ATOM_FF *atom_ff) {

  int sumN,Z,id;
  double A,sumA,Zo,Ye,minY;
  HISTOGRAM *fhist;
  HISTOGRAM_POINT *point,*min_point,*points;

  // normalise and smoothe histogram:

  histogram->normalisation_method = ABSOLUTE_NORMALISATION;

  histogram->start_s = 0.2;
  histogram->end_s = 2.0;
  histogram->step_s = 0.2;

  process_histogram(histogram);

  // clone histogram to store Z frequency distribution in separate histogram:

  id = histogram->id;

  fhist = add_histogram("Zf",histogram->n_points);

  histogram = get_histogram(id);

  id = fhist->id;

  points = fhist->points;

  memcpy(fhist,histogram,sizeof(HISTOGRAM));

  fhist->points = points;

  memcpy(fhist->points,histogram->points,(histogram->n_points)*sizeof(HISTOGRAM_POINT));

  fhist->id = id;

  // then scale histogram relative to expected fractions:

  A = atom_ff->max_polar_area;

  sumN = histogram->sumN;
  sumA = histogram->sumY;

  if (sumA < 1.0E-10) {

    error_fn("process_atom_ff_histogram: sumA is (nearly) zero");
  }

  Zo = (A*((double) sumN))/sumA;

  min_point = NULL;

  for (Z=0,point=histogram->points;Z<histogram->n_points;Z++,point++) {

    Ye = ((pow(Zo,(double) Z))*exp(-(Zo)))/factorial(Z);

    point->Y /= Ye;

    if ((point->Y > 1.0E-10) && ((min_point == NULL) || (point->Y < min_point->Y))) {

      min_point = point;
    }
  }

  if (min_point == NULL) {

    error_fn("process_atom_ff_histogram: min_point undefined");
  }

  // log the propensities:

  minY = 0.5*(min_point->Y);

  for (Z=0,point=histogram->points;Z<histogram->n_points;Z++,point++) {

    point->Y = (point->Y > 1.0E-10) ? log(point->Y) : log(minY);
  }

  histogram->extrapolation_method = EXTRAPOLATE_TO_NEAREST;

  set_histogram_extrema(histogram);

  return(fhist);
}


static void process_atom_ff_histogram_old(HISTOGRAM *histogram,ATOM_FF *atom_ff) {

  int Z,Zo,n_points,i,hist_n_points,hist_id;
  long int sumN;
  double sumA,A,maxA,C,X,Zop,x,fe,fo,minY;
  HISTOGRAM *A_hist;
  HISTOGRAM_POINT *A_point,*Z_point,*min_point;

  // normalise and smoothe histogram:

  histogram->normalisation_method = ABSOLUTE_NORMALISATION;

  histogram->start_s = 0.2;
  histogram->end_s = 2.0;
  histogram->step_s = 0.2;

  process_histogram(histogram);

  // then scale histogram relative to expected fractions:

  maxA = atom_ff->max_polar_area;

  sumN = histogram->sumN;
  sumA = histogram->sumY;

  if (sumA < 1.0E-10) {

    error_fn("process_atom_ff_histogram: sumA is (nearly) zero");
  }

  C = ((double) sumN)/sumA;

  n_points = ceil(maxA/0.01);

  min_point = NULL;

  // first add the "matrix" histograms, because otherwise memory can get screwed up:

  hist_id = histogram->id;
  hist_n_points = histogram->n_points;

  for (Z=0;Z<hist_n_points;Z++) {

    A_hist = add_histogram("AZ",n_points);

    atom_ff->A_histogram_ids[Z] = A_hist->id;
  }

  histogram = get_histogram(hist_id);

  for (Z=0;Z<histogram->n_points;Z++) {

    A_hist = get_histogram(atom_ff->A_histogram_ids[Z]);

    A_hist->startX = 0.0;
    A_hist->endX = maxA;
    A_hist->stepX = maxA/((double) A_hist->n_points);

    for (i=0,A_point=A_hist->points;i<A_hist->n_points;i++,A_point++) {

      A = (((double) i) + 0.5)*(A_hist->stepX);
      X = A/maxA;

      fe = ((pow(A*C,(double) Z))*exp(-(A*C)))/factorial(Z);

      fo = 0.0;

      for (Zo=0,Z_point=histogram->points;Zo<histogram->n_points;Zo++,Z_point++) {

	Zop = ((double) Zo)*X;

	x = fabs(Zop-((double) Z));

	// maybe need something smoother than this:

	if (x < 1.0) {

	  fo += (1-x)*(Z_point->Y);
	}
      }

      A_point->X = A;
      A_point->Y = fo/fe;

      if ((A_point->Y > 1.0E-10) && ((min_point == NULL) || (A_point->Y < min_point->Y))) {

	min_point = A_point;
      }
    }
  }


  if (min_point == NULL) {

    error_fn("process_atom_ff_histogram: min_point undefined");
  }

  minY = 0.5*(min_point->Y);

  for (Z=0;Z<histogram->n_points;Z++) {


    A_hist = get_histogram(atom_ff->A_histogram_ids[Z]);

    for (i=0,A_point=A_hist->points;i<A_hist->n_points;i++,A_point++) {

      A_point->Y = (A_point->Y > 1.0E-10) ? log(A_point->Y) : log(minY);

    }


    A_hist->extrapolation_method = EXTRAPOLATE_TO_NEAREST;

    set_histogram_extrema(A_hist);
  }
}



static void process_nonbonded_ff_histogram(HISTOGRAM *histogram,NONBONDED_FF *nonbonded_ff) {

  ATOM_GEOMETRY *geom1,*geom2;

  geom1 = nonbonded_ff->type1->geometry;
  geom2 = nonbonded_ff->type2->geometry;

  if (nonbonded_ff->R_histogram_id == histogram->id) {

    histogram->limits = geom1->r_limits;

    if ((params_get_parameter("pliff_r_correction"))->value.i) {

      histogram->geometric_correction = PLIFF_R_CORRECTION;
    }

  } else if (nonbonded_ff->ALPHA1_histogram_id == histogram->id) {

    histogram->limits = geom1->alpha_limits;

    histogram->geometric_correction = CONICAL_CORRECTION;

  } else if (nonbonded_ff->BETA1_histogram_id == histogram->id) {

    histogram->limits = geom1->beta_limits;

  } else if (nonbonded_ff->ALPHA2_histogram_id == histogram->id) {

    histogram->limits = geom2->alpha_limits;

    histogram->geometric_correction = CONICAL_CORRECTION;

  } else if (nonbonded_ff->BETA2_histogram_id == histogram->id) {

    histogram->limits = geom2->beta_limits;
  }

  histogram->extrapolation_method = EXTRAPOLATE_TO_SMALLEST;
  histogram->normalisation_method = RELATIVE_NORMALISATION;
  histogram->log_scale = 1;

  process_histogram(histogram);

  set_histogram_extrema(histogram);
}



static int process_R_potential(NONBONDED_FF *nbff) {

  char *style;

  if (nbff->R_histogram_id != -1) {

    init_R_potential(nbff);

    style = (params_get_parameter("pliff_potentials"))->value.s;

    if (!strcmp(style,"smoothed")) {
      
      if (smoothe_R_potential(nbff)) {

	nbff->surrogate = 1;
      }

    } else {
            
      if (fit_R_potential(nbff)) {
	
	nbff->surrogate = 1;
      }
    }

  } else {

    nbff->surrogate = 2;
  }

  if (nbff->surrogate) {

    set_surrogate_R_potential(nbff);
  }

  finish_R_potential(nbff);

  return(0);
}



static void init_R_potential(NONBONDED_FF *nbff) {

  HISTOGRAM *histogram;
  ATOM_GEOMETRY *geom1;

  histogram = get_histogram(nbff->R_histogram_id);
  
  geom1 = nbff->type1->geometry;
  
  histogram->extrapolation_method = EXTRAPOLATE_TO_SMALLEST;
  histogram->normalisation_method = RELATIVE_NORMALISATION;
  histogram->log_scale = 1;
  histogram->limits = geom1->r_limits;
  histogram->max_peaks = 1;
  
  if ((params_get_parameter("pliff_r_correction"))->value.i) {
    
    histogram->geometric_correction = PLIFF_R_CORRECTION;
  }
}



static void finish_R_potential(NONBONDED_FF *nbff) {

  HISTOGRAM *histogram;

  histogram = get_histogram(nbff->R_histogram_id);

  set_histogram_extrema(histogram);

  calc_R_potential_fit_quality(histogram,nbff);
}



static int smoothe_R_potential(NONBONDED_FF *nbff) {

  int i,j,flag;
  double D,slope,intercept,Xo,Yo;
  ATOM_TYPE *type1,*type2;
  HISTOGRAM *histogram;
  HISTOGRAM_POINT *point,*point1,*point2,*point3;

  histogram = get_histogram(nbff->R_histogram_id);

  // old style smoothing:

  process_histogram(histogram);
  
  // find two points to calculate clash distance from:

  i = 0;
  point1 = histogram->points;

  if (point1->Y > 0.0) {

    warning_fn("%s: potential starts positive for %s - %s",__func__,nbff->type1->name,nbff->type2->name);

    point2 = point1 + 1;

  } else {

    point2 = histogram->points;

    do {

      i++;
      point2++;

    } while ((i < histogram->n_points) && (point2->Y < 0.0));

    if (i == histogram->n_points) {

      error_fn("%s: corrupt potential (1)",__func__);
    }
  
    point1 = point2 - 1;
  }

  // intrapolate (or extrapolate if potential starts positive) to get clash distance:

  flag = solve_line(point1->X,point1->Y,point2->X,point2->Y,&slope,&intercept);

  if ((flag) || (slope < 0.0)) {

    error_fn("%s: corrupt potential (2)",__func__);
  }

  nbff->Dc = (fabs(slope) < 1.0E-30) ? nbff->Dvdw : nbff->Dvdw - (intercept/slope);

  // now find three points around the optimum:

  point3 = point2;
  point2 = point1;

  do {

    i++;
    point2++;
    point3++;

  } while ((i < histogram->n_points) && (point3->Y > point2->Y));

  if (i == histogram->n_points) {

    error_fn("%s: corrupt potential (3)",__func__);
  }

  point1 = point2 - 1;
 
  // fit parabola to three points to find optimum distance:

  flag = calc_parabola_vertex(point1->X,point1->Y,point2->X,point2->Y,point3->X,point3->Y,&Xo,&Yo);

  if (flag) {

    error_fn("%s: unable to derive optimum distance from R histogram",__func__);
  }

  nbff->Do = nbff->Dvdw + Xo;
  nbff->Eo = Yo;
  nbff->Po = exp(Yo);

  nbff->LJ_n = log(2.0)/log(nbff->Do/nbff->Dc);

  if ((params_get_parameter("pliff_lj_clash"))->value.i) {

    for (j=0,point=histogram->points;j<i;j++,point++) {

      if (point->X + nbff->Dvdw < nbff->Dc) {

	point->Y = -lennard_jones_energy(point->X + nbff->Dvdw,nbff->Do,nbff->Eo,nbff->LJ_n);
      }
    }
  }

  return(0);
}



static void calc_R_potential_fit_quality(HISTOGRAM *histogram,NONBONDED_FF *nbff) {

  int i,n1,n2;
  double D,y1[1000],sy1[1000],yfit1[100],y2[1000],sy2[1000],yfit2[1000];
  HISTOGRAM_POINT *point;

  n1 = n2 = 0;

  for (i=0,point=histogram->points;i<histogram->n_points;i++,point++) {

    D = point->X + nbff->Dvdw;

    if ((D > nbff->Dc) && (D < nbff->LR_D)) {

      y1[n1] = point->Yraw;
      sy1[n1] = point->sYraw;
      yfit1[n1++] = exp(point->Y);

    } else if (D > nbff->LR_D) {

      y2[n2] = point->Yraw;
      sy2[n2] = point->sYraw;
      yfit2[n2++] = exp(point->Y);
    }
  }

  nbff->Z1 = normalised_fit_error(y1,yfit1,sy1,n1);
  nbff->Z2 = normalised_fit_error(y2,yfit2,sy2,n2);

  nbff->sY1 = weighted_fit_error(y1,yfit1,sy1,n1);
  nbff->sY2 = weighted_fit_error(y2,yfit2,sy2,n2);
}



static void normalise_R_histogram(HISTOGRAM *histogram,NONBONDED_FF *nonbonded_ff) {

  int i,n,pliff_r_correction;
  double stepX,intY,intF,C,F,D,Dmax;
  HISTOGRAM_POINT *point;

  Dmax = nonbonded_ff->Dvdw + 2.0*((get_atom_type("H2O",(get_settings())->atom_typing_scheme))->united_atom->vdw_radius);

  intY = 0.0;
  intF = 0.0;

  pliff_r_correction = ((params_get_parameter("pliff_r_correction"))->value.i) ? 1 : 0;

  stepX = histogram->stepX;

  n = 0;

  for (i=0,point=histogram->points;i<histogram->n_points;i++,point++) {

    D = point->X + nonbonded_ff->Dvdw;

    if (D < Dmax) {

      if ((histogram->limits == NULL) || ((point->X > histogram->limits[0]) && (point->X < histogram->limits[1]))) {

	intY += (point->Y)*stepX;
	
	intF += (pliff_r_correction) ? ((sqr(D) - (pow(D,3.0)/Dmax))*stepX) : stepX;
      }

      n++;
    }
  }

  histogram->n_points = n;

  if (intF < 1.0E-10) {

    error_fn("normalise_histogram: intF (nearly) zero (%8e) for histogram %s",intF,histogram->name);
  }

  C = (histogram->normalisation_method == RELATIVE_NORMALISATION) ? (intF/intY) : (1.0/intY);

  for (i=0,point=histogram->points;i<histogram->n_points;i++,point++) {

    D = point->X + nonbonded_ff->Dvdw;

    F = (pliff_r_correction) ? (sqr(D) - (pow(D,3.0)/Dmax)) : 1.0;

    point->Y *= (C/F);
    point->sY *= (C/F);

    point->Yraw *= (C/F);
    point->sYraw *= (C/F);
  }
}



static int fit_R_potential(NONBONDED_FF *nbff) {

  int i;
  double D;
  char *style;
  HISTOGRAM *histogram;
  HISTOGRAM_POINT *point;

  histogram = get_histogram(nbff->R_histogram_id);

  // re-calculate error estimates:

  recalc_R_histogram_errors(histogram);

  // normalise and log histogram:

  normalise_R_histogram(histogram,nbff);

  // generate rough initial fit for entire potential:

  if (fit_R_potential_full(histogram,nbff)) {

    warning_fn("%s: unable to fit full R potential for %s (%d) - %s (%d)",__func__,
	       nbff->type1->name,nbff->type1->id,
	       nbff->type2->name,nbff->type2->id);

    return(1);
  }

  // refine the position and height of the peak:

  refine_R_optimum(histogram,nbff);

  // fit short-range potential:

  fit_short_range_potential(histogram,nbff);

  // fit long-range potential:

  if (fit_long_range_potential(histogram,nbff)) {

    warning_fn("%s: unable to fit long-range R potential for %s (%d) - %s (%d)",__func__,
	       nbff->type1->name,nbff->type1->id,
	       nbff->type2->name,nbff->type2->id);

    return(1);
  }

  apply_R_potential_softness(nbff);

  // calculate potential:

  style = (params_get_parameter("pliff_potentials"))->value.s;

  for (i=0,point=histogram->points;i<histogram->n_points;i++,point++) {
    
    D = point->X + nbff->Dvdw;
    
    if (!strcmp(style,"LJ")) {

      if (D < nbff->Do) {

	point->Y = -lennard_jones_energy(D,nbff->Do,nbff->Eo,nbff->LJ_n);

      } else {

	point->Y = -lennard_jones_energy(D,nbff->Do,nbff->Eo,nbff->LR_n);
      }

    } else if (!strcmp(style,"fitted")) {

      point->Y = calc_fitted_R_potential(D,nbff);

    } else {

      error_fn("%s: unknown potential style '%s'",__func__,style);
    }
  }
  
  return(0);
}



static void apply_R_potential_softness(NONBONDED_FF *nbff) {
  
  double Xo,Yo,X1,x1,x2;
  
  Xo = nbff->Do - nbff->Dvdw;
  Yo = nbff->Eo;

  // find point where Y=0:

  if (!calc_parabola_x_intercepts(-(nbff->k1),2.0*(nbff->k1)*Xo,Yo-(nbff->k1)*sqr(Xo),&x1,&x2)) {

    error_fn("%s: unexpected error occurred",__func__);
  }

  // apply softness and calculate LJ order:

  X1 = fmin(x1,x2) - (params_get_parameter("pliff_softness"))->value.d;

  nbff->k1 = Yo/sqr(X1-Xo);

  nbff->Dc = X1 + nbff->Dvdw;

  nbff->LJ_n = log(2.0)/log(nbff->Do/nbff->Dc);

}



static double calc_fitted_R_potential(double D,NONBONDED_FF *nbff) {

  if (D < nbff->Do) {
    
    return((nbff->Eo) - (nbff->k1)*sqr((D-nbff->Do)));
  }

  if (D < nbff->LR_D) {
    
    return((nbff->Eo) - (nbff->k2)*sqr((D-nbff->Do)));
  }
    
  return(log(nbff->LR_B) - (nbff->LR_A)*(D-(nbff->LR_D)));
}



static int fit_R_potential_full(HISTOGRAM *histogram,NONBONDED_FF *nbff) {

  int i,j,maxi,flag;
  double *s,*int1,int_tot,k1,k2,y,w,sumw;
  HISTOGRAM_POINT *point,*point1,*point2;

  // if there is lots of data, we may be able to get the optimum straight from the raw points:

  if (!get_R_optimum_from_points(histogram,nbff)) {

    return(0);
  }

  int1 = (double*) calloc(histogram->n_points,sizeof(double));
  s = (double*) calloc(histogram->n_points,sizeof(double));

  if ((int1 == NULL) || (s == NULL)) {

    error_fn("%s: out of memory allocating int1",__func__);
  }

  // get total integral, and left-side integral for each point:

  int_tot = 0.0;
 
  for (i=0,point=histogram->points;i<histogram->n_points;i++,point++) {

    int_tot += (point->Y)*(histogram->stepX);

    int1[i] = int_tot - 0.5*(point->Y)*(histogram->stepX);
  }

  maxi = 0;

  for (i=0,point=histogram->points;i<histogram->n_points;i++,point++) {

    s[i] = 1.0E30;

    // attempt fit if propensity is above 1.0

    if (point->Y > 1.0) {

      // calculate k1 and k2 from integral and Y value:

      k1 = 0.25*PI*sqr((point->Y)/int1[i]);

      k2 = (point->Y)/(int_tot - int1[i]);

      if ((k1 > MIN_K1) && (k2 < MAX_LR_A)) {

	// calculate fit error:

	s[i] = 0.0;
	sumw = 0.0;
	
	for (j=0,point1=histogram->points;j<histogram->n_points;j++,point1++) {
	
	  y = (j <= i) ? (point->Y)*exp(-k1*sqr((point1->X)-(point->X))) : (point->Y)*exp(-k2*((point1->X)-(point->X)));

	  w = 1.0/sqr(point->sY);
	  
	  s[i] += w*sqr(y - (point1->Y));
	
	  sumw += w;
	}
	
	s[i] = sqrt(s[i]/sumw);

	// keep track of best fit:

	if (s[i] < s[maxi]) {
	  
	  maxi = i;
	}
      }
    }
  }

  if (maxi == 0) {

    warning_fn("%s: no meaningful full R potential fit for %s (%d) - %s (%d)",__func__,
	       nbff->type1->name,nbff->type1->id,
	       nbff->type2->name,nbff->type2->id);

    return(1);
  }

  point = histogram->points + maxi;

  // check if the optimum is ambiguous:

  flag = 0;

  for (i=0,point1=histogram->points;i<histogram->n_points;i++,point1++) {

    if ((fabs((point1->X)-(point->X)) > MIN_AMBIGUOUS_PEAK_DIST) && ((s[maxi]/s[i]) > MIN_AMBIGUOUS_PEAK_RATIO)) {

      warning_fn("%s: ambiguous full R potential fit for %s (%d) - %s (%d)",__func__,
		 nbff->type1->name,nbff->type1->id,
		 nbff->type2->name,nbff->type2->id);

      return(1);
 
    } else if ((fabs((point1->X)-(point->X)) > WARNING_AMBIGUOUS_PEAK_DIST) && ((s[maxi]/s[i]) > WARNING_AMBIGUOUS_PEAK_RATIO)) {

      flag = 1;
    }
  }

  if (flag) {

    warning_fn("%s: potentially ambiguous full R potential fit for %s (%d) - %s (%d)",__func__,
	       nbff->type1->name,nbff->type1->id,
	       nbff->type2->name,nbff->type2->id);
  }

  // assign potential parameters:

  nbff->Do = point->X + nbff->Dvdw;
  nbff->Po = point->Y;
  nbff->Eo = log(point->Y);

  free(int1);
  free(s);

  return(0);
}



static int get_R_optimum_from_points(HISTOGRAM *histogram,NONBONDED_FF *nbff) {

  int i,n,maxi,fluct;
  double x[1000],y[1000],sy[1000],dy[1000],yfit[1000],c[3],sumsy,avgsy;

  histogram2xys(histogram,-10.0,10.0,1.0,999.9,x,y,sy,&n);

  // find maximum and average error:

  maxi = 0;
  
  sumsy = 0.0;

  for (i=0;i<n;i++) {

    if (y[i] > y[maxi]) {

      maxi = i;
    }

    sumsy += sy[i];
  }

  avgsy = sumsy/((double) n);

  if ((maxi < 2) || (maxi > n-3) || (avgsy > 0.05)) {

    return(1);
  }

  // check if potential rises and falls continuously:

  for (i=1;i<n;i++) {

    dy[i] = y[i] - y[i-1];

    if (i <= maxi) {

      if (dy[i] < 0.0) {

	return(2);
      }

    } else {

      if (dy[i] > 0.0) {

	return(2);
      }
    }
  }

  // check fluctuations in first derivative:

  fluct = 0;

  for (i=2;i<n-1;i++) {

    if (((dy[i] > dy[i-1]) && (dy[i] > dy[i+1])) || ((dy[i] < dy[i-1]) && (dy[i] < dy[i+1]))) {

      fluct++;
    }
  }

  if (fluct > 5) {

    return(3);
  }

  // set optimum:

  nbff->Do = x[maxi] + nbff->Dvdw;
  nbff->Po = y[maxi];
  nbff->Eo = log(nbff->Po);

  return(0);
}




static int refine_R_optimum(HISTOGRAM *histogram,NONBONDED_FF *nbff) {

  int i,n;
  double Xo,Yo,X,Y,x[1000],y[1000],sy[1000],yfit[1000],c[3],Z,avg_sy,k1,k2,fk;

  Xo = nbff->Do - nbff->Dvdw;
  Yo = log(nbff->Po);

  histogram2xys(histogram,Xo-1.5*MAX_REFINE_SHIFT_X,Xo+1.5*MAX_REFINE_SHIFT_X,1.001,999.99,x,y,sy,&n);

  avg_sy = 0.0;

  for (i=0;i<n;i++) {

    sy[i] /= y[i];
    y[i] = log(y[i]);

    avg_sy += sy[i];
  }

  avg_sy /= (double) n;

  nbff->opt_Z = nbff->opt_sY = 1.0E30;

  for (X=Xo-MAX_REFINE_SHIFT_X;X<Xo+MAX_REFINE_SHIFT_X+0.00001;X+=REFINE_DX) {

    for (Y=Yo-(MAX_REFINE_SHIFT_Y*avg_sy);Y<Yo+(MAX_REFINE_SHIFT_Y*avg_sy)+0.00001;Y+=REFINE_DY*Yo) {
     
      if (!weighted_asymmetric_parabola_fit(x,y,sy,n,X,Y,c,yfit)) {
	
	k1 = -c[0];
	k2 = -c[1];

	// check here against the k1 value that would be obtained for the
	// short-range potential:

	fk = (calc_short_range_k(histogram,X,exp(Y))/k2);

	Z = normalised_fit_error(y,yfit,sy,n);

	if ((k1 > MIN_K1) && (k2 > MIN_K2) && (fk > MIN_K_RATIO) && (fk < MAX_K_RATIO)) {

	  if (Z < nbff->opt_Z) {
  
	    //printf("Z %10.4lf %10.4lf %10.4lf %10.4lf %10.4lf %10.4lf %10.4lf\n",Z,k1,calc_short_range_k(histogram,X,exp(Y)),k2,fk,X + nbff->Dvdw,exp(Y));

	    nbff->opt_Z = Z;
	    
	    nbff->opt_sY = weighted_fit_error(y,yfit,sy,n);
	    
	    nbff->Do = X + nbff->Dvdw;
	    nbff->Eo = Y;
	    
	    nbff->k1 = k1;
	    nbff->k2 = k2;
	  }
	}
      }
    }
  }
  
  nbff->Po = exp(nbff->Eo);
}



static int fit_short_range_potential(HISTOGRAM *histogram,NONBONDED_FF *nbff) {
  
  double Xo,Yo,x1,x2;
  
  Xo = nbff->Do - nbff->Dvdw;
  Yo = nbff->Po;
  
  nbff->k1 = calc_short_range_k(histogram,Xo,Yo);

  // find point where Y=0:

  if (!calc_parabola_x_intercepts(-(nbff->k1),2.0*(nbff->k1)*Xo,log(Yo)-(nbff->k1)*sqr(Xo),&x1,&x2)) {

    error_fn("%s: unexpected error occurred",__func__);
  }

  // calculateDc and  LJ order:

  nbff->Dc = fmin(x1,x2) + nbff->Dvdw;

  nbff->LJ_n = log(2.0)/log(nbff->Do/nbff->Dc);
}



static double calc_short_range_k(HISTOGRAM *histogram,double Xo,double Yo) {
  
  double int1;
  HISTOGRAM_POINT *point;

  // integrate left-hand side of peak:
  
  int1 = 0.0;
  
  for (point=histogram->points;(point->X)<Xo;point++) {
 
    int1 += (point->Y)*(histogram->stepX);
  }
  
  point--;
  
  // add final bit between last point and top of the peak:
  
  int1 -= 0.5*(point->Y)*(histogram->stepX);
  int1 += 0.5*(Yo + (point->Y))*(Xo - (point->X));

  return((PI/4.0)*sqr(Yo/int1));
}



static int fit_long_range_potential(HISTOGRAM *histogram,NONBONDED_FF *nbff) {

  int i,n;
  double Xo,Yo,X1,Y1,X2,Y2,x[1000],y[1000],sy[1000],w[1000],yfit[1000],c[3],sumx2,sumx,sumy,sumw,a,b,k2,Z,Zmin,x1,x2,n1,n2,fk;
  HISTOGRAM *hist;
  HISTOGRAM_POINT *point;

  Xo = nbff->Do - nbff->Dvdw;
  Yo = log(nbff->Po);

  hist = clone_histogram(histogram,0);
  
  for (i=0,point=hist->points;i<hist->n_points;i++,point++) {

    if (point->N) {

      point->sY = point->sY/point->Y;
      point->Y = log(point->Y);
    }
  }

  // find optimum distance to change from parabola to straight line:

  Zmin = 1.0E10;

  for (X1=Xo+0.2;X1<Xo+2.0;X1+=0.01) {

    // fit parabola to first section after optimum:

    histogram2xys(hist,Xo,X1,-999.9,999.99,x,y,sy,&n);

    sumx2 = sumy = 0.0;

    for (i=0;i<n;i++) {

      w[i] = 1.0/sqr(sy[i]);

      sumx2 += w[i]*sqr(x[i]-Xo);
      sumy += w[i]*(y[i] - Yo);
    }

    k2 = -(sumy/sumx2);

    fk = (nbff->k1/k2);

    // if outside range, revert to previous (valid) k2 value:

    if ((k2 < MIN_K2) || (fk < MIN_K_RATIO) || (fk > MAX_K_RATIO)) {

      k2 = nbff->k2;

      fk = (nbff->k1/k2);
    }

    if ((k2 > MIN_K2) && (fk > MIN_K_RATIO) && (fk < MAX_K_RATIO)) {

      Y1 = (nbff->Po)*exp(-k2*sqr(X1-Xo));
      
      // fit straight line through remainder:
      
      histogram2xys(histogram,X1+0.0001,999.9,0.0,999.99,x,y,sy,&n);
      
      sumy = 0.5*(x[0]-X1)*(y[0]+Y1) + 0.5*y[0]*(histogram->stepX);
      
      for (i=1;i<n;i++) {
	
	sumy += y[i]*(histogram->stepX);
      }
      
      a = Y1/sumy;
      b = Y1;
      
      // test overall long-range fit:
      
      histogram2xys(histogram,Xo,999.9,0.0,999.99,x,y,sy,&n);
      
      for (i=0;i<n;i++) {
	
	sy[i] = 1.0;
	yfit[i] = (x[i] < X1) ? (nbff->Po)*exp(-k2*(sqr(x[i]-Xo))) : b*exp(-a*(x[i]-X1));
      }
      
      Z = normalised_fit_error(y,yfit,sy,n);
      
      if (Z < Zmin) {
	
	//printf("Z %10.4lf %10.4lf %10.4lf %10.4lf %10.4lf %10.4lf %10.4lf\n",Z,nbff->k1,k2,fk,MIN_K2,MIN_K_RATIO,MAX_K_RATIO);

	nbff->k2 = k2;
	
	nbff->LR_D = X1 + nbff->Dvdw;
	nbff->LR_A = a;
	nbff->LR_B = b;
	
	Zmin = Z;
      }
    }
  }

  free_histogram(hist);

  if (Zmin > 1000.0) {

    return(1);
  }

  // fit long-range Lennard-Jones:

  Y1 = 0.50*Yo;
  X1 = (Y1 > log(nbff->LR_B)) ? sqrt((Yo-Y1)/(nbff->k2)) : ((log(nbff->LR_B)-Y1)/(nbff->LR_A)) + (nbff->LR_D) - (nbff->Do);

  if (!calc_parabola_x_intercepts(-Yo,2.0*Yo,-Y1,&x1,&x2)) {

    error_fn("%s: couldn't fit Lennard-Jones to long-range section of potential (1)",__func__);
  }

  if ((x1 < 0.0) || (x2 < 0.0)) {

    error_fn("%s: couldn't fit Lennard-Jones to long-range section of potential (2)",__func__);
  }

  n1 = log(x1)/log((nbff->Do)/(X1+(nbff->Do)));
  n2 = log(x2)/log((nbff->Do)/(X1+(nbff->Do)));

  if ((n1 < 0.0) && (n2 < 0.0)) {

    error_fn("%s: couldn't fit Lennard-Jones to long-range section of potential (2)",__func__);
  }

  nbff->LR_n = (n1 > 0.0) ? n1 : n2;

  return(0);
}



static void recalc_R_histogram_errors(HISTOGRAM *histogram) {

  int i,j,n_pts,di;
  double x[100],y[100],sy[100],yfit[100],c[3],af[100],a,suma;
  HISTOGRAM_POINT *point,*spoint;

  for (i=0,point=histogram->points;i<histogram->n_points;i++,point++) {

    n_pts = 0;
      
    for (j=0,spoint=histogram->points;j<histogram->n_points;j++,spoint++) {
      
      if ((fabs(point->X - spoint->X) < 0.5) && (spoint->N > 0)) {
	
	x[n_pts] = spoint->X;
	y[n_pts] = spoint->Y/((double) spoint->N);
	
	sy[n_pts] = 1.0;
	
	n_pts++;
      }
    }
 
    if (n_pts > 3) {

      weighted_polynomial_fit(x,y,sy,n_pts,1,c,yfit);
    
      af[i] = polynomial(point->X,c,1);

    } else {

      af[i] = -1.0;
    }
  }

  for (i=0,point=histogram->points;i<histogram->n_points;i++,point++) {

    if (af[i] < 0.0) {

      di = 1;

      do {

	n_pts = 0;
	suma = 0.0;

	for (j=i-di;j<=i+di;j++) {

	  if ((j >= 0) && (j < histogram->n_points) && (af[j] > 0.0)) {

	    suma += af[j];
	    
	    n_pts++;
	  }
	}

	di++;

      } while (!n_pts);

      a = suma/((double) n_pts);

    } else {

      a = af[i];
    }

    if (point->N) {
      
      point->sY = a*sqrt((double) point->N);
      
    } else {

      point->sY = a*((double) histogram->sumN)*(1.0 - pow(0.32,1.0/((double) histogram->sumN)));
    }

    point->sYraw = point->sY;
  }
}



static double calc_2d_correction(int id,double (*g)(double)) {

  int i,conical;
  double sumF,sumGF,Y,F,G;
  HISTOGRAM *histogram;
  HISTOGRAM_POINT *point;

  if (id == -1) {

    return(1.0);
  }

  histogram = get_histogram(id);

  if (histogram == NULL) {

    error_fn("calc_2d_correction: histogram %d undefined",id);
  }

  conical = (!strncmp(histogram->name,"ALPHA",5)) ? 1 : 0;

  sumF = sumGF = 0.0;

  for (i=0,point=histogram->points;i<histogram->n_points;i++,point++) {

    Y = (histogram->log_scale) ? exp(point->Y) : point->Y;

    F = (conical) ? Y*sin((PI/180.0)*(point->X)) : Y;
    G = (*g)(point->X);

    sumF += F;
    sumGF += (G*F);
  }

  return(sumF/sumGF);
}



static void apply_2d_correction(int id,double f) {

  int i;
  HISTOGRAM *histogram;
  HISTOGRAM_POINT *point;

  if (id == -1) {

    return;
  }

  histogram = get_histogram(id);

  if (histogram == NULL) {

    error_fn("apply_2d_correction: histogram %d undefined",id);
  }

  for (i=0,point=histogram->points;i<histogram->n_points;i++,point++) {

    point->Y *= f;
  }
}



static void init_ami_ff(FORCE_FIELD *ff) {

  int flag;
  char *pli_dir,filename[1000],line[MAX_LINE_LEN],word[MAX_LINE_LEN];
  HISTOGRAM *histogram;
  PLI_FILE *file;

  pli_dir = get_pli_dir();

  sprintf(filename,"%s/ff/nonbonded/u.pliff",pli_dir);

  file = open_file(filename,"r");

  if (file == NULL) {

    error_fn("init_ami_ff: failed to open ff file '%s'",filename);
  }

  ff->AMI_HISTOGRAM_ID = -1;

  while (!end_of_file(file)) {

    if (read_line(line,MAX_LINE_LEN,file) == NULL)
      break;

    flag = sscanf(line,"%s",word);

    if (flag != EOF) {

      if (!strcmp(word,"histogram")) {

	histogram = read_histogram(file,line);

	ff->AMI_HISTOGRAM_ID = histogram->id;
      }
    }
  }

  close_file(file);

  if (ff->AMI_HISTOGRAM_ID == -1) {

    error_fn("init_ami_ff: failed to read ff file '%s'",filename);
  }
}



static void init_atom_ff_list(ATOM_FF *atom_ff_list,ATOM_TYPING_SCHEME *scheme,FORCE_FIELD *ff) {

  int i;
  ATOM_TYPE *type;

  for (i=0,type=scheme->atom_types;i<scheme->n_atom_types;i++,type++) {

    init_atom_ff(atom_ff_list + i,type,ff);
  }
}



static void init_nonbonded_ff_matrix(NONBONDED_FF **nonbonded_ff_matrix,ATOM_TYPING_SCHEME *scheme,FORCE_FIELD *ff) {

  int i,j,n_types;
  ATOM_TYPE *type1,*type2;

  n_types = scheme->n_atom_types;

  for (i=0,type1=scheme->atom_types;i<n_types;i++,type1++) {

    for (j=0,type2=scheme->atom_types;j<n_types;j++,type2++) {

      init_nonbonded_ff(nonbonded_ff_matrix[i] + j,type1,type2,ff);
    }
  }
}



static void init_atom_ff(ATOM_FF *atom_ff,ATOM_TYPE *type,FORCE_FIELD *ff) {

  int Z;

  atom_ff->ready = 0;
  atom_ff->type = type;

  atom_ff->avg_total_area = 0.0;
  atom_ff->std_total_area = 0.0;

  atom_ff->avg_inter_area = 0.0;
  atom_ff->std_inter_area = 0.0;

  atom_ff->max_polar_area = 0.0;

  atom_ff->Zp_histogram_id = -1;
  atom_ff->Zf_histogram_id = -1;

  atom_ff->force_field = ff;

  for (Z=0;Z<MAX_Z;Z++) {

    atom_ff->A_histogram_ids[Z] = -1;
  }
}



static void init_nonbonded_ff(NONBONDED_FF *nonbonded_ff,ATOM_TYPE *type1,ATOM_TYPE *type2,FORCE_FIELD *ff) {

  nonbonded_ff->surrogate = 0;
  nonbonded_ff->ready = 0;
  nonbonded_ff->usage_count = 0;

  nonbonded_ff->type1 = type1;
  nonbonded_ff->type2 = type2;

  nonbonded_ff->P = -1.0;
  nonbonded_ff->sP = 0.0;
  nonbonded_ff->lnP = 0.0;
  nonbonded_ff->slnP = 0.0;

  strcpy(nonbonded_ff->precision,"undefined");
  strcpy(nonbonded_ff->shape,"undefined");

  nonbonded_ff->Dvdw = get_atom_type_vdw_radius(type1) + get_atom_type_vdw_radius(type2);

  nonbonded_ff->Do = 0.0;
  nonbonded_ff->Eo = 0.0;
  nonbonded_ff->Po = 1.0;
  nonbonded_ff->Dc = 0.0;

  nonbonded_ff->k1 = 0.0;
  nonbonded_ff->k2 = 0.0;

  nonbonded_ff->opt_Z = 0.0;
  nonbonded_ff->opt_sY = 0.0;

  nonbonded_ff->LJ_n = 6.0;

  nonbonded_ff->LJ_Z = 0.0;
  nonbonded_ff->LJ_sY = 0.0;

  nonbonded_ff->LR_A = 0.0;
  nonbonded_ff->LR_B = 0.0;
  nonbonded_ff->LR_D = 0.0;
  nonbonded_ff->LR_n = 6.0;

  nonbonded_ff->R_histogram_id = -1;
  nonbonded_ff->ALPHA1_histogram_id = -1;
  nonbonded_ff->BETA1_histogram_id = -1;
  nonbonded_ff->ALPHA2_histogram_id = -1;
  nonbonded_ff->BETA2_histogram_id = -1;
  nonbonded_ff->force_field = ff;
}



static FORCE_FIELD* alloc_ff(ATOM_TYPING_SCHEME *scheme) {

  FORCE_FIELD *ff;
  NONBONDED_FF **nonbonded;

  ff = (FORCE_FIELD*) malloc(sizeof(FORCE_FIELD));

  if (ff == NULL) {

    error_fn("alloc_ff: out of memory allocating force field");
  }

  ff->nonbonded = alloc_nonbonded_ff_matrix(scheme);

  ff->atom = alloc_atom_ff_list(scheme);

  return(ff);
}



static NONBONDED_FF** alloc_nonbonded_ff_matrix(ATOM_TYPING_SCHEME *scheme) {

  int i,n_types;
  NONBONDED_FF **nonbonded_ff_matrix;

  n_types = scheme->n_atom_types;

  nonbonded_ff_matrix = (NONBONDED_FF**) calloc(n_types,sizeof(NONBONDED_FF*));

  if (nonbonded_ff_matrix == NULL) {

    error_fn("alloc_nonbonded_ff_matrix: out of memory allocating nonbonded (1)");
  }

  nonbonded_ff_matrix[0] = (NONBONDED_FF*) calloc(n_types*n_types,sizeof(NONBONDED_FF));

  if (nonbonded_ff_matrix[0] == NULL) {

    error_fn("alloc_nonbonded_ff_matrix: out of memory allocating nonbonded (2)");
  }

  for (i=1;i<n_types;i++) {

    nonbonded_ff_matrix[i] = nonbonded_ff_matrix[i-1] + n_types;
  }

  return(nonbonded_ff_matrix);
}


static void reset_nonbonded_ff(NONBONDED_FF *nonbonded_ff) {

  ATOM_TYPE *type1,*type2;
  NONBONDED_FF *paired_ff,**nonbonded_ff_matrix;
  FORCE_FIELD *ff;

  type1 = nonbonded_ff->type1;
  type2 = nonbonded_ff->type2;

  ff = nonbonded_ff->force_field;

  nonbonded_ff_matrix = ff->nonbonded;

  reset_histogram(nonbonded_ff->R_histogram_id);
  reset_histogram(nonbonded_ff->ALPHA1_histogram_id);
  reset_histogram(nonbonded_ff->BETA1_histogram_id);
  reset_histogram(nonbonded_ff->ALPHA2_histogram_id);
  reset_histogram(nonbonded_ff->BETA2_histogram_id);

  init_nonbonded_ff(nonbonded_ff,type1,type2,ff);
  
  paired_ff = (type1 == type2) ? (NULL) : nonbonded_ff_matrix[type2->i] + type1->i;
  
  if (paired_ff) {

    init_nonbonded_ff(paired_ff,type2,type1,ff);
  }

  n_stored_potentials--;
}



static NONBONDED_FF* get_least_used_nonbonded_ff(FORCE_FIELD *ff) {

  int i,j,n_types;
  ATOM_TYPE *type1,*type2;
  ATOM_TYPING_SCHEME *scheme;
  NONBONDED_FF **nonbonded_ff_matrix,*nb_ff,*least_used_nb_ff;

  scheme = ff->settings->atom_typing_scheme;

  n_types = scheme->n_atom_types;

  least_used_nb_ff = NULL;

  nonbonded_ff_matrix = ff->nonbonded;

  for (i=0,type1=scheme->atom_types;i<n_types;i++,type1++) {

    for (j=i,type2=scheme->atom_types;j<n_types;j++,type2++) {

      nb_ff = nonbonded_ff_matrix[i] + j;

      if ((nb_ff->ready) && ((least_used_nb_ff == NULL) || (nb_ff->usage_count < least_used_nb_ff->usage_count))) {
	  
	least_used_nb_ff = nb_ff;
      }
    }
  }

  return(least_used_nb_ff);
}



static ATOM_FF* alloc_atom_ff_list(ATOM_TYPING_SCHEME *scheme) {

  int i,n_types;
  ATOM_FF *atom_ff;

  n_types = scheme->n_atom_types;

  atom_ff = (ATOM_FF*) calloc(n_types,sizeof(ATOM_FF));

  if (atom_ff == NULL) {

    error_fn("alloc_atom_ff_list: out of memory allocating nonbonded");
  }

  return(atom_ff);
}



static void free_ff(FORCE_FIELD *ff) {

  if (ff == NULL) {

    return;
  }

  free_nonbonded_ff_matrix(ff->nonbonded);

  if (ff->atom) {

    free(ff->atom);
  }

  free(ff);
}



static void free_nonbonded_ff_matrix(NONBONDED_FF **nonbonded_ff_matrix) {

  if (nonbonded_ff_matrix == NULL) {

    return;
  }

  free(*nonbonded_ff_matrix);
  free(nonbonded_ff_matrix);
}
