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



static struct FFSettings {
  int use_alpha_potentials;
  int use_beta_potentials;
  int use_clash_potentials;
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
static NONBONDED_FF *methyl_methyl_nonbonded_ff = NULL;
static NONBONDED_FF *last_read_ff = NULL;
static n_stored_potentials = 0;



static void assign_methyl_methyl_R_potential(NONBONDED_FF*);
static void reset_nonbonded_ff(NONBONDED_FF*);
static NONBONDED_FF* get_least_used_nonbonded_ff(FORCE_FIELD*);
static HISTOGRAM* process_atom_ff_histogram(HISTOGRAM*,ATOM_FF*);
static void process_atom_ff_histogram_old(HISTOGRAM*,ATOM_FF*);
static void process_nonbonded_ff_histogram(HISTOGRAM*,NONBONDED_FF*);
static void process_R_potential(HISTOGRAM*,NONBONDED_FF*);
static double derive_vdw_distance(NONBONDED_FF*);
static double calc_2d_correction(int,double (*g)(double));
static void apply_2d_correction(int,double);
static void read_atom_ff(ATOM_FF*);
static void read_nonbonded_ff(NONBONDED_FF*);
static void copy_nonbonded_ff(NONBONDED_FF*,NONBONDED_FF*);
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
  ff_settings->use_clash_potentials = (params_get_parameter("use_clash_potentials"))->value.i;

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

  ATOM_TYPE *methyl;
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

  // pre-load methyl-methyl nonbonded_ff for cases where the R distribution isnt known:

  methyl = get_atom_type("-CH3",settings->atom_typing_scheme);

  methyl_methyl_nonbonded_ff = get_nonbonded_ff(ff->nonbonded,methyl,methyl);

  if ((methyl_methyl_nonbonded_ff == NULL) || (methyl_methyl_nonbonded_ff->R_histogram_id == -1)) {

    error_fn("init_ff: methyl-methyl R distribution undefined");
  }

  histogram = get_histogram(methyl_methyl_nonbonded_ff->R_histogram_id);

  histogram->protected = 1;
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

  NONBONDED_FF *nonbonded_ff,*paired_ff;

  if ((type1 == NULL) || (type2 == NULL)) {

    error_fn("get_nonbonded_ff: type(s) undefined");
  }

  nonbonded_ff = nonbonded_ff_matrix[type1->i] + type2->i;
  paired_ff = (type1 == type2) ? (NULL) : nonbonded_ff_matrix[type2->i] + type1->i;

  if (nonbonded_ff->ready == 0) {

    read_nonbonded_ff(nonbonded_ff);

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



static void read_nonbonded_ff(NONBONDED_FF *nonbonded_ff) {

  int i,i1,i2,flag,id;
  char *pli_dir,name1[10],name2[10],filename[1000],line[MAX_LINE_LEN],word[MAX_LINE_LEN];
  double f_r_alpha,f_r_beta,f_alpha1_beta,f_alpha2_beta;
  ATOM_TYPE *type1,*type2;
  HISTOGRAM *histogram,*methyl_methyl_R_hist;
  HISTOGRAM_POINT *point;
  NONBONDED_FF *nb_ff;
  PLI_FILE *file;

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

    error_fn("read_nonbonded_ff: failed to open ff file '%s'",filename);
  }

  while (!end_of_file(file)) {

    if (read_line(line,MAX_LINE_LEN,file) == NULL)
      break;

    flag = sscanf(line,"%s",word);

    if (flag != EOF) {

      if (!strcmp(word,"P")) {

	sscanf(line,"%*s %lf %lf",&(nonbonded_ff->P),&(nonbonded_ff->sP));

      } else if (!strcmp(word,"lnP")) {

	sscanf(line,"%*s %lf %lf",&(nonbonded_ff->lnP),&(nonbonded_ff->slnP));

      } else if (!strcmp(word,"Precision")) {

	sscanf(line,"%*s %s",nonbonded_ff->precision);

      } else if (!strcmp(word,"histogram")) {

	histogram = read_histogram(file,line);

	if (!strcmp(histogram->name,"R")) {

	  histogram->max_peaks = 1;

	  nonbonded_ff->R_histogram_id = histogram->id;

	} else if ((!strcmp(histogram->name,"ALPHA1")) && (ff_settings->use_alpha_potentials)) {

	  histogram->max_peaks = 2;

	  if (i1 < i2) {

	    nonbonded_ff->ALPHA1_histogram_id = histogram->id;

	  } else {

	    nonbonded_ff->ALPHA2_histogram_id = histogram->id;
	  }

	} else if ((!strcmp(histogram->name,"BETA1")) && (ff_settings->use_beta_potentials)) {

	  histogram->max_peaks = 2;

	  if (i1 < i2) {

	    nonbonded_ff->BETA1_histogram_id = histogram->id;

	  } else {

	    nonbonded_ff->BETA2_histogram_id = histogram->id;
	  }

	} else if ((!strcmp(histogram->name,"ALPHA2")) && (ff_settings->use_alpha_potentials)) {

	  histogram->max_peaks = 2;

	  if (i1 < i2) {

	    nonbonded_ff->ALPHA2_histogram_id = histogram->id;

	  } else {

	    nonbonded_ff->ALPHA1_histogram_id = histogram->id;
	  }

	} else if ((!strcmp(histogram->name,"BETA2")) && (ff_settings->use_beta_potentials)) {

	  histogram->max_peaks = 2;

	  if (i1 < i2) {

	    nonbonded_ff->BETA2_histogram_id = histogram->id;

	  } else {

	    nonbonded_ff->BETA1_histogram_id = histogram->id;
	  }
	}

	process_nonbonded_ff_histogram(histogram,nonbonded_ff);
      }
    }
  }
  
  close_file(file);

  f_r_alpha = calc_2d_correction(nonbonded_ff->R_histogram_id,ff_settings->g_r_alpha);
  f_r_beta = calc_2d_correction(nonbonded_ff->R_histogram_id,ff_settings->g_r_beta);

  f_alpha1_beta = calc_2d_correction(nonbonded_ff->ALPHA1_histogram_id,ff_settings->g_alpha_beta);
  f_alpha2_beta = calc_2d_correction(nonbonded_ff->ALPHA2_histogram_id,ff_settings->g_alpha_beta);

  apply_2d_correction(nonbonded_ff->ALPHA1_histogram_id,f_r_alpha);
  apply_2d_correction(nonbonded_ff->ALPHA2_histogram_id,f_r_alpha);

  apply_2d_correction(nonbonded_ff->BETA1_histogram_id,f_r_beta*f_alpha1_beta);
  apply_2d_correction(nonbonded_ff->BETA2_histogram_id,f_r_beta*f_alpha2_beta);

  if (nonbonded_ff->R_histogram_id == -1) {

    // no distance histogram found - use methyl-methyl one instead:

    assign_methyl_methyl_R_potential(nonbonded_ff);
  }

  nonbonded_ff->ready = 1;

  last_read_ff = nonbonded_ff;

  n_stored_potentials++;
}



void copy_nonbonded_ff(NONBONDED_FF *nonbonded_ff1,NONBONDED_FF *nonbonded_ff2) {

  if (nonbonded_ff1 == nonbonded_ff2) {

    return;
  }

  if ((nonbonded_ff1->type1 != nonbonded_ff2->type2) || (nonbonded_ff1->type2 != nonbonded_ff2->type1)) {

    error_fn("copy_nonbonded_ff: atom types are not matched");
  }

  nonbonded_ff2->ready = 1;
  nonbonded_ff2->P = nonbonded_ff1->P;
  nonbonded_ff2->sP = nonbonded_ff1->sP;
  nonbonded_ff2->lnP = nonbonded_ff1->lnP;
  nonbonded_ff2->slnP = nonbonded_ff1->slnP;

  nonbonded_ff2->Do = nonbonded_ff1->Do;
  nonbonded_ff2->Dc = nonbonded_ff1->Dc;

  nonbonded_ff2->use_clash_potential = nonbonded_ff1->use_clash_potential;

  nonbonded_ff2->LJ_A = nonbonded_ff1->LJ_A;
  nonbonded_ff2->LJ_B = nonbonded_ff1->LJ_B;

  nonbonded_ff2->R_histogram_id = nonbonded_ff1->R_histogram_id;
  nonbonded_ff2->ALPHA1_histogram_id = nonbonded_ff1->ALPHA2_histogram_id;
  nonbonded_ff2->BETA1_histogram_id = nonbonded_ff1->BETA2_histogram_id;
  nonbonded_ff2->ALPHA2_histogram_id = nonbonded_ff1->ALPHA1_histogram_id;
  nonbonded_ff2->BETA2_histogram_id = nonbonded_ff1->BETA1_histogram_id;
  nonbonded_ff2->force_field = nonbonded_ff2->force_field;
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



static void assign_methyl_methyl_R_potential(NONBONDED_FF *nonbonded_ff) {

  HISTOGRAM *histogram;

  nonbonded_ff->R_histogram_id = methyl_methyl_nonbonded_ff->R_histogram_id;

  histogram = get_histogram(nonbonded_ff->R_histogram_id);

  process_R_potential(histogram,nonbonded_ff);
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
    histogram->sharpen = 1;

  } else if (nonbonded_ff->ALPHA1_histogram_id == histogram->id) {

    histogram->limits = geom1->alpha_limits;

    histogram->geometric_correction = CONICAL_CORRECTION;
    histogram->sharpen = 1;

  } else if (nonbonded_ff->BETA1_histogram_id == histogram->id) {

    histogram->limits = geom1->beta_limits;
    histogram->sharpen = 1;

  } else if (nonbonded_ff->ALPHA2_histogram_id == histogram->id) {

    histogram->limits = geom2->alpha_limits;
    histogram->sharpen = 1;

    histogram->geometric_correction = CONICAL_CORRECTION;

  } else if (nonbonded_ff->BETA2_histogram_id == histogram->id) {

    histogram->limits = geom2->beta_limits;
    histogram->sharpen = 1;
  }

  histogram->extrapolation_method = EXTRAPOLATE_TO_SMALLEST;
  histogram->normalisation_method = RELATIVE_NORMALISATION;
  histogram->log_scale = 1;

  process_histogram(histogram);

  if (nonbonded_ff->R_histogram_id == histogram->id) {

    process_R_potential(histogram,nonbonded_ff);
  }

  set_histogram_extrema(histogram);
}



static void process_R_potential(HISTOGRAM *histogram,NONBONDED_FF *nonbonded_ff) {

  int i,j,flag;
  double D,slope,intercept,Dc,Xo,Yo,a,b;
  ATOM_TYPE *type1,*type2;
  HISTOGRAM_POINT *point,*point1,*point2,*point3;

  type1 = nonbonded_ff->type1;
  type2 = nonbonded_ff->type2;

  D = get_atom_type_vdw_radius(type1) + get_atom_type_vdw_radius(type2);

  // find two points to calculate clash distance from:

  i = 0;
  point1 = histogram->points;

  if (point1->Y > 0.0) {

    warning_fn("process_R_potential: potential starts positive for %s - %s",type1->name,type2->name);

    point2 = point1 + 1;

  } else {

    point2 = histogram->points;

    do {

      i++;
      point2++;

    } while ((i < histogram->n_points) && (point2->Y < 0.0));

    if (i == histogram->n_points) {

      error_fn("process_R_potential: corrupt potential (1)");
    }
  
    point1 = point2 - 1;
  }

  // intrapolate (or extrapolate if potential starts positive) to get clash distance:

  flag = solve_line(point1->X,point1->Y,point2->X,point2->Y,&slope,&intercept);

  if ((flag) || (slope < 0.0)) {

    error_fn("process_R_potential: corrupt potential (2)");
  }

  Dc = (fabs(slope) < 1.0E-30) ? D : D - (intercept/slope);

  // fit LJ potential to clash distance and slope:

  b = -(slope*pow(Dc,7))/6;

  a = b*pow(Dc,6);

  if (nonbonded_ff->use_clash_potential) {

    for (j=0,point=histogram->points;j<i;j++,point++) {

      point->Y = lennard_jones_energy(point->X + D,a,b);
    }
  }

  // now find three points around the optimum:

  point3 = point2;
  point2 = point1;

  do {

    i++;
    point2++;
    point3++;

  } while ((i < histogram->n_points) && (point3->Y > point2->Y));

  if (i == histogram->n_points) {

    error_fn("process_R_potential: corrupt potential (3)");
  }

  point1 = point2 - 1;
 
  // fit parabola to three points to find optimum distance:

  flag = calc_parabola_vertex(point1->X,point1->Y,point2->X,point2->Y,point3->X,point3->Y,&Xo,&Yo);

  if (flag) {

    error_fn("read_nonbonded_ff: unable to derive optimum distance from R histogram");
  }

  nonbonded_ff->Dc = Dc;
  nonbonded_ff->Do = D + Xo;

  nonbonded_ff->LJ_A = a;
  nonbonded_ff->LJ_B = b;
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

  nonbonded_ff->ready = 0;
  nonbonded_ff->usage_count = 0;

  nonbonded_ff->type1 = type1;
  nonbonded_ff->type2 = type2;

  nonbonded_ff->P = 0.0;
  nonbonded_ff->sP = 0.0;
  nonbonded_ff->lnP = 1.0;
  nonbonded_ff->slnP = 0.0;

  strcpy(nonbonded_ff->precision,"undefined");

  nonbonded_ff->Do = 0.0;
  nonbonded_ff->Dc = 0.0;

  nonbonded_ff->use_clash_potential = ff_settings->use_clash_potentials;

  nonbonded_ff->LJ_A = 0.0;
  nonbonded_ff->LJ_B = 0.0;

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
