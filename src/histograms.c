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



static struct HistogramSettings {
  double start_s;
  double end_s;
  double step_s;
  double max_point_zscore;
  double max_point_error;
} *histogram_settings = NULL;



static HISTOGRAM *histograms = NULL;
static int n_histograms = 0;
static int n_alloc_histograms = 0;
static int histogram_batch_size = 2;
static int histogram_point_batch_size = 20;



static double histogram_intrapolate(HISTOGRAM*,double);
static int count_histogram_maxima(HISTOGRAM*);
static int local_histogram_maximum(HISTOGRAM_POINT*,int,int);
static int local_histogram_minimum(HISTOGRAM_POINT*,int,int);
static void init_histogram(HISTOGRAM*);
static void normalise_histogram(HISTOGRAM*,int);
static void smoothe_histogram(HISTOGRAM*,HISTOGRAM*,double);
static void sharpen_histogram(HISTOGRAM*,HISTOGRAM*,double);
static int remove_local_extreme(HISTOGRAM*);
static double point_intrapolate(HISTOGRAM_POINT*,int,int,double*);
static double two_point_extrapolate(HISTOGRAM_POINT*,HISTOGRAM_POINT*,double);
static void subtract_histograms(HISTOGRAM*,HISTOGRAM*,HISTOGRAM*);
static void log_histogram(HISTOGRAM*);
static double avg_histogram_point_error(HISTOGRAM*,HISTOGRAM*);
static double avg_histogram_point_zscore(HISTOGRAM*,HISTOGRAM*);
static void free_histogram(HISTOGRAM*);
static void alloc_histograms(void);
static void alloc_histogram_points(HISTOGRAM*);



void init_histogram_settings(void) {

  if (histogram_settings == NULL) {

    histogram_settings = (struct HistogramSettings*) malloc(sizeof(struct HistogramSettings));
  }

  if (histogram_settings == NULL) {

    error_fn("init_histogram_settings: out of memory allocating settings");
  }

  histogram_settings->start_s = (params_get_parameter("hist_start_s"))->value.d;
  histogram_settings->end_s = (params_get_parameter("hist_end_s"))->value.d;
  histogram_settings->step_s = (params_get_parameter("hist_step_s"))->value.d;

  histogram_settings->max_point_zscore = (params_get_parameter("hist_max_point_zscore"))->value.d;
  histogram_settings->max_point_error = (params_get_parameter("hist_max_point_error"))->value.d;
}



HISTOGRAM* get_histogram(int histogram_id) {

  if (histogram_id == -1) {

    return(NULL);
  }

  if (histogram_id >= n_histograms) {

    error_fn("get_histogram: histogram %d does not exist",histogram_id);
  }

  return(histograms + histogram_id);
}



void run_histogram(SETTINGS *settings) {

  int hist_id;
  ATOM_TYPE *type1,*type2;
  ATOM_TYPING_SCHEME *scheme;
  NONBONDED_FF *nonbonded_ff;
  HISTOGRAM *raw_histogram,*histogram;

  scheme = settings->atom_typing_scheme;

  if ((settings->id1 == -1) || (settings->id2 == -1)) {

    error_fn("histogram: id1 or id2 undefined");
  }

  if (!strcmp(settings->hist_name,"UNKNOWN")) {

    error_fn("histogram: hist_name undefined");
  }

  type1 = get_atom_type_by_id(settings->id1,scheme);
  type2 = get_atom_type_by_id(settings->id2,scheme);

  nonbonded_ff = get_nonbonded_ff(settings->force_field->nonbonded,type1,type2);

  hist_id = -1;

  if (!strcmp(settings->hist_name,"R")) {

    hist_id = nonbonded_ff->R_histogram_id;

  } else if (!strcmp(settings->hist_name,"ALPHA1")) {

    hist_id = nonbonded_ff->ALPHA1_histogram_id;

  } else if (!strcmp(settings->hist_name,"BETA1")) {

    hist_id = nonbonded_ff->BETA1_histogram_id;
  }


  if (hist_id != -1) {

    histogram = histograms + hist_id;
  }

  if (!strcmp(settings->mode->name,"histogram")) {

    write_histogram(histogram,type1,type2);
  }
}



HISTOGRAM *read_histogram(PLI_FILE *file,char *headline) {

  int flag;
  char line[MAX_LINE_LEN],word[MAX_LINE_LEN];
  HISTOGRAM *histogram;
  HISTOGRAM_POINT *point;

  histogram = new_histogram();

  sscanf(headline,"%s %s %*s %d %lf %lf %lf %lf",
	 word,histogram->name,
	 &(histogram->sumN),&(histogram->sumY),
	 &(histogram->startX),&(histogram->endX),&(histogram->stepX));

  while ((!end_of_file(file)) && (strcmp(word,"end_histogram"))) {
  
    if (read_line(line,MAX_LINE_LEN,file) == NULL)
      break;

    flag = sscanf(line,"%s",word);

    if (flag != EOF) {

      if (strcmp(word,"end_histogram")) {

	alloc_histogram_points(histogram);

	point = histogram->points + histogram->n_points;

	sscanf(line,"%lf %lf %lf",&(point->X),&(point->Y),&(point->sY));

	histogram->n_points++;
      }
    }
  }

  return(histogram);
}



HISTOGRAM *new_histogram(void) {

  int i;
  HISTOGRAM *histogram;

  for (i=0,histogram=histograms;i<n_histograms;i++,histogram++) {

    if (!histogram->used) {

      histogram->used = 1;

      return(histogram);
    }
  }

  alloc_histograms();

  histogram = histograms + n_histograms;

  init_histogram(histogram);

  histogram->used = 1;

  n_histograms++;

  return(histogram);
}



HISTOGRAM* add_histogram(char *name,int n_points) {

  int i;
  HISTOGRAM *histogram;

  alloc_histograms();

  histogram = histograms + n_histograms;

  init_histogram(histogram);

  strcpy(histogram->name,name);

  for (i=0;i<n_points;i++) {

    alloc_histogram_points(histogram);

    histogram->n_points++;
  }

  n_histograms++;

  return(histogram);
}



void process_histogram(HISTOGRAM *histogram) {

  int i;
  double sX,start_sX,end_sX,step_sX,avg_point_error,avg_point_zscore,R1,R2,X;
  ATOM_TYPE *type1,*type2;
  HISTOGRAM *new_histogram;
  HISTOGRAM_POINT *point1,*point2;

  // normalise histogram (including conical correction where needed):

  normalise_histogram(histogram,histogram->geometric_correction);

  // initially try smoothing using bin width as sigma:

  new_histogram = clone_histogram(histogram);

  start_sX = (histogram->start_s)*(histogram->stepX);
  end_sX = 1.0001*(histogram->end_s)*(histogram->stepX);
  step_sX = (histogram->step_s)*(histogram->stepX);

  sX = start_sX;

  smoothe_histogram(histogram,new_histogram,sX);

  avg_point_zscore = avg_histogram_point_zscore(new_histogram,histogram);

  // if that gives too much deviation from original, take half the bin width as sigma:

  if (avg_point_zscore >= histogram_settings->max_point_zscore) {

    sX *= 0.5;

    smoothe_histogram(histogram,new_histogram,sX);

    avg_point_zscore = avg_histogram_point_zscore(new_histogram,histogram);
  }

  avg_point_error = avg_histogram_point_error(new_histogram,histogram);

  // if required try increasingly higher values for sigma until error is below zero
  // or deviation from original histogram becomes too large

  if (((avg_point_error > histogram_settings->max_point_error) && (avg_point_zscore < histogram_settings->max_point_zscore)) ||
      (count_histogram_maxima(new_histogram) > histogram->max_peaks)) {

    while ((sX+step_sX < end_sX) && 
	   (((avg_point_error > histogram_settings->max_point_error) &&
	     (avg_point_zscore < histogram_settings->max_point_zscore)) ||
	    (count_histogram_maxima(new_histogram) > histogram->max_peaks))) {

      sX += step_sX;

      smoothe_histogram(histogram,new_histogram,sX);

      avg_point_error = avg_histogram_point_error(new_histogram,histogram);

      avg_point_zscore = avg_histogram_point_zscore(new_histogram,histogram);
    }

    // back track if deviation from original was too large:

    if (avg_point_zscore >= histogram_settings->max_point_zscore) {

      sX -= step_sX;

      smoothe_histogram(histogram,new_histogram,sX);    

      avg_point_error = avg_histogram_point_error(new_histogram,histogram);

      avg_point_zscore = avg_histogram_point_zscore(new_histogram,histogram);
    }
  }

  memcpy(histogram->points,new_histogram->points,(histogram->n_points)*sizeof(HISTOGRAM_POINT));

  while (remove_local_extreme(histogram)) { }
  
  free_histogram(new_histogram);

  // normalise histogram again - smoothing etc will have changed it - dont do conical correction again though:

  normalise_histogram(histogram,0);

  if (histogram->log_scale) {

    log_histogram(histogram);
  }
}



void set_histogram_extrema(HISTOGRAM *histogram) {

  int i;
  HISTOGRAM_POINT *point,*min_point,*max_point;

  min_point = NULL;
  max_point = NULL;

  for (i=0,point=histogram->points;i<histogram->n_points;i++,point++) {

    if ((min_point == NULL) || (point->Y < min_point->Y)) {

      min_point = point;
    }

    if ((max_point == NULL) || (point->Y > max_point->Y)) {

      max_point = point;
    }
  }

  if ((min_point == NULL) || (max_point == NULL)) {

    error_fn("set_histogram_extrema: corrupt histogram '%s' (%d)",histogram->name,histogram->id);
  }

  histogram->minX = min_point->X;
  histogram->minY = min_point->Y;

  histogram->maxX = max_point->X;
  histogram->maxY = max_point->Y;

  if (histogram->extrapolation_method == EXTRAPOLATE_TO_SMALLEST) {

    histogram->startY = histogram->minY;
    histogram->endY = histogram->minY;

  } else if (histogram->extrapolation_method == EXTRAPOLATE_TO_NEAREST) {

    histogram->startY = histogram_intrapolate(histogram,histogram->startX);
    histogram->endY = histogram_intrapolate(histogram,histogram->endX);
  }
}



double histogram_X2Y(HISTOGRAM *histogram,double X) {

  double Y;

  if (X < histogram->startX) {

    return(histogram->startY);
  }

  if (X > histogram->endX) {

    return(histogram->endY);
  }

  Y = histogram_intrapolate(histogram,X);

  return(Y);
}



static double histogram_intrapolate(HISTOGRAM *histogram,double X) {

  int bin;
  double X1,X2,Y1,Y2,Y,f1,f2,f;
  HISTOGRAM_POINT *points;

  points = histogram->points;

  f1 = X-(histogram->startX)-0.5*(histogram->stepX);
  f2 = histogram->stepX;

  f = (X-(histogram->startX)-0.5*(histogram->stepX))/(histogram->stepX);
  bin = floor(f);

  X1 = (histogram->startX) + (1.0*bin + 0.5)*(histogram->stepX);
  X2 = X1 + (histogram->stepX);

  if (bin < 0) {

    Y1 = Y2 = points->Y;

  } else if (bin > (histogram->n_points)-2) {

    Y1 = Y2 = (points + (histogram->n_points) - 1)->Y;

  } else {

    Y1 = (points + bin)->Y;
    Y2 = (points + bin + 1)->Y;
  }


  Y = Y1 + (Y2-Y1)*((X-X1)/(X2-X1));

  return(Y);
}



static int count_histogram_maxima(HISTOGRAM *histogram) {

  int i,n_points,n_max;
  HISTOGRAM_POINT *point,*point1,*point2;

  if (histogram->n_points < 3) {

    return(0);
  }

  n_points = histogram->n_points;

  n_max = 0;

  for (i=0,point=histogram->points;i<n_points;i++,point++) {

    if (local_histogram_maximum(point,i,n_points)) {

      n_max++;
    }
  }

  return(n_max);
}



static int local_histogram_maximum(HISTOGRAM_POINT *point,int i,int n_points) {

  double x;
  HISTOGRAM_POINT *point1,*point2;

  if ((i < 0) || (i >= n_points)) {

    error_fn("local_histogram_maximum: point out of range");
  }

  if (i == 0) {

    point1 = point + 1;

    if (point->Y > point1->Y) {

      return(1);
    }

  } else if (i == n_points - 1) {

    point1 = point - 1;

    if (point->Y > point1->Y) {

      return(1);
    }

  } else {

    if (point->Y > 1) { x=1;}

    point1 = point - 1;
    point2 = point + 1;
    
    if (point1->Y > 1) { x=1;}
    if (point2->Y > 1) { x=1;}

    if ((point->Y > point1->Y) && (point->Y > point2->Y)) {
      
      return(1);
    }
  }

  return(0);
}



static int local_histogram_minimum(HISTOGRAM_POINT *point,int i,int n_points) {

  HISTOGRAM_POINT *point1,*point2;

  if (i == 0) {

    point1 = point + 1;

    if (point->Y < point1->Y) {

      return(1);
    }

  } else if (i == n_points - 1) {

    point1 = point - 1;

    if (point->Y < point1->Y) {

      return(1);
    }

  } else {
    
    point1 = point - 1;
    point2 = point + 1;
    
    if ((point->Y < point1->Y) && (point->Y < point2->Y)) {
      
      return(1);
    }
  }

  return(0);
}



void write_histogram(HISTOGRAM *histogram,ATOM_TYPE *type1,ATOM_TYPE *type2) {

  int i;
  HISTOGRAM_POINT *point;

  for (i=0,point=histogram->points;i<histogram->n_points;i++,point++) {

    if (!strcmp(histogram->name,"R")) {

      if ((type1 == NULL) || (type2 == NULL)) {

	printf("%10.4lf %10.4lf %10.4lf\n",point->X,point->Y,point->sY);

      } else {

	printf("%10.4lf %10.4lf %10.4lf\n",point->X + type1->united_atom->vdw_radius + type2->united_atom->vdw_radius,point->Y,point->sY);
      }

    } else {

      printf("%10.4lf %10.4lf %10.4lf\n",point->X,point->Y,point->sY);
    }
  }
}



void reset_histogram(int id) {

  HISTOGRAM *histogram;

  histogram = get_histogram(id);

  if ((!histogram) || (histogram->protected)) {

    return;
  }

  if (histogram->n_points) {

    free(histogram->points);
  }

  init_histogram(histogram);

  histogram->id = id;
}



static void init_histogram(HISTOGRAM *histogram) {

  histogram->id = n_histograms;

  histogram->used = 0;
  histogram->protected = 0;

  strcpy(histogram->name,"");

  histogram->startX = 0.0;
  histogram->endX = 0.0;
  histogram->stepX = 0.0;

  histogram->startY = 0.0;
  histogram->endY = 0.0;

  histogram->minX = 0.0;
  histogram->minY = 0.0;
  histogram->maxX = 0.0;
  histogram->maxY = 0.0;

  histogram->n_alloc_points = 0;
  histogram->n_points = 0;
  histogram->points = NULL;

  histogram->start_s = histogram_settings->start_s;
  histogram->end_s = histogram_settings->end_s;
  histogram->step_s = histogram_settings->step_s;

  histogram->limits = NULL;

  histogram->max_peaks = 9999;
  histogram->extrapolation_method = EXTRAPOLATE_TO_SMALLEST;
  histogram->geometric_correction = NO_GEOMETRIC_CORRECTION;
  histogram->normalisation_method = ABSOLUTE_NORMALISATION;

  histogram->sharpen = 0;
  histogram->sharpen_sX = 0.0;

  histogram->log_scale = 0;
}



HISTOGRAM* clone_histogram(HISTOGRAM *histogram) {

  HISTOGRAM *new_histogram;

  new_histogram = (HISTOGRAM*) malloc(sizeof(HISTOGRAM));

  if (new_histogram == NULL) {

    error_fn("clone_histogram: out of memory allocating histogram");
  }

  memcpy(new_histogram,histogram,sizeof(HISTOGRAM));

  if (new_histogram->n_points) {

    new_histogram->points = (HISTOGRAM_POINT*) calloc(new_histogram->n_points,sizeof(HISTOGRAM_POINT));

    if (new_histogram->points == NULL) {

      error_fn("clone_histogram: out of memory allocating histogram points");
    }

    memcpy(new_histogram->points,histogram->points,(new_histogram->n_points)*sizeof(HISTOGRAM_POINT));
  }

  return(new_histogram);
}



static void normalise_histogram(HISTOGRAM *histogram,int geometric_correction) {

  int i,conical;
  double stepX,intY,intF,C,F;
  HISTOGRAM_POINT *point;

  intY = 0.0;
  intF = 0.0;

  conical = (geometric_correction == CONICAL_CORRECTION) ? 1 : 0;
  stepX = histogram->stepX;

  for (i=0,point=histogram->points;i<histogram->n_points;i++,point++) {

    if ((histogram->limits == NULL) || ((point->X > histogram->limits[0]) && (point->X < histogram->limits[1]))) {

      intY += (point->Y)*stepX;

      intF += (conical) ? (sin((PI/180.0)*(point->X))*stepX) : stepX;
    }
  }

  if (intF < 1.0E-10) {

    error_fn("normalise_histogram: intF (nearly) zero (%8e) for histogram %s",intF,histogram->name);
  }

  C = (histogram->normalisation_method == RELATIVE_NORMALISATION) ? (intF/intY) : (1.0/intY);

  for (i=0,point=histogram->points;i<histogram->n_points;i++,point++) {

    F = (conical) ? sin((PI/180.0)*(point->X)) : 1.0;

    point->Y *= (C/F);
    point->sY *= (C/F);
  }
}



static void smoothe_histogram(HISTOGRAM *histogram,HISTOGRAM *smooth_histogram,double sX) {

  int n_bins,i,j,k,l,periodic;
  double dX,Y,sY,f,sumY,sumSY2,sumf;
  HISTOGRAM_POINT *point1,*point2;

  periodic = (!strcmp(histogram->name,"R")) ? 0 : 1;

  n_bins = ceil((4.0*sX)/histogram->stepX);

  for (i=0,point1=histogram->points,point2=smooth_histogram->points;i<histogram->n_points;i++,point1++,point2++) {

    sumY = sumf = sumSY2 = 0.0;

    for (j=-n_bins;j<=n_bins;j++) {

      dX = (histogram->stepX)*j;

      k = i + j;

      l = -1;

      if (k < 0) {

	if ((-k < histogram->n_points) && (periodic)) {

	  l = -k;
	}

      } else if (k >= histogram->n_points) { 

	if ((2*(histogram->n_points)-k-1 > 0) && (periodic)) {

	  l = 2*(histogram->n_points)-k-1;
	}

      } else {

	l = k;
      }

      if (l != -1) {

	Y = (histogram->points + l)->Y;
	sY = (histogram->points + l)->sY;

      } else {

	Y = sY = 0.0;
      }

      f = exp(-0.5*sqr(dX/sX));

      sumY += f*Y;
      sumSY2 += sqr(f*sY);
      sumf += f;
    }

    point2->X = point1->X;

    point2->Y = (sumf > 1.0E-10) ? sumY/sumf : point1->Y;

    point2->sY = (sumf > 1.0E-10) ? sqrt(sumSY2)/sumf : point1->sY;
  }
}



static void sharpen_histogram(HISTOGRAM *histogram,HISTOGRAM *sharpened_histogram,double sX) {

  int i;
  double f,z,f_opt,z_max;
  HISTOGRAM *smooth_histogram,*diff_histogram;
  HISTOGRAM_POINT *ph,*pd,*psh,*psm;

  smooth_histogram = clone_histogram(histogram);

  smoothe_histogram(histogram,smooth_histogram,sX);

  diff_histogram = clone_histogram(histogram);

  subtract_histograms(histogram,smooth_histogram,diff_histogram);

  z_max = 999.999;
  f_opt = -1.0;

  for (f=0.05;f<4.0;f+=0.05) {
    
    for (i=0,ph=histogram->points,pd=diff_histogram->points,psh=sharpened_histogram->points;i<histogram->n_points;i++,ph++,pd++,psh++) {

      psh->Y = ph->Y - f*(pd->Y);

      if (psh->Y < 0.0) {

	psh->Y = 0.0 ;
      }
    }

    smoothe_histogram(sharpened_histogram,smooth_histogram,sX);

    z = avg_histogram_point_zscore(smooth_histogram,histogram);

    if (z < z_max) {

      z_max = z;
      f_opt = f;
    }
  }

  if (f_opt > 0.0) {
    
    for (i=0,ph=histogram->points,pd=diff_histogram->points,psh=sharpened_histogram->points;i<histogram->n_points;i++,ph++,pd++,psh++) {

      psh->Y = ph->Y - f_opt*(pd->Y);

      if (psh->Y < 0.0) {

	psh->Y = 0.0 ;
      }
    }
  }

  free_histogram(smooth_histogram);
  free_histogram(diff_histogram);
}



static int remove_local_extreme(HISTOGRAM *histogram) {

  int i,n_points,best_i;
  double Y,sY,min_sY;
  HISTOGRAM_POINT *point,*best_point;
  HISTOGRAM_POINT *test;

  best_point = NULL;
  best_i = -1;

  min_sY = 9999.999;
  sY = 9999.999;

  n_points = histogram->n_points;

  for (i=0,point=histogram->points;i<n_points;i++,point++) {

    test = point;

    if (local_histogram_maximum(point,i,n_points)) {

      Y = point->Y;

      point->Y = point_intrapolate(point,i,n_points,&sY);

      if (!local_histogram_maximum(point,i,n_points)) {

	if (sY < min_sY) {

	  best_point = point;
	  best_i = i;

	  min_sY = sY;
	}
      }

      point->Y = Y;

    } else if (local_histogram_minimum(point,i,n_points)) {

      Y = point->Y;

      point->Y = point_intrapolate(point,i,n_points,&sY);

      if (!local_histogram_minimum(point,i,n_points)) {

	if (sY < min_sY) {

	  best_point = point;
	  best_i = i;

	  min_sY = sY;
	}
      }

      point->Y = Y;
    }
  }

  return(0);

  if (best_point) {

    best_point->Y = point_intrapolate(best_point,best_i,n_points,&sY);

    return(1);
  }

  return(0);
}



static double point_intrapolate(HISTOGRAM_POINT *point,int i,int n_points,double *sY) {

  double X,Y,Y1,Y2;

  X = point->X;

  if (i < 2) {

    *sY = 999.999;

    return(two_point_extrapolate(point+2,point+1,X));

  } else if (i > n_points - 3) {

    *sY = 999.999;

    return(two_point_extrapolate(point-2,point-1,X));
  }

  Y1 = two_point_extrapolate(point-2,point-1,X);
  Y2 = two_point_extrapolate(point+2,point+1,X);

  Y = 0.5*(Y1+Y2);

  *sY = fabs((Y1-Y2)/Y);

  return(Y);
}



static double two_point_extrapolate(HISTOGRAM_POINT *point1,HISTOGRAM_POINT *point2,double X) {

  double X1,X2,Y1,Y2;

  X1 = point1->X;
  X2 = point2->X;
  Y1 = point1->Y;
  Y2 = point2->Y;

  return(Y1+((X-X1)/(X2-X1))*(Y2-Y1));
}



static void subtract_histograms(HISTOGRAM *hist1,HISTOGRAM *hist2,HISTOGRAM *diff_hist) {

  int i;
  HISTOGRAM_POINT *point1,*point2,*diff_point;

  for (i=0,point1=hist1->points,point2=hist2->points,diff_point=diff_hist->points;i<hist1->n_points;i++,point1++,point2++,diff_point++) {

    diff_point->Y = point2->Y - point1->Y;
  }
}



static void log_histogram(HISTOGRAM *histogram) {

  int i;
  double minY,Y,sY;
  HISTOGRAM_POINT *point,*min_point;

  min_point = NULL;

  for (i=0,point=histogram->points;i<histogram->n_points;i++,point++) {

    if ((point->Y > 1.0E-10) && ((min_point == NULL) || (point->Y < min_point->Y))) {

      min_point = point;
    }
  }

  if (min_point == NULL) {

    error_fn("log_histogram: min_point undefined");
  }

  minY = 0.5*(min_point->Y);

  for (i=0,point=histogram->points;i<histogram->n_points;i++,point++) {

    Y = (point->Y > 1.0E-10) ? log(point->Y) : log(minY);

    sY = (point->Y > 1.0E-10) ? (point->sY)/(point->Y) : 1.0;

    point->Y = Y;
    point->sY = sY;
 }
}



static double avg_histogram_point_error(HISTOGRAM *histogram,HISTOGRAM *ref_histogram) {

  int i,n_points;
  double sum,w,sumw;
  HISTOGRAM_POINT *point,*ref_point;

  if (histogram->n_points == 0) {

    return(0.0);
  }

  sum = sumw = 0.0;
  n_points = 0;

  for (i=0,point=histogram->points,ref_point=ref_histogram->points;i<histogram->n_points;i++,point++,ref_point++) {

    if ((ref_point->sY > 1.0E-10) && (point->Y > 1.0E-10)) {

      w = sqr((ref_point->Y)/(ref_point->sY));

      sum += w*(point->sY/point->Y);

      sumw += w;
    }
  }

  if (sumw < 1.0E-10) {

    return(0.0);
  }

  return(sum/sumw);
}



static double avg_histogram_point_zscore(HISTOGRAM *histogram,HISTOGRAM *ref_histogram) {

  int i,n_points;
  double sum,w,sumw;
  HISTOGRAM_POINT *point,*ref_point;

  if (histogram->n_points == 0) {

    return(0.0);
  }

  sum = sumw = 0.0;
  n_points = 0;

  for (i=0,point=histogram->points,ref_point=ref_histogram->points;i<histogram->n_points;i++,point++,ref_point++) {

    if ((ref_point->sY > 1.0E-10) && (point->Y > 1.0E-10)) {

      w = sqr((ref_point->Y)/(ref_point->sY));

      sum += w*fabs(((point->Y) -(ref_point->Y))/ref_point->sY);

      sumw += w;
    }
  }

  if (sumw < 1.0E-10) {

    return(0.0);
  }

  return(sum/sumw);
}



static void free_histogram(HISTOGRAM *histogram) {

  if (histogram->n_points) {

    free(histogram->points);
  }

  free(histogram);
}



static void alloc_histograms(void) {

  if (n_alloc_histograms == 0) {

    n_alloc_histograms += histogram_batch_size;

    histograms = (HISTOGRAM*) calloc(n_alloc_histograms,sizeof(HISTOGRAM));

    if (histograms == NULL) {

      error_fn("alloc_histograms: out of memory allocating histograms");
    }

  } else if (n_histograms == n_alloc_histograms) {

    n_alloc_histograms += histogram_batch_size;

    histograms = (HISTOGRAM*) realloc(histograms,n_alloc_histograms*sizeof(HISTOGRAM));

    if (histograms == NULL) {

      error_fn("alloc_histograms: out of memory reallocating histograms");
    }
  }
}



static void alloc_histogram_points(HISTOGRAM *histogram) {

  if (histogram->n_alloc_points == 0) {

    histogram->n_alloc_points = histogram_point_batch_size;

    histogram->points = (HISTOGRAM_POINT*) calloc(histogram->n_alloc_points,sizeof(HISTOGRAM_POINT));

    if (histogram->points == NULL) {

      error_fn("alloc_histogram_points: out of memory allocating histogram_points");
    }

  } else if (histogram->n_alloc_points == histogram->n_points) {

    histogram->n_alloc_points += histogram_point_batch_size;

    histogram->points = (HISTOGRAM_POINT*) realloc(histogram->points,(histogram->n_alloc_points)*sizeof(HISTOGRAM_POINT));

    if (histogram->points == NULL) {

      error_fn("alloc_histogram_points: out of memory reallocating histogram_points");
    }
  }
}
