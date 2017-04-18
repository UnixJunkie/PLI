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



static int frand_seeded = 0;
static int irand_seeded = 0;



static void seed_frand(void);
static void seed_irand(void);


void use_fixed_seed (void) {
  frand_seeded = 1;
  irand_seeded = 1;
}

double frand(void) {

  if (!frand_seeded) {

    seed_frand();
  }

  return(drand48());
}



int irand(void) {

  if (!irand_seeded) {

    seed_irand();
  }

  return(rand());
}



static void seed_frand(void) {

  struct timeval tp;
  struct timezone tzp;
  long lseed;

  do {

    gettimeofday(&tp, &tzp);

  } while (tp.tv_usec == 0);

  lseed = tp.tv_usec;
  srand48(lseed);

  frand_seeded = 1;
}

double uniform_rand(double min, double max) {
  return(min + frand() * (max-min));
}

int uniform_rand_int(int min, int max){
  return(irand() % (max + 1 - min) + min);
}

double normal_rand(double mean, double sd) {
  double U1, U2, W, multiplier;
  static double X1, X2;
  static int call = 0;

  if (call == 1) {
    call = !call;
    return (mean + sd * (double) X2);
  }

  do {
    U1 = -1 + (frand()) * 2;
    U2 = -1 + (frand()) * 2;
    W = pow (U1, 2) + pow (U2, 2);
  } while (W >= 1 || W == 0);

  multiplier = sqrt((-2 * log (W)) / W);
  X1 = U1 * multiplier;
  X2 = U2 * multiplier;

  call = !call;
  return (mean + sd * (double) X1);
}


static void seed_irand(void) {

  struct timeval tp;
  struct timezone tzp;
  long lseed;

  do {

    gettimeofday(&tp, &tzp);

  } while (tp.tv_usec == 0);

  lseed = tp.tv_usec;
  srand(lseed);

  irand_seeded = 1;
}



double distance(double *v1,double *v2) {

  return(sqrt(sqr(v2[0]-v1[0]) + sqr(v2[1]-v1[1]) + sqr(v2[2]-v1[2])));
}



double sqr_distance(double *v1,double *v2) {

  return(sqr(v2[0]-v1[0]) + sqr(v2[1]-v1[1]) + sqr(v2[2]-v1[2]));
}



int points_within_distance(double *p1,double *p2,double dist) {

  double dx,dy,dz;

  dx = fabs(p2[0] - p1[0]);

  if (dx > dist)
    return(0);

  dy = fabs(p2[1] - p1[1]);

  if (dy > dist)
    return(0);

  dz = fabs(p2[2] - p1[2]);

  if (dz > dist)
    return(0);

  if ((dx*dx + dy*dy + dz*dz) < (dist*dist))
    return(1);

  return(0);
}



void copy_vector(double *v1,double *v2) {

  v2[0] = v1[0];
  v2[1] = v1[1];
  v2[2] = v1[2];

  v2[3] = 1.0;
}



void set_vector(double *v,double x,double y,double z) {

  v[0] = x;
  v[1] = y;
  v[2] = z;

  v[3] = 1.0;
}



void null_vector(double *v) {

  set_vector(v,0.0,0.0,0.0);
}



void calc_vector(double *p1,double *p2,double *v) {

  int i;

  for (i=0;i<3;i++) {

    v[i] = p2[i] - p1[i];
  }

  v[3] = 1.0;
}



void invert_vector(double *v1,double *v2) {

  int i;

  for (i=0;i<3;i++) {

    v2[i] = -v1[i];
  }

  v2[3] = 1.0;
}



double vector_length(double *v) {

  return(sqrt(sqr(v[0])+sqr(v[1])+sqr(v[2])));
}



double sqr_vector_length(double *v) {

  return(sqr(v[0])+sqr(v[1])+sqr(v[2]));
}



double dotproduct(double *v1,double *v2) { 

  return(v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]);
}



void calc_crossproduct(double* v1,double* v2,double *v3) {

  v3[0] = v1[1]*v2[2] - v1[2]*v2[1];
  v3[1] = v1[2]*v2[0] - v1[0]*v2[2];
  v3[2] = v1[0]*v2[1] - v1[1]*v2[0];

  v3[3] = 1.0;
}




double vector_angle(double *v1,double *v2) {
  double a, dp;
  a = vector_length(v1)*vector_length(v2);
  if (fabs(a) < 1.0E-6) {
    return (0.0);
  }

  dp = dotproduct(v1,v2);
  if (fabs(dp/a + 1.0) < 1.0E-6) {   // dp == -1
    return(180.0);
  }
  return((180.0/PI)*acos(dp/a));
}


void scale_vector(double *v,double len) {
 
  double a;
  int i;

  a = len/vector_length(v);

  for (i=0; i<3;i++) {

    v[i] *= a;
  }
}



void sum_vector(double *v1,double *v2,double *v3) {

  int i;

  for (i=0;i<3;i++) {

    v3[i] = v1[i] + v2[i];
  }

  v3[3] = 1.0;
}


void shift_vector(double *v,double *shift) {

  int i;

  for (i=0;i<3;i++) {

    v[i] = v[i] + shift[i];
  }

  v[3] = 1.0;
}


void transform_vector(double *v1,double m[4][4]) {
 
  int i,j;
  double v2[4];

  for (i=0;i<4;i++) {

    v2[i] = 0.0;

    for (j=0;j<4;j++) {

      v2[i] += v1[j]*m[i][j];
    }
  }

  for (i=0;i<3;i++) {

    v1[i] = v2[i];
  }
}


void calc_transformed_vector(double *v1,double m[4][4],double *v2) {
 
  int i,j;

  for (i=0;i<4;i++) {

    v2[i] = 0.0;

    for (j=0;j<4;j++) {

      v2[i] += v1[j]*m[i][j];
    }
  }
}



void write_matrix(char *name,double m[4][4]) {

  int i,j;

  printf("matrix %s\n",name);

  for (i=0;i<4;i++) {

    for (j=0;j<4;j++) {

      printf("[%1d%1d]%10.4lf",i,j,m[i][j]);
    }

    printf("\n");
  }
}


void unit_matrix(double m[4][4]) {

  int i,j;

  for (i=0;i<4;i++) {

    for (j=0;j<4;j++) {

      m[i][j] = (i == j) ? 1.0 : 0.0;
    }
  }
}

void scalar_matrix(double value, double m[4][4]) {

  int i,j;

  for (i=0;i<4;i++) {

    for (j=0;j<4;j++) {

      m[i][j] = (i == j) ? value : 0.0;
    }
  }
}

void copy_matrix(double m1[4][4], double m2[4][4]) {
  int i, j;
  for (i=0; i<4; i++) {
    for (j=0; j<4; j++) {
      m2[i][j] = m1[i][j];
    }
  }
}




void translation_matrix(double *v,double m[4][4]) {

  int i,j;

  for (i=0;i<4;i++) {

    for (j=0;j<3;j++) {

      m[i][j] = (i == j) ? 1.0 : 0.0;
    }

    m[i][3] = v[i];
  }
}



void euler_matrix(double ax,double ay,double az,double m[4][4]) {
  
  ax *= (PI/180.0);
  ay *= (PI/180.0);
  az *= (PI/180.0);

  m[0][0] =  cos(ay)*cos(az);
  m[0][1] = -cos(ax)*sin(az) + sin(ax)*sin(ay)*cos(az);
  m[0][2] =  sin(ax)*sin(az) + cos(ax)*sin(ay)*cos(az);
  m[0][3] = 0.0;

  m[1][0] =  cos(ay)*sin(az);
  m[1][1] =  cos(ax)*cos(az) + sin(ax)*sin(ay)*sin(az);
  m[1][2] = -sin(ax)*cos(az) + cos(ax)*sin(ay)*sin(az);
  m[1][3] = 0.0;

  m[2][0] = -sin(ay);
  m[2][1] =  sin(ax)*cos(ay);
  m[2][2] =  cos(ax)*cos(ay);
  m[2][3] = 0.0;

  m[3][0] = 0.0;
  m[3][1] = 0.0;
  m[3][2] = 0.0;
  m[3][3] = 1.0;
}



void rotation_matrix(double *p1,double *p2,double q,double m[4][4]) {
  // FIXME: when p1==p2!

  double axis[4],a,b,c,u,v,w,u2,v2,w2,uv,uw,vw,sinq,cosq;

  q *= (PI/180.0);

  calc_vector(p1,p2,axis);

  scale_vector(axis,1.0);

  a = p1[0]; b = p1[1]; c = p1[2];
  u = axis[0]; v = axis[1]; w = axis[2];

  u2 = sqr(u); v2 = sqr(v); w2 = sqr(w);

  uv = u*v; uw = u*w; vw = v*w;

  sinq = sin(q); cosq = cos(q);

  m[0][0] = u2 + (v2+w2)*cosq;
  m[0][1] = uv*(1-cosq) - w*sinq;
  m[0][2] = uw*(1-cosq) + v*sinq;
  m[0][3] = (a*(v2+w2)-u*(b*v+c*w))*(1-cosq) + (b*w-c*v)*sinq;

  m[1][0] = uv*(1-cosq) + w*sinq;
  m[1][1] = v2 + (u2+w2)*cosq;
  m[1][2] = vw*(1-cosq) - u*sinq;
  m[1][3] = (b*(u2+w2)-v*(a*u+c*w))*(1-cosq) + (c*u-a*w)*sinq;

  m[2][0] = uw*(1-cosq) - v*sinq;
  m[2][1] = vw*(1-cosq) + u*sinq;
  m[2][2] = w2 + (u2+v2)*cosq;
  m[2][3] = (c*(u2+v2)-w*(a*u+b*v))*(1-cosq) + (a*v-b*u)*sinq;

  m[3][0] = 0.0;
  m[3][1] = 0.0;
  m[3][2] = 0.0;
  m[3][3] = 1.0;  
}



void calc_matrix_product(double m1[4][4],double m2[4][4],double m3[4][4]) {

  int i,j,k;
  
  for (i=0;i<4;i++) {

    for (j=0;j<4;j++) {

      m3[i][j] = 0.0;

      for (k=0;k<4;k++) {

	m3[i][j] += m2[i][k]*m1[k][j];
      }
    }
  }
}



double matrix_determinant(double **m,int n) {

  int i,j,k,subi,subj,subn,sign;
  double **subm,det;

  if (n == 2) {

    return((m[0][0]*m[1][1])-(m[0][1]*m[1][0]));
  }

  subn = n - 1;

  subm = alloc_2d_fmatrix(subn,subn);

  sign = 1;

  det = 0.0;

  for(k=0;k<n;k++) {

    for (i=1,subi=0;i<n;i++,subi++) {

      for(j=0,subj=0;j<n;j++) {

	if (j != k) {

	  subm[subi][subj] = m[i][j];

	  subj++;
	}
      }
    }

    det += sign*m[0][k]*matrix_determinant(subm,subn);

    sign *= -1;
  }

  free(subm);

  return(det);
}



double three_point_angle(double *v1,double *v2,double *v3) {

  double v21[4],v23[4];

  calc_vector(v2,v1,v21);
  calc_vector(v2,v3,v23);

  return(vector_angle(v21,v23));
}



double dihedral_angle(double *v1,double *v2,double *v3,double *v4) {

  int i;
  double v12[4],v23[4],v34[4],vn1[4],vn2[4],vtp[4];
  double dp1,dp2,l12,ratio,angle;

  for (i=0;i<3;i++) {

    v12[i] = v2[i] - v1[i];
    v23[i] = v3[i] - v2[i];
    v34[i] = v4[i] - v3[i];
  }
    
  calc_crossproduct(v12,v23,vn1);
  calc_crossproduct(v23,v34,vn2);
  calc_crossproduct(vn1,vn2,vtp);

  dp1 = dotproduct(vn1,vn2);
  dp2 = dotproduct(vtp,v23);

  l12 = vector_length(vn1)*vector_length(vn2);

  ratio = dp1/l12;

  if (ratio > 1.0) {

    ratio = 1.0;

  } else if (ratio < -1.0) {

    ratio = -1.0;
  }

  angle = (180.0/PI)*acos(ratio);

  if (dp2 < 0.0) {

    return(-angle);
  }

  return(angle);
}



double triangle_area(double *v1,double *v2,double *v3) {

  double v12[4],v13[4];

  calc_vector(v1,v2,v12);
  calc_vector(v1,v3,v13);

  return(0.5*sqrt(sqr_vector_length(v12)*sqr_vector_length(v13) - sqr(dotproduct(v12,v13))));
}



int solve_line(double x1,double y1,double x2,double y2,double *a,double *b) {

  double denom;

  denom = x2 - x1;

  if (fabs(denom) < 1.0E-30) {

    return(1);
  }

  *a = (y2-y1)/denom;

  *b = y1 - (*a)*x1;

  return(0);
}



int line_line_intersection(double a1,double b1,double a2,double b2,double *x,double *y) {

  double da;

  da = a2 - a1;

  if (fabs(da) < 1.0E-30) {

    return(1);
  }

  *x = (b1-b2)/da;

  *y = a1*(*x) + b1;

  return(0);
}



int solve_parabola(double x1,double y1,double x2,double y2,double x3,double y3,double *a,double *b,double *c) {

  double denom;

  denom = (x1 - x2) * (x1 - x3) * (x2 - x3);

  if (fabs(denom) < 1.0E-30) {

    return(1);
  }

  *a = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom;
  *b = (x3*x3 * (y1 - y2) + x2*x2 * (y3 - y1) + x1*x1 * (y2 - y3)) / denom;
  *c = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 + x1 * x2 * (x1 - x2) * y3) / denom;

  return(0);
}



int calc_parabola_vertex(double x1,double y1,double x2,double y2,double x3,double y3,double *xv,double *yv) {

  double a,b,c;

  if (solve_parabola(x1,y1,x2,y2,x3,y3,&a,&b,&c)) {

    return(1);
  }

  if (fabs(a) < 1.0E-30) {

    return(2);
  }

  *xv = -b/(2*a);
  *yv = c - (b*b)/(4*a);

  return(0);
}



int calc_parabola_x_intercepts(double a,double b,double c,double *x1,double *x2) {

  double F,sF;

  F = sqr(b) - 4.0*a*c;

  if (F < 0.0) {

    return(0);
  }

  sF = sqrt(F);

  *x1 = (-b - sF)/(2.0*a);
  *x2 = (-b + sF)/(2.0*a);

  return(2);
}



long int factorial(int n) {

  int i;
  long int f;

  if (n < 2) {

    return(1);
  }

  f = 1;

  for (i=2;i<=n;i++) {

    f *= i;
  }

  return(f);
}



double ramp_function(double x,double x1,double x2,double y1,double y2) {

  if (x < x1) {

    return(y1);

  } else if (x > x2) {

    return(y2);

  }

  return(y1 + (x-x1)*((y2-y1)/(x2-x1)));
}



double linear_interpolate(double v1, double v2, double d) {

  return (v1*(1.0-d) + v2*d);
}
