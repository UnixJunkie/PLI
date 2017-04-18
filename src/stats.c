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


#define MAX_POLY_ORDER 10



static double power_fit_S(double*,double*,double*,int,double);



double polynomial(double x,double *c,int n) {

  int i;
  double y;

  y = c[0];

  for (i=1;i<=n;i++) {
    
    y += c[i]*pow(x,(double) i);
  }

  return(y);
}



int weighted_asymmetric_parabola_fit(double *x,double *y,double *s,int n,double xo,double yo,double *a,double *yfit) {

  int i,n1,n2;
  double xi,yi,wi,sumw,sumy,sum1x2,sum2x2,sum1x4,sum2x4,sum1x2y,sum2x2y;

  /*
  n1 = n2 = 0;
  sumw = sumy = sum1x2 = sum2x2 = sum1x4 = sum2x4 = sum1x2y = sum2x2y = 0.0;

  for (i=0;i<n;i++) {

    wi = (1.0/sqr(s[i]));

    xi = x[i] - xo;

    if (xi <= 0.0) {

      sum1x2 += wi*sqr(xi);
      sum1x4 += wi*pow(xi,4.0);
      sum1x2y += wi*sqr(xi)*y[i];

      n1++;

    } else {

      sum2x2 += wi*sqr(xi);
      sum2x4 += wi*pow(xi,4.0);
      sum2x2y += wi*sqr(xi)*y[i];

      n2++;
    }

    sumw += wi;
    sumy += wi*y[i];
  }

  if ((n1 < 10) || (n2 < 10)) {

    return(1);
  }

  a[0] = ((sumy*sum1x4*sum2x4) - (sum1x2*sum2x4*sum1x2y) - (sum2x2*sum1x4*sum2x2y)) / ((sumw*sum1x4*sum2x4) - (sum2x4*sqr(sum1x2)) - (sum1x4*sqr(sum2x2)));

  a[1] = (sum1x2y - (a[0]*sum1x2)) / sum1x4;
  a[2] = (sum2x2y - (a[0]*sum2x2)) / sum2x4;

  printf("test a[0] = %10.4lf, a[1] = %10.4lf, a[2] = %10.4lf\n",a[0],a[1],a[2]);

  for (i=0;i<n;i++) {

    xi = x[i] - xo;

    yfit[i] = (xi <= 0.0) ? a[0] + a[1]*sqr(xi) : a[0] + a[2]*sqr(xi);

    //printf("test %10.4lf %10.4lf\n",y[i],yfit[i]);
  }
*/
  
  n1 = n2 = 0;
  sum1x4 = sum2x4 = sum1x2y = sum2x2y = 0.0;

  for (i=0;i<n;i++) {

   wi = (1.0/sqr(s[i]));

    xi = x[i] - xo;
    yi = y[i] - yo;

    if (xi <= 0.0) {

      sum1x4 += wi*pow(xi,4.0);
      sum1x2y += wi*sqr(xi)*yi;

      n1++;

    } else {

      sum2x4 += wi*pow(xi,4.0);
      sum2x2y += wi*sqr(xi)*yi;

      n2++;
   }
  }

  if ((n1 < 3) || (n2 < 3)) {
  
    return(1);
  }

  a[0] = sum1x2y/sum1x4;
  a[1] = sum2x2y/sum2x4;

  for (i=0;i<n;i++) {

    xi = x[i] - xo;

    yfit[i] = (xi <= 0.0) ? yo + a[0]*sqr(xi) : yo + a[1]*sqr(xi);
  }

  return(0);
}


int weighted_polynomial_fit(double *x,double *y,double *s,int n,int order,double *a,double *yfit) {

  int i,j,k;
  double **mx,**mc;
  double sumx[2*(MAX_POLY_ORDER+1)],sumy[MAX_POLY_ORDER+1];
  double w[1000],sumw,detmx;

  if ((order < 1) || (order > MAX_POLY_ORDER)) {

    return(1);
  }

  mx = alloc_2d_fmatrix(order+1,order+1);
  mc = alloc_2d_fmatrix(order+1,order+1);

  // calculate sums:

  set_double_array(sumx,(order+1)*2,0.0);
  set_double_array(sumy,order+1,0.0);

  sumw = 0.0;

  for (i=0;i<n;i++) {

    // TODO: weighting the points seems to make the polynomial unstable
    // not sure if it is to do with the implementation?

    w[i] = (1.0/sqr(s[i]));

    //w[i] = 1.0;

    sumw += w[i];
  }

  for (i=0;i<n;i++) {

    w[i] = (w[i]/sumw);
  }

  for (i=0;i<n;i++) {

    for (j=0;j<=2*order;j++) {

      sumx[j] += (w[i]*pow(x[i],(double) j));
    }

    for (j=0;j<=order;j++) {

      sumy[j] += (w[i]*pow(x[i],(double) j)*y[i]);
    }
  }

  // generate matrix:

  for (i=0;i<=order;i++) {

    for (j=0;j<=order;j++) {

      mx[i][j] = sumx[i+j];
    }
  }

  // apply Cramer's rule:

  detmx = matrix_determinant(mx,order+1);

  for (k=0;k<=order;k++) {

    for (i=0;i<=order;i++) {
      
      for (j=0;j<=order;j++) {

	mc[i][j] = (j == k) ? sumy[i] : mx[i][j];
      }
    }
   
    a[k] = matrix_determinant(mc,order+1)/detmx;
  }

  for (i=0;i<n;i++) {

    yfit[i] = a[0];

    for (j=1;j<=order;j++) {

      yfit[i] += a[j]*pow(x[i],(double) j);
    }
  }

  return(0);
}



int weighted_power_fit(double *x,double *y,double *s,int n,double *a,double *yfit) {

  int i,nl;
  double xl[100],yl[100],sl[100],yfitl[100],c[2],S1,S2,S3,b1,b2,b3,f,db;
  double w,sumwxby,sumwx2b;

  // estimate coefficients from linear regression first:

  nl = 0;

  for (i=0;i<n;i++) {

    if (x[i] < 0.0) {

      return(1);
    }

    if (y[i] > 0.0) {

      xl[nl] = log(x[i]);
      yl[nl] = log(y[i]);

      sl[nl] = s[i]/y[i];

      nl++;
    }
  }

  if (nl < 2) {

    return(2);
  }

  if (weighted_polynomial_fit(xl,yl,sl,nl,1,c,yfitl)) {

    return(3);
  }

  // bracket minimum for S(b):

  f = 1.1;
  b2 = c[1];
  S2 = power_fit_S(x,y,s,n,b2);

  if (S2 > 0.0) {

    // optimum b is lower:

    b3 = b2;
    S3 = S2;

    do {

      b1 = (b3 < 0.0) ? b3*f : b3/f;

      S1 = power_fit_S(x,y,s,n,b1);

      f *= 1.1;
   
    } while ((S1 > 0.0) && (f < 4.0));

  } else {

    // optimum b is higher:

    b1 = b2;
    S1 = S2;

    do {

      b3 = (b1 < 0.0) ? b1/f : b1*f;

      S3 = power_fit_S(x,y,s,n,b3);

      f *= 1.1;
   
    } while ((S3 < 0.0) && (f < 4.0));
  }

  b2 = 0.5*(b1+b3);
  S2 = power_fit_S(x,y,s,n,b2);

  if (((S2/S1) > 0.0) && ((S2/S3) > 0.0)) {

    return(4);
  }

  // optimise b:
    
  do {
    
    if ((S1/S2) < 0.0) {
      
      b3 = b2;
      S3 = S2;
      
    } else {
      
      b1 = b2;
      S1 = S2;
    }
    
    b2 = 0.5*(b1+b3);
    
    db = b2-b1;
    
    S2 = power_fit_S(x,y,s,n,b2);
    
  } while (fabs(db/b2) > 1.0E-10);
 
  a[1] = b2;

  // calculate a:

  sumwxby = sumwx2b = 0.0;

  for (i=0;i<n;i++) {

    w = (1.0/sqr(s[i]));

    sumwxby += (w*pow(x[i],a[1])*y[i]);
    sumwx2b += (w*pow(x[i],2.0*a[1]));
  }

  a[0] = sumwxby/sumwx2b;

  for (i=0;i<n;i++) {
    
    yfit[i] = a[0]*pow(x[i],a[1]);
  }
 
  if (fabs(S2) > 0.1) {

    for (i=0;i<n;i++) {
            
      printf("powerfit %10.4lf %10.4lf %10.4lf\n",x[i],y[i],yfit[i]);
    }
    
    exit(0);

    return(5);
  }

  return(0);
}



static double power_fit_S(double *x,double *y,double *s,int n,double b) {

  int i;
  double w,sumw,sumxby,sumlnxx2b,sumx2b,sumlnxxby;

  sumw = sumxby =  sumlnxx2b = sumx2b = sumlnxxby = 0.0;

  for (i=0;i<n;i++) {

    w = (1.0/sqr(s[i]));

    sumw += w;

    sumxby += (w*pow(x[i],b)*y[i]);
    sumlnxx2b += (w*log(x[i])*pow(x[i],2.0*b));
    sumx2b += (w*pow(x[i],2.0*b));
    sumlnxxby += (w*log(x[i])*pow(x[i],b)*y[i]);
  }

  return(((sumxby*sumlnxx2b)-(sumx2b*sumlnxxby))/sumw);
}



int weighted_linear_least_squares(double *x,double *y,double *s,int n,double *a,double *b,double *yfit) {

  int i;
  double sumw,sumwx,sumwy,sumwxx,sumwxy,Xw,Yw,*w;

  w = (double*) calloc(n,sizeof(double));

  if (w == NULL) {

    error_fn("%s: out of memory allocating weights",__func__);
  }

  for (i=0;i<n;i++) {

    w[i] = 1.0/(sqr(s[i]));
  }

  sumw = sumwx = sumwy = 0.0;

  for (i=0;i<n;i++) {

    sumw += w[i];
 
    sumwx += (w[i]*x[i]);
    sumwy += (w[i]*y[i]);
  }

  if (sumw < 1.0E-30) {

    free(w);

    return(1);
  }

  Xw = sumwx/sumw;
  Yw = sumwy/sumw;

  sumwxx = sumwxy = 0.0;

  for (i=0;i<n;i++) {

    sumwxx += (w[i]*sqr(x[i]-Xw));
    sumwxy += (w[i]*(x[i]-Xw)*(y[i]-Yw));
  }

  if (sumwxx < 1.0E-30) {

    free(w);

    return(2);
  }

  *a = sumwxy/sumwxx;
  *b = Yw - (*a)*Xw;

  for (i=0;i<n;i++) {

    yfit[i] = ((*a)*x[i]) + (*b);
  }

  free(w);

  return(0);
}



int weighted_LJ_Eo_fit(double *x,double *y,double *s,int n,int order,double *Eo,double *yfit) {

  int i;
  double w,sumwx2,sumwx3,sumwx4,sumwx1y,sumwx2y;

  sumwx2 = sumwx3= sumwx4= sumwx1y = sumwx2y = 0.0;
  
  for (i=0;i<n;i++) {

    w = 1.0/(sqr(s[i]));

    sumwx2 += (w*pow(x[i],-2.0*((double) order)));
    sumwx3 += (w*pow(x[i],-3.0*((double) order)));
    sumwx4 += (w*pow(x[i],-4.0*((double) order)));

    sumwx1y += (w*y[i]*pow(x[i],-1.0*((double) order)));
    sumwx2y += (w*y[i]*pow(x[i],-2.0*((double) order)));
  }

  *Eo = ((2.0*sumwx1y) - sumwx2y)/(sumwx4+(4.0*sumwx2)-(4.0*sumwx3));

  for (i=0;i<n;i++) {

    yfit[i] = (*Eo)*(pow(x[i],2.0*((double) order)) -2.0*pow(x[i],(double) order));

    //printf("%10.4lf %10.4lf %10.4lf %10.4lf\n",x[i],y[i],s[i],yfit[i]);
  }

  return(0);
}



int weighted_LJ_fit(double *x,double *y,double *w,int n,int order,double *a,double *b,double *yfit) {

  int i;
  double sumwx2,sumwx3,sumwx4,sumwx1y,sumwx2y;

  sumwx2 = sumwx3= sumwx4= sumwx1y = sumwx2y = 0.0;
  
  for (i=0;i<n;i++) {
 
    sumwx2 += (w[i]*pow(x[i],-2.0*((double) order)));
    sumwx3 += (w[i]*pow(x[i],-3.0*((double) order)));
    sumwx4 += (w[i]*pow(x[i],-4.0*((double) order)));

    sumwx1y += (w[i]*y[i]*pow(x[i],-1.0*((double) order)));
    sumwx2y += (w[i]*y[i]*pow(x[i],-2.0*((double) order)));
  }

  *a = ((sumwx2y*sumwx2) - (sumwx1y*sumwx3))/((sumwx4*sumwx2) - (sumwx3*sumwx3));

  *b = (((*a)*sumwx4) - sumwx2y)/sumwx3;

  for (i=0;i<n;i++) {

    yfit[i] = ((*a)*pow(x[i],-2.0*((double) order))) - ((*b)*pow(x[i],-1.0*((double) order)));
  }

  return(0);
}



double normalised_fit_error(double *y,double *yfit,double *sy,int n) {

  int i;
  double sumE2;

  sumE2 = 0.0;

  for (i=0;i<n;i++) {

    sumE2 += sqr((yfit[i]-y[i])/sy[i]);
  }

  return(sqrt((sumE2)/((double) n)));
}


double weighted_fit_error( double *y,double *yfit,double *sy,int n) {

  int i;
  double sumE2,sumw;

  sumE2 = sumw = 0.0;

  for (i=0;i<n;i++) {

    sumE2 += sqr((yfit[i]-y[i])/sy[i]);

    sumw += sqr(1.0/sy[i]);
  }

  return(sqrt(sumE2/sumw));
}



