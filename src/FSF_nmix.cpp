// In order to do this I need to draw from a discrete random variable.
// I can use my own hacked version of this and see what happens.

#include "FSF_nmix.h"
#include <vector>
#include "RNG.h"
#include <math.h>
#include <stdio.h>

#ifdef USE_R
#include "R.h"
#include "Rmath.h"
#include <R_ext/Utils.h>
#endif

using std::vector;

void draw_indicators_logistic(int *r, double *z, double *lambda, int *n, 
			      double *w, double *s, int *J)
{
  vector<double> lpmf(*J);
  vector<double> pmf(*J);
  vector<double> cdf(*J);
  
  vector<double> lcoef(*J);
  for (int j = 0; j < *J; j++)
    lcoef[j] = log(w[j]) - log(s[j]);
  
  #ifdef USE_R
  GetRNGstate();
  #endif

  RNG rng;

  for (int i = 0; i < *n; i++) {
    #ifdef USE_R
    if (i%10==0) R_CheckUserInterrupt();
    #endif

    // Construct CDF.
    cdf[0] = 0.0;
    for (int j = 0; j < *J; j++) {
      lpmf[j] = lcoef[j] - 0.5 * pow((z[i] - log(lambda[i])) / s[j], 2);
      pmf[j] = exp(lpmf[j]);
      cdf[j] = pmf[j] + (j > 0 ? cdf[j-1] : 0.0);
    }
    // Draw from discrete.
    double U = rng.flat(0.0, cdf[*J-1]);
    int k = 0;
    while(cdf[k] < U)
      k++;
    r[i] = k+1; // k in C indexing, k+1 in R indexing.
  }

  #ifdef USE_R
  PutRNGstate();
  #endif

} // indicators_logistic

void draw_indicators_generic(int *r, double *res, int *n, 
			     double *w, double *m, double *s, int *J)
{
  // res_i = \ep_i, \ep_i \sim N(m_{r_i}, v_{r_i}).
  vector<double> lpmf(*J);
  vector<double> pmf(*J);
  vector<double> cdf(*J);
  
  vector<double> lcoef(*J);
  for (int j = 0; j < *J; j++)
    lcoef[j] = log(w[j]) - log(s[j]);
  
  #ifdef USE_R
  GetRNGstate();
  #endif

  RNG rng;

  for (int i = 0; i < *n; i++) {
    #ifdef USE_R
    if (i%10==0) R_CheckUserInterrupt();
    #endif

    // Construct CDF.
    cdf[0] = 0.0;
    for (int j = 0; j < *J; j++) {
      lpmf[j] = lcoef[j] - 0.5 * pow((res[i] - m[j]) / s[j], 2);
      pmf[j] = exp(lpmf[j]);
      cdf[j] = pmf[j] + (j > 0 ? cdf[j-1] : 0.0);
    }
    
    // Draw from discrete.
    double U = rng.flat(0.0, cdf[*J-1]);
    int k = 0;
    while(cdf[k] < U)
      k++;
    r[i] = k+1; // k in C indexing, k+1 in R indexing.
  }

  #ifdef USE_R
  PutRNGstate();
  #endif

} // indicators_generic

