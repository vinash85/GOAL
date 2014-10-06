////////////////////////////////////////////////////////////////////////////////

// Copyright 2012 Nick Polson, James Scott, and Jesse Windle.

// This file is part of BayesLogit.

// BayesLogit is free software: you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or any later version.

// BayesLogit is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with
// BayesLogit.  If not, see <http://www.gnu.org/licenses/>.

////////////////////////////////////////////////////////////////////////////////

#ifdef USE_R
#include "R.h"
#include "Rmath.h"
#endif

#include "LogitWrapper.h"
#include "Logit.h"
#include "MultLogit.h"
#include "RNG.h"
#include "PolyaGamma.h"
#include "PolyaGammaAlt.h"
#include "PolyaGammaSP.h"
#include <exception>
#include <stdio.h>

////////////////////////////////////////////////////////////////////////////////
				// PolyaGamma //
////////////////////////////////////////////////////////////////////////////////

void rpg_gamma(double *x, double *n, double *z, int *num, int *trunc)
{
  RNG r;
  PolyaGamma pg(*trunc);

  #ifdef USE_R
  GetRNGstate();
  #endif

  for(int i=0; i < *num; ++i){
    #ifdef USE_R
      if (i % 1000 == 0) R_CheckUserInterrupt();
    #endif
      if (n[i]!=0.0) 
	x[i] = pg.draw_sum_of_gammas(n[i], z[i], r);
      else 
	x[i] = 0.0;

  }

  #ifdef USE_R
  PutRNGstate();
  #endif
} // rpg



void rpg_devroye(double *x, int *n, double *z, int *num)
{
  RNG r;
  PolyaGamma pg(1);

  #ifdef USE_R
  GetRNGstate();
  #endif

  for(int i=0; i < *num; ++i){
    if (n[i]!=0)
      x[i] = pg.draw(n[i], z[i], r);
    else
      x[i] = 0.0;
  }

  #ifdef USE_R
  PutRNGstate();
  #endif
} // rpg

void rpg_alt(double *x, double *h, double *z, int* num)
{
  RNG r;
  PolyaGammaAlt pg;

  #ifdef USE_R
  GetRNGstate();
  #endif

  for(int i=0; i < *num; ++i){
    if (h[i]!=0)
      x[i] = pg.draw(h[i], z[i], r);
    else
      x[i] = 0.0;
  }

  #ifdef USE_R
  PutRNGstate();
  #endif
}

void rpg_sp(double *x, double *h, double *z, int* num, int *iter)
{
  RNG r;
  PolyaGammaSP pg;
  
  #ifdef USE_R
  GetRNGstate();
  #endif

  for(int i=0; i < *num; ++i){
    if (h[i]!=0)
      iter[i] = pg.draw(x[i], h[i], z[i], r);
    else
      x[i] = 0.0;
  }

  #ifdef USE_R
  PutRNGstate();
  #endif
}

void rpg_hybrid(double *x, double *h, double *z, int* num)
{
  RNG r;
  PolyaGamma dv;
  PolyaGammaAlt alt;
  PolyaGammaSP sp;
  
  #ifdef USE_R
  GetRNGstate();
  #endif

  for(int i=0; i < *num; ++i){
    double b = h[i];
    if (b > 170) {
      double m = dv.pg_m1(b,z[i]);
      double v = dv.pg_m2(b,z[i]) - m*m;
      x[i] = r.norm(m, sqrt(v));
    }
    else if (b > 13) {
      sp.draw(x[i], b, z[i], r);
    }
    else if (b==1 || b==2) {
      x[i] = dv.draw((int)b, z[i], r);
    }
    else if (b > 0) {
      x[i] = alt.draw(b, z[i], r);
    }
    else {
      x[i] = 0.0;
    }
  }

  #ifdef USE_R
  PutRNGstate();
  #endif
}


////////////////////////////////////////////////////////////////////////////////
			   // POSTERIOR INFERENCE //
////////////////////////////////////////////////////////////////////////////////

// Posterior by Gibbs.
//------------------------------------------------------------------------------
void gibbs(double *wp, double *betap,                            // Posterior
	   double *yp, double *tXp, double *np,                  // Data
	   double *m0p, double *P0p,                             // Prior
	   int *N, int *P,                                       // Dim
	   int *samp, int *burn)                                 // MCMC
{

  // Set up data.
  Matrix y (yp,  *N,  1);
  Matrix tX(tXp, *P, *N);
  Matrix n (np,  *N,  1);
  Matrix m0(m0p, *P,  1);
  Matrix P0(P0p, *P, *P);

  // Declare posteriors.
  Matrix w, beta;

  // Random number generator.
  RNG r;

  #ifdef USE_R
  GetRNGstate();
  #endif

  // Logit Gibbs
  try{
    Logit logit(y, tX, n);
    logit.set_prior(m0, P0);
    // logit.compress();

    // Set the correct dimensions after combining data.
    w.resize(logit.get_N(), 1, *samp);
    beta.resize(logit.get_P(), 1, *samp);
    MatrixFrame w_mf    (wp   , w.rows()   , w.cols()   , w.mats());
    MatrixFrame beta_mf (betap, beta.rows(), beta.cols(), beta.mats());

    // Copy values to test code.  Must be using known value with correct dim.
    // w.copy(w_mf);
    // beta.copy(beta_mf);

    // Run gibbs.
    logit.gibbs(w, beta, *samp, *burn, r);

    // Copy values to return.
    w_mf.copy(w);
    beta_mf.copy(beta);

    // Adjust for combined data.
    *N = w.rows();
  }
  catch (std::exception& e) {
    Rprintf("Error: %s\n", e.what());
    Rprintf("Aborting Gibbs sampler.\n");
  }

  #ifdef USE_R
  PutRNGstate();
  #endif
} // gibbs

// Posterior Mode by EM
//------------------------------------------------------------------------------
void EM(double *betap,
	double *yp, double *tXp, double *np,
	int *Np, int *Pp,
	double *tolp, int *max_iterp)
{
  int N = *Np;
  int P = *Pp;

  // Set up data.
  Matrix y (yp,  N, 1);
  Matrix tX(tXp, P, N);
  Matrix n (np,  N, 1);

  // Declare posteriors.
  Matrix beta(P);

  // Logit EM
  try{
    Logit logit(y, tX, n);
    // logit.compress();
    *max_iterp = logit.EM(beta, *tolp, *max_iterp);

    // Copy.
    MatrixFrame beta_mf (betap, beta.rows(), beta.cols(), beta.mats());
    beta_mf.copy(beta);

  }
  catch (std::exception& e) {
    Rprintf("Error: %s\n", e.what());
    Rprintf("Aborting EM.\n");
  }

  #ifdef USE_R
  PutRNGstate();
  #endif
}

////////////////////////////////////////////////////////////////////////////////

// combine_data
//------------------------------------------------------------------------------
void combine(double *yp, double *tXp, double *np,                  // Data
	     int *N, int *P)
{

  // Set up data.
  Matrix y (yp,  *N,  1);
  Matrix tX(tXp, *P, *N);
  Matrix n (np,  *N,  1);

  // Logit Gibbs
  try{
    Logit logit(y, tX, n);
    logit.compress();
    logit.get_data(y, tX, n);

    // Copy.
    MatrixFrame y_mf    (yp ,  y.rows(),  y.cols(),  y.mats());
    MatrixFrame tX_mf   (tXp, tX.rows(), tX.cols(), tX.mats());
    MatrixFrame n_mf    (np ,  n.rows(),  n.cols(),  n.mats());

    y_mf.copy(y);
    tX_mf.copy(tX);
    n_mf.copy(n);

    *N = tX.cols();
  }
  catch (std::exception& e) {
    Rprintf("Error: %s\n", e.what());
    Rprintf("Aborting combine.\n");
  }

} // combine_data

////////////////////////////////////////////////////////////////////////////////
			    // Multinomial Gibbs //
////////////////////////////////////////////////////////////////////////////////

void mult_gibbs(double *wp, double *betap,
		double *typ, double *tXp, double *np,
		double *m0p, double *P0p,
		int *N, int *P, int *J,
		int *samp, int *burn)
{
  int unk = *J-1; // number unknown

  // Set up data.
  Matrix ty(typ, unk, *N);
  Matrix tX(tXp, *P , *N);
  Matrix n (np,  *N ,  1);

  Matrix m0(m0p, *P,  1, unk);
  Matrix P0(P0p, *P, *P, unk);

  // Declare posteriors.
  Matrix w, beta;

  // Random number generator.
  RNG r;

  #ifdef USE_R
  GetRNGstate();
  #endif

  // Multinomial Logistic Regression via Gibbs
  try{
    MultLogit logit(ty, tX, n);

    // Set the correct dimensions after combining data.
    w.resize   (logit.get_N(), unk, *samp);
    beta.resize(logit.get_P(), unk, *samp);
    MatrixFrame w_mf    (wp   , w.rows()   , w.cols()   , w.mats());
    MatrixFrame beta_mf (betap, beta.rows(), beta.cols(), beta.mats());

    // Copy values to test code.  Must be using known value with correct dim.
    // w.copy(w_mf);
    // beta.copy(beta_mf);

    // Run gibbs.
    logit.gibbs(w, beta, m0, P0, *samp, *burn, r);

    // Copy values to return.
    w_mf.copy(w);
    beta_mf.copy(beta);

    // Adjust for combined data.
    *N = w.rows();
  }
  catch (std::exception& e) {
    Rprintf("Error: %s\n", e.what());
    Rprintf("Aborting Gibbs sampler.\n");
  }

  #ifdef USE_R
  PutRNGstate();
  #endif
} // mult_gibbs

void mult_combine(double *typ, double *tXp, double *np,
		  int *N, int *P, int *J)
{
  int unk = *J-1; // number unknown

  // Set up data.
  Matrix ty(typ, unk, *N);
  Matrix tX(tXp,  *P, *N);
  Matrix n (np,   *N,  1);

  // Logit Gibbs
  try{
    MultLogit logit(ty, tX, n);

    logit.get_data(ty, tX, n);

    // Copy.
    MatrixFrame ty_mf   (typ, ty.rows(), ty.cols(), ty.mats());
    MatrixFrame tX_mf   (tXp, tX.rows(), tX.cols(), tX.mats());
    MatrixFrame n_mf    (np ,  n.rows(),  n.cols(),  n.mats());

    ty_mf.copy(ty);
    tX_mf.copy(tX);
    n_mf.copy(n);

    *N = tX.cols();

  }
  catch (std::exception& e) {
    Rprintf("Error: %s\n", e.what());
    Rprintf("Aborting combine.\n");
  }

} // mult_combine
