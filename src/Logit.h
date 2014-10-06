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

// This class possess the functions needed to implement Polson and Scott's
// "Default Bayesian Logistic Regression."  See <http://arxiv.org/pdf/1109.4180>
// for details for the algorithm.

#ifndef __LOGIT__
#define __LOGIT__

#include "Matrix.h"
#include "RNG.h"
#include "PolyaGamma.h"
#include "Normal.h"
#include <stdexcept>
#include <list>
#include <time.h>
#include <algorithm>
#include <stdio.h>

// I need this to interrupt the Gibbs sampler from R.
#ifdef USE_R
#include <R_ext/Utils.h>
#endif

using std::list;

////////////////////////////////////////////////////////////////////////////////
				  // LOGIT //
////////////////////////////////////////////////////////////////////////////////

class Logit{

  // Variables.
  uint P;
  uint N;

  // Data
  Matrix tX; // Transpose of design matrix.
  Matrix n;
  Matrix y; // Not a sufficient statistic, but need when combining data.

  // Prior precision, mean, b = Pm.
  Matrix P0;
  Matrix m0;
  Matrix b0;

  // Posterior precision, b = Pm.
  Matrix PP;
  Matrix bP;

  // Random variates.
  Normal mvnorm;
  PolyaGamma pg;

public:

  // Constructors.
  Logit();
  // tX_data is transpose of design matrix!
  Logit(const Matrix& y_data, const Matrix& tX_data, const Matrix& n_data);

  // Utilities.
  void set_data (const Matrix& y_data, const Matrix& tX_data, const Matrix& n_data);
  void set_prior(const Matrix& m0_, const Matrix& P0_);

  void compress();

  bool data_conforms(const Matrix& y_data, const Matrix& tX_data, const Matrix& n_data);

  void get_data(Matrix& y_data, Matrix& tX_data, Matrix& n_data);

  uint get_N() { return N; }
  uint get_P() { return P; }

  // Gibbs sampling -- default
  inline void draw_w   (MF w, MF psi, RNG& r);
  inline void draw_beta(MF beta, MF w, MF beta_prev, RNG& r);
  inline void draw_beta(MF beta, MF w, RNG& r);
  void gibbs(Matrix& w, Matrix& beta, int samp, int burn, RNG& r);

  double gibbs_block(MF beta_space, MF w_space,
		     MF beta_init, MF w_init,
		     int samp, int period, RNG& r);

  // Exepectation Maximization.
  int EM(Matrix& beta, double tol=1e-9, int max_iter=1000);

protected:

  void set_bP();

}; // Logit

////////////////////////////////////////////////////////////////////////////////
			       // CONSTRUCTORS //
////////////////////////////////////////////////////////////////////////////////

Logit::Logit()
{
  // We do not want to call this constructor.
  // Rprintf("You must add data.");
} // Logit

Logit::Logit(const Matrix& y_data , const Matrix& tX_data, const Matrix& n_data)
  : mvnorm(tX_data.rows())
{
  set_data(y_data, tX_data, n_data);
} // Logit

////////////////////////////////////////////////////////////////////////////////
				// UTILITIES //
////////////////////////////////////////////////////////////////////////////////

bool Logit::data_conforms(const Matrix& y_data, const Matrix& tX_data, const Matrix& n_data)
{
  bool ok = true;
  bool check[2];

  ok *= check[0] = y_data.area()  == tX_data.cols();
  ok *= check[1] = y_data.area()  == n_data.area();

  for(int i = 0; i < 2; i++)
    if (!check[i]) Rprintf("Problem with check %i .\n", i);

  return ok;
}

void Logit::set_data(const Matrix& y_data, const Matrix& tX_data, const Matrix& n_data)
{
  // Check that the data is valid.
  if (!data_conforms(y_data, tX_data, n_data)) {
    throw std::invalid_argument("set_data: data does not conform.");
  }

  P = tX_data.rows();
  N = tX_data.cols();

  // Set data.
  y.clone(y_data);
  tX.clone(tX_data);
  n.clone(n_data);

  // Set default prior.
  P0.resize(P, P); P0.fill(0.0);
  m0.resize(P)   ; m0.fill(0.0);
  b0.resize(P)   ; b0.fill(0.0);

  // Set posterior dim.
  PP.resize(P, P);
  bP.resize(P, 1);

}

void Logit::set_bP()
{
  // Set up bP and alpha.
  Matrix alpha(N);
  bP.clone(b0);

  for(uint i = 0; i < N; ++i)
    alpha(i) = n(i) * (y(i) - 0.5);
  gemm(bP, tX, alpha, 'N', 'N', 1.0, 1.0);
}

void Logit::set_prior(const Matrix& m0_, const Matrix& P0_)
{
  m0.clone(m0_);
  P0.clone(P0_);
  mult(b0, P0, m0);
}

void Logit::compress()
{
  // Push everything into a list.
  list<double> ylist;
  list<Matrix> xlist;
  list<double> nlist;

  // Our data should not have n_data(i)=0.
  for(uint i=0; i < N; ++i){
      ylist.push_back(y(i));
      xlist.push_back(tX.col(i));
      nlist.push_back(n(i));
  }

  // Merge data.
  list<double>::iterator y_i;
  list<Matrix>::iterator x_i;
  list<double>::iterator n_i;

  list<double>::iterator y_j;
  list<Matrix>::iterator x_j;
  list<double>::iterator n_j;

  x_i = xlist.begin();
  y_i = ylist.begin();
  n_i = nlist.begin();

  while(x_i != xlist.end()){
    x_j = x_i; y_j = y_i; n_j = n_i;
    ++x_j; ++y_j; ++n_j;
    while(x_j != xlist.end()){
      if (*x_i == *x_j) {
  	double sum = *n_i + *n_j;
   	*y_i = (*n_i/sum) * *y_i + (*n_j/sum) * *y_j;
  	*n_i = sum;
  	// Increment THEN erase!
	// Actually, send pointer to erase, then increment, then erase.
  	xlist.erase(x_j++);
  	ylist.erase(y_j++);
  	nlist.erase(n_j++);
      }
      else {
  	++x_j; ++y_j; ++n_j;
      }

    }
    ++x_i; ++y_i; ++n_i;
  }

  uint old_N = N; // Record to see if we have changed data.

  // Set up y, tX, n.
  N = xlist.size();

  // cout << "Old N: " << old_N << " N: " << N << "\n";
  // Warning...
  if (old_N != N) {
    Rprintf("Warning: data was combined!\n");
    Rprintf("N: %i, P: %i \n", N, P);
  }

  // Matrix y(N);
  y.resize(N);
  tX.resize(P, N);
  n.resize(N);

  x_i = xlist.begin();
  y_i = ylist.begin();
  n_i = nlist.begin();
  for(uint i = 0; i < N; ++i){
    y(i)      = *y_i++;
    tX.col(i) = *x_i++;
    n(i)      = *n_i++;
  }

  // cout << "y:\n" << y;
  // cout << "tX:\n" << tX;
  // cout << "n:\n" << n;
}

void Logit::get_data(Matrix& y_data, Matrix& tX_data, Matrix& n_data)
{
  y_data  = y;
  tX_data = tX;
  n_data  = n;
}

////////////////////////////////////////////////////////////////////////////////
			    // POSTERIOR BY GIBBS //
////////////////////////////////////////////////////////////////////////////////

inline void Logit::draw_w(MF w, MF psi, RNG& r)
{
  for(int i = 0; i < (int)N; ++i) {
    // w(i) = pg.draw(n(i), psi(i), r);
    w(i) = pg.draw((int)n(i), psi(i), r);
  }
}

inline void Logit::draw_beta(MF beta, MF w, RNG& r)
{
  // tXRtOm = tX sqrt(Om)
  Matrix tXRtOm(P, N);
  for(unsigned int j=0; j<tX.cols(); j++)
    for(unsigned int i=0; i<tX.rows(); i++)
      tXRtOm(i,j) = tX(i,j) * sqrt(w(j));

  // PP = X' Om X + P0.
  PP = P0; syrk(PP, tXRtOm, 'N', 1.0, 1.0);

  // PP = U'U.
  // U'U mP = bP --> U' y = bP; U mP = y;
  // ep = U z
  // dr = mP + z.

  Matrix U; chol(U, PP, 'U');

  r.norm(beta, 1.0);
  Matrix mP(bP);

  trsm(U, mP,   'U', 'L', 'T'); // U' y = bP
  trsm(U, mP,   'U', 'L', 'N'); // U mP = y
  trsm(U, beta, 'U', 'L', 'N'); // U z = ep.

  for(uint i=0; i<P; i++)
    beta(i) += mP(i);
}

inline void Logit::draw_beta(MF beta, MF w, MF beta_prev, RNG& r)
{
  // Maybe have a special case for when P0 = 0.
  Matrix tXOmega(tX);
  prodonrow(tXOmega, w);
  Matrix tXOmX; mult(tXOmX, tXOmega, tX, 'N', 'T');
  PP = P0; PP += tXOmX;
  // Joint draw.
  mvnorm.set_from_likelihood(bP, PP);
  mvnorm.draw(beta[0], r);

  // // Gibbs sampling.
  // beta.copy(beta_prev);
  // for(uint i = 0; i < P; ++i){
  //   double m_i = bP(i);
  //   for(uint j = 0; j < P; ++j){
  // 	m_i -= j != i ? tXOmX(i,j) * beta(j) : 0.0;
  //   }
  //   m_i /= tXOmX(i,i);
  //   beta(i) = r.norm(m_i, sqrt(1/tXOmX(i,i)));
  // }
}

double Logit::gibbs_block(MF beta_space, MF w_space,
			  MF beta_init, MF w_init,
			  int samp, int period, RNG& r)
{
  // Basically, implementing a simple iterator here.

  double *beta_curr = &beta_space(0);
  double *beta_prev = beta_curr;
  double *w_curr = &w_space(0);

  MF beta_curr_mf(beta_curr, (int)P, 1, 1);
  MF beta_prev_mf(beta_prev, (int)P, 1, 1);
  MF w_curr_mf(w_curr, (int)N, 1, 1);
  Matrix psi((int)N, 1, 1);

  // Initialize.
  beta_curr_mf.copy(beta_init);
  beta_prev_mf.copy(beta_curr_mf);
  w_curr_mf.copy(w_init);
  gemm(psi, tX, beta_curr_mf, 'T');

  clock_t start, end;
  start = clock();

  for (int m=1; m<=samp*period; m++) {

    draw_w   (w_curr_mf,  psi, r);
    // draw_beta(beta_curr_mf, w_curr_mf, beta_prev_mf, r);
    draw_beta(beta_curr_mf, w_curr_mf, r);
    gemm(psi, tX, beta_curr_mf, 'T');

    // Increment on i-th iteration % period.
    if (m % period == 0) {
      beta_prev = beta_curr;
      beta_curr += P;
      w_curr += N;
      // Could also do something like:
      // beta_curr += P * sizeof(double) * (double)(m % thin);

      beta_curr_mf.setp(beta_curr);
      beta_prev_mf.setp(beta_prev);
      w_curr_mf.setp(w_curr);
    } // Increment.

    #ifdef USE_R
    if (m%100==0) R_CheckUserInterrupt();
    #endif

  }

  end = clock();

  double total_time = (double)(end - start) / CLOCKS_PER_SEC;

  return total_time;
} // gibbs_block

// Gibbs sampling -- Default Logit.
void Logit::gibbs(Matrix& w, Matrix& beta, int samp, int burn, RNG& r)
{
  set_bP();

  uint M = (uint)samp;

  w.resize(N, 1, M);
  beta.resize(P, 1, M);
  // beta.fill(-1.0);

  double total_time;

  // Burn.
  total_time = gibbs_block(beta, w, beta[0], w[0], 1, burn, r);
  Rprintf("Burn-in complete: %g sec. for %i iterations.\n", total_time, burn);
  Rprintf("Expect approx. %g sec. for %i samples.\n", total_time * samp / burn, samp);

  // Sample.
  total_time = gibbs_block(beta, w, beta[0], w[0], samp, 1, r);
  Rprintf("Sampling complete: %g sec. for %i iterations.\n", total_time, samp);

} // gibbs

////////////////////////////////////////////////////////////////////////////////
			   // POSTERIOR MODE BY EM //
////////////////////////////////////////////////////////////////////////////////

// Exepectation Maximization.
int Logit::EM(Matrix& beta, double tol, int max_iter)
{
  // Preprocess.
  set_bP();

  Matrix U; 
  Matrix psi(N);
  Matrix w(N);
  double dist = tol + 1.0;

  // Proper size.
  beta.resize(P);
  beta.fill(0.0);

  int  iter = 0;
  while (dist > tol && iter < max_iter) {

    // std::cout << beta;

    // w: posterior mean
    gemm(psi, tX, beta, 'T', 'N');
    for (int i = 0; i < (int)N; ++i) {
      double hpsi = psi(i) * 0.5;
      // n * tanh(psi/2) / (psi/2) * 0.5
      if ( fabs(hpsi) < 0.01 ) {
	w(i) = n(i) / cosh(hpsi)
	  * (1 + hpsi*hpsi / 6.0 + pow(hpsi, 4.0) / 120.0 + pow(hpsi, 6) / 5040.0)
	  * 0.25 ;
      }
      else
	w(i) = n(i) * tanh(hpsi) / hpsi * 0.25;
    }

    // beta: posterior mode
    Matrix old_beta(beta);

    // tXRtOm = tX sqrt(Om)
    Matrix tXRtOm(P, N);
    for(unsigned int j=0; j<tX.cols(); j++) {
      double rtw = sqrt(w(j));
      for(unsigned int i=0; i<tX.rows(); i++) {
	tXRtOm(i,j) = tX(i,j) * rtw;
      }
    }

    // PP = X' Om X + P0.
    PP = P0;
    // syrk(PP, tXRtOm, 'N', 1.0, 1.0);
    gemm(PP, tXRtOm, tXRtOm, 'N', 'T', 1.0, 1.0);

    chol(U, PP, 'U');
    beta.clone(bP);
    trsm(U, beta,   'U', 'L', 'T'); // U' y = bP
    trsm(U, beta,   'U', 'L', 'N'); // U beta = y
    // symsolve(PP, beta);

    // Check how much we improved.
    // Matrix diff = beta - old_beta;
    // dist = sqrt( dot(diff, diff) );
    Matrix diff = fabs(beta - old_beta);
    dist = *std::max_element(diff.begin(), diff.end());

    ++iter;
  }

  return iter;
}

////////////////////////////////////////////////////////////////////////////////
			       // END OF CODE //
////////////////////////////////////////////////////////////////////////////////

#endif

////////////////////////////////////////////////////////////////////////////////
				 // APPENDIX //
////////////////////////////////////////////////////////////////////////////////

/*

  The biggest bottleneck in our algorithm appears to be drawing from the Poyla
  Gamma distribution.  There may not be anyway around this.  I tried sampling
  all the gamma variables at once, but this didn't seem to improve things.  In
  fact, it seemed to make things take longer.  Basically, it seems that you the
  overhead to create an array to to sample a bunch of gammas at "one time" does
  not overcome any optential speedup from preenting repeatedly calling the
  function.  Of course it may be that the compiler is optimizing things and
  inlining things.

  Basically, it appears to just be a problem to sample a Poyla Gamma.

  UPDATE: You can sample a PG more quickly following the method of Devroye.

  ---

  The algorithm specifies drawing a block of omegas and then a block of psi's.
  But due to the independence of each we can sample them in an alternatating
  fashion.  But this doesn't work.  Why?  Too many degrees of freedom when using
  psi?  We could parallelize this... if this worked.

 */
