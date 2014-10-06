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

#ifndef __MULTLOGIT__
#define __MULTLOGIT__

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

class MultLogit{

  // Variables.
  uint P; // Dimension of beta_j.
  uint N; // Number of observations.
  uint J; // Total number of categories

  // Sufficient Statistics
  Matrix tX;
  Matrix n;
  Matrix Z;
  Matrix ty; // Not a sufficient statistic, but need when combining data.

  // Random variates.
  Normal mvnorm;
  PolyaGamma pg;

public:

  // Constructors.
  MultLogit();
  MultLogit(const Matrix& ty_data, const Matrix& tX_data, const Matrix& n_data);

  // Utilities.
  void set_data(const Matrix& ty_data, const Matrix& tX_data, const Matrix& n_data);

  bool data_conforms(const Matrix& ty_data, const Matrix& tX_data, const Matrix& n_data);

  void get_data(Matrix& ty_data, Matrix& tX_data, Matrix& n_data);

  uint get_N() { return N; }
  uint get_P() { return P; }

  // Gibbs samplign -- Normal prior.
  inline void draw_w(MF w, MF psi, RNG& r);
  inline void draw_beta(MF beta, MF w, MF Z_j, MF C_j, MF b0, MF P0, MF beta_prev, RNG& r);
  void gibbs(Matrix& w, Matrix& beta, Matrix& m0, Matrices &P0, int samp, int burn, RNG& r);

  // Exepectation Maximization.
  // int EM(Matrix& beta, double tol, int max_iter);

}; // MultLogit

////////////////////////////////////////////////////////////////////////////////
			       // CONSTRUCTORS //
////////////////////////////////////////////////////////////////////////////////

MultLogit::MultLogit()
{
  // We do not want to call this constructor.
  throw std::invalid_argument("MultLogit: default constructor called.");
} // MultLogit

MultLogit::MultLogit(const Matrix& ty_data, const Matrix& tX_data, const Matrix& n_data)
  : mvnorm(tX_data.rows())
{
  set_data(ty_data, tX_data, n_data);
} // MultLogit

////////////////////////////////////////////////////////////////////////////////
				// UTILITIES //
////////////////////////////////////////////////////////////////////////////////

bool MultLogit::data_conforms(const Matrix& ty_data, const Matrix& tX_data,
			      const Matrix& n_data)
{
  bool ok = true;
  bool check[2];

  ok *= check[0] = ty_data.cols() == tX_data.cols();
  ok *= check[1] = ty_data.cols() == n_data.rows();

  for(int i = 0; i < 2; i++)
    if (!check[i]) Rprintf("Problem with check %i .\n", i);

  return ok;
}

void MultLogit::set_data(const Matrix& ty_data, const Matrix& tX_data, const Matrix& n_data)
{

  // Check that the data is valid.
  if (!data_conforms(ty_data, tX_data, n_data)) {
    throw std::invalid_argument("set_data: data does not conform.");
  }

  Matrix tX_temp(tX_data);
  Matrix ty_temp(ty_data);

  P = tX_data.rows();
  N = tX_data.cols();
  J = ty_data.rows() + 1;

  // Push everything into a list.
  list<Matrix> ylist;
  list<Matrix> xlist;
  list<double> nlist;

  // Our data should not have n_data(i)=0.
  for(uint i=0; i < N; ++i){
      ylist.push_back(ty_temp.col(i));
      xlist.push_back(tX_temp.col(i));
      nlist.push_back(n_data(i));
  }

  // Merge data.
  list<Matrix>::iterator y_i;
  list<Matrix>::iterator x_i;
  list<double>::iterator n_i;

  list<Matrix>::iterator y_j;
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

  // Set up ty, tX, n.
  N = xlist.size();

  // cout << "Old N: " << old_N << " N: " << N << "\n";
  // Warning...
  if (old_N != N) {
    Rprintf("Warning: data was combined!\n");
    Rprintf("N: %i, P: %i \n", N, P);
  }

  // Matrix y(N);
  ty.resize(J-1, N);
  tX.resize(P, N);
  n.resize(N);

  x_i = xlist.begin();
  y_i = ylist.begin();
  n_i = nlist.begin();
  for(uint i = 0; i < N; ++i){
    ty.col(i) = *y_i++;
    tX.col(i) = *x_i++;
    n(i)      = *n_i++;
  }

  // cout << "y:\n" << y;
  // cout << "tX:\n" << tX;
  // cout << "n:\n" << n;

  // Set up Z and kappa.
  Matrix tkappa(J-1, N);
  Z.resize(P, J-1);
  for(uint i = 0; i < N; ++i)
    tkappa.col(i) = n(i) * (ty.col(i) - 0.5);
  gemm(Z, tX, tkappa, 'N', 'T');

}

void MultLogit::get_data(Matrix& ty_data, Matrix& tX_data, Matrix& n_data)
{
  ty_data = ty;
  tX_data = tX;
  n_data  = n;
}

////////////////////////////////////////////////////////////////////////////////
			    // POSTERIOR BY GIBBS //
////////////////////////////////////////////////////////////////////////////////

inline void MultLogit::draw_w(MF w, MF psi, RNG& r)
{
  for(int i = 0; i < (int)N; ++i) {
    // w(i) = pg.draw(n(i), psi(i), r);
    w(i) = pg.draw((int)n(i), psi(i), r);
  }
}

inline void MultLogit::draw_beta(MF beta, 
				 MF w, MF Z_j, MF C_j, MF b0, MF P0, 
				 MF beta_prev, RNG& r)
{
    Matrix tXOmega(tX);
    prodonrow(tXOmega, w);
    Matrix tXOmX; mult(tXOmX, tXOmega, tX, 'N', 'T');
    Matrix tXOmC; mult(tXOmC, tXOmega, C_j);

    // Joint draw.
    Matrix P1 = tXOmX + P0;
    Matrix b1 = Z_j + tXOmC + b0;

    mvnorm.set_from_likelihood(b1, P1);
    mvnorm.draw(beta[0], r);

}

// Gibbs sampling -- Default MultLogit.
void MultLogit::gibbs(Matrix& w, Matrix& beta, Matrix& m0, Matrices &P0, 
		      int samp, int burn, RNG& r)
{
  uint M = (uint)samp;

  w.resize(N, J-1, M);
  beta.resize(P, J-1, M);
  // beta.fill(-1.0);

  // Set up prior.
  Matrix b0(P, J-1);
  for (uint j=0; j < J-1; j++)
    gemm(b0.col(j), P0[j], m0[j]);

  // Could put in the loop to reset every iteration.
  Matrix XB(N, J);
  gemm(XB.col(0,J-1), tX, beta[0], 'T'); // Mult into first J-1 columns.
  Matrix XB_no_j(N, J-1);
  Matrix A(N);

  // Keep track of time.
  clock_t start, end;

  start = clock();
  // Burn-in - Take an extra for first sample of MCMC.
  for(int m = 0; m < burn+1; ++m){

    XB_no_j.copy(XB, seq((uint)0,N-1), seq((uint)1,J-1));

    // Gibbs for beta_j, w_j.
    for (uint j = 0; j < J-1; j++) {
      
      A = rowSums(exp(XB_no_j));
      
      // Worried abotu instability.
      // A -= exp(XB.col(j));

      Matrix c_j   = log(A);
      Matrix eta_j = XB.col(j) - c_j;

      draw_w(w[0].col(j), eta_j, r);

      draw_beta(beta[0].col(j), 
		w[0].col(j), Z.col(j), c_j, b0.col(j), P0[j], 
		beta[0].col(j), r);

      gemm(XB.col(j), tX, beta[0].col(j), 'T');
      
      // Worried about instabiltiy.
      // A += exp(XB.col(j));

      // If using the XB_no_j method.
      if (j < J-2) XB_no_j.col(j) = XB.col(j);

    }

    // In case we are using R.
    #ifdef USE_R
    if (m%1==0) R_CheckUserInterrupt();
    #endif

  }
  end = clock();

  double total_time = (double)(end - start) / CLOCKS_PER_SEC;
  Rprintf("Burn-in complete: %g sec. for %i iterations.\n", total_time, burn);
  Rprintf("Expect approx. %g sec. for %i samples.\n", total_time * samp / burn, samp);

  start = clock();
  // Sample - Already took one sample from burn-in.
  for(int m = 1; m < samp; ++m){

    XB_no_j.copy(XB, seq((uint)0,N-1), seq((uint)1,J-1));

    // Gibbs for beta_j, w_j.
    for (uint j = 0; j < J-1; j++) {
      
      A = rowSums(exp(XB_no_j));
      
      // Worried abotu instability.
      // A -= exp(XB.col(j));

      Matrix c_j   = log(A);
      Matrix eta_j = XB.col(j) - c_j;

      draw_w(w[m].col(j), eta_j, r);

      draw_beta(beta[m].col(j), 
		w[m].col(j), Z.col(j), c_j, b0.col(j), P0[j], 
		beta[m-1].col(j), r);

      gemm(XB.col(j), tX, beta[m].col(j), 'T');
      
      // Worried about instabiltiy.
      // A += exp(XB.col(j));

      // If using the XB_no_j method.
      if (j < J-2) XB_no_j.col(j) = XB.col(j);

    }

    // In case we are using R.
    #ifdef USE_R
    if (m%1==0) R_CheckUserInterrupt();
    #endif
  }
  end = clock();

  total_time = (double)(end - start) / CLOCKS_PER_SEC;
  Rprintf("Sampling complete: %g sec. for %i iterations.\n", total_time, samp);

}

////////////////////////////////////////////////////////////////////////////////
			   // POSTERIOR MODE BY EM //
////////////////////////////////////////////////////////////////////////////////

// // Exepectation Maximization.
// int MultLogit::EM(Matrix& beta, double tol, int max_iter)
// {
//   Matrix psi(N);
//   Matrix w(N);
//   double dist = tol + 1.0;

//   // Proper size.
//   beta.resize(P);

//   int  iter = 0;
//   while (dist > tol && iter < max_iter) {

//     // w: posterior mean
//     gemm(psi, tX, beta, 'T', 'N');
//     for (int i = 0; i < (int)N; ++i) {
//       double hpsi = psi(i) * 0.5;
//       // n * tanh(psi/2) / (psi/2) * 0.5
//       if ( hpsi < 0.01 ) {
// 	w(i) = n(i) / cosh(hpsi)
// 	  * (1 + hpsi*hpsi / 6.0 + pow(hpsi, 4.0) / 120.0 + pow(hpsi, 6) / 5040.0)
// 	  * 0.25 ;
//       }
//       else
// 	w(i) = n(i) * tanh(hpsi) / hpsi * 0.25;
//     }

//     // beta: posterior mode
//     Matrix old_beta(beta);

//     Matrix tXOmega(tX);
//     prodonrow(tXOmega, w);
//     Matrix tXOmX(tXOmega, tX, 'N', 'T');
//     beta.clone(Z); symsolve(tXOmX, beta);

//     // Check how much we improved.
//     // Matrix diff = beta - old_beta;
//     // dist = sqrt( dot(diff, diff) );
//     Matrix diff = fabs(beta - old_beta);
//     dist = *std::max_element(diff.begin(), diff.end());

//     ++iter;
//   }

//   return iter;
// }

////////////////////////////////////////////////////////////////////////////////
			       // END OF CODE //
////////////////////////////////////////////////////////////////////////////////

#endif

////////////////////////////////////////////////////////////////////////////////
				 // APPENDIX //
////////////////////////////////////////////////////////////////////////////////

