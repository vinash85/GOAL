// -*- mode: c++; -*-
#ifndef __AR1__
#define __AR1__

#include <vector>
#include <Eigen/Dense>
#include <stdio.h>

#ifdef USE_R
#include <R_ext/Utils.h>
#endif

using std::vector;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::MatrixBase;
using Eigen::Matrix;
using Eigen::Dynamic;
using Eigen::LLT;


////////////////////////////////////////////////////////////////////////////////

extern "C" {

  void ar1_llh(double *alpha_, double *beta_,
	       double *mu_, double *phi_, double *W_, 
	       double *m0_, double *C0_, int *N_b_, int *N_, int *T_,
	       double *log_dens);

}

////////////////////////////////////////////////////////////////////////////////

// dM and dV let us switch between Matrix and Map.
template <typename dM, typename dV>
double ar1_llh(MatrixBase<dV> &alpha, MatrixBase<dM> &beta,
	       MatrixBase<dV> &mu, MatrixBase<dV> &phi, MatrixBase<dM> &W, 
	       MatrixBase<dV> &m0, MatrixBase<dM> &C0)
{
  int T   = beta.cols() - 1;
  int N   = m0.size();
  int N_b = beta.rows();
  int N_a = N - N_b;

  double llh = 0.0;

  // t = 1..T
  // LLT<MatrixXd> llt = W.llt();
  MatrixXd      L   = W.llt().matrixL();
  double        CL  = L.diagonal().array().log().sum();

  VectorXd m_i(N_b);
  VectorXd e_i(N_b);

  for (int i=1; i < (T+1); i++) {
    m_i = mu.array() + phi.array() * (beta.col(i-1).array() - mu.array());
    e_i = beta.col(i) - m_i;
    // llh += -0.5 * llt.solve(e_i).squaredNorm() - CL;
    llh += -0.5 * L.triangularView<Eigen::Lower>().solve(e_i).squaredNorm() - CL;
  }

  // t = 0
  L   = C0.llt().matrixL();;
  CL  = L.diagonal().array().log().sum();

  VectorXd theta(N); theta.segment(N_a, N_b) = beta.col(0);
  if (N_a > 0)       theta.segment(0  , N_a) = alpha;

  e_i = theta - m0;
  llh += -0.5 * L.triangularView<Eigen::Lower>().solve(e_i).squaredNorm() - CL;

  // return
  return llh;
}

#endif
