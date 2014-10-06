/*
  This class represents a normal random variable.  It allows you to
  draw from the distribution of a normal random variable with the mean
  and covariance chosen by the user.
*/

// We only want to read this file once.
// Begin if statement for compiler.
#ifndef _NORMAL_
#define _NORMAL_

// Include the following.
#include <stdexcept>
#include <iostream>

// We want to use dgrid and rng.
#include "Matrix.h"
#include "RNG.h"

class Normal
{
 public:
  //========== CONSTRUCTORS ==========//

  //! Construct an n-dimensional standard Normal RV.
  Normal(uint n=1);

  //! Construct a N(m, V) random variate.
  Normal(MF& m, MF& V);

  //========== METHODS ==========/

  // Sample from this random variable.
  void draw(MF d, RNG & r);

  // Sample using a pointer.  Error prone.
  void draw(double *d, RNG& r);

  // Update the mean and variance.
  void set(MF& m, MF& V);

  //
  void set_from_likelihood(MF& b, MF& Prec);

  //========== DATA ==========/

  //! The dimension of the normal random variable.
  uint dim;

  //! The mean of the normal random variable.
  Matrix mean;

  //! The Cholesky decomposition of the normal random variable.
  Matrix lower;

  //! Identity;
  Matrix Id;

};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

//Construct an n-dimensional standard normal RV.
Normal::Normal(uint n)
  : dim(n)
  , mean(dim)
  , lower("I", dim)
  , Id("I", dim) {}

Normal::Normal(MF & m, MF & V)
  : dim(m.rows())
  , mean(dim)
  , lower(dim, dim)
  , Id("I", dim)
{
  set(m, V);
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void Normal::set(MF& m, MF& V)
{
  // Check dimensions.
  if (m.rows() != dim || V.rows() != dim || V.cols() != dim) {
    throw std::invalid_argument("set: data does not conform.");
  }

  // Mean.
  mean.copy(m);

  // Cholesky.
  lower.copy(Id);
  chol(lower, V, 'L');
}

void Normal::set_from_likelihood(MF& b, MF& R)
{
  // Check dimensions.
  if (b.rows() != dim || R.rows() != dim || R.cols() != dim) {
    throw std::invalid_argument("set: data does not conform.");
  }

  // Variance
  Matrix V("I", dim);
  symsolve(R, V);

  // Mean
  gemm(mean, V, b);

  // Cholesky
  lower.copy(Id);
  chol(lower, V, 'L');
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void Normal::draw(MF d, RNG& r)
{
  // Draw uniform normal.
  r.norm(d, 0.0, 1.0);

  // draw = Chol * draw.
  trmm(lower, d, 'L');
  // trmm beats out gemm.  See Appendix B.

  // d = d + mean;
  hsumeq(d, mean);
}

void Normal::draw(double *d, RNG& r)
{
  MatrixFrame mf(d, dim);
  draw(mf, r);
}



#endif

