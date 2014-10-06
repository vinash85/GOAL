// -*- mode: c++; -*-
#ifndef __CUBS_UPDATE__
#define __CUBS_UPDATE__

#include <gsl/gsl_multiroots.h>

////////////////////////////////////////////////////////////////////////////////

typedef void (*CUBS_transform) (const double* rs, const double* fq, double* out);
typedef int  (*CUBS_gsl_call)  (const gsl_vector* x, void* param, gsl_vector* f);
typedef void (*CUBS_update)    (const double* prior, double* post, double y, double n, const double* ival, double epsrel, int max_iter); 

////////////////////////////////////////////////////////////////////////////////

void binom_transform (const double* rs, const double* fq, double* out);
void utest_binom_transform(double r, double s);

void binom_approx_transform(const double* rs, const double* fq, double* out);
void utest_binom_approx_transform(double r, double s);

int  binom_transform_gsl (const gsl_vector* x, void* p, gsl_vector* f);
void utest_binom_transform_gsl(double r, double s);

int solver(const double* fq, double* rs, double epsabs, double epsrel, int max_iter,
	   int (*gsl_transform) (const gsl_vector*, void*, gsl_vector*));

// "CUBS_gsl_call gls_transofrm" is the same as:
// int (*gsl_transform) (const gsl_vector*, void*, gsl_vector*)

void binom_post (const double* prior, double* post, double y, double n, const double* ival, double epsrel, int max_iter);
void nbinom_post(const double* prior, double* post, double y, double n, const double* ival, double epsrel, int max_iter);
void norm_post  (const double* prior, double* post, double y, double n, const double* ival, double epsrel, int max_iter);

void binom_post_paranoid (const double* prior, double* post, double y, double n, const double* ival, double epsrel, int max_iter);
void nbinom_post_paranoid(const double* prior, double* post, double y, double n, const double* ival, double epsrel, int max_iter);

////////////////////////////////////////////////////////////////////////////////
			       // CLASS BASED //
////////////////////////////////////////////////////////////////////////////////

//------------------------------------------------------------------------------
class CUBSSolver
{
protected:
  const gsl_multiroot_fsolver_type * T;
  gsl_multiroot_fsolver            * s;
  CUBS_gsl_call gsl_transform;

public:

  CUBSSolver();
  CUBSSolver(CUBS_gsl_call gsl_transform_);
  virtual ~CUBSSolver();

  void solve(const double* fq, double* rs, double epsabs, double epsrel, int max_iter);

};

//------------------------------------------------------------------------------
class CUBSUpdate 
{
public:
  virtual void update(const double* prior, double* post, double y, double n, double epsrel, int max_iter) = 0; 
};

// DERIVED CLASSES
class NormUpdate : public CUBSUpdate
{
  void update(const double* prior, double* post, double y, double n, double epsrel, int max_iter);
};

class BinomUpdate : public CUBSUpdate
{
  CUBSSolver cs;
public:
  BinomUpdate() : cs(&binom_transform_gsl) {};
  void update(const double* prior, double* post, double y, double n, double epsrel, int max_iter);
};

class NBinomUpdate : public CUBSUpdate
{
  CUBSSolver cs;
public:
  NBinomUpdate() : cs(&binom_transform_gsl) {};
  void update(const double* prior, double* post, double y, double n, double epsrel, int max_iter);
};

#endif
