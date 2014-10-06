#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_sf.h>
#include <stdio.h>
#include <cstring>
#include "CUBS_update.h"

#ifdef USE_R
#include "R.h"
#endif

void binom_transform (const double* rs, const double* fq, double* out)
{
  double r=rs[0]; double s=rs[1];
  double f=fq[0]; double q=fq[1];
  double E = gsl_sf_psi(r) - gsl_sf_psi(s);
  double V = gsl_sf_psi_1(r) + gsl_sf_psi_1(s);
  out[0] = E - f;
  out[1] = V - q;
  // Rprintf("r=%g, s=%g, f=%g, q=%g, out=%g, %g\n", r, s, f, q, out[0], out[1]);
}

void utest_binom_transform(double r, double s)
{
  double rs[2]; rs[0] = r; rs[1] = s;
  double fq[2] = {0, 0};
  double out[2];

  binom_transform(rs, fq, out);
  Rprintf("r=%g, s=%g, f=%g, q=%g\n", rs[0], rs[1], out[0], out[1]);
}

int binom_transform_gsl (const gsl_vector* x, void* p, gsl_vector* f) {
  double rs[2];
  double out[2];
  double* fq = (double *)p;

  rs[0] = gsl_vector_get(x, 0);
  rs[1] = gsl_vector_get(x, 1);

  binom_transform(rs, fq, out);

  gsl_vector_set (f, 0, out[0]);
  gsl_vector_set (f, 1, out[1]);

  return GSL_SUCCESS;
}

void utest_binom_transform_gsl(double r, double s)
{
  gsl_vector* x = gsl_vector_alloc(2);
  x->data[0] = r;
  x->data[1] = s;

  double fq[2] = {0, 0};
  void* fqv = (void *)fq;

  gsl_vector* out = gsl_vector_alloc(2);
  binom_transform_gsl(x, fqv, out);

  Rprintf("r=%g, s=%g, f=%g, q=%g\n", r, s, out->data[0], out->data[1]);
}

int binom_transform_df(const gsl_vector* x, void* p, gsl_matrix* J)
{
  double r = gsl_vector_get(x, 0);
  double s = gsl_vector_get(x, 1);

  double df00 = gsl_sf_psi_1(r);
  double df01 = -1.0 * gsl_sf_psi_1(s);
  double df10 = gsl_sf_psi_n(2, r);
  double df11 = gsl_sf_psi_n(2, s);

  gsl_matrix_set(J, 0, 0, df00);
  gsl_matrix_set(J, 0, 1, df01);
  gsl_matrix_set(J, 1, 0, df10);
  gsl_matrix_set(J, 1, 1, df11);

  return GSL_SUCCESS;
}

int binom_transform_fdf(const gsl_vector * x, void *params, gsl_vector * f, gsl_matrix * J)
{
  binom_transform_gsl(x, params, f);
  binom_transform_df(x, params, J);

  return GSL_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
			  // FUNCTION BASED SOLVER //
////////////////////////////////////////////////////////////////////////////////

int solver(const double* fq, double* rs, const double* ival, double epsabs, double epsrel, int max_iter,
	   int (*gsl_transform) (const gsl_vector*, void*, gsl_vector*))
{
  #ifdef USE_R
  gsl_set_error_handler_off ();
  #endif

  double params[2]; memmove(params, fq, 2 * sizeof(double));
  // fq[0] = prior[0]; fq[1] = prior[1];

  const gsl_multiroot_fsolver_type * T = gsl_multiroot_fsolver_hybrid;
  gsl_multiroot_fsolver            * s = gsl_multiroot_fsolver_alloc(T, 2);

  gsl_multiroot_function F;

  // Set up F.
  F.f = gsl_transform;
  F.n = 2;
  F.params = (void *)params;

  // Set up initial vector.
  gsl_vector* x = gsl_vector_alloc(2);
  memcpy(x->data, ival, 2 * sizeof(double));

  gsl_multiroot_fsolver_set(s, &F, x);
  // Rprintf("x: %g, %g \t f: %g, %g\n", s->x->data[0], s->x->data[1], s->f->data[0], s->f->data[0]);

  int i = 0;
  int msg = GSL_CONTINUE;
  for(i = 0; i < max_iter && msg != GSL_SUCCESS; i++) {
    msg = gsl_multiroot_fsolver_iterate(s);
    if (msg == GSL_EBADFUNC || msg == GSL_ENOPROG) break;
    // Rprintf("x: %g, %g \t f: %g, %g \t dx: %g, %g\n", s->x->data[0], s->x->data[1],
    // s->f->data[0], s->f->data[0], s->dx->data[0], s->dx->data[1]);
    // check |dx| < epsabs + epsrel * |x|
    msg = gsl_multiroot_test_delta(s->dx, s->x, epsabs, epsrel);
  }

  // You can turn off GSL error handling so it doesn't crash things.
  if (msg != GSL_SUCCESS) {
    Rprintf( "CUBS_udpate.cpp::solver Error %i.  Break on %i.\n", msg, i);
    Rprintf( "error: %s\n", gsl_strerror (msg));
    Rprintf( "Init: r=%g, s=%g, f=%g, q=%g\n", ival[0], ival[1], fq[0], fq[1]);
    Rprintf( "Exit: r=%g, s=%g, ", s->x->data[0], s->x->data[1]);
    Rprintf( "F0=%g, F1=%g, ", s->f->data[0], s->f->data[1]);
    Rprintf( "D0=%g, D1=%g\n", s->dx->data[0], s->dx->data[1]);
  }

  memmove(rs, s->x->data, 2 * sizeof(double));

  // Free mem.
  gsl_multiroot_fsolver_free (s);
  gsl_vector_free (x);

  return msg;
}

////////////////////////////////////////////////////////////////////////////////
		       // BINOM SOLVER USING GRADIENT //
////////////////////////////////////////////////////////////////////////////////

int binom_solver(const double* fq, double* rs, const double* ival, double epsabs, double epsrel, int max_iter)
{
  #ifdef USE_R
  gsl_set_error_handler_off ();
  #endif

  double params[2]; memmove(params, fq, 2 * sizeof(double));
  // fq[0] = prior[0]; fq[1] = prior[1];

  const gsl_multiroot_fdfsolver_type *T;
  gsl_multiroot_fdfsolver *s;

  const size_t n = 2;
  // Set up F.
  gsl_multiroot_function_fdf F = {&binom_transform_gsl,
				  &binom_transform_df,
				  &binom_transform_fdf,
				  n, (void *)params};

  // Set up initial vector.
  gsl_vector* x = gsl_vector_alloc(2);
  memcpy(x->data, ival, 2 * sizeof(double));

  T = gsl_multiroot_fdfsolver_gnewton;
  s = gsl_multiroot_fdfsolver_alloc (T, n);
  gsl_multiroot_fdfsolver_set (s, &F, x);

  // Rprintf("x: %g, %g \t f: %g, %g\n", s->x->data[0], s->x->data[1], s->f->data[0], s->f->data[0]);

  int i = 0;
  int msg = GSL_CONTINUE;
  for(i = 0; i < max_iter && msg != GSL_SUCCESS; i++) {
    msg = gsl_multiroot_fdfsolver_iterate(s);
    if (msg == GSL_EBADFUNC || msg == GSL_ENOPROG) break;
    // Rprintf("x: %g, %g \t f: %g, %g \t dx: %g, %g\n", s->x->data[0], s->x->data[1],
    // s->f->data[0], s->f->data[0], s->dx->data[0], s->dx->data[1]);
    // check |dx| < epsabs + epsrel * |x|
    msg = gsl_multiroot_test_delta(s->dx, s->x, epsabs, epsrel);
  }

  // You can turn off GSL error handling so it doesn't crash things.
  if (msg != GSL_SUCCESS) {
    Rprintf( "CUBS_udpate.cpp::solver Error %i.  Break on %i.\n", msg, i);
    Rprintf( "error: %s\n", gsl_strerror (msg));
    Rprintf( "Init: r=%g, s=%g, f=%g, q=%g\n", ival[0], ival[1], fq[0], fq[1]);
    Rprintf( "Exit: r=%g, s=%g, ", s->x->data[0], s->x->data[1]);
    Rprintf( "F0=%g, F1=%g, ", s->f->data[0], s->f->data[1]);
    Rprintf( "D0=%g, D1=%g\n", s->dx->data[0], s->dx->data[1]);
  }

  memmove(rs, s->x->data, 2 * sizeof(double));

  // Free mem.
  gsl_multiroot_fdfsolver_free (s);
  gsl_vector_free (x);

  return msg;
}

////////////////////////////////////////////////////////////////////////////////
			    // POSTERIOR ROUTINES //
////////////////////////////////////////////////////////////////////////////////

void binom_post(const double* prior, double* post, double y, double n, const double* ival, double epsrel, int max_iter)
{
  double rs[2];
  double zero[2] = {0, 0};
  // solver(prior, rs, ival, 0.0, epsrel, max_iter, &binom_transform_gsl);
  binom_solver(prior, rs, ival, 0.0, epsrel, max_iter);
  rs[0] = rs[0] + y;
  rs[1] = rs[1] + n - y;
  binom_transform(rs, zero, post);
  // Rprintf("fq: %g, %g; rs: %g, %g\n", post[0], post[1], rs[0], rs[1]);
  // Rprintf("epsrel: %g, max_iter, %i\n", epsrel, max_iter);
}

void nbinom_post(const double* prior, double* post, double y, double n, const double* ival, double epsrel, int max_iter)
{
  double rs[2];
  double zero[2] = {0, 0};
  double fq[2]; memmove(fq, prior, 2 * sizeof(double));
  fq[0] = fq[0] - log(n);
  // solver(fq, rs, ival, 0.0, epsrel, max_iter, &binom_transform_gsl);
  binom_solver(fq, rs, ival, 0.0, epsrel, max_iter);
  rs[0] = rs[0] + y;
  rs[1] = rs[1] + n;
  binom_transform(rs, zero, post);
  post[0] = post[0] + log(n);
}

void norm_post  (const double* prior, double* post, double y, double n, const double* ival, double epsrel, int max_iter)
{
  double f = prior[0];
  double q = prior[1];
  post[1] =  (q * n) / (n + q);
  post[0] = ((f / q) + (y / n)) * post[1];
}

////////////////////////////////////////////////////////////////////////////////
			       // CLASS BASED //
////////////////////////////////////////////////////////////////////////////////

CUBSSolver::CUBSSolver(CUBS_gsl_call gsl_transform_)
  : T(gsl_multiroot_fsolver_hybrid)
  , s(NULL)
  , gsl_transform(gsl_transform_)
{
  s = gsl_multiroot_fsolver_alloc(T, 2);
}

CUBSSolver::~CUBSSolver()
{
  gsl_multiroot_fsolver_free (s);
}

void CUBSSolver::solve(const double* fq, double* rs, double epsabs, double epsrel, int max_iter)
{
  #ifdef USE_R
  gsl_set_error_handler_off ();
  #endif

  double params[2]; memmove(params, fq, 2 * sizeof(double));
  // fq[0] = prior[0]; fq[1] = prior[1];

  gsl_multiroot_function F;

  // Set up F.
  F.f = gsl_transform;
  F.n = 2;
  F.params = (void *)params;

  // Set up initial vector.
  gsl_vector* x = gsl_vector_alloc(2);

  gsl_vector_set_all(x, 0.01);
  gsl_multiroot_fsolver_set(s, &F, x);
  // Rprintf("x: %g, %g \t f: %g, %g\n", s->x->data[0], s->x->data[1], s->f->data[0], s->f->data[0]);

  int i = 0;
  int msg = GSL_CONTINUE;
  for(i = 0; i < max_iter && msg != GSL_SUCCESS; i++) {
    msg = gsl_multiroot_fsolver_iterate(s);
    if (msg == GSL_EBADFUNC || msg == GSL_ENOPROG) break;
    // Rprintf("x: %g, %g \t f: %g, %g\n", s->x->data[0], s->x->data[1], s->f->data[0], s->f->data[0]);
    // check |dx| < epsabs + epsrel * |x|
    msg = gsl_multiroot_test_delta(s->dx, s->x, epsabs, epsrel);
  }

  // You can turn off GSL error handling so it doesn't crash things.
  if (msg != GSL_SUCCESS) {
    Rprintf( "CUBSSolver::solve: Error %i.  Break on %i.\n", msg, i);
    Rprintf( "error: %s\n", gsl_strerror (msg));
    Rprintf( "r=%g, s=%g\n", s->x->data[0], s->x->data[1]);
    Rprintf( "f=%g, q=%g\n", s->f->data[0], s->f->data[1]);
  }

  memmove(rs, s->x->data, 2 * sizeof(double));

  // Free mem.
  gsl_vector_free (x);
}

//------------------------------------------------------------------------------

// Update
void BinomUpdate::update(const double* prior, double* post, double y, double n, double epsrel, int max_iter)
{
  double rs[2];
  double zero[2] = {0, 0};
  cs.solve(prior, rs, 0.0, epsrel, max_iter);
  rs[0] = rs[0] + y;
  rs[1] = rs[1] + n - y;
  binom_transform(rs, zero, post);
  // Rprintf("fq: %g, %g; rs: %g, %g\n", post[0], post[1], rs[0], rs[1]);
  // Rprintf("epsrel: %g, max_iter, %i\n", epsrel, max_iter);
}

void NBinomUpdate::update(const double* prior, double* post, double y, double n, double epsrel, int max_iter)
{
  double rs[2];
  double zero[2] = {0, 0};
  double fq[2]; memmove(fq, prior, 2 * sizeof(double));
  fq[0] = fq[0] - log(n);
  cs.solve(fq, rs, 0.0, epsrel, max_iter);
  rs[0] = rs[0] + y;
  rs[1] = rs[1] + n;
  binom_transform(rs, zero, post);
  post[0] = post[0] + log(n);
}

void NormUpdate::update  (const double* prior, double* post, double y, double n, double epsrel, int max_iter)
{
  double f = prior[0];
  double q = prior[1];
  post[1] =  (q * n) / (n + q);
  post[0] = ((f / q) + (y / n)) * post[1];
}

////////////////////////////////////////////////////////////////////////////////

void test_time_fast(unsigned int N)
{
  double prior[2] = {0.0, 3.3};
  double rs[2];
  CUBSSolver cs(&binom_transform_gsl);
  for(unsigned int i = 0; i < N; i++)
    cs.solve(prior, rs, 0.0, 1e-8, 100);
}

void test_time_slow(unsigned int N)
{
  double prior[2] = {0.0, 3.3};
  double rs[2];
  for(unsigned int i = 0; i < N; i++) {
    CUBSSolver cs(&binom_transform_gsl);
    cs.solve(prior, rs, 0.0, 1e-8, 100);
  }
}

#ifdef CUBS_UPDATE_MAIN
int main(int argc, char** argv)
{
  // Function tests
  utest_binom_transform_gsl(1, 1);
  utest_binom_transform(1, 1);

  double prior[2] = {0.0, 3.3};
  double post[2];
  double rs[2];
  double ival[2] = {0.1, 0.1};

  // Exact
  Rprintf("Exact...\n");
  solver(prior, rs, ival, 0.0, 1e-8, 100, &binom_transform_gsl);
  Rprintf("rs = %g, %g\n", rs[0], rs[1]);

  binom_post(prior, post, 1, 1, ival, 1e-8, 100);
  Rprintf("post = %g, %g\n", post[0], post[1]);

  // Newton
  Rprintf("Exact with Jacobian...\n");
  binom_solver(prior, rs, ival, 0.0, 1e-8, 100);
  Rprintf("rs = %g, %g\n", rs[0], rs[1]);

  binom_post(prior, post, 1, 1, ival, 1e-8, 100);
  Rprintf("post = %g, %g\n", post[0], post[1]);

  // Class tests
  Rprintf("Exact (by class)...\n");
  CUBSSolver cs(&binom_transform_gsl);
  cs.solve(prior, rs, 0.0, 1e-8, 100);
  Rprintf("rs = %g, %g\n", rs[0], rs[1]);

  BinomUpdate binom;
  binom.update(prior, post, 1, 1, 1e-8, 100);
  Rprintf("post = %g, %g\n", post[0], post[1]);

  return 0;
}
#endif

////////////////////////////////////////////////////////////////////////////////
				 // APPENDIX //
////////////////////////////////////////////////////////////////////////////////

// I was running to into errors saying the "solver was not making progress" or
// something like that.  If you think about the form of the function, when r=s,
// then you do not change the f coordinate.  Thus, if you move in the (1,1)
// direction, the f coordinate will not change.  Thus, if you guess to move in
// that direction it will seem like nothing is happening, no matter the step
// size.  Not sure if this was precisely what was happening, but something like
// that may have been occuring.  Things seems to work better after havign
// specified the Jacobian.

////////////////////////////////////////////////////////////////////////////////

// #define GSLTRANSFORM(NAME, CALL)			
//   int NAME (const gsl_vector* x, void* p, gsl_vector* f) {	
//     double rs[2];						
//     double out[2];						
//     double* fq = (double *)p;				
//     								
//     rs[0] = gsl_vector_get(x, 0);				
//     rs[1] = gsl_vector_get(x, 1);				
//     								
//     CALL (rs, fq, out);					
//     								
//     gsl_vector_set (f, 0, out[0]);				
//     gsl_vector_set (f, 1, out[1]);				
//     								
//     return GSL_SUCCESS;					
//   }								

// GSLTRANSFORM(binom_transform_gsl, binom_transform)


////////////////////////////////////////////////////////////////////////////////
				   // OLD //
////////////////////////////////////////////////////////////////////////////////

// Prior to adding the Jacobian I thought it might be useful to approximate
// digamma and trigamma functions, find the root, then use that root as an
// initial value.  Just including Jacobian appears to be a better idea.

// int binom_approx_transform_gsl (const gsl_vector* x, void* p, gsl_vector* f) {
//   double rs[2];
//   double out[2];
//   double* fq = (double *)p;

//   rs[0] = gsl_vector_get(x, 0);
//   rs[1] = gsl_vector_get(x, 1);

//   binom_approx_transform(rs, fq, out);

//   gsl_vector_set (f, 0, out[0]);
//   gsl_vector_set (f, 1, out[1]);

//   return GSL_SUCCESS;
// }

// void binom_approx_transform(const double* rs, const double* fq, double* out)
// {
//   double r=rs[0]; double s=rs[1];
//   double f=fq[0]; double q=fq[1];
//   double E = log(r) - log(s);
//   double V = 1.0 / r + 1.0 / s;
//   out[0] = E - f;
//   out[1] = V - q;
// }

// void utest_binom_approx_transform(double r, double s)
// {
//   double rs[2]; rs[0] = r; rs[1] = s;
//   double fq[2] = {0, 0};
//   double out[2];

//   binom_approx_transform(rs, fq, out);
//   Rprintf("r=%g, s=%g, f=%g, q=%g\n", rs[0], rs[1], out[0], out[1]);
// }

// void nbinom_post_paranoid(const double* prior, double* post, double y, double n, const double* ival, double epsrel, int max_iter)
// {
//   double rs[2];
//   double zero[2] = {0, 0};
//   double fq[2]; memmove(fq, prior, 2 * sizeof(double));
//   fq[0] = fq[0] - log(n);
//   int msg;
//   msg = solver(fq, rs, ival, 0.0, epsrel, max_iter, &binom_approx_transform_gsl);
//   if (msg != 0) Rprintf( "approx: solver problem %i\n", msg);
//   double seed[2]; seed[0] = rs[0]; seed[1]=rs[1];
//   msg = solver(fq, rs, seed, 0.0, epsrel, max_iter, &binom_transform_gsl);
//   if (msg != 0) Rprintf( "nbinom: solver problem %i\n", msg);
//   rs[0] = rs[0] + y;
//   rs[1] = rs[1] + n;
//   binom_transform(rs, zero, post);
//   post[0] = post[0] + log(n);
// }

// void binom_approx_post(const double* prior, double* post, double y, double n, const double* ival, double epsrel, int max_iter)
// {
//   double rs[2];
//   double zero[2] = {0, 0};
//   solver(prior, rs, ival, 0.0, epsrel, max_iter, &binom_approx_transform_gsl);
//   rs[0] = rs[0] + y;
//   rs[1] = rs[1] + n - y;
//   binom_approx_transform(rs, zero, post);
//   // Rprintf("fq: %g, %g; rs: %g, %g\n", post[0], post[1], rs[0], rs[1]);
//   // Rprintf("epsrel: %g, max_iter, %i\n", epsrel, max_iter);
// }

// void binom_post_paranoid(const double* prior, double* post, double y, double n, const double* ival, double epsrel, int max_iter)
// {
//   double rs[2];
//   double zero[2] = {0, 0};
//   solver(prior, rs, ival, 0.0, epsrel, max_iter, &binom_approx_transform_gsl);
//   double seed[2]; seed[0] = rs[0]; seed[1]=rs[1];
//   solver(prior, rs, seed, 0.0, epsrel, max_iter, &binom_transform_gsl);
//   rs[0] = rs[0] + y;
//   rs[1] = rs[1] + n - y;
//   binom_transform(rs, zero, post);
// }
