#include "InvertY2.h"

#ifdef USE_R
#include "R.h"
#endif

//------------------------------------------------------------------------------

YV::YV() : T(gsl_root_fdfsolver_newton), s(NULL)
{
  s = gsl_root_fdfsolver_alloc (T);

  FDF.f   = &f_eval;
  FDF.df  = &df_eval;
  FDF.fdf = &fdf_eval;
  FDF.params = 0;
  
  ylower = ygrid[0];
  yupper = ygrid[grid_size-1];
}

YV::~YV() { 
  if (s != NULL) { 
    gsl_root_fdfsolver_free(s); 
    s = NULL; 
  } 
}

double YV::y_func(double v) {
  return y_eval(v);
}

double YV::v_func(double y, int maxiter) 
{
  double v = 1.0;

  double ycopy = y;
  FDF.params = &ycopy;

  if (y < ylower) {
    return -1. / (y*y);
  } else if (y > yupper) {
    v = atan(0.5 * y * IYPI);
    return v*v;
  }
  else if (y==1) return 0.0;
    
  double id = (log(y) / log(2) + 4.0) / 0.1;
  // Rprintf("y, id, y[id], v[id]: %g, %g, %g, %g\n", y, id, ygrid[(int)id], vgrid[(int)id]);
  
  gsl_root_fdfsolver_set(s, &FDF, vgrid[(int)id]);

  int iter = 0;
  int status = 0;
  double vp = 0.0;
  // double fval, dfval;

  do {
    iter++;
    status = gsl_root_fdfsolver_iterate (s);
    vp = v;
    v  = gsl_root_fdfsolver_root(s);
    status = gsl_root_test_delta(v, vp, 0, 1e-8);

    // ydy_eval(v, &fval, &dfval);
    // Rprintf("yval, dyval, v: %g, %g, %g\n", fval, dfval, v);

    // fdf_eval(v, FDF.params, &fval, &dfval);
    // Rprintf("fval, dfval: %g, %g\n", fval, dfval);

  } while (status == GSL_CONTINUE && iter < maxiter);

  if (iter >= maxiter) {
    Rprintf( "YV: v reached maxiter.  ");
    Rprintf( "y: %g; init v: %g; cur. v: %g\n", y, vgrid[(int)id], v);
  }

  return v;
}

double YV::upperIncompleteGamma(double x, double shape, double rate)
{
  double t = rate * x;
  return gsl_sf_gamma_inc(shape, t);
}
