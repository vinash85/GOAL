#ifndef __INVERTY2__
#define __INVERTY2__

#include "InvertY.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf.h>

//------------------------------------------------------------------------------

class YV {

private:
  double yupper;
  double ylower;

public:

  const gsl_root_fdfsolver_type * T;
  gsl_root_fdfsolver            * s;
  gsl_function_fdf              FDF;
  
  YV();
  virtual ~YV();

  double y_func(double v);
  double v_func(double y, int maxiter=100);

  double upperIncompleteGamma(double x, double shape, double rate);

};

//------------------------------------------------------------------------------

////////////////////////////////////////////////////////////////////////////////

#endif
