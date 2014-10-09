#include <Eigen/Dense>
#include <stdexcept>
#include "RNG.h"
#include "CUBS_update.h"
#include "CUBS.h"

using Eigen::Map;

// We could use #define to simplify this, but that makes it harder to bebug.

////////////////////////////////////////////////////////////////////////////////
				  // EXTERN //
////////////////////////////////////////////////////////////////////////////////

void cubs_norm(double *alpha_, double *beta_,
	       double *z_, double *X_, double *V_,
	       double *mu_, double *phi_, double *W_, 
	       double *m0_, double *C0_, int *N_b_, int *N_, int *T_,
	       double *log_dens, double *eps_rel_, int* max_iter_)
{
  RNG r;
  
  int T = *T_;
  int N = *N_;
  int N_b = *N_b_;
  int N_a = N - N_b;

  if (T > 10000 || N > 1000) {
    Rprintf( "cubs: T=%i or N=%i is very large.  Aborting.\n", T, N);
    // throw std::runtime_error("T or N is very large.  Aborting.\n");
    return;
  }

  Map<VectorXd> alpha(alpha_, N_a > 1 ? N_a : 1);
  Map<MatrixXd> beta (beta_ , N_b, T+1);
  Map<VectorXd> z    (z_    , T       );
  Map<MatrixXd> X    (X_    , T  , N  );
  Map<VectorXd> V    (V_    , T       );
  Map<VectorXd> mu   (mu_   , N_b     );
  Map<VectorXd> phi  (phi_  , N_b     );
  Map<MatrixXd> W    (W_    , N_b, N_b);
  Map<VectorXd> m0   (m0_   , N       );
  Map<MatrixXd> C0   (C0_   , N  , N  );

  // NormUpdate norm;
  // cubs(alpha, beta, z, X, V, mu, phi, W, m0, C0, log_dens, norm, *eps_rel_, *max_iter_, r);
  cubs(alpha, beta, z, X, V, mu, phi, W, m0, C0, log_dens, &norm_post, *eps_rel_, *max_iter_, r);
}

void cubs_binom(double *alpha_, double *beta_,
		double *z_, double *X_, double *n_,
		double *mu_, double *phi_, double *W_, 
		double *m0_, double *C0_, int *N_b_, int *N_, int *T_,
		double *log_dens, double *eps_rel_, int* max_iter_)
{
  RNG r;
  
  int T = *T_;
  int N = *N_;
  int N_b = *N_b_;
  int N_a = N - N_b;
  
  if (T > 10000 || N > 1000) {
    Rprintf( "cubs: T=%i or N=%i is very large.  Aborting.\n", T, N);
    // throw std::runtime_error("T or N is very large.  Aborting.\n");
    return;
  }

  Map<VectorXd> alpha(alpha_, N_a > 1 ? N_a : 1);
  Map<MatrixXd> beta (beta_ , N_b, T+1);
  Map<VectorXd> z    (z_    , T       );
  Map<MatrixXd> X    (X_    , T  , N  );
  Map<VectorXd> n    (n_    , T       );
  Map<VectorXd> mu   (mu_   , N_b     );
  Map<VectorXd> phi  (phi_  , N_b     );
  Map<MatrixXd> W    (W_    , N_b, N_b);
  Map<VectorXd> m0   (m0_   , N       );
  Map<MatrixXd> C0   (C0_   , N  , N  );

  // BinomUpdate binom;
  // cubs(alpha, beta, z, X, n, mu, phi, W, m0, C0, log_dens, binom, *eps_rel_, *max_iter_, r);
  cubs(alpha, beta, z, X, n, mu, phi, W, m0, C0, log_dens, &binom_post, *eps_rel_, *max_iter_, r);
}

void cubs_nbinom(double *alpha_, double *beta_,
		 double *z_, double *X_, double *n_,
		 double *mu_, double *phi_, double *W_, 
		 double *m0_, double *C0_, int *N_b_, int *N_, int *T_,
		 double *log_dens, double* eps_rel_, int* max_iter_)
{
  RNG r;
  
  int T = *T_;
  int N = *N_;
  int N_b = *N_b_;
  int N_a = N - N_b;
  
  if (T > 10000 || N > 1000) {
    Rprintf( "cubs: T=%i or N=%i is very large.  Aborting.\n", T, N);
    // throw std::runtime_error("T or N is very large.  Aborting.\n");
    return;
  }

  Map<VectorXd> alpha(alpha_, N_a > 1 ? N_a : 1);
  Map<MatrixXd> beta (beta_ , N_b, T+1);
  Map<VectorXd> z    (z_    , T       );
  Map<MatrixXd> X    (X_    , T  , N  );
  Map<VectorXd> n    (n_    , T       );
  Map<VectorXd> mu   (mu_   , N_b     );
  Map<VectorXd> phi  (phi_  , N_b     );
  Map<MatrixXd> W    (W_    , N_b, N_b);
  Map<VectorXd> m0   (m0_   , N       );
  Map<MatrixXd> C0   (C0_   , N  , N  );

  try {
    // NBinomUpdate nbinom;
    // cubs(alpha, beta, z, X, n, mu, phi, W, m0, C0, log_dens, nbinom, *eps_rel_, *max_iter_, r);
    cubs(alpha, beta, z, X, n, mu, phi, W, m0, C0, log_dens, &nbinom_post, *eps_rel_, *max_iter_, r);
  }
  catch (std::exception& e) {
    Rprintf( "%s", e.what());
  }
}
