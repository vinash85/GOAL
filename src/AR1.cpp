#include "AR1.h"

using Eigen::Map;

void ar1_llh(double *alpha_, double *beta_,
	     double *mu_, double *phi_, double *W_, 
	     double *m0_, double *C0_, int *N_b_, int *N_, int *T_,
	     double *llh)
{
  int T   = *T_;
  int N   = *N_;
  int N_b = *N_b_;
  int N_a =  N - N_b;

  Map<VectorXd> alpha(alpha_, N_a > 1 ? N_a : 1);
  Map<MatrixXd> beta (beta_ , N_b, T+1);
  Map<VectorXd> mu   (mu_   , N_b     );
  Map<VectorXd> phi  (phi_  , N_b     );
  Map<MatrixXd> W    (W_    , N_b, N_b);
  Map<VectorXd> m0   (m0_   , N       );
  Map<MatrixXd> C0   (C0_   , N  , N  );

  *llh = ar1_llh(alpha, beta, mu, phi, W, m0, C0);

}
