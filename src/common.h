#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include "LogitWrapper.h"
#define sqr(x) x*x
using namespace arma;		// shorthand
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma); 
vec logistic(vec x);
double logistic(double x);
vec logit(vec x);
double logit(double x);
