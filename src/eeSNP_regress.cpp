#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include "LogitWrapper.h"



using namespace arma;		// shorthand
typedef  std::vector<int> index_t;
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
    int ncols = sigma.n_cols;
    arma::mat Y = arma::randn(n, ncols);
    return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}

vec logistic(vec x)
{
    return(1/(1 +exp(-(x))));
}
double logistic(double x)
{
    return(1/(1 +exp(-(x))));
}

double logit(double x)
{
    return(log(x/(1-x)));
}
vec logit(vec x)
{
    return(log(x/(1-x)));
}


using namespace arma;		// shorthand
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
Rcpp::List eeSNP_regress(arma::mat x, arma::vec y, arma::mat feature, arma::colvec alpHa,  arma::colvec beTa, bool estimate_alpha, bool estimate_beta, arma::colvec B_inv_alpHa, int itermax, int  thin, int burnIn)
{
    //initial checking;
    bool oversample = false;
    int num_out_samples = ceil((itermax - burnIn)/thin);
    int p=x.n_cols,  n=x.n_rows; //, numSamp_logit = n/2, ;
    mat feature1 = join_rows( ones(p,1), feature);
    int f = feature1.n_cols;
    mat beTa_mat(p, num_out_samples), alpHa_mat(f, num_out_samples), gamMa_mat(p, num_out_samples);
    colvec w(p); 
    y = y - (ones(n) * mean(y)) ;

    int count  = 0; // double tau2 = .003; int c2 = 10000;
    colvec logit_freq = ones(p), logit_freq_nonzero = ones(p), gamMa_inclusion=zeros(p), zz=zeros(p);
    Rcpp::NumericVector sample_prob(p);
    mat diag_w, V_w; vec m_w(p);
    double sigma2 = var(y);
    if(beTa.size() !=p) beTa = zeros(p);
    //double v0=0, v1 =  (max(y) - min(y))/2*3, delta = .01*sqrt(sigma2)/2;
    double v0=0, v1 =  sigma2/2, delta = .1*sqrt(sigma2)/2;
    //#eqn 8 of Geroge et. al. (approaches for BVS)
    double minR = 10, maxR = 1e12, ratio=0, mindelc = 1e5, delc =0;
    double bb, maxbinsearch = 1e3;
    for (bb = 0; bb < maxbinsearch; bb++) {
      ratio = exp( log(minR) + ( bb/maxbinsearch) *log(maxR/minR)); 
      delc = abs(log( v1*log(ratio) /(ratio -1)) - 2*log(delta));
      if(mindelc> delc){
	mindelc = delc; v0 = v1/ratio;
      } 
    }

    int localmaxiter = 1;//TODO put it as parameter
    uvec nonzero; int k;
    mat omega(16,16); colvec omegaxy(p);
    vec residue=zeros(p);
    int bootstrap=1000; vec ret(bootstrap); uvec bootstrapS(p); 
    bootstrapS = linspace<uvec>(0, p-1, p);
    mat feature_nonzero; 
    vec fTalpHa(p);
    double rate, post0, post1, prior0, prior1;
    int p_nonzero;
    //return logit_freq
    vec gamMa = zeros(p);
    colvec gamMa_prob=logistic(feature1 * alpHa); 
    for (int snp = 0; snp < p; snp++){ 
	gamMa(snp)  =  (Rcpp::rbinom(1, 1, gamMa_prob(snp)))(0);
	//if (gamMa_prob(snp) >= .5) gamMa(snp) =1;
    }
    gamMa_inclusion = gamMa;
    //std::cout<<"gamMa" <<  gamMa << std::endl;
    for (int iter = 0; iter < itermax; iter++) {
	//std::cout << "iter" << iter<< std::endl;
	if(estimate_beta){
	    //#gibbs sampling;
	    if (oversample) logit_freq.zeros();
	    gamMa_inclusion.zeros();
	    nonzero.clear();
	    nonzero = find(gamMa > 0);
	    k = nonzero.size();
	    //for (int snp = 0; snp < p; snp++){ 
		//if(gamMa >0) nonzero.push_back(snp);
	    //}
	    for (int localiter = 0; localiter < localmaxiter; localiter++) {
		mat x_nonzero(k,p);
		x_nonzero = x.cols(nonzero);
		//#Chol = chol(t(x_nonzero) %*% x_nonzero + diag(1/sigma2A, k,k));
		//#xty =  t(x_nonzero) %*% y;
		//#Ctxty = solve( t(Chol), xty) ;
		//beTa = zeros(p);
		//std::cout <<"here95"<< std::endl;
		beTa = randn(p) * pow(v0, .5);
		//std::cout <<"here97"<< beTa<< std::endl;
		if(k >0){ 
    if(k > omega.n_cols) omega.set_size(k, k);
    //std::cout << omega.n_cols  << std::endl;
		    omega(span(0,k-1), span(0, k-1)) = inv(trans(x_nonzero) * x_nonzero + v1*eye(k,k)) ;
		    omegaxy(span(0,k-1)) = omega(span(0,k-1), span(0, k-1)) * trans(x_nonzero) * y;
		    //std::cout <<"sigma2:"<< sigma2<< std::endl;
		    ////std::cout << omega(span(1,k), span(1, k))<< std::endl;
		    ////std::cout << x_nonzero<< std::endl;
		    ////std::cout << omega(span(1,k), span(1, k))  << std::endl;
    //Rcpp::List out1; out1["alpHa"] = x_nonzero; out1["mult"]=trans(x_nonzero) * x_nonzero + c2*tau2*eye(k,k);
    //return out1;
    //std::cout << "nonzero"<<  trans(nonzero)<< std::endl;
    //std::cout << "omegaxy"<<  trans(omegaxy)<< std::endl;
    ////std::cout << "beTa"<< mvrnormArma(1, omegaxy(span(0,k-1)),  omega(span(0,k-1), span(0, k-1))*sigma2) << std::endl;
    			beTa(nonzero)  = mvrnormArma(1, omegaxy(span(0,k-1)),  omega(span(0,k-1), span(0, k-1))*sigma2);
    //beTa(nonzero)  = omegaxy(span(0,k-1)); 
    //std::cout << "beTa"<<  trans(beTa(nonzero)) << std::endl;
		    
		}
		residue = y - x_nonzero * beTa(nonzero) ;
		double temp1=R::rgamma( n/2 , 2/sum(square(residue)) );
		sigma2 = 1/temp1;
		for (int snp = 0; snp < p; snp++) {
		    prior1 = R::dnorm(beTa(snp), 0.0, pow(v1, .5), true);
		    post1 =  prior1 + log(gamMa_prob(snp));
		    prior0 = R::dnorm(beTa(snp), 0.0, pow(v0, .5), true);
		    post0 =  prior0 + log(1 - gamMa_prob(snp));
		    rate =  logistic( post1 - post0);
		    gamMa(snp)  =  (Rcpp::rbinom(1, 1, rate))(0);
		}
    //std::cout << "gamMa"<< trans(gamMa) << std::endl;
		//# sample for logistic regression;
		nonzero.clear();
		nonzero = find(gamMa > 0);
		 gamMa_inclusion = gamMa_inclusion + gamMa;
		k = nonzero.size();
	    }
	    if(oversample){
		double sum_inclusion = sum(gamMa_inclusion);
		if(localmaxiter*p  > sum_inclusion) 
		    sample_prob = (localmaxiter*p - sum_inclusion) * gamMa_inclusion +  sum_inclusion *(localmaxiter- gamMa_inclusion);
		else std::fill(sample_prob.begin(), sample_prob.end(), 1.0);

		//std::cout << (localmaxiter*p - sum(gamMa_inclusion)) * gamMa_inclusion +  sum(gamMa_inclusion) *(localmaxiter- gamMa_inclusion) << std::endl;
		//std::cout << sample_prob.size()<< std::endl;
		//std::cout << "bootstrapS:"<< bootstrapS.size() << std::endl;
		//std::cout << sample_prob<< std::endl;

		uvec ret = Rcpp::RcppArmadillo::sample(bootstrapS, bootstrap, true, sample_prob);
		//logit_freq = as_numeric(tabulate(sample_int(n=p, size=1000, replace=T, prob=sample_prob), p));
		for (int ii = 0; ii < bootstrap; ii++) {
		    logit_freq(ret(ii)) = 	logit_freq(ret(ii))+1; 
		}
	    }else logit_freq = ones(p);
	    gamMa_inclusion = gamMa_inclusion/localmaxiter;
	}
	if(estimate_alpha){
	    //#sample wi;
	    uvec freq_nonzero = find(logit_freq > 0); p_nonzero = freq_nonzero.size(); zz = logit_freq % (gamMa_inclusion - ones(p)/2);
	    feature_nonzero = feature1.rows(freq_nonzero); 
    //std::cout << logit_freq(freq_nonzero) << std::endl;
    //std::cout << logit_freq_nonzero.size()<< ":" << p_nonzero << std::endl;
	    logit_freq_nonzero(span(0,p_nonzero-1)) = logit_freq(freq_nonzero);
	    fTalpHa(span(0,p_nonzero-1)) = feature_nonzero * alpHa;
	    ////std::cout << w << std::endl;
    w = zeros(p);
	int temp =1;
    //rpg_hybrid( &w[0], &logit_freq_nonzero[0], &fTalpHa[0], &p_nonzero);
    for (int wc = 0; wc < p_nonzero; wc++) {
	rpg_hybrid( &w[wc], &logit_freq_nonzero[wc], &fTalpHa[wc], &temp);
    }	    
    //rpg_hybrid( &w[0], &logit_freq_nonzero[0], &fTalpHa[0], &p_nonzero);
    //std::cout << "w"<< trans(w) << std::endl;
    //std::cout << "fTalpHa"<< trans(fTalpHa) << std::endl;
    //std::cout << "logit_freq"<< trans(logit_freq) << std::endl;
    //std::cout << "zz"<< trans(zz) << std::endl;
	    //calculate m_w V_w;
	    mat diag_w  = diagmat(w(span(0, p_nonzero-1 )));
	    mat V_w = ((trans(feature_nonzero) * diag_w * feature_nonzero)+  diagmat( B_inv_alpHa)).i();
	    vec m_w = V_w * (trans(feature_nonzero) * zz(freq_nonzero));
	    //std::cout<< "V_w" << V_w << "m_w" << m_w << std::endl;
	    //std::cout<< "feature_nonzero" << feature_nonzero << "fTz" << trans(feature_nonzero) * zz(freq_nonzero) << std::endl;
	    ////std::cout<< mvrnormArma(1, m_w, V_w) << std::endl;
	    alpHa = trans(mvrnormArma(1, m_w, V_w));
	    gamMa_prob = logistic(feature1 * alpHa);
	}

	if((((iter + 1) %thin)==0)  && (iter >= burnIn)){
	    gamMa_mat.col(count) = gamMa;
	    beTa_mat.col(count) = beTa;
	    alpHa_mat.col(count) = alpHa;
	    
	    count = count +1;;
	}

    }
    Rcpp::List out; out["alpHa.mat"] = trans(alpHa_mat); out["beTa.mat"] = trans(beTa_mat); out["gamMa"] = trans(gamMa_mat); 
    
    return out;
}
