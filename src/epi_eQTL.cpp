#ifdef _OPENMP
#include <omp.h>
#endif
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include "LogitWrapper.h"
#include "common.h"
#include "ransampl.h"
#include "dbg.h"
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <time.h>
//#include <trng/yarn2.hpp>
//#include <trng/gamma_dist.hpp>
#include <math.h>

#define min(a,b) a<b?a:b 
#define max(a,b) a>b?a:b 
#define ABS(a) (a>0.0) ? a:-a
#define bounds(a, b ,c) min(max(a,b),c)
double square(double x)
{
    return(x*x);
}
void center_norm(mat & xx)
{
    int n = xx.n_rows; rowvec meanx = mean(xx), sdx = stddev(xx);
    xx = (xx - ones(n)*meanx) /(ones(n) * sdx); 
}

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
int bernoulli(double xx)
{
    vec temp = randu(1);
    int out = (xx > temp(0))? 1:0;
    return(out);
}
void center(mat & xx)
{
    int n = xx.n_rows; rowvec meanx = mean(xx);
    xx = (xx - ones(n)*meanx) ; 
}
using namespace arma;		// shorthand
typedef  std::vector<int> index_t;
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec find_correlation(Rcpp::List mask, arma::mat x, arma::mat y) {
    int yLen =  mask.size(), totalMaskLen = 0;
    int yy;
    for (yy = 0; yy < yLen; yy++){
	SEXP ll = mask[yy];
	Rcpp::NumericVector maskCurr(ll); 
	totalMaskLen += maskCurr.size();
    }
    arma::vec corr=zeros(totalMaskLen);
    int maskStart =0; 
    for (yy = 0; yy < yLen; yy++) {
	//cout << 38 << "\t" <<  yy <<endl;
	SEXP ll = mask[yy];
	Rcpp::NumericVector maskCurr(ll); 
	//cout << 41 << "\t" <<  yy <<endl;
	arma::rowvec yCurr=y.row(yy);
	//cout << 43 << "\t" <<  yy <<endl;
	for (int i = 0; i < maskCurr.size() ; i++) {
	    arma::mat temp = arma::cor(x.row(maskCurr[i]-1), yCurr);
	    //cout << 47 << "\t" <<  i <<endl;
	    corr(maskStart + i) = temp(0,0); 
	}
	maskStart += maskCurr.size();
    }
    return corr; 
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
Rcpp::List filterSNP(Rcpp::List obj, bool parseFeature, double threshold){

    //#ifdef _OPENMP
    //if ( threads > 0 )
    //omp_set_num_threads( threads );
    //REprintf("Number of threads=%i\n", omp_get_max_threads());
    //#endif
    mat x = obj["x"], y = obj["y"];
    imat mask2 = obj["mask2"];
    //mat feature, pairFeature; 
    //if(parseFeature) {
	mat feature = obj["feature"];
	mat pairFeature = obj["pairFeature.mat"];
    //}
    int n = x.n_rows;
    //Environment base("package:stats");
    Function qt("qt");
    Rcpp::NumericVector tt = qt(1.0 - threshold, n-2);
    double r_thres = tt[0]/ abs(sqrt(n -2 + tt[0] *tt[0])); 
    center_norm(y); center_norm(x);
    double numGenes=y.n_cols; 
    field<uvec> mask(numGenes);
    uvec maskCol = conv_to<umat>::from(mask2.col(1));
    field<uvec> inx(numGenes); uvec inxall;
    ivec outmaskvec2;
    //#pragma omp parallel for schedule(dynamic)   
    for (int gene = 0; gene < numGenes; gene++) {
	uvec mask = maskCol(find(mask2.col(0) == (gene +1) )) - 1; 	
	mat r = abs( y.col(gene).t() * x.cols(mask))/(n-1);
	//if(gene == 0) cout << r << endl;
	inx(gene) = mask(find(r.row(0) > r_thres));
	//p_curr = mask(gene).size();
    }
    //creating new mask
    for (int gene = 0; gene < numGenes; gene++) {
	inxall = join_cols(inxall, inx(gene)+1);
	double inxSize= inx(gene).size();
	ivec temp(inxSize); temp.fill(gene + 1);
	if(gene ==0 ) outmaskvec2 = temp;
	else outmaskvec2 = join_cols(outmaskvec2 , temp);
    }
    Rcpp::List out;
    int outn =  outmaskvec2.size();
    int maskColnum = mask2.n_cols;
    imat outmask2( outn, maskColnum);  
    outmask2.col(0) = outmaskvec2;
    outmask2.col(1) = linspace<ivec>(1, outn, outn);
    ivec mt;
    if(maskColnum >2){
	    for (int mm = 2; mm < maskColnum; mm++) {
		    mt = mask2.col( mm);
		    outmask2.col(mm) = mt(inxall -1); 
	    }
    }
    out["mask2"] = outmask2;
    mat outx	= x.cols(inxall-1); 
    out["x"] = outx;
    mat outF, outpairF;
    if(parseFeature){
	outF = feature.rows(inxall-1);  out["feature"] = outF;
	outpairF = pairFeature.rows(inxall-1); out["pairFeature.mat"] = outpairF;
    }
    out["inx"] = inx; out["inxall"]= inxall;
    return out;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec filter_inconsequential(Rcpp::List mask, arma::mat x, arma::mat y, int  threshold) {
    int yLen =  mask.size(), totalMaskLen = 0;
    int yy;
    for (yy = 0; yy < yLen; yy++){
	SEXP ll = mask[yy];
	Rcpp::NumericVector maskCurr(ll); 
	totalMaskLen += maskCurr.size();
    }
    arma::vec corr=zeros(totalMaskLen);
    int maskStart =0; 
    for (yy = 0; yy < yLen; yy++) {
	SEXP ll = mask[yy];
	Rcpp::NumericVector maskCurr(ll); 
	arma::rowvec yCurr=y.row(yy);
	for (int i = 0; i < maskCurr.size() ; i++) {
	    arma::rowvec xCurr=x.row(maskCurr[i]-1);
	    arma::mat temp = arma::cor(x.row(maskCurr[i]-1), yCurr);
	    corr(maskStart + i) = temp(0,0); 
	}
	maskStart += maskCurr.size();
    }
    return corr; 
}
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec univariate_eQTL(arma::umat mask2, arma::mat x, arma::mat y, arma::vec numRegulators, int threads =1) {
    #ifdef _OPENMP
    if ( threads > 0 )
    omp_set_num_threads( threads );
    REprintf("Number of threads=%i\n", omp_get_max_threads());
    #endif
    int yLen =  y.n_cols; 
    //int yLen =1;
    int yy;
    vec ev=zeros(yLen);
    uvec maskCol = mask2.col(1);
    #pragma omp parallel for schedule(dynamic)   
    for (int yy = 0; yy < yLen; yy++) {
	int regulators =  numRegulators(yy);
	if(regulators > 0){
	    uvec mask = maskCol(find(mask2.col(0) == (yy +1) )) - 1; 	
	    arma::mat xCurr = x.cols(mask);
	    center(xCurr);
	    arma::colvec yCurr= y.col(yy);
	    yCurr = yCurr - mean(yCurr);
	    arma::mat corr = arma::cor(xCurr, yCurr); // the output are in column
	    arma::uvec nz = sort_index(corr.col(0), "descend");
	    nz = nz(span(0, regulators -1));
	    xCurr = xCurr.cols(nz);
	    rowvec small(regulators);
	    small.fill(.001); 
	    mat V_w = ((trans(xCurr)  * xCurr)+  diagmat(small)).i();
	    vec m_w = V_w * (trans(xCurr) * yCurr);
	    vec residue = yCurr - ( xCurr * m_w);
	    double varresidue = var(residue);
	    double vary = var(yCurr);
	    ev(yy) = 1 - varresidue/vary;
	}
    }
    return ev; 
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec multivariate_eQTL(arma::umat mask2, arma::mat x, arma::mat y, arma::vec numRegulators, int threads =1) {
    #ifdef _OPENMP
    if ( threads > 0 )
    omp_set_num_threads( threads );
    REprintf("Number of threads=%i\n", omp_get_max_threads());
    #endif
    int yLen =  y.n_cols; 
    //int yLen =1;
    int yy;
    vec ev=zeros(yLen);
    uvec maskCol = mask2.col(1);
    #pragma omp parallel for schedule(dynamic)   
    for (int yy = 0; yy < yLen; yy++) {
	int regulators =  numRegulators(yy);
	if(regulators > 0){
	    uvec mask = maskCol(find(mask2.col(0) == (yy +1) )) - 1; 	
	    arma::mat xCurr = x.cols(mask);
	    center(xCurr);
	    arma::colvec yCurr= y.col(yy);
	    yCurr = yCurr - mean(yCurr);
	    //arma::mat corr = arma::cor(xCurr, yCurr); // the output are in column
	    vec small(xCurr.n_cols);
	    small.fill(.001); 
	    mat V_w = ((trans(xCurr)  * xCurr)+  diagmat(small)).i();
	    vec m_w = V_w * (trans(xCurr) * yCurr);
	    arma::uvec nz = sort_index(abs(m_w.col(0)), "descend");
	    nz = nz(span(0, regulators -1));
	    xCurr = xCurr.cols(nz);
	    small =zeros(regulators);
	    small.fill(.001); 
	    V_w = ((trans(xCurr)  * xCurr)+  diagmat(small)).i();
	    m_w = V_w * (trans(xCurr) * yCurr);
	    vec residue = yCurr - ( xCurr * m_w);
	    double varresidue = var(residue);
	    double vary = var(yCurr);
	    ev(yy) = 1 - varresidue/vary;
	}
    }
    return ev; 
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::mat cholGivensRotateLeft(arma::mat TT, arma::mat inp = arma::mat() ) 
{
    mat TTtemp = TT; 
    int n = TT.n_cols;
    mat u_in = eye(n, n); mat u1 = u_in;
    u1.cols(span(0,n-2)) = u_in.cols(span(1, n-1)); u1.col(n-1) = u_in.col(0); 
    TT = TT * u1;
    //uvec inx(n); inx(0) = n-1; inx(span(1,n-1)) = linspace<uvec>(0, n-2 , n-1);
    //TT = TT.cols(inx); 
    //TT = TT.rows(inx); 
    mat::fixed<2, 2>  ui;
    double r, a, b, t, c, s;
    double u;
    for (int k = 0; k <= n-2; k++) { //givens rotation wiki
	a = TT(k,k); b = TT(k+1, k);
	/* if (!(b==0)){ 
	    if(a==0) {
		c = 0; s = copysign(1.0, b);  r = abs(b);
	    }else if(abs(b) > abs(a)){
		t = a/b;
		u = copysign(sqrt(1+t*t), b);
		s = -1.0/u;
		c = -s*t;
		r = b*u;
	    }else{
		t = b/a;
		u = copysign(sqrt(1+t*t), a);
		c = 1.0/u;
		s = -c*t;
		r = a*u;
	    } */
	    r = hypot(TT(k,k), TT(k+1, k)); 
	    c =  TT(k,k)/r; 
	    s = - TT(k+1,k)/r;
	    ui(0,0) = c; ui(1,1) = c;
	    ui(1,0) =  s; ui(0,1) = -s;
	    TT(k,k) = r; TT(k+1, k) =0;
	    TT(span(k,k+1), span(k+1,n-1) )  =   ui * TT(span(k,k+1), span(k+1,n-1) ) ; 
	    //d = ABS(sqrt(TT(k, k) * TT(k, k) + TT(k+1, k)* TT(k+1, k)));
	    //if(d < 0){
		//cout << "d is negative \t " <<d <<"\t" <<  TT(k, k) * TT(k, k) + TT(k+1, k)* TT(k+1, k) << "\t"<< sqrt(TT(k, k) * TT(k, k) + TT(k+1, k)* TT(k+1, k))<< endl;
		//exit(1);
	    //}
	    //ui(0,0) = ui(1,1) = TT(k, k)/d ; //cos=TT(k, k)/d;
	    //ui(0,1) =  TT(k+1,k)/d; ui(1,0) = -ui(0,1); //sin= TT(k+1,k)/d;
	    //TT(k,k) = d; TT(k+1, k) =0;
	    //TT(span(k,k+1), span(k+1,n-1) )  =   ui * TT(span(k,k+1), span(k+1,n-1) ) ; 
	//}
    }
    TT(n-1, n-1) = ABS(TT(n-1, n-1));
    uvec yy = find( TT.diag() < 0);
    /* if (yy.size() > 0) {
	cout << "TT  while add is negative \t " << TT << endl;
	cout << TTtemp << endl;
	exit(-1);
    }
    TTtemp = chol(inp);
    if(sum(sum(TTtemp - TT)) > .05){
	cout << sum(sum(abs(TTtemp)))<< "\t" << sum(sum(abs(TT)))  << endl;
	cout << "mismatch in  rotation TT \t "  << endl;
	cout << TT << endl;
	cout << TTtemp << endl;
	exit(-1);
    } */
    return TT;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::mat cholAddCol(arma::mat TT, arma::vec xx, arma::mat inp = arma::mat()) 
{
    //cout << xx.size() << " inside " << TT.n_cols<<  endl;
    mat TTtemp = TT; 
    int n = TT.n_cols+1;
    TT.resize(n,n);
    //TT(span(1, n-1), n-1) = zeros(n-1); 
    double d;
    for (int i = 0; i <= n-1; i++) {
	if (i==0) {
	    TT(0, n-1) = xx(0)/ TT(0,0);
	}else{
	    TT(i, n-1) = xx(i) - sum(TT( span(0, i-1), n-1).t() * TT( span(0,i-1), i));
	    if( i < n-1){
		TT(i, n-1 ) = TT(i, n-1)/TT(i,i);
	    }else{
		d =   sqrt(TT(n-1,n-1));  TT(n-1,n-1) = (d>0) ? d:(-1*d);
		//if(TT(n-1 , n-1) < 0 ) cout << "erorr ******************* "<< TT(n-1, n-1) << endl; 
	    }
	}
    }
    /* uvec yy = find( TT.diag() < 0);
    if (yy.size() > 0) {
	cout << "TT is negative \t " << TT << endl;
	cout << TTtemp << endl;
	exit(-1);
    }
    TTtemp = chol(inp);
    if(sum(sum(TTtemp - TT)) > .01){
	cout << sum(sum(abs(TTtemp)))<< "\t" << sum(sum(abs(TT)))  << endl;
	cout << "mismatch in  TT at Add col \t " << endl;
	cout << TT << endl;
	cout << TTtemp << endl;
	exit(-1);
    } */
    //cout<< "complete" <<endl;
    return TT;
}

//arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
    //int ncols = sigma.n_cols;
    //arma::mat Y = arma::randn(n, ncols);
    //return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
//}
//
//vec logistic(vec x)
//{
    //return(1/(1 +exp(-(x))));
//}
//double logistic(double x)
//{
    //return(1/(1 +exp(-(x))));
//}
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double pinvgammaC(double q, double shape,double scale, int lower_tail = 1, int log_p = 0) 
{
    int lower_tailCurr = lower_tail?0:1;
    return(R::pgamma(1/q, shape, 1/scale, lower_tailCurr, log_p));
}
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double rinvgammaC(double shape, double scale) 
{
    return(1/R::rgamma(shape, 1/scale));
}

//double rinvgammaTRNG(double shape, double scale) 
//{
    //return(1/R::rgamma(shape, 1/scale));
//}


//double dnorm(double x){
    //return((exp(-(sqr(x))/2.0))/sqrt(2.0*M_PI));
//}
double dnorm(double x, double mu=0.0, double sigma=1.0){
    return((exp(-(sqr(x - mu))/(2.0*sqr(sigma))))/( sigma* ABS(sqrt(2.0*M_PI))));
}

void logistic_variable_selection_mcmc( arma::mat x, arma::vec y, arma::vec& gamMa, double& cc, double& sigma2 , double lambdal, double  nu)
{
    //x = num_samples x num_feature, y = num_samples x 1
    int  p_curr, k;
    double gamMa_snpCurr, tempd;
    p_curr = gamMa.size();
    int n = x.n_cols; 
    double small = .01;
    uvec nonzero, zero, nonzero_curr;
    vec  gamMa_rate, randDraw;
    mat xTx; colvec xTy;
    xTy = x.t() * y;
    xTx = x.t()* x + small * eye(p_curr, p_curr) ;
    nonzero_curr = find(gamMa > 0);
    nonzero = nonzero_curr;	
    int qq_curr  = nonzero_curr.size(); int qq = qq_curr;
    double gg;
    vec mm, WW;
    mat TT;
    double y2 = sum(square(y));
    vec  beTaCurr(p_curr);
    if(qq_curr == 0 ){
	TT.set_size(1,1); TT.zeros();
	gg  =  -(n + nu)/2.0 * log(nu * lambdal + y2);
    }else {
	TT = chol( x.cols(nonzero_curr).t() * x.cols(nonzero_curr)+ small * eye(qq_curr, qq_curr));
	WW = solve(trimatl(TT.t()), xTy(nonzero_curr)); 
	gg = -qq_curr/2.0 * log(1 + cc) - (n + nu)/2.0 * log (nu * lambdal + y2  - sum(square(WW))/(1 + 1/cc) );
    }

    mat TT1, TT_inv,VV, TT_curr = TT;
    colvec xbeTa(n), residue(n);
    vec xTx_snp;  
    vec gamMa_inclusion = zeros(p_curr);
    vec logit_freqCurr = zeros(p_curr);
    vec gamMaChange(p_curr); gamMaChange = gamMa;
    vec gamMa_prob(p_curr); gamMa_prob.fill(.5);
    double gg1, gg0;
    // nonzero need to be stored and when updating with unsorted nonzero and 1 chol at a time 
    zero = find(gamMa==0); 
    gamMa_rate = zeros(p_curr);
    double bf;

    for(int covar = 0;  covar < p_curr; covar++){
	int TTSize = TT_curr.n_cols;
	if((TTSize!= nonzero.size() )  && (TTSize > 1) ){
	    exit(0);
	}
	if(gamMa(covar) ==1){
	    qq_curr = qq-1;
	    gamMaChange(covar) = 0;
	    if(qq_curr > 0 ){
		nonzero_curr= nonzero(span(1, qq_curr));  
		nonzero(qq-1) = nonzero(0); nonzero(span(0, qq-2)) = nonzero_curr; 
		//TT_curr = cholGivensRotateLeft(TT_curr); 
		TT_curr = cholGivensRotateLeft(TT_curr, xTx(nonzero, nonzero)); 
		gg1 = gg;
		WW = solve(trimatl(TT_curr(span(0,qq_curr-1), span(0, qq_curr-1) ).t()), xTy(nonzero_curr)); 
		gg0 = -qq_curr/2.0 * log(1 + cc) - (n + nu)/2.0 * log (nu * lambdal + y2  - sum(square(WW))/(1 + 1/cc) );
	    }else{
		nonzero_curr.clear();
		gg0 =  - (n + nu)/2 * log(nu * lambdal  + y2 );
	    }	    
	}else{
	    gamMaChange(covar) = 1;
	    qq_curr = qq+1;
	    gg0 = gg;
	    nonzero_curr = nonzero;
	    nonzero_curr.resize(qq_curr); nonzero_curr(qq_curr-1) = covar;//  = find(gamMaChange > 0);
	    if(qq >  0){
		xTx_snp = xTx.col(covar);
		TT1 = cholAddCol(TT_curr, xTx_snp(nonzero_curr), xTx(nonzero_curr, nonzero_curr));
	    } else{
		TT1.set_size(1,1); tempd = sqrt(xTx(covar,covar)); TT1(0,0) = (tempd > 0) ? tempd : (-1*tempd);
	    }
	    WW =  solve(trimatl(TT1.t()), xTy(nonzero_curr)); //check there is no transponse
	    gg1 = -qq_curr/2.0 * log(1 + cc) - (n + nu)/2.0 * log (nu * lambdal + y2  -  sum(square(WW))/(1 + 1/cc) );
	}

	bf = logistic(gg1 - gg0); 
	randDraw = randu(2);
	gamMa_rate(covar) = logistic(gg1 - gg0 + log(gamMa_prob(covar)) - log(1 - gamMa_prob(covar)));	
	if( (bf > .5 ) && (bf > randDraw(0))) {
	    //randDraw = randu(1);
	    gamMa_snpCurr = (gamMa_rate(covar) > randDraw(1)) ?1:0;
	}else gamMa_snpCurr = 0;

	if(gamMa_snpCurr != gamMa(covar) ) {
	    gamMa(covar)   = gamMa_snpCurr; 
	    if(gamMa_snpCurr==1){
		nonzero = nonzero_curr; 
		gg = gg1;
		qq = qq_curr;
		TT_curr = TT1;
	    }else{
		nonzero = nonzero_curr;
		gg = gg0;
		qq = qq_curr;
		if(qq_curr >0){
		    TT_curr = TT_curr(span(0,qq-1), span(0,qq-1));
		}else TT_curr = zeros<mat>(1,1);
	    }
	}else gamMaChange(covar)   = gamMa_snpCurr;
    }
    if(qq > 0) {
	TT_inv = solve(trimatu(TT_curr), eye(qq,qq)); 
	VV = cc*sigma2/(cc +1) * (TT_inv * TT_inv.t() );
	mm = VV * xTy(nonzero)/sigma2;
	beTaCurr.zeros();
	beTaCurr(nonzero) = mvrnormArma(1, mm, VV);
	xbeTa  = x.cols(nonzero) * beTaCurr(nonzero);
	cout << 318 << endl;
	if(y.size() != xbeTa.size() ) {
	    cout << " error y_size not equal  to xbeTa size"  << "\t " << y.size() << "\t" << xbeTa.size() << "\t" << nonzero.t() << endl;
	}
	residue = y - xbeTa; 
	cout << 320 << endl;
    }else{
	beTaCurr.zeros(); xbeTa.zeros();
	residue = y;
    }
    double 	  residue2sum = sum(square(residue)); double xbeta2sum = sum(square(xbeTa));   

    sigma2 = rinvgammaC((qq + nu+ n +1)/2.0, 1/2.0 * (residue2sum + nu * lambdal));     
    cc = rinvgammaC((qq+1)/2.0, (xbeta2sum + n) /2.0);
}


using namespace arma;		// shorthand
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
Rcpp::List GOAL(arma::mat x, arma::mat y, arma::mat feature, arma::mat pairFeature, arma::umat mask2, arma::colvec alpHa,   arma::colvec gamMa, bool estimate_alpha, bool estimate_beta, arma::colvec B_inv_alpHa, int itermax, int  thin, int burnIn,  int threads, bool verbose=true, bool balance = true, bool logistic_variable_selection = true, bool oversample = false,bool  sample_logit = false, bool use_raoblackwell =true, int num_logit_train = 1000, int negtrain_num =1000, double ratio = 1, double gamMa_thres = .9, double beTa_thres=.1, double regulator_prior = 1.0/300.0, int accIter = 1, double rho = 100.0, double prior = .1)
{
    //regulator_prior = 0 : implies no re-estimation of alpha. num_logit_train = -1: no threshold is applied.
    if(gamMa_thres ==0) {
  cout << "gamMa_thres =0" << endl;
	exit(0);
    }
    #ifdef _OPENMP
    if ( threads > 0 )
    omp_set_num_threads( threads );
    REprintf("Number of threads=%i\n", omp_get_max_threads());
    #endif
    double eQTL_threshold = 1; // 1- implies eQTL_threshold is ignored 
    double regulator_threshold = 0; // 0- implies regulator_threshold is ignored
    bool dual_filter = false;
    bool usebeTaprob = beTa_thres > 0 ? true : false; 
    bool Berger=true;
    std::clock_t t1, t2, t3; 
    t2 = std::clock();
    float startOfepi_eQTL = float(t2);
    int localmaxiter =1; int bootstrap=50/localmaxiter; 
    int num_out_samples = floor(((double)(itermax - burnIn))/((double) thin));
    int p=mask2.n_rows,  n=y.n_rows, numGenes=y.n_cols; //, numSamp_logit = n/2, ;
    //cout << 454 << endl;
    mat feature1 = join_rows( ones(p,1), feature);
    //cout << 456 << endl;
    //mat feature1 = feature; // ones is added in logistic part
    int pairf = pairFeature.n_cols; int f1 = feature1.n_cols ; int f = f1 + pairf;
    int num_logit_postrain = (numGenes <  40)?  5*numGenes: 1000;
    mat beTa_mat(p, num_out_samples), alpHa_mat(f, num_out_samples), tau_mat(f, num_out_samples), gamMa_mat(p, num_out_samples), gamMa_prob_mat(p, num_out_samples), bf_mat(p, num_out_samples),  gamMa_rate_mat(p, num_out_samples),  gamMa_rateAcc(p, accIter), logistic_select_mat(f-1, num_out_samples), 
    ev_mat(numGenes, num_out_samples), cc_mat(numGenes, num_out_samples), sigma2_mat(numGenes, num_out_samples);
    colvec w(p), beTa(p), beTaAcc(p), beTaABS(p), gamMaAcc(p), gamMaABS(p), gamMa_rateAll(p);
    colvec alpHaAdj= alpHa; 
    uvec nonzero_logitInx(p); 
    gamMa_rateAcc.fill(.5);
    vec gamMa_prob, theta, aTf;
    vec ptheta(p);
    //y = y - (ones(n) * mean(y)) ;
    center_norm(y);
    center_norm(x);
    double LDthres = 0.9 * n-1; //LD threshold in terms of Pearson correlation coefficient. 
    //x = x - (ones(n) * mean(x)) ;
    colvec xTy(p); 
    field<mat> TT(numGenes); 
    field<mat> xTx(numGenes); 
    field<uvec> mask(numGenes); 
    field<field<uvec> > LD(numGenes); 
    field<vec> bf_LD(numGenes); 
    vec gg(numGenes); gg.zeros(); 
    uvec gene1_start(numGenes) , gene1_end(numGenes); gene1_start.zeros(); gene1_end.zeros();
    vec cc(numGenes); cc.fill(200); double cc_epi =10; // g-priors
    vec evAll(numGenes); evAll.zeros();
    int count  = 0;
    colvec logit_freq, geneEpi_freq, logit_freqtop, logit_freqbottom,  logit_freq_nonzero,  zz, bf(p), bf_single(p), tau(4);
    tau << 0 << 1 << 1<< 0;
    bool useGibbs = false;



    Rcpp::NumericVector	 ret1(num_logit_postrain), ret2(4*num_logit_postrain);
    logit_freq = ones(p); geneEpi_freq = ones(p); logit_freqtop = ones(p); logit_freqbottom = ones(p);  logit_freq_nonzero = ones(p);  zz=zeros(p); 
    zz = gamMa - .5; 
    vec logistic_select(f-1); 
    if( logistic_variable_selection) {
	logistic_select.zeros();
	for (int finx = 0; finx < f-1; finx++) {
	    if(alpHa(finx + 1) !=0)  logistic_select(finx) = 1;
	}
	//uvec temp = find(alpHa  !=0.0);
    }else logistic_select.ones();
    mat featureCat(p, f); 
    //cout << 493 << endl;
    if(pairf > 0) {
	//if( pairf ==1) featureCat.col(f) = pairFeature;
	//else 
	featureCat.cols(span(f1, f-1 )) = pairFeature;
    }
    //cout << 499 << endl;
    mat diag_w, V_w; vec m_w(p);
    if(beTa.size() !=p) beTa = zeros(p);

    mat feature_nonzero; 
    vec fTalpHa(p), geneEpiTtau(p);
    int p_nonzero;
    double nu1 = log( regulator_prior/ ( 1 - regulator_prior));
    // using regulatory prior 
    vec B_inv_m_w_prior= zeros(f);  // B_inv_m_w_prior(0) = B_inv_alpHa(0)  * logit(regulator_prior); // B-1 b page 6. polson et. al 	
    //alpHa(0) = logit(regulator_prior);


    // epi hyper-parameter 
    double sigma2_epi = .01, lambdal_epi = .01, nu_epi = 100;
    double a11, a10; // a10 = p(gamMa = 1| theta =0)
    a11 = prior;
    a10 = prior/rho;
    cout << " a11 \t"  << a11 << "\t a10 \t" << a10 << endl;  
    gamMa_prob =  zeros(p); gamMa_prob.fill(0.5);//this is theta
    //double bf_max = -logit(exp(regulator_prior * ( log(a11))  + (1-regulator_prior) * log(a10)));
    //double avg_a1 = exp(a11 * log(a11) + (1-a10) * log(a10)) ; 
    //double bf_max = max(logit(.75) - logit(avg_a1), min(logit(.9) - logit(avg_a1) , -logit(a10)/2)) ;
    double avg_gamMa_prob = 0.75;
    double avg_enh_a1 = exp(avg_gamMa_prob * log(a11) + (1-avg_gamMa_prob) * log(a10)) ;
    double bf_max = 19 * (1-avg_enh_a1)/(avg_enh_a1); /// at enhancer probability = .75 and bf_max, posterior probability = 0.95
    cout << "bf_max" << bf_max <<endl;
    double bf_min = -.01; double bf_thres = -12;
    cout << "bf_min" << bf_min <<endl;

    bf_max = log(bf_max); //bf_max is compared at log stage


    //# needed only  once
    RNGScope scope ;
    Environment base("package:base");
    Function sample = base["sample"];
    Rcpp::List out;
    bool display_progress = true;
    Progress progress(itermax, display_progress); 
    double small = .01;
    rowvec varY; varY =  var(y);
    rowvec sigma2,  lambdal;
    sigma2 = lambdal = 1*varY/9; // estimate of sigma2 variance of residue
    rowvec y2 = sum(square(y), 0);
    double minR = 0, maxR = 50, by=.5; 
    vec nu =  zeros(numGenes);
    double seqT, seqTmin, nuCurr,  bb;
    uvec maskCol = mask2.col(1);
    uvec LDCol = mask2.col(2);
    // calculation of bias of logistic regression
    mat featureMean(f,numGenes);

    // for theta regression 
    double sigma2_l =1,  sigma_l =1, lambdal_l = 0, nu_l =0, cc_l =10.0;
    int gamMa_initialized =0; 
    if( sum(gamMa) > 0) gamMa_initialized =1; 


    // initialize

    gamMa_rateAll =  zeros(p); //this is theta
    for (int gene = 0; gene < numGenes; gene++) {
	int p_curr, qq, qq_curr, k;
	uvec temp=maskCol(find(mask2.col(0) == (gene +1) ));  
	temp = temp-1;
	mask(gene) = temp; 
	p_curr = mask(gene).size();
	if(gene >0)
	gene1_start(gene) = gene1_end(gene-1) + 1; 
	gene1_end(gene) = gene1_start(gene) + p_curr - 1;
    }
    #pragma omp parallel for  schedule(dynamic)   
    for (int gene = 0; gene < numGenes; gene++) {
	seqTmin = 1; nuCurr = 0; 
	for (int ii = 0; ii < (maxR/by + 9000); ii++) {
	    if(ii < maxR/by) bb = minR + by*ii;
	    else bb = exp( -30.0 + 60.0/9000.0 *   (double) (ii - maxR/by));
	    seqT = pinvgammaC(varY(gene), bb/2.0, lambdal(gene)*bb/2.0) - .9;
	    if(( seqT< seqTmin) && (seqT > 0)){
		seqTmin = seqT; nuCurr = bb;
	    }
	}
	if(!Berger) nu(gene) = nuCurr;
	else nu(gene) =0;

	int p_curr, qq, qq_curr, k;
	uvec  maskCurr =  mask(gene);
	//cout << "maskCurr "<< maskCurr.t() << endl;
	//uvec temp=maskCol(find(mask2.col(0) == (gene +1) ));  
	//mask(gene) = temp-1;

	//if(gene >0)
	//gene1_start(gene) = gene1_end(gene-1) + 1; 
	//gene1_end(gene) = gene1_start(gene) + p_curr - 1;

	p_curr = mask(gene).size();
	if(p_curr > 0){
	    uvec nonzero, zero, nonzero_curr, maxinx;
	    //cout << 551<<endl;
	    featureCat(span(gene1_start(gene), gene1_end(gene)), span(0, f1-1)) = feature1.rows(mask(gene));
	    //cout << 553<<endl;
	    featureMean.col(gene) = trans(mean(featureCat.rows(span(gene1_start(gene), gene1_end(gene)) ))); 
	    //cout << 555<<endl;
	    xTy(span(gene1_start(gene), gene1_end(gene))) = x.cols(mask(gene)).t() * y.col(gene);
	    //cout << 557<<endl;
	    xTx(gene) = x.cols(mask(gene)).t()* x.cols(mask(gene)) + small*eye(p_curr, p_curr);
	    vec xTxdiag = xTx(gene).diag();
	    //cout << 558<<endl;
	    vec gamMaCurr = gamMa(span(gene1_start(gene),gene1_end(gene)));
	    vec xTy_gene = xTy(span(gene1_start(gene),gene1_end(gene)));
	    //cout << 559<<endl;
	    nonzero_curr = find(gamMaCurr > 0); 
	    int qq_curr  = nonzero_curr.size();
	    //cout << 550<<endl;
	    double ccCurr = cc(gene); vec WW; mat xTxcurr;


	    //double gg0Curr =  -(n + nu(gene))/2.0 * log(nu(gene) * lambdal(gene) + y2(gene));
	    double gg0Curr =  0;
	    double R2;
	    if(qq_curr == 0 ){
		TT(gene).set_size(1,1); TT(gene).zeros();
		//cout << 552 << "\t" << gg.size()  << "\t" << lambdal.size() << "\t" << gene << "\t" << nu.size() << "\t" << y2.size()<<endl;
		gg(gene)  = gg0Curr;
		//cout << 554<<endl;
	    }else {
		//TT(gene) = chol( x.cols(maskCurr(nonzero_curr)).t() * x.cols(maskCurr(nonzero_curr)) + small * eye(qq_curr, qq_curr)); 
		xTxcurr = (xTx(gene))(nonzero_curr, nonzero_curr); 
		TT(gene) = chol( xTxcurr);
		//cout << 556 << nonzero_curr.t() <<endl;
		WW = solve(trimatl(TT(gene).t()), xTy(maskCurr(nonzero_curr))); 
		//cout << 558<<endl;
		R2 = sum(square(WW))/y2(gene); 
		gg(gene) = (n - qq_curr -1.0)/2.0 * log(1 + ccCurr) - (n -1.0 )/2.0 * log(1.0+ccCurr*(1.0-R2));
		//gg(gene) = -qq_curr/2.0 * log(1 + ccCurr) - (n + nu(gene))/2.0 * log (nu(gene) * lambdal(gene) + y2(gene)  - sum(square(WW))/(1 + 1/ccCurr) );
		//cout << 561<<endl;
	    }

	    int ldsize = p_curr;
	    field<uvec> LDgene(ldsize,1);
	    vec bf_LDgene(ldsize);
	    int ldStart=0, ldEnd = 0, ldCnt = 0;
	    uvec tempLD = LDCol(maskCurr);
	    int lastld = tempLD(0);
	    uvec maskLD;
	    vec bf_temp;
	    int ldSize;
	    double ggA, fac;
	    mat TTtemp;
	    for (int ldinx = 0; ldinx < p_curr; ldinx++) {
		//ldsize = LDgene.size(); 
		//if(ldCnt > ( ldsize - 2) ) LDgene.resize(2*ldsize());
		ldEnd = (ldinx < p_curr -1) ? (ldinx - 1) : ldinx;
		if((lastld !=tempLD(ldinx))|| (ldinx == p_curr -1)|| ((tempLD(ldinx) < 0)&& ((ldEnd - ldStart) >= 2  ))){ // negative values imply singletons // enter if block only at end of LD block.
		    maskLD = maskCurr(span(ldStart, ldEnd)) -gene1_start(gene);
		    //cout << "maskLD " << maskLD.t()<< " lasld"<<lastld << "tempLD(ldinx)) " << tempLD(ldinx)   << endl;
		    //cout << xTx(gene).n_cols << " " << maskCurr(ldStart) << " " << maskCurr(ldEnd) <<  " " << ldStart << " " << ldEnd << endl;
		    xTxcurr = (xTx(gene))(maskLD, maskLD); 
		    TTtemp = chol( xTxcurr);
		    ldSize = ldEnd -ldStart + 1;
		    WW = solve(trimatl(TTtemp.t()), xTy(maskLD + gene1_start(gene))); 
		    R2 = sum(square(WW))/y2(gene);
		    //non-penalized for length version of bayesian factor qq is removed
		    ggA = (n  - 1.0)/2.0 * log(1.0 + ccCurr) - (n -1.0 )/2.0 * log(1+ccCurr*(1.0-R2));
		    bf_LDgene(ldinx) = exp(ggA);
		    //ggA = -gg0Curr - (ldEnd -ldStart + 1)/2.0 * log(1 + ccCurr) - (n + nu(gene))/2.0 * log (nu(gene) * lambdal(gene) + y2(gene)  - sum(square(WW))/(1 + 1/ccCurr) );
		    //R2_vec = square(xTy(maskLD))/(xTxdiag(maskLD) * y2(gene));
		    bf_temp =  (n - 2.0)/2.0 * log(1.0 + ccCurr) - (n -1.0 )/2.0 * log(1+ccCurr*(1.0- (square(xTy_gene(maskLD))/(xTxdiag(maskLD) * y2(gene))) )); //previously it was not exponential  
		    //bf_temp =   -gg0Curr - 1.0/2.0* log(1 + ccCurr) - (n + nu(gene))/2.0 * log (nu(gene) * lambdal(gene) + y2(gene)  - (square(xTy( maskLD))/((1 + 1/ccCurr)*xTxdiag(maskLD) )) );
		    //cout << "bf_temp" << bf_temp <<endl;
		    // previously bf_temp was normalized
		    /* bf_temp = bf_temp/sum(bf_temp) * ggA; */
		    maxinx = find(bf_temp > bf_max);
		    if(maxinx.size() >0)
		    bf_temp(maxinx) = ones<vec>(maxinx.size()) * bf_max;
		    bf_temp = exp(bf_temp);
		    bf_single(maskLD + gene1_start(gene)) = bf_temp;
		    //cout << "bf_single ld " << bf_single(maskLD + gene1_start(gene)) <<endl;
		    LDgene(ldCnt, 0) = maskLD; 
		    ldStart = ldinx;
		    ldCnt++;

		}
		lastld = tempLD(ldinx);

	    }
	    //cout << ldCnt << " gene "  << endl;
	    LD(gene) = LDgene.rows(0, ldCnt -1);
	    bf_LD(gene) = bf_LDgene;
	    //cout << LDgene(10) << endl; 

	    //bf_single(span(gene1_start(gene), gene1_end(gene)))  =   - gg0Curr -1.0/2.0* log(1 + ccCurr) - (n + nu(gene))/2.0 * log (nu(gene) * lambdal(gene) + y2(gene)  - (square(xTy(span(gene1_start(gene), gene1_end(gene))))/((1 + 1/ccCurr)*xTxdiag )) );
	    //cout << gene << "bf_single "  << bf_single(span(gene1_start(gene), gene1_end(gene))).t() << endl;
	    //cout << gene << "bf_LD "  << bf_LD(gene).t() << endl;

	    for (int i = gene1_start(gene); i < gene1_end(gene); i++) {
		if(bf_single(i) != bf_single(i) ){
		    cout << gene << "nan generated " << i << endl;
		    cout << gg0Curr <<  " " << ccCurr << "  " << xTy(i) <<endl;
		}
	    }
	}
    }
    //cout<< "bf many lines " << bf <<endl;
    //cout << "gg" << gg.t()<< endl;
    //cout << "nu" << nu.t() << endl;
    //cout << "gene_start " << gene1_start.t() <<endl;
    //cout << "gene_end " << gene1_end.t() <<(endl;
    aTf = featureCat * alpHaAdj;
    double sumAlpHa = sum(alpHa);
    //double prior =  (p > (40* numGenes))?  ((20.0*numGenes)/p ): 0.05;
    //double prior =  (p > (100* numGenes))?  ((50.0*numGenes)/p ): 0.01;

    // switching model

    double theta_bias = 0.0;
    double theta_odds = regulator_prior/(1-regulator_prior); // sum(gamMa_prob)/(p - sum(gamMa_prob));
    for (int iter = 0; iter < itermax; iter++) {
	if(verbose) 
	{ 
	    t1 = std::clock();
	    cout << "start of beta cycle \t" << iter << "\t" << ((float)t1 - (float) t2) / CLOCKS_PER_SEC << "s"<< endl;
	}
	if (Progress::check_abort() )
	return out;
	progress.increment();
	////cout << "running iteration " << iter << endl;
	if(estimate_beta){
	    logit_freq.zeros(); logit_freqtop.zeros(); logit_freqbottom.zeros();
	    if((iter %  accIter) ==0) {
		beTaAcc.zeros();
		gamMaAcc.zeros();
	    }
	    // find phi (std of theta)
	    //double phi = stddev(gamMa_prob - aTf);
	    //cout << "phi \t" << phi << endl;
	    ////phi = .5;
	    //vec randDraw1; randDraw1 = randu(1);
	    //dual_filter = (randDraw1(0) < .1) ? true : false; 
	    //double sigma2_eff = 1.0/ (1.0/(sigma2_l) + 1.0/square(sigma_0) ); double sigma_eff = sqrt(sigma2_eff);
	    /*if(false) 
	    {
		if(iter == 0) {
		    a11 = .8; a10 = 0.01;
		}else{
		    a11 = sum(gamMa_rateAll % gamMa_prob)/sum(gamMa_prob); 
		    a10 = sum(gamMa_rateAll % ( ones(p) - gamMa_prob) )/(p - sum(gamMa_prob)); 
		}
	    }else{
		a11 =  min(.99, .99 * numGenes * 5/sum(gamMa_prob));  
		a10 =  min(.99, .01 * numGenes * 5/(p - sum(gamMa_prob)));
	    }*/

	    #pragma omp parallel for schedule(dynamic)   
	    for (int gene = 0; gene < numGenes; gene++) {
		int  p_curr, qq, qq_curr, k;
		double gamMa_snpCurr, gamMa_probsnpCurr, tempd;
		double ccCurr = cc(gene), sigma2Curr = sigma2(gene);
		uvec  maskCurr =  mask(gene);
		bool update_bf=true;

		//#gibbs sampling;
		p_curr = maskCurr.size();
		if(p_curr > 0){
		    int topp, bottomm;
		    uvec sort_inx;
		    //Rcpp::NumericVector
		    ransampl_ws* ws = ransampl_alloc( p_curr );
		    vec sample_prob(p_curr);
		    double ret;
		    uvec bootstrapS(p_curr);
		    vec a1, a2; 
		    //vec ret(bootstrap); 
		    bootstrapS = linspace<uvec>(0, p_curr-1, p_curr);
		    double epiDiscount = 0; 
		    uvec nonzero, zero, nonzero_curr;
		    vec  bfCurr(p_curr), gamMa_rate, randDraw;
		    double bf_snpCurr, prob_proposal, proposal_ratio, accept;
		    int accept_count = 0, proposal_count =0;
		    vec mm, WW;
		    mat TT1, TT_inv,VV, TT_proposal, TT_curr = TT(gene); // chol( x.cols(maskCurr(nonzero_curr)).t() * x.cols(maskCurr(nonzero_curr)) + small * eye(qq_curr, qq_curr))
		    int gene_start = gene1_start(gene), gene_end = gene1_end(gene); 
		    colvec xbeTa(n), xbeTaM(n), residue(n);
		    mat xTx_gene = xTx(gene); vec xTx_snp_gene;  
		    vec  gamMaCurr = gamMa(span(gene_start, gene_end)); vec gamMa_inclusion = zeros(p_curr);
		    vec logit_freqCurr = zeros(p_curr);
		    vec gamMaChange(p_curr); gamMaChange = gamMaCurr;
		    //cout << "here" << gene << 181 << gene_start << "endl" << gene_end << endl;
		    vec  beTaCurr(p_curr), geneEpiCurr(4);
		    //vec alpHaCurr = alpHa;
		    vec featureMeanCurr = featureMean.col(gene);
		    //double expalpHaF = sum(trans(alpHa(span(1,f-1))) * featureMeanCurr(span(1,f-1)));
		    //if ((p_curr > (3*regulator_prior)) && (regulator_prior  > 0) ) alpHaCurr[0] = -log(p_curr/regulator_prior -1) - expalpHaF;
		    //vec gamMa_probCurr = logistic(featureCat.rows(span(gene_start,gene_end)) * alpHaCurr);
		    //alpHaCurr(0) = 0;
		    //cout << "alpHa" << alpHa.t() << endl; 
		    //cout << "gamMa" << sum(gamMa) << endl;
		    //cout << "initialized" << gamMa_initialized << endl;
		    //return out;


		    double gg1, gg0, gg_proposal, R2, gg_curr = gg(gene);
		    // nonzero need to be stored and when updating with unsorted nonzero and 1 chol at a time 
		    k = gamMaCurr.size(); zero = find(gamMaCurr==0); nonzero = find(gamMaCurr > 0); qq  = nonzero.size();
		    //cout << 763 << endl;
		    gamMa_rate = gamMa_rateAll(span(gene_start, gene_end));
		    //cout << 764 << endl;
		    nonzero_curr = nonzero;
		    // sample gamMa_prob 
		    vec gamMa_probCurr(p_curr);
		    gamMa_probCurr= gamMa_prob(span(gene_start, gene_end));
		    vec aTfCurr = aTf(span(gene_start, gene_end));
		    //cout << 763 << endl;
		    //vec Randn = randn<vec>(p_curr)*sigma_eff; // N(0, sigma_eff)
		    //vec Randn_0 = randn<vec>(p_curr)*sigma_0; Randn_0 = Randn_0 + thetamean * ones(p_curr); //  N(mu0, sigma0) 
		    double aTfbin;
		    //gamMa_probCurr  = gamMa_probCurr * (1-2*epiDiscount) + epiDiscount;

		    //double phi1;

		    if( qq > 0) {
			WW = solve(trimatl(TT_curr.t()), xTy(maskCurr(nonzero)));
			R2 = sum(square(WW))/y2(gene); 
			gg_curr = (n - qq -1.0)/2.0 * log(1 + ccCurr) - (n -1.0 )/2.0 * log(1.0+ccCurr*(1.0-R2));
		    }else gg_curr = 0;
		    for (int localiter = 0; localiter < localmaxiter; localiter++) {
			if((useGibbs) ||  ((iter ==0) && (gamMa_initialized ==0) )){
			    for(int snp = 0;  snp < p_curr; snp++){
				//cout << qq_curr<< " here "  << 243 << gene << endl;
				if(gamMaCurr(snp) ==1){
				    qq_curr = qq-1;
				    //TT1 = TT_curr;
				    gamMaChange(snp) = 0;
				    if(qq_curr > 0 ){
					nonzero_curr= nonzero(span(1, qq_curr));  //find(gamMaChange > 0);
					nonzero(qq-1) = nonzero(0); nonzero(span(0, qq-2)) = nonzero_curr; 
					//cout << "here" << gene << 198 << endl;
					TT_curr = cholGivensRotateLeft(TT_curr, xTx_gene(nonzero, nonzero)); //chol( x.cols(maskCurr(nonzero_curr)).t() * x.cols(maskCurr(nonzero_curr)) + small * eye(qq_curr, qq_curr));
					//cout << "here" << gene << 192 << endl;
					gg1 = gg_curr;
					WW = solve(trimatl(TT_curr(span(0,qq_curr-1), span(0, qq_curr-1) ).t()), xTy(maskCurr(nonzero_curr))); 
					//cout << "here" << gene << 19206 << endl;
					R2 = sum(square(WW))/y2(gene); 
					gg0 = (n - qq_curr -1.0)/2.0 * log(1 + ccCurr) - (n -1.0 )/2.0 * log(1.0+ccCurr*(1.0-R2));
					//cout << "here" << gene << 19208<< endl;
				    }else{
					//TT0 = zeros<mat>(1,1);
					nonzero_curr.clear();
					//cout << "here" << gene << 19210<< endl;
					gg0 = 0;
					//gg0 =  - (n + nu(gene))/2 * log(nu(gene) * lambdal(gene)  + y2(gene) );
					//cout << "here" << gene << 19213 << endl;
				    }	    
				}else{
				    //TT0 = TT_curr; 
				    //cout << "here" << gene << 19217 << endl;
				    gamMaChange(snp) = 1;
				    qq_curr = qq+1;
				    gg0 = gg_curr;
				    //cout <<"here 265 " <<  qq_curr <<nonzero.t() << endl;
				    nonzero_curr = nonzero;
				    nonzero_curr.resize(qq_curr); nonzero_curr(qq_curr-1) = snp;//  = find(gamMaChange > 0);
				    //cout << "here" << nonzero_curr.t() << " "  << gene << 19220 << " " <<qq_curr << "  " << TT1.n_cols<< endl;
				    if(qq >  0){
					xTx_snp_gene = xTx_gene.col(snp);
					TT1 = cholAddCol(TT_curr, xTx_snp_gene(nonzero_curr), xTx_gene(nonzero_curr, nonzero_curr));//chol( x.cols(maskCurr(nonzero_curr)).t() * x.cols(maskCurr(nonzero_curr)) + small * eye(qq_curr, qq_curr));
				    } else{
					TT1.set_size(1,1); tempd = sqrt(xTx_gene(snp,snp)); TT1(0,0) = (tempd > 0) ? tempd : (-1*tempd);
				    }
				    //cout << "here" << gene << 19222 << endl;
				    WW =  solve(trimatl(TT1.t()), xTy(maskCurr(nonzero_curr))); //check there is no transponse
				    R2 = sum(square(WW))/y2(gene); 
				    gg1 = (n - qq_curr -1.0)/2.0 * log(1 + ccCurr) - (n -1.0 )/2.0 * log(1.0+ccCurr*(1.0-R2));
				    //gg1 = -qq_curr/2.0 * log(1 + ccCurr) - (n + nu(gene))/2.0 * log (nu(gene) * lambdal(gene) + y2(gene)  -  sum(square(WW))/(1 + 1/ccCurr) );
				    ////cout << "here" << gene << 19226 << endl;
				}
				//cout << gg0 << " here " << gg1 << sum(square(WW)) << endl;
				//randDraw = randu(1);
				// sample gamMa_prob 

				aTfbin = logistic(aTfCurr(snp));
				randDraw = randu(3);
				//phi1 = ABS(aTfCurr(snp) - nu1);
				//gamMa_probsnpCurr = gamMa_probCurr(snp);
				gamMa_probsnpCurr = logistic(aTfCurr(snp));

				bfCurr(snp) =  logistic(gg1  - gg0  + logit(.05)); 
				//cout << 863 << "\t" << gg1 <<  "\t" << gg0  << "\t" << gamMa_probCurr(snp)  << "\t" <<  prior <<  "\t"<<  gg1 - gg0 + gamMa_probCurr(snp)+ log(prior/(1-prior)) << endl;
				//gamMa_rate(snp) = logistic(gg1 - gg0 + gamMa_probsnpCurr);
				//double pa1x = min(a11 * gamMa_probsnpCurr + a10 * (1 - gamMa_probsnpCurr) , prior);
				double pa1x;
				pa1x = (iter==0)?  a10: a11 * gamMa_probsnpCurr + a10 * (1 - gamMa_probsnpCurr);
				//cout << " pa1x \t" << pa1x << endl;
				gamMa_rate(snp) = logistic(gg1 - gg0 +   logit(pa1x));
				double max_rate = .99*aTfbin + .8*(1-aTfbin);
				//gamMa_rate(snp) = (aTfCurr(snp) > 0)?gamMa_rate(snp): min(gamMa_rate(snp), 0.9);

				//cout << gamMa_rate(snp) << endl;
				gamMa_snpCurr = 0;
				//gamMa_snpCurr = ((gamMa_rate(snp) >= randDraw(1)) & (bfCurr(snp) > randDraw(2)) ) ? 1:0;
				if(dual_filter){
				    if ((gamMa_rate(snp) >= randDraw(1)) && (aTfbin > randDraw(2)) ) gamMa_snpCurr = 1;
				}else{
				    //if((gamMa_rate(snp) >= randDraw(1)) && (bfCurr(snp) > randDraw(2)) ) gamMa_snpCurr = 1;
				    if((gamMa_rate(snp) * max_rate) >= randDraw(1)  ) gamMa_snpCurr = 1;
				}
				double pax1 =  gamMa_snpCurr* log(a11/a10)  + (1-gamMa_snpCurr ) * log( (1 - a11)/(1 - a10)) ;  
				// this may not be correct, should it be arthmetic mean or geometric mean, double pax1 =  gamMa_rate(snp)* log(a11/a10)  + (1-gamMa_rate(snp) ) * log( (1 - a11)/(1 - a10)) ;  
				//double temp1 = log(theta_odds) + pax1 + aTfCurr(snp); //  - logit(.99); // a prior can be introduced in gamMa_prob/theta
				double temp1 = - theta_bias + pax1 + aTfCurr(snp); //  - logit(.99); // a prior can be introduced in gamMa_prob/theta
				//if(temp1 > 8.0 ) temp1 = 8.0;
				//if(temp1 <  -8.0 ) temp1 = -8.0;
				gamMa_probCurr(snp) = (iter == 0)? gamMa_snpCurr : temp1;
				//cout  << gene << " snp " << snp << " prior " <<  prior <<  " pa1x "<<  pa1x << " pax1 "<<  pax1 << "gamMa_rate(snp)" << gamMa_rate(snp) << " gamMa " <<  gamMa_snpCurr << " gamMa_probCurr(snp) " << gamMa_probCurr(snp) << " logit(gamMa_rate) " << gg1 - gg0 +  logit(prior) + logit(pa1x) <<" bf "  << bfCurr(snp) << " aTf " <<  aTfCurr(snp) << " theta_bias "<< theta_bias<<  endl;
				//if( (gamMa_probCurr(snp)  == 1.0) |(gamMa_probCurr(snp) == 1.0)  ){
				    //cout << "error probCurr out of range " << endl;
				    //exit(-1);
				//}



				// sampling of theta
				/*double mu = aTfCurr(snp)* sigma2_eff/sigma2_l;
				if(ABS(aTfbin -  gamMaCurr(snp)) > 0.5 ){ // proposal N(0, s) 
				    prob_proposal = Randn(snp);
				    proposal_ratio = (gamMaCurr(snp) > 0) ? (logistic(prob_proposal )/logistic(gamMa_probsnpCurr)) : ((1 - logistic(prob_proposal))/(1- logistic(gamMa_probsnpCurr))); 
				    //proposal_ratio	= exp( (1 - gamMa_rate(snp)) * ( log(1- logistic(prob_proposal)) - log(1- logistic(gamMa_probsnpCurr)   )) +
				    //gamMa_rate(snp)  *( log(logistic(prob_proposal )) -  log(logistic(gamMa_probsnpCurr)) ));
				    accept = proposal_ratio * dnorm(prob_proposal, mu, sigma_eff )/dnorm(gamMa_probsnpCurr, mu, sigma_eff )*\
				    dnorm(gamMa_probsnpCurr, 0.0, sigma_eff)/dnorm(prob_proposal, 0.0, sigma_eff);
				    if((accept > randDraw(0))&& (prob_proposal > thetamin) && (prob_proposal <=  thetamax)){
					accept_count++;
					gamMa_probCurr(snp) = prob_proposal; //accept 
				    }
				}else{ //proposal N(aTf, s)
				    prob_proposal = Randn(snp) + mu;
				    proposal_ratio = (gamMaCurr(snp) > 0) ? (logistic(prob_proposal )/logistic(gamMa_probsnpCurr)) : ((1 - logistic(prob_proposal))/(1- logistic(gamMa_probsnpCurr))); 
				    //proposal_ratio	= exp( (1 - gamMa_rate(snp)) * ( log(1- logistic(prob_proposal)) - log(1- logistic(gamMa_probsnpCurr)   )) +
				    //gamMa_rate(snp)  *( log(logistic(prob_proposal )) -  log(logistic(gamMa_probsnpCurr)) ));
				    accept = proposal_ratio;
				    if((accept > randDraw(0))&& (prob_proposal > thetamin) && (prob_proposal <=  thetamax)){
					accept_count++;
					gamMa_probCurr(snp) = prob_proposal; //accept
				    }
				}
				//cout << "prob_proposal \t" << prob_proposal << "\t" <<  gamMa_probCurr(snp) << endl;
				proposal_count++; * / sampling of theta 

				//gamMa_probCurr(snp) = aTfCurr(snp); //accept
				//gamMa_snpCurr = (gamMa_rate(snp) >= randDraw(1))  ? 1:0;
				//cout << gene_start + snp <<"\t bfCurr\t" <<  bfCurr(snp) <<  "\t bfCurr\t " <<  gg0 - gg1 << "\t gamMa_prob \t " << gamMa_probCurr(snp) << "\t gamMa_rate\t "  << 
				//gamMa_rate(snp) <<" \t gamMa \t"   << gamMa_snpCurr<< "\t prior_prob \t"<<log(prior/(1-prior))<<  endl;
				/*if(((gene_start+ snp) > 32000)  && ((gene_start+ snp) < 35000)){
				    double GG1, GG0 ; int QQ_curr; mat tt1;
				    uvec NONzero; uvec temp12 = find(nonzero != snp);
				    NONzero = nonzero(temp12);
				    QQ_curr= NONzero.size();
				    cout << NONzero.t() << "\tsnp \t" << snp<<endl; 
				    tt1 = chol( x.cols(maskCurr(NONzero)).t() * x.cols(maskCurr(NONzero)) + small * eye(QQ_curr, QQ_curr)); 
				    WW =  solve(trimatl(tt1.t()), xTy(maskCurr(NONzero))); //check there is no transponse
				    GG1 = -QQ_curr/2.0 * log(1 + ccCurr) - (n + nu(gene))/2.0 * log (nu(gene) * lambdal(gene) + y2(gene)  -  sum(square(WW))/(1 + 1/ccCurr) ); 
				    QQ_curr = QQ_curr + 1;	
				    NONzero.resize(QQ_curr); NONzero(QQ_curr-1) = snp;//  = find(gamMaChange > 0);
				    tt1 = chol( x.cols(maskCurr(NONzero)).t() * x.cols(maskCurr(NONzero)) + small * eye(QQ_curr, QQ_curr)); 
				    WW =  solve(trimatl(tt1.t()), xTy(maskCurr(NONzero))); //check there is no transponse
				    GG0 = -QQ_curr/2.0 * log(1 + ccCurr) - (n + nu(gene))/2.0 * log (nu(gene) * lambdal(gene) + y2(gene)  -  sum(square(WW))/(1 + 1/ccCurr) );
				    cout << "check bfCurr \t"<< GG1 -GG0 <<endl;
				}*/

				//gamMa_rate(snp) = logistic(sum(geneEpiCurr.t() * tau) );
				//cout << gamMa_rate(snp) << endl;
				//gamMa_rate(snp) = logistic(gg1 - gg0 + log(gamMa_probCurr(snp)) - log(1 - gamMa_probCurr(snp)));	
				//bf_snpCurr = logistic( bfCurr(snp)); 
				//cout << "randDraw \t" << randDraw.t() << "gamMa_rate \t" << gamMa_rate(snp) << "\t bf_snpCurr" << endl;

				/*if( (bf_snpCurr >= eQTL_threshold ) || (bf_snpCurr >= randDraw(0))) {
				    if(gamMa_probCurr(snp) >= regulator_threshold){
					//cout << "inside" << endl;
					gamMa_snpCurr = (gamMa_rate(snp) >= randDraw(1)) ?1:0;
					//if(gamMa_snpCurr > 0) cout << "success" << endl;
				    }else  gamMa_snpCurr = 0;
				    if(dual_filter){
					gamMa_snpCurr = (gamMa_probCurr(snp) >= randDraw(1)) ?1:0;
				    }
				}else gamMa_snpCurr = 0;*/

				//cout <<gamMa_snpCurr << endl;

				//gamMa_snpCurr = (Rcpp::rbinom(1, 1, gamMa_rate(snp)))(0);
				//gamMa_snpCurr = (gamMa_rate(snp) <.5 )? 0:1;
				////cout << gamMa_rate(snp) << " here " <<  gamMa_snpCurr << " old " << gamMaCurr(snp)  <<  endl;
				if(gamMa_snpCurr != gamMaCurr(snp) ) {
				    gamMaCurr(snp)   = gamMa_snpCurr; 
				    if(gamMa_snpCurr==1){
					nonzero = nonzero_curr; 
					gg_curr = gg1;
					qq = qq_curr;
					TT_curr = TT1;
				    }else{
					nonzero = nonzero_curr;
					gg_curr = gg0;
					qq = qq_curr;
					if(qq_curr >0){
					    TT_curr = TT_curr(span(0,qq-1), span(0,qq-1));
					}else TT_curr = zeros<mat>(1,1);
				    }
				}else gamMaChange(snp)   = gamMa_snpCurr;
			    }
			    if(use_raoblackwell) gamMa_inclusion = gamMa_inclusion + gamMa_rate;
			    else gamMa_inclusion = gamMa_inclusion + gamMaCurr;
			}else{
			    int qqLD, qqLD_proposal, qq_proposal, qqLDLD, qqLDLD_proposal;
			    uvec nonzero_proposal, maxinx, mininx;
			    field<uvec> LDgene = LD(gene);
			    int ldcnt = LDgene.size();
			    uvec LDCurr, snps_inx;
			    vec aTf_bin = logistic(aTfCurr);
			    //vec pa1x_vec =   a11 * gamMa_probCurr + a10 * (1 - gamMa_probCurr);
			    vec pa1x_vec =   exp(log(a11) * aTf_bin + log(a10) * (1 - aTf_bin));
			    vec pa1x_vecneg = 1 - pa1x_vec;
			    vec bf_singleCurr = bf_single(span(gene_start, gene_end));
			    if(update_bf){
				vec xTxdiag = xTx_gene.diag();
				vec xTy_gene = xTy(span(gene_start, gene_end));
				vec bf_temp;
				for(int ldinx = 0;  ldinx < ldcnt; ldinx++){
				    snps_inx = LDgene(ldinx);
				    bf_temp =  (n - 2.0)/2.0 * log(1.0 + ccCurr) - (n -1.0 )/2.0 * log(1+ccCurr*(1.0- (square(xTy_gene(snps_inx))/(xTxdiag(snps_inx) * y2(gene))) )); 
				    maxinx = find(bf_temp > bf_max);
				    if(maxinx.size() >0)
				    bf_temp(maxinx) = ones<vec>(maxinx.size()) * bf_max;
				    mininx = find(bf_temp < bf_min);
				    if(mininx.size() >0)
				    bf_temp(mininx) = ones<vec>(mininx.size()) * bf_thres;
				    bf_temp = exp(bf_temp);
				    bf_singleCurr(snps_inx) = bf_temp;
				}
				 bf_single(span(gene_start, gene_end)) = bf_singleCurr;
				//cout << "bf_single ld " << bf_single(snps_inx + gene1_start(gene)) <<endl;
			    }
			    //cout << 1030 << endl;
			    //vec proposal_gamMa = logistic(logit(pa1x_vec) + bf_singleCurr);
			    vec proposal_gamMa = bf_singleCurr % pa1x_vec;
			    vec proposal_gamMaCurr;
			    vec bf_LDgene = bf_LD(gene);
			    //cout << 1031 << endl;
			    vec Randu = randu<vec>(p_curr); // N(0, sigma_eff)
			    double prior_proposal, prior_curr, trans_proposal, trans_curr, prior_changeNotLD, trans_changeNotLD, pa1x_neg, pa1xNotLD_neg;
			    //bfCurr = gamMaCurr;
			    bfCurr = bf_singleCurr;
			    uvec LDLDInx, nonzeroNotLD, nonzeroNotLD_proposal, zeroNotLD_proposal;
			    uvec LDLD, zeroLDLDInx; 
			    vec proposal_gamMaNotLD;
				//uvec nonzeroLDLD, nonzeroLDLDInx, nonzeroLDLD_proposal, zeroLD_proposal, zeroNotLD_proposal;
				//int cholEstimate;
			    //gamMa proposal
			    for(int ldinx = 0;  ldinx < ldcnt; ldinx++){
				//cout << 1118 << endl;
				//cholEstimate = 0;
				snps_inx = LDgene(ldinx);
				//cout <<"LDCurr " <<  LDCurr.t() << endl;
				//cout << 1032 << " ldinx "<<  ldinx<<  " dfds "<< snps_inx.t()<< endl;
				proposal_gamMaCurr = proposal_gamMa(snps_inx);
				pa1x_neg = prod(pa1x_vecneg(snps_inx));
				pa1xNotLD_neg = 1;
				double proposal_gamMaSum = sum(proposal_gamMaCurr);
				//cout << trans(bf_singleCurr(snps_inx)) << " " <<trans( pa1x_vec(snps_inx)) << endl;
				//cout << 1033 << " nonzeroLD_proposal "<<  nonzeroLD_proposal.t() << endl;
				uvec nonzeroLDInx =find( gamMaCurr(snps_inx) > 0.5); 
				uvec nonzeroLD = snps_inx(nonzeroLDInx); 
				uvec zeroLDInx =find( gamMaCurr(snps_inx) <= 0.5); 
				uvec zeroLD = snps_inx(zeroLDInx); 
				
				//cout << 1034 << " nonzeroLD_ "<<  nonzeroLD.t() << endl;
				//if gamMa proposal is same then accept
				accept = 0; qqLD = nonzeroLD.size() ;
				//cout << 1035 << endl;

				// if so include that in LD. 
				// indices of included variables not in LD
				int notqqLDLD = 0, cholEstimate =0;

				//cout << " before  rotate " <<  nonzero.t() << endl;
				//rotate all nonzeroLD 
				if((qq > 0) && (qq > qqLD) && (qqLD > 0 )){ // if qqLD == qq don't need to rotate
				    for (int rot = 0; rot < qqLD; rot++) TT_curr = cholGivensRotateLeft(TT_curr);
				    nonzero_curr.set_size(qq);
				    nonzero_curr(span(0,qq - qqLD - 1)) = nonzero(span( qqLD  , qq  -1 ));
				    //cout << 1038 << endl;
				    nonzero_curr(span(qq - qqLD , qq -1 )) = nonzero(span( 0  , qqLD  -1 ));
				    //cout << 1039 << endl;
				    nonzero = nonzero_curr; //rotated nonzero

				    nonzeroNotLD = nonzero(span(0,qq - qqLD - 1));
				}else{
				    nonzeroNotLD = nonzero;
				}

				int qqNotLD = qq - qqLD;
				if(  qqNotLD > 0){
				    //cout << 1150 << endl;
				    umat xTxTemp = ( abs(xTx_gene(nonzeroNotLD, snps_inx)) >= (LDthres ) );
				    //cout << 1151 << endl;
				    uvec LDLDInx = find(sum(xTxTemp, 1) > 0);
				    LDLD = nonzeroNotLD( LDLDInx);
				    //nonzero_remaining = nonzeroNotLD(find(1- LDLDInx ));
				    //cout << 1153 << endl;

				    qqLDLD = LDLD.size();
				    if(  qqLDLD > 0){
					//cout << 1157 << endl;
					vec RanduTemp = randu<vec>(LDLD.size()); // N(0, sigma_eff)
					//cout << 1159 << endl;
					pa1xNotLD_neg = prod(pa1x_vecneg(LDLD));
					proposal_gamMaSum = proposal_gamMaSum + sum(proposal_gamMa(LDLD) );
					proposal_gamMaNotLD = ( proposal_gamMa(LDLD))/( proposal_gamMaSum + (pa1xNotLD_neg * pa1x_neg));
				        zeroLDLDInx = find((proposal_gamMaNotLD - RanduTemp) < 0);// this indicies are with respect to LDLD
					//cout << "zeroLDLDInx " << zeroLDLDInx << endl;
					zeroNotLD_proposal = LDLD(zeroLDLDInx);  // nonzero snps in LD with current LDBlock that will be zero 
					//cout << "pa1xNotLD_neg " <<pa1xNotLD_neg << " proposal_gamMaNotLD " << proposal_gamMaNotLD  << "sum"<<   ( proposal_gamMaSum + (pa1xNotLD_neg * pa1x_neg)) << "zeroNotLD_proposal" << zeroNotLD_proposal <<  endl;
					notqqLDLD = zeroNotLD_proposal.size();
					if(notqqLDLD > 0){ 
					    uvec nonzeroNotLD_proposalInx = ones<uvec>(qqNotLD);
					    //cout << "LDLDInx " << LDLDInx.t() << " zeroLDLDInx \t " << zeroLDLDInx.t() <<  " nonzeroNotLD \t " << nonzeroNotLD.t() << " noqqLDLD " << notqqLDLD << " inx "<< LDLDInx(zeroLDLDInx) << endl;
					    nonzeroNotLD_proposalInx(LDLDInx(zeroLDLDInx)) = zeros<uvec>(notqqLDLD);
					    //cout << "	nonzeroNotLD_proposalInx" << nonzeroNotLD_proposalInx << endl;
					    nonzeroNotLD_proposal = nonzeroNotLD(find( nonzeroNotLD_proposalInx));
					    //cout << 1160 << endl;
					    cholEstimate  = 1;
					    accept = 0;
					}
				    }
				}
				proposal_gamMaCurr = proposal_gamMaCurr /( proposal_gamMaSum +  (pa1xNotLD_neg * pa1x_neg)); 
				uvec nonzeroLD_proposalInx = find((proposal_gamMaCurr - Randu(snps_inx)) >= 0); 
				uvec nonzeroLD_proposal = snps_inx(nonzeroLD_proposalInx); //sorted indices
				uvec zeroLD_proposalInx = find((proposal_gamMaCurr - Randu(snps_inx)) < 0); 
				uvec zeroLD_proposal = snps_inx(zeroLD_proposalInx); //sorted indices
				qqLD_proposal = nonzeroLD_proposal.size();
				//cout << gene << " snps "<< snps_inx.t() << " LD " <<  nonzeroLD.t() << " proposal " << nonzeroLD_proposal.t() << " trans_curr " << proposal_gamMaCurr.t() << " pa1x " << pa1x_vec(snps_inx).t() << "bf_single " << bf_singleCurr(snps_inx).t()   << endl; 
				//cout << "Randu(snps_inx)" << Randu(snps_inx).t() << " proposal_gamMaCurr" << proposal_gamMaCurr.t() << endl; 
				//cout << 1181 << endl;
				//if(qqLDLD > 0){
				    //if( (qqLDLD == qqLDLD_proposal) ){
					//cholEstimate=0; // no need check each bit because they were all one anyone change will change qqLDLD_proposal
					//cout << 1036 << endl;
				    //}
				//}
				if(cholEstimate ==0){
				    if(qqLD == qqLD_proposal ){
					if( sum(nonzeroLD_proposal == nonzeroLD) == qqLD) accept = 1; 
					//cout << 1037 << endl;
				    }
				}
				//cout << 101138<< " accept "<<  accept << endl;
				if(accept !=1){
				    if(cholEstimate==1){
					//cout << 1216 <<endl;
					// nonzero_proposal
					qq_proposal = qq - qqLD - notqqLDLD;
					qq_proposal = qq_proposal + qqLD_proposal;
					//cout << 1199 << endl;
					if(qq_proposal > 0){
					    nonzero_proposal = join_cols( nonzeroNotLD_proposal,  nonzeroLD_proposal);
					    //cout << " nonzeroLDLD_proposal " <<  nonzero_proposal.t() << endl; 
					    //cout << "chol "<< xTx_gene(nonzero_proposal, nonzero_proposal)  << endl;
					    TT_proposal = chol(xTx_gene(nonzero_proposal, nonzero_proposal));
					}else{
					    nonzero_proposal.clear(); 
					    TT_proposal.set_size(1,1); TT_proposal.zeros();
					}
					prior_changeNotLD =  prod(1- pa1x_vec( zeroNotLD_proposal))/prod( pa1x_vec( zeroNotLD_proposal));
					trans_changeNotLD  = prod( 1- proposal_gamMaNotLD(zeroLDLDInx)) / prod(  proposal_gamMaNotLD(zeroLDLDInx)); 

				    }else{
					//cout << 1235 <<endl;
					prior_changeNotLD = 1; trans_changeNotLD = 1;

					qq_proposal = qq - qqLD;
					if (qq_proposal > 0){
					    TT_proposal = TT_curr(span(0,qq_proposal -1), span(0, qq_proposal -1) );
					    //cout << 1031 << endl;
					}else{
					    TT_proposal.set_size(1,1); TT_proposal.zeros();
					    //TT_proposal = zeros<mat>(1,1);
					    //cout << 1032 << endl;
					}
					//cout << 1248 <<endl;
					qq_proposal = qq_proposal + qqLD_proposal;
					nonzero_proposal =  nonzero; 
					//if(qqLD_proposal ==1) cout << "this 1 "  <<  nonzero_proposal.t() << endl;
					nonzero_proposal.resize(qq_proposal);
					//if(qqLD_proposal ==1) cout << "this  2 "  <<  nonzero_proposal.t() << endl;
					if(qqLD_proposal > 0){
					    //cout << 1033 << " " << qq_proposal << " " << qqLD_proposal << " " << nonzeroLD_proposal << endl;
					    nonzero_proposal(span(qq_proposal - qqLD_proposal , qq_proposal -1 )) = nonzeroLD_proposal;
					    //if(qqLD_proposal ==1) cout << "this 3 "  <<  nonzero_proposal.t() << endl;
					    //cout << 1034 << endl;
					    for (int rot = 0; rot < qqLD_proposal; rot++){
						xTx_snp_gene = xTx_gene.col(nonzeroLD_proposal(rot));
						//cout << 1055 << endl;
						//cout << qq - qqLD + rot <<" nz prop "  << nonzero_proposal.size() << endl;  
						//cout <<  xTx_snp_gene.size() <<" nz prop "  << trans(nonzero_proposal( span(0, qq - qqLD +rot) )) << endl;  
						if((qq - qqLD + rot) > 0){
						    TT_proposal = cholAddCol(TT_proposal, xTx_snp_gene(nonzero_proposal(span(0, qq - qqLD + rot))));
						}else{
						    TT_proposal.set_size(1,1); 
						    tempd = sqrt(xTx_snp_gene(nonzero_proposal(0)));
						    TT_proposal(0,0) = (tempd > 0) ? tempd : (-1*tempd);
						}
						//cout << 1056 << endl;
					    }
					}
				    }
				    //find gg_proposal
				    //cout << nonzero_proposal.t() << "nz " << maskCurr.size() << endl;
				    if(qq_proposal > 0) {
					WW =  solve(trimatl(TT_proposal.t()), xTy(maskCurr(nonzero_proposal)));
					//cout << 1257 << endl;
					R2 = sum(square(WW))/y2(gene); 
					gg_proposal = (n - qq_proposal - 1.0)/2.0 * log(1 + ccCurr) - (n -1.0 )/2.0 * log(1.0+ccCurr*(1.0-R2));
				    }else{
					gg_proposal = 0;
				    }
				    //gg_proposal = -qq_proposal/2.0 * log(1 + ccCurr) - (n + nu(gene))/2.0 * log (nu(gene) * lambdal(gene) + y2(gene)  -  sum(square(WW))/(1 + 1/ccCurr) );
				    //cout << 1038 << endl;
				    //cout << "bhow bhow!" <<endl;
				    //cout << zeroLD_proposal.t() <<  " nonzeroLD" << endl;
				    //cout << nonzeroLD.t() << "zeroLD " << endl;
				    //cout << zeroLD.t() << endl;
				    //cout << nonzeroLD_proposal.t() << " zeroLD_proposal " << endl; 
				    ////find acceptance probability
				    ////arma prod function is 1 if the matrix is empty and even prod(1-empty.matrix) is 1.        
				    prior_proposal = prod( pa1x_vec(nonzeroLD_proposal)) * prod(1 - pa1x_vec(zeroLD_proposal));
				    //cout << 1039 << endl;
				    trans_proposal = prod(proposal_gamMaCurr(nonzeroLD_proposalInx)) * prod(1 - proposal_gamMaCurr(zeroLD_proposalInx)); 
				    //cout << 1040 << endl;
				    prior_curr = prod( pa1x_vec(nonzeroLD)) * prod(1 - pa1x_vec(zeroLD));
				    //cout << 1041 << endl;
				    trans_curr = prod( proposal_gamMaCurr(nonzeroLDInx)) * prod(1 - proposal_gamMaCurr(zeroLDInx)); 
				    //cout << 1174 << endl;
				    accept = exp( gg_proposal - gg_curr) * (prior_proposal/prior_curr) * (trans_curr/trans_proposal) * (prior_changeNotLD /trans_changeNotLD);
				    //cout << gene << " snps "<< snps_inx.t() << " accept " <<  accept <<  " gg_proposal "<< gg_proposal << " gg_curr  "<< gg_curr <<  " prior_proposal "<< prior_proposal *prior_changeNotLD << " prior_curr  "<< prior_curr  << " trans_proposal " <<   trans_proposal *trans_changeNotLD <<  " trans_curr "<< trans_curr<<  endl;
				    //cout << nonzero.t() << "nz " << nonzero_proposal.t() << endl;
				    randDraw = randu(1);
				    if(verbose){
					//bfCurr(nonzeroLD) = zeros<vec>(qqLD);
					//bfCurr(nonzeroLD_proposal) = ones<vec>(qqLD_proposal);
					gamMa_rate(snps_inx) = zeros<vec>(snps_inx.size());
				    }

				    if( accept > randDraw(0) ){
					gamMaCurr(nonzeroLD) = zeros<vec>(qqLD);
					if(cholEstimate ==1) gamMaCurr(zeroNotLD_proposal) = zeros<vec>(zeroNotLD_proposal.size());
					gamMaCurr(nonzeroLD_proposal) = ones<vec>(qqLD_proposal);
					nonzero = nonzero_proposal; 
					gg_curr = gg_proposal;
					qq = qq_proposal;
					TT_curr = TT_proposal;
					if(verbose){	
					    gamMa_rate(snps_inx) = ones<vec>(snps_inx.size());
					}
				    }
				    //cout << nonzero.t() << "nz " << maskCurr.size() << endl;
				}
				//cout <<  ldinx  << " ld inx " << gene << "gene " << endl;
			    }
			    //find gamMa_prob
			    vec pax1_vec =  gamMaCurr* log(a11/a10)  + (1-gamMaCurr ) * log( (1 - a11)/(1 - a10));  
			    //cout << 1040 << endl;
			    gamMa_probCurr = (iter == 0)? gamMaCurr : (- theta_bias + pax1_vec + aTfCurr);
			    //cout << 1040 << endl;
			    //gamMa_rate = gamMaCurr;
			}
			   gamMa_inclusion = gamMa_inclusion + gamMaCurr;
		    }
		    //cout << 932<<  "\t" << gene << endl;
		    if(sample_logit){
			if(oversample){
			    double sum_inclusion = sum(gamMa_inclusion);
			    if((localmaxiter*p_curr  > sum_inclusion) && (sum_inclusion > 0)){ 
				sample_prob = (localmaxiter*p_curr - sum_inclusion) * gamMa_inclusion +  sum_inclusion *(localmaxiter- gamMa_inclusion);
			    } else sample_prob.ones(); 
			    //cout << localmaxiter*p_curr<< "sample_prob " << sum_inclusion << endl;
			    //uvec ret = Rcpp::RcppArmadillo::sample(bootstrapS, bootstrap, true, sample_prob);
			    //ret = sample(p_curr,Named("size",bootstrap),Named("prob",sample_prob), Named("replace",true));
			    //sample.prob.ones(); // removed the oversampling and balancing of enhancers.
			    ransampl_set( ws, sample_prob.memptr() );
			    a1.randu(bootstrap); a2.randu(bootstrap);
			    for (int ii = 0; ii < bootstrap; ii++) {
				ret = ransampl_draw( ws, a1(ii), a2(ii));
				logit_freqCurr(ret) = 	logit_freqCurr(ret)+1; 
			    }

			}else if(!balance){
			    //logit_freqCurr = ones(p_curr);
			    sort_inx = arma::sort_index(gamMa_inclusion);
			    if(p_curr >= 50) {
				topp = 5; bottomm = 45; 
			    }else{
				//topp = std::max(1.0, .10 * (double) p_curr);
				//bottomm = std::max(1.0, .50 * (double) p_curr);
				topp = p_curr; bottomm=0;
			    }
			}
			//logit_freqCurr(sort_inx(span(p_curr - topp, p_curr-1) )) = ones<vec>(topp); 
			//logit_freqCurr(sort_inx(span(0, bottomm -1) )) = ones<vec>(bottomm); 
		    }
		    zero = find(gamMaCurr==0);
		    //VV =  1/(1/sigma2 + 1/cc)(X'X)^-1,  mm = VV X'Y    
		    //cout << 962<<  "\t" << gene << endl;
		    if(qq > 0) {
			TT_inv = solve(trimatu(TT_curr), eye(qq,qq)); 
			VV = ccCurr*sigma2Curr/(1+ccCurr) * (TT_inv * TT_inv.t() );
			//cout << "here" << gene << 19263 << endl;
			mm = VV * xTy(maskCurr(nonzero))/sigma2Curr;
			//cout << nonzero.t() << endl;
			//cout << "here beTa" << gene << " " << VV(0,0)  << endl;
			beTaCurr.zeros();
			beTaCurr(nonzero) = mvrnormArma(1, mm, VV);
			//cout << "here" << gene << 19267 << endl;
			xbeTa  = x.cols(maskCurr(nonzero)) * beTaCurr(nonzero);
			xbeTaM =  x.cols(maskCurr(nonzero)) * mm;
		    }else{
			beTaCurr.zeros(); xbeTa.zeros(); xbeTaM.zeros();

		    }
		    //cout << 979<<  "\t" << gene << endl;
		    if(qq > 1){
			mat TTtemp = chol(xTx_gene(nonzero, nonzero));
			if(sum(sum(TTtemp - TT_curr)) > .001){
			    cout << sum(sum(TTtemp))<< "\t" << sum(sum(TT_curr))  << endl;
			    cout << "mismatch TT \t " << nonzero.t() << endl;
			    cout << TT_curr << endl;
			    cout << TTtemp << endl;
			} 
		    }


		    residue = y.col(gene) - xbeTa; 
		    double ev = var(y.col(gene) - xbeTaM)/varY(gene); 
		    double weightCurr = 10* pow(2, -10*(1-ev)) + 1;
		    double 	  residue2sum = sum(square(residue)); double xbeta2sum = sum(square(xbeTa));
		    //cout << 996<<  "\t" << gene << endl;
		    //cout << "gamMa" << sum(gamMaCurr) <<endl;
		    //cout << ((double) accept_count)/((double) proposal_count) << endl; 

		    //#pragma omp critical(dataupdate)
		    {
			//sigma2Curr = rinvgammaC((qq + nu(gene)+ n )/2.0, 1/2.0 * (residue2sum + nu(gene) * lambdal(gene)));     
			//ccCurr = rinvgammaC((qq+1)/2.0, (xbeta2sum + n) /2.0);
			sigma2Curr = rinvgammaC(( nu(gene)+ n + qq)/2.0, 1/2.0 * (residue2sum + xbeta2sum/ccCurr +  nu(gene) * lambdal(gene)));     
			ccCurr = rinvgammaC((qq+1)/2.0, (xbeta2sum/sigma2Curr + n) /2.0);
			sigma2(gene) = sigma2Curr;
			cc(gene) = ccCurr;
			//cout << 1007<<  "\t" << gene << "\t"<<zz.size() <<"\t"<< gene_end << endl;
			evAll(gene) =1-ev;
			if(evAll(gene) < 0) { 
			    cout << ev << " ev neg for gene" <<  gene << endl;
			    cout <<   beTaCurr(nonzero) << endl;
			    cout <<  mm  << endl;
			    cout << nonzero.t() << endl;
			    cout << TT_curr << endl;
			}

			zz(span(gene_start,gene_end)) = gamMa_inclusion/localmaxiter; //gamMa_rate - 0.5; //# TODO Instead of gamMa, gamMa_rate can be used 
			//cout << 1017<<  "\t" << gene << endl;
			//zz(span(gene_start,gene_end)) = logistic(gamMa_probCurr); //gamMa_rate - 0.5; //# TODO Instead of gamMa, gamMa_rate can be used 
			if(sample_logit){
			    if( oversample){
				logit_freq(span(gene_start,gene_end)) = logit_freqCurr; //TODO change behaviour for overlapping snps 
			    }else if (!balance) {
				logit_freqtop(gene_start + sort_inx(span(p_curr - topp, p_curr-1))) = weightCurr* ones<vec>(topp); 
				logit_freq(gene_start + sort_inx(span(p_curr - topp, p_curr-1))) = weightCurr* ones<vec>(topp); 
				if(bottomm > 0) {
				    logit_freqbottom(gene_start + sort_inx(span(0, bottomm -1))) = weightCurr*ones<vec>(bottomm);
				    logit_freq(gene_start + sort_inx(span(0, bottomm -1))) = weightCurr*ones<vec>(bottomm);
				}
			    }
			}else{
			    //cout << p_curr<< "here " << gene_end - gene_start +1<<endl;
			    logit_freq(span(gene_start, gene_end)) = ones<vec>(gene_end - gene_start + 1);
			}
			//TODO instead of condition apply pre-processor conditional.

			gamMa(span(gene_start, gene_end)) = gamMaCurr; 
			gamMa_prob(span(gene_start, gene_end)) = gamMa_probCurr; 
			bf(span(gene_start, gene_end)) = bfCurr; 
			gamMa_rateAll(span(gene_start, gene_end)) = gamMa_rate; 
			beTa(span(gene_start, gene_end)) = beTaCurr; 
			//gg(gene) = gg_curr; // not required only when qq is zero it is not recalcualted 
			TT(gene) = TT_curr;
			if((((iter -burnIn +1) %thin)==0)  && (iter >= burnIn)){
			    gamMa_rate_mat(span(gene_start, gene_end),count)= gamMa_rate;
			    bf_mat(span(gene_start, gene_end),count)= bfCurr;
			    gamMa_mat(span(gene_start, gene_end),count)= gamMaCurr;
			    gamMa_prob_mat(span(gene_start, gene_end),count)= gamMa_probCurr;
			    beTa_mat( span(gene_start, gene_end), count) = beTaCurr;
			    cc_mat(gene, count) = cc(gene);
			    ev_mat(gene, count) = evAll(gene);
			    sigma2_mat(gene, count) = sigma2(gene);
			}
		    }

		    ransampl_free( ws );
		}
		//cout << gene << "gene done" << endl;
	    }
	    beTaAcc = beTaAcc + beTa;
	    //gamMaAcc = gamMaAcc + gamMa;
	    gamMa_rateAcc.col(iter%accIter) = gamMa;
	}
	//estimate_alpha = false;
	if(verbose) 
	{ 
	    t2 = std::clock();
	    cout << "start of alpha cycle \t " << iter << "\t" << ((float)t2 - (float) t1) / CLOCKS_PER_SEC << "s"<<   endl;
	    cout << "number of eeSNP \t" << sum(gamMa) << endl;
	}

	if(estimate_alpha){
	    //sampling of snp of genes.
	    if(sample_logit){
		uvec nonzero_logittop = find(logit_freqtop > 0 );
		uvec nonzero_logitbottom = find(logit_freqbottom > 0 );
		bool canSample = ((nonzero_logittop.size()  > (num_logit_postrain)) &&  (nonzero_logitbottom.size()  > (4*num_logit_postrain)))? true :  false;
		if((!oversample) && canSample){
		    logit_freq.zeros();
		    ret1 = sample(p, Named("size", num_logit_postrain), Named("prob", logit_freqtop), Named("replace", false));
		    for (int tt = 0; tt < num_logit_postrain; tt++) 
		    logit_freq(ret1[tt] -1) = 1;
		    ret2 = sample(p, Named("size", 4*num_logit_postrain), Named("prob", logit_freqbottom), Named("replace", false));
		    for (int tt = 0; tt < 4*num_logit_postrain; tt++) 
		    logit_freq(ret2[tt] -1) = 1;
		}	
	    }
	    gamMaABS = mean(gamMa_rateAcc,1);
	    uvec nonzero_logit;

	    //adjusting gamMa_prob
	    bool use_logisitic = true;
	    double bias_adjusted = 0;
	    if (use_logisitic) {
		if ((regulator_prior > 0) && (iter > 0) ){
		    uvec temp1= find(gamMa_prob > 0.0);
		    int ynum =  (temp1).size();
		    cout<< "number of regulatory enhancers \t " << ynum << endl;
		    vec gamMa_sort = sort(gamMa_prob);
		    int reg_inx = regulator_prior * p;
		    bias_adjusted = gamMa_sort(p - reg_inx);
		    //alpHa(0) = alpHa(0) - bias_adjusted;
		    gamMa_prob = gamMa_prob - bias_adjusted;
		    temp1=find(gamMa_prob > 0);
		    cout<< "number of regulatory enhancers after adjustment \t " << temp1.size() << " \t bias adjusted =" << bias_adjusted << endl;
		}
		    gamMa_prob = logistic(gamMa_prob ); 
	    }

	    int logitIter = 2;
	    for (int logitCnt = 0; logitCnt < logitIter; logitCnt++) {
		if(balance){
		    //cout<< "here" << gamMa_rateAcc.n_cols<<endl;
		    if(usebeTaprob) beTaABS = mean(beTaAcc,1);
		    //cout<< "here1" <<gamMaABS.size() << endl;
		    //if((iter+1) >= accIter){
			//if (((iter+1) % accIter)==0) { 
			    //beTaABS = ABS(beTaAcc)/accIter;
			    //gamMaABS = sum(gamMa_rateAcc);
			//}
		    //} else {
			//beTaABS = ABS(beTa);
			//gamMaABS = gamMa;
		    //}
		    colvec nonzero_prob;
		    if (usebeTaprob) {
			nonzero_logitInx = ((gamMa > gamMa_thres) % (beTaABS > beTa_thres)); 
			nonzero_logit = find(nonzero_logitInx);
		    }else nonzero_logit = find(gamMa > gamMa_thres);
		    nonzero_prob = gamMaABS(nonzero_logit);
		    uvec zero_logit = find(gamMa < 0.1);
		    colvec zero_prob = ones(zero_logit.size());
		    zero_prob = zero_prob - gamMaABS(zero_logit);
		    int postrain_num = (num_logit_train > -1) ? min( nonzero_logit.size()*ratio , num_logit_train/2):  nonzero_logit.size();
		    //int negtrain_num = (num_logit_train > -1) ? min( zero_logit.size()*ratio , num_logit_train/2):  zero_logit.size();
		    //int train_num = min(postrain_num, negtrain_num);
		    int train_num = postrain_num;
		    //negtrain_num = train_num*2;
		    //negtrain_num = 10000;
		    //train_num = min(num_logit_train/2, train_num);
		    //num_logit_postrain = min( num_logit_train/2,  nonzero_logit * ratio) ;
		    //int num_logit_negtrain = min( num_logit_train/2,  (p - nonzero_logit)*ratio) ;
		    logit_freq.zeros();
		    geneEpi_freq.zeros();
		    cout <<"train_num \t"<< train_num << "\t" << negtrain_num << "\t" << train_num << "\t"<< gamMa_thres    << endl;
		    if(train_num > (f/2 +1) ){
			ret1 = sample(nonzero_logit.size(), Named("size", train_num), Named("prob", nonzero_prob),Named("replace", false));
			for (int tt = 0; tt < ret1.size(); tt++) 
			logit_freq(nonzero_logit(ret1[tt] -1)) = logit_freq(nonzero_logit(ret1[tt] -1)) + 1;
		    }else{
			cout << "oops very few confident regulators to train" << endl;
			ret1 = sample(nonzero_logit.size(), Named("size", f/2.0 + 1), Named("prob", nonzero_prob),Named("replace", true));
			for (int tt = 0; tt < ret1.size(); tt++) 
			logit_freq(nonzero_logit(ret1[tt] -1)) = logit_freq(nonzero_logit(ret1[tt] -1)) + 1;
			//ret1 = sample(p, Named("size", f/2 + 1), Named("prob", nonzero_prob), Named("replace", true));
			//for (int tt = 0; tt < ret1.size(); tt++) 
			//logit_freq(ret1[tt] -1) = logit_freq(ret1[tt] -1) + 1;
		    }
		    ret2 = sample(p , Named("size", negtrain_num),  Named("replace", false));
		    for (int tt = 0; tt < ret2.size(); tt++) 
		    logit_freq(ret2[tt] -1) = logit_freq(ret2[tt] -1)+ 1;
		}


		if(verbose){ 
		    cout << " end of balance \t " << iter << "\t" << ((float)(std::clock())  - (float) t2) / CLOCKS_PER_SEC << "s"<<   endl;
		    t2 = std::clock();
		}
		//zz = logistic(gamMa_prob);
		//zz = gamMaABS;
		uvec freq_nonzero = find(logit_freq > 0); p_nonzero = freq_nonzero.size();
		//cout << " p_nonzero " <<  p_nonzero << endl;
		mat feature_train1 = featureCat.rows(freq_nonzero);
		if (use_logisitic) {

		    zz = gamMa_prob;
		    zz = zz - 0.5;
		    zz(freq_nonzero) = logit_freq(freq_nonzero) % zz(freq_nonzero);

		    // variable selection of alpha
		    uvec alpha_select;
		    int f_curr=2;
		    if(logistic_variable_selection){
			mat feature_train = feature_train1.cols(span(1,(f-1))); // ones are not added anymore 
			logistic_variable_selection_mcmc( feature_train, zz(freq_nonzero),  logistic_select, cc_epi, sigma2_epi, lambdal_epi, nu_epi);
			uvec alpha_select = find(logistic_select > 0);
			int f_curr = alpha_select.size() + 1;
			//else logistic_select.ones(); 
		    }

		    //#sample wi
		    if(f_curr > 1){
			uvec alpha_select1 = ones<uvec>(f_curr);
			if(logistic_variable_selection) { 
			    alpha_select1(0) = 0; alpha_select1(span(1, f_curr-1)) =  alpha_select + 1; 
			    feature_nonzero =  feature_train1.cols(alpha_select1); // first column (bias) is always part of the model
			    //cout << "feature_nonzero sd " << trans(stddev(feature_nonzero)) << endl;  
			    logit_freq_nonzero(span(0,p_nonzero-1)) = logit_freq(freq_nonzero);
			    fTalpHa = feature_nonzero * alpHa(alpha_select1);
			    if(verbose){ 
				cout << " start of rpg \t " << iter << "\t" << ((float)(std::clock())  - (float) t2) / CLOCKS_PER_SEC << "s"<<   endl;
				t2 = std::clock();
			    }
			    rpg_hybrid( &w(0), &logit_freq_nonzero(0), &fTalpHa(0), &p_nonzero);
			    if(verbose){ 
				cout << " end of rpg \t " << iter << "\t" << ((float)(std::clock())  - (float) t2) / CLOCKS_PER_SEC << "s"<<   endl;
				t2 = std::clock();
			    }
			    //calculate m_w V_w;
			    mat diag_w  = diagmat(w(span(0, p_nonzero-1 )));
			    mat V_w = ((trans(feature_nonzero) * diag_w * feature_nonzero)+  diagmat( B_inv_alpHa(alpha_select1) )).i();
			    vec m_w = V_w * ( (trans(feature_nonzero) * zz(freq_nonzero)) +   B_inv_m_w_prior(alpha_select1));
			    //cout << "m_w " << m_w.t() << endl;
			    alpHa.zeros();
			    alpHa(alpha_select1)  = trans(mvrnormArma(1, m_w, V_w));

			}else{
			    //uvec freq_nonzero = find(logit_freq > 0); p_nonzero = freq_nonzero.size();
			    feature_nonzero =  feature_train1; // first column (bias) is always part of the model
			    //cout << "feature_nonzero sd " << trans(stddev(feature_nonzero)) << endl;  
			    logit_freq_nonzero(span(0,p_nonzero-1)) = logit_freq(freq_nonzero);
			    fTalpHa = feature_nonzero * alpHa;
			    if(verbose){ 
				cout << " start of rpg \t " << iter << "\t" << ((float)(std::clock())  - (float) t2) / CLOCKS_PER_SEC << "s"<<   endl;
				t2 = std::clock();
			    }
			    rpg_hybrid( &w(0), &logit_freq_nonzero(0), &fTalpHa(0), &p_nonzero);
			    if(verbose){ 
				cout << " end of rpg \t " << iter << "\t" << ((float)(std::clock())  - (float) t2) / CLOCKS_PER_SEC << "s"<<   endl;
				t2 = std::clock();
			    }
			    //calculate m_w V_w;
			    mat diag_w  = diagmat(w(span(0, p_nonzero-1 )));
			    mat V_w = ((trans(feature_nonzero) * diag_w * feature_nonzero)+  diagmat( B_inv_alpHa )).i();
			    vec m_w = V_w * ( (trans(feature_nonzero) * zz(freq_nonzero)) +   B_inv_m_w_prior);
			    //cout << "m_w " << m_w.t() << endl;
			    alpHa.zeros();
			    alpHa  = trans(mvrnormArma(1, m_w, V_w));
			}
			/*if( logitCnt ==0) {
			    aTf = featureCat * alpHa;
			    double aTfthres = mean(aTf); aTfthres = max(aTfthres, 0);
			    uvec sortInx = sort_index(aTf);
			    logit_freq.zeros();
			    int train_num = min(num_logit_train/2, p/10);
			    uvec tempInx = sortInx(span(0, train_num-1));
			    uvec labelInx = tempInx( find(zz(tempInx) < 0)); 
			    logit_freq(labelInx) = ones(labelInx.size());
			    cout<<"negtrain size \t " << labelInx.size() ;
			    tempInx = sortInx(sortInx(span(p - train_num, p-1)));
			    labelInx = tempInx( find(zz(tempInx) > 0)); 
			    cout<<"\t postrain size \t " << labelInx.size() << endl;
			    logit_freq(labelInx) = ones(labelInx.size());
			    labelInx = nonzero_logit( find(aTf(nonzero_logit) > aTfthres));
			    cout<<"\t postrain size  \t " << labelInx.size() << endl;
			    logit_freq(labelInx) = ones(labelInx.size());
			    freq_nonzero = find(logit_freq > 0); p_nonzero = freq_nonzero.size(); 
			    zz(freq_nonzero) = logit_freq(freq_nonzero) % zz(freq_nonzero);
			    feature_train1 = featureCat.rows(freq_nonzero);
			}*/
			//gamMa_prob = logistic(featureCat * alpHa)
		    }else{
			alpHa.zeros();
		    }
		}else{ // theta regression
		    // input gamMa_probs, featureCat
		    // output alpHa
		    vec alpHa_nobias = alpHa(span(1,f-1));
		    //	zz = gamMa_prob(freq_nonzero);
		    //cout << "theta_mean \t" << endl;
		    //cout<< "gamMa_prob\t " << gamMa_prob(span(1,10)) << endl; 
		    vec y_l = logit(gamMa_prob(freq_nonzero));
		    double aTfbias = mean(y_l);
		    //cout << "y"<< (y_l(span(1,10))).t() << endl;
		    y_l = y_l - ones(p_nonzero) * aTfbias;
		    //cout << "y"<< (y_l(span(1,10))).t() << endl;
		    mat xx = feature_train1.cols(span(1,(f-1))); // ones are not added anymore
		    //cout << "xx"<< (xx.col(0)).t() << endl;
		    mat V_w =  cc_l*sigma2_l/(1+cc_l) * ((xx.t()* xx + small*eye(f-1, f-1)).i());
		    //cout << "V_w"<< V_w.row(0) << endl;
		    vec m_w = V_w * xx.t() * y_l/sigma2_l;
		    //cout << "m_w"<< m_w.t() << endl;
		    alpHa_nobias = trans(mvrnormArma(1, m_w, V_w));
		    vec xalpHa  = xx * alpHa_nobias;
		    vec residue = y_l - xalpHa; 
		    double 	  residue2sum = sum(square(residue)); double xalpHa2sum = sum(square(xalpHa));
		    sigma2_l = rinvgammaC(( nu_l+ f - 1 + p_nonzero)/2.0, 1/2.0 * (residue2sum + xalpHa2sum/cc_l +  nu_l * lambdal_l));     
		    cc_l = rinvgammaC(f/2.0, (xalpHa2sum/sigma2_l + p_nonzero) /2.0);
		    sigma_l = sqrt(sigma2_l);
		    alpHa(span(1,f-1)) = alpHa_nobias;
		    alpHa(0) = aTfbias;
		    aTf = featureCat * alpHa;
		}

		if(verbose){ 
		    cout << " end of alpHa estimate \t " << iter << "\t" << ((float)(std::clock())  - (float) t2) / CLOCKS_PER_SEC << "s"<<   endl;
		    t2 = std::clock();
		}

	    }

	    //genetic epigenetic weight
	    /* cout << 1011 << endl;
	    freq_nonzero = find(geneEpi_freq > 0); p_nonzero = freq_nonzero.size(); 
	    cout << 1012 << endl;
	    zz(freq_nonzero) = geneEpi_freq(freq_nonzero) % zz(freq_nonzero);
	    cout << 1014) << endl;
	    feature_train1 = featureCat.rows(freq_nonzero);
	    cout << 1014 << endl;
	    mat geneEpi(p_nonzero, 4);
	    cout << 1015 << endl;
	    if(f_curr > 1){
		feature_nonzero =  feature_train1.cols(alpha_select + 1); // alpHa select don't contain bias. First column  not part of the model   
		cout << 1016 << endl;
		logit_freq_nonzero(span(0,p_nonzero-1)) = geneEpi_freq(freq_nonzero);
		cout << 1017 << endl;
		fTalpHa = feature_nonzero * alpHa(alpha_select + 1 );
		cout << 1018 << endl;
		geneEpi.col(0) = ones(p_nonzero,1);  geneEpi.col(1) =  fTalpHa;geneEpi.col(2) =  bf(span(1,p_nonzero)); geneEpi.col(3) =  bf(span(1,p_nonzero))%fTalpHa;
		cout << 1019 << endl;
		geneEpiTtau = geneEpi * tau;
		cout << 1010 << endl;
		rpg_hybrid( &w(0), &logit_freq_nonzero(0), &geneEpiTtau(0), &p_nonzero);
		cout << 1011 << endl;
		diag_w  = diagmat(w(span(0, p_nonzero-1 )));
		cout << 1012 << endl;
		V_w = ((trans(geneEpi) * diag_w * geneEpi)+  diagmat( B_inv_alpHa(span(1,4)) )).i();
		cout << 1013 << endl;
		m_w = V_w * (trans(geneEpi) * zz(freq_nonzero));
		cout << 1014 << endl;
		tau.zeros();
		cout << 1015 << endl;
		tau  = trans(mvrnormArma(1, m_w, V_w));
		cout << 1016 << endl;
	    }
	    mat tempGeneEpi(p, 4);
	    alpHa(0) = 0;
	    fTalpHa = featureCat * alpHa;
	    tempGeneEpi.col(0) = ones(p,1);  tempGeneEpi.col(1) =   fTalpHa;tempGeneEpi.col(2) =  bf; tempGeneEpi.col(3) =  bf%fTalpHa;
	    geneEpiTtau = tempGeneEpi * tau;
	    freq_nonzero = find(geneEpiTtau > 0);
	    cout << tau.t() << endl;
	    cout << freq_nonzero.size() << endl; 	    
	    */


	    //weight readjustment
	    //if(regulator_prior > 0){
		//aTf = featureCat * alpHa;
		//uvec temp1= find((logistic(aTf)) > 0.5);
		//int ynum =  (temp1).size();
		//cout<< "number of enhancers \t " << ynum << endl;
		//double ym = ((double) ynum)/p; // King et. al eqn. 7.
		//alpHaAdj = alpHa;
		//alpHaAdj(0) =  alpHa(0) - log( ym*(1-regulator_prior)/(regulator_prior*(1 -ym))); 
		//double aTfthres = mean(aTf); 
		/* {vec aTf_sort = sort(aTf);
		    int reg_inx = regulator_prior * p;
		    //theta_bias =  aTf_sort(p - reg_inx); 
		    //theta_bias = 0.0; //
		    //alpHa(0) = alpHa(0) - log( ym*(1-regulator_prior)/(regulator_prior*(1 -ym)));
		    //vec alpHa1 = alpHa;
		    alpHa(0) = alpHa(0) - aTf_sort(p - reg_inx); 
		    temp1=find((logistic(featureCat * alpHa)) > 0.5);
		cout<< "number of enhancers after adjustment \t " << temp1.size() << endl; }*/
		//alpHa(0) = alpHaAdj(0);
	    //}
	    //uvec temp1=find((logistic(featureCat * alpHa)) > 0.5);
	    //cout<< "number of enhancers \t " << temp1.size() << endl;
	    alpHaAdj = alpHa;
	    aTf = featureCat * alpHaAdj;
		if ((regulator_prior > 0)  ){
		    uvec temp1= find(aTf > 0.0);
		    int ynum =  (temp1).size();
		    cout<< "number of regulatory enhancers \t " << ynum << endl;
		    vec gamMa_sort = sort(aTf);
		    int reg_inx = regulator_prior * p;
		    bias_adjusted = gamMa_sort(p - reg_inx);
		    alpHaAdj(0) = alpHa(0) - bias_adjusted;
		    aTf = aTf - bias_adjusted;
		    temp1=find(aTf > 0);
		    cout<< "number of regulatory enhancers after adjustment \t " << temp1.size() << " \t bias adjusted =" << bias_adjusted << endl;
		}
		    gamMa_prob = logistic(gamMa_prob ); 
	    cout<<"alpHa \t"<< alpHa.t() <<endl;

	}
	if((((iter - burnIn + 1) %thin)==0)  && (iter >= burnIn)){
	    alpHa_mat.col(count) = alpHa;
	    //tau_mat.col(count) = tau;
	    //logistic_select_mat.col(count) = logistic_select;
	    count = count +1;
	}

    }
    out["alpHa.mat"] = trans(alpHa_mat); out["beTa.mat"] = trans(beTa_mat); out["gamMa.mat"] = trans(gamMa_mat); out["gamMa_rate.mat"] = trans(gamMa_rate_mat); out["gamMa.prob.mat"] = trans(gamMa_prob_mat); out["cc.mat"] = trans(cc_mat);out["sigma2.mat"] = trans(sigma2_mat); out["ev"] = trans(ev_mat);
    out["bf.mat"] = trans(bf_mat); out["bf.single"] = bf_single;	
    if(logistic_variable_selection)  out["logistic_select"] = logistic_select_mat.t(); 
    if(verbose) 
    { 
	t1 = std::clock();
	cout << "end of program \t " << "\t" << ((float)t1 - startOfepi_eQTL) / CLOCKS_PER_SEC << "s"<< endl;
    }

    return out;
}








using namespace arma;		// shorthand
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
Rcpp::List eQTL_empirical_prior(arma::mat x, arma::mat y, arma::mat feature, arma::mat pairFeature, arma::umat mask2, arma::colvec alpHa,   arma::colvec gamMa, bool estimate_alpha, bool estimate_beta, arma::colvec B_inv_alpHa, int itermax, int  thin, int burnIn,  int threads, bool verbose=true, bool balance = true, bool logistic_variable_selection = true, bool oversample = false,bool  sample_logit = false, bool use_raoblackwell =true, int num_logit_train = 1000, double ratio = .5, double gamMa_thres = .9, double beTa_thres=.1, double regulator_prior = 1.0/300.0, int accIter = 1, double rho = 1, double prior = .1)
{
    //regulator_prior = 0 : implies no re-estimation of alpha. num_logit_train = -1: no threshold is applied.
    #ifdef _OPENMP
    if ( threads > 0 )
    omp_set_num_threads( threads );
    REprintf("Number of threads=%i\n", omp_get_max_threads());
    #endif
    double eQTL_threshold = 1; // 1- implies eQTL_threshold is ignored 
    double regulator_threshold = 0; // 0- implies regulator_threshold is ignored 
    bool dual_filter = false;
    bool usebeTaprob = false, Berger=true;
    std::clock_t t1, t2, t3; 
    t2 = std::clock();
    float startOfepi_eQTL = float(t2);
    int localmaxiter =20; int bootstrap=50/localmaxiter; 
    int num_out_samples = ceil(((double)(itermax - burnIn))/((double) thin));
    int p=mask2.n_rows,  n=y.n_rows, numGenes=y.n_cols; //, numSamp_logit = n/2, ;
    //cout << 454 << endl;
    mat feature1 = join_rows( ones(p,1), feature);
    //cout << 456 << endl;
    //mat feature1 = feature; // ones is added in logistic part
    int pairf = pairFeature.n_cols; int f1 = feature1.n_cols ; int f = f1 + pairf;
    int num_logit_postrain = (numGenes <  40)?  5*numGenes: 1000;
    mat beTa_mat(p, num_out_samples), alpHa_mat(f, num_out_samples), tau_mat(f, num_out_samples), gamMa_mat(p, num_out_samples), gamMa_prob_mat(p, num_out_samples), bf_mat(p, num_out_samples),  gamMa_rate_mat(p, num_out_samples),  gamMa_rateAcc(p, accIter), logistic_select_mat(f-1, num_out_samples), 
    ev_mat(numGenes, num_out_samples), cc_mat(numGenes, num_out_samples), sigma2_mat(numGenes, num_out_samples);
    colvec w(p), beTa(p), beTaAcc(p), beTaABS(p), gamMaAcc(p), gamMaABS(p), gamMa_rateAll(p);
    colvec alpHaAdj= alpHa;
    uvec nonzero_logitInx(p); 
    gamMa_rateAcc.fill(.5);
    vec gamMa_prob, theta, aTf;
    vec ptheta(p);
    //y = y - (ones(n) * mean(y)) ;
    center_norm(y);
    center_norm(x);
    //x = x - (ones(n) * mean(x)) ;
    colvec xTy(p); 
    field<mat> TT(numGenes); 
    field<mat> xTx(numGenes); 
    field<uvec> mask(numGenes); 
    vec gg(numGenes); gg.zeros(); 
    uvec gene1_start(numGenes) , gene1_end(numGenes); gene1_start.zeros(); gene1_end.zeros();
    vec cc(numGenes); cc.fill(10.0); double cc_epi =10; // g-priors
    vec evAll(numGenes); evAll.zeros();
    int count  = 0;
    colvec logit_freq, geneEpi_freq, logit_freqtop, logit_freqbottom,  logit_freq_nonzero,  zz, bf(p), tau(4);
    tau << 0 << 1 << 1<< 0;



    Rcpp::NumericVector	 ret1(num_logit_postrain), ret2(4*num_logit_postrain);
    logit_freq = ones(p); geneEpi_freq = ones(p); logit_freqtop = ones(p); logit_freqbottom = ones(p);  logit_freq_nonzero = ones(p);  zz=zeros(p); 
    zz = gamMa - .5; 
    vec logistic_select(f-1); 
    if( logistic_variable_selection) {
	logistic_select.zeros();
	for (int finx = 0; finx < f-1; finx++) {
	    if(alpHa(finx + 1) !=0)  logistic_select(finx) = 1;
	}
	//uvec temp = find(alpHa  !=0.0);
    }else logistic_select.ones();
    mat featureCat(p, f); 
    //cout << 493 << endl;
    if(pairf > 0) {
	//if( pairf ==1) featureCat.col(f) = pairFeature;
	//else 
	featureCat.cols(span(f1, f-1 )) = pairFeature;
    }
    //cout << 499 << endl;
    mat diag_w, V_w; vec m_w(p);
    if(beTa.size() !=p) beTa = zeros(p);

    mat feature_nonzero; 
    vec fTalpHa(p), geneEpiTtau(p);
    int p_nonzero;
    double nu1 = log( regulator_prior/ ( 1 - regulator_prior));
    // using regulatory prior 
    vec B_inv_m_w_prior= zeros(f);  // B_inv_m_w_prior(0) = B_inv_alpHa(0)  * logit(regulator_prior); // B-1 b page 6. polson et. al 	
    alpHa(0) = logit(regulator_prior);


    // epi hyper-parameter 
    double sigma2_epi = .01, lambdal_epi = .01, nu_epi = 100;  

    //# needed only  once
    RNGScope scope ;
    Environment base("package:base");
    Function sample = base["sample"];
    Rcpp::List out;
    bool display_progress = true;
    Progress progress(itermax, display_progress); 
    double small = .001;
    rowvec varY; varY =  var(y);
    rowvec sigma2,  lambdal;
    sigma2 = lambdal = 1*varY/9; // estimate of sigma2 variance of residue
    rowvec y2 = sum(square(y), 0);
    double minR = 0, maxR = 50, by=.5; 
    vec nu =  zeros(numGenes);
    double seqT, seqTmin, nuCurr,  bb;
    uvec maskCol = mask2.col(1);
    // calculation of bias of logistic regression
    mat featureMean(f,numGenes);

    // for theta regression 
    double sigma2_l =1,  sigma_l =1, lambdal_l = 0, nu_l =0, cc_l =10.0;


    // initialize

    gamMa_rateAll =  zeros(p); //this is theta
    for (int gene = 0; gene < numGenes; gene++) {
	int p_curr, qq, qq_curr, k;
	uvec temp=maskCol(find(mask2.col(0) == (gene +1) ));  
	mask(gene) = temp-1;
	p_curr = mask(gene).size();

	if(gene >0)
	gene1_start(gene) = gene1_end(gene-1) + 1; 
	gene1_end(gene) = gene1_start(gene) + p_curr - 1;
    }
    #pragma omp parallel for  schedule(dynamic)   
    for (int gene = 0; gene < numGenes; gene++) {
	seqTmin = 1; nuCurr = 0; 
	for (int ii = 0; ii < (maxR/by + 9000); ii++) {
	    if(ii < maxR/by) bb = minR + by*ii;
	    else bb = exp( -30.0 + 60.0/9000.0 *   (double) (ii - maxR/by));
	    seqT = pinvgammaC(varY(gene), bb/2.0, lambdal(gene)*bb/2.0) - .9;
	    if(( seqT< seqTmin) && (seqT > 0)){
		seqTmin = seqT; nuCurr = bb;
	    }
	}
	if(!Berger) nu(gene) = nuCurr;
	else nu(gene) =0;

	int p_curr, qq, qq_curr, k;
	//uvec temp=maskCol(find(mask2.col(0) == (gene +1) ));  
	//mask(gene) = temp-1;

	//if(gene >0)
	//gene1_start(gene) = gene1_end(gene-1) + 1; 
	//gene1_end(gene) = gene1_start(gene) + p_curr - 1;

	p_curr = mask(gene).size();
	if(p_curr > 0){
	    uvec nonzero, zero, nonzero_curr;
	    //cout << 551<<endl;
	    featureCat(span(gene1_start(gene), gene1_end(gene)), span(0, f1-1)) = feature1.rows(mask(gene));
	    //cout << 553<<endl;
	    featureMean.col(gene) = trans(mean(featureCat.rows(span(gene1_start(gene), gene1_end(gene)) ))); 
	    //cout << 555<<endl;
	    xTy(span(gene1_start(gene), gene1_end(gene))) = x.cols(mask(gene)).t() * y.col(gene);
	    //cout << 557<<endl;
	    xTx(gene) = x.cols(mask(gene)).t()* x.cols(mask(gene)) + small*eye(p_curr, p_curr);
	    //cout << 558<<endl;
	    vec gamMaCurr = gamMa(span(gene1_start(gene),gene1_end(gene)));
	    //cout << 559<<endl;
	    nonzero_curr = find(gamMaCurr > 0); 
	    int qq_curr  = nonzero_curr.size();
	    uvec  maskCurr =  mask(gene);
	    //cout << 550<<endl;
	    double ccCurr = cc(gene);
	    double gg0Curr =  -(n + nu(gene))/2.0 * log(nu(gene) * lambdal(gene) + y2(gene));
	    if(qq_curr == 0 ){
		TT(gene).set_size(1,1); TT(gene).zeros();
		//cout << 552 << "\t" << gg.size()  << "\t" << lambdal.size() << "\t" << gene << "\t" << nu.size() << "\t" << y2.size()<<endl;
		gg(gene)  = gg0Curr;
		//cout << 554<<endl;
	    }else {
		TT(gene) = chol( x.cols(maskCurr(nonzero_curr)).t() * x.cols(maskCurr(nonzero_curr)) + small * eye(qq_curr, qq_curr)); 
		//cout << 556 << nonzero_curr.t() <<endl;
		vec WW = solve(trimatl(TT(gene).t()), xTy(maskCurr(nonzero_curr))); 
		//cout << 558<<endl;
		gg(gene) = -qq_curr/2.0 * log(1 + ccCurr) - (n + nu(gene))/2.0 * log (nu(gene) * lambdal(gene) + y2(gene)  - sum(square(WW))/(1 + 1/ccCurr) );
		//cout << 561<<endl;
	    }
	}
    }
    //cout << "gg" << gg.t()<< endl;
    //cout << "nu" << nu.t() << endl;
    //cout << "gene_start " << gene1_start.t() <<endl;
    //cout << "gene_end " << gene1_end.t() <<(endl;

    double sumAlpHa = sum(alpHa);
    //double prior =  (p > (40* numGenes))?  ((20.0*numGenes)/p ): 0.05;
    //double prior =  (p > (100* numGenes))?  ((50.0*numGenes)/p ): 0.01;
    double a11, a10; // a10 = p(gamMa = 1| theta =0)
    a11 = prior;
    a10 = prior/rho;
    cout << " a11 \t"  << a11 << "\t a10 \t" << a10 << endl;  
    //double thetamean = logit(prior); double thetamax = logit(0.9); double thetamin = 2.0*thetamean - thetamax;
    //double sigma_0 = (thetamax - thetamin)/6.0;
    gamMa_prob =  zeros(p); gamMa_prob.fill(0.5);//this is theta
    //cout << "thetamean \t" << thetamean << endl;
    // switching model

    double theta_bias = 0.0;
    double theta_odds = regulator_prior/(1-regulator_prior); // sum(gamMa_prob)/(p - sum(gamMa_prob));
    for (int iter = 0; iter < itermax; iter++) {
	if(verbose) 
	{ 
	    t1 = std::clock();
	    cout << "start of beta cycle \t" << iter << "\t" << ((float)t1 - (float) t2) / CLOCKS_PER_SEC << "s"<< endl;
	}
	if (Progress::check_abort() )
	return out;
	progress.increment();
	    aTf = featureCat * alpHa;
	if(prior > 0.5){ 
		uvec temp10 = find(aTf > 5);
		if(temp10.size() > 0)
		aTf(temp10) = ones(temp10.size()) * 5;
		aTf = exp(aTf);
	}
	////cout << "running iteration " << iter << endl;
	if(estimate_beta){
	    logit_freq.zeros(); logit_freqtop.zeros(); logit_freqbottom.zeros();
	    if((iter %  accIter) ==0) {
		beTaAcc.zeros();
		gamMaAcc.zeros();
	    }
	    // find phi (std of theta)
	    //double phi = stddev(gamMa_prob - aTf);
	    //cout << "phi \t" << phi << endl;
	    ////phi = .5;
	    //vec randDraw1; randDraw1 = randu(1);
	    //dual_filter = (randDraw1(0) < .1) ? true : false; 
	    //double sigma2_eff = 1.0/ (1.0/(sigma2_l) + 1.0/square(sigma_0) ); double sigma_eff = sqrt(sigma2_eff);
	    /*if(false) 
	    {
		if(iter == 0) {
		    a11 = .8; a10 = 0.01;
		}else{
		    a11 = sum(gamMa_rateAll % gamMa_prob)/sum(gamMa_prob); 
		    a10 = sum(gamMa_rateAll % ( ones(p) - gamMa_prob) )/(p - sum(gamMa_prob)); 
		}
	    }else{
		a11 =  min(.99, .99 * numGenes * 5/sum(gamMa_prob));  
		a10 =  min(.99, .01 * numGenes * 5/(p - sum(gamMa_prob)));
	    }*/

	    #pragma omp parallel for schedule(dynamic)   
	    for (int gene = 0; gene < numGenes; gene++) {
		int  p_curr, qq, qq_curr, k;
		double gamMa_snpCurr, gamMa_probsnpCurr, tempd;
		double ccCurr = cc(gene), sigma2Curr = sigma2(gene);
		uvec  maskCurr =  mask(gene);

		//#gibbs sampling;
		p_curr = maskCurr.size();
		if(p_curr > 0){
		    int topp, bottomm;
		    uvec sort_inx;
		    //Rcpp::NumericVector
		    ransampl_ws* ws = ransampl_alloc( p_curr );
		    vec sample_prob(p_curr);
		    double ret;
		    uvec bootstrapS(p_curr);
		    vec a1, a2; 
		    //vec ret(bootstrap); 
		    bootstrapS = linspace<uvec>(0, p_curr-1, p_curr);
		    double epiDiscount = 0; 
		    uvec nonzero, zero, nonzero_curr;
		    vec  bfCurr(p_curr), gamMa_rate, randDraw;
		    double bf_snpCurr, prob_proposal, proposal_ratio, accept;
		    int accept_count = 0, proposal_count =0;
		    vec mm, WW;
		    mat TT1, TT_inv,VV, TT_curr = TT(gene); // chol( x.cols(maskCurr(nonzero_curr)).t() * x.cols(maskCurr(nonzero_curr)) + small * eye(qq_curr, qq_curr))
		    int gene_start = gene1_start(gene), gene_end = gene1_end(gene); 
		    colvec xbeTa(n), xbeTaM(n), residue(n);
		    mat xTx_gene = xTx(gene); vec xTx_snp_gene;  
		    vec  gamMaCurr = gamMa(span(gene_start, gene_end)); vec gamMa_inclusion = zeros(p_curr);
		    vec logit_freqCurr = zeros(p_curr);
		    vec gamMaChange(p_curr); gamMaChange = gamMaCurr;
		    //cout << "here" << gene << 181 << gene_start << "endl" << gene_end << endl;
		    vec  beTaCurr(p_curr), geneEpiCurr(4);
		    //vec alpHaCurr = alpHa;
		    vec featureMeanCurr = featureMean.col(gene);
		    //double expalpHaF = sum(trans(alpHa(span(1,f-1))) * featureMeanCurr(span(1,f-1)));
		    //if ((p_curr > (3*regulator_prior)) && (regulator_prior  > 0) ) alpHaCurr[0] = -log(p_curr/regulator_prior -1) - expalpHaF;
		    //vec gamMa_probCurr = logistic(featureCat.rows(span(gene_start,gene_end)) * alpHaCurr);
		    //alpHaCurr(0) = 0; 


		    double gg1, gg0, gg_curr = gg(gene);
		    // nonzero need to be stored and when updating with unsorted nonzero and 1 chol at a time 
		    k = gamMaCurr.size(); zero = find(gamMaCurr==0); nonzero = find(gamMaCurr > 0); qq  = nonzero.size();
		    //cout << 763 << endl;
		    gamMa_rate = gamMa_rateAll(span(gene_start, gene_end));
		    //cout << 764 << endl;
		    nonzero_curr = nonzero;
		    // sample gamMa_prob 
		    vec gamMa_probCurr(p_curr);
		    gamMa_probCurr= gamMa_prob(span(gene_start, gene_end));
		    vec aTfCurr = aTf(span(gene_start, gene_end)) ;
		    if(prior >0.5){
			//aTfCurr = exp(aTfCurr);
			double sumTemp = sum(aTfCurr);
		        aTfCurr = aTfCurr/sumTemp;
			if((iter ==0) && (gene==0)) cout << "aTfCurr \t" << aTfCurr.t() << endl;
		        aTfCurr = logit(aTfCurr);
		    }
		    //cout << 763 << endl;
		    //cout << 763 << endl;
		    //vec Randn = randn<vec>(p_curr)*sigma_eff; // N(0, sigma_eff)
		    //vec Randn_0 = randn<vec>(p_curr)*sigma_0; Randn_0 = Randn_0 + thetamean * ones(p_curr); //  N(mu0, sigma0) 
		    double aTfbin;
		    //gamMa_probCurr  = gamMa_probCurr * (1-2*epiDiscount) + epiDiscount;

		    //double phi1;

		    if( qq > 0) {
			WW = solve(trimatl(TT_curr.t()), xTy(maskCurr(nonzero)));
			gg_curr = -qq/2.0 * log(1 + ccCurr) - (n + nu(gene))/2.0 * log (nu(gene) * lambdal(gene) + y2(gene)  - sum(square(WW))/(1 + 1/ccCurr) ) ; 
		    }else gg_curr = gg(gene);
		    for (int localiter = 0; localiter < localmaxiter; localiter++) {
			for(int snp = 0;  snp < p_curr; snp++){
			    //cout << qq_curr<< " here "  << 243 << gene << endl;
			    if(gamMaCurr(snp) ==1){
				qq_curr = qq-1;
				//TT1 = TT_curr;
				gamMaChange(snp) = 0;
				if(qq_curr > 0 ){
				    nonzero_curr= nonzero(span(1, qq_curr));  //find(gamMaChange > 0);
				    nonzero(qq-1) = nonzero(0); nonzero(span(0, qq-2)) = nonzero_curr; 
				    //cout << "here" << gene << 198 << endl;
				    TT_curr = cholGivensRotateLeft(TT_curr, xTx_gene(nonzero, nonzero)); //chol( x.cols(maskCurr(nonzero_curr)).t() * x.cols(maskCurr(nonzero_curr)) + small * eye(qq_curr, qq_curr));
				    //cout << "here" << gene << 192 << endl;
				    gg1 = gg_curr;
				    WW = solve(trimatl(TT_curr(span(0,qq_curr-1), span(0, qq_curr-1) ).t()), xTy(maskCurr(nonzero_curr))); 
				    //cout << "here" << gene << 19206 << endl;
				    gg0 = -qq_curr/2.0 * log(1 + ccCurr) - (n + nu(gene))/2.0 * log (nu(gene) * lambdal(gene) + y2(gene)  - sum(square(WW))/(1 + 1/ccCurr) );
				    //cout << "here" << gene << 19208<< endl;
				}else{
				    //TT0 = zeros<mat>(1,1);
				    nonzero_curr.clear();
				    //cout << "here" << gene << 19210<< endl;
				    gg0 =  - (n + nu(gene))/2 * log(nu(gene) * lambdal(gene)  + y2(gene) );
				    //cout << "here" << gene << 19213 << endl;
				}	    
			    }else{
				//TT0 = TT_curr; 
				//cout << "here" << gene << 19217 << endl;
				gamMaChange(snp) = 1;
				qq_curr = qq+1;
				gg0 = gg_curr;
				//cout <<"here 265 " <<  qq_curr <<nonzero.t() << endl;
				nonzero_curr = nonzero;
				nonzero_curr.resize(qq_curr); nonzero_curr(qq_curr-1) = snp;//  = find(gamMaChange > 0);
				//cout << "here" << nonzero_curr.t() << " "  << gene << 19220 << " " <<qq_curr << "  " << TT1.n_cols<< endl;
				if(qq >  0){
				    xTx_snp_gene = xTx_gene.col(snp);
				    TT1 = cholAddCol(TT_curr, xTx_snp_gene(nonzero_curr), xTx_gene(nonzero_curr, nonzero_curr));//chol( x.cols(maskCurr(nonzero_curr)).t() * x.cols(maskCurr(nonzero_curr)) + small * eye(qq_curr, qq_curr));
				} else{
				    TT1.set_size(1,1); tempd = sqrt(xTx_gene(snp,snp)); TT1(0,0) = (tempd > 0) ? tempd : (-1*tempd);
				}
				//cout << "here" << gene << 19222 << endl;
				WW =  solve(trimatl(TT1.t()), xTy(maskCurr(nonzero_curr))); //check there is no transponse
				gg1 = -qq_curr/2.0 * log(1 + ccCurr) - (n + nu(gene))/2.0 * log (nu(gene) * lambdal(gene) + y2(gene)  -  sum(square(WW))/(1 + 1/ccCurr) );
				//cout << "here" << gene << 19226 << endl;
			    }
			    //cout << gg0 << " here " << gg1 << sum(square(WW)) << endl;
			    //randDraw = randu(1);
			    // sample gamMa_prob 

			    aTfbin = logistic(aTfCurr(snp));
			    randDraw = randu(3);
			    //phi1 = ABS(aTfCurr(snp) - nu1);
			    gamMa_probsnpCurr = gamMa_probCurr(snp);

			    bfCurr(snp) =  logistic(gg1  - gg0  + logit(.05)); 
			    //cout << 863 << "\t" << gg1 <<  "\t" << gg0  << "\t" << gamMa_probCurr(snp)  << "\t" <<  prior <<  "\t"<<  gg1 - gg0 + gamMa_probCurr(snp)+ log(prior/(1-prior)) << endl;
			    //gamMa_rate(snp) = logistic(gg1 - gg0 + gamMa_probsnpCurr);
			    //double pa1x = min(a11 * gamMa_probsnpCurr + a10 * (1 - gamMa_probsnpCurr) , prior);
			    double pa1x = a11 * gamMa_probsnpCurr + a10 * (1 - gamMa_probsnpCurr);
			    //cout << " pa1x \t" << pa1x << endl;
			    if(prior < 0.5){
				    gamMa_rate(snp) = logistic(gg1 - gg0 +  aTfCurr(snp) + logit(prior));
			    }else{
				    gamMa_rate(snp) = logistic(gg1 - gg0 +  aTfCurr(snp));
			    }
			    //cout << gamMa_rate(snp) << endl;
			    gamMa_snpCurr = 0;
			    //gamMa_snpCurr = ((gamMa_rate(snp) >= randDraw(1)) & (bfCurr(snp) > randDraw(2)) ) ? 1:0;
			    if(dual_filter){
				if ((gamMa_rate(snp) >= randDraw(1)) && (aTfbin > randDraw(2)) ) gamMa_snpCurr = 1;
			    }else{
				//if((gamMa_rate(snp) >= randDraw(1)) && (bfCurr(snp) > randDraw(2)) ) gamMa_snpCurr = 1;
				if(gamMa_rate(snp) >= randDraw(1)  ) gamMa_snpCurr = 1;
			    }

			    //if((iter ==0) && (gene==0)) cout << snp << "\t" << gamMa_rate_snp<< aTfCurr.t() << endl;
			    double pax1 =  gamMa_snpCurr* log(a11/a10)  + (1-gamMa_snpCurr ) * log( (1 - a11)/(1 - a10)) ;  
			    // this may not be correct, should it be arthmetic mean or geometric mean, double pax1 =  gamMa_rate(snp)* log(a11/a10)  + (1-gamMa_rate(snp) ) * log( (1 - a11)/(1 - a10)) ;  
			    //double temp1 = log(theta_odds) + pax1 + aTfCurr(snp); //  - logit(.99); // a prior can be introduced in gamMa_prob/theta
			    double temp1 = - theta_bias + aTfCurr(snp); //  - logit(.99); // a prior can be introduced in gamMa_prob/theta
			    //if(temp1 > 8.0 ) temp1 = 8.0;
			    //if(temp1 <  -8.0 ) temp1 = -8.0;
			    gamMa_probCurr(snp) = (iter == 0)? gamMa_snpCurr : temp1;
			    //if((localiter ==0) && (gene==0) && (snp > 994 ) && (snp < 999)) 
			    //cout  << gene << " snp " << snp << " prior " <<  prior <<  " pa1x "<<  pa1x << " pax1 "<<  pax1 << "gamMa_rate(snp)" << gamMa_rate(snp) << " gamMa " <<  gamMa_snpCurr << " gamMa_probCurr(snp) " << gamMa_probCurr(snp) << " logit(gamMa_rate) " << gg1 - gg0 +  logit(prior) + logit(pa1x) <<" bf "  << bfCurr(snp) << " aTf " <<  aTfCurr(snp) << " theta_bias "<< theta_bias<<  endl;
			    //if( (gamMa_probCurr(snp)  == 1.0) |(gamMa_probCurr(snp) == 1.0)  ){
				//cout << "error probCurr out of range " << endl;
				//exit(-1);
			    //}



			    // sampling of theta
			    /*double mu = aTfCurr(snp)* sigma2_eff/sigma2_l;
			    if(ABS(aTfbin -  gamMaCurr(snp)) > 0.5 ){ // proposal N(0, s) 
				prob_proposal = Randn(snp);
				proposal_ratio = (gamMaCurr(snp) > 0) ? (logistic(prob_proposal )/logistic(gamMa_probsnpCurr)) : ((1 - logistic(prob_proposal))/(1- logistic(gamMa_probsnpCurr))); 
				//proposal_ratio	= exp( (1 - gamMa_rate(snp)) * ( log(1- logistic(prob_proposal)) - log(1- logistic(gamMa_probsnpCurr)   )) +
				//gamMa_rate(snp)  *( log(logistic(prob_proposal )) -  log(logistic(gamMa_probsnpCurr)) ));
				accept = proposal_ratio * dnorm(prob_proposal, mu, sigma_eff )/dnorm(gamMa_probsnpCurr, mu, sigma_eff )*\
				dnorm(gamMa_probsnpCurr, 0.0, sigma_eff)/dnorm(prob_proposal, 0.0, sigma_eff);
				if((accept > randDraw(0))&& (prob_proposal > thetamin) && (prob_proposal <=  thetamax)){
				    accept_count++;
				    gamMa_probCurr(snp) = prob_proposal; //accept 
				}
			    }else{ //proposal N(aTf, s)
				prob_proposal = Randn(snp) + mu;
				proposal_ratio = (gamMaCurr(snp) > 0) ? (logistic(prob_proposal )/logistic(gamMa_probsnpCurr)) : ((1 - logistic(prob_proposal))/(1- logistic(gamMa_probsnpCurr))); 
				//proposal_ratio	= exp( (1 - gamMa_rate(snp)) * ( log(1- logistic(prob_proposal)) - log(1- logistic(gamMa_probsnpCurr)   )) +
				//gamMa_rate(snp)  *( log(logistic(prob_proposal )) -  log(logistic(gamMa_probsnpCurr)) ));
				accept = proposal_ratio;
				if((accept > randDraw(0))&& (prob_proposal > thetamin) && (prob_proposal <=  thetamax)){
				    accept_count++;
				    gamMa_probCurr(snp) = prob_proposal; //accept
				}
			    }
			    //cout << "prob_proposal \t" << prob_proposal << "\t" <<  gamMa_probCurr(snp) << endl;
			    proposal_count++; * / sampling of theta 

			    //gamMa_probCurr(snp) = aTfCurr(snp); //accept
			    //gamMa_snpCurr = (gamMa_rate(snp) >= randDraw(1))  ? 1:0;
			    //cout << gene_start + snp <<"\t bfCurr\t" <<  bfCurr(snp) <<  "\t bfCurr\t " <<  gg0 - gg1 << "\t gamMa_prob \t " << gamMa_probCurr(snp) << "\t gamMa_rate\t "  << 
			    //gamMa_rate(snp) <<" \t gamMa \t"   << gamMa_snpCurr<< "\t prior_prob \t"<<log(prior/(1-prior))<<  endl;
			    /*if(((gene_start+ snp) > 32000)  && ((gene_start+ snp) < 35000)){
				double GG1, GG0 ; int QQ_curr; mat tt1;
				uvec NONzero; uvec temp12 = find(nonzero != snp);
				NONzero = nonzero(temp12);
				QQ_curr= NONzero.size();
				cout << NONzero.t() << "\tsnp \t" << snp<<endl; 
				tt1 = chol( x.cols(maskCurr(NONzero)).t() * x.cols(maskCurr(NONzero)) + small * eye(QQ_curr, QQ_curr)); 
				WW =  solve(trimatl(tt1.t()), xTy(maskCurr(NONzero))); //check there is no transponse
				GG1 = -QQ_curr/2.0 * log(1 + ccCurr) - (n + nu(gene))/2.0 * log (nu(gene) * lambdal(gene) + y2(gene)  -  sum(square(WW))/(1 + 1/ccCurr) ); 
				QQ_curr = QQ_curr + 1;	
				NONzero.resize(QQ_curr); NONzero(QQ_curr-1) = snp;//  = find(gamMaChange > 0);
				tt1 = chol( x.cols(maskCurr(NONzero)).t() * x.cols(maskCurr(NONzero)) + small * eye(QQ_curr, QQ_curr)); 
				WW =  solve(trimatl(tt1.t()), xTy(maskCurr(NONzero))); //check there is no transponse
				GG0 = -QQ_curr/2.0 * log(1 + ccCurr) - (n + nu(gene))/2.0 * log (nu(gene) * lambdal(gene) + y2(gene)  -  sum(square(WW))/(1 + 1/ccCurr) );
				cout << "check bfCurr \t"<< GG1 -GG0 <<endl;
			    }*/

			    //gamMa_rate(snp) = logistic(sum(geneEpiCurr.t() * tau) );
			    //cout << gamMa_rate(snp) << endl;
			    //gamMa_rate(snp) = logistic(gg1 - gg0 + log(gamMa_probCurr(snp)) - log(1 - gamMa_probCurr(snp)));	
			    //bf_snpCurr = logistic( bfCurr(snp)); 
			    //cout << "randDraw \t" << randDraw.t() << "gamMa_rate \t" << gamMa_rate(snp) << "\t bf_snpCurr" << endl;

			    /*if( (bf_snpCurr >= eQTL_threshold ) || (bf_snpCurr >= randDraw(0))) {
				if(gamMa_probCurr(snp) >= regulator_threshold){
				    //cout << "inside" << endl;
				    gamMa_snpCurr = (gamMa_rate(snp) >= randDraw(1)) ?1:0;
				    //if(gamMa_snpCurr > 0) cout << "success" << endl;
				}else  gamMa_snpCurr = 0;
				if(dual_filter){
				    gamMa_snpCurr = (gamMa_probCurr(snp) >= randDraw(1)) ?1:0;
				}
			    }else gamMa_snpCurr = 0;*/

			    //cout <<gamMa_snpCurr << endl;

			    //gamMa_snpCurr = (Rcpp::rbinom(1, 1, gamMa_rate(snp)))(0);
			    //gamMa_snpCurr = (gamMa_rate(snp) <.5 )? 0:1;
			    ////cout << gamMa_rate(snp) << " here " <<  gamMa_snpCurr << " old " << gamMaCurr(snp)  <<  endl;
			    if(gamMa_snpCurr != gamMaCurr(snp) ) {
				gamMaCurr(snp)   = gamMa_snpCurr; 
				if(gamMa_snpCurr==1){
				    nonzero = nonzero_curr; 
				    gg_curr = gg1;
				    qq = qq_curr;
				    TT_curr = TT1;
				}else{
				    nonzero = nonzero_curr;
				    gg_curr = gg0;
				    qq = qq_curr;
				    if(qq_curr >0){
					TT_curr = TT_curr(span(0,qq-1), span(0,qq-1));
				    }else TT_curr = zeros<mat>(1,1);
				}
			    }else gamMaChange(snp)   = gamMa_snpCurr;
			}
			if(use_raoblackwell) gamMa_inclusion = gamMa_inclusion + gamMa_rate;
			else gamMa_inclusion = gamMa_inclusion + gamMaCurr;
		    }
		    //cout << 932<<  "\t" << gene << endl;
		    if(sample_logit){
			if(oversample){
			    double sum_inclusion = sum(gamMa_inclusion);
			    if((localmaxiter*p_curr  > sum_inclusion) && (sum_inclusion > 0)){ 
				sample_prob = (localmaxiter*p_curr - sum_inclusion) * gamMa_inclusion +  sum_inclusion *(localmaxiter- gamMa_inclusion);
			    } else sample_prob.ones(); 
			    //cout << localmaxiter*p_curr<< "sample_prob " << sum_inclusion << endl;
			    //uvec ret = Rcpp::RcppArmadillo::sample(bootstrapS, bootstrap, true, sample_prob);
			    //ret = sample(p_curr,Named("size",bootstrap),Named("prob",sample_prob), Named("replace",true));
			    //sample.prob.ones(); // removed the oversampling and balancing of enhancers.
			    ransampl_set( ws, sample_prob.memptr() );
			    a1.randu(bootstrap); a2.randu(bootstrap);
			    for (int ii = 0; ii < bootstrap; ii++) {
				ret = ransampl_draw( ws, a1(ii), a2(ii));
				logit_freqCurr(ret) = 	logit_freqCurr(ret)+1; 
			    }

			}else if(!balance){
			    //logit_freqCurr = ones(p_curr);
			    sort_inx = arma::sort_index(gamMa_inclusion);
			    if(p_curr >= 50) {
				topp = 5; bottomm = 45; 
			    }else{
				//topp = std::max(1.0, .10 * (double) p_curr);
				//bottomm = std::max(1.0, .50 * (double) p_curr);
				topp = p_curr; bottomm=0;
			    }
			}
			//logit_freqCurr(sort_inx(span(p_curr - topp, p_curr-1) )) = ones<vec>(topp); 
			//logit_freqCurr(sort_inx(span(0, bottomm -1) )) = ones<vec>(bottomm); 
		    }
		    zero = find(gamMaCurr==0);
		    //VV =  1/(1/sigma2 + 1/cc)(X'X)^-1,  mm = VV X'Y    
		    //cout << 962<<  "\t" << gene << endl;
		    if(qq > 0) {
			TT_inv = solve(trimatu(TT_curr), eye(qq,qq)); 
			VV = ccCurr*sigma2Curr/(1+ccCurr) * (TT_inv * TT_inv.t() );
			//cout << "here" << gene << 19263 << endl;
			mm = VV * xTy(maskCurr(nonzero))/sigma2Curr;
			//cout << nonzero.t() << endl;
			//cout << "here beTa" << gene << " " << VV(0,0)  << endl;
			beTaCurr.zeros();
			beTaCurr(nonzero) = mvrnormArma(1, mm, VV);
			//cout << "here" << gene << 19267 << endl;
			xbeTa  = x.cols(maskCurr(nonzero)) * beTaCurr(nonzero);
			xbeTaM =  x.cols(maskCurr(nonzero)) * mm;
		    }else{
			beTaCurr.zeros(); xbeTa.zeros(); xbeTaM.zeros();

		    }
		    //cout << 979<<  "\t" << gene << endl;
		    if(qq > 1){
			mat TTtemp = chol(xTx_gene(nonzero, nonzero));
			if(sum(sum(TTtemp - TT_curr)) > .001){
			    cout << sum(sum(TTtemp))<< "\t" << sum(sum(TT_curr))  << endl;
			    cout << "mismatch TT \t " << nonzero.t() << endl;
			    cout << TT_curr << endl;
			    cout << TTtemp << endl;
			} 
		    }


		    residue = y.col(gene) - xbeTa; 
		    double ev = var(y.col(gene) - xbeTaM)/varY(gene); 
		    double weightCurr = 10* pow(2, -10*(1-ev)) + 1;
		    double 	  residue2sum = sum(square(residue)); double xbeta2sum = sum(square(xbeTa));
		    //cout << 996<<  "\t" << gene << endl;
		    //cout << "gamMa" << sum(gamMaCurr) <<endl;
		    //cout << ((double) accept_count)/((double) proposal_count) << endl; 

		    //#pragma omp critical(dataupdate)
		    {
			//sigma2Curr = rinvgammaC((qq + nu(gene)+ n )/2.0, 1/2.0 * (residue2sum + nu(gene) * lambdal(gene)));     
			//ccCurr = rinvgammaC((qq+1)/2.0, (xbeta2sum + n) /2.0);
			sigma2Curr = rinvgammaC(( nu(gene)+ n + qq)/2.0, 1/2.0 * (residue2sum + xbeta2sum/ccCurr +  nu(gene) * lambdal(gene)));     
			ccCurr = rinvgammaC((qq+1)/2.0, (xbeta2sum/sigma2Curr + n) /2.0);
			sigma2(gene) = sigma2Curr;
			cc(gene) = ccCurr;
			//cout << 1007<<  "\t" << gene << "\t"<<zz.size() <<"\t"<< gene_end << endl;
			evAll(gene) =1-ev;
			if(evAll(gene) < 0) { 
			    cout << ev << " ev neg for gene" <<  gene << endl;
			    cout <<   beTaCurr(nonzero) << endl;
			    cout <<  mm  << endl;
			    cout << nonzero.t() << endl;
			    cout << TT_curr << endl;
			}

			zz(span(gene_start,gene_end)) = gamMa_inclusion/localmaxiter; //gamMa_rate - 0.5; //# TODO Instead of gamMa, gamMa_rate can be used 
			//cout << 1017<<  "\t" << gene << endl;
			//zz(span(gene_start,gene_end)) = logistic(gamMa_probCurr); //gamMa_rate - 0.5; //# TODO Instead of gamMa, gamMa_rate can be used 
			if(sample_logit){
			    if( oversample){
				logit_freq(span(gene_start,gene_end)) = logit_freqCurr; //TODO change behaviour for overlapping snps 
			    }else if (!balance) {
				logit_freqtop(gene_start + sort_inx(span(p_curr - topp, p_curr-1))) = weightCurr* ones<vec>(topp); 
				logit_freq(gene_start + sort_inx(span(p_curr - topp, p_curr-1))) = weightCurr* ones<vec>(topp); 
				if(bottomm > 0) {
				    logit_freqbottom(gene_start + sort_inx(span(0, bottomm -1))) = weightCurr*ones<vec>(bottomm);
				    logit_freq(gene_start + sort_inx(span(0, bottomm -1))) = weightCurr*ones<vec>(bottomm);
				}
			    }
			}else{
			    //cout << p_curr<< "here " << gene_end - gene_start +1<<endl;
			    logit_freq(span(gene_start, gene_end)) = ones<vec>(gene_end - gene_start + 1);
			}
			//TODO instead of condition apply pre-processor conditional.

			gamMa(span(gene_start, gene_end)) = gamMaCurr; 
			gamMa_prob(span(gene_start, gene_end)) = gamMa_probCurr; 
			bf(span(gene_start, gene_end)) = bfCurr; 
			gamMa_rateAll(span(gene_start, gene_end)) = gamMa_rate; 
			beTa(span(gene_start, gene_end)) = beTaCurr; 
			//gg(gene) = gg_curr; // not required only when qq is zero it is not recalcualted 
			TT(gene) = TT_curr;
			if((((iter -burnIn +1) %thin)==0)  && (iter >= burnIn)){
			    gamMa_rate_mat(span(gene_start, gene_end),count)= gamMa_rate;
			    bf_mat(span(gene_start, gene_end),count)= bfCurr;
			    gamMa_mat(span(gene_start, gene_end),count)= gamMaCurr;
			    gamMa_prob_mat(span(gene_start, gene_end),count)= gamMa_probCurr;
			    beTa_mat( span(gene_start, gene_end), count) = beTaCurr;
			    cc_mat(gene, count) = cc(gene);
			    ev_mat(gene, count) = evAll(gene);
			    sigma2_mat(gene, count) = sigma2(gene);
			}
		    }

		    ransampl_free( ws );
		}
	    }
	    beTaAcc = beTaAcc + beTa;
	    //gamMaAcc = gamMaAcc + gamMa;
	    gamMa_rateAcc.col(iter%accIter) = gamMa;
	}
	//estimate_alpha = false;
	if(verbose) 
	{ 
	    t2 = std::clock();
	    cout << "start of alpha cycle \t " << iter << "\t" << ((float)t2 - (float) t1) / CLOCKS_PER_SEC << "s"<<   endl;
	    cout << "number of eeSNP \t" << sum(gamMa) << endl;
	}

	if(estimate_alpha){
	    //sampling of snp of genes.
	    if(sample_logit){
		uvec nonzero_logittop = find(logit_freqtop > 0 );
		uvec nonzero_logitbottom = find(logit_freqbottom > 0 );
		bool canSample = ((nonzero_logittop.size()  > (num_logit_postrain)) &&  (nonzero_logitbottom.size()  > (4*num_logit_postrain)))? true :  false;
		if((!oversample) && canSample){
		    logit_freq.zeros();
		    ret1 = sample(p, Named("size", num_logit_postrain), Named("prob", logit_freqtop), Named("replace", false));
		    for (int tt = 0; tt < num_logit_postrain; tt++) 
		    logit_freq(ret1[tt] -1) = 1;
		    ret2 = sample(p, Named("size", 4*num_logit_postrain), Named("prob", logit_freqbottom), Named("replace", false));
		    for (int tt = 0; tt < 4*num_logit_postrain; tt++) 
		    logit_freq(ret2[tt] -1) = 1;
		}	
	    }
	    gamMaABS = mean(gamMa_rateAcc,1);
	    uvec nonzero_logit;
	    if(balance){
		//cout<< "here" << gamMa_rateAcc.n_cols<<endl;
		if(usebeTaprob) beTaABS = mean(beTaAcc,1);
		//cout<< "here1" <<gamMaABS.size() << endl;
		//if((iter+1) >= accIter){
		    //if (((iter+1) % accIter)==0) { 
			//beTaABS = ABS(beTaAcc)/accIter;
			//gamMaABS = sum(gamMa_rateAcc);
		    //}
		//} else {
		    //beTaABS = ABS(beTa);
		    //gamMaABS = gamMa;
		//}
		colvec nonzero_prob;
		if (usebeTaprob) {
		    nonzero_logitInx = (( zz > gamMa_thres) % (beTaABS > beTa_thres)); 
		    nonzero_logit = find(nonzero_logitInx);
		}else nonzero_logit = find(zz > gamMa_thres);
		nonzero_prob = gamMaABS(nonzero_logit);
		uvec zero_logit = find(zz < 0.1);
		colvec zero_prob = ones(zero_logit.size());
		zero_prob = zero_prob - gamMaABS(zero_logit);
		int postrain_num = (num_logit_train > -1) ? min( nonzero_logit.size()*ratio , num_logit_train/2):  nonzero_logit.size();
		int negtrain_num = (num_logit_train > -1) ? min( zero_logit.size()*ratio , num_logit_train/2):  zero_logit.size();
		int train_num = min(postrain_num, negtrain_num);
		//negtrain_num = train_num*2;
		negtrain_num = 2000;
		//train_num = min(num_logit_train/2, train_num);
		//num_logit_postrain = min( num_logit_train/2,  nonzero_logit * ratio) ;
		//int num_logit_negtrain = min( num_logit_train/2,  (p - nonzero_logit)*ratio) ;
		logit_freq.zeros();
		geneEpi_freq.zeros();
		cout <<"train_num \t"<< postrain_num << "\t" << negtrain_num << "\t" << train_num << "\t"<< gamMa_thres    << endl;
		if(train_num > (f/2 +1) ){
		    ret1 = sample(nonzero_logit.size(), Named("size", train_num), Named("prob", nonzero_prob),Named("replace", false));
		    for (int tt = 0; tt < ret1.size(); tt++) 
		    logit_freq(nonzero_logit(ret1[tt] -1)) = logit_freq(nonzero_logit(ret1[tt] -1)) + 1;
		}else{
		    cout << "oops very few confident regulators to train" << endl;
		    ret1 = sample(nonzero_logit.size(), Named("size", f/2.0 + 1), Named("prob", nonzero_prob),Named("replace", true));
		    for (int tt = 0; tt < ret1.size(); tt++) 
		    logit_freq(nonzero_logit(ret1[tt] -1)) = logit_freq(nonzero_logit(ret1[tt] -1)) + 1;
		    //ret1 = sample(p, Named("size", f/2 + 1), Named("prob", nonzero_prob), Named("replace", true));
		    //for (int tt = 0; tt < ret1.size(); tt++) 
		    //logit_freq(ret1[tt] -1) = logit_freq(ret1[tt] -1) + 1;
		}
		ret2 = sample(p , Named("size", negtrain_num),  Named("replace", false));
		for (int tt = 0; tt < ret2.size(); tt++) 
		logit_freq(ret2[tt] -1) = logit_freq(ret2[tt] -1)+ 1;
		/*if(negtrain_num > (f/2 +1) ){
		    ret2 = sample(zero_logit.size(), Named("size", negtrain_num), Named("prob", zero_prob), Named("replace", false));
		    for (int tt = 0; tt < ret2.size(); tt++) 
		    logit_freq(zero_logit(ret2[tt] -1)) = logit_freq(zero_logit(ret2[tt] -1))+ 1;
		}else{
		    ret2 = sample(p, Named("size", (f/2 + 1)), Named("prob", 1-zz), Named("replace", false));
		    for (int tt = 0; tt < ret2.size(); tt++) 
		    logit_freq(ret2[tt] -1) = logit_freq(ret2[tt] -1)+ 1;
		}*/
		//genetic and epigenetic weight
		//ret1 = sample(nonzero_logit.size(), Named("size", train_num/2), Named("replace", false));
		//for (int tt = 0; tt < ret1.size(); tt++) 
		//geneEpi_freq(ret1[tt] -1) = geneEpi_freq(ret1[tt] -1) + 1;
		//ret2 = sample(p , Named("size", negtrain_num*2),  Named("replace", false));
		//for (int tt = 0; tt < ret2.size(); tt++) 
		//geneEpi_freq(ret2[tt] -1) = geneEpi_freq(ret2[tt] -1)+ 1;
	    }


	    //zz = logistic(gamMa_prob);
	    //zz = gamMaABS;
	    uvec freq_nonzero = find(logit_freq > 0); p_nonzero = freq_nonzero.size();
	    cout << " p_nonzero " <<  p_nonzero << endl;
	    mat feature_train1 = featureCat.rows(freq_nonzero);
	    bool use_logisitic = true;
	    if (use_logisitic) {


		//adjusting gamMa_prob
		/*if ((regulator_prior > 0) && (iter > 0)){
		    uvec temp1= find(gamMa_prob > 0.0);
		    int ynum =  (temp1).size();
		    cout<< "number of regulatory enhancers \t " << ynum << endl;
		    vec gamMa_sort = sort(gamMa_prob);
		    int reg_inx = regulator_prior * p;
		    gamMa_prob = logistic(gamMa_prob - gamMa_sort(p - reg_inx)); 
		    temp1=find(gamMa_prob > 0.5);
		    cout<< "number of regulatory enhancers after adjustment \t " << temp1.size() << endl;
		}*/

		//zz = gamMa_prob;
		zz = zz - 0.5;
		zz(freq_nonzero) = logit_freq(freq_nonzero) % zz(freq_nonzero);

		// variable selection of alpha
		if(logistic_variable_selection){
		    mat feature_train = feature_train1.cols(span(1,(f-1))); // ones are not added anymore 
		    logistic_variable_selection_mcmc( feature_train, zz(freq_nonzero),  logistic_select, cc_epi, sigma2_epi, lambdal_epi, nu_epi);
		    //else logistic_select.ones(); 
		}

		//#sample wi
		cout <<  "logisitic regression " <<endl;
		int logitIter = 20;
		vec alpHaMean=zeros(f);
		uvec alpha_select = find(logistic_select > 0);
		int f_curr = alpha_select.size() + 1;
		if(f_curr > 1){
		    uvec alpha_select1(f_curr); alpha_select1(0) = 0; alpha_select1(span(1, f_curr-1)) =  alpha_select + 1;   
		    //uvec freq_nonzero = find(logit_freq > 0); p_nonzero = freq_nonzero.size();
		    for (int logitCnt = 0; logitCnt < logitIter; logitCnt++) {
			feature_nonzero =  feature_train1.cols(alpha_select1); // first column (bias) is always part of the model   
			logit_freq_nonzero(span(0,p_nonzero-1)) = logit_freq(freq_nonzero);
			fTalpHa = feature_nonzero * alpHa(alpha_select1);
			rpg_hybrid( &w(0), &logit_freq_nonzero(0), &fTalpHa(0), &p_nonzero);
			//calculate m_w V_w;
			mat diag_w  = diagmat(w(span(0, p_nonzero-1 )));
			mat V_w = ((trans(feature_nonzero) * diag_w * feature_nonzero)+  diagmat( B_inv_alpHa(alpha_select1) )).i();
			vec m_w = V_w * ( (trans(feature_nonzero) * zz(freq_nonzero)) +   B_inv_m_w_prior(alpha_select1));
			alpHa.zeros();
			alpHa(alpha_select1)  = trans(mvrnormArma(1, m_w, V_w));
			alpHaMean = alpHaMean + alpHa;
			/*if( logitCnt ==0) {
			    aTf = featureCat * alpHa;
			    double aTfthres = mean(aTf); aTfthres = max(aTfthres, 0);
			    uvec sortInx = sort_index(aTf);
			    logit_freq.zeros();
			    int train_num = min(num_logit_train/2, p/10);
			    uvec tempInx = sortInx(span(0, train_num-1));
			    uvec labelInx = tempInx( find(zz(tempInx) < 0)); 
			    logit_freq(labelInx) = ones(labelInx.size());
			    cout<<"negtrain size \t " << labelInx.size() ;
			    tempInx = sortInx(sortInx(span(p - train_num, p-1)));
			    labelInx = tempInx( find(zz(tempInx) > 0)); 
			    cout<<"\t postrain size \t " << labelInx.size() << endl;
			    logit_freq(labelInx) = ones(labelInx.size());
			    labelInx = nonzero_logit( find(aTf(nonzero_logit) > aTfthres));
			    cout<<"\t postrain size  \t " << labelInx.size() << endl;
			    logit_freq(labelInx) = ones(labelInx.size());
			    freq_nonzero = find(logit_freq > 0); p_nonzero = freq_nonzero.size(); 
			    zz(freq_nonzero) = logit_freq(freq_nonzero) % zz(freq_nonzero);
			    feature_train1 = featureCat.rows(freq_nonzero);
			}*/
		    }
		    alpHa = alpHaMean/logitIter;
		    //gamMa_prob = logistic(featureCat * alpHa)
		}else{
		    alpHa.zeros();
		}
	    }else{ // theta regression
		// input gamMa_probs, featureCat
		// output alpHa
		vec alpHa_nobias = alpHa(span(1,f-1));
		//	zz = gamMa_prob(freq_nonzero);
		//cout << "theta_mean \t" << endl;
		//cout<< "gamMa_prob\t " << gamMa_prob(span(1,10)) << endl; 
		vec y_l = logit(gamMa_prob(freq_nonzero));
		double aTfbias = mean(y_l);
		//cout << "y"<< (y_l(span(1,10))).t() << endl;
		y_l = y_l - ones(p_nonzero) * aTfbias;
		//cout << "y"<< (y_l(span(1,10))).t() << endl;
		mat xx = feature_train1.cols(span(1,(f-1))); // ones are not added anymore
		//cout << "xx"<< (xx.col(0)).t() << endl;
		mat V_w =  cc_l*sigma2_l/(1+cc_l) * ((xx.t()* xx + small*eye(f-1, f-1)).i());
		//cout << "V_w"<< V_w.row(0) << endl;
		vec m_w = V_w * xx.t() * y_l/sigma2_l;
		//cout << "m_w"<< m_w.t() << endl;
		alpHa_nobias = trans(mvrnormArma(1, m_w, V_w));
		vec xalpHa  = xx * alpHa_nobias;
		vec residue = y_l - xalpHa; 
		double 	  residue2sum = sum(square(residue)); double xalpHa2sum = sum(square(xalpHa));
		sigma2_l = rinvgammaC(( nu_l+ f - 1 + p_nonzero)/2.0, 1/2.0 * (residue2sum + xalpHa2sum/cc_l +  nu_l * lambdal_l));     
		cc_l = rinvgammaC(f/2.0, (xalpHa2sum/sigma2_l + p_nonzero) /2.0);
		sigma_l = sqrt(sigma2_l);
		alpHa(span(1,f-1)) = alpHa_nobias;
		alpHa(0) = aTfbias;
		aTf = featureCat * alpHa;
	    }






	    //genetic epigenetic weight
	    /* cout << 1011 << endl;
	    freq_nonzero = find(geneEpi_freq > 0); p_nonzero = freq_nonzero.size(); 
	    cout << 1012 << endl;
	    zz(freq_nonzero) = geneEpi_freq(freq_nonzero) % zz(freq_nonzero);
	    cout << 1014) << endl;
	    feature_train1 = featureCat.rows(freq_nonzero);
	    cout << 1014 << endl;
	    mat geneEpi(p_nonzero, 4);
	    cout << 1015 << endl;
	    if(f_curr > 1){
		feature_nonzero =  feature_train1.cols(alpha_select + 1); // alpHa select don't contain bias. First column  not part of the model   
		cout << 1016 << endl;
		logit_freq_nonzero(span(0,p_nonzero-1)) = geneEpi_freq(freq_nonzero);
		cout << 1017 << endl;
		fTalpHa = feature_nonzero * alpHa(alpha_select + 1 );
		cout << 1018 << endl;
		geneEpi.col(0) = ones(p_nonzero,1);  geneEpi.col(1) =  fTalpHa;geneEpi.col(2) =  bf(span(1,p_nonzero)); geneEpi.col(3) =  bf(span(1,p_nonzero))%fTalpHa;
		cout << 1019 << endl;
		geneEpiTtau = geneEpi * tau;
		cout << 1010 << endl;
		rpg_hybrid( &w(0), &logit_freq_nonzero(0), &geneEpiTtau(0), &p_nonzero);
		cout << 1011 << endl;
		diag_w  = diagmat(w(span(0, p_nonzero-1 )));
		cout << 1012 << endl;
		V_w = ((trans(geneEpi) * diag_w * geneEpi)+  diagmat( B_inv_alpHa(span(1,4)) )).i();
		cout << 1013 << endl;
		m_w = V_w * (trans(geneEpi) * zz(freq_nonzero));
		cout << 1014 << endl;
		tau.zeros();
		cout << 1015 << endl;
		tau  = trans(mvrnormArma(1, m_w, V_w));
		cout << 1016 << endl;
	    }
	    mat tempGeneEpi(p, 4);
	    alpHa(0) = 0;
	    fTalpHa = featureCat * alpHa;
	    tempGeneEpi.col(0) = ones(p,1);  tempGeneEpi.col(1) =   fTalpHa;tempGeneEpi.col(2) =  bf; tempGeneEpi.col(3) =  bf%fTalpHa;
	    geneEpiTtau = tempGeneEpi * tau;
	    freq_nonzero = find(geneEpiTtau > 0);
	    cout << tau.t() << endl;
	    cout << freq_nonzero.size() << endl; 	    
	    */


	    //weight readjustment
	    if(regulator_prior > 0){
		aTf = featureCat * alpHa;
		uvec temp1= find((logistic(aTf)) > 0.5);
		int ynum =  (temp1).size();
		cout<< "number of enhancers \t " << ynum << endl;
		//double ym = ((double) ynum)/p; // King et. al eqn. 7.
		//alpHaAdj = alpHa;
		//alpHaAdj(0) =  alpHa(0) - log( ym*(1-regulator_prior)/(regulator_prior*(1 -ym))); 
		//double aTfthres = mean(aTf); 
		/* {vec aTf_sort = sort(aTf);
		    int reg_inx = regulator_prior * p;
		    //theta_bias =  aTf_sort(p - reg_inx); 
		    //theta_bias = 0.0; //
		    //alpHa(0) = alpHa(0) - log( ym*(1-regulator_prior)/(regulator_prior*(1 -ym)));
		    //vec alpHa1 = alpHa;
		    alpHa(0) = alpHa(0) - aTf_sort(p - reg_inx); 
		    temp1=find((logistic(featureCat * alpHa)) > 0.5);
		cout<< "number of enhancers after adjustment \t " << temp1.size() << endl; }*/
		//alpHa(0) = alpHaAdj(0);
	    }
	    //uvec temp1=find((logistic(featureCat * alpHa)) > 0.5);
	    //cout<< "number of enhancers \t " << temp1.size() << endl;
	    cout<<"alpHa \t"<< alpHa.t() <<endl;

	}
	if((((iter - burnIn + 1) %thin)==0)  && (iter >= burnIn)){
	    alpHa_mat.col(count) = alpHa;
	    //tau_mat.col(count) = tau;
	    //logistic_select_mat.col(count) = logistic_select;
	    count = count +1;
	}

	out["alpHa.mat"] = trans(alpHa_mat); out["beTa.mat"] = trans(beTa_mat); out["gamMa.mat"] = trans(gamMa_mat); out["gamMa_rate.mat"] = trans(gamMa_rate_mat); out["gamMa.prob.mat"] = trans(gamMa_prob_mat); out["cc.mat"] = trans(cc_mat);out["sigma2.mat"] = trans(sigma2_mat); out["ev"] = trans(ev_mat);
	out["bf.mat"] = trans(bf_mat);	
	if(logistic_variable_selection)  out["logistic_select"] = logistic_select_mat.t(); 
    }
    if(verbose) 
    { 
	t1 = std::clock();
	cout << "end of program \t " << "\t" << ((float)t1 - startOfepi_eQTL) / CLOCKS_PER_SEC << "s"<< endl;
    }

    return out;
}


