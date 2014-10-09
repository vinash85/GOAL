\name{generate.epieQTL.logistic}
\alias{generate.epieQTL.logistic}
\alias{generate.epieQTL.logistic}
\title{Simulated data for GOAL}
\description{
  Generates simulated data for GOAL.
}

\usage{
  generate.epieQTL.logistic ( numSNP = 512, avg.eSNP=20,numGenes = 1, numSamp=20, numFeature=9,
  featureStat= list(enhancer=list(mean=rep(2,numFeature), cov= 2* diag(numFeature)),
  random=list(mean=rep(1,numFeature), cov=4* diag(numFeature))), seed=NULL,   evSNP =0.025)
}
\arguments{
  \item{numSNP}{The number of SNP per gene P/G.}
  \item{avg.eSNP}{The number of expression regulators per gene.}
  \item{numGenes }{The number of genes (G).}
  \item{numSamps }{The number of samples (N).}
  \item{numFeature }{The number of epigenetic features F.}
  \item{seed}{Value of seed for randomization.}
  \item{ev}{Explained variance by eeSNP for each gene}
  \item{LD}{Should I generate SNPs in LD (perfect LD)?}
}
\value{
  \code{generate.epieQTL.logistic ()} Returns an object:
  \item{x}{An N x P SNP matrix; with N samples and P SNPs.}
  \item{y}{An N x G gene matrix; with N samples and G genes.}
  \item{y.predict}{An N x G predicted gene matrix, component of expression predicted by SNP.}
  \item{feature}{A P x FALSE feature matrix corresponding to P SNPs; FALSE is number of features.}
  \item{alpHa}{An FALSE dimensional  feature importance parameters.}
  \item{gamMa}{An P dimensional inclusion for each SNPs.}
  \item{beTa}{An P dimensional regression weights SNPs.}
  \item{gamMa.prob}{An P dimensional interacting potentials vector.}
}
\details{
  generate.epieQTL.logistic generates a simulated data for eQTeL model. It generates
  SNPs matrix, gene matrix, feature matrix and also real parameters of eQTeL model (alpHa,beTa, gamMa and gamMa.prob). 
 }
\references{
  \enumerate{
    \item Avinash Das, MagNet Consortium, Shane Jensen, Sridhar Hannenhalli. Bayesian integration of genetics and epigenetics detects causal regulatory SNPs underlying expression variability. (Unpublished).
  }

}

\examples{
featureStat= list(enhancer=list(mean=rep(2, 7), cov=diag(7)) ,
		  random= list(mean=rep(0, 7), cov=2*diag(7)))
obj <- generate.epieQTL.logistic( numSNP = 512, avg.eSNP=5,numGenes = 100, 
numSamp=1000, numFeature=7, featureStat= featureStat, seed=NULL, evSNP =0.15)
}
\author{Avinash Das }

