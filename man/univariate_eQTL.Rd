\name{univariate_eQTL}
\alias{univariate_eQTL}
\alias{eQTL}
\title{Univariate eQTL}
\description{
Correlation based univariate eQTL. Calculate explained variation. 
}

\usage{
univariate_eQTL(mask2, x, y, numRegulators, threads = 1L)
}
\arguments{
  \item{x}{An N x P SNP matrix; with N samples and P SNPs.}
  \item{y}{An N x G gene matrix; with N samples and G genes.}
  \item{mask2}{A P x 2 matrix that maps which SNPs are considered for which genes. Column 1: genes, column 2: SNPs}
  \item{alpHa}{An F+F' dimensional initial feature importance parameters.}
  \item{gamMa}{An F+F' dimensional initial inclusion for each SNPs.}
  \item{numRegulators}{A vector of number of regulator for each gene.}
  \item{threads}{The number of maximum threads.}
}
\value{
  \code{univariate_eQTL()} Returns a vector of explained variance  for each gene.  
}
\details{
  univariate_eQTL is standard univariate eQTL analysis. It calculates correlation between SNPs and gene pair and estimates explained variance.    
}
\references{
  \enumerate{
    \item Avinash Das, MagNet Consortium, Shane Jensen, Sridhar Hannenhalli. Bayesian integration of genetics and epigenetics detects causal regulatory SNPs underlying expression variability. (Unpublished).
  }

}

\seealso{
 \code{\link[GOAL:multivariate_eQTL]{multivariate_eQTL}} 
}
\examples{
featureStat= list(enhancer=list(mean=rep(2, 7), cov=diag(7)) ,
		  random= list(mean=rep(0, 7), cov=2*diag(7)))
obj <- generate.epieQTL.logistic( numSNP = 512, avg.eSNP=5,numGenes = 100, numSamp=1000, numFeature=7, featureStat= featureStat, seed=NULL, evSNP =0.15)
obj$pairFeature.mat <- matrix(0,0,0 )
x1 <- obj$x
x1[is.na(obj$x)] <- 0
ev <- univariate_eQTL(x=x1, y=obj$y, mask=obj$mask2, numRegulators=rep(3, ncol(obj$y)), threads=1)
}
\author{Avinash Das }

