\name{GOAL-package}
\alias{GOAL-package}
\alias{GOAL-package}
\docType{package}
\title{
  Genetics Omics Assciation resolve Linkage (GOAL)
}
\description{
GOAL implements eQTeL model that integrates epigenetic and genetic to infer causal
SNP that regulate expression. GOAL estimates regulatory potential and interaction 
potential of a SNP, and use this as prior of SNP to be expression regulators. 
}
\details{
\tabular{ll}{
Package: \tab GOAL\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-10-08\cr
License: \tab GPL (>=3.0)\cr
}
The main function of GOAL, the input is 1) SNP matrix 2) gene matrix 3) epigenetic feature, 4) pair feature, and 5) LD block. A simulate data can be generate using generate.epieQTL.logistic function.
}
\author{
Avinash Das Sahu <vinash85@umiacs.umd.edu>, Sridhar Hannenhalli <sridhar@umiacs.umd.edu>
}
\references{
     Avinash Das, MagNet Consortium, Shane Jensen, Sridhar Hannenhalli. Bayesian integration of genetics and epigenetics detects causal regulatory SNPs underlying expression variability. (Unpublished).

}
\keyword{ eQTL }
\seealso{
 \code{\link[GOAL:GOAL]{GOAL}} 
}

\examples{

}
