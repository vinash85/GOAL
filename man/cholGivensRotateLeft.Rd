\name{cholGivensRotateLeft}
\alias{cholGivensRotateLeft}
\title{  Updating Cholesky decomposition when design matrix is rotated left}
\description{
  The function updates Cholesky decomposition when design matrix is rotated left. This is useful  in linear regression when covariate is already in design matrix.
}

\usage{
cholGivensRotateLeft(TT, xx)
}
\arguments{
  \item{TT}{Old cholesky matrix.}
  \item{xx}{Column to added to design matrix.}
}
\value{
  \code{cholGivensRotateLeft()}   returns updated Cholesky decomposition. 
}
\details{
 Using Givens rotation cholGivensRotateLeft efficiently updates a Cholesky decomposition of a design matrix when its columns are rotated left. 
}
\references{
  \enumerate{
    \item Dongarra, J. J., Bunch, J. R., Moler, C. B. & Stewart, G. W. LINPACK users’ guide, vol. 8 (Siam, 1979).
  }

}

\examples{
}
\seealso{
 \code{\link[GOAL:cholAddCol]{cholAddCol}} 
}
\author{Avinash Das }

