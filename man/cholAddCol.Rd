\name{cholAddCol}
\alias{cholAddCol}
\title{  Updating Cholesky decomposition on addition of column to design matrix}
\description{

 The function updates Cholesky decomposition on addition of column to design matrix
}

\usage{
cholAddCol(TT, xx)
}
\arguments{
  \item{TT}{Old cholesky matrix.}
  \item{xx}{Column to added to design matrix.}
}
\value{
  \code{cholAddCol()}   returns updated Cholesky decomposition. 
}
\details{
  Using Givens rotation cholAddCol efficiently updatea a Cholesky decomposition of a design matrix when a column is added to the design matrix. 
}
\references{
  \enumerate{
    \item Dongarra, J. J., Bunch, J. R., Moler, C. B. & Stewart, G. W. LINPACK users’ guide, vol. 8 (Siam, 1979).
  }

}

\examples{
}
\seealso{
 \code{\link[GOAL:cholGivensRotateLeft]{cholGivensRotateLeft}} 
}
\author{Avinash Das }

