\name{vcov.npsf}
\alias{vcov.npsf}
\title{
 'vcov' method for class 'npsf'
}

\description{
 Extracts the variance-covariance matrix of the ML parameters of a stochastic frontier model estimated by \code{\link{sf}}.
}

\usage{
 \method{vcov}{npsf}( object, \dots )
}

\arguments{
\item{object}{
 an object of class \code{npsf} returned by the function \code{\link{sf}}.
}
\item{\dots}{
 currently unused.
}
}

\details{
 The estimated variance-covariance matrix of the ML parameters 
 is the inverse of the negative Hessian evaluated at the MLE.
}

\value{
 \code{vcov.npsf} returns the estimated variance-covariance matrix of the ML parameters.
}

\seealso{
 \code{\link{coef.npsf}}, \code{\link{nobs.npsf}}, \code{\link{summary.npsf}}, and \code{\link{sf}}.
}

\author{
 Oleg Badunenko <obadunen@uni-koeln.de>
}

\examples{
require( npsf )
 
# Load Penn World Tables 5.6 dataset
 
data( pwt56 )
head( pwt56 )
 
# Create some missing values
 
pwt56 [4, "K"] <- NA 
 
# Stochastic production frontier model with 
# homoskedastic error components (half-normal)
 
# Use subset of observations - for year 1965
 
m1 <- sf(log(Y) ~ log(L) + log(K), data = pwt56, 
 subset = year == 1965, distribution = "h")
vcov( m1 )
}

\keyword{ Stochastic Frontier Analysis }
\keyword{ Model }
\keyword{ Truncated normal }
\keyword{ Parametric efficiency analysis }
\keyword{ Variance-covariance matrix }
\keyword{ summary }