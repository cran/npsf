\name{nobs.npsf}
\alias{nobs.npsf}
\title{
 'nobs' method for class 'npsf'
}

\description{
 Extracts the number of observations for which efficiencies are estimated 
 by SF or DEA model estimated by \code{\link{sf}}, \code{\link{teradial}}, 
 \code{\link{tenonradial}}, or \code{\link{teradialbc}}.
}

\usage{
 \method{nobs}{npsf}( object, \dots )
}

\arguments{
\item{object}{
 an object of class \code{npsf} returned by the function \code{\link{sf}}).
}
\item{\dots}{
 currently unused.
}
}

\value{
 \code{nobs.npsf} returns the number of observations for which efficiencies are estimated by SF or DEA model.
}

\seealso{
 \code{\link{vcov.npsf}}, \code{\link{coef.npsf}}, \code{\link{summary.npsf}}, and \code{\link{sf}}.
}

\author{
 Oleg Badunenko <oleg.badunenko@brunel.ac.uk>
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
nobs( m1 )

# DEA

t1 <- teradialbc ( Y ~ K + L, data = pwt56, subset = Nu < 10)
nobs(t1)
}

\keyword{ Stochastic Frontier Analysis }
\keyword{ Model }
\keyword{ Truncated normal }
\keyword{ Parametric efficiency analysis }
\keyword{ summary }