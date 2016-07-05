\name{summary.npsf}
\alias{summary.npsf}
\alias{print.summary.npsf}
\title{
 'summary' method for class 'npsf'
}

\description{
 Prints summary of SF or DEA model estimated by \code{\link{sf}}, 
 \code{\link{teradial}}, \code{\link{tenonradial}}, and 
 \code{\link{teradialbc}}, or 
 testing procedures \code{\link{nptestrts}} and \code{\link{nptestind}}.
}

\usage{
 \method{summary}{npsf}( object, \dots )
 \method{print}{summary.npsf}( x, digits = NULL, print.level = NULL, \dots )
}

\arguments{
\item{object}{
 an object of class \code{npsf} returned by one of the functions 
 \code{\link{sf}}, \code{\link{teradial}}, \code{\link{tenonradial}}, 
 \code{\link{teradialbc}}, \code{\link{nptestrts}} or \code{\link{nptestind}}.
}
\item{x}{
 an object of class \code{npsf} returned by one of the functions 
 \code{\link{sf}}, \code{\link{teradial}}, \code{\link{tenonradial}}, 
 \code{\link{teradialbc}}, \code{\link{nptestrts}} or \code{\link{nptestind}}.
}
\item{digits}{
 numeric. Number of digits to be displayed in estimation results and for efficiency estimates. Default is 4.
}
\item{print.level}{
numeric. 0 - nothing is printed; 1 - print summary of the model and data. 2 - print summary of technical efficiency measures. 3 - print estimation results observation by observation (for DEA models). Default is 1.
}
\item{\dots}{
 currently unused.
}
}

\details{
 The summary depends on the model or testing procedure that is being estimated
}

\value{
 Currently no value is returned
}

\author{
 Oleg Badunenko <obadunen@uni-koeln.de>
}

\seealso{
 \code{\link{sf}}, \code{\link{teradial}}, \code{\link{tenonradial}}, 
 \code{\link{teradialbc}}, \code{\link{nptestrts}}, and \code{\link{nptestind}}
}

\examples{
require( npsf )
 
# Load Penn World Tables 5.6 dataset
 
data( pwt56 )
 
# Stochastic production frontier model with 
# homoskedastic error components (half-normal)
 
# Use subset of observations - for year 1965

# DEA

t1 <- teradialbc ( Y ~ K + L, data = pwt56, subset = Nu < 10, 
 reps = 199, print.level = 0)
summary(t1)

# SFA
 
m1 <- sf(log(Y) ~ log(L) + log(K), data = pwt56, 
 subset = year == 1965, distribution = "h",
 print.level = 0)
summary( m1 )

# Load U.S. commercial banks dataset
 
data(banks05)

m3 <- sf(lnC ~ lnw1 + lnw2 + lny1 + lny2, uhet = ~ ER, 
 vhet = ~ LA, data = banks05, distribution = "t", 
 prod = FALSE, print.level = 0)
 
summary(m3)

}

\keyword{ Stochastic Frontier Analysis }
\keyword{ Heteroskedasticity }
\keyword{ Truncated normal }
\keyword{ Parametric efficiency analysis }
\keyword{ Data Envelopment Analysis }
\keyword{ Nonradial measure }
\keyword{ Radial measure }
\keyword{ Nonparametric efficiency analysis }
\keyword{ Linear programming }
\keyword{ Bootstrap }
\keyword{ Smoothed bootstrap }
\keyword{ Subsampling bootstrap }
\keyword{ Parallel computing }
\keyword{ summary }