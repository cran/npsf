\name{rescale}
\alias{rescale}
\title{
  'rescale' method for class 'npsf'
}

\description{
  rescales a vector.
}

\usage{
  rescale(x, lb = min(x), ub = max(x))
}

\arguments{
  \item{x}{
    a numeric vector.
  }
  \item{lb}{
    numeric. lower bound.
  }
  \item{ub}{
  numeric. upper bound.
 }
}

\value{
  \code{rescale} returns rescaled vector.
}

\seealso{
  \code{\link{sf}}, \code{\link{primes}}, and \code{\link{halton}}.
}

\author{
  Oleg Badunenko <oleg.badunenko@brunel.ac.uk>
}

\examples{
  require( npsf )
  
  # obtain first 30 prime numbers
  set.seed(8265897)
  t1 <- runif(10, min = 1, max = 2)
  summary(t1)
  summary(rescale(t1, 0, 10))
  
}

\keyword{ Halton draws }
\keyword{ Prime numbers }