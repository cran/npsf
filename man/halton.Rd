\name{halton}
\alias{halton}
\title{
  'halton' method for class 'npsf'
}

\description{
  Provides Halton draws, deviates from a uniform distribution.
}

\usage{
  halton(n = 1, n.bases = 1, bases = NULL, 
                   start = 0, random.primes = FALSE, seed = 7, 
                   scale.coverage = FALSE, shuffle = FALSE)
}

\arguments{
  \item{n}{
    number of prime numbers to be returned (the row number in the value).
  }
  \item{n.bases}{
    numeric. number of bases used.  (the column number in the value).
  }
  \item{bases}{
    numeric. Supply specific order numbers for getting primes, see \code{\link{primes}}.  See examples.
  }
  \item{start}{
    numeric. from which value in the halton sequence to start.  Default is 0, which is actually 0.
  }
  \item{random.primes}{
    logical. if \code{TRUE}, the \code{n.bases} primes are chosen on a random basis from 100008 available prime numbers. See \code{\link{primes}}.
  }
  \item{seed}{
  set seed for replicability. Default is 17345168.
 }
 \item{scale.coverage}{
    logical. at larger primes not whole [0,1] interval is covered. if \code{TRUE}, \code{rescale} is used to fill the coverage.
  }
  \item{shuffle}{
    logical. if \code{TRUE}, each column in the value is randomly reshuffled (\code{seed} is used).
  }
}

\value{
  \code{halton} returns Halton draws.
}

\seealso{
  \code{\link{sf}} and \code{\link{primes}}.
}

\author{
  Oleg Badunenko <oleg.badunenko@brunel.ac.uk>
}

\examples{
  require( npsf )
  
  # obtain first 10 x 7 matrix with the first 7 primes as bases
  npsf::halton(10, n.bases = 7)
  
  # obtain first 10 x 7 matrix with the randomly chosen 7 primes as bases
  npsf::halton(10, n.bases = 7, random.primes = TRUE, seed = 17345168)
  
  # just one column with desired prime
  npsf::halton(10, bases = 1)
  
  # or 2 columns
  npsf::halton(10, bases = c(1,7))
  
  # if bases are large
  npsf::halton(10, bases = c(1,7)*1000)
  
  # the coverage is not great
  npsf::halton(10, bases = c(1,7)*1000, scale.coverage = TRUE)
  
  # reshuffle them, use seed for replicability
  npsf::halton(10, bases = c(1,7)*1000, scale.coverage = TRUE, shuffle = TRUE, seed = 17345168)
  
}

\keyword{ Stochastic Frontier Analysis }
\keyword{ Model }
\keyword{ Truncated normal }
\keyword{ Halton draws }
\keyword{ Prime numbers }
\keyword{ Parametric efficiency analysis }
\keyword{ summary }