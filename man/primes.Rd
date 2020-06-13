\name{primes}
\alias{primes}
\title{
  'primes' method for class 'npsf'
}

\description{
  Provides prime numbers.
}

\usage{
  primes(n = NULL, which = NULL, random.primes = FALSE, seed = 7)
}

\arguments{
  \item{n}{
    number of prime numbers to be returned. Should be smaller than 100008.
  }
  \item{which}{
    numeric. if specific prime numbers are requred.  See examples.
  }
  \item{random.primes}{
    logical. if \code{n} is supplied and \code{random.primes = TRUE}, the  \code{n} primes are chosen on a random basis from 100008 available prime numbers.
  }
  \item{seed}{
  set seed for replicability. Default is 17345168.
 }
}

\details{
  \code{primes} just returns prime numbers, which come from \url{https://primes.utm.edu/lists/small/100000.txt}, see \url{ https://primes.utm.edu}
}

\source{
\url{https://primes.utm.edu/lists/small/100000.txt} and \url{ https://primes.utm.edu}
}

\value{
  \code{primes} returns prime numbers.
}

\seealso{
  \code{\link{sf}} and \code{\link{halton}}.
}

\author{
  Oleg Badunenko <oleg.badunenko@brunel.ac.uk>
}

\examples{
  require( npsf )
  
  # obtain first 30 prime numbers
  npsf::primes( 30 )
  
  # the same as
  npsf::primes( n = 30 )
  
  # the result in both case above are 30 prime numbers
  
  # if we use
  npsf::primes( which = 30 )
  
  # the 30th prime is returns, just a scalar
  
  # both cannot be used
  # npsf::primes(n = 30, which = 30, random.primes = FALSE, seed = 17345168)
  # will give a mistake
  
  # you can get random 30 primes, use seed for replicability
  npsf::primes(n = 30, random.primes = TRUE, seed = 17345168)
  
  # obtain specific primes: which take order number(s)
  npsf::primes(which = c(3,67,30, 100008))

}

\keyword{ Stochastic Frontier Analysis }
\keyword{ Model }
\keyword{ Halton draws }
\keyword{ Prime numbers }
\keyword{ Truncated normal }
\keyword{ Parametric efficiency analysis }
\keyword{ summary }