\name{pwt56}
\docType{data}
\alias{pwt56}

\title{Penn World Tables 5.6 (compiled in 1995)}

\description{
 The data set is from Penn World Tables (PWT) 5.6.  This data set provides only selected variables for years 1965 and 1990.
}

\usage{data( pwt56 )}

\format{
	This data frame contains the following variables (columns):
  \describe{
    \item{\code{Nu}}{Order Number}
    \item{\code{Country}}{Country Name}
    \item{\code{year}}{1965 or 1990}
		\item{\code{Y}}{Real GDP chain, international prices of 1985}
    \item{\code{K}}{Capital stock, international prices of 1985}
    \item{\code{L}}{Number of workers, in thousands}
  }
}


\details{
	The Penn World Table was developed by Robert Summers and Alan Heston (and others) to facilitate consistent national accounts comparisons across countries as well as over time. The data can be used to evaluate the efficiency of economies of various countries in years 1965 and 1990.
}

\source{
\url{http://www.rug.nl/research/ggdc/data/pwt/pwt-5.6}. These data were originally hosted on the website of the Center for International Comparisons at the University of Pennsylvania.
}

\references{
Heston, A. and Summers, R. (1991), The Penn World Table (Mark 5): An Expanded Set of International Comparisons, 1950-1988, \emph{The Quarterly Journal of Economics}, \bold{106}, 327--368.

}

\keyword{datasets}
