\name{usmanuf}
\docType{data}
\alias{usmanuf}

\title{US Manufacturing Industry Data}

\description{
 This data come from the National Bureau of Economic Research Center for Economic Studies manufacturing industry database.  This data set provides only selected variables.
}

\usage{data( usmanuf )}

\format{
	This data frame contains the following variables (columns):
  \describe{
    \item{\code{naics}}{NAICS 6-digit Codes}
    \item{\code{year}}{Year ranges from 1990 to 2009}
    \item{\code{Y}}{Total value of shipments in $1m divided by the deflator for VSHIP 1997=1.000 (vship/piship)}
    \item{\code{K}}{Total real capital stock in $1m (cap)}
    \item{\code{L}}{Total employment in thousands (emp)}
		\item{\code{M}}{Total cost of materials in $1m divided by the deflator for MATCOST 1997=1.000 plus oost of electric & fuels in $1m divided by the deflator for ENERGY 1997=1.000 (matcost/pimat +  energy/pien)}
  }
}


\details{
These data come from the National Bureau of Economic Research Center for Economic Studies manufacturing industry database, which was compiled by Randy A. Becker and Wayne B. Gray. This database is a joint effort between the National Bureau of Economic Research (NBER) and U.S. Census Bureau's Center for Economic Studies (CES), containing annual industry-level data from 1958-2009 on output, employment, payroll and other input costs, investment, capital stocks, TFP, and various industry-specific price indexes. Because of the change from SIC to NAICS industry definitions in 1997, the database is provided in two versions: one with 459 four-digit 1987 SIC industries and the other with 473 six-digit 1997 NAICS industries.
}

\source{
\url{https://www.nber.org/nberces/}.
}

\references{
Bartelsman, E.J. and Gray, W. (1996), The NBER Manufacturing Productivity Database, \emph{National Bureau of Economic Research}, Technical Working Paper Series, \doi{10.3386/t0205}

}


\keyword{datasets}
