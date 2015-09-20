\name{nptestind}
\alias{nptestind}
\title{
Nonparametric Test of Independence
}
\description{
In output based efficiency measurement, routine \code{nptestind} perform test that radial (Debreu-Farrell) output-based measure of technical efficiency under chosen assumption about the technology and mix of outputs are independent. In input-based efficiency measurement, routine \code{nptestind} perform test that radial (Debreu-Farrell) input-based measure of technical efficiency under chosen assumption about the technology and mix of inputs are independent. Testing is performed using bootstrap technique.
}
\usage{
nptestind(formula, data, subset,
 rts = c("C", "NI", "V"), base = c("output", "input"),
 reps = 999, alpha = 0.05,
 print.level = 1, dots = TRUE)
}
\arguments{
  \item{formula}{
an object of class ``formula'' (or one that can be coerced to that class): a symbolic description of the model. The details of model specification are given under `Details'.
}
  \item{data}{
an optional data frame containing the variables in the model. If not found in data, the variables are taken from environment (\code{formula}), typically the environment from which \code{teradial} is called.
}
  \item{subset}{
an optional vector specifying a subset of observations for which technical efficiency is to be computed.
}
  \item{rts}{
character or numeric. string: first letter of the word ``c'' for constant, ``n'' for non-increasing, or ``v'' for variable returns to scale assumption. numeric: 3 for constant, 2 for non-increasing, or 3 for variable returns to scale assumption.  
}
  \item{base}{
character or numeric. string: first letter of the word ``o'' for computing output-based or ``i'' for computing input-based technical efficiency measure. string: 2 for computing output-based or 1 for computing input-based technical efficiency measure
}
  \item{reps}{
specifies the number of bootstrap replications to be performed.  The default is 999.  The minimum is 100.  Adequate estimates of confidence intervals using bias-corrected methods typically require 1,000 or more replications.
}
  \item{alpha}{
sets significance level; default is \code{alpha=0.05}.
}
  \item{dots}{
logical. Relevant if \code{print.level>=1}. If TRUE, one dot character is displayed for each successful replication; if FALSE,  display of the replication dots is suppressed.
}
  \item{print.level}{
numeric. 0 - nothing is printed; 1 - print summary of the model and data. 2 - print summary of technical efficiency measures. 3 - print estimation results observation by observation. Default is 1.
}
}
\details{
In output based efficiency measurement, routine \code{nptestind} perform test that radial (Debreu-Farrell) output-based measure of technical efficiency under chosen assumption about the technology and mix of outputs are independent. In input-based efficiency measurement, routine \code{nptestind} perform test that radial (Debreu-Farrell) input-based measure of technical efficiency under chosen assumption about the technology and mix of inputs are independent.

Testing is performed using bootstrap technique.
}

\value{
\code{nptestrts} returns a list of class \code{npsf} containing the following elements:
  \item{K}{
numeric: number of data points.
}
  \item{M}{
numeric: number of inputs.
}
  \item{N}{
numeric: number of outputs.
}
  \item{rts}{
string: RTS assumption.
}
  \item{base}{
string: base for efficiency measurement.
}
  \item{reps}{
numeric: number of bootstrap replications.
}
  \item{alpha}{
numeric: significance level.
}
  \item{t4n}{
numeric: value of the T4n statistic.
}

  \item{pval}{
numeric: p-value of the test of independence.
}
}

\references{
Färe, R. and Lovell, C. A. K. (1994), Measuring the technical efficiency of production, \emph{Journal of Economic Theory}, \bold{19}, 150--162 

Färe, R., Grosskopf, S. and Lovell, C. A. K. (1994), \emph{Production Frontiers}, Cambridge U.K.: Cambridge University Press

Wilson P.W.  (2003), Testing Independence in Models of Productive Efficiency, \emph{Journal of Productivity Analysis}, \bold{20}, 361--390
}

\author{
Oleg Badunenko <obadunen@uni-koeln.de>, Pavlo Mozharovskyi <mozharovskyi@statistik.uni-koeln.de>

}

\seealso{
\code{\link{teradial}}, \code{\link{tenonradial}}, \code{\link{teradialbc}}, \code{\link{nptestrts}}, \code{\link{sf}}
}

\examples{

\dontrun{

require( npsf )

# Prepare data and matrices

data( ccr81 )
head( ccr81 )

# Create some missing values

ccr81 [64, "x4"] <- NA # just to create missing
ccr81 [68, "y2"] <- NA # just to create missing

Y2 <- as.matrix( ccr81[ , c("y1", "y2", "y3"), drop = FALSE] )
X2 <- as.matrix( ccr81[ , c("x1", "x2", "x3", "x4", "x5"), drop = FALSE] )

# Perform nonparametric test that radial (Debreu-Farrell) 
# output-based measure of technical efficiency under assumption of 
# NIRS technology and mix of outputs are independent. Test is 
# performed based on 999 replications at the 5% significance level.

t1 <- nptestind ( y1 + y2 + y3 ~ x1 + x2 + x3 + x4 + x5, 
	data = ccr81, base = "o", rts = "n", 
	reps = 999, dots = TRUE)


# Really large data-set

data(usmanuf)
head(usmanuf)

nrow(usmanuf)
table(usmanuf$year)

# This will take some time depending on computer power

data(usmanuf)
head(usmanuf)

t2 <- nptestind ( Y ~ K + L + M, data = usmanuf, 
	subset = year >= 1999 & year <= 2000, 
	reps = 999, dots = TRUE, base = "i", rts = "v")

}

}

\keyword{ Data Envelopment Analysis }
\keyword{ Radial measure }
\keyword{ Nonparametric efficiency analysis }
\keyword{ Linear programming }
\keyword{ Nonparametric test }
\keyword{ Bootstrap }

