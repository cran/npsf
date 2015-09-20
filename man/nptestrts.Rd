\name{nptestrts}
\alias{nptestrts}
\title{
Nonparametric Test of Returns to Scale
}
\description{
Routine \code{nptestrts} performs nonparametric tests the returns to scale of the underlying technology via bootstrapping techniques.
}
\usage{
nptestrts(formula, data, subset,
 base = c("output", "input"),
 homogeneous = TRUE, test.two = TRUE,
 reps = 999, alpha = 0.05,
 core.count = 1, cl.type = c("SOCK", "MPI"),
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
  \item{base}{
character or numeric. string: first letter of the word ``o'' for computing output-based or ``i'' for computing input-based technical efficiency measure. string: 2 for computing output-based or 1 for computing input-based technical efficiency measure
}
  \item{homogeneous}{
logical. If TRUE, the reference set is bootstrapped with homogeneous smoothing; if FALSE, the reference set is bootstrapped with heterogeneous smoothing.
}
  \item{test.two}{
logical. If TRUE, test 2, where efficiency measures under assumption of non-increasing and variable returns to scale technology are compared; if FALSE, \code{nptestrts} stops after test 1 is completed.
}
  \item{reps}{
specifies the number of bootstrap replications to be performed.  The default is 999.  The minimum is 100.  Adequate estimates of confidence intervals using bias-corrected methods typically require 1,000 or more replications.
}
  \item{alpha}{
sets significance level; default is \code{alpha=0.05}.
}
  \item{core.count}{
positive integer. Number of cluster nodes. If \code{core.count=1}, the process runs sequentially. See  \code{\link[snowFT]{performParallel}} for more details.
}
  \item{cl.type}{
Character string that specifies cluster type (see \code{\link[snowFT]{makeClusterFT}}). Possible values are 'MPI' and 'SOCK' ('PVM' is currently not available). See  \code{\link[snowFT]{performParallel}} for more details.
}
  \item{dots}{
logical. Relevant if \code{print.level>=1}. If TRUE, one dot character is displayed for each successful replication; if FALSE,  display of the replication dots is suppressed.
}
  \item{print.level}{
numeric. 0 - nothing is printed; 1 - print summary of the model and data. 2 - print summary of technical efficiency measures. 3 - print estimation results observation by observation. Default is 1.
}
}
\details{
Routine \code{nptestrts} \code{nptestrts} performs nonparametric tests the returns to scale of the underlying technology. 

If \code{test.two} is not specified, \code{nptestrts} performs only Test #1, which consists of two parts. First, the null hypothesis that the technology is globally CRS (vs VRS) is tested. Second, the null hypothesis that the data point is scale efficient is tested.

If \code{test.two} is specified, \code{nptestrts} may perform Test #2. If the null hypothesis that the technology is CRS is rejected, \code{test.two} requests that \code{nptestrts} tests the null hypothesis that the technology is NIRS (vs VRS). If not all data points are scale efficient, \code{nptestrts} tests that the reason for scale inefficiency is DRS. If the null hypothesis that the technology is CRS is not rejected and all data points are scale efficient, \code{nptestrts} will not perform Test #2 even if \code{test.two} is specified.

Models for \code{nptestrts} are specified symbolically. A typical model has the form \code{outputs ~ inputs}, where \code{outputs} (\code{inputs}) is a series of (numeric) terms which specifies outputs (inputs).  Refer to the examples.

If \code{core.count>=1}, \code{nptestrts} will perform bootstrap on multiple cores.  Parallel computing requires package \code{snowFT}. By the default cluster type is defined by option \code{cl.type="SOCK"}.  Specifying \code{cl.type="MPI"} requires package \code{Rmpi}.

On some systems, specifying option \code{cl.type="SOCK"} results in much quicker execution than specifying option \code{cl.type="MPI"}.  Option \code{cl.type="SOCK"} might be problematic on Mac system.

Parallel computing make a difference for large data sets.  Specifying option \code{dots=TRUE} will indicate at what speed the bootstrap actually proceeds.  Specify \code{reps=100} and compare two runs with option \code{core.count=1} and \code{core.count>1} to see if parallel computing speeds up the bootstrap.  For small samples, parallel computing may actually slow down the \code{nptestrts}.

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
  \item{teCrs}{
numeric: measures of technical efficiency under the assumption of CRS.
}
  \item{teNrs}{
numeric: measures of technical efficiency under the assumption of NiRS.
}
  \item{teVrs}{
numeric: measures of technical efficiency under the assumption of VRS.
}
  \item{sefficiency}{
numeric: scale efficiency.
}
  \item{sefficiencyMean}{
numeric: ratio of means of technical efficiency measures under CRS and VRS.
}
  \item{pGlobalCRS}{
numeric: p-value of the test that the technology is globally CRS.
}

  \item{psefficient}{
numeric: p-value of the test that data point is statistically scale efficient.
}
  \item{sefficient}{
logical: returns \code{TRUE}, if statistically scale efficient; \code{FALSE} otherwise.
}
  \item{nsefficient}{
numeric: number of statistically scale efficient.
}

  \item{nrsOVERvrsMean}{
numeric: ratio of means of technical efficiency measures under NIRS and VRS (if \code{test.two=TRUE}).
}
  \item{pGlobalNRS}{
numeric: p-value of the test the technology is globally NIRS (if \code{test.two=TRUE}).
}

  \item{sineffdrs}{
logical: returns \code{TRUE} if statistically scale inefficient due to DRS and \code{FALSE} if statistically scale inefficient due to IRS (if \code{test.two=TRUE} and not all data points are statistically scale efficient \code{nsefficient<K}).
}
  \item{pineffdrs}{
numeric: p-value of the test that data point is scale inefficient due to DRS (if \code{test.two=TRUE} and not all data points are statistically scale efficient \code{nsefficient<K}).
}
  \item{nrsOVERvrs}{
numeric: ratio of measures of technical efficiency under NiRS and VRS (if \code{test.two=TRUE} and not all data points are statistically scale efficient \code{nsefficient<K}).
}

  \item{esample}{
logical: returns TRUE if the observation in user supplied data is in the estimation subsample and FALSE otherwise.
}

}

\note{
Before specifying option \code{homogeneous} it is advised to preform the test of independence (see \code{\link{nptestind}}).
}

\references{
Färe, R. and Lovell, C. A. K. (1994), Measuring the technical efficiency of production, \emph{Journal of Economic Theory}, \bold{19}, 150--162 

Färe, R., Grosskopf, S. and Lovell, C. A. K. (1994), \emph{Production Frontiers}, Cambridge U.K.: Cambridge University Press

Simar L. and P.W. Wilson (2002), Nonparametric Tests of Return to Scale, \emph{European Journal of Operational Research}, \bold{139}, 115--132
}

\author{
Oleg Badunenko <obadunen@uni-koeln.de>, Pavlo Mozharovskyi <mozharovskyi@statistik.uni-koeln.de>

}

\seealso{
\code{\link{teradial}}, \code{\link{tenonradial}}, \code{\link{teradialbc}}, \code{\link{nptestind}}, \code{\link{sf}}
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

# Perform output-based test of returns to scale smoothed 
# homogeneous bootstrap with 999 replications at the 5% 
# significance level.  Also perform Test #2

t1 <- nptestrts(y1 + y2 + y3 ~ x1 + x2 + x3 + x4 + x5,
	data = ccr81, homogeneous = TRUE,
	reps = 999, dots = TRUE, base = "o")

# suppress printing replication dots
t2 <- nptestrts(Y2 ~ X2,
	homogeneous = TRUE,
	reps = 100, dots = FALSE, base = "o")


# heterogeneous
t3 <- nptestrts(Y2 ~ X2,
	homogeneous = FALSE,
	reps = 100, dots = TRUE, base = "o")


# ===========================
# ===  Parallel computing ===
# ===========================

# Perform previous test but use 4 cores and
# cluster type MPI

t3 <-  nptestrts(y1 + y2 + y3 ~ x1 + x2 + x3 + x4 + x5,
	data = ccr81, homogeneous = FALSE,
	reps = 100, dots = TRUE, base = "o",
	core.count = 4, cl.type = "MPI")


# Really large data-set

data(usmanuf)
head(usmanuf)

nrow(usmanuf)
table(usmanuf$year)

# Figure industries to include in the sample (first quarter)

summary(usmanuf[usmanuf$year >= 1999 & usmanuf$year < 2000, "naics"])

# This test is quite demanding and it will take some time
# depending on computer power

t4 <- nptestrts(Y ~ K + L + M, data = usmanuf,
	subset = year >= 1999 & year < 2000 & naics < 321900,
	homogeneous = FALSE, reps = 100, dots = TRUE, base = "o",
	core.count = 4, cl.type = "MPI")

# This is very computer intensive task

t5 <- nptestrts(Y ~ K + L + M, data = usmanuf,
	subset = year >= 1999 & year < 2000,
	homogeneous = FALSE, reps = 100, dots = TRUE, base = "o",
	core.count = 4, cl.type = "MPI")

}

}

\keyword{ Data Envelopment Analysis }
\keyword{ Radial measure }
\keyword{ Nonparametric efficiency analysis }
\keyword{ Linear programming }
\keyword{ Nonparametric test }
\keyword{ Bootstrap }
\keyword{ Smoothed bootstrap }
\keyword{ Subsampling bootstrap }
\keyword{ Parallel computing }
