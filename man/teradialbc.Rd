\name{teradialbc}
\alias{teradialbc}
\title{
 Statistical Inference Regarding the Radial Measure of Technical Efficiency
}
\description{
Routine \code{teradialbc} performs bias correction of the radial Debrue-Farrell input- or output-based measure of technical efficiency, computes bias and constructs confidence intervals via bootstrapping techniques.
}
\usage{
teradialbc(formula, data, subset,
 ref = NULL, data.ref = NULL, subset.ref = NULL,
 rts = c("C", "NI", "V"), base = c("output", "input"),
 homogeneous = TRUE, smoothed = TRUE, kappa = NULL,
 reps = 999, level = 95,
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
\item{rts}{
character or numeric. string: first letter of the word ``c'' for constant, ``n'' for non-increasing, or ``v'' for variable returns to scale assumption. numeric: 3 for constant, 2 for non-increasing, or 1 for variable returns to scale assumption. 
}
\item{base}{
character or numeric. string: first letter of the word ``o'' for computing output-based or ``i'' for computing input-based technical efficiency measure. string: 2 for computing output-based or 1 for computing input-based technical efficiency measure
}
\item{ref}{
an object of class ``formula'' (or one that can be coerced to that class): a symbolic description of inputs and outputs that are used to define the technology reference set. The details of technology reference set specification are given under `Details'. If reference is not provided, the technical efficiency measures for data points are computed relative to technology based on data points themselves.
}
\item{data.ref}{
an optional data frame containing the variables in the technology reference set. If not found in \code{data.ref}, the variables are taken from environment(ref), typically the environment from which \code{teradial} is called.
}
\item{subset.ref}{
an optional vector specifying a subset of observations to define the technology reference set.
}
\item{smoothed}{
logical. If TRUE, the reference set is bootstrapped with smoothing; if FALSE, the reference set is bootstrapped with subsampling.
}
\item{homogeneous}{
logical. Relevant if \code{smoothed=TRUE}. If TRUE, the reference set is bootstrapped with homogeneous smoothing; if FALSE, the reference set is bootstrapped with heterogeneous smoothing.
}
\item{kappa}{
relevant if \code{smoothed=TRUE}. 'kappa' sets the size of the subsample as K^kappa, where K is the number of data points in the original reference set. The default value is 0.7. 'kappa' may be between 0.5 and 1.
}
\item{reps}{
specifies the number of bootstrap replications to be performed.  The default is 999.  The minimum is 100.  Adequate estimates of confidence intervals using bias-corrected methods typically require 1,000 or more replications.
}
\item{level}{
sets confidence level for confidence intervals; default is \code{level = 95}.
}
\item{core.count}{
positive integer. Number of cluster nodes. If \code{core.count=1}, the process runs sequentially. See  \code{performParallel} in package \code{snowFT} for more details.
}
\item{cl.type}{
Character string that specifies cluster type (see \code{makeClusterFT} in package \code{snowFT}). Possible values are 'MPI' and 'SOCK' ('PVM' is currently not available). See \code{performParallel} in package \code{snowFT} for more details.
}
\item{dots}{
logical. Relevant if \code{print.level>=1}. If TRUE, one dot character is displayed for each successful replication; if FALSE,  display of the replication dots is suppressed.
}
\item{print.level}{
numeric. 0 - nothing is printed; 1 - print summary of the model and data. 2 - print summary of technical efficiency measures. 3 - print estimation results observation by observation. Default is 1.
}
}
\details{
Routine \code{teradialbc} performs bias correction of the radial Debrue-Farrell input- or output-based measure of technical efficiency, computes bias and constructs confidence intervals via bootstrapping techniques. See Simar and Wilson (1998) \doi{10.1287/mnsc.44.1.49}, Simar and Wilson (2000) \doi{10.1080/02664760050081951}, Kneip, Simar, and Wilson (2008) \doi{10.1017/S0266466608080651}, and references with links below.

Models for \code{teradialbc} are specified symbolically. A typical model has the form \code{outputs ~ inputs}, where \code{outputs} (\code{inputs}) is a series of (numeric) terms which specifies outputs (inputs). The same goes for reference set.  Refer to the examples.

If \code{core.count>=1}, \code{teradialbc} will perform bootstrap on multiple cores.  Parallel computing requires package \code{snowFT}. By the default cluster type is defined by option \code{cl.type="SOCK"}.  Specifying \code{cl.type="MPI"} requires package \code{Rmpi}.

On some systems, specifying option \code{cl.type="SOCK"} results in much quicker execution than specifying option \code{cl.type="MPI"}.  Option \code{cl.type="SOCK"} might be problematic on Mac system.

Parallel computing make a difference for large data sets.  Specifying option \code{dots=TRUE} will indicate at what speed the bootstrap actually proceeds.  Specify \code{reps=100} and compare two runs with option \code{core.count=1} and \code{core.count>1} to see if parallel computing speeds up the bootstrap.  For small samples, parallel computing may actually slow down the \code{teradialbc}.

Results can be summarized using \code{\link{summary.npsf}}.

}

\value{
\code{teradialbc} returns a list of class \code{npsf} containing the following elements:
\item{K}{
numeric: number of data points.
}
\item{M}{
numeric: number of outputs.
}
\item{N}{
numeric: number of inputs.
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
\item{level}{
numeric: confidence level for confidence intervals.
}
\item{te}{
numeric: radial measure (Russell) of technical efficiency.
}
\item{tebc}{
numeric: bias-corrected radial measures of technical efficiency.
}
\item{biasboot}{
numeric: bootstrap bias estimate for the original radial measures of technical efficiency.
}
\item{varboot}{
numeric: bootstrap variance estimate for the radial measures of technical efficiency.
}
\item{biassqvar}{
numeric:  one-third of the ratio  of bias squared to variance for radial measures of technical efficiency.
}
\item{realreps}{
numeric: actual number of replications used for statistical inference.
}
\item{telow}{
numeric: lower bound estimate for radial measures of technical efficiency.
}
\item{teupp}{
numeric: upper bound estimate for radial measures of technical efficiency.
}
\item{teboot}{
numeric: \code{reps x K} matrix containing bootstrapped measures of technical efficiency from each of \code{reps} bootstrap replications.
}
\item{esample}{
logical: returns TRUE if the observation in user supplied data is in the estimation subsample and FALSE otherwise.
}

}

\note{
Before specifying option \code{homogeneous} it is advised to preform the test of independence (see \code{\link{nptestind}}). Routine \code{\link{nptestrts}} may help deciding regarding option \code{rts}.

Results can be summarized using \code{\link{summary.npsf}}.
}

\references{
Färe, R. and Lovell, C. A. K. (1994), Measuring the technical efficiency of production, \emph{Journal of Economic Theory}, \bold{19}, 150--162 

Färe, R., Grosskopf, S. and Lovell, C. A. K. (1994), \emph{Production Frontiers}, Cambridge U.K.: Cambridge University Press

Kneip, A., Simar L., and P.W. Wilson (2008), Asymptotics and Consistent Bootstraps for DEA Estimators in Nonparametric Frontier Models, \emph{Econometric Theory}, \bold{24}, 1663--1697, \doi{10.1017/S0266466608080651}

Simar, L. and P.W. Wilson (1998), Sensitivity Analysis of Efficiency Scores: How to Bootstrap in Nonparametric Frontier Models, \emph{Management Science}, \bold{44}, 49--61, \doi{10.1287/mnsc.44.1.49}

Simar, L. and P.W. Wilson (2000), A General Methodology for Bootstrapping in Nonparametric Frontier Models, \emph{Journal of Applied Statistics}, \bold{27}, 779--802, \doi{10.1080/02664760050081951}
}

\author{
Oleg Badunenko <oleg.badunenko@brunel.ac.uk>, Pavlo Mozharovskyi <pavlo.mozharovskyi@telecom-paris.fr>

}

\seealso{
\code{\link{teradial}}, \code{\link{tenonradial}}, \code{\link{tenonradialbc}}, \code{\link{nptestrts}}, \code{\link{nptestind}}, \code{\link{sf}}
}

\examples{

\dontrun{

require( npsf )

# Prepare data and matrices

data( pwt56 )
head( pwt56 )

# Create some missing values

pwt56 [49, "K"] <- NA # just to create missing

Y1 <- as.matrix ( pwt56[ pwt56$year == 1965, c("Y"), drop = FALSE] )
X1 <- as.matrix ( pwt56[ pwt56$year == 1965, c("K", "L"), drop = FALSE] )

X1 [51, 2] <- NA # just to create missing
X1 [49, 1] <- NA # just to create missing

data( ccr81 )
head( ccr81 )

# Create some missing values

ccr81 [64, "x4"] <- NA # just to create missing
ccr81 [68, "y2"] <- NA # just to create missing

Y2 <- as.matrix( ccr81[ , c("y1", "y2", "y3"), drop = FALSE] )
X2 <- as.matrix( ccr81[ , c("x1", "x2", "x3", "x4", "x5"), drop = FALSE] )

# Compute output-based measures of technical efficiency under 
# the assumption of CRS (the default) and perform bias-correctiion
# using smoothed homogeneous bootstrap (the default) with 999
# replications (the default).

t1 <- teradialbc ( y1 + y2 + y3 ~ x1 + x2 + x3 + x4 + x5, 
	data = ccr81)

# or just

t2 <- teradialbc ( Y2 ~ X2)

# Combined formula and matrix

t3 <- teradialbc ( Y ~ K + L, data = pwt56, subset = Nu < 10, 
	ref = Y1[-2,] ~ X1[-1,] )

# Compute input-based measures of technical efficiency under 
# the assumption of VRS and perform bias-correctiion using
# subsampling heterogenous bootstrap with 1999 replications.
# Choose to report 99% confidence interval. The reference set
# formed by data points where x5 is not equal 10. 
# Suppress printing dots.

t4 <- teradialbc ( y1 + y2 + y3 ~ x1 + x2 + x3 + x4 + x5, 
	data = ccr81, ref = y1 + y2 + y3 ~ x1 + x2 + x3 + x4 + x5, 
	subset.ref = x5 != 10, data.ref = ccr81, reps = 1999, 
	smoothed = FALSE, kappa = 0.7, dots = FALSE, 
	base = "i", rts = "v", level = 99)

# Compute input-based measures of technical efficiency under
# the assumption of NRS and perform bias-correctiion using 
# smoothed heterogenous bootstrap with 499 replications for 
# all data points. The reference set formed by data points 
# where x5 is not equal 10.

t5 <- teradialbc ( y1 + y2 + y3 ~ x1 + x2 + x3 + x4 + x5, 
	data = ccr81, ref = y1 + y2 + y3 ~ x1 + x2 + x3 + x4 + x5, 
	subset.ref = x5 != 10, data.ref = ccr81, homogeneous = FALSE, 
	reps = 999, smoothed = TRUE, dots = TRUE, base = "i", rts = "n")


# ===========================
# ===  Parallel computing ===
# ===========================

# Perform previous bias-correction but use 8 cores and 
# cluster type SOCK

t51 <-  teradialbc ( y1 + y2 + y3 ~ x1 + x2 + x3 + x4 + x5, 
	data = ccr81, ref = y1 + y2 + y3 ~ x1 + x2 + x3 + x4 + x5, 
	subset.ref = x5 != 10, data.ref = ccr81, homogeneous = FALSE, 
	reps = 999, smoothed = TRUE, dots = TRUE, base = "i", rts = "n", 
	core.count = 8, cl.type = "SOCK")


# Really large data-set

data(usmanuf)
head(usmanuf)

nrow(usmanuf)
table(usmanuf$year)

# This will take some time depending on computer power

data(usmanuf)
head(usmanuf)

t6 <- teradialbc ( Y ~ K + L + M, data = usmanuf, 
	subset = year >= 1999 & year <= 2000, homogeneous = FALSE, 
	base = "o", reps = 100, 
	core.count = 8, cl.type = "SOCK")

}

}

\keyword{ Data Envelopment Analysis }
\keyword{ Radial measure }
\keyword{ Nonparametric efficiency analysis }
\keyword{ Linear programming }
\keyword{ Bias-correction }
\keyword{ Bootstrap }
\keyword{ Smoothed bootstrap }
\keyword{ Subsampling bootstrap }
\keyword{ Parallel computing }
\keyword{ summary }