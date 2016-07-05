\name{tenonradial}
\alias{tenonradial}
\title{
 Nonradial Measure of Technical Efficiency, the Russell Measure
}
\description{
Routine \code{tenonradial} uses reduced linear programming to compute the nonradial output- or input-based measure of technical efficiency, which is known as the Russell measure. In input-based nonradial efficiency measurement, this measure allows for non-proportional/different reductions in each positive input, and this is what permits it to shrink an input vector all the way back to the efficient subset. In output-based nonradial efficiency measurement, the Russell measure allows for non-proportional/different expansions of each positive output. 
}
\usage{
tenonradial(formula, data, subset,
 rts = c("C", "NI", "V"), 
 base = c("output", "input"), 
 ref = NULL, data.ref = NULL, subset.ref = NULL,
 print.level = 1)
}
\arguments{
\item{formula}{
an object of class ``formula'' (or one that can be coerced to that class): a symbolic description of the model. The details of model specification are given under `Details'.
}
\item{data}{
an optional data frame containing the variables in the model. If not found in data, the variables are taken from environment (\code{formula}), typically the environment from which \code{tenonradial} is called.
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
an optional data frame containing the variables in the technology reference set. If not found in \code{data.ref}, the variables are taken from environment(ref), typically the environment from which \code{tenonradial} is called.
}
\item{subset.ref}{
an optional vector specifying a subset of observations to define the technology reference set.
}
\item{print.level}{
numeric. 0 - nothing is printed; 1 - print summary of the model and data. 2 - print summary of technical efficiency measures. 3 - print estimation results observation by observation. Default is 1.
}
}
\details{
Routine \code{tenonradial} computes the nonradial output- or input-based measure of technical efficiency under assumption of constant, non-increasing, or variable returns to scale technology.

Models for \code{tenonradial} are specified symbolically. A typical model has the form \code{outputs ~ inputs}, where \code{outputs} (\code{inputs}) is a series of (numeric) terms which specifies outputs (inputs). The same goes for reference set.  Refer to the examples.

Results can be summarized using \code{\link{summary.npsf}}.
}

\value{
\code{tenonradial} returns a list of class \code{npsf} containing the following elements:
\item{model}{
string: model name.
}
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
\item{te}{
numeric: nonradial measure (Russell) of technical efficiency.
}
\item{esample}{
logical: returns TRUE if the observation in user supplied data is in the estimation subsample and FALSE otherwise.
}
\item{esample.ref}{
logical: returns TRUE if the observation in the user supplied reference is in the reference subsample and FALSE otherwise.
}
}

\note{
In case of one input (output), the input (output)-based Russell measure is equal to Debrue-Farrell (\code{\link{teradial}}) measure of technical efficiency.

Results can be summarized using \code{\link{summary.npsf}}.
}

\references{
Färe, R. and Lovell, C. A. K. (1994), Measuring the technical efficiency of production, \emph{Journal of Economic Theory}, \bold{19}, 150--162 

Färe, R., Grosskopf, S. and Lovell, C. A. K. (1994), \emph{Production Frontiers}, Cambridge U.K.: Cambridge University Press
}

\author{
Oleg Badunenko <obadunen@uni-koeln.de>, Pavlo Mozharovskyi <pavlo.mozharovskyi@univ-rennes1.fr>
}

\seealso{
\code{\link{teradial}}, \code{\link{teradialbc}}, \code{\link{nptestrts}}, \code{\link{nptestind}}, \code{\link{sf}}, \code{\link{summary.npsf}} for printing summary results
}
\examples{

require( npsf )

# Prepare data and matrices

data( pwt56 )
head( pwt56 )

# Create some missing values

pwt56 [49, "K"] <- NA # create missing

Y1 <- as.matrix ( pwt56[ pwt56$year == 1965, c("Y"), drop = FALSE] )
X1 <- as.matrix ( pwt56[ pwt56$year == 1965, c("K", "L"), drop = FALSE] )

X1 [51, 2] <- NA # create missing
X1 [49, 1] <- NA # create missing

data( ccr81 )
head( ccr81 )

# Create some missing values

ccr81 [64, "x4"] <- NA # create missing
ccr81 [68, "y2"] <- NA # create missing

Y2 <- as.matrix( ccr81[ , c("y1", "y2", "y3"), drop = FALSE] )
X2 <- as.matrix( ccr81[ , c("x1", "x2", "x3", "x4", "x5"), drop = FALSE] )

# Computing without reference set

# Using formula

# Country is a categorical variable, so nonradial gives error message

# t1 <- tenonradial ( Country ~ K + L, data = pwt56 )

# for computing the efficiencies of countries in 1965 
# with technology reference set is defined by observations in 1965
# (that same sample of countries)

t2 <- tenonradial ( Y ~ K + L, data = pwt56, rts = "v", 
base = "in", print.level = 2)

# Using a subset

t3 <- tenonradial ( Y ~ K + L, data = pwt56, subset = year == 1965,
	rts = "VRS", base = "in", print.level = 3 )

t4 <- tenonradial ( Y ~ K + L, data = pwt56, subset = Nu < 10,
	rts = "vrs", base = "I" )

t5 <- tenonradial ( Y ~ L, data = pwt56, subset = Nu < 10, rts = "v" )

# Multiple outputs

t8 <- tenonradial ( y1 + y2 + y3 ~ x1 + x2 + x3 + x4 + x5, data = ccr81,
	rts = "v", base = "i" )

# Using a subset

t7 <- tenonradial ( y1 + y2 + y3 ~ x1 + x2 + x3 + x4 + x5, data = ccr81,
	subset = x5 != 22, rts = "n", base = "o" )

# Computation using matrices

t9 <- tenonradial ( Y1 ~ X1, rts = "v", base = "i" )

# Define subsets on a fly

t10 <- tenonradial ( Y1[-1,] ~ X1[-2,1] )
t11 <- tenonradial ( Y1[-3,] ~ X1[-1,], rts = "v", base = "o" )

# Multiple outputs

t12 <- tenonradial ( Y2 ~ X2 )
t13 <- tenonradial ( Y2[-66,] ~ X2[-1, -c(1,3)] )


# Computing with reference set

# Using formula

# For computing the efficiencies of countries with order number
# less than 10 with technology reference set defined by countries
# with order number larger than 10 and smaller than 11 (in effect 
# no reference set, hence warning) type

t14 <- tenonradial ( Y ~ K + L, data = pwt56, subset = Nu < 10, 
	ref = Y ~ K + L, data.ref = pwt56,
	subset.ref = Nu > 10 & Nu < 11 ) # warning

# For computing the efficiencies of countries with order number
# less than 10 with technology reference set defined by countries 
# with order number larger than 10 and smaller than 15 type

t15 <- tenonradial ( Y ~ K + L, data = pwt56, subset = Nu < 10, ref = Y ~ K + L, 
	data.ref = pwt56, subset.ref = Nu > 10 & Nu < 15 )

# For computing the efficiencies of countries in 1965 
# with technology reference set is defined by observations in both
# 1965 and 1990 (all) type
	
t16 <- tenonradial ( Y ~ K + L, data = pwt56, subset = year == 1965,
	rts = "v", base = "i", 
	ref = Y ~ K + L, data.ref = pwt56 )

# For computing the efficiencies of countries in 1990
# with technology reference set is defined by observations in 1965
# type

t17 <- tenonradial ( Y ~ K + L, data = pwt56, subset = year == 1990, 
	ref = Y ~ K + L, data.ref = pwt56, subset.ref = year == 1965 )

# Using matrices

t18 <- tenonradial ( Y1[-1,] ~ X1[-2,], ref = Y1[-2,] ~ X1[-1,] )

# error: not equal number of observations in outputs and inputs

# t19 <- tenonradial ( Y1[-1,] ~ X1[-(1:2),], 
# ref = Y1[-2,] ~ X1[-1,1] )

# Combined formula and matrix

# error: not equal number of inputs in data and reference set

# t20 <- tenonradial ( Y ~ K + L, data = pwt56, subset = Nu < 10,
# ref = Y1[-2,] ~ X1[-1,1] )

t21 <- tenonradial ( Y ~ K + L, data = pwt56, subset = Nu < 10, 
	ref = Y1[-2,] ~ X1[-1,] )

\dontrun{

# Really large data-set

data(usmanuf)
head(usmanuf)

nrow(usmanuf)
table(usmanuf$year)

# This will take some time depending on computer power

t22 <- tenonradial ( Y ~ K + L + M, data = usmanuf, 
	subset = year >= 1995 & year <= 2000 ) 

# Summary

summary ( t22$te )

# Write efficiencies to the data frame:

usmanuf$te_nonrad_crs_out[ t22$esample ] <- t22$te

head(usmanuf, 17)

}


}

\keyword{ Data Envelopment Analysis }
\keyword{ Nonradial measure }
\keyword{ Nonparametric efficiency analysis }
\keyword{ Linear programming }
\keyword{ summary }