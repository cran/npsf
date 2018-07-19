\name{npsf-package}
\alias{npsf-package}
\alias{npsf}
\docType{package}
\title{
 Introduction to Nonparametric and Stochastic Frontier Analysis
}

\description{
This package provides a variety of tools for nonparametric and parametric efficiency measurement. 
% A vignette containing many of the examples found in the help files accompanying the np package that is intended to serve as a gentle introduction to this package can be accessed via vignette("np", package="np").
}

\details{
The nonparametric models in \code{npsf} comprise nonradial efficiency measurement (\code{\link{tenonradial}}), where non-proportional reductions (expansions) in each positive input (output) are allowed, as well as popular radial efficiency measurement (\code{\link{teradial}}), where movements to the frontier are proportional.

Using bootstrapping techniques,  \code{\link{teradialbc}}, \code{\link{nptestrts}}, \code{\link{nptestind}} deal with statistical inference about the radial efficiency measurement. \code{\link{nptestind}} helps in deciding which type of the bootstrap to employ. Global return to scale and individual scale efficiency is tested by \code{\link{nptestrts}}.  Finally, \code{\link{teradialbc}} performs bias correction of the radial Debrue-Farrell input- or output-based measure of technical efficiency, computes bias and constructs confidence intervals.

Computer intensive functions \code{\link{teradialbc}} and \code{\link{nptestrts}} allow making use of parallel computing, even on a single machine with multiple cores. Help files contain examples that are intended to introduce the usage.

We use own LP solver for solving linear programming problems.

The parametric stochastic frontier models in \code{npsf} can be estimated by \code{\link{sf}}, which performs maximum likelihood estimation of the frontier parameters and technical or cost efficiencies.  Inefficiency error component can be assumed to be have either half-normal or truncated normal distribution.  \code{\link{sf}} allows modelling multiplicative heteroskedasticity of either inefficiency or random noise component, or both.  Additionally, marginal effects of exogenous variable(s) on the expected value of inefficiency term can be computed.

For details of the respective method please see the reference at the end of this introduction and of the respective help file.

All function in \code{npsf} accept formula with either names of variables in the data set, or names of the matrices.  Except for \code{\link{nptestind}}, all function return \code{esample}, a logical vector length of which is determined by \code{data} and \code{subset} (if specified) or number of rows in matrix \code{outputs}. \code{esample} equals \code{TRUE} if this data point parted in estimation procedure, and \code{FALSE} otherwise.

Results can be summarized using \code{\link{summary.npsf}}.

}

\author{
Oleg Badunenko, <oleg.badunenko@port.ac.uk>

Yaryna Kolomiytseva, <kolomiytseva@wiso.uni-koeln.de>

Pavlo Mozharovskyi, <pavlo.mozharovskyi@ensai.fr>

Maintainer: Oleg Badunenko <oleg.badunenko@port.ac.uk>
}

\references{

Bartelsman, E.J. and Gray, W. (1996), The NBER Manufacturing Productivity Database, \emph{National Bureau of Economic Research}, Technical Working Paper Series, \url{http://www.nber.org/papers/t0205}.

Battese, G., Coelli, T. (1988), Prediction of firm-level technical effiiencies with a generalized frontier production function and panel data. \emph{Journal of Econometrics}, \bold{38}, 387--399.

Battese, G., Coelli, T. (1992), Frontier production functions, technical efficiency and panel data: With application to paddy farmers in India. \emph{Journal of Productivity Analysis}, \bold{3}, 153--169.

Charnes, A., W. W. Cooper, and E. Rhodes. 1981. Evaluating Program and Managerial Efficiency: An Application of Data Envelopment Analysis to Program Follow Through. \emph{Management Science} \bold{27}: 668--697.

Caudill, S., Ford, J., Gropper, D. (1995), Frontier estimation and firm-specific inefficiency measures in the presence of heteroscedasticity. \emph{Journal of Business and Economic Statistics}, \bold{13}, 105--111.

Debreu, G. 1951. The Coefficient of Resource Utilization. \emph{Econometrica} \bold{19}: 273--292.

Färe, R. and Lovell, C. A. K. (1994), Measuring the technical efficiency of production, \emph{Journal of Economic Theory}, \bold{19}, 150--162.

Färe, R., Grosskopf, S. and Lovell, C. A. K. (1994), \emph{Production Frontiers}, Cambridge U.K.: Cambridge University Press.

Farrell, M. J. 1957. The Measurement of Productive Efficiency. \emph{Journal of the Royal Statistical Society. Series A (General)} \bold{120}(3): 253--290.

Heston, A., and R. Summers. 1991. The Penn World Table (Mark 5): An Expanded Set of International Comparisons, 1950-1988. \emph{The Quarterly Journal of Economics} \bold{106}: 327--368.

Horrace, W. and Schmidt, P. (1996), On ranking and selection from independent truncated normal distributions. \emph{Journal of Productivity Analysis}, \bold{7}, 257--282.

Jondrow, J., Lovell, C., Materov, I., Schmidt, P. (1982), On estimation of technical inefficiency in the stochastic frontier production function model. \emph{Journal of Econometrics}, \bold{19}, 233--238.

Kneip, A., Simar L., and P.W. Wilson (2008), Asymptotics and Consistent Bootstraps for DEA Estimators in Nonparametric Frontier Models, \emph{Econometric Theory}, \bold{24}, 1663--1697

Koetter, M., Kolari, J., and Spierdijk, L. (2012), Enjoying the quiet life under deregulation? Evidence from adjusted Lerner indices for U.S. banks. \emph{Review of Economics and Statistics}, \bold{94}, 2, 462--480.

Kumbhakar, S. (1990),  Production Frontiers, Panel Data, and Time-varying Technical Inefficiency. \emph{Journal of Econometrics}, \bold{46}, 201--211.

Kumbhakar, S. and Lovell, C. (2003), \emph{Stochastic Frontier Analysis.} Cambridge: Cambridge University Press.

Restrepo-Tobon, D. and Kumbhakar, S. (2014), Enjoying the quiet life under deregulation? Not Quite. \emph{Journal of Applied Econometrics}, \bold{29}, 2, 333--343.

Simar, L. and P.W. Wilson (1998), Sensitivity Analysis of Efficiency Scores: How to Bootstrap in Nonparametric Frontier Models, \emph{Management Science}, \bold{44}, 49--61

Simar, L. and P.W. Wilson (2000), A General Methodology for Bootstrapping in Nonparametric Frontier Models, \emph{Journal of Applied Statistics}, \bold{27}, 779--802

Simar, L. and P.W. Wilson (2002), Nonparametric Tests of Return to Scale, \emph{European Journal of Operational Research}, \bold{139}, 115--132

Wang, H.-J. (2002), Heteroskedasticity and non-monotonic efficiency effects of a stochastic frontier model. \emph{Journal of Productivity Analysis}, \bold{18}, 241--253.

Wilson P.W.  (2003), Testing Independence in Models of Productive Efficiency, \emph{Journal of Productivity Analysis}, \bold{20}, 361--390

}

\keyword{ Stochastic Frontier Analysis }
\keyword{ Panel Data Models }
\keyword{ Panel Data Analysis }
\keyword{ Data Envelopment Analysis }
\keyword{ Radial measure }
\keyword{ Nonradial measure }
\keyword{ Nonparametric efficiency analysis }
\keyword{ Linear programming }
\keyword{ Bootstrap }
\keyword{ Bias-correction }
\keyword{ Smoothed bootstrap }
\keyword{ Subsampling bootstrap }
\keyword{ Parallel computing }
\keyword{ Variance-covariance matrix }
\keyword{ summary }