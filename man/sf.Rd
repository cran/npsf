\name{sf}
\alias{sf}
\title{
 Stochastic Frontier Models Using Cross-Sectional Data
}
\description{
 \code{sf} performs maximum likelihood estimation of the parameters and technical or cost efficiencies in cross-sectional stochastic (production or cost) frontier models with half-normal or truncated normal distributional assumption imposed on inefficiency error component. 
}
\usage{
sf(formula, uhet = NULL, vhet = NULL,
 tmean = NULL, prod = TRUE, data, subset, 
 distribution = c("h", "t"), start.val = NULL,
 lmtol = sqrt(.Machine$double.eps),
 level = 95, marg.eff = FALSE, digits = 4, 
 print.level = 2)
}
\arguments{
 \item{formula}{
  an object of class ``formula'' (or one that can be coerced to that class): a symbolic description of the model. The details of model specification are given under `Details'.
 }
 \item{uhet}{
 one-sided formula; e.g. \code{uhet ~ z1 + z2}. Specifies exogenous variables entering the expression for variance of inefficiency error component. If \code{NULL}, inefficiency term is assumed to be homoskedastic, i.e. \eqn{\sigma_u^2 = exp(\gamma[0])}.
 }
 \item{vhet}{
 one-sided formula; e.g. \code{vhet ~ z1 + z2}. Specifies exogenous variables entering the expression for variance of random noise error component. If \code{NULL}, random noise component is assumed to be homoskedastic, i.e. \eqn{\sigma_v^2 = exp(\gamma[0])}.
 }
 \item{tmean}{
 one-sided formula; e.g. \code{tmean ~ z1 + z2}. Specifies whether the mean of pre-truncated normal distribution of inefficiency term is a linear function of exogenous variables. Used only when \code{distribution = "t"}. If \code{NULL}, mean is assumed to be constant. 
 }
 \item{prod}{
 logical. If \code{TRUE}, the estimates of parameters of stochastic production frontier model and of technical efficiencies are returned; if \code{FALSE}, the estimates of parameters of stochastic cost frontier model and of cost efficiencies are returned.
 }
  \item{data}{
  an optional data frame containing variables in the model. If not found in data, the variables are taken from environment (\code{formula}), typically the environment from which \code{sf} is called.
  }
  \item{subset}{
  an optional vector specifying a subset of observations for which technical or cost efficiencies are to be computed.
  }
 \item{distribution}{
 either \code{"h"} (half-normal) or \code{"t"} (truncated normal), specifying the distribution of inefficiency term.
 }
 \item{start.val}{
 numeric. Starting values to be supplied to the optimization routine. If \code{NULL}, OLS and method of moments estimates are used (see Kumbhakar and Lovell 2000).
 }
 \item{lmtol}{
 numeric. Tolerance for the scaled gradient in ML optimization. Default is sqrt(.Machine$double.eps).
 }
  
 \item{level}{
 numeric. Defines level\% two-sided prediction interval for technical or cost efficiencies (see Horrace and Schmidt 1996). Default is 95.
 }
 \item{marg.eff}{
 logical. If \code{TRUE}, unit-specific marginal effects of exogenous variables on the mean of distribution of inefficiency term are returned.
 }
 \item{digits}{
 numeric. Number of digits to be displayed in estimation results and for efficiency estimates. Default is 4.
 }
 \item{print.level}{
  numeric. 1 - print estimation results. 2 - print optimization details. 3 - print summary of point estimates of technical or cost efficiencies. 4 - print unit-specific point and interval estimates of technical or cost efficiencies. Default is 2.
 }
  }
  
  
\details{
 Models for \code{sf} are specified symbolically. A typical model has the form \code{y ~ x1 + ...}, where \code{y} represents the logarithm of outputs or total costs and \code{\{x1,...\}} is a series of inputs or outputs and input prices (in logs).  

Options \code{uhet} and \code{vhet} can be used if multiplicative heteroskedasticity of either inefficiency or random noise component (or both) is assumed; i.e. if their variances can be expressed as exponential functions of (e.g. size-related) exogenous variables (including intercept) (see Caudill et al. 1995).

If \code{marg.eff = TRUE} and \code{distribution = "h"}, the marginal effect of \emph{k}th exogenous variable on the expected value of inefficiency term of unit \emph{i} is computed as: \eqn{\gamma[k]\sigma[i]/\sqrt2\pi}, where \eqn{\sigma_u[i] = \sqrt exp(z[i]'\gamma)}. If \code{distribution = "t"}, marginal effects are returned if either \code{tmean} or \code{uhet} are not \code{NULL}. If the same exogenous variables are specified under both options, (non-monotonic) marginal effects are computed as explained in Wang (2002).

}

\value{
 \code{sf} returns a list of class \code{npsf} containing the following elements:
  \item{coef}{
   numeric. Named vector of ML parameter estimates.
  }
 \item{vcov}{
  matrix. Estimated covariance matrix of ML estimator. 
 }
 \item{loglik}{
  numeric. Value of log-likelihood at ML estimates.
 }
 \item{efficiencies}{
  data frame. Contains point estimates of unit-specific technical or cost efficiencies:  exp(-E(u|e)) of Jondrow et al. (1982), E(exp(-u)|e) of Battese and Coelli (1988), and exp(-M(u|e)), where M(u|e) is the mode of conditional distribution of inefficiency term. In addition, estimated lower and upper bounds of (1-\eqn{\alpha})100\% two-sided prediction intervals are returned. 
 }
 \item{marg.effects}{
  data frame. Contains unit-specific marginal effects of exogenous variables on the expected value of inefficiency term. 
 }
 \item{sigmas_u}{
  matrix. Estimated unit-specific variances of inefficiency term. Returned if \code{uhet} is not \code{NULL}.
 }
  \item{sigmas_v}{
   matrix. Estimated unit-specific variances of random noise component. Returned if \code{vhet} is not \code{NULL}.
 }
  \item{mu}{
   matrix. Estimated unit-specific means of pre-truncated normal distribution of inefficiency term. Returned if \code{tmean} is not \code{NULL}.
 }
 \item{esample}{
  logical. Returns TRUE if the observation in user supplied data is in the estimation subsample and FALSE otherwise.
 }
}


\references{
Battese, G., Coelli, T. (1988), Prediction of firm-level technical effiiencies with a generalized frontier production function and panel data. \emph{Journal of Econometrics}, \bold{38}, 387--399.

Caudill, S., Ford, J., Gropper, D. (1995), Frontier estimation and firm-specific inefficiency measures in the presence of heteroscedasticity. \emph{Journal of Business and Economic Statistics}, \bold{13}, 105--111.

Horrace, W. and Schmidt, P. (1996), On ranking and selection from independent truncated normal distributions. \emph{Journal of Productivity Analysis}, \bold{7}, 257--282.

Jondrow, J., Lovell, C., Materov, I., Schmidt, P. (1982), On estimation of technical inefficiency in the stochastic frontier production function model. \emph{Journal of Econometrics}, \bold{19}, 233--238.

Kumbhakar, S. and Lovell, C. (2003), \emph{Stochastic Frontier Analysis.} Cambridge: Cambridge University Press.

Wang, H.-J. (2002), Heteroskedasticity and non-monotonic efficiency effects of a stochastic frontier model. \emph{Journal of Productivity Analysis}, \bold{18}, 241--253.
}

\author{
 Oleg Badunenko <obadunen@uni-koeln.de>, Yaryna Kolomiytseva <yaryna.kolomiytseva@uni-koeln.de>
}

\seealso{
\code{\link{teradial}}, \code{\link{tenonradial}}, \code{\link{teradialbc}}, \code{\link{nptestrts}}, \code{\link{nptestind}}
}

\examples{ 
require( npsf )
 
# Load Penn World Tables 5.6 dataset
 
data( pwt56 )
head( pwt56 )
 
# Create some missing values
 
pwt56 [4, "K"] <- NA 
 
# Stochastic production frontier model with 
# homoskedastic error components (half-normal)
 
# Use subset of observations - for year 1965
 
m1 <- sf(log(Y) ~ log(L) + log(K), data = pwt56, 
 subset = year == 1965, distribution = "h")

# Test CRS: 'car' package in required for that
\dontrun{linearHypothesis(m1, "log(L) + log(K) = 1")}
 
# Write efficiencies to the data frame using 'esample':

pwt56$BC[ m1$esample ] <- m1$efficiencies$BC
\dontrun{View(pwt56)}
  
# Computation using matrices
 
Y1 <- as.matrix(log(pwt56[pwt56$year == 1965, 
c("Y"), drop = FALSE]))
X1 <- as.matrix(log(pwt56[pwt56$year == 1965,
c("K", "L"), drop = FALSE]))

X1 [51, 2] <- NA # create missing
X1 [49, 1] <- NA # create missing
 
m2 <- sf(Y1 ~ X1, distribution = "h")
  
# Load U.S. commercial banks dataset
 
data(banks05)
head(banks05)
 
# Doubly heteroskedastic stochastic cost frontier 
# model (truncated normal)
 
# Print summaries of cost efficiencies' estimates
 
m3 <- sf(lnC ~ lnw1 + lnw2 + lny1 + lny2, uhet = ~ ER, 
 vhet = ~ LA, data = banks05, distribution = "t", 
 prod = FALSE, print.level = 3)

# Non-monotonic marginal effects of equity ratio on 
# the mean of distribution of inefficiency term
 
m4 <- sf(lnC ~ lnw1 + lnw2 + lny1 + lny2, uhet = ~ ER,
 tmean = ~ ER, data = banks05, distribution = "t", 
 prod = FALSE, marg.eff = TRUE)
 
summary(m4$marg.effects)

}

\keyword{ Stochastic Frontier Analysis }
\keyword{ Heteroskedasticity }
\keyword{ Truncated normal }
\keyword{ Parametric efficiency analysis }
\keyword{ summary }
