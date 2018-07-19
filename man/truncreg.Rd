\name{truncreg}
\alias{truncreg}
\title{
 Parametric truncated regression for cross-sectional data
}
\description{
 \code{truncreg} performs maximum likelihood estimation of the parameters in cross-sectional truncated regression. 
}
\usage{
truncreg(formula, data, subset, 
 ll = -Inf, ul = Inf, 
 lmtol = .Machine$double.eps, maxiter = 150, 
 marg.eff = FALSE, 
 print.level = 1)
}
\arguments{
 \item{formula}{
  an object of class ``formula'' (or one that can be coerced to that class): a symbolic description of the model. The details of model specification are given under `Details'.
 }
  \item{data}{
  an optional data frame containing variables in the model. If not found in data, the variables are taken from environment (\code{formula}), typically the environment from which \code{sf} is called.
  }
  \item{subset}{
  an optional vector specifying a subset of observations for which technical or cost efficiencies are to be computed.
  }
 \item{ll}{
 scalar or one-sided formula for lower limit for left-truncation; e.g. \code{ll = ~ 1} or \code{ll = ~ z1}.
 }
 \item{ul}{
 scalar or one-sided formula for upper limit for right-truncation; e.g. \code{ul = ~ 800} or \code{ul = ~ z1}.
 }
 \item{lmtol}{
 numeric. Tolerance for the scaled gradient in ML optimization. Default is .Machine$double.eps.
 }
  
 \item{maxiter}{
 numeric. maximum number of iteration for maximization of the log likelihood function.
 }
 \item{marg.eff}{
 logical. If \code{TRUE}, unit-specific marginal effects of exogenous variables on the mean of distribution of inefficiency term are returned.
 }
 \item{print.level}{
  numeric. 0 - nothing is printed. 1 - optimization steps and print estimation results. 2 - detailed optimization steps and print estimation results. Default is 1.
 }
}
  
\details{
 \code{truncreg} performs a regression from a sample drawn from a restricted part of the population.  Under the assumption that the error term of the whole population is normal, the error terms in the truncated regression model have a truncated normal distribution.
 
 Both lower limit for left-truncation and upper limit for right-truncation can be specified simultaneously.

 Models for \code{truncreg} are specified symbolically. A typical model has the form \code{y ~ x1 + ...}, where \code{y} represents the left hand side variable and \code{\{x1,...\}} right hand side variables.  

 If \code{marg.eff = TRUE}, the marginal effects are computed.

}

\value{
 \code{truncreg} returns a list of class \code{npsf} containing the following elements:
 \item{call}{
   call. 'truncreg.cs'.
  }
  \item{model}{
   character. Unevaluated call to function \code{truncreg}.
  }
  \item{coef}{
   numeric. Named vector of ML parameter estimates.
  }
  \item{table}{
  matrix. Table with results.
  }
 \item{vcov}{
  matrix. Estimated covariance matrix of ML estimator. 
 }
 \item{ll}{
  numeric. Value of log-likelihood at ML estimates.
 }
 \item{lmtol}{
 numeric. Convergence criterion: tolerance for the scaled gradient.
 }
 \item{LM}{
 numeric. Value of the scaled gradient.
 }
 \item{esttime}{
 numeric. Estimation time.
 }
 \item{marg.effects}{
  data frame. Contains unit-specific marginal effects of exogenous variables. 
 }
 \item{sigma}{
  numeric. estimate of sigma.
 }
 \item{LL}{
 numeric. The lower limit for left-truncation
 }
 \item{UL}{
 numeric. The upper limit for left-truncation
 }
 \item{n}{
 numeric. Number of observations (used in regression).
 }
 \item{n.full}{
 numeric. Number of observations (used and not used in regression).
 }
 \item{nontruncsample}{
   logical. Returns TRUE if the observation in user supplied data is in the estimation subsample and in non-truncated part of the sample, and FALSE otherwise.
 }
 \item{esample}{
  logical. Returns TRUE if the observation in user supplied data is in the estimation subsample and FALSE otherwise.
 }
}

\author{
 Oleg Badunenko <oleg.badunenko@port.ac.uk>
}

\seealso{
\code{\link{teradial}}, \code{\link{tenonradial}}, \code{\link{teradialbc}}, \code{\link{nptestrts}}, \code{\link{nptestind}}, \code{\link{sf}}
}

\examples{ 
require( npsf )
 
# Female labor force participation dataset
 
data(mroz)
head(mroz)

t1 <- npsf::truncreg(hours ~ kidslt6 + kidsge6 + age*educ, 
 data = mroz, ll = 0, lmtol = 1e-16, print.level = 2)
 
# matrices also can be used
myY <- mroz$hours
myX <- cbind(mroz$kidslt6, mroz$kidsge6, mroz$age, mroz$educ, mroz$age*mroz$educ)

t1.m <- truncreg(myY ~ myX, ll = 0)

# gives identical result to `t1':
# compare summary(t1) and summary(t1.m)

# using variable for limits is admissible
# we obtain the same result as before

mroz$myll <- 0
t11 <- npsf::truncreg(hours ~ kidslt6 + kidsge6 + age*educ, 
 data = mroz, ll = ~ myll, lmtol = 1e-16, print.level = 0)
summary(t11)

# if you believe that the sample is additionally truncted from above at say 3500

t12 <- update(t1, ul = 3500, print.level = 1)

# for obtaining marginal effects

t13 <- update(t12, marg.eff = TRUE)

summary(t13$marg.effects)

}

\keyword{ Truncated regression }
\keyword{ summary }
