\name{sf}
\alias{sf}
\title{
 Stochastic Frontier Models Using Cross-Sectional and Panel Data
}
\description{
 \code{sf} performs maximum likelihood estimation of the parameters and technical or cost efficiencies in cross-sectional stochastic (production or cost) frontier models with half-normal or truncated normal distributional assumption imposed on inefficiency error component. 
}
\usage{
sf(formula, data, it = NULL, subset,
 prod = TRUE, model = "K1990", distribution = c("h"),
 eff.time.invariant = FALSE, 
 mean.u.0i.zero     = FALSE,
 mean.u.0i          = NULL,
 ln.var.u.0i        = NULL,
 ln.var.v.0i        = NULL,
 ln.var.v.it        = NULL,  
 simtype = c("halton", "rnorm"), halton.base = 2, R = 500,
 initialize.halton = TRUE, my.leap = NULL,
 simtype_GHK = c("halton", "runif"), R_GHK = 1000,
 cost.eff.less.one  = FALSE, level = 95, marg.eff = FALSE,
 start.val = NULL, maxit = 199, report.ll.optim = 10, 
 reltol = 1e-8, lmtol = sqrt(.Machine$double.eps),
 digits = 4, print.level = 4)
}
\arguments{
 \item{formula}{
  an object of class ``formula'' (or one that can be coerced to that class): a symbolic description of the model. The details of model specification are given under `Details'.
 }
 \item{data}{
  an optional data frame containing variables in the model. If not found in data, the variables are taken from environment (\code{formula}), typically the environment from which \code{sf} is called.
  }
 \item{it}{
 vector with two character entries.  E.g., c("ID", "TIME"), where "ID" defines individuals that are observed in time periods defined by "TIME". The default is \code{NULL}.  At default, cross-sectional model will be estimated. 
 }
 \item{subset}{
  an optional vector specifying a subset of observations for which technical or cost efficiencies are to be computed.
  }
 \item{prod}{
 logical. If \code{TRUE}, the estimates of parameters of stochastic production frontier model and of technical efficiencies are returned; if \code{FALSE}, the estimates of parameters of stochastic cost frontier model and of cost efficiencies are returned.
 }
 \item{model}{
 character. Five panel data models are estimated for now. "K1990" and "K1990modified" (see Kumbhakar, 1990), "BC1992" (see Battese and Coelli, 1992), "4comp" (see Badunenko and Kumbhakar (2016) and Filippini and Greene, 2016). They specify common evolution of inefficiency.  Deffault is "K1990". The time functions are "( 1 + exp(b*t + c*t^2) )^-1", "1 + d*(t-T_i) + e*(t-T_i)^2", and "exp( -f*(t-T_i) )", respectively.
 }
 \item{distribution}{
 either \code{"h"} (half-normal) or \code{"t"} (truncated normal), specifying the distribution of inefficiency term.
 }
 \item{eff.time.invariant}{
 logical. If \code{TRUE}, the 1st generation of panel data models is estimated, otherwise, the 2nd generation or 4 components panel data model is estimated.
 }
 \item{mean.u.0i.zero}{
 logical. If \code{TRUE}, normal-half normal model is estimated, otherwise, normal-truncated model is estimated.
 }
 \item{mean.u.0i}{
 one-sided formula; e.g. \code{mean.u.0i ~ z1 + z2}. Specifies whether the mean of pre-truncated normal distribution of inefficiency term is a linear function of exogenous variables. In cross-sectional models, used only when \code{distribution = "t"}. If \code{NULL}, mean is assumed to be constant for all ids. 
 }
 \item{ln.var.u.0i}{
 one-sided formula; e.g. \code{ln.var.u.0i ~ z1 + z2}. Specifies exogenous variables entering the expression for the log of variance of inefficiency error component. If \code{NULL}, inefficiency term is assumed to be homoskedastic, i.e. \eqn{\sigma_u^2 = exp(\gamma[0])}. Time invariant variables are expected.
 }
 \item{ln.var.v.0i}{
 one-sided formula; e.g. \code{ln.var.v.0i ~ z1 + z2}. Specifies exogenous variables entering the expression for variance of random noise error component. If \code{NULL}, random noise component is assumed to be homoskedastic, i.e. \eqn{\sigma_v^2 = exp(\gamma[0])}. Time invariant variables are expected.
 }
 \item{ln.var.v.it}{
 one-sided formula; e.g. \code{ln.var.v.it ~ z1 + z2}. Specifies exogenous variables entering the expression for variance of random noise error component. If \code{NULL}, random noise component is assumed to be homoskedastic, i.e. \eqn{\sigma_v^2 = exp(\gamma[0])}. Time invariant variables are expected.
 }
 \item{simtype}{
 character. Type of random deviates for the 4 components model. 'halton' draws are default.  One can specify 'rnorm.'
 }
 \item{halton.base}{
 numeric.  The prime number which is the base for the Halton draws. Should be less then 102.  Default is 2.
 }
 \item{R}{
 numeric.  Number of draws.  Default is 500.  Can be time consuming.
 }
 \item{initialize.halton}{
 logical.  If TRUE, the Halton draws use different base for each id.  Then \code{halton.base} is ignored.  Function \code{\link{halton}} from \code{\link{randtoolbox}} package is used.  Otherwise, N*R Halton draws come from the same base. Function \code{\link{sHalton}} from \code{sfsmisc} package is used.
 }
 \item{my.leap}{
 if \code{initialize.halton} is FALSE, \code{my.leap} draws are skipped. See \code{\link{sHalton}}.  By default none are skipped.
 }
 \item{simtype_GHK}{
 character. Type of random deviates for use in GHK for efficiency estimating by approximation. 'halton' draws are default.  One can specify 'runif.'
 }
 \item{R_GHK}{
 numeric.  Number of draws for GHK.  Default is 1000.  Can be time consuming.
 }
 \item{cost.eff.less.one}{
 logical. If \code{TRUE}, and \code{prod = FALSE}, reported cost efficiencies are larger than one. Interpretation: by what factor is total cost larger than the potential total cost.
 }
 \item{level}{
 numeric. Defines level\% two-sided prediction interval for technical or cost efficiencies (see Horrace and Schmidt 1996). Default is 95.
 }
 \item{marg.eff}{
 logical. If \code{TRUE}, unit-specific marginal effects of exogenous variables on the mean of distribution of inefficiency term are returned.
 }
 \item{start.val}{
 numeric. Starting values to be supplied to the optimization routine. If \code{NULL}, OLS and method of moments estimates are used (see Kumbhakar and Lovell 2000).
 }
 \item{maxit}{
 numeric. Maximum number of iterations. Default is 199.
 }
 \item{report.ll.optim}{
 numeric. Not used for now.
 } 
 \item{reltol}{
 numeric. One of convergence criteria. Not used for now.
 }
 \item{lmtol}{
 numeric. Tolerance for the scaled gradient in ML optimization. Default is sqrt(.Machine$double.eps).
 }
 \item{digits}{
 numeric. Number of digits to be displayed in estimation results and for efficiency estimates. Default is 4.
 }
 \item{print.level}{
  numeric. 1 - print estimation results. 2 - print optimization details. 3 - print summary of point estimates of technical or cost efficiencies. 7 - print unit-specific point and interval estimates of technical or cost efficiencies. Default is 4.
 }
}
  
  
\details{
 Models for \code{sf} are specified symbolically. A typical model has the form \code{y ~ x1 + ...}, where \code{y} represents the logarithm of outputs or total costs and \code{\{x1,...\}} is a series of inputs or outputs and input prices (in logs).  

Options \code{ln.var.u.0i} and \code{ln.var.v.0i} can be used if multiplicative heteroskedasticity of either inefficiency or random noise component (or both) is assumed; i.e. if their variances can be expressed as exponential functions of (e.g. size-related) exogenous variables (including intercept) (see Caudill et al. 1995).

If \code{marg.eff = TRUE} and \code{distribution = "h"}, the marginal effect of \emph{k}th exogenous variable on the expected value of inefficiency term of unit \emph{i} is computed as: \eqn{\gamma[k]\sigma[i]/\sqrt2\pi}, where \eqn{\sigma_u[i] = \sqrt exp(z[i]'\gamma)}. If \code{distribution = "t"}, marginal effects are returned if either \code{mean.u.0i} or \code{ln.var.u.0i} are not \code{NULL}. If the same exogenous variables are specified under both options, (non-monotonic) marginal effects are computed as explained in Wang (2002).

}

\value{
 \code{sf} returns a list of class \code{npsf} containing the following elements:
  \item{coef}{
   numeric. Named vector of ML parameter estimates.
  }
 \item{vcov}{
  matrix. Estimated covariance matrix of ML estimator. 
 }
 \item{ll}{
  numeric. Value of log-likelihood at ML estimates.
 }
 \item{efficiencies}{
  data frame. Contains point estimates of unit-specific technical or cost efficiencies:  exp(-E(u|e)) of Jondrow et al. (1982), E(exp(-u)|e) of Battese and Coelli (1988), and exp(-M(u|e)), where M(u|e) is the mode of conditional distribution of inefficiency term. In addition, estimated lower and upper bounds of (1-\eqn{\alpha})100\% two-sided prediction intervals are returned. 
 }
 \item{marg.effects}{
  data frame. Contains unit-specific marginal effects of exogenous variables on the expected value of inefficiency term. 
 }
 \item{sigmas_u}{
  matrix. Estimated unit-specific variances of inefficiency term. Returned if \code{ln.var.u.0i} is not \code{NULL}.
 }
  \item{sigmas_v}{
   matrix. Estimated unit-specific variances of random noise component. Returned if \code{ln.var.v.0i} is not \code{NULL}.
 }
  \item{mu}{
   matrix. Estimated unit-specific means of pre-truncated normal distribution of inefficiency term. Returned if \code{mean.u.0i} is not \code{NULL}.
 }
 \item{esample}{
  logical. Returns TRUE if the observation in user supplied data is in the estimation subsample and FALSE otherwise.
 }
}


\references{

Badunenko, O. and Kumbhakar, S.C. (2016), When, Where and How to Estimate Persistent and Transient Efficiency in Stochastic Frontier Panel Data Models, \emph{European Journal of Operational Research}, \bold{255}(1), 272--287.

Battese, G., Coelli, T. (1988), Prediction of firm-level technical effiiencies with a generalized frontier production function and panel data. \emph{Journal of Econometrics}, \bold{38}, 387--399.

Battese, G., Coelli, T. (1992), Frontier production functions, technical efficiency and panel data: With application to paddy farmers in India. \emph{Journal of Productivity Analysis}, \bold{3}, 153--169.

Caudill, S., Ford, J., Gropper, D. (1995), Frontier estimation and firm-specific inefficiency measures in the presence of heteroscedasticity. \emph{Journal of Business and Economic Statistics}, \bold{13}, 105--111.

Filippini, M. and Greene, W.H. (2016), Persistent and transient productive inefficiency: A maximum simulated likelihood approach. \emph{Journal of Productivity Analysis}, \bold{45} (2), 187--196.

Horrace, W. and Schmidt, P. (1996), On ranking and selection from independent truncated normal distributions. \emph{Journal of Productivity Analysis}, \bold{7}, 257--282.

Jondrow, J., Lovell, C., Materov, I., Schmidt, P. (1982), On estimation of technical inefficiency in the stochastic frontier production function model. \emph{Journal of Econometrics}, \bold{19}, 233--238.

Kumbhakar, S. (1990),  Production Frontiers, Panel Data, and Time-varying Technical Inefficiency. \emph{Journal of Econometrics}, \bold{46}, 201--211.

Kumbhakar, S. and Lovell, C. (2003), \emph{Stochastic Frontier Analysis.} Cambridge: Cambridge University Press.

Wang, H.-J. (2002), Heteroskedasticity and non-monotonic efficiency effects of a stochastic frontier model. \emph{Journal of Productivity Analysis}, \bold{18}, 241--253.
}

\author{
 Oleg Badunenko <oleg.badunenko@port.ac.uk>, Yaryna Kolomiytseva <yaryna.kolomiytseva@uni-koeln.de>
}

\seealso{
\code{\link{teradial}}, \code{\link{tenonradial}}, \code{\link{teradialbc}}, \code{\link{nptestrts}}, \code{\link{nptestind}}
}

\examples{ 
require( npsf )

# Cross-sectional examples begin ------------------------------------------

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
## Not run: linearHypothesis(m1, "log(L) + log(K) = 1")

# Write efficiencies to the data frame using 'esample':

pwt56$BC[ m1$esample ] <- m1$efficiencies$BC
## Not run: View(pwt56)

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

m3 <- sf(lnC ~ lnw1 + lnw2 + lny1 + lny2, ln.var.u.0i = ~ ER, 
         ln.var.v.0i = ~ LA, data = banks05, distribution = "t", 
         prod = FALSE, print.level = 3)

# Non-monotonic marginal effects of equity ratio on 
# the mean of distribution of inefficiency term

m4 <- sf(lnC ~ lnw1 + lnw2 + lny1 + lny2, ln.var.u.0i = ~ ER,
         mean.u.0i = ~ ER, data = banks05, distribution = "t", 
         prod = FALSE, marg.eff = TRUE)

summary(m4$marg.effects)


# Cross-sectional examples end --------------------------------------------

\dontrun{

# Panel data examples begin -----------------------------------------------

# The only way to differentiate between cross-sectional and panel-data
# models is by specifying "it".
# If "it" is not specified, cross-sectional model will be estimated.
# Example is below.

# Read data ---------------------------------------------------------------

# Load U.S. commercial banks dataset

data(banks00_07)
head(banks00_07, 5)

banks00_07$trend <- banks00_07$year - min(banks00_07$year) + 1

# Model specification -----------------------------------------------------

my.prod     <- FALSE
my.it       <- c("id","year")

# my.model = "BC1992"
# my.model = "K1990modified"
# my.model = "K1990"

# translog ----------------------------------------------------------------
formu <- log(TC) ~ (log(Y1) + log(Y2) + log(W1) + log(W2) + trend)^2 +
 I(0.5*log(Y1)^2) + I(0.5*log(Y2)^2) + I(0.5*log(W1)^2) + I(0.5*log(W2)^2) +
 trend + I(0.5*trend^2)

# Cobb-Douglas ------------------------------------------------------------
# formu <- log(TC) ~ log(Y1) + log(Y2) + log(W1) + log(W2) + trend + I(trend^2)

ols <- lm(formu, data = banks00_07)
# just to mark the results of the OLS model
summary(ols)

# Components specifications -----------------------------------------------

ln.var.v.it <- ~ log(TA)
ln.var.u.0i <- ~ ER_ave
mean.u.0i_1 <- ~ LLP_ave + LA_ave
mean.u.0i_2 <- ~ LLP_ave + LA_ave - 1

# Suppose "it" is not specified
# Then it is a cross-sectional model

t0_1st_0 <- sf(formu, data = banks00_07, subset = year > 2000,
               prod = my.prod,
               ln.var.v.it = ln.var.v.it,
               ln.var.u.0i = ln.var.u.0i,
               eff.time.invariant = TRUE,
               mean.u.0i.zero = TRUE)

# 1st generation models ---------------------------------------------------

# normal-half normal ------------------------------------------------------

# the same as above but "it" is specified

t0_1st_0 <- sf(formu, data = banks00_07, it = my.it, subset = year > 2000,
               prod = my.prod,
               ln.var.v.it = ln.var.v.it,
               ln.var.u.0i = ln.var.u.0i, 
               eff.time.invariant = TRUE, 
               mean.u.0i.zero = TRUE)

# Note efficiencies are time-invariant

# confidence intervals for efficiencies -----------------------------------

head(t0_1st_0$efficiencies, 20)


# normal-truncated normal -------------------------------------------------

# truncation point is constant (for all ids) ------------------------------

t0_1st_1 <- sf(formu, data = banks00_07, it = my.it, subset = year > 2000,
               prod = my.prod,
               eff.time.invariant = TRUE,
               mean.u.0i.zero = FALSE,
               ln.var.v.it = ln.var.v.it,
               ln.var.u.0i = ln.var.u.0i,
               mean.u.0i = NULL,
               cost.eff.less.one = TRUE)



# truncation point is determined by variables -----------------------------

t0_1st_2 <- sf(formu, data = banks00_07, it = my.it, subset = year > 2000,
               prod = my.prod,
               eff.time.invariant = TRUE,
               mean.u.0i.zero = FALSE,
               mean.u.0i = mean.u.0i_1,
               ln.var.v.it = ln.var.v.it,
               ln.var.u.0i = ln.var.u.0i,
               cost.eff.less.one = TRUE)



# the same, but without intercept in mean.u.0i

t0_1st_3 <- sf(formu, data = banks00_07, it = my.it, subset = year > 2000,
               prod = my.prod,
               eff.time.invariant = TRUE,
               mean.u.0i.zero = FALSE,
               mean.u.0i = mean.u.0i_2,
               ln.var.v.it = ln.var.v.it,
               ln.var.u.0i = ln.var.u.0i,
               cost.eff.less.one = TRUE)

# 2nd generation models ---------------------------------------------------

# normal-half normal ------------------------------------------------------

# Kumbhakar (1990) model --------------------------------------------------

t_nhn_K1990 <- sf(formu, data = banks00_07, it = my.it, subset = year > 2000,
                  prod = my.prod,
                  eff.time.invariant = FALSE,
                  mean.u.0i.zero = TRUE, 
                  ln.var.v.it = ln.var.v.it,
                  ln.var.u.0i = ln.var.u.0i, 
                  cost.eff.less.one = FALSE)


# Kumbhakar (1990) modified model -----------------------------------------

t_nhn_K1990modified <- sf(formu, data = banks00_07, it = my.it, subset = year > 2000,
                          prod = my.prod, model = "K1990modified",
                          eff.time.invariant = FALSE,
                          mean.u.0i.zero = TRUE, 
                          ln.var.v.it = ln.var.v.it,
                          ln.var.u.0i = ln.var.u.0i, 
                          cost.eff.less.one = FALSE)


# Battese and Coelli (1992) model -----------------------------------------

t_nhn_BC1992 <- sf(formu, data = banks00_07, it = my.it, subset = year > 2000,
                   prod = my.prod, model = "BC1992",
                   eff.time.invariant = FALSE,
                   mean.u.0i.zero = TRUE, 
                   ln.var.v.it = ln.var.v.it,
                   ln.var.u.0i = ln.var.u.0i, 
                   cost.eff.less.one = FALSE)

# normal-truncated normal -------------------------------------------------

# Kumbhakar (1990) model --------------------------------------------------

# truncation point is constant (for all ids) ------------------------------

t_ntn_K1990_0 <- sf(formu, data = banks00_07, it = my.it, subset = year > 2000,
                    prod = my.prod,
                    eff.time.invariant = FALSE, 
                    mean.u.0i.zero = FALSE,
                    ln.var.v.it = ln.var.v.it, 
                    ln.var.u.0i = ln.var.u.0i,
                    cost.eff.less.one = FALSE)


# truncation point is determined by variables -----------------------------

t_ntn_K1990_1 <- sf(formu, data = banks00_07, it = my.it, subset = year > 2000,
                    prod = my.prod,
                    eff.time.invariant = FALSE, 
                    mean.u.0i.zero = FALSE, 
                    mean.u.0i = mean.u.0i_1,
                    ln.var.v.it = ln.var.v.it, 
                    ln.var.u.0i = ln.var.u.0i,
                    cost.eff.less.one = FALSE)

# Efficiencies are tiny, since empirically truncation points are quite big.
# Try withouth an intercept in conditional mean f-n

t_ntn_K1990_2 <- sf(formu, data = banks00_07, it = my.it, subset = year > 2000,
                    prod = my.prod,
                    eff.time.invariant = FALSE, 
                    mean.u.0i.zero = FALSE, 
                    mean.u.0i = mean.u.0i_2,
                    ln.var.v.it = ln.var.v.it, 
                    ln.var.u.0i = ln.var.u.0i,
                    cost.eff.less.one = FALSE)

# Kumbhakar (1990) modified model -----------------------------------------

t_ntn_K1990modified <- sf(formu, data = banks00_07, it = my.it, subset = year > 2000,
                          prod = my.prod, model = "K1990modified",
                          eff.time.invariant = FALSE, 
                          mean.u.0i.zero = FALSE, 
                          mean.u.0i = mean.u.0i_1,
                          ln.var.v.it = ln.var.v.it, 
                          ln.var.u.0i = ln.var.u.0i,
                          cost.eff.less.one = FALSE)

# Battese and Coelli (1992) model -----------------------------------------


t_ntn_BC1992 <- sf(formu, data = banks00_07, it = my.it, subset = year > 2000,
                   prod = my.prod, model = "BC1992",
                   eff.time.invariant = FALSE, 
                   mean.u.0i.zero = FALSE, 
                   mean.u.0i = mean.u.0i_1,
                   ln.var.v.it = ln.var.v.it, 
                   ln.var.u.0i = ln.var.u.0i,
                   cost.eff.less.one = FALSE)

# The next one (without "subset = year > 2000" option) converges OK

t_ntn_BC1992 <- sf(formu, data = banks00_07, it = my.it,
                   prod = my.prod, model = "BC1992",
                   eff.time.invariant = FALSE, 
                   mean.u.0i.zero = FALSE, 
                   mean.u.0i = mean.u.0i_1,
                   ln.var.v.it = ln.var.v.it, 
                   ln.var.u.0i = ln.var.u.0i,
                   cost.eff.less.one = FALSE)

# 4 component model ------------------------------------------------------

# Note, R should better be more than 200, this is just for illustration.
# This is the model that takes long to be estimated.  
# For the following example, 'mlmaximize' required 357 iterations and
# took 8 minutes.  
# The time will increase with more data and more parameters.

formu <- log(TC) ~ log(Y1) + log(Y2) + log(W1) + log(W2) + trend

t_4comp <- sf(formu, data = banks00_07, it = my.it, 
              subset = year >= 2001 & year < 2006,
              prod = my.prod, model = "4comp",
              R = 200, initialize.halton = TRUE, 
              lmtol = 1e-5, maxit = 500, print.level = 4)

# @e_i0, @e_it, and @e_over give efficiencies, 
# where @ is either 'c' or 't' for cost or production function.
# e.g., t_ntn_4comp$ce_i0 from last model, gives persistent cost efficiencies

# Panel data examples end -------------------------------------------------

}

}

\keyword{ Stochastic Frontier Analysis }
\keyword{ Panel Data Models }
\keyword{ Panel Data Analysis }
\keyword{ Heteroskedasticity }
\keyword{ Truncated normal }
\keyword{ Half normal }
\keyword{ Halton draws }
\keyword{ Parametric efficiency analysis }
\keyword{ summary }
