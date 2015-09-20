\name{banks05}
\docType{data}
\alias{banks05}

\title{U.S. Commercial Banks Data}

\description{
\code{banks05} data frame contains selected variables from the U.S. commercial banks data of Koetter et al. (2012) for year 2005 and 500 banks randomly sampled from around 5000. Dependent variable was randomly generated, as described under 'Details', to satisfy the assumptions of doubly heteroskedastic stochastic cost frontier model. This data is, therefore, not suitable for research purposes. 
}

\usage{data(banks05)}

\format{
 This data frame contains the following variables (columns):
  \describe{
   \item{\code{lnC}}{Randomly generated total operating costs.}
   \item{\code{lnw1}}{Logarithm of the cost of fixed assets divided by the cost of borrowed funds.}
   \item{\code{lnw2}}{Logarithm of the cost of labor (in thousands of US dollars) divided by the cost of borrowed funds.}
   \item{\code{lny1}}{Logarithm of total securities (in thousands of US dollars).}
   \item{\code{lny2}}{Logarithm of total loans and leases (in thousands of US dollars).}
   \item{\code{ER}}{Gross total equity to gross total assets ratio.}
   \item{\code{LA}}{Total loans and leases to gross total assets ratio.}
  }
}


\details{
The variable representing total operating costs was generated as follows: \deqn{lnC = \beta0 + \beta1*lnw1 + \beta2*lnw2 + \gamma1*lny1 + \gamma2*lny2 + \nu + u,} where \eqn{\nu ~ N(0, exp(\alpha0 + \alpha1*LA))} and \eqn{u ~ N+(\delta0 + \delta1*ER, exp(\omega0 + \omega1*ER)).} More detailed description of input prices, outputs, and exogenous variables is provided in Koetter et al. (2012). See also related study of Restrepo-Tobon and Kumbhakar (2014).
}
 
 \source{
 \url{http://qed.econ.queensu.ca/jae/2014-v29.2/restrepo-tobon-kumbhakar/}.
 }
 
 \references{
Koetter, M., Kolari, J., and Spierdijk, L. (2012), Enjoying the quiet life under deregulation? Evidence from adjusted Lerner indices for U.S. banks. \emph{Review of Economics and Statistics}, \bold{94}, 2, 462--480.

Restrepo-Tobon, D. and Kumbhakar, S. (2014), Enjoying the quiet life under deregulation? Not Quite. \emph{Journal of Applied Econometrics}, \bold{29}, 2, 333--343.
 
 }
 
 
\keyword{datasets}
