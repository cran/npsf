\name{tenonradialbc}
\alias{tenonradialbc}
\title{
 Statistical Inference Regarding the Russell Measure of Technical Efficiency
}
\description{
 Routine \code{tenonradialbc} performs bias correction of the nonradial Russell input- or output-based measure of technical efficiency, computes bias and constructs confidence intervals via bootstrapping techniques.
}
\usage{
 tenonradialbc(formula, data, subset,
            ref = NULL, data.ref = NULL, subset.ref = NULL,
            rts = c("C", "NI", "V"), base = c("output", "input"),
            homogeneous = TRUE, smoothed = TRUE, kappa = NULL,
            reps = 999, level = 95,
            print.level = 1, show.progress = TRUE, seed = NULL)
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
 logical. Relevant if \code{smoothed=TRUE}. If TRUE, the reference set is bootstrapped with homogeneous smoothing; if FALSE, the reference set is bootstrapped with heterogeneous subsampling.
}
\item{kappa}{
 relevant if \code{smoothed=TRUE}. 'kappa' sets the size of the subsample as K^kappa, where K is the number of data points in the original reference set. The default value is 0.7. 'kappa' may be between 0.5 and 1.
}
\item{reps}{
 specifies the number of bootstrap replications to be performed.  The default is 999.  The minimum is 100.  Adequate estimates of confidence intervals using bias-corrected methods typically require 1,000 or more replications.
}
\item{level}{
 sets confidence level for confidence intervals; default is \code{level=95}.
}
\item{show.progress}{
 logical. Relevant if \code{print.level>=1}. If TRUE, progress of the bootstrap is displayed; if FALSE, display of the bootstrap progress is suppressed.
}
\item{print.level}{
 numeric. 0 - nothing is printed; 1 - print summary of the model and data. 2 - print summary of technical efficiency measures. 3 - print estimation results observation by observation. Default is 1.
}
\item{seed}{
 numeric. The seed (for replication purposes).
}
}
\details{
 Routine \code{tenonradialbc} performs bias correction of the nonradial Russell input- or output-based measure of technical efficiency, computes bias and constructs confidence intervals via bootstrapping techniques.
 
 Models for \code{tenonradialbc} are specified symbolically. A typical model has the form \code{outputs ~ inputs}, where \code{outputs} (\code{inputs}) is a series of (numeric) terms which specifies outputs (inputs). The same goes for reference set.  Refer to the examples.
 
 Results can be summarized using \code{\link{summary.npsf}}.
 
}

\value{
 \code{tenonradialbc} returns a list of class \code{npsf} containing the following elements:
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
 Badunenko, Oleg and Mozharovskyi, Pavlo (2020), Statistical inference for the Russell measure of technical efficiency, \emph{Journal of the Operational Research Society}, \bold{71}3, 517--527, https://doi.org/10.1080/01605682.2019.1599778
}

\author{
 Oleg Badunenko <oleg.badunenko@brunel.ac.uk>, Pavlo Mozharovskyi <pavlo.mozharovskyi@telecom-paris.fr>
  
}

\seealso{
 \code{\link{teradial}}, \code{\link{tenonradial}}, \code{\link{teradialbc}}, \code{\link{nptestrts}}, \code{\link{nptestind}}, \code{\link{sf}}
}

\examples{
 
 \dontrun{
  
  data( ccr81 )
  head( ccr81 )
  
  # Subsampling

  t9 <- tenonradialbc(y1 + y2 + y3 ~ x1 + x2 + x3 + x4 + x5, data = ccr81,
                          ref = NULL, data.ref = NULL, subset.ref = NULL,
                          rts = "v", base = "i",
                          homogeneous = FALSE, smoothed = TRUE, kappa = .6,
                          reps = 999, level = 95,
                          print.level = 1, show.progress = TRUE, seed = NULL)
  # display the results

  cbind(te = t9$te, telow = t9$telow, tebc = t9$tebc, teupp = t9$teupp, 
        biasboot = t9$biasboot, varboot = t9$varboot, biassqvar = t9$biassqvar)
  
  # Smoothing

  t10 <- tenonradialbc(y1 + y2 + y3 ~ x1 + x2 + x3 + x4 + x5, data = ccr81,
                          ref = NULL, data.ref = NULL, subset.ref = NULL,
                          rts = "v", base = "i",
                          homogeneous = TRUE, smoothed = TRUE, kappa = .6,
                          reps = 999, level = 95,
                          print.level = 1, show.progress = TRUE, seed = NULL)
  # display the results

  cbind(te = t10$te, telow = t10$telow, tebc = t10$tebc, teupp = t10$teupp, 
        biasboot = t10$biasboot, varboot = t10$varboot, biassqvar = t10$biassqvar)
  
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