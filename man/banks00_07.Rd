\name{banks00_07}
\docType{data}
\alias{banks00_07}

\title{U.S. Commercial Banks Data}

\description{
\code{banks00_07} data frame contains selected variables for 500 (randomly sampled from around 5000) U.S. commercial banks from data of Koetter et al. (2012) for years 2000-2007.  This data are used for illustrution purposes and are not suitable for research purposes. 
}

\usage{data(banks00_07)}

\format{
 This data frame contains the following variables (columns):
  \describe{
  \item{\code{year}}{Year.}
  \item{\code{id}}{Entity ID.}
  \item{\code{TA}}{Gross total assets.}
  \item{\code{LLP}}{Loan loss provisions.}
  \item{\code{Y1}}{Total securities (in thousands of US dollars).}
  \item{\code{Y2}}{Total loans and leases (in thousands of US dollars).}
  \item{\code{W1}}{Cost of fixed assets divided by the cost of borrowed funds.}
  \item{\code{W2}}{Cost of labor (in thousands of US dollars) divided by the cost of borrowed funds.}
  \item{\code{ER}}{Gross total equity to gross total assets ratio.}
  \item{\code{TC}}{Total operating cost.}
  \item{\code{LA}}{Total loans and leases to gross total assets ratio.}
  \item{\code{Ti}}{Times bank is observed.}
  \item{\code{TA_ave}}{Mean value of TA.}
  \item{\code{TA_initial}}{Value of TA in the first observed year.}
  \item{\code{LLP_ave}}{Mean value of LLP.}
  \item{\code{LLP_initial}}{Value of LLP in the first observed year.}
  \item{\code{ER_ave}}{Mean value of ER.}
  \item{\code{ER_initial}}{Value of ER in the first observed year.}
  \item{\code{LA_ave}}{Mean value of LA.}
  \item{\code{LA_initial}}{Value of LA in the first observed year.}
}
}

\details{
 The data were sampled and generated as shown in section "Examples".
 
 }
 
\examples{

\dontrun{

# Download data from the link in "Source"

banks00_07 <- read.delim("2b_QLH.txt")

# rename 'entity' to 'id'
colnames(banks00_07) [colnames(banks00_07) == "entity"] <- "id"

table(banks00_07$year)
# keep if 2000 -- 2007

banks00_07 <- 
 banks00_07[(banks00_07$year >= 2000 & banks00_07$year <= 2007),]
dim(banks00_07)

q1q3 <- quantile(banks00_07$TA, probs = c(.25,.75))

banks00_07 <- 
 banks00_07[(banks00_07$TA >= q1q3[1] & banks00_07$TA <= q1q3[2]),]
dim(banks00_07)


# generate required variables
banks00_07$TC <-banks00_07$TOC
banks00_07$ER <- banks00_07$Z / banks00_07$TA
banks00_07$LA <- banks00_07$Y2 / banks00_07$TA

banks00_07 <- 
 banks00_07[, colnames(banks00_07) %in% 
 c("id", "year", "Ti", "TC", "Y1", "Y2", "W1","W2", "ER", "LA", "TA", "LLP")]
dim(banks00_07)

t0 <- as.vector( by(data = banks00_07$id, 
                    INDICES = banks00_07$id, 
                    FUN = function(qq) length(qq)) )
banks00_07$Ti <- rep(t0, times = t0)
banks00_07 <- banks00_07[banks00_07$Ti > 4,]

# complete observations
banks00_07 <- banks00_07[complete.cases(banks00_07),]
dim(banks00_07)

id_names <- unique(banks00_07$id)
N_total <- length(id_names)
set.seed(816376586)
ids_n2choose <- sample(1:N_total, 500)
ids2choose <- id_names[ids_n2choose]
banks00_07 <- banks00_07[banks00_07$id %in% ids2choose,]
dim(banks00_07)

t0 <- as.vector( by(data = banks00_07$id, 
                    INDICES = banks00_07$id, 
                    FUN = function(qq) length(qq)) )
length(rep(t0, times = t0))

banks00_07$Ti <- rep(t0, times = t0)

banks00_07[1:50,c("id","year","Ti")]

# keep if Ti > 4

banks00_07 <- banks00_07[banks00_07$Ti > 4,]
dim(banks00_07)

# sort
banks00_07 <- banks00_07[order(banks00_07$id, banks00_07$year),]

# TC = TOC
#
# ER = Z / TA
# Gross total equity to gross total assets ratio.
#
# LA = Y2 / TA
# Total loans and leases to gross total assets ratio.

banks00_07$TA_ave <- 
 rep(as.vector( by(data = banks00_07$TA,
                   INDICES = banks00_07$id, 
                   FUN = function(qq) mean(qq))), times = t0)

banks00_07$TA_initial <- 
 rep(as.vector( by(data = banks00_07$TA, 
                   INDICES = banks00_07$id,
                   FUN = function(qq) qq[1])), times = t0)

banks00_07$LLP_ave <- 
 rep(as.vector( by(data = banks00_07$LLP,
                   INDICES = banks00_07$id,
                   FUN = function(qq) mean(qq))), times = t0)

banks00_07$LLP_initial <- 
 rep(as.vector( by(data = banks00_07$LLP, 
                   INDICES = banks00_07$id, 
                   FUN = function(qq) qq[1])), times = t0)

banks00_07$ER_ave <- 
 rep(as.vector( by(data = banks00_07$ER, 
                   INDICES = banks00_07$id, 
                   FUN = function(qq) mean(qq))), times = t0)

banks00_07$ER_initial <- 
 rep(as.vector( by(data = banks00_07$ER, 
                   INDICES = banks00_07$id, 
                   FUN = function(qq) qq[1])), times = t0)

banks00_07$LA_ave <- 
 rep(as.vector( by(data = banks00_07$LA, 
                   INDICES = banks00_07$id, 
                   FUN = function(qq) mean(qq))), times = t0)

banks00_07$LA_initial <- 
 rep(as.vector( by(data = banks00_07$LA, 
                   INDICES = banks00_07$id, 
                   FUN = function(qq) qq[1])), times = t0)

cols2export <- c("id","year","Ti","TA","TA_ave",
                 "TA_initial","LLP","LLP_ave",
                 "LLP_initial","ER_ave","ER_initial","LA_ave","LA_initial")

write.table(x = banks00_07, file = "banks00_07.txt", row.names = FALSE)

}

}
 
 
 \source{
 \url{http://qed.econ.queensu.ca/jae/2014-v29.2/restrepo-tobon-kumbhakar/}.
 }
 
 \references{
Koetter, M., Kolari, J., and Spierdijk, L. (2012), Enjoying the quiet life under deregulation? Evidence from adjusted Lerner indices for U.S. banks. \emph{Review of Economics and Statistics}, \bold{94}, 2, 462--480.

Restrepo-Tobon, D. and Kumbhakar, S. (2014), Enjoying the quiet life under deregulation? Not Quite. \emph{Journal of Applied Econometrics}, \bold{29}, 2, 333--343.
 
 }
 

 
\keyword{datasets}
