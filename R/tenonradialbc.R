tenonradialbc <- function(formula, data, subset,
                       ref = NULL, data.ref = NULL, subset.ref = NULL,
                       rts = c("C", "NI", "V"), base = c("output", "input"),
                       homogeneous = TRUE, smoothed = TRUE, kappa = NULL,
                       reps = 999, level = 95,
                       print.level = 1, show.progress = TRUE, seed = NULL){
 if( !is.null(ref) & is.null(data.ref) ){
  warning("If you use variable names in 'ref', 'data.ref' is required", call. = FALSE)
 }

 if( is.null(ref) & !is.null(data.ref) ){
  stop("If you use 'data.ref', 'ref' is required", call. = FALSE)
 }
 
 if( !homogeneous ){
  subsampling <- TRUE
  warning("Subsampling will be used")
  if(!is.null(kappa)){
   if (kappa <= 0.5 | kappa >= 1) {
    stop("'kappa' must be between 0.5 and 1")
   }
  }
  else {
   if(!smoothed){
    stop("'kappa' must be provided for subsampling bootstrap")
   }
  }
 } else {
  subsampling <- FALSE
 }
 
 if (level < 10 | level > 99.99) {
  stop("'level' must be between 10 and 99.99 inclusive")
 }
 
 
 
 if (reps < 100) {
  stop("'reps' must be at least 100")
 }
 
 if(print.level == 0){
  dots <- FALSE
 }
 
 my.warnings <- NULL
 
 if (reps < 200) {
  warning(" Statistical inference may be unreliable \n          for small number of bootstrap replications; \n          consider setting 'reps' larger than 200", call. = FALSE, immediate. = TRUE)
  warning(" Statistical inference may be unreliable for small number of bootstrap replications; consider setting 'reps' larger than 200\n", call. = FALSE, immediate. = FALSE)
  my.warnings <- c(my.warnings, " Statistical inference may be unreliable for small number of bootstrap replications; consider setting 'reps' larger than 200.")
 }
 
 if (reps > 2000) {
  warning(" Unnecessary too many bootstrap replications; \n          consider setting 'reps' smaller than 2000", call. = FALSE, immediate. = TRUE)
  warning(" Unnecessary too many bootstrap replications; consider setting 'reps' smaller than 2000\n", call. = FALSE, immediate. = FALSE)
  my.warnings <- c(my.warnings, " Unnecessary too many bootstrap replications; consider setting 'reps' smaller than 2000.")
 }
 
 # # begin require for parallel computing
 # if (!is.numeric(core.count) || floor(core.count) != core.count || 
 #     core.count < 1){
 #  stop("'core.count' must be a positive integer", call. = FALSE)
 # }
 # 
 # if (core.count > 1){
 #  if (!(is.character(cl.type)) || !(cl.type %in% c("SOCK", "MPI"))){
 #   stop("invalid cluster type in 'cl.type'; 'cl.type' must be \"SOCK\" or \"MPI\"", call. = FALSE)
 #  }
 #  if (!("snowFT" %in% rownames(installed.packages()))){
 #   # mymessage <- unlist(strsplit("Package 'snowFT' required for parallel computing is not installed; type -install.packages(''snowFT'')- to install it", split = " "))
 #   stop("Package 'snowFT' required for parallel computing is not installed; \nuse -install.packages(",paste(dQuote("snowFT")), ")- to install it\n")
 #  }else{
 #   # require("snowFT")
 #   requireNamespace("snowFT", quietly = TRUE)
 #  }
 #  if (cl.type == "MPI"){
 #   if (!("Rmpi" %in% rownames(installed.packages()))){
 #    stop("Package 'Rmpi' required for parallel computing with option 'MPI' is not installed; \nuse -install.packages(",paste(dQuote("Rmpi")), ")- to install it\n")
 #   }else{
 #    # require("Rmpi")
 #    requireNamespace("Rmpi", quietly = TRUE)
 #   }
 #  }
 # }
 # # end require for parallel computing
 
 winw <- getOption("width")
 if (winw > 100+5){
  winw <- 100
 }
 else if (winw > 90+5) {
  winw <- 90
 }
 else if (winw > 80+5) {
  winw <- 80
 }
 else if (winw > 70+5) {
  winw <- 70
 }
 else if (winw > 60+5) {
  winw <- 60
 }
 else if (winw > 50+5) {
  winw <- 50
 }
 else {
  winw <- 0
 }
 
 # get the data in matrices
 
 YX <- .prepareYX(formula = formula, data = data, subset = subset, rts = rts,
                  base = base, ref = ref,	data.ref = data.ref, subset.ref = subset.ref,
                  print.level = print.level, type = "RM", winw = winw, 
                  sysnframe = sys.nframe())
 Y  <- t(YX$y)
 X  <- t(YX$x)
 M  <- nrow(Y) # number of Ys
 N  <- nrow(X) # number of Xs
 K  <- ncol(Y) # number of obs
 Yr <- t(YX$y.ref)
 Xr <- t(YX$x.ref)
 Kr <- ncol(Yr)
 rt <- YX$myrts
 ba <- my.base <- YX$mybase
 esample <- YX$esample
 
 # cat.print( dim(Y) )
 # cat.print( dim(X) )
 # cat.print( dim(Yr) )
 # cat.print( dim(Xr) )
 # cat.print( M )
 # cat.print( N )
 # cat.print( K )
 # cat.print( Kr )
 

 # original Russell measures
 
 # te0 <- .teRad(t(Y),t(X),M,N,K,t(Yr),t(Xr),Kr,rt,ba,0,print.level=print.level)
 # eff <- nonradial.c(YOBS = t(Y), XOBS = t(X), YREF = t(Yr), XREF = t(Xr), RTS = rt, base = ba, lmdConstr = TRUE, print.level = print.level)
 
 eff <- .teNonrad2(Y, X, M, N, K, Yr, Xr, Kr, rt, ba,
                  ifqh = FALSE, print.level = 0,
                  lmdConstr = TRUE, full.solution = TRUE)
 
 
 # te <- .teRad(t(t1$y), t(t1$x), ncol(t1$y), ncol(t1$x), nrow(t1$y),
 #              t(t1$y.ref), t(t1$x.ref), nrow(t1$y.ref), t1$myrts, t1$mybase,
 #              0, print.level = print.level)
 
 
 # redefine if some Farrell measures are not computed
 
 te.good <- !is.na(eff$te)
 K  <- sum(te.good) # number of obs
 if(K == 0){
  # cat.print(eff$te)
  stop("Could not compute measure of technical efficiency for a single data point")
 }
 eff$te <- eff$te[te.good]
 Y  <- Y[, te.good, drop = FALSE]
 X  <- X[, te.good, drop = FALSE]
 esample[!te.good] <- FALSE
 
 # return("done")
 
 # cat.print(subsampling)
 
 # return("done")
 
 if(subsampling){

  msub  <- round(Kr^kappa)
  
  # B <- reps
  
  if(!is.null(seed)) set.seed(seed)
  # set.seed(871635)
  
  eff.boot <- matrix(NA, nrow = K, ncol = reps)
  
  if(show.progress)
  {
    cat("\n")
    pb <- utils::txtProgressBar(min=0, max=reps, style=3)
  }
  
  for(b in seq_len(reps)){
   
   # gnerate Data
   newSample <- sample( seq_len(Kr), msub, replace = TRUE)
   # if(b<3){
   #   if(b==1){
   #     cat("\n")
   #   }
   #   cat.print(newSample)
   # }
   ystar <- Yr[, newSample, drop = FALSE]
   xstar <- Xr[, newSample, drop = FALSE]
   
   # Efficiency relative to the bootstrapped frontier
   
   eff.boot[,b] <- .teNonrad2(Y, X, M, N, K, ystar, xstar, msub, rt, ba,
                              ifqh = FALSE, print.level = 0,
                              lmdConstr = TRUE, full.solution = FALSE)
   # if(b<3) cat.print(eff.boot[,b])
   # if(b==3) break
   # if(my.base == 1){
   #  # input base
   #  # eff.boot[,b] <- nonradial.c(YOBS = t(Y), XOBS = t(X), YREF = ystar, XREF = xstar, RTS = rt, base = ba, lmdConstr = TRUE, print.level = print.level)$te
   #  
   #  eff.boot[,b] <- .teNonrad2(Y, X, M, N, K, ystar, xstar, msub, rt, ba,
   #             ifqh = FALSE, print.level = 0,
   #             lmdConstr = TRUE, full.solution = FALSE)$te
   # } else {
   #  # output base
   #  eff.boot[,b] <- nonradial.c(YOBS = t(Y), XOBS = t(X), YREF = ystar, XREF = xstar, RTS = rt, base = ba, lmdConstr = TRUE, print.level = print.level)$te
   # }
   
   if(show.progress)  utils::setTxtProgressBar(pb, b)
   }
 } else {
   # cat.print(101)
  # my.eff <- cbind( eff$te, eff$te.detail )
  # colnames(my.eff) <- c("RM",paste0("lambda_",1:ncol(eff$te.detail)))
  # head(my.eff)
  # cat.print(102)
  
  # cor(t(eff$te.all))
  
  # 0
  E <- eff$te.detail
  # K <- ncol(Xobs) # number of obs
  # N <- nrow(Xobs) # number of Xs
  # M <- nrow(Yobs) # number of Ys
  # M <- ncol(E)
  Q <- ifelse(my.base == 1, N, M)
  
  # cat.print(103)
  
  if(Q == 1){
   stop(paste0(ifelse(my.base == 1, "Input", "Output"),"-based nonradial efficiency measurerement is the same as ",ifelse(my.base == 1, "input", "output"),"-based radial efficiency measurerement: use 'teradial' instead"), call. = FALSE)
  }
  
  seq_len_Q <- seq_len(Q)
  
  # cat.print(104)
  
  E <- cbind(E, 2-E)
  # cat("dim(E) = ", dim(E),"\n")
  
  # Create grid of combinations of eff and reflected eff
  e.names <- list()
  for(i in seq_len(Q)){
   e.names[[paste0("e",i)]] <- c(1,2)
  }
  # @
  # cat("\nprint(e.names) \n")
  # print(e.names)
  # e.names
  my.grid <- expand.grid(e.names)
  # @
  # cat("\nprint(my.grid) \n")
  # print(my.grid)
  n.comb <- nrow(my.grid)
  
  # cat.print(105)
  
  #  Big matrix
  # define columns in E matrix:
  # if my.grid == 2, get reflected, otherwise normal eff
  Er <- NULL
  for(j in seq_len( n.comb )){
   # print(my.grid[j,])
   Er <- rbind(Er, E[,ifelse( my.grid[j,] - 1 == 1, Q + seq_len_Q, seq_len_Q) ] )
   # cat("j = ", j, "\n")
   # print(ifelse( my.grid[j,] - 1 == 1, Q + seq_len_Q, seq_len_Q))
  }
  
  # cat.print(106)
  
  # dim(Er)
  # print(dim(Er))
  # print(Er)
  # return(my.grid)
  
  #  Create list with VC and Chol for each row of "my.grid"
  VCs <- Cs <- list()
  for(j in seq_len( n.comb )){
   VCs[[paste0("vc",j)]] <- cov( E[,ifelse( my.grid[j,] - 1 == 1, Q + seq_len_Q, seq_len_Q), drop = FALSE ] )
   Cs[[paste0("c",j)]] <- t( chol( VCs[[paste0("vc",j)]]) )
  }
  
  # cat.print(107)
  
  #  work with intervals
  
  my.intervals <- seq_len( n.comb) * K
  
  # .findInterval(x = 120, my.intervals)
  
  h <- (4/(5 * K))^(1/6)
  h <- (4/(Q+2)) ^ (1/(4+Q))
  h <- (4/((Q+2) * K)) ^ (1/(Q+4))
  # B <- reps
  
  if(!is.null(seed)) set.seed(seed)
  # set.seed(871635)
  
  eff.boot <- matrix(NA, nrow = K, ncol = reps)
  
  # cat.print(108)
  
  if(show.progress)
  {
    cat("\n")
    pb <- utils::txtProgressBar(min=0, max=reps, style=3)
  }
  
  for(b in seq_len(reps)) {
   # Sample
   ii <- sample(seq_len( n.comb*Kr ) , Kr)
   from.inte <- .findInterval(x = ii, my.intervals)
   
   # From Big matrix Er
   E1 <- Er[ii,]
   E1bar <- matrix( rep(colMeans(E1), each = Kr), ncol = Q)
   
   # Draw from multivariate normal with known VC
   e <- matrix(rnorm(Q * Kr), nrow = Kr, ncol = Q)
   eps <- matrix(NA, nrow = Kr, ncol = Q)
   for(i in seq_len(n.comb)){
    from.inte.i <- from.inte == i
    eps[from.inte.i, ] <- e[from.inte.i, ] %*% Cs[[paste0("c",j)]]
   }
   
   E2. <- sqrt(1 + h^2) * ( E1bar + (E1 - E1bar + h * eps) )
   
   # Reflect back if necessary
   use.Farrell <- FALSE
   if(my.base == 1){
    # input base
    E2 <- ifelse(E2. < 1, E2., 2 - E2.)
    if( min(E2) < 0){
     use.Farrell <- TRUE
     E2 <- ifelse(E2. > 1, E2., 2 - E2.)
    }
   } else {
    # output base
    E2 <- ifelse(E2. > 1, E2., 2 - E2.)
   }
   
   # transform Data
   if(my.base == 1){
    # input base
     # cat.print(dim(Xr))
     # cat.print(dim(eff$te.detail))
     # cat.print(dim(E2))
     
    if(use.Farrell){
     Xobs.star <- Xr * t( eff$te.detail * E2 )
    } else {
     Xobs.star <- Xr * t( eff$te.detail / E2 )
    }
     
   } else {
    # output base
    Yobs.star <- Yr * t(eff$te.detail / E2)
   }
   
   # if(b < 3) print(summary(eff$te.detail * E2))
   
   # Efficiency relative to the bootstrapped frontier
   if(my.base == 1){
    # input base
    # eff.boot[,b] <- nonradial.c(YOBS = t(Y), XOBS = t(X), YREF = t(Yr), XREF = Xobs.star, RTS = rt, base = ba, lmdConstr = TRUE, print.level = print.level)$te
    eff.boot[,b] <- .teNonrad2(Y, X, M, N, K, Yr, Xobs.star, Kr, rt, ba,
                                ifqh = FALSE, print.level = 0,
                                lmdConstr = TRUE, full.solution = FALSE)
   } else {
    # output base
    # eff.boot[,b] <- nonradial.c(YOBS = t(Y), XOBS = t(X), YREF = Yobs.star, XREF = t(Xr), RTS = rt, base = ba, lmdConstr = TRUE, print.level = print.level)$te
    eff.boot[,b] <- .teNonrad2(Y, X, M, N, K, Yobs.star, Xr, Kr, rt, ba,
                               ifqh = FALSE, print.level = 0,
                               lmdConstr = TRUE, full.solution = FALSE)
   }
   
   if(show.progress)  utils::setTxtProgressBar(pb, b)
  }
 }
 cat("\n")
 
 # cat.print(109)
 
 # return("done")
 
 
 # CI
 
 # Working with large BOOT matrix
 # Step 3: bias correction and CIs
 bci <- .biasAndCI(te=eff$te, teboot=t(eff.boot), msub = msub, K=K, Kr=K, M=M, N=N,
                   level=level, smoothed=!subsampling, forceLargerOne = FALSE)
 # cat.print(110)
 # cat.print(21)
 
 # Check if any bias corrected measure is negative;
 #	If yes, do inference in terms of Shephard
 #	distance functions (must be input oriented)
 if( min(bci$tebc, na.rm = TRUE) < 0 ){
  bci <- .biasAndCI(te=eff$te, t(eff.boot), msub = msub, K, Kr, M, N,
                    level, smoothed=!subsampling, forceLargerOne = TRUE)
  eff$te <- 1/eff$te
  if(print.level >= 1){
   warning(" One or more bias-corrected Russell ",YX$base.string,"-based measures of technical efficiency is negative.  Analysis will be done in terms of Shephard distance function, a reciprocal of the Farrell measure. \n", call. = FALSE, immediate. = FALSE)
   my.warnings <- c(my.warnings, paste0(" One or more bias-corrected Russell ",YX$base.string,"-based measures of technical efficiency is negative.  Analysis will be done in terms of Shephard distance function, a reciprocal of the Farrell measure."))
  }
 }
 
 # cat.print(111)
 
 # Check if bias squared over var is larger for all
 if( min(bci$BovV, na.rm = TRUE) < 1){
  warning(" For one or more data points statistic [(3*(bias)^2/var] is smaller than 1; bias-correction should not be used.\n", call. = FALSE, immediate. = FALSE)
  my.warnings <- c(my.warnings, " For one or more data points statistic [(3*(bias)^2/var] is smaller than 1; bias-correction should not be used.")
 }
 
 # cat.print(112)
 
 tymch <- list(call = match.call(), model ="tenonradialbc", K = K, M = M, N = N, reps = reps, level = level,
                rts = YX$rts.string, base = YX$base.string, 
                Kr = Kr, ref = ifelse(is.null(ref), "FALSE", "TRUE"),
                te = eff$te, tebc = bci$tebc, biasboot = bci$bias, varboot = bci$vari,
                biassqvar = bci$BovV, realreps = bci$reaB,
                telow = bci$LL, teupp = bci$UU, teboot = t(eff.boot),
                esample = esample, esample.ref = YX$esample.ref, warnings = my.warnings)
 # cat.print(113)
 class(tymch) <- "npsf"
 return(tymch)
 
}



#