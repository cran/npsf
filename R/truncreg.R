truncreg <- function(formula, data, subset, 
                     ll = -Inf, ul = Inf, 
                     lmtol = .Machine$double.eps, maxiter = 150, 
                     marg.eff = FALSE, 
                     print.level = 1) {

 YX <- .prepareYXllul(formula = formula, ll = ll, ul = ul, data, subset, sysnframe = sys.nframe())
 y <- YX$Y
 x <- YX$X
 a <- YX$LL
 b <- YX$UL
 # sample size
 n <- YX$n
 n.full <- YX$n.full
 if(print.level >= 1){
  cat("\n(Note: ",n.full-n," (out of ",n.full,") observations truncated) \n\n", sep = "")
 }
 
 # initial values
 beta0 <- 0.999999*solve(t(x) %*% x) %*% t(x) %*% y
 sigma0 <- 1.000001*sqrt( sum((y - x %*% beta0)^2)/(n-ncol(x)) )
 names(sigma0) <- "/sigma"
 theta0 <- rbind( beta0, sigma0)
 # initial LL
 l_init <- .ll.trunc (y, x, beta0, sigma0, a, b)
 a0 <- a
 b0 <- b
 
 # print(summary(a))
 # print(summary(b))
 
 Delta <- 1
 mylm <- 1
 step <- 0
 time.05 <- proc.time()
 while ( abs(mylm) > lmtol & step <= maxiter){
  # print(step)
  step <- step + 1
  # print(step)
  XBeta <- x %*% beta0
  aXBeta <- a - XBeta
  bXBeta <- b - XBeta
  yXBeta <- y - XBeta
  h0 <- .h.trunc(al = aXBeta/sigma0, be = bXBeta/sigma0, a1 = 1, b1 = 1, a = a0, b = b0 )
  h1 <- .h.trunc(al = aXBeta/sigma0, be = bXBeta/sigma0, a1 = aXBeta, b1 = bXBeta, a = a0, b = b0 )
  h2 <- .h.trunc(al = aXBeta/sigma0, be = bXBeta/sigma0, a1 = aXBeta^2, b1 = bXBeta^2, a = a0, b = b0 )
  h3 <- .h.trunc(al = aXBeta/sigma0, be = bXBeta/sigma0, a1 = aXBeta^3, b1 = bXBeta^3, a = a0, b = b0 )
  h0 <- as.matrix(h0)
  h1 <- as.matrix(h1)
  h2 <- as.matrix(h2)
  h3 <- as.matrix(h3)
  # gradient
  dl <- rbind(
   sigma0 ^ (-2) * t(x) %*% yXBeta + sigma0 ^ (-1) * t(x) %*% h0,
   -n * sigma0 ^ (-1) + sigma0 ^ (-3) * sum(yXBeta^2) + sigma0 ^
    (-2) * sum(h1)
  )
  # colnames(dl) <- "gradient"
  # print(dl)
  # Hessian
  # ddlA11 <- -sigma0^(-2)*t(x) %*% (diag( 1 - as.vector(sigma0^(-1) * h1 + h0^2) ) ) %*% x
  ddlA11 <- (-sigma0^(-2)*t(x * (1 - as.vector(sigma0^(-1) * h1 + h0^2) ) )  ) %*% x
  ddlA12 <- t(x) %*% (-2*sigma0^(-3)* yXBeta - sigma0^(-2)*h0 + sigma0^(-4)*h2 + sigma0^(-3)*h0*h1)
  ddlA21 <- t(-2*sigma0^(-3)* yXBeta - sigma0^(-2)*h0 + sigma0^(-4)*h2 + sigma0^(-3)*h0*h1) %*% x
  ddlA22 <- n*sigma0^(-2) - 3*sigma0^(-4)* sum(yXBeta^2) - 2*sigma0^(-3) * sum(h1) + sigma0^(-5) * sum(h3) + sigma0^(-4) * sum(h1^2)
  ddl <- rbind( cbind(ddlA11, ddlA12), cbind(ddlA21, ddlA22) )
  #print(ddl)
  iddl <- solve(ddl)
  theta0 <- theta0 - iddl %*% dl
  #print(theta0)
  beta0 <- theta0[-length(theta0), , drop = FALSE]
  sigma0 <- theta0[length(theta0)]
  # sig0 <- sqrt(-diag(iddl))
  #print(cbind(dl, theta0))
  likelihood <- .ll.trunc (y,x, beta0, sigma0,a,b)
  # print(likelihood, digits = 4)
  mylm <- as.vector( t(dl)%*%-iddl%*%dl )
  if ( is.na(mylm) ) stop("cannot compute an improvement")
  # cat(" Iteration ",step,": LM criterion = ",mylm,", LL = ",l_init," \n", sep = "")
  
  if ( print.level >= 1) {
   if ( print.level >= 2){
    cat(paste("Iteration ",formatC(step, width = 2)," (hessian is analytical, g*inv(H)*g' = ",formatC(mylm, format = "E", digits = 1),"): log likelihood = ",format(l_init, nsmall = 5),"\n", sep = ""), sep = "")
   } else {
    cat(paste("Iteration ",formatC(step, width = 2),": log likelihood = ",format(l_init, nsmall = 5),"\n", sep = ""), sep = "")
   }
  }
  Delta <- abs(likelihood - l_init)
  l_init <- likelihood
 } # the end of while
 
 time.06 <- proc.time()
 est.time.sec <- (time.06-time.05)[1]
 names(est.time.sec) <- "sec"
 if(print.level >= 1){
  .timing(est.time.sec, "\nLog likelihood maximization completed in\n")
  # cat("___________________________________________________\n")
 }
 
 if(step > maxiter) cat(" \n Not converged in ",maxiter," iterations \n\n", sep = "")
 # cat(" \n", sep = "")

 if ( print.level >= 1){
  cat("\nConvergence given g*inv(H)*g' = ",formatC(abs(mylm), format = "e", digits = 3)," < lmtol (",lmtol,")\n", sep = "")
  cat(paste("\nFinal log likelihood = ",formatC(l_init,digits=7,format="f"),"\n\n", sep = ""), sep = "")
  cat("Truncated regression for cross-sectional data\n", sep = "")
  cat("Limits:\n", sep = "")
  #  check if infinite
  ll.inf <- sum(is.infinite(YX$LL))
  ul.inf <- sum(is.infinite(YX$UL))
  
  if (ll.inf == n) {
   #  if all are infinite
   cat("lower limit for left-truncation  = ",YX$LL[1],"\n", sep = "")
  } else {
   if (sd(YX$LL[is.finite(YX$LL)]) != 0) {
    cat("lower limit for left-truncation:", sep = "")
    print(summary(YX$LL))
   } else {
    cat("lower limit for left-truncation  = ",YX$LL[1],"\n", sep = "")
   }
  }
  
  if (ul.inf == n) {
   #  if all are infinite
   cat("upper limit for right-truncation = ",YX$UL[1],"\n", sep = "")
  } else {
   if (sd(YX$UL[is.finite(YX$UL)]) != 0) {
    cat("upper limit for right-truncation:", sep = "")
    print(summary(YX$UL))
   } else {
    cat("upper limit for right-truncation = ",YX$UL[1],"\n", sep = "")
   }
  }
  
  cat("Number of observations (used in regression)               = ",n,"\n", sep = "")
  cat("Number of truncated observations (not used in regression) = ",n.full-n,"\n", sep = "")
  cat("\nEstimation results:\n\n", sep = "")
 }
 

 sig0 <- sqrt(-diag(iddl))

 # marginal effects
 # lmd <- NULL
 if(marg.eff){
  m.eff <- matrix(nrow = n, ncol = length(beta0))
  for(qq in 1:n)
  {
   x0 <- x[qq, ,drop = FALSE]
   a0 <- ifelse(length(a)==1, a, a[qq])
   b0 <- ifelse(length(b)==1, b, b[qq])

   if(b0 == Inf)
   {
    ai <- (a0 - x0 %*% beta0) / sigma0
    lambdai <- dnorm(ai) / (1 - pnorm(ai))
    # print(c(ai,lambdai))
    m.eff[qq,] <- as.matrix(1 - lambdai^2 + ai*lambdai) %*% t(beta0)
   } else if (a0 == -Inf)
   {
    ai <- (b0 - x0 %*% beta0) / sigma0
    lambdai <- -dnorm(ai) / pnorm(ai)
    m.eff[qq,] <- as.matrix(1 - lambdai^2 + ai*lambdai) %*% t(beta0)
   } else
   {
    ai1 <- (a0 - x0 %*% beta0) / sigma0
    ai2 <- (b0 - x0 %*% beta0) / sigma0
    lambdai <- (dnorm(ai2)-dnorm(ai1)) / (pnorm(ai2) - pnorm(ai1))
    m.eff[qq,] <- as.matrix(1 - lambdai^2 + ai2*lambdai - (b0-a0)*dnorm(ai1) / sigma0 / ( pnorm(ai2) - pnorm(ai1) )) %*% t(beta0)
   }
   # lmd <- c(lmd, lambdai)
  } # end of qq loop
  colnames(m.eff) <- rownames(beta0)
 }

 coeffs <- cbind(theta0, sig0,theta0/sig0,(1-pnorm(abs(theta0/sig0),lower.tail =TRUE))*2)
 rownames(coeffs)[length(beta0)+1] <- "/sigma"
 colnames(coeffs) <- c("Estimate","Std. Error","z","Pr(>|z|)")

 if (print.level >= 1){
  printCoefmat(coeffs)
 }

 tymch <- list(call = match.call(), model = "truncreg.cs", coef = coeffs[,1], table = coeffs, vcov = -iddl, ll = l_init, lmtol = lmtol, LM = mylm, esttime = est.time.sec, sigma = theta0[length(theta0)], LL = a, UL = b, n = n, n.full = n.full, nontruncsample = as.vector(YX$nontruncsample), esample = YX$esample)
 
 if(marg.eff){
  tymch$marg.effects <- m.eff
 }

 class(tymch) <- "npsf"
 return(tymch)
}

