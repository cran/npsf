sf <- function(formula, data, it = NULL, subset,
               prod = TRUE, model = "K1990", distribution = c("h"),
               eff.time.invariant = FALSE, 
               mean.u.0i.zero     = FALSE,
               mean.u.0i          = NULL,
               ln.var.u.0i        = NULL,
               ln.var.v.0i        = NULL,
               ln.var.v.it        = NULL,  
               simtype = c("halton", "rnorm"), halton.base = NULL, R = 500,
               simtype_GHK = c("halton", "runif"), R_GHK = 500,
               random.primes = FALSE, 
               cost.eff.less.one  = FALSE, level = 95, marg.eff = FALSE,
               start.val = NULL, maxit = 199, report.ll.optim = 10, 
               reltol = 1e-8, lmtol = sqrt(.Machine$double.eps),
               digits = 4, print.level = 4, seed = 17345168,
               only.data = FALSE) {
 
 # prefixes for the names will be:
 # lnVARu0i_
 # lnVARuit_
 # lnVARv0i_
 # lnVARvit_
 # mean_u0i_
  
 myprod <- ifelse(prod, 1, -1)
  
 if(!is.null(halton.base)){
   if(length(halton.base) != 1){
     stop("Length of 'halton.base' must be one")
   }
 }
  
 if (level < 0 | level > 99.99) {
  stop("'level' must be between 0 and 99.99 inclusive", call. = FALSE)
 }
 alpha <- 1 - level / 100
 
 if(is.null(it)){
  
  # Cross-sectional models begin --------------------------------------------
  
  warning("This is a cross-sectional model", immediate. = TRUE, call. = FALSE)
  
  if(length(distribution) != 1){
   stop("Distribution of inefficiency term should be specified.")
  } else {
   distribution <- tolower(substr(distribution, 1,1 ))
  }
  
  if( !distribution %in% c("t","h","e") ){
   stop("'distribution' is invalid")
  }
  
  if(is.null(mean.u.0i) == FALSE & distribution != "t"){
   stop("Option 'mean.u.0i' can be used only when distribution of inefficiency term is truncated normal.")
  }
  
  YXZ <- .prepareYXZ.cs(formula = formula, ln.var.u.0i = ln.var.u.0i, ln.var.v.0i = ln.var.v.0i, mean.u.0i = mean.u.0i, data, subset, sysnframe = sys.nframe())
  
  y <- YXZ$Y
  X <- YXZ$X
  Zu <- YXZ$Zu
  Zv <- YXZ$Zv
  Zdel <- YXZ$Zdel
  n <- YXZ$n
  k <- YXZ$k
  ku <- YXZ$ku
  kv <- YXZ$kv
  kdel <- YXZ$kdel
  esample <- YXZ$esample
  
  abbr.length <- 34
  names_x <- abbreviate(colnames(X), abbr.length+6, strict = TRUE, dot = FALSE, method = "both.sides")
  # names_zu <-  abbreviate(colnames(Zu), 9, strict = TRUE, dot = FALSE)
  # names_zv <-  abbreviate(colnames(Zv), 9, strict = TRUE, dot = FALSE)
  Zu_colnames <- abbreviate(colnames(Zu), abbr.length, strict = TRUE, dot = FALSE, method = "both.sides")
  names_zu <- paste0("lnVARu0i_", abbreviate(colnames(Zu), abbr.length, strict = TRUE, dot = FALSE, method = "both.sides"))
  names_zv <- paste0("lnVARv0i_", abbreviate(colnames(Zv), abbr.length, strict = TRUE, dot = FALSE, method = "both.sides"))
  
  if(distribution == "t"){
   # names_del <- abbreviate(colnames(Zdel), 9, strict = T, dot = F)
   kdel <- ncol(Zdel)
   names_del <- paste0("mean_u0i_", abbreviate(colnames(Zdel), 9, strict = TRUE, dot = FALSE))
   if(kdel == 1){
    # if(names_del == "mu_intercept"){
    #  names_del <- "mu"
    # }
   }
  } else {
   names_del <- NULL
   kdel <- 0
  }
  
  if(is.null(names_del)){
   coef.names.full <- c(
    names_x,
    paste("lnVARv0i_",c("(Intercept)", names_zv[-1]),"", sep = ""),
    paste("lnVARu0i_",c("(Intercept)", names_zu[-1]),"", sep = "")
   )
  } else {
   coef.names.full <- c(
    names_x,
    paste("lnVARv0i_",c("(Intercept)", names_zv[-1]),"", sep = ""),
    paste("lnVARu0i_",c("(Intercept)", names_zu[-1]),"", sep = ""),
    paste("mean_u0i_",c("(Intercept)", names_del[-1]),"", sep = "")
   )
  }
  
  coef.names.full.output <- c(names_x, names_zv, names_zu, names_del)
  
  # starting values
  ols <- lm(y ~ 0 + X)
  # print(summary(ols))
  beta0 <- as.matrix(coef(ols), nrow = length(coef(ols)), ncol = 1)
  olsResid <- y - X%*%beta0
  olsSkewness <- .skewness(olsResid)
  
  # Moment estimators for sigma_u/v squared
  if(is.null(start.val) == TRUE){
   m3 <- sum(olsResid^3)/n
   m3 <- ifelse(m3*myprod < 0, m3, -0.0001*myprod)
   m2 <- sum(olsResid^2)/n
   su2init <- ((m3/(sqrt(2/pi)*(1-4/pi)))^2)^(1/3)
   if(is.nan(su2init) | su2init < 0) su2init <- 0.1
   sv2init <- m2 - (1 - 2/pi)*su2init
   if(is.nan(sv2init) | sv2init < 0) sv2init <- 0.1
   ou2 <- (-myprod*m3/2)^(2/3)
   beta0[1] <- beta0[1] + sqrt(2/pi)*sqrt(su2init)
   # theta0[1] <- theta0[1] + myprod*(sqrt(2/pi))*sqrt(ou2)
   y1 <- 0.5*log(((olsResid^2 - sv2init)/(1 - 2/pi))^2)
   reg_hetu <- lm(y1 ~ 0 + Zu)
   gu0 <- as.matrix(coef(reg_hetu), nrow = length(coef(reg_hetu)), ncol = 1)
   y2 <- 0.5*log((olsResid^2 - (1 - 2/pi)* su2init)^2)
   reg_hetv <- lm(y2 ~ 0 + Zv)
   gv0 <- as.matrix(coef(reg_hetv), nrow = length(coef(reg_hetv)), ncol = 1)
   if((is.null(mean.u.0i)) & (distribution == "t")){
    theta0 <- rbind( beta0, gv0, gu0, 0)
   } else if (!is.null(mean.u.0i)) {
    theta0 <- rbind( beta0, gv0, gu0, as.matrix(rep(0, kdel)))
   } else {
    theta0 <- rbind( beta0, gv0, gu0)
   }
   if(distribution == "e"){
     olsres2 <- olsResid^2; m2 = mean(olsres2)
     olsres3 <- olsResid^3; m3 = mean(olsres3)
     m3t     <- m3/sqrt(6*m2^3/n)
     # /* negative skewness, so one-side test */
     p_m3t   <- pnorm(myprod*m3t)
     # /* correct 3rd moment to be negative */
     m3      <- ifelse(m3*myprod < 0, m3, -0.0001*myprod)
     theta0  <- coef(ols)
     ou2 <- (-myprod*m3/2)^(2/3)
     lnou2 <- log(ou2)
     ov2 <- m2 - ou2
     lnov2 <- ifelse(ov2>0, log(ov2), log(.0001))
     theta0[1] <- theta0[1] + myprod*(sqrt(2/pi))*sqrt(ou2)
     
     # log SV2
     if(kv == 1){
       theta0 <- c(theta0, lnov2)
     } else {
       # cat.print(cbind(Zv[,-1,drop = FALSE],1))
       theta0 <- c(theta0, coef(lm( rep(lnov2,n) ~ Zv[,-1,drop = FALSE])))
       # cat.print(coef(lm( rep(lnov2,n) ~ Zv[,-1,drop = FALSE])))
       # cat.print(theta0)
     }
     # log SU2
     if(ku == 1){
       theta0 <- c(theta0, lnou2)
     } else {
       theta0 <- c(theta0, coef(lm( rep(lnou2,n) ~ Zu[,-1,drop = FALSE] )))
       # cat.print(coef(lm( rep(lnou2,n) ~ Zu[,-1,drop = FALSE] )))
       # cat.print(theta0)
     }
     theta0 <- as.matrix(theta0)
   }
  } else {
   theta0 <- start.val
  }
  rownames(theta0) <- coef.names.full#c(names_x, names_zv, names_zu, names_del)
  
  # cat.print(theta0)

  max.name.length <- max(nchar(rownames(theta0)))
  
  # print(theta0)

  # Optimization -----------------------------------------------------------

  time.05 <- proc.time()
  if(print.level >= 2){
   # cat("__________________________________________________")
   cat("\nOptimization using 'mlmaximize': \n\n", sep = "")
  }
  obj <- eval(parse(text = paste(" tryCatch(.mlmaximize(theta0, .ll.",noquote(distribution),"n, gr = .g.",noquote(distribution),"n, hess = .hess.",noquote(distribution),"n, prod = prod, k = k, kv = kv, ku = ku, y = y, Zv = Zv, Zu = Zu, kdel = kdel, Zdel = Zdel, X = X, print.level = print.level, lmtol = lmtol), error = function(e) e )", sep = "")))
  # print(obj)
  
  if(inherits(obj, "error")){
    print(obj)
   if(print.level >= 2){
    cat("",rep("_", max.name.length+42-1),"", "\n", sep = "")
    cat(" Failed, trying to optimize using 'optim':\n", sep = "")
    cat(" ('optim' minimizes, so 'll' is\n", sep = "")
    cat("  negative of what is printed) \n\n", sep = "")
   }
   obj <- eval(parse(text = paste(" tryCatch( optim(par = theta0, fn = .ll.",noquote(distribution),"n, gr = .g.",noquote(distribution),"n, prod = prod, k = k, kv = kv, ku = ku, y = y, Zv = Zv, Zu = Zu, kdel = kdel, Zdel = Zdel, X = X, method = c('BFGS'), control = list(reltol  = .Machine$double.eps, maxit = 1000, trace = ifelse(print.level >=2, 1, 0), REPORT = ifelse(print.level >=2, 1, 1), fnscale = -1), hessian = TRUE), error = function(e) e )", sep = "")))
   # print(obj)
   cannot.est.model <- inherits(obj, "error")
   if(cannot.est.model){
    # what is the error?
    # print(obj)
    stop("convergence was not reached or 'optim' cannot estimate the model at provided starting values")
   }
   if(obj$convergence == 0){
    convergence <- 0
    obj$gHg <- 17
    obj$ll <- obj$value
    obj$vcov = solve(-obj$hessian)
   } else if(cannot.est.model){
    # what is the error?
    # print(obj)
    warning("convergence was not reached or 'optim' cannot estimate the model at provided starting values")
   } else {
    theta0 <- obj$par
    cat("\n")
    if(print.level >= 2){
     cat("\n__________________________________________________")
     cat("\nOptimization using 'mlmaximize': \n\n", sep = "")
    }
    obj <- eval(parse(text = paste(" tryCatch(.mlmaximize(theta0, .ll.",noquote(distribution),"n, gr = .g.",noquote(distribution),"n, hess = .hess.",noquote(distribution),"n, prod = prod, k = k, kv = kv, ku = ku, y = y, Zv = Zv, Zu = Zu,kdel = kdel, Zdel = Zdel, X = X, print.level = print.level), error = function(e) e )", sep = "")))
    cannot.est.model <- inherits(obj, "error")
    if ( cannot.est.model ){
     # what is the error?
     # print(obj)
     stop("Cannot estimate the model at provided starting values")
    } else {
     convergence <- obj$conv.crite
    }
   }} else {
    convergence <- obj$conv.crite
   }
  
  time.06 <- proc.time()
  est.time.sec <- (time.06-time.05)[1]
  names(est.time.sec) <- "sec"
  if(print.level >= 2){
   .timing(est.time.sec, "\nLog likelihood maximization completed in\n")
   # cat("___________________________________________________\n")
  }
  
  # SE for sigmau and sigmav
  sigv <- sqrt(exp(obj$par[k+1]))
  sigu <- sqrt(exp(obj$par[k+kv+1]))
  se_sigv <- (0.5*sigv) * sqrt(obj$vcov[k+1,k+1])
  se_sigu <- (0.5*sigu) * sqrt(obj$vcov[k+kv+1,k+kv+1])
  
  # SE for lambda
  lmd <- sqrt(exp(obj$par[k+kv+1] - obj$par[k+1]))
  glmd <- c(0.5*lmd, -0.5*lmd)
  se_lmd <- as.vector( sqrt( t(glmd) %*% obj$vcov[c((k+1),(k+kv+1)), c((k+1),(k+kv+1))] %*% glmd))
  
  # SE for gamma
  gam <- exp(obj$par[k+kv+1])/(exp(obj$par[k+kv+1]) + exp(obj$par[k+1]))
  ggam <- c(gam*(1 - gam), -gam^2/lmd^2)
  se_gam <- as.vector( sqrt( t(ggam) %*% obj$vcov[c((k+1),(k+kv+1)), c((k+1),(k+kv+1))] %*% ggam))
  
  if((is.null(ln.var.u.0i) == TRUE) & (is.null(ln.var.v.0i) == TRUE)){
   par <- c(obj$par, sigv, sigu, lmd, gam); se <- c(sqrt(diag(obj$vcov)), se_sigv, se_sigu, se_lmd, se_gam)
   names(par) = c(names_x, names_zv, names_zu, names_del, "sigma_v  ", "sigma_u  ", "lambda   ", "gamma   ")
  } else if((is.null(ln.var.u.0i) == T) & (is.null(ln.var.v.0i) == F)){
   par <- c(obj$par, sigu); se <- c(sqrt(diag(obj$vcov)), se_sigu)
   names(par) = c(names_x, names_zv, names_zu, names_del, "sigma_u  ")
  } else if((is.null(ln.var.u.0i) == F) & (is.null(ln.var.v.0i) == T)){
   par <- c(obj$par, sigv); se <- c(sqrt(diag(obj$vcov)), se_sigv)
   names(par) = c(names_x, names_zv, names_zu, names_del, "sigma_v  ")
  } else {
   par <- c(obj$par); se <- c(sqrt(diag(obj$vcov)))
   names(par) = c(names_x, names_zv, names_zu, names_del)
  }
  
  output <- cbind(round(par, digits = digits), round(se, digits = digits), round(par/se,digits = 2), round(pnorm(abs(par/se), lower.tail = FALSE)*2, digits = digits))
  colnames(output) <- c("Coef.", "SE ", "z ",  "P>|z|")
  
  # estimation results
  max.name.length <- max(nchar(row.names(output)))
  
  if(print.level >= 1){
   cat("",rep("_", max.name.length+42-1),"", "\n", sep = "")
   if(prod){
    cat("\nCross-sectional stochastic (production) frontier model\n",  sep = "")}
   else {
    cat("\nCross-sectional stochastic (cost) frontier model\n",  sep = "")
   }
   cat("\nDistributional assumptions\n\n", sep = "")
   Assumptions <- rep("heteroskedastic",2)
   if(kv==1){
    Assumptions[1] <- "homoskedastic"
   }
   if(ku==1){
    Assumptions[2] <- "homoskedastic"
   }
   Distribution = c("normal ", "half-normal ")
   if(distribution == "t"){
    Distribution[2] <- "truncated-normal "
   }   
   if(distribution == "e"){
    Distribution[2] <- "exponential "
   }
   a1 <- data.frame(
    Component = c("Random noise: ","Inefficiency: "),
    Distribution = Distribution,
    Assumption = Assumptions
   )
   print(a1, quote = FALSE, right = FALSE)
   cat("\nNumber of observations = ", n, "\n", sep = "")
   max.name.length <- max(nchar(row.names(output)))
   est.rez.left <- floor( (max.name.length+42-22) / 2 )
   est.rez.right <- max.name.length+42-22 - est.rez.left
   cat("\n",rep("-", est.rez.left)," Estimation results: ",rep("-", est.rez.right),"\n\n", sep ="")
   # cat("\n--------------- Estimation results: --------------\n\n", sep = "")
   .printoutcs(output, digits = digits, k = k, kv = kv, ku = ku, kdel = kdel, na.print = "NA", dist = distribution, max.name.length = max.name.length)
  }
  
  # Technical efficiencies  -------------------------------
  if(distribution == "e"){
    sn1 <- ifelse(prod, 1, -1)
    if(!prod & cost.eff.less.one) sn1 <- 1
    
    xbeta <- X %*%    obj$par[1                        :k]
    eps0 <-       y - xbeta
    sv   <- sigmas_v <- exp( 0.5 * Zv %*% obj$par[(k + 1)                  :(k + kv)])
    su   <- sigmas_u <- exp( 0.5 * Zu %*% obj$par[(k + kv  + 1)     :(k + kv + ku)] )
    
    mu1  <- -myprod *eps0 - (sv^2/su)
    s1   <- sv
    z    <- mu1/s1
    
    # u = E(u|e)
    u_JMLS  <- mu1 + s1 * dnorm(-z)/pnorm(z)
    # colnames(u_JMLS) <- 'u_JMLS'
    te_jlms_mean <- exp(-sn1*u_JMLS)
    # colnames(te_JMLS) <- 'te_JMLS'
    
    # u = M(u|e)
    m_JMLS  <- ifelse(mu1 >= 0, mu1, 0)
    # colnames(m_JMLS) <- 'm_JMLS'
    te_jlms_mode <- exp(-sn1*m_JMLS)
    # colnames(te_JMLS_m) <- 'te_JMLS_m'
    
    # TE
    te_bc   <- (pnorm(-sn1*s1+z)) /pnorm(z)  * exp(-sn1*mu1+1/2*s1^2)
    colnames(te_bc) <- 'te_bc'
    
    zl    <- qnorm( 1 - alpha / 2 * pnorm(z) )
    zu    <- qnorm( 1 - ( 1 - alpha/2 ) * pnorm(z) )
    te_l  <- exp( -sn1*mu1 - sn1* zl * s1 )
    te_u  <- exp( -sn1*mu1 - sn1* zu * s1 )
    eff <- data.frame(te_u, te_jlms_mean, te_jlms_mode,te_bc,te_l)
    colnames(eff) <- c("Lower bound","JLMS", "Mode", "BC","Upper bound" )
    
  } else {
    xbeta <- X%*%obj$par[1:k,]
    e <- y - xbeta
    sigmas_v <- sqrt(exp(Zv%*%obj$par[(k+1):(k+kv)]))
    sigmas_u <- sqrt(exp(Zu%*%obj$par[(k+kv+1):(k+kv+ku)]))
    
    if(distribution == "h"){
      # eff <- round(.u2efftnm(e, sigmas_u, sigmas_v, mu = 0, alpha = alpha, prod = prod), digits = digits)
      eff <- .u2efftnm(e, sigmas_u, sigmas_v, mu = 0, alpha = alpha, prod = prod, cost.eff.less.one = cost.eff.less.one)
      mu = NULL
    } else {
      mu <- Zdel%*%obj$par[-c(1:(k+kv+ku))]
      # eff <- round(.u2efftnm(e, sigmas_u, sigmas_v, mu = mu, alpha = alpha, prod = prod), digits = digits)
      eff <- .u2efftnm(e, sigmas_u, sigmas_v, mu = mu, alpha = alpha, prod = prod, cost.eff.less.one = cost.eff.less.one)
    }
  }
  # cat.print(eff)
  
  # Marginal effects
  if(marg.eff & distribution != "e"){
   meff = tryCatch(.me(obj$par[-c(1:(k+kv)),1,drop = FALSE], Zu = Zu, Zdel = Zdel, ku = ku, kdel = kdel, n = n, dist = distribution), error = function(e) { cat('In error handler\n'); print(e); e })
   if(inherits(meff, "error"))  {
    meff <- NULL
   } else {
    meff <- as.data.frame(meff)
    colnames(meff) <- colnames(Zu)[-1]
   }
  } else {
   meff <- NULL
  }
  
  myeff <- ifelse(prod, "technical", "cost")
  
  if(print.level >= 3){
   eff.name <- paste0("Summary of ",myeff," efficiencies")
   len.eff.name <- nchar(eff.name)
   est.eff.left <- floor( (max.name.length+42-len.eff.name-4) / 2 )
   est.eff.right <- max.name.length+42-len.eff.name-4 - est.eff.left
   cat("\n",rep("-", est.eff.left)," ",eff.name,": ",rep("-", est.eff.right),"\n\n", sep ="")
   eff1 <- eff[,2:4]
   colnames(eff1) <- formatC(colnames(eff1), width = 4, flag = "-")
   cat("",rep(" ", est.eff.left+1),"JLMS:= exp(-E[ui|ei])\n", sep = "")
   cat("",rep(" ", est.eff.left+1),"Mode:= exp(-M[ui|ei])\n", sep = "")
   cat("",rep(" ", est.eff.left+3),"BC:= E[exp(-ui)|ei]\n", sep = "")
   cat("\n")
   .su1(eff1, transpose = TRUE)
   
   # cat("\n=================")
   # cat(" Summary of ",myeff," efficiencies, exp(-E[ui|ei]): \n\n", sep = "")
   # .su1(eff1[,1, drop = FALSE], transpose = TRUE)
   #
   # cat("\n=================")
   # cat(" Summary of ",myeff," efficiencies, exp(-M[ui|ei]): \n\n", sep = "")
   # .su1(eff1[,2, drop = FALSE])
   #
   # cat("\n=================")
   # cat(" Summary of ",myeff," efficiencies, E[exp(-ui)|ei]: \n\n", sep = "")
   # .su1(eff1[,3, drop = FALSE])
   # cat("\n\n")
  }
  
  if(print.level > 6){
   cat("Point and interval estimates of unit-specific ",myeff," efficiencies: \n\n", sep = "")
   print(eff1)
  }
  
  if(length(unique(sigmas_u)) == 1) sigmas_u = NULL
  if(length(unique(sigmas_v)) == 1) sigmas_v = NULL
  if(distribution == "e"){
    mu = NULL
  } else {
    if(length(unique(mu)) == 1) mu = NULL
  }
  
  if(is.null(ln.var.u.0i) & is.null(ln.var.v.0i)){
   if (prod & olsSkewness > 0) {
    warning("OLS residuals are positively skewed, leading to zero estimate of inefficiency variance. Remaining ML estimates are equal to OLS. All units are estimated to be fully technically efficient.")
   } else if(!prod & olsSkewness < 0) {
    warning("OLS residuals are negatively skewed, leading to zero estimate of inefficiency variance. Remaining ML estimates are equal to OLS. All units are estimated to be fully cost efficient.")
   }
  }
  
  colnames(obj$vcov) <- rownames(obj$vcov) <- c(names_x, names_zv, names_zu, names_del)
  
  temp <- list(call = match.call(), model = "sf_sc", coef = par[names(par) %in% colnames(obj$vcov)], table = output, vcov = obj$vcov, loglik = obj$ll, lmtol = lmtol, LM = obj$gHg, convergence = convergence, olsSkewness = olsSkewness, esttime = est.time.sec, prod = prod, efficiencies = eff, marg.effects = meff, sigmas_u = sigmas_u, sigmas_v = sigmas_v, n = n, mu = mu, k = k, kv = kv, ku = ku, kdel = kdel, distribution = distribution, fitted = as.vector(unname(xbeta)), yobs = y, xobs = X, esample = esample)
  if (convergence == 2){
   temp$delta_rel <- obj$delta_rel
   temp$ltol <- obj$ltol
   
  }
  if (convergence == 3){
   temp$theta_rel_ch <- obj$theta_rel_ch
   temp$steptol <- sqrt(.Machine$double.eps)
  }
  
  class(temp) <- "npsf"
  
  return(temp)
  
  # Cross-sectional models end ----------------------------------------------
  
  
 } else {
  
  # Panel data models 1st and 2nd gen begin ---------------------------------
  
  warning("This is a panel data model", immediate. = TRUE, call. = FALSE)
  
  if(eff.time.invariant) model <- "K1990"
  
  model.correct <- model %in% c("K1990", "K1990modified", "BC1992", "4comp")
  if( !model.correct ){
   stop("'model' is not supported")
  }
  
  steptol = .Machine$double.eps
  
  # get the data and other inputs to mlmaximize -----------------------------
  # cat("0\n")
  
  if(mean.u.0i.zero) mean.u.0i <- NULL
  if(model == "4comp"){
   ln.var.v.it <- ln.var.u.0i <- mean.u.0i <- NULL
  }
  
  YXZ <- .prepareYXZ.panel.simple(formula = formula, data, it, subset, ln.var.v.it = ln.var.v.it, ln.var.u.0i = ln.var.u.0i, mean.u.0i = mean.u.0i, sysnframe = sys.nframe())
  
  if(only.data) return(YXZ)
  
  y <- YXZ$Y
  X <- YXZ$X
  n <- YXZ$n
  nt <- YXZ$nt
  
  timevar <- YXZ$timevar
  # maxT <- rep(t0, times = t0) #YXZ$dat.descr[5]
  
  idvar <- YXZ$idvar
  ids <- YXZ$ids
  t0 <- YXZ$t0
  dat.descr <- YXZ$dat.descr
  
  my.t <- timevar - min(timevar) + 1
  
  # print(summary(my.t))
  
  maxT <- as.vector( unlist(by(data = my.t, INDICES = idvar, FUN = function(qq) rep(max(qq), length(qq))) ))
  
  # print(summary(maxT))
  
  # print(maxT)
  
  esample <- YXZ$esample
  
  Zvi <- YXZ$Zvi
  Zu0 <- YXZ$Zu0
  Zdeli <- YXZ$Zdeli
  # print(head(Zdeli))
  
  Kb <- YXZ$Kb
  Kvi <- YXZ$Kvi
  Ku0 <- YXZ$Ku0
  Kdeli <- YXZ$Kdeli
  
  abbr.length <- 121
  names_x  <- abbreviate(colnames(X), abbr.length, strict = TRUE, dot = FALSE, method = "both.sides")
  names_vi <- abbreviate(colnames(Zvi), abbr.length, strict = TRUE, dot = FALSE, method = "both.sides")
  names_u0 <- abbreviate(colnames(Zu0), abbr.length, strict = TRUE, dot = FALSE, method = "both.sides")
  names_udeli <- abbreviate(colnames(Zdeli), abbr.length, strict = TRUE, dot = FALSE, method = "both.sides")
  max.name.length <- max(nchar(c(names_x, names_vi, names_u0, names_udeli)))
  
  # if(mean.u.0i.zero){
  #  coef.names <- c(
  #   names_x,
  #   paste("ln_vit_var_",c("Intercept", names_vi[-1]),"", sep = ""),
  #   paste("ln_u0_var_",c("Intercept", names_u0[-1]),"", sep = "")
  #  )
  # } else {
  #  coef.names <- c(
  #   names_x,
  #   paste("ln_vit_var_",c("Intercept", names_vi[-1]),"", sep = ""),
  #   paste("ln_u0_var_",c("Intercept", names_u0[-1]),"", sep = ""),
  #   paste("u0_mean_",c("Intercept", names_udeli[-1]),"", sep = "")
  #  )
  # }
  
  coef.names.full.no.eta <- c(
   names_x,
   paste("lnVARvit_",c("(Intercept)", names_vi[-1]),"", sep = ""),
   paste("lnVARu0i_",c("(Intercept)", names_u0[-1]),"", sep = ""),
   paste("mean_u0i_",c("(Intercept)", names_udeli[-1]),"", sep = "")
  )
  
  coef.names.full.no.eta.output <- c(
   names_x,names_vi,names_u0,names_udeli
  )
  
  if(model == "4comp"){
   coef.names.full <- coef.names.full.no.eta <- c(
    names_x, "ln_lmd", "ln_sig", "lnVARv0","lnVARu0")
  }
  
  coef.fixed <- rep( FALSE, length(coef.names.full.no.eta) )
  
  # print(coef.names)
  
  # cat("1\n")

  # model: 1st and 2nd gen --------------------------------------------------
  if( model %in% c("K1990", "K1990modified", "BC1992")){
   
  # starting values 1st and 2nd gen -----------------------------------------
  
   # run OLS anyways
   ols <- lm(y ~ X - 1)
   ols.r1 <- resid(ols)
   ols.r1_mean <- as.vector(by(ols.r1, idvar, mean))
   ols.coef <- c(coef(ols))
   names(ols.coef)[1] <- "(Intercept)"
   # print(model)
   # cat("eff.time.invariant = ", eff.time.invariant,"\n")
   theta0 <- c(1.0*ols.coef, .75*coef(lm(log(ols.r1^2) ~ Zvi - 1)), .75*coef(lm(log(ols.r1_mean^2) ~ Zu0 - 1)), 0, rep(0, (ncol(Zdeli)-1)))
   names(theta0) <- c(coef.names.full.no.eta)
   
   
   if(model == "K1990"){
    if(eff.time.invariant) {
     coef.fixed <- c(coef.fixed, rep(TRUE, 2))
     model.no <- "0"
     coef.names.output <- c(coef.names.full.no.eta.output)
     time.fn <- "Constant over time"
    } else {
     coef.fixed <- c(coef.fixed, rep(FALSE, 2))
     theta0 <- c(theta0, rep(0, 2))
     model.no <- "`"
     names(theta0) <- c(coef.names.full.no.eta, paste0(c("b","c"),model.no))
     coef.names.output <- c(coef.names.full.no.eta.output,paste0(c("b","c"),model.no) )
     time.fn <- "( 1 + exp(b*t + c*t^2) )^-1"
    }
   }
   
   if(model == "K1990modified"){
    if(eff.time.invariant) {
     coef.fixed <- c(coef.fixed, rep(TRUE, 2))
     model.no <- "0"
     coef.names.output <- c(coef.names.full.no.eta.output)
     time.fn <- "Constant over time"
    } else {
     coef.fixed <- c(coef.fixed, rep(FALSE, 2))
     theta0 <- c(theta0, rep(0, 2))
     model.no <- "``"
     names(theta0) <- c(coef.names.full.no.eta, paste0(c("d","e"),model.no))
     coef.names.output <- c(coef.names.full.no.eta.output,paste0(c("d","e"),model.no) )
     time.fn <- "1 + d*(t-T_i) + e*(t-T_i)^2"
    }
   }
   
   if(model == "BC1992"){
    if(eff.time.invariant) {
     coef.fixed <- c(coef.fixed, rep(TRUE, 1))
     model.no <- "0"
     coef.names.output <- c(coef.names.full.no.eta.output)
     time.fn <- "Constant over time"
    } else {
     coef.fixed <- c(coef.fixed, rep(FALSE, 1))
     theta0 <- c(theta0, rep(0, 1))
     model.no <- "```"
     names(theta0) <- c(coef.names.full.no.eta, paste0("f",model.no))
     coef.names.output <- c(coef.names.full.no.eta.output,paste0(c("f"),model.no) )
     time.fn <- "exp( -f*(t-T_i) )"
    }
   }
   
   if(mean.u.0i.zero) coef.fixed[Kb+Kvi+Ku0+Kdeli] <- TRUE
   
   theta0 <- theta0[!coef.fixed]
   coef.names.output <- coef.names.output[!coef.fixed]
   
   # cat("3\n")
   # print(theta0)
   
   Ktheta <- length(theta0)
   
   # cat("4\n")
   # print(theta0)
   # cat("2\n")
   # print(theta0)
   # print(rownames(theta0))
   
   if(!is.null(start.val)){
    # if start.val partially available
    if(length(theta0) != length(start.val)){
     warning("\nLength of 'start.val' is not appropriate: trying to deal with it", call. = FALSE,immediate. = TRUE)
     theta2 <- rep(0, length(theta0))
     # print(theta2)
     # print(as.vector(rownames(theta0)))
     # print(as.vector(rownames(start.val)))
     # print(as.vector(rownames(theta0)) %in% as.vector(rownames(start.val)))
     # theta2[as.vector(rownames(theta0)) %in% as.vector(rownames(start.val))] <- start.val
     # print(names(theta0))
     # print(names(start.val))
     # cat("2.5\n")
     # print(theta0)
     # print(start.val)
     # cat("2.6\n")
     theta2[names(theta0) %in% names(start.val)] <- start.val
     # print(theta2)
     names(theta2) <- names(theta0)
     theta0 <- theta2
     # if(model == "K1990" | model == "K1990modified"){
     #  if(theta0["eta1*"] == 0) theta0["eta1*"] <- -1e-3
     #  if(theta0["eta2*"] == 0) theta0["eta2*"] <- 1e-3
     # }
     # if(model == "BC1992"){
     #  if(theta0["eta1*"] == 0) theta0["eta1*"] <- 1e-3
     # }
    } else {
     theta0 <- start.val
    }
   }
   
   # cat("3\n")
   # print(theta0)
   
   # obj10 <- eval(parse(text = paste(" optim(par = theta0, fn = .ll.panel, gr = .gr.panel, prod = prod, my.n = n, k = Kb, kv = Kvi, ku = Ku0, kdel = Kdeli, yit = y, zvit = Zvi, zui = Zu0, xit = X, zdeli = Zdeli, eff.time.invariant = eff.time.invariant, model = model, timevar = my.t, maxT = maxT, ids = ids, idvar = idvar, t0 = t0, method = c('Nelder-Mead'), control = list(reltol  = 1e-16, maxit = maxit^2, trace = ifelse(print.level >=2, 1, -1), REPORT = ifelse(print.level >=2, 7, -1), fnscale = -1), hessian = FALSE)", sep = "")))
   # print(obj10$par)
   
   
   # time.05 <- proc.time()
   # obj20 <- eval(parse(text = paste(" optim(par = theta0, fn = .ll.panel, gr = .gr.panel, coef.fixed = coef.fixed, prod = prod, my.n = n, Ktheta = Ktheta, k = Kb, kv = Kvi, ku = Ku0, kdel = Kdeli, yit = y, zvit = Zvi, zui = Zu0, xit = X, zdeli = Zdeli, eff.time.invariant = eff.time.invariant, mean.u.0i.zero = mean.u.0i.zero, model = model, timevar = my.t, maxT = maxT, ids = ids, idvar = idvar, t0 = t0, method = c('BFGS'), control = list(reltol  = 1e-16, maxit = maxit, trace = ifelse(print.level >=2, 1, -1), REPORT = ifelse(print.level >=2, report.ll.optim, -1), fnscale = -1), hessian = TRUE)", sep = "")))
   # co_ <- obj20$par
   # sd_ <- suppressWarnings( sqrt(diag(solve(-obj20$hessian))) )
   # ts_ <- co_ / sd_
   # pv_ <- pt(abs(ts_), n-length(theta0), lower.tail = FALSE) * 2
   # printCoefmat(cbind(co_,sd_,ts_,pv_), digits = digits)           
   # # print(cbind(obj20$par, 1))
   # time.06 <- proc.time()
   # est.time.sec <- (time.06-time.05)[1]
   # names(est.time.sec) <- "sec"
   # if(print.level >= 2){
   #  .timing(est.time.sec, "\n Log likelihood maximization completed in ")
   #  # cat("___________________________________________________\n")
   # }
   
   
   # obj30 <- eval(parse(text = paste(" optim(par = obj20$par, fn = .ll.panel, gr = .gr.panel, prod = prod, my.n = n, k = Kb, kv = Kvi, ku = Ku0, kdel = Kdeli, yit = y, zvit = Zvi, zui = Zu0, xit = X, zdeli = Zdeli, eff.time.invariant = eff.time.invariant, model = model, timevar = timevar, maxT = maxT, ids = ids, idvar = idvar, t0 = t0, method = c('Nelder-Mead'), control = list(reltol  = 1e-16, maxit = maxit, trace = ifelse(print.level >=2, 1, -1), REPORT = ifelse(print.level >=2, 7, -1), fnscale = -1), hessian = FALSE)", sep = "")))
   # 
   # obj40 <- eval(parse(text = paste(" optim(par = obj30$par, fn = .ll.panel, gr = .gr.panel, prod = prod, my.n = n, k = Kb, kv = Kvi, ku = Ku0, kdel = Kdeli, yit = y, zvit = Zvi, zui = Zu0, xit = X, zdeli = Zdeli, eff.time.invariant = eff.time.invariant, model = model, timevar = timevar, maxT = maxT, ids = ids, idvar = idvar, t0 = t0, method = c('BFGS'), control = list(reltol  = 1e-32, maxit = maxit, trace = ifelse(print.level >=2, 1, -1), REPORT = ifelse(print.level >=2, 7, -1), fnscale = -1), hessian = FALSE)", sep = "")))
   # 
   # obj50 <- eval(parse(text = paste(" optim(par = obj40$par, fn = .ll.panel, gr = .gr.panel, prod = prod, my.n = n, k = Kb, kv = Kvi, ku = Ku0, kdel = Kdeli, yit = y, zvit = Zvi, zui = Zu0, xit = X, zdeli = Zdeli, eff.time.invariant = eff.time.invariant, model = model, timevar = timevar, maxT = maxT, ids = ids, idvar = idvar, t0 = t0, method = c('Nelder-Mead'), control = list(reltol  = 1e-16, maxit = maxit, trace = ifelse(print.level >=2, 1, -1), REPORT = ifelse(print.level >=2, 7, -1), fnscale = -1), hessian = FALSE)", sep = "")))
   # 
   # print(cbind(obj10$par,obj20$par,obj30$par,obj40$par,obj50$par))
   # print(c(obj10$value,obj20$value,obj30$value,obj40$value,obj50$value))
   
   # obj  <- obj10

  # optimization 1st and 2nd gen -------------------------------------------------------------
  
  if(print.level >= 2){
   my.message <- paste0(" Optimization using 'mlmaximize': ")
   m.left <- floor( (max.name.length+42-nchar(my.message)-1) / 2 )
   m.right <- max.name.length+42-nchar(my.message) -1- m.left
   if(m.left <= 0) m.left <- 2
   if(m.right <= 0) m.right <- 2
   cat("\n",rep("-", m.left), my.message, rep("-", m.right),"\n\n", sep ="")
  }
  
  time.05 <- proc.time()
  
  # print(theta0)
  # print(length(theta0))
  
  obj <- eval(parse(text = paste(".mlmaximize.panel(theta0, .ll.panel, gr = .gr.panel, hess = .hess.panel, gr.hess = .gr.hess.panel, coef.fixed = coef.fixed, prod = prod, my.n = n, Ktheta = Ktheta, k = Kb, kv = Kvi, ku = Ku0, kdel = Kdeli, yit = y, zvit = Zvi, zui = Zu0, xit = X, zdeli = Zdeli, eff.time.invariant = eff.time.invariant, mean.u.0i.zero = mean.u.0i.zero, model = model, timevar = my.t, maxT = maxT, ids = ids, idvar = idvar, print.level = print.level, t0 = t0, reltol  = reltol, lmtol = lmtol, steptol = steptol, max.backedup = 17, maxit = maxit, when.backedup = .Machine$double.eps)", sep = "")))
  
  
  
  time.06 <- proc.time()
  est.time.sec <- (time.06-time.05)[1]
  names(est.time.sec) <- "sec"
  if(print.level >= 2){
   .timing(est.time.sec, "\n Log likelihood maximization completed in ")
   # cat("___________________________________________________\n")
  }
  # print(obj)
  if(obj$converged == 0) return(obj)
  
  # print(obj$table)

  # 1st and 2nd: auxiliary parameters ---------------------------------------
  
  if(is.null(ln.var.v.it) | is.null(ln.var.u.0i)){
   # auxiliary parameters
   # SE for sigmau and sigmav
   sigv <- sqrt(exp(obj$par[Kb+1]))
   sigu <- sqrt(exp(obj$par[Kb+Kvi+1]))
   se_sigv <- (0.5*sigv) * sqrt(obj$vcov[Kb+1,Kb+1])
   se_sigu <- (0.5*sigu) * sqrt(obj$vcov[Kb+Kvi+1,Kb+Kvi+1])
   
   # cat("sigv = ",sigv,"se_sigv = ",se_sigv,"\n")
   # cat("sigu = ",sigu,"se_sigu = ",se_sigu,"\n")
   
   # SE for lambda
   lmd <- sqrt(exp(obj$par[Kb+Kvi+1] - obj$par[Kb+1]))
   glmd <- c(0.5*lmd, -0.5*lmd)
   se_lmd <- as.vector( sqrt( t(glmd) %*% obj$vcov[c((Kb+1),(Kb+Kvi+1)), c((Kb+1),(Kb+Kvi+1))] %*% glmd))
   
   # cat("lmd = ",lmd,"se_lmd = ",se_lmd,"\n")
   
   # SE for gamma
   gam <- exp(obj$par[Kb+Kvi+1])/(exp(obj$par[Kb+Kvi+1]) + exp(obj$par[Kb+1]))
   ggam <- c(gam*(1 - gam), -gam^2/lmd^2)
   se_gam <- as.vector( sqrt( t(ggam) %*% obj$vcov[c((Kb+1),(Kb+Kvi+1)), c((Kb+1),(Kb+Kvi+1))] %*% ggam))
   
   # cat("gam = ",gam,"se_gam = ",se_gam,"\n")
  }
  
  par.aux <- output <- NULL
  if(is.null(ln.var.v.it) & is.null(ln.var.u.0i)){
   par.aux <- c(sigv, sigu, lmd, gam); se.aux <- c(se_sigv, se_sigu, se_lmd, se_gam)
   names(par.aux) = c("sigma_vit  ", "sigma_u0i  ", "lambda   ", "gamma   ")
  } else if(is.null(ln.var.u.0i) & !is.null(ln.var.v.it)){
   par.aux <- c(sigu); se.aux <- c(se_sigu)
   names(par.aux) = c("sigma_u0i  ")
  } else if(!is.null(ln.var.v.it) & is.null(ln.var.u.0i)){
   par.aux <- c(sigv); se.aux <- c(se_sigv)
   names(par.aux) = c("sigma_vit  ")
  } else {
   my.message = "het"
  }
  
  # output estimation results -----------------------------------------------
  
  my.table <- obj$table[,1:4]
  # rownames(my.table)[1:(Kb+Kvi+Ku0+Kdeli)] <- coef.names.output
  rownames(my.table) <- coef.names.output
  if(!is.null(par.aux)){
   output <- rbind(my.table, 
                   cbind(par.aux, se.aux, par.aux/se.aux, pnorm(abs(par.aux/se.aux), lower.tail = FALSE)*2))
  } else {
   output <- my.table
  }
  
  myeff <- ifelse(prod, "technical", "cost")
  
  if(print.level >= 1){
   my.gen <- "2nd"
   if(eff.time.invariant) my.gen <- "1st"
   cat("",rep("_", max.name.length+42-1),"", "\n", sep = "")
   cat("\nStochastic ",myeff," frontier panel data (",my.gen," gen.) model\n",  sep = "")
   cat("\nDistributional assumptions\n\n", sep = "")
   Assumptions <- rep("heteroskedastic",2)
   if(is.null(ln.var.v.it)) Assumptions[1] <- "homoskedastic"
   if(is.null(ln.var.u.0i)) Assumptions[2] <- "homoskedastic"
   Distribution = c("normal ", "half-normal ")
   if(!is.null(mean.u.0i)){
    Distribution[2] <- "truncated-normal "
   }
   a1 <- data.frame(
    Component = c("Random noise:","Inefficiency:"),
    Distribution = Distribution,
    Assumption = Assumptions
   )
   print(a1, quote = FALSE, right = FALSE)
  }
  
  if(print.level >= 1){
   est.spd.left <- floor( (max.name.length+42-29) / 2 )
   est.spd.right <- max.name.length+42-29 - est.spd.left
   # cat("\n------------")
   # cat(" Summary of the panel data: ------------\n\n", sep = "")
   # cat("\n------------ Summary of the panel data: ----------\n\n", sep = "")
   cat("\n",rep("-", est.spd.left)," Summary of the panel data: ",rep("-", est.spd.right),"\n\n", sep ="")
   cat("   Number of obs       (NT) =",nt,"", "\n")
   cat("   Number of groups     (N) =",n,"", "\n")
   cat("   Obs per group: (T_i) min =",dat.descr[3],"", "\n")
   cat("                        avg =",dat.descr[4],"", "\n")
   cat("                        max =",dat.descr[5],"", "\n")
  }
  
  if(print.level >= 1){
   # cat("\n---------------")
   # max.name.length <- max(nchar(row.names(output)))
   est.rez.left <- floor( (max.name.length+42-22) / 2 )
   est.rez.right <- max.name.length+42-22 - est.rez.left
   cat("\n",rep("-", est.rez.left)," Estimation results: ",rep("-", est.rez.right),"\n\n", sep ="")
   # cat("\n--------------- Estimation results: --------------\n\n", sep = "")
   # .printgtresfhet(output, digits = digits, Kb, Kv0, Ku0, Kvi, Kui, na.print = "NA", max.name.length)
   .printpanel2nd(x = output, digits = 4, kb = Kb, kvi  = Kvi, ku0 = Ku0, kdeli= Kdeli, eff.time.invariant = eff.time.invariant, mean.u.0i.zero = mean.u.0i.zero, model = model, na.print = "NA", max.name.length = max.name.length, mycutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), mysymbols = c("***", "**", "*", ".", " "))
   
   # cat.print(output)
  }
  
  if(!eff.time.invariant){
   if(model == "BC1992"){
    my.etas <- "f"
   } else if(model == "K1990modified"){
    my.etas <- "d and e"
   } else if(model == "K1990"){
    my.etas <- "b and c"
   }
   
   cat("\n",model.no,"Pay attention to interpretation of",my.etas,"\n")
   cat(" Function of inefficiency change over time is\n ",time.fn," \n")
  }
  
  
  # Efficiency calculation --------------------------------------------------
  
  eit <- as.vector( y - X%*%obj$par[1:Kb]) 
  # if(model == "K1990") eit <- eit + obj$par[1]
  sigv2it <- as.vector(  exp(Zvi%*%obj$par[(Kb+1):(Kb+Kvi)]))
  gu <- obj$par[(Kb+Kvi+1):(Kb+Kvi+Ku0)]
  sigu2i <- as.vector(  exp(Zu0%*%gu) )
  delta <- obj$par[(Kb+Kvi+Ku0+1):(Kb+Kvi+Ku0+Kdeli)]
  if(mean.u.0i.zero){
   mui <- rep(0, n)
  } else {
   mui <- as.vector( Zdeli%*%delta )
   # Print truncated point
   if(print.level >= 2){
    
    sum.te <- paste0(" Summary of mean(s) of the trancated distribution: ")
    est.cle.left <- floor( (max.name.length+42-nchar(sum.te)-1) / 2 )
    est.cle.right <- max.name.length+42-nchar(sum.te) -1- est.cle.left
    if(est.cle.left <= 0) est.cle.left <- 1
    if(est.cle.right <= 0) est.cle.right <- 1
    cat("\n",rep("-", est.cle.left),sum.te,rep("-", est.cle.right),"\n", sep ="")
    
    # eff1 <- data.frame(eff[,4:6])
    # colnames(eff1) <- c("exp(-E[ui|ei])", "exp(-M[ui|ei])", "E[exp(-ui)|ei]")
    cat("(id-specific time-invariant)", "\n\n", sep = "")
    .su(mui, print = TRUE, width = 5, format = "fg", drop0trailing = FALSE, names = c("mu_i"))
   }
  }
  
  
  if(eff.time.invariant){
   Git <- rep(1, length(my.t))
  } else {
   # eta <- obj$par[-c(1:(Kb+Kvi+Ku0+Kdeli))]
   if(model == "BC1992"){
    eta <- obj$par[c( Ktheta )]
    Git = exp(-eta*(my.t - maxT))
   } else if(model == "K1990modified"){
    eta <- obj$par[c( (Ktheta-1) : (Ktheta) )]
    Git = 1 + eta[1]*(my.t - maxT) + eta[2]*(my.t - maxT)^2
   } else if(model == "K1990"){
    eta <- obj$par[c( (Ktheta-1) : (Ktheta) )]
    Git = (1 + exp(eta[1]*my.t + eta[2]*my.t^2))^(-1)
   } else {
    stop("Unknown model\n")
   }
  }
  
  # prod1 <- prod
  # if(!prod & cost.eff.less.one) prod1 <- TRUE 
  eff <- .u2efftnm.panel(eit, sigv2it, sigu2i, mui, Git, alpha = alpha, prod = prod, cost.eff.less.one, ids, idvar, timevar, t0, it.names = it)
  
  if(print.level >= 2){
   
   sum.te <- paste0(" Summary of ",myeff," efficiencies: ")
   est.cle.left <- floor( (max.name.length+42-nchar(sum.te)-1) / 2 )
   est.cle.right <- max.name.length+42-nchar(sum.te) -1- est.cle.left
   if(est.cle.left <= 0) est.cle.left <- 1
   if(est.cle.right <= 0) est.cle.right <- 1
   cat("\n",rep("-", est.cle.left),sum.te,rep("-", est.cle.right),"\n\n", sep ="")
   
   # cat("",rep(" ", est.eff.left+1),"JLMS:= exp(-E[ui|ei])\n", sep = "")
   # cat("",rep(" ", est.eff.left+1),"Mode:= exp(-M[ui|ei])\n", sep = "")
   # cat("",rep(" ", est.eff.left+3),"BC:= E[exp(-ui)|ei]\n", sep = "")
   
   eff1 <- data.frame(eff[,4:6])
   colnames(eff1) <- c("exp(-E[ui|ei])", "exp(-M[ui|ei])", "E[exp(-ui)|ei]")
   
   .su(eff1, print = TRUE, width = 5, format = "fg", drop0trailing = FALSE, names = c("exp(-E[ui|ei])", "exp(-M[ui|ei])", "E[exp(-ui)|ei]"))
  }
  
  # print(summary(eff))
  
  # Marginal effects
  if(marg.eff){
   meff = tryCatch(.me(as.matrix(obj$par)[-c(1:(Kb+Kvi)),1,drop = FALSE], Zu = Zu0, Zdel = Zdeli, ku = Ku0, kdel = Kdeli, n = n, dist = "t"), error = function(e) { cat('In error handler\n'); print(e); e })
   if(inherits(meff, "error"))  {
    meff <- NULL
   } else if (!is.null(meff)){
    meff <- as.data.frame(meff)
    colnames(meff) <- colnames(Zu0)[-1]
   }
  } else {
   meff <- NULL
  }
  
  if (!is.null(meff)) .su(meff)
  
  # colnames(obj$vcov)[1:(Kb+Kvi+Ku0+Kdeli)] <- rownames(obj$vcov)[1:(Kb+Kvi+Ku0+Kdeli)] <- names(theta0)#c(names_x, names_zv, names_zu, names_del)
  # colnames(obj$vcov) <- rownames(obj$vcov) <- names(theta0)#c(names_x, names_zv, names_zu, names_del)
  
  mymodel <- paste0("sf_p_2_", model)
  if(eff.time.invariant)  mymodel <- "sf_p_1"
  
  temp <- list(call = match.call(), model = mymodel, coef = obj$par, table = output, assumptions = a1, vcov = obj$vcov, ll = obj$ll, lmtol = lmtol, LM = obj$gHg, convergence = obj$converged, esttime = est.time.sec, prod = prod, efficiencies = eff, marg.effects = meff, Kb = Kb, Kvi = Kvi, Ku0 = Ku0, Kdeli = Kdeli, eff.time.invariant = eff.time.invariant,  mean.u.0i.zero = mean.u.0i.zero, ln.var.v.it = ln.var.v.it, ln.var.u.0i = ln.var.u.0i, model.no = model.no, time.fn = time.fn, fitted = as.vector(unname(X%*%obj$par[1:Kb])), yobs = as.vector(y), xobs = X, n = n, nt = nt, dat.descr = dat.descr, esample = esample)

  # print(temp)
  if (obj$converged == 2){
   temp$delta_rel <- obj$delta_rel
   temp$ltol <- obj$ltol
  }
  if (obj$converged == 3){
   temp$theta_rel_ch <- obj$theta_rel_ch
   temp$steptol <- sqrt(.Machine$double.eps)
  }
  
  } # end of if( model %in% c("K1990", "K1990modified", "BC1992")){

 # model: 4comp ------------------------------------------------------------
  
  if(model == "4comp"){
   
   # drawing random deviates -------------------------------------------------
   
   # myprimes <- c(2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 
                 # 59, 61, 67, 71, 73, 79, 83, 89, 97, 101)
   if(length(simtype) != 1) simtype <- "halton"
   if(length(simtype_GHK) != 1) simtype_GHK <- "halton"
   if(simtype == "rnorm"){
    V0 <- matrix(rnorm(R*n), nrow = n, ncol = R)
    U0 <- abs( matrix(rnorm(R*n), nrow = n, ncol = R) )
   } else {
     set.seed(seed)
     if(is.null(halton.base)){
       # deviates for each ID come from own base, but then scaled
       myrand0 <- t( npsf::halton(R, n.bases = 2*n, random.primes = random.primes, seed = seed, start = 1, scale.coverage = TRUE, shuffle = TRUE) )
     } else {
       if(halton.base > 100008){
         stop("halton.base should be smaller than 100,008")
       } else {
         # deviates for all ID come from a sequence with said base
         myrand0 <- t( matrix( npsf::halton(R*2*n, bases = halton.base, random.primes = random.primes, seed = seed, start = 1, scale.coverage = TRUE, shuffle = TRUE), nrow = R, ncol = 2*n) )
       }
     }
     # apply qnorm
     V0 <-      qnorm(myrand0[1:n,])
     U0 <- abs( qnorm(myrand0[(n+1):(2*n),]) )
     rm(myrand0)
   }
     
    # if(initialize.halton){
    #   
    #   
    #   
    #  # begin initiated sequence
    #  # total number generated dimentions
    #  tngdim <- floor(n*2)
    #  # initialize Halton algorithm
    #  # myrand0 <- randtoolbox::halton(R, dim = tngdim, init = TRUE, normal = TRUE)
    #  myrand0 <- t( randtoolbox::halton(R, dim = tngdim, init = TRUE, normal = TRUE) )
    #  # myrand1 <- t( randtoolbox::halton(R, dim = tngdim, init = FALSE, normal = TRUE) )
    #  # set.seed(seed)
    #  # smpl.v <- sample(seq_len(tngdim), n)
    #  # smpl.u0 <- sample(seq_len(tngdim), n)
    #  # smpl.u <- smpl.u0
    #  # my.repeats <- 0
    #  # while (any(smpl.v - smpl.u) == 0) {
    #  #   smpl.u <- sample(smpl.u0)
    #  #   my.repeats <- my.repeats + 1
    #  #   print(my.repeats)
    #  # }
    #  # rm(my.repeats)
    #  # V0 <-      myrand0[smpl.v,]
    #  # U0 <- abs( myrand1[smpl.u,] )
    #  # V0 <-      myrand0[(1:n),]
    #  # U0 <- abs( myrand0[((n+1):(2*n)),] )
    #  # rm(myrand0)
    #  # another way to do this
    #  myrand0 <- t( randtoolbox::halton(R, dim = 2*n, init = TRUE, normal = TRUE) )
    #  V0 <-      myrand0[1:n,]
    #  U0 <- abs( myrand0[(n+1):(2*n),] )
    #  rm(myrand0)
    #  # end initiated sequence
    # } else {
    #  # begin noninitiated sequence
    #  if(halton.base < 2 ){
    #   stop("base for Halton draws should be >= 2")
    #  } else {
    #   if(halton.base > 101){
    #    stop("pick a prime smaller or equal 101")
    #   }
    #   is.halton.base.prime <- sum(halton.base == myprimes) == 1
    #   if(!is.halton.base.prime){
    #    mydiff <- ifelse(halton.base-myprimes < 0, 102, halton.base-myprimes)
    #    closest.prime <- which(mydiff == min(mydiff))
    #    warning(paste0("halton.base must be a prime: picked ",myprimes[closest.prime]," as base"))
    #   } else {
    #    closest.prime <- which(myprimes == halton.base)
    #   }
    #   if(is.null(my.leap)) my.leap <- 1
    #   V0 <- qnorm(matrix(sfsmisc::sHalton(R*n*my.leap, leap = my.leap, base = myprimes[closest.prime]), nrow = n, ncol = R, byrow = TRUE))
    #   U0 <- abs( qnorm(matrix(sfsmisc::sHalton(R*n*my.leap, leap = my.leap, base = myprimes[closest.prime+7]), nrow = n, ncol = R, byrow = TRUE)) )
    #  }
    #  # end noninitiated sequence
    # }
   # }
   
   # starting values 4comp ---------------------------------------------------
   
   # run OLS anyways
   ols <- lm(y ~ X - 1)
   ols.r1 <- resid(ols)
   ols.r1_mean <- as.vector(by(ols.r1, idvar, mean))
   ols.r1_mean.sd <- sd(ols.r1_mean)
   ols.coef <- c(coef(ols))
   names(ols.coef)[1] <- "(Intercept)"
   # print(model)
   # cat("eff.time.invariant = ", eff.time.invariant,"\n")
   theta0 <- c(1.0*ols.coef, ln_lmd = log(1), ln_sig = log(sqrt(2*ols.r1_mean.sd^2)), ln_sigv0 = log(ols.r1_mean.sd), ln_sigu0 = log(ols.r1_mean.sd))
   names(theta0) <- c(coef.names.full)
   
   Ktheta <- length(theta0)
   k <- Ktheta - 4
   
   if(!is.null(start.val)){
    # if start.val partially available
    if(length(theta0) != length(start.val)){
     warning("\nLength of 'start.val' is not appropriate: trying to deal with it", call. = FALSE,immediate. = TRUE)
     theta2 <- rep(0, length(theta0))
     # print(theta2)
     # print(as.vector(rownames(theta0)))
     # print(as.vector(rownames(start.val)))
     # print(as.vector(rownames(theta0)) %in% as.vector(rownames(start.val)))
     # theta2[as.vector(rownames(theta0)) %in% as.vector(rownames(start.val))] <- start.val
     # print(names(theta0))
     # print(names(start.val))
     # cat("2.5\n")
     # print(theta0)
     # print(start.val)
     # cat("2.6\n")
     theta2[names(theta0) %in% names(start.val)] <- start.val
     # print(theta2)
     names(theta2) <- names(theta0)
     theta0 <- theta2
     # if(model == "K1990" | model == "K1990modified"){
     #  if(theta0["eta1*"] == 0) theta0["eta1*"] <- -1e-3
     #  if(theta0["eta2*"] == 0) theta0["eta2*"] <- 1e-3
     # }
     # if(model == "BC1992"){
     #  if(theta0["eta1*"] == 0) theta0["eta1*"] <- 1e-3
     # }
    } else {
     theta0 <- start.val
    }
    
    # optimize 4comp ----------------------------------------------------------

    # optimize 4comp; start val provided --------------------------------------
    
    time.05 <- proc.time()

    # obj <- tryCatch( mlmaximize(theta0 = theta0, ll = C_gtrelnls, gr = C_gtregrad, hess = C_gtrehess, BHHH = FALSE, alternate = FALSE, prod = prod, log.lmd = log.lmd, log.sig = log.sig, log.sigv0 = log.sigv0, log.sigu0 = log.sigu0, level = 0.99, step.back = step.back, reltol = reltol, lmtol = lmtol, steptol = steptol, digits = digits, when.backedup = when.backedup, max.backedup = max.backedup, only.maximize = only.maximize, print.level = print.level, n = nt, maxit = maxit), error = function(e) e )
    obj <- eval(parse(text = paste("tryCatch(.mlmaximize.panel(theta0 = theta0, ll = .fourC_gtrelnls, gr = .fourC_gtregrad, hess = .fourC_gtrehess, gr.hess = .fourC_gtregr_gtrehess, alternate = NULL, BHHH = FALSE, prod = prod, V0=V0,U0=U0,yit=y,xit=X,ids=ids,idvar=idvar,nt=nt,my.n=n,n=n,R=R,idlenmax=dat.descr[5],Ktheta=Ktheta,n_threads=1, step.back = .Machine$double.eps^.5, reltol =  sqrt(.Machine$double.eps), lmtol =  lmtol, steptol =  .Machine$double.eps, digits = 4, when.backedup = sqrt(.Machine$double.eps), max.backedup = 17, print.level = print.level, only.maximize = FALSE, maxit = maxit), error = function(e) e )", sep = "")))
    cannot.est.model1 <- inherits(obj, "error")
    
   } else {
    # starting values are not provided
    
    time.05 <- proc.time()
    
    # optimize 4comp; start val not provided ----------------------------------

    if(print.level >= 2){
     my.message <- paste0(" Optimization using 'mlmaximize': ")
     m.left <- floor( (max.name.length+42-nchar(my.message)-1) / 2 )
     m.right <- max.name.length+42-nchar(my.message) -1- m.left
     if(m.left <= 0) m.left <- 2
     if(m.right <= 0) m.right <- 2
     cat("\n",rep("-", m.left), my.message, rep("-", m.right),"\n\n", sep ="")
    }
    
    
    
    # print(theta0)
    # print(length(theta0))
    
    
    # compare times in c and numDeriv: begin ----------------------------------
    
    
    # t123 <- .fourC_gtrelnls(theta=theta0, prod = prod, V0=V0,U0=U0,yit=y,xit=X,ids=ids,idvar=idvar,nt=nt,my.n=n,R=R,idlenmax=dat.descr[5],Ktheta=Ktheta,n_threads=1)
    # cat.print(t123)
    
    # cat("\n\n\nGradient\n\n\n")
    # 
    # q10 <- proc.time()
    # t1234 <- .fourC_gtregrad(theta=theta0, prod = prod, V0=V0,U0=U0,yit=y,xit=X,ids=ids,idvar=idvar,nt=nt,my.n=n,R=R,idlenmax=dat.descr[5],Ktheta=Ktheta,n_threads=1)
    # q11 <- proc.time()
    # cat.print(q11 - q10)
    # cat.print(t1234)
    # q10 <- proc.time()
    # t1234_num <- numDeriv::grad(func=.fourC_gtrelnls, x=theta0, prod = prod, V0=V0,U0=U0,yit=y,xit=X,ids=ids,idvar=idvar,nt=nt,my.n=n,R=R,idlenmax=dat.descr[5],Ktheta=Ktheta,n_threads=1)
    # q11 <- proc.time()
    # cat.print(q11 - q10)
    # cat.print(t1234_num)
    # cat.print(t1234 - t1234_num)
    # 
    # cat("\n\n\nHessian\n\n\n")
    # 
    # q10 <- proc.time()
    # t1234 <- .fourC_gtregr_gtrehess(theta=theta0, prod = prod, V0=V0,U0=U0,yit=y,xit=X,ids=ids,idvar=idvar,nt=nt,my.n=n,R=R,idlenmax=dat.descr[5],Ktheta=Ktheta,n_threads=1)$hessian1
    # q11 <- proc.time()
    # cat.print(q11 - q10)
    # cat.print(t1234)
    # q10 <- proc.time()
    # t1234_num <- hessian1(func=.fourC_gtrelnls, at=theta0, prod = prod, V0=V0,U0=U0,yit=y,xit=X,ids=ids,idvar=idvar,nt=nt,my.n=n,R=R,idlenmax=dat.descr[5],Ktheta=Ktheta,n_threads=1)
    # q11 <- proc.time()
    # cat.print(q11 - q10)
    # cat.print(t1234_num)
    # cat.print(t1234 - t1234_num)
    
    # stop("this")
    
    # compare times in c and numDeriv: end ------------------------------------
    
    # obj <- eval(parse(text = paste(".mlmaximize.panel(theta0 = theta0, ll = .fourC_gtrelnls, gr = .fourC_gtregrad, hess = .fourC_gtrehess, gr.hess = .fourC_gtregr_gtrehess, alternate = NULL, BHHH = FALSE, prod = prod, V0=V0,U0=U0,yit=y,xit=X,ids=ids,idvar=idvar,nt=nt,my.n=n,n=n,R=R,idlenmax=dat.descr[5],Ktheta=Ktheta,n_threads=1, step.back = .Machine$double.eps^.5, reltol =  sqrt(.Machine$double.eps) , lmtol =  sqrt(.Machine$double.eps) , steptol =  .Machine$double.eps, digits = 4, when.backedup = sqrt(.Machine$double.eps), max.backedup = 17, print.level = print.level, only.maximize = FALSE, maxit = maxit)", sep = "")))
    # cat("I am here")
    obj <- eval(parse(text = paste("tryCatch( .mlmaximize.panel(theta0 = theta0, ll = .fourC_gtrelnls, gr = .fourC_gtregrad, hess = .fourC_gtrehess, gr.hess = .fourC_gtregr_gtrehess, alternate = NULL, BHHH = FALSE, prod = prod, V0 = V0, U0 = U0, yit = y, xit = X, ids = ids, idvar = idvar, nt = nt, my.n = n, n = n, R = R, idlenmax = dat.descr[5], Ktheta = Ktheta, n_threads = 1, step.back = .Machine$double.eps^.5, reltol = sqrt(.Machine$double.eps), lmtol = lmtol, steptol = .Machine$double.eps, digits = 4, when.backedup = sqrt(.Machine$double.eps), max.backedup = 17, print.level = print.level, only.maximize = FALSE, maxit = maxit), error = function(e) e)", sep = "")))
    
    
   } # end '# starting values are not provided'
   
   time.06 <- proc.time()
   est.time.sec <- (time.06-time.05)[1]
   names(est.time.sec) <- "sec"
   if(print.level >= 2){
    .timing(est.time.sec, "\n Log likelihood maximization completed in ")
    # cat("___________________________________________________\n")
   }
   # print(obj)
   if(obj$converged == 0) return(obj)


   # output estimation results -----------------------------------------------

   Kv0 <- Ku0 <- Kvi <- Kui <- 1
   
   if(print.level >= 1){
    cat("",rep("_", max.name.length+42-1),"", "\n", sep = "")
    cat("\nStochastic ",ifelse(prod, "production", "cost")," frontier 4 components model\n",  sep = "")
    cat("\nDistributional assumptions\n\n", sep = "")
    Assumptions <- rep("heteroskedastic",4)
    if(Kv0==1) Assumptions[1] <- "homoskedastic"
    if(Ku0==1) Assumptions[2] <- "homoskedastic"
    if(Kvi==1) Assumptions[3] <- "homoskedastic"
    if(Kui==1) Assumptions[4] <- "homoskedastic"
    a1 <- data.frame(
     Component = c("Random effects:","Persistent ineff.: ","Random noise:","Transient ineff.: "),
     Distribution = c("normal", "half-normal  ","normal", "half-normal  "),
     Assumption = Assumptions
    )
    toinclude <- rep(TRUE,4)
    print(a1[toinclude,], quote = FALSE, right = FALSE)
   }
   
   if(print.level >= 1){
    est.spd.left <- floor( (max.name.length+42-29) / 2 )
    est.spd.right <- max.name.length+42-29 - est.spd.left
    # cat("\n------------")
    # cat(" Summary of the panel data: ------------\n\n", sep = "")
    # cat("\n------------ Summary of the panel data: ----------\n\n", sep = "")
    cat("\n",rep("-", est.spd.left)," Summary of the panel data: ",rep("-", est.spd.right),"\n\n", sep ="")
    cat("   Number of obs       (NT) =",nt,"", "\n")
    cat("   Number of groups     (N) =",n,"", "\n")
    cat("   Obs per group: (T_i) min =",dat.descr[3],"", "\n")
    cat("                        avg =",dat.descr[4],"", "\n")
    cat("                        max =",dat.descr[5],"", "\n")
   }
   
   par <- c(obj$par)
   names(par) <- coef.names.full
   se <- c(sqrt(diag(obj$vcov)))
   output <- cbind(round(par, digits = digits),
                   round(se, digits = digits),
                   round(par/se,digits = 2),
                   round(pnorm(abs(par/se), lower.tail = FALSE)*2, digits = digits))
   colnames(output) <- c("Coef.", "SE ", "z ",  "P>|z|")
   
   if(print.level >= 1){
    # cat("\n---------------")
    # max.name.length <- max(nchar(row.names(output)))
    est.rez.left <- floor( (max.name.length+42-22) / 2 )
    est.rez.right <- max.name.length+42-22 - est.rez.left
    cat("\n",rep("-", est.rez.left)," Estimation results: ",rep("-", est.rez.right),"\n\n", sep ="")
    # cat("\n--------------- Estimation results: --------------\n\n", sep = "")
    .printgtresfhet(output, digits = digits, Kb, Kv0, Ku0, Kvi, Kui, na.print = "NA", max.name.length)
   }
   obj$results <- output
   colnames(obj$vcov) <- rownames(obj$vcov) <- names(obj$par) <- names(obj$gr) <- rownames(obj$results) <- coef.names.full#[toincludetheta]
   
   # aux parms: 4comp --------------------------------------------------------
   
   my.sig   <- unname( exp( obj$par["ln_sig"] ) )
   my.lmd   <- unname( exp( obj$par["ln_lmd"] ) )
   my.sigv0 <- unname( exp( obj$par["lnVARv0"] ) )
   my.sigu0 <- unname( exp( obj$par["lnVARu0"] ) )

   sv <- my.sig / sqrt(1 + my.lmd * my.lmd)
   su <- sv * my.lmd
   
   # prepare
   # eit : NT x 1
   # svi2 : Ti x 1
   # sv02 :  1 x 1
   # sui2 : Ti x 1
   # su02 :  1 x 1
   
   # define indices
   m.begin <- c(1,(cumsum(YXZ$t0)+1)[-n])
   m.end <- cumsum(YXZ$t0)
   
   my.rez <- as.vector(y - X %*% obj$par[1:Kb])
   obj$resid <- my.rez
   
   sui2 <- rep(su, nt)
   svi2 <- rep(sv, nt)
   sv02 <- rep(my.sigv0, n)
   su02 <- rep(my.sigu0, n)
   
   # Calculate efficiencies --------------------------------------------------
   
   myeff <- ifelse(prod, "technical", "cost")
   # ndots <- ifelse(prod, "....", ".........")
   clceff <- paste0(" Calculating ",myeff," efficiencies: ")
   # print(nchar(clceff))
   if(print.level >= 1){
    # cat("\n-----------")
    # cat(" Calculating ",myeff," efficiencies: ",ndots,"", sep = "")
    est.cle.left <- floor( (max.name.length+42-nchar(clceff)-1) / 2 )
    est.cle.right <- max.name.length+42-nchar(clceff) -1- est.cle.left
    # print(c(max.name.length+42, nchar(clceff),est.cle.left),est.cle.right)
    cat("\n",rep("-", est.cle.left),clceff,rep("-", est.cle.right),"\n\n", sep ="")
   }
   
   # cat.print(length(su02))
   # cat.print(length(sv02))

   time.15 <- proc.time()
   
   # te_i0 <- te_it <- NULL
   # R1 <- ifelse(R < 500, 500, R)
   te_i0 <- ghkReLmd <- rep(NA, n)
   te_it <- rep(NA, nt)
   for(i in seq_len(n) ){
    sample.i <- (m.begin[i]):(m.end[i])
    # cat("id is ",idvar[ (m.begin[i])],"\n", sep = "")
    #    print(idvar[ (m.begin[i]):(m.end[i]) ])
    #    if(i == 1){
    #      print(svi2[ (m.begin[i]):(m.end[i]) ])
    #    }
    t1 <- tryCatch( .e_exp_tu (ei = my.rez[ sample.i ], id = ids[i], sui2=sui2[ sample.i ], svi2 = svi2[ sample.i ], sv02 = sv02[i], su02 = su02[i], prod = prod, simtype = simtype_GHK, R = R_GHK, inv.tol = .Machine$double.eps, print.level = print.level, seed = seed, random.primes = random.primes, halton.base = halton.base ), error = function(e) e )
    if (inherits(t1, "error")){
     print(t1)
     # te_i0 <- c(te_i0, NA )
     # te_it <- c(te_it, rep(NA, YXZ$t0[i]) )
     # nothing to do
    } else {
     # te_i0 <- c(te_i0, t1$te_pers )
     # te_it <- c(te_it, t1$te_resi )
     # if(!is.na(t1$range_star) & !is.na(t1$ghkReLmd)){
     # if(t1$range_star < 0.5 & t1$ghkReLmd > 1e-7){
     te_i0[i] <- t1$te_pers
     ghkReLmd[i] <- t1$ghkReLmd
     te_it[ sample.i ] <- t1$te_resi
     # }
     # }
     # else remains NA
    }
    rm(t1)
   }
   # te_i0 <- ifelse(te_i0 > 1, NA, te_i0)
   # te_it <- ifelse(te_it > 1, NA, te_it)
   te_cs_long <- rep(te_i0, YXZ$t0)
   ghkReLmd_long <- rep(ghkReLmd, YXZ$t0)
   te_over <- te_cs_long * te_it
   #  print(YXZ$t0)
   #  print(sum(YXZ$t0))
   #  print(length(te_i0))
   #  print(te_i0)
   #  print(length(te_it))
   
   
   time.16 <- proc.time()
   # timing((time.16-time.15)[1])
   
   if(print.level >= 1){
    .timing((time.16-time.15)[3], paste("",clceff," completed in\n", sep = ""))
    # cat("____________________________________________________\n")
   }
   
   if(print.level >= 2){
    
    sum.te <- paste0(" Summary of ",myeff," efficiencies: ")
    est.cle.left <- floor( (max.name.length+42-nchar(sum.te)-1) / 2 )
    est.cle.right <- max.name.length+42-nchar(sum.te) -1- est.cle.left
    if(est.cle.left <= 0) est.cle.left <- 1
    if(est.cle.right <= 0) est.cle.right <- 1
    cat("\n",rep("-", est.cle.left),sum.te,rep("-", est.cle.right),"\n\n", sep ="")
    .su(list(te_i0,te_it,te_over), print = TRUE, width = 5, format = "fg", drop0trailing = FALSE, names = c("Persistent", "Transient", "Overall"))
    
   }
   
   # return eff --------------------------------------------------------------
   
   cond1 <- YXZ$old.sorted
   
   temp <- list(call = match.call(), model = paste0("sf_p_4comp_homosk"), coef = obj$par, table = output, vcov = obj$vcov, ll = obj$ll, lmtol = lmtol, LM = obj$gHg, convergence = obj$converged, esttime = est.time.sec, prod = prod, Kb = Kb, fitted = as.vector(unname(X%*%obj$par[1:Kb])), yobs = as.vector(y)[cond1], xobs = X[cond1,], n = n, nt = nt, dat.descr = dat.descr, esample = esample)
   
   
   if(prod){
    temp$efficiencies_i <- data.frame(id = ids, te_i = te_i0, su02, sv02)
    temp$efficiencies_it <- data.frame(id = YXZ$itvar[,1], t = YXZ$itvar[,2], te_i = te_cs_long, te_it = te_it, te_over = te_over, svi2, sv02 = rep(sv02, YXZ$t0), sui2, su02 = rep(su02, YXZ$t0))
    temp$te_i0 <- te_cs_long[cond1]
    temp$te_it <- te_it[cond1]
    temp$te_over <- te_over[cond1]
   }
   else {
    temp$efficiencies_i <- data.frame(id = ids, ce_i = te_i0, su02, sv02)
    temp$efficiencies_it <- data.frame(id = YXZ$itvar[,1], t = YXZ$itvar[,2], ce_i = te_cs_long, ce_it = te_it, ce_over = te_over, svi2, sv02 = rep(sv02, YXZ$t0), sui2, su02 = rep(su02, YXZ$t0))
    temp$ce_i0 <- te_cs_long[cond1]
    temp$ce_it <- te_it[cond1]
    temp$ce_over <- te_over[cond1]
   }
   
   # print(temp)
   
   if (obj$converged == 2){
    temp$delta_rel <- obj$delta_rel
    temp$ltol <- obj$ltol
   }
   if (obj$converged == 3){
    temp$theta_rel_ch <- obj$theta_rel_ch
    temp$steptol <- sqrt(.Machine$double.eps)
   }
   
  }
  
  class(temp) <- "npsf"
  
  return(temp)
  
 }
 
}
