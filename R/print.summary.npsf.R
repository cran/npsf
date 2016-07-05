print.summary.npsf <- function( x, digits = NULL, print.level = NULL, ... ) {
 # digits
 if( is.null( digits ) ) {
  digits <- 4
 }
 # print level
 if( is.null( print.level ) ) {
  print.level <- 3
 }
 # width of the window
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
 # model
 if(is.null(x$model) | is.null(x)){
  stop( gettextf("There is no model in '%s'.", deparse(substitute(x))) )
 }
 mymodel <- x$model
 # possible.models <- c("teradial", "tenonradial", "teradialbc", "nptestrts", "nptestind", "sfsc")
 valid.model <- mymodel %in% c("teradial", "tenonradial", "teradialbc", "nptestrts", "nptestind", "sfsc") # possible.models
 if(!valid.model){
  stop( gettextf("Not valid model in '%s'.", deparse(substitute(x))) )
 }
 # select the model
 # tenoradial or teradial -------------------------------------------------------
 if(mymodel == "tenonradial" | mymodel == "teradial"){
  # begin tenoradial or teradial
  if(print.level >= 1 & winw > 50){
   # DEA
   mymesage <- paste("",ifelse(mymodel == "tenonradial", "Nonradial (Russell)", "Radial (Debrue-Farrell)")," ",x$base,"-based measures of technical efficiency under assumption of ",x$rts," technology are computed for the following data:\n", sep = "")
   cat("",unlist(strsplit(mymesage, " ")),"", sep = " ", fill = winw-10 )
   cat("  Number of data points (K) = ",x$K,"\n", sep = "")
   cat("  Number of outputs     (M) = ",x$M,"\n", sep = "")
   cat("  Number of inputs      (N) = ",x$N,"\n", sep = "")
   # Reference
   if(is.null(x$esample.ref)){
    mymesage <- paste("\nData for reference set are not provided. Reference set is formed by ",x$K," data ",ngettext(x$K, "point", "point(s)"), " for which measures of technical efficiency are computed", sep = "")
    cat("",unlist(strsplit(mymesage, " ")),"", sep = " ", fill = winw-10 )
   }
   else {
    mymesage <- paste("\nReference set is formed by ",sum(x$esample.ref)," provided reference data ",ngettext(sum(x$esample.ref), "point", "point(s)"), "", sep = "")
    cat("",unlist(strsplit(mymesage, " ")),"", sep = " ", fill = winw-10 )
   }
  }
  if(print.level >= 2){
   cat(paste("",rep("_", (winw-10)/1),"", sep = ""), "\n", sep = "")
   cat("Summary of efficiencies:\n\n", sep = "")
   if(mymodel == "tenonradial"){
    RM <- x$te
    .su(RM, print = TRUE)
   } else {
    DF <- x$te
    .su(DF, print = TRUE)
   }
  }
  # end tenoradial or teradial
 }
 # nptestind --------------------------------------------------------------------
 else if(mymodel == "nptestind"){
  # begin nptestind
  if(print.level >= 1 & winw > 50){
   # DEA
   mymesage <- paste("Radial (Debrue-Farrell) ",x$base,"-based measures of technical efficiency under assumption of ",x$rts," technology are computed for the following data:\n", sep = "")
   cat("",unlist(strsplit(mymesage, " ")),"", sep = " ", fill = winw-10 )
   cat("  Number of data points (K) = ",x$K,"\n", sep = "")
   cat("  Number of outputs     (M) = ",x$M,"\n", sep = "")
   cat("  Number of inputs      (N) = ",x$N,"\n", sep = "")
   mymesage <- paste("\nReference set is formed by ",x$K," data point(s), for which measures of technical efficiency are computed", sep = "")
   cat("",unlist(strsplit(mymesage, " ")),"", sep = " ", fill = winw-10 )
   # Null
   if ( x$base == "input" ){
    inps <- ngettext(x$N, "input", "mix of inputs")
    mymesage <- paste("Ho: input-based measure of technical efficiency and ",inps," are independent\n", sep = "")
   }
   else {
    outs <- ngettext(x$M, "output", "mix of outputs")
    # 	  cat(" Test: Ho: output-based measure of technical efficiency and\n", sep = "")
    # 	  cat("           ",outs," are independent \n\n", sep = "")
    mymesage <- paste("\nHo: input-based measure of technical efficiency and ",outs," are independent\n", sep = "")
   }
   cat("", paste("",rep("_", (winw-10)/1),"", sep = ""), sep = "")
   cat("\n          Test of independence\n", sep = "")
   cat("",unlist(strsplit(mymesage, " ")),"", sep = " ", fill = winw-10 )
  }
  # Bootstrap
  cat("Bootstrapping test statistic T4n (",x$reps," replications)\n", sep = "")
  # Result of the test
  if( x$base == "input" ){
   mymesage <- paste("\np-value of the Ho that input-based measure of technical efficiency and ",inps," are independent = ",formatC(x$pval, digits = 4, format = "f"),":", sep = "")
  }
  else {
   mymesage <- paste("\np-value of the Ho that output-based measure of technical efficiency and ",outs," are independent = ",formatC(x$pval, digits = 4, format = "f"),":", sep = "")
  }
  cat("",unlist(strsplit(mymesage, " ")),"", sep = " ", fill = winw-10 )
  
  mymesage <- paste("\n",ifelse(x$pval <= x$alpha, "Heterogeneous", "Homogeneous")," bootstrap ",ifelse(x$pval <= x$alpha, "should", "can")," be used when performing ",x$base,"-based technical efficiency measurement under assumption of ",x$rts," technology", sep = "")
  cat("",unlist(strsplit(mymesage, " ")),"", sep = " ", fill = winw-10 )
  # end nptestind
 }
 # teradialbc -------------------------------------------------------------------
 else if(mymodel == "teradialbc"){
  # begin teradialbc
  if(print.level >= 1 & winw > 50){
   # DEA
   mymesage <- paste("",ifelse(mymodel == "tenonradial", "Nonradial (Russell)", "Radial (Debrue-Farrell)")," ",x$base,"-based measures of technical efficiency under assumption of ",x$rts," technology are computed for the following data:\n", sep = "")
   cat("",unlist(strsplit(mymesage, " ")),"", sep = " ", fill = winw-10 )
   cat("  Number of data points (K) = ",x$K,"\n", sep = "")
   cat("  Number of outputs     (M) = ",x$M,"\n", sep = "")
   cat("  Number of inputs      (N) = ",x$N,"\n", sep = "")
   # Reference
   if(is.null(x$esample.ref)){
    mymesage <- paste("\nData for reference set are not provided. Reference set is formed by ",x$K," data ",ngettext(x$K, "point", "point(s)"), " for which measures of technical efficiency are computed", sep = "")
    cat("",unlist(strsplit(mymesage, " ")),"", sep = " ", fill = winw-10 )
   }
   else {
    mymesage <- paste("\nReference set is formed by ",sum(x$esample.ref)," provided reference data ",ngettext(sum(x$esample.ref), "point", "point(s)"), "", sep = "")
    cat("",unlist(strsplit(mymesage, " ")),"", sep = " ", fill = winw-10 )
   }
  }
  # bootstrap
  cat("", paste("",rep("_", (winw-10)/1),"", sep = ""), sep = "")
  mymesage <- paste("\nBootstrapping reference set formed by ",x$Kr," ",ifelse(is.null(x$ref), "reference data", "data")," point(s) and computing radial (Debreu-Farrell) ",x$base,"-based measures of technical efficiency under assumption of ",x$rts," technology for each of ",x$K," data point(s) realtive to the bootstrapped reference set\n", sep = "")
  cat("",unlist(strsplit(mymesage, " ")),"", sep = " ", fill = winw-10 )
  cat(x$boot.type)
  # warnings
  if(!is.null(x$warnings)){
   for(warn in x$warnings){
    warning(warn, call. = FALSE)
   }
  }
  # end teradialbc
 }
 # nptestrts --------------------------------------------------------------------
 if(mymodel == "nptestrts"){
  # begin nptestrts
  if(print.level >= 1 & winw > 50){
   mymesage <- paste("",ifelse(mymodel == "tenonradial", "Nonradial (Russell)", "Radial (Debrue-Farrell)")," ",x$base,"-based measures of technical efficiency under assumption of CRS, NiRS, and VRS technology are computed for the following data:\n", sep = "")
   cat("",unlist(strsplit(mymesage, " ")),"", sep = " ", fill = winw-10 )
   cat("  Number of data points (K) = ",x$K,"\n", sep = "")
   cat("  Number of outputs     (M) = ",x$M,"\n", sep = "")
   cat("  Number of inputs      (N) = ",x$N,"\n", sep = "")
   
   mymesage <- paste("\nReference set is formed by ",x$K," data point(s), for which measures of technical efficiency are computed", sep = "")
   cat("",unlist(strsplit(mymesage, " ")),"", sep = " ", fill = winw-10 )
   # test 1
   cat("", paste("",rep("_", (winw-10)/1),"", sep = ""), sep = "")
   cat("\n          Test #1\n\n", sep = "")
   cat(" Ho: mean(F_i^CRS)/mean(F_i^VRS) = 1\n", sep = "")
   cat("   and\n", sep = "")
   cat(" Ho: F_i^CRS/F_i^VRS = 1 for each of ",x$K," data ", ngettext(x$K, "point", "point(s)"), "\n", sep = "")
   mymesage <- paste("\nBootstrapping reference set formed by ",x$K," data point(s) and computing radial (Debreu-Farrell) ",x$base,"-based measures of technical efficiency under assumptions of CRS and VRS technology for each of ",x$K," data point(s) realtive to the bootstrapped reference set\n", sep = "")
   cat("",unlist(strsplit(mymesage, " ")),"", sep = " ", fill = winw-10 )
   cat(x$boot.type)
   mymesage <- paste("\np-value of the Ho that mean(F_i^CRS)/mean(F_i^VRS) = 1 (Ho that the global technology is CRS) = ",formatC(x$pGlobalCRS, digits = 4, format = "f"),":\n\nmean(hat{F_i^CRS})/mean(hat{F_i^VRS}) = ",formatC(x$sefficiencyMean, digits = 4, format = "f")," ",ifelse(x$pGlobalCRS>x$alpha, "is not", "is")," statistically greater than 1 at the ",x$alpha*100,"% significance level", sep = "")
   cat("",unlist(strsplit(mymesage, " ")),"", sep = " ", fill = winw-10 )
  }
  if(!is.null(x$nsefficient)){
   if(x$nsefficient == x$K ){
    cat("\nAll data points are scale efficient\n", sep = "")
   }
  }
  # test 2
  if(!is.null(x$pGlobalNRS)){
   if(x$pGlobalCRS > x$alpha){
    # begin CRS is not rejected so perform only individual tests
    if ( !is.null(x$sineffdrs) ){
     cat("", paste("",rep("_", (winw-10)/1),"", sep = ""), sep = "")
     cat("\n          Test #2\n\n", sep = "")
     cat(" Ho: F_i^CRS/F_i^VRS = 1 for each of ",x$K-x$nsefficient," data ", ngettext(x$K, "point", "point(s)"), "\n", sep = "")
     mymesage <- paste("\nBootstrapping reference set formed by ",x$K," data point(s) and computing radial (Debreu-Farrell) ",x$base,"-based measures of technical efficiency under assumptions of CRS and VRS technology for each of ",x$K-x$nsefficient," data point(s) realtive to the bootstrapped reference set\n", sep = "")
     cat("",unlist(strsplit(mymesage, " ")),"", sep = " ", fill = winw-10 )
    }
    # end # CRS is not rejected so perform only individual tests
   } else {
    # begin CRS is rejected, perform NiRS global and individual tests
    cat("", paste("",rep("_", (winw-10)/1),"", sep = ""), sep = "")
    cat("\n          Test #2\n\n", sep = "")
    cat(" Ho: mean(F_i^NiRS)/mean(F_i^VRS) = 1\n", sep = "")
    if (x$K-x$nsefficient > 0) {
     cat("   and\n", sep = "")
     cat(" Ho: F_i^NiRS/F_i^VRS = 1 for each of ",x$K-x$nsefficient," data ", ngettext(x$K, "point", "point(s)"), "\n", sep = "")
    }
    mymesage <- paste("\nBootstrapping reference set formed by ",x$K," data point(s) and computing radial (Debreu-Farrell) ",x$base,"-based measures of technical efficiency under assumptions of CRS and VRS technology for each of ",x$K," data point(s) realtive to the bootstrapped reference set\n", sep = "")
    cat("",unlist(strsplit(mymesage, " ")),"", sep = " ", fill = winw-10 )
    # next
    if(!is.null(x$nrsOVERvrsMean)){
     mymesage <- paste("p-value of the Ho that mean(F_i^NiRS)/mean(F_i^VRS) = 1 (Ho that the global technology is NiRS) = ",formatC(x$pGlobalNRS, digits = 4, format = "f"),":\n\nmean(hat{F_i^NiRS})/mean(hat{F_i^VRS}) = ",formatC(x$nrsOVERvrsMean, digits = 4, format = "f")," ",ifelse(x$pGlobalNRS>x$alpha, "is not", "is")," statistically greater than 1 at the ",x$alpha*100,"% significance level", sep = "")
     cat("",unlist(strsplit(mymesage, " ")),"", sep = " ", fill = winw-10 )
    }
    # end CRS is rejected, perform NiRS global and individual tests
   }
  }
  # warnings
  if(!is.null(x$warnings)){
   for(warn in x$warnings){
    warning(warn, call. = FALSE)
   }
  }
  # end nptestrts
 }
 # sf sross-sectional -----------------------------------------------------------
 if(mymodel == "sfsc"){
  # begin sf sross-sectional
  if(print.level >= 1 & winw > 50){
   if(x$LM < x$lmtol){
    cat("Convergence given g*inv(H)*g' = ",formatC(x$LM,digits=1,format="e")," < lmtol(",formatC(x$lmtol,digits=1,format="e"),")\n", sep = "")
   }
   else {
    cat("'optim' did it\n", sep = "")
   }
   .timing(x$esttime, "Log likelihood maximization completed in ")
   cat("Final log likelihood = ",formatC(x$loglik,digits=7,format="f"),"\n",  sep = "")
   # model
   max.name.length <- max(nchar(row.names(x$table)))
   cat("",rep("_", max.name.length+42-1),"", "\n", sep = "")
   if(x$prod){
    cat("\nCross-sectional stochastic (production) frontier model\n",  sep = "")} 
   else {
    cat("\nCross-sectional stochastic (cost) frontier model\n",  sep = "")
   }
   cat("Distributional assumptions:\n\n", sep = "")
   Assumptions <- rep("heteroskedastic",2)
   if(x$kv==1){
    Assumptions[1] <- "homoskedastic"
   } 
   if(x$ku==1){
    Assumptions[2] <- "homoskedastic"
   }
   Distribution = c("normal ", "half-normal ")
   if(x$distribution == "t"){
    Distribution[2] <- "truncated-normal "
   }
   a1 <- data.frame(
    Component = c("Random noise: ","Inefficiency: "),
    Distribution = Distribution,
    Assumption = Assumptions
   )
   print(a1, quote = FALSE, right = FALSE)
   cat("\nNumber of observations = ", nrow(x$eff), "\n", sep = "")
   # estimation results
   est.rez.left <- floor( (max.name.length+42-22) / 2 )
   est.rez.right <- max.name.length+42-22 - est.rez.left
   cat("\n",rep("-", est.rez.left)," Estimation results: ",rep("-", est.rez.right),"\n\n", sep ="")
   # cat("\n--------------- Estimation results: --------------\n\n", sep = "")
   .printoutcs(x$table, digits = digits, k = x$k, kv = x$kv, ku = x$ku, kdel = x$kdel, na.print = "NA", dist = x$distribution,max.name.length = max.name.length)
   # efficiencies
   myeff <- ifelse(x$prod, "technical", "cost")
   eff.name <- paste0("Summary of ",myeff," efficiencies")
   len.eff.name <- nchar(eff.name)
   est.eff.left <- floor( (max.name.length+42-len.eff.name-4) / 2 )
   est.eff.right <- max.name.length+42-len.eff.name-4 - est.eff.left
   cat("\n",rep("-", est.eff.left)," ",eff.name,": ",rep("-", est.eff.right),"\n\n", sep ="")
   eff1 <- x$eff[,2:4]
   colnames(eff1) <- formatC(colnames(eff1), width = 4, flag = "-")
   cat("",rep(" ", est.eff.left+1),"JLMS:= exp(-E[ui|ei])\n", sep = "")
   cat("",rep(" ", est.eff.left+1),"Mode:= exp(-M[ui|ei])\n", sep = "")
   cat("",rep(" ", est.eff.left+3),"BC:= E[exp(-ui)|ei]\n", sep = "")
   cat("\n")
   .su1(eff1, transpose = TRUE)
  }
  # end sf sross-sectional
 }
 
 
 # if(x$LM < x$lmtol){
 #    cat("\nConvergence given g*inv(H)*g' = ",formatC(x$LM,digits=1,format="e")," < lmtol(",x$lmtol,")\n", sep = "")
 #  }
 #  else {
 #    cat("\nCriterion g*inv(H)*g' = ",formatC(x$LM,digits=1,format="e")," > lmtol(",x$lmtol,")\n", sep = "")
 #    if(x$bhhh){
 #      cat("Note that Hessian is computed as outer product (BHHH)\n", sep = "")
 #      cat("Criterion g'g = ",obg$gg,"\n", sep = "")
 #    }
 #    warning("Convergence given g*inv(H)*g' is still not reached; one of optim's convergence criteria is used", call. = FALSE)
 #  }
 #  .timing(x$esttime, "Log likelihood maximization completed in ")
 #  cat("Log likelihood = ",formatC(x$ll,digits=4,format="f"),"\n",  sep = "")
 #  cat("____________________________________________________\n")

  invisible( x )
}