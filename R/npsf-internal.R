.prepareYX <- function(formula, data, subset, rts = c("C", "NI", "V"),
                       base = c("output", "input"), ref = NULL,
                       data.ref = NULL, subset.ref = NULL, print.level = 1,
                       type = "RM", winw = 50, rts.subst = NULL, sysnframe = 1, ...)
{
 # get y and x matrices

 # mf0 <- match.call(expand.dots = FALSE, call = sys.call(which = 1))
 # needed.frame <- sys.nframe() - 1
 needed.frame <- sys.nframe() - sysnframe
 # cat("sys.nframe() = ",sys.nframe(),"\n")
 # cat("needed.frame = ",needed.frame,"\n")
 # cat("sysnframe = ",sysnframe,"\n")
 # this is the entire call
 # this one is from previous call
 # mf0 <- match.call(expand.dots = FALSE, call = sys.call(sys.parent(n = needed.frame)))
 # print(mf0)
 # this one is from the initial call (where tenonradial might have been called)
 mf0 <- match.call(expand.dots = FALSE, call = sys.call(sys.parent(n = needed.frame)))
 # print(mf0)
 # print( match.call(expand.dots = FALSE, call = sys.call(sys.parent(n = sysnframe+2))) )
 # for(i in 1:1){
 #  cat("i=",i,"\n")
 #  mf0 <- match.call(expand.dots = FALSE, call = sys.call(sys.parent(n = i)))
 #  mf <- mf0
 #  m <- match(c("formula", "data", "subset"), names(mf), 0L)
 #  mf <- mf[c(1L, m)]
 #  mf[[1L]] <- as.name("model.frame")
 #  # mf$formula <- Formula( formula )
 #  # mf <- eval(mf, parent.frame(n = 1))
 #  print(mf)
 # }
 # cat("end\n")

 # if(print.level >= 1) print(mf0)

 # check if it is a matrix
 datasupplied <- !(match("data", names(mf0), 0) == 0)
 # subsetsupplied <- !(match("subset", names(mf0), 0) == 0)
 # print(subsetsupplied)
 # cat("1\n")

 # if data are supplied
 if(datasupplied){
  # begin get a logical vector equal TRUE if !missing

  # first using data and subset to get x without NA
  mf <- mf0
  m <- match(c("formula", "data", "subset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- as.name("model.frame")
  # mf$formula <- Formula( formula )
  # print(mf)
  mf <- eval(mf, sys.frame(sys.parent(n = needed.frame)))
  # mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  x <- as.matrix(model.matrix(mt, mf))
  # print(x)
  # now get the names in the entire data
  # esample <- seq_len( nrow(data) ) %in% as.numeric(rownames(x))
  esample <- rownames(data) %in% rownames(x)

  # print(rownames(data) )
  # print(rownames(x))
  #
  # print(table(esample))
  # end get a logical vector equal TRUE if !missing

  # get the data
  y <- as.matrix(model.part(Formula(formula), data = data[esample,], lhs = 1))
  x <- as.matrix(model.matrix(Formula(formula), data = data[esample,], rhs = 1)[,-1])

  # print(y)
  # print(x)
 }
 # if data are not supplied
 else {
  mf <- mf0
  subsetsupplied <- !(match("subset", names(mf), 0) == 0)
  if(subsetsupplied) stop("Subset with matrices is not allowed")
  m <- match(c("formula"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- as.name("model.frame")
  # mf <- eval(mf, parent.frame())
  mf <- eval(mf, sys.frame(sys.parent(n = needed.frame)))
  mt <- attr(mf, "terms")
  y <- as.matrix(model.response(mf))
  x <- as.matrix(model.matrix(mt, mf)[,-1])
  # get a logical vector equal TRUE if !missing
  with.na <- model.frame(mt, na.action = na.pass)
  esample <- rownames(with.na) %in% rownames(x)
  # print(table(esample))
  # print(y)
  # print(x)
 }

 if( !is.numeric(y) ) stop("Some of the outputs are not numeric.")
 if( !is.numeric(x) ) stop("Some of the inputs are not numeric.")
 
 if(type == "RM"){
  # check negative
  x.negative <- apply(x, MARGIN = 1, FUN = function(q) sum(q<0)) >= 1
  if(sum(x.negative) == nrow(x)){
   stop("Negative values in at least one of the inputs for each data point", call. = FALSE)
  }
  y.negative <- apply(y, MARGIN = 1, FUN = function(q) sum(q<0)) >= 1
  if(sum(y.negative) == nrow(y)){
   stop("Negative values in at least one of the outputs for each data point", call. = FALSE)
  }
  # deal with negative
  if(sum(x.negative) > 0){
   warning(paste0("There are negative values in inputs for ",sum(x.negative)," data points; these data points will not be considered"), call. = FALSE, immediate. = TRUE)
   cat("\n")
  }
  if(sum(y.negative) > 0){
   warning(paste0("There are negative values in outputs for ",sum(y.negative)," data points; these data points will not be considered"), call. = FALSE, immediate. = TRUE)
   cat("\n")
  }
  x <- x[!x.negative & !y.negative,,drop = FALSE]
  y <- y[!x.negative & !y.negative,,drop = FALSE]
 } else {
  # check nonpositive
  x.nonpositive <- apply(x, MARGIN = 1, FUN = function(q) sum(q<=0)) >= 1
  if(sum(x.nonpositive) == nrow(x)){
   stop("Nonpositive values in at least one of the inputs for each data point", call. = FALSE)
  }
  y.nonpositive <- apply(y, MARGIN = 1, FUN = function(q) sum(q<00)) >= 1
  if(sum(y.nonpositive) == nrow(y)){
   stop("Nonpositive values in at least one of the outputs for each data point", call. = FALSE)
  }
  # deal with negative
  if(sum(x.nonpositive) > 0){
   warning(paste0("There are nonpositive values in inputs for ",sum(x.nonpositive)," data points; these data points will not be considered"), call. = FALSE, immediate. = TRUE)
   cat("\n")
  }
  if(sum(y.nonpositive) > 0){
   warning(paste0("There are nonpositive values in outputs for ",sum(y.nonpositive)," data points; these data points will not be considered"), call. = FALSE, immediate. = TRUE)
   cat("\n")
  }
  x <- x[(!x.nonpositive & !y.nonpositive),,drop = FALSE]
  y <- y[(!x.nonpositive & !y.nonpositive),,drop = FALSE]
 }
 
 rts <- rts[1]
 base <- base[1]

 rts1 <- tolower(substr(rts, 1,1 ))
 base1 <- tolower(substr(base, 1,1 ))

 if(is.null(rts.subst)){
  if(rts1 == "c" | rts == 3){
   myrts <- 3
   myrts1 <- "CRS"
  } else if (rts1 == "n" | rts == 2){
   myrts <- 2
   myrts1 <- "NIRS"
  } else if (rts1 == "v" | rts == 1){
   myrts <- 1
   myrts1 <- "VRS"
  } else {
   stop("invalid 'rts'; 'rts' must be 'CRS', 'NIRS', or 'VRS'")
  }
 }
 else {
  myrts1 <- rts.subst
 }

 if (base1 == "o" | base == 2){
  mybase <- 2
  mybase1 <- "output"
 } else if (base1 == "i" | base == 1){
  mybase <- 1
  mybase1 <- "input"
 } else {
  stop("invalid 'base'; 'base' must be 'output' or 'input'")
 }

 # for printing
 if(print.level >= 1 & winw > 50){
  mymesage <- paste("\n",ifelse(type == "RM", "Nonradial (Russell)", "Radial (Debrue-Farrell)")," ",mybase1,"-based measures of technical efficiency under assumption of ",myrts1," technology are computed for the following data:\n", sep = "")
  cat("",unlist(strsplit(mymesage, " ")),"", sep = " ", fill = winw-10 )
  # 		cat("\n ",ifelse(type == "RM", "Nonradial (Russell)", "Radial Debrue-Farrell")," ",mybase1,"-based measures of technical efficiency \n", sep = "")
  # 		cat(" under assumption of ",myrts1," technology are computed for the \n", sep = "")
  # 		cat(" following data:\n\n", sep = "")
  cat("  Number of data points (K) = ",nrow(y),"\n", sep = "")
  cat("  Number of outputs     (M) = ",ncol(y),"\n", sep = "")
  cat("  Number of inputs      (N) = ",ncol(x),"\n", sep = "")
 }

 # get y_ref and x_ref matrices

 if(!is.null(ref)){
  mf <- mf0

  # check if it is a matrix
  datasupplied <- !(match("data.ref", names(mf), 0) == 0)

  # if data are supplied
  if(datasupplied){
   N_all_ref <- nrow(data.ref)
   if(N_all_ref == 0) warning("Provided data for reference set does not have a signle data point", call. = FALSE)

   # begin get a logical vector equal TRUE if !missing

   # first using data and subset to get x without NA
   mf <- mf0
   m <- match(c("ref", "data.ref", "subset.ref"), names(mf), 0L)
   mf <- mf[c(1L, m)]
   # change the names for eval
   names(mf)[which(names(mf) == "ref")] <- "formula"
   names(mf)[which(names(mf) == "data.ref")] <- "data"
   names(mf)[which(names(mf) == "subset.ref")] <- "subset"
   mf[[1L]] <- as.name("model.frame")
   # mf$formula <- Formula( formula )
   mf <- eval(mf, parent.frame())
   mt <- attr(mf, "terms")
   x_ref <- as.matrix(model.matrix(mt, mf))
   if(length(x_ref) == 0){
    warning(" Given 'subset' reference set contains zero data points", call. = FALSE)
    # warning(" Reference set will be based on observations,")
    # warning(" for which efficiency is to be computed")
    y_ref <- y
    x_ref <- x
    esample_ref <- NULL
   } else {
    # now get the names in the entire data
    esample_ref <- seq_len( N_all_ref ) %in% as.numeric(rownames(x_ref))
    # print(table(esample_ref))
    # end get a logical vector equal TRUE if !missing
    y_ref <- as.matrix(model.part(Formula(ref), data = data.ref[esample_ref,], lhs = 1))
    x_ref <- as.matrix(model.matrix(Formula(ref), data = data.ref[esample_ref,], rhs = 1)[,-1])
   }
   # print(y_ref)
   # print(x_ref)
  }
  # if data are not supplied
  else {
   mf <- mf0
   subsetsupplied <- !(match("subset.ref", names(mf), 0) == 0)
   if(subsetsupplied) stop("Subset for reference set with matrices is not allowed. \n   Or: \nOption 'data.ref' must be provided if option 'ref' is provided", call. = FALSE)
   m <- match(c("ref", "data.ref", "subset.ref"), names(mf), 0L)
   mf <- mf[c(1L, m)]
   names(mf)[which(names(mf) == "ref")] <- "formula"
   mf[[1L]] <- as.name("model.frame")
   mf <- eval(mf, parent.frame())
   mt <- attr(mf, "terms")
   y_ref <- as.matrix(model.response(mf))
   x_ref <- as.matrix(model.matrix(mt, mf)[,-1])
   # print(y_ref)
   # print(x_ref)
   # get a logical vector equal TRUE if !missing
   with.na <- model.frame(mt, na.action = na.pass)
   esample_ref <- rownames(with.na) %in% rownames(x_ref)
   # print(table(esample_ref))
  }
  # if reference set not provided
  # for printing
  # 		if(print.level >= 1){
  # 			if(!is.null(esample_ref)){
  # 				cat("\n Reference set is composed as follows:\n\n", sep = "")
  # 				cat(" Number of observations (K) = ",nrow(y_ref),"\n", sep = "")
  # 				cat(" Number of outputs      (M) = ",ncol(y_ref),"\n", sep = "")
  # 				cat(" Number of inputs       (N) = ",ncol(x_ref),"\n\n", sep = "")
  # 			}
  # 		}
  if(!is.null(y_ref)){
   if(ncol(y_ref) != ncol(y)) stop("Number of outputs for data points and reference set must be the same.")
   if(ncol(x_ref) != ncol(x)) stop("Number of inputs for data points and reference set must be the same.")
  }
  
  if(type == "RM"){
   # check negative
   x.negative <- apply(x_ref, MARGIN = 1, FUN = function(q) sum(q<0)) >= 1
   if(sum(x.negative) == nrow(x_ref)){
    stop("Negative values in at least one of the reference inputs for each data point", call. = FALSE)
   }
   y.negative <- apply(y_ref, MARGIN = 1, FUN = function(q) sum(q<0)) >= 1
   if(sum(y.negative) == nrow(y_ref)){
    stop("Negative values in at least one of the reference outputs for each data point", call. = FALSE)
   }
   # deal with negative
   if(sum(x.negative) > 0){
    warning(paste0("There are negative values in reference inputs for ",sum(x.negative)," data points; these data points will not be considered"), call. = FALSE, immediate. = TRUE)
    cat("\n")
   }
   if(sum(y.negative) > 0){
    warning(paste0("There are negative values in reference outputs for ",sum(y.negative)," data points; these data points will not be considered"), call. = FALSE, immediate. = TRUE)
    cat("\n")
   }
   x_ref <- x_ref[!x.negative & !y.negative,,drop = FALSE]
   y_ref <- y_ref[!x.negative & !y.negative,,drop = FALSE]
  } else {
   # check nonpositive
   x.nonpositive <- apply(x_ref, MARGIN = 1, FUN = function(q) sum(q<=0)) >= 1
   if(sum(x.nonpositive) == nrow(x)){
    stop("Nonpositive values in at least one of the reference inputs for each data point", call. = FALSE)
   }
   y.nonpositive <- apply(y_ref, MARGIN = 1, FUN = function(q) sum(q<00)) >= 1
   if(sum(y.nonpositive) == nrow(y_ref)){
    stop("Nonpositive values in at least one of the reference outputs for each data point", call. = FALSE)
   }
   # deal with negative
   if(sum(x.nonpositive) > 0){
    warning(paste0("There are nonpositive values in reference inputs for ",sum(x.nonpositive)," data points; these data points will not be considered"), call. = FALSE, immediate. = TRUE)
    cat("\n")
   }
   if(sum(y.nonpositive) > 0){
    warning(paste0("There are nonpositive values in outputs for ",sum(y.nonpositive)," data points; these data points will not be considered"), call. = FALSE, immediate. = TRUE)
    cat("\n")
   }
   x_ref <- x_ref[!x.nonpositive & !y.nonpositive,,drop = FALSE]
   y_ref <- y_ref[!x.nonpositive & !y.nonpositive,,drop = FALSE]
  }
  
 } else {
  y_ref <- y
  x_ref <- x
  esample_ref <- NULL
 }

 if(print.level >= 1  & winw > 50){
  if(is.null(esample_ref)){
   mymesage <- paste("\nData for reference set are not provided. Reference set is formed by ",nrow(y)," data ",ngettext(nrow(y), "point", "point(s)"), " for which measures of technical efficiency are computed", sep = "")
   cat("",unlist(strsplit(mymesage, " ")),"", sep = " ", fill = winw-10 )
   #     cat("\n Data for reference set are not provided. Reference set is formed by ",nrow(y),"\n", sep = "")
   #     cat(" data points, for which measures of technical efficiency are computed.\n\n", sep = "")
  }
  else {
   mymesage <- paste("\nReference set is formed by ",nrow(y_ref)," provided reference data ",ngettext(nrow(y_ref), "point", "point(s)"), "", sep = "")
   cat("",unlist(strsplit(mymesage, " ")),"", sep = " ", fill = winw-10 )
   # cat("\n Reference set is formed by ",nrow(y_ref)," provided reference data points.\n\n", sep = "")
  }
 }

 tymch <- list(y = y, x = x, y.ref = y_ref, x.ref = x_ref, esample = esample, esample.ref = esample_ref, myrts = myrts, rts.string = myrts1, mybase = mybase, base.string = mybase1)
 class(tymch) <- "npsf"
 return(tymch)
}

.prepareYXnoRef <- function(formula, data, subset,
                            rts = c("C", "NI", "V"), rts.subst = NULL,
                            base = c("output", "input"), print.level = 1,
                            type = "RM", winw = 50, ...)
{
 # get y and x matrices

 # mf0 <- match.call(expand.dots = FALSE, call = sys.call(which = 1))
 needed.frame <- sys.nframe() - 1
 mf0 <- match.call(expand.dots = FALSE, call = sys.call(sys.parent(n = needed.frame)))

 # if(print.level >= 1) print(mf0)

 # check if it is a matrix
 datasupplied <- !(match("data", names(mf0), 0) == 0)
 # subsetsupplied <- !(match("subset", names(mf0), 0) == 0)
 # print(subsetsupplied)
 # cat("1\n")

 # if data are supplied
 if(datasupplied){
  # begin get a logical vector equal TRUE if !missing

  # first using data and subset to get x without NA
  mf <- mf0
  m <- match(c("formula", "data", "subset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- as.name("model.frame")
  # mf$formula <- Formula( formula )
  # mf <- eval(mf, parent.frame(n = 1))
  mf <- eval(mf, sys.frame(sys.parent(n = needed.frame)))
  mt <- attr(mf, "terms")
  x <- as.matrix(model.matrix(mt, mf))
  # now get the names in the entire data
  esample <- seq_len( nrow(data) ) %in% as.numeric(rownames(x))
  # print(table(esample))
  # end get a logical vector equal TRUE if !missing

  # get the data
  y <- as.matrix(model.part(Formula(formula), data = data[esample,], lhs = 1))
  x <- as.matrix(model.matrix(Formula(formula), data = data[esample,], rhs = 1)[,-1])

  # print(y)
  # print(x)
 }
 # if data are not supplied
 else {
  mf <- mf0
  subsetsupplied <- !(match("subset", names(mf), 0) == 0)
  if(subsetsupplied) stop("Subset with matrices is not allowed")
  m <- match(c("formula", "data", "subset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- as.name("model.frame")
  # mf <- eval(mf, parent.frame())
  mf <- eval(mf, sys.frame(sys.parent(n = needed.frame)))
  mt <- attr(mf, "terms")
  y <- as.matrix(model.response(mf))
  x <- as.matrix(model.matrix(mt, mf)[,-1])
  # get a logical vector equal TRUE if !missing
  with.na <- model.frame(mt, na.action = na.pass)
  esample <- rownames(with.na) %in% rownames(x)
  # print(table(esample))
  # print(y)
  # print(x)
 }

 if( !is.numeric(y) ) stop("Some of the outputs are not numeric.")
 if( !is.numeric(x) ) stop("Some of the inputs are not numeric.")
 
 if(type == "RM"){
  # check negative
  x.negative <- apply(x, MARGIN = 1, FUN = function(q) sum(q<0)) >= 1
  if(sum(x.negative) == nrow(x)){
   stop("Negative values in at least one of the inputs for each data point", call. = FALSE)
  }
  y.negative <- apply(y, MARGIN = 1, FUN = function(q) sum(q<0)) >= 1
  if(sum(y.negative) == nrow(y)){
   stop("Negative values in at least one of the outputs for each data point", call. = FALSE)
  }
  # deal with negative
  if(sum(x.negative) > 0){
   warning(paste0("There are negative values in inputs for ",sum(x.negative)," data points; these data points will not be considered"), call. = FALSE, immediate. = TRUE)
   cat("\n")
  }
  if(sum(y.negative) > 0){
   warning(paste0("There are negative values in outputs for ",sum(y.negative)," data points; these data points will not be considered"), call. = FALSE, immediate. = TRUE)
   cat("\n")
  }
  x <- x[!x.negative & !y.negative,,drop = FALSE]
  y <- y[!x.negative & !y.negative,,drop = FALSE]
 } else {
  # check nonpositive
  x.nonpositive <- apply(x, MARGIN = 1, FUN = function(q) sum(q<=0)) >= 1
  if(sum(x.nonpositive) == nrow(x)){
   stop("Nonpositive values in at least one of the inputs for each data point", call. = FALSE)
  }
  y.nonpositive <- apply(y, MARGIN = 1, FUN = function(q) sum(q<00)) >= 1
  if(sum(y.nonpositive) == nrow(y)){
   stop("Nonpositive values in at least one of the outputs for each data point", call. = FALSE)
  }
  # deal with negative
  if(sum(x.nonpositive) > 0){
   warning(paste0("There are nonpositive values in inputs for ",sum(x.nonpositive)," data points; these data points will not be considered"), call. = FALSE, immediate. = TRUE)
   cat("\n")
  }
  if(sum(y.nonpositive) > 0){
   warning(paste0("There are nonpositive values in outputs for ",sum(y.nonpositive)," data points; these data points will not be considered"), call. = FALSE, immediate. = TRUE)
   cat("\n")
  }
  x <- x[!x.nonpositive & !y.nonpositive,,drop = FALSE]
  y <- y[!x.nonpositive & !y.nonpositive,,drop = FALSE]
 }

 rts <- rts[1]
 base <- base[1]

 rts1 <- tolower(substr(rts, 1,1 ))
 base1 <- tolower(substr(base, 1,1 ))

 if(is.null(rts.subst)){
  if(rts1 == "c" | rts == 1){
   myrts <- 3
   myrts1 <- "CRS"
  } else if (rts1 == "n" | rts == 2){
   myrts <- 2
   myrts1 <- "NIRS"
  } else if (rts1 == "v" | rts == 3){
   myrts <- 1
   myrts1 <- "VRS"
  } else {
   stop("invalid 'rts'; 'rts' must be 'CRS', 'NIRS', or 'VRS'")
  }
 }
 else {
  myrts1 <- rts.subst
 }

 if (base1 == "o" | base == 2){
  mybase <- 2
  mybase1 <- "output"
 } else if (base1 == "i" | base == 1){
  mybase <- 1
  mybase1 <- "input"
 } else {
  stop("invalid 'base'; 'base' must be 'output' or 'input'")
 }

 # for printing
 #  if(print.level >= 1){
 #   cat("\n Radial (Debreu-Farrell) ",mybase1,"-based measures of technical efficiency \n", sep = "")
 #   cat(" under assumption of CRS, NIRS, and VRS technology are computed for the\n", sep = "")
 #   cat(" following data:\n\n", sep = "")
 #   cat("  Number of data points (K) = ",nrow(y),"\n", sep = "")
 #   cat("  Number of outputs     (M) = ",ncol(y),"\n", sep = "")
 #   cat("  Number of inputs      (N) = ",ncol(x),"\n", sep = "")
 #  }

 if(print.level >= 1 & winw > 50){
  mymesage <- paste("\n",ifelse(type == "RM", "Nonradial (Russell)", "Radial (Debrue-Farrell)")," ",mybase1,"-based measures of technical efficiency under assumption of ",myrts1," technology are computed for the following data:\n", sep = "")
  cat("",unlist(strsplit(mymesage, " ")),"", sep = " ", fill = winw-10 )
  # 		cat("\n ",ifelse(type == "RM", "Nonradial (Russell)", "Radial Debrue-Farrell")," ",mybase1,"-based measures of technical efficiency \n", sep = "")
  # 		cat(" under assumption of ",myrts1," technology are computed for the \n", sep = "")
  # 		cat(" following data:\n\n", sep = "")
  cat("  Number of data points (K) = ",nrow(y),"\n", sep = "")
  cat("  Number of outputs     (M) = ",ncol(y),"\n", sep = "")
  cat("  Number of inputs      (N) = ",ncol(x),"\n", sep = "")

  mymesage <- paste("\nReference set is formed by ",nrow(y)," data point(s), for which measures of technical efficiency are computed", sep = "")
  cat("",unlist(strsplit(mymesage, " ")),"", sep = " ", fill = winw-10 )
 }

 #  if(print.level >= 1){
 #    cat("\n Reference set is formed by ",nrow(y)," data points,\n", sep = "")
 #    cat(" for which measures of technical efficiency are computed.\n\n", sep = "")
 #  }

 tymch <- list(y = y, x = x, esample = esample, mybase = mybase, base.string = mybase1)
 if(is.null(rts.subst)){
  tymch$myrts <- myrts
  tymch$rts.string <- myrts1
 }
 class(tymch) <- "npsf"
 return(tymch)
}

# .su <- function(x, mat.var.in.col = TRUE, digits = 4, probs = c(0.1, 0.25, 0.5, 0.75, 0.9), print = FALSE){
# 
#  xvec2 <- xvec1 <- FALSE
# 
#  if(is.matrix(x)){
#   if(min(dim(x)) == 1){
#    # xvec1 <- TRUE
#    if(which(dim(x) == 1) == 2){
#     mynames <- colnames(x)
#    } else {
#     mynames <- rownames(x)
#     x <- t(x)
#    }
#    # x <- as.vector(x)
#   } else {
#    if(!mat.var.in.col){
#     x <- t(x)
#    }
#    mynames <- colnames(x)
#   }
#   # print(mynames)
#   if(is.null(mynames)) mynames <- paste("Var", seq_len(ncol(x)), sep = "")
#   x <- as.data.frame(x)
#   # print(x)
#   # mynames <- colnames(x)
#  } # end if matrix
# 
#  if(is.vector(x)){
#   xvec2 <- TRUE
#   mynames <- deparse(substitute(x))
#   x <- data.frame(Var1 = x)
#  } # end if vector
# 
#  # cat("nymanes", sep ="")
#  # print(mynames)
# 
#  if(!is.vector(x) & !is.matrix(x) & !is.data.frame(x)){
#   stop("Provide vector, matrix, or data.frame")
#  } else {
#   t1 <- apply(x, 2, function(x) c(Obs = length(x), NAs = sum(is.na(x)), Mean = mean(x, na.rm = TRUE), StDev = sd(x, na.rm = TRUE), IQR = IQR(x, na.rm = TRUE), Min = min(x, na.rm = TRUE), quantile( x, probs = probs, na.rm = TRUE ), Max = max(x, na.rm = TRUE)))
#   # print(t1)
#   # print(mynames)
#   # print(class(t1))
#   # print(dim(t1))
#   if(xvec2 & !xvec1) colnames(t1) <- mynames
#   if(print){
#    tymch <- formatC(t(t1), digits = 4, format = "f", width = 4+1)
#    tymch <- gsub(".0000","",tymch, fixed = TRUE)
#    print(tymch, quote = FALSE)
#   }
#  }
#  class(t1) <- "npsf"
#  return(t(t1))
# }

.my.prettyNum <- function(xx2){
 n <- length(xx2)
 xx2.pn <- prettyNum(xx2)
 my.integers <- !(1:n %in%  grep(".", xx2.pn, fixed = TRUE))
 xx3 <- ifelse(abs(xx2) < 1e-04 & xx2 != 0, formatC(xx2, digits = 1, format = "e", width = 5), ifelse(abs(xx2) >= 1e-04 & abs(xx2) < 10, formatC(xx2, format = "f", digits = 4), ifelse(abs(xx2) >= 10 & abs(xx2) < 100, formatC(xx2, format = "f", digits = 3), ifelse(abs(xx2) >= 100 & abs(xx2) < 1000, formatC(xx2, format = "f", digits = 2), ifelse(abs(xx2) >= 1000, formatC(xx2, format = "e", digits = 1), formatC(xx2, format = "f", digits = 1))))))
 xx3[my.integers] <- xx2.pn[my.integers]
 xx3
}

.su <- function(x, mat.var.in.col = TRUE, digits = 4, probs = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9), print = FALSE, names = NULL, ...){
 
 xvec2 <- xvec1 <- FALSE
 # print(class(x))
 
 if(is.list(x)){
  # cat("this is a list\n")
  mynames <- paste0("Var",1:length(x))
 }
 
 if(is.matrix(x)){
  # cat("this is a matrix\n")
  if(min(dim(x)) == 1){
   # xvec1 <- TRUE
   if(which(dim(x) == 1) == 2){
    mynames <- colnames(x)
    
   } else {
    mynames <- rownames(x)
    x <- t(x)
   }
   # x <- as.vector(x)
  } else {
   if(!mat.var.in.col){
    x <- t(x)
   }
   mynames <- colnames(x)
  }
  
  # print(mynames)
  if(is.null(mynames)) mynames <- paste("Var", seq_len(ncol(x)), sep = "")
  x <- as.data.frame(x)
  # print(x)
  # mynames <- colnames(x)
 } # end if matrix
 
 if(is.vector(x) & !is.list(x)){
  # cat("this is a vector\n")
  xvec2 <- TRUE
  mynames <- deparse(substitute(x))
  x <- data.frame(Var1 = x)
 } # end if vector
 
 # print(mynames)
 # print(names)
 
 if(!is.null(names)){
  if(length(mynames) == length(names)){
   mynames <- names
  }
 }
 
 # cat("nynames\n", sep ="")
 # print(mynames)
 if(is.list(x)){
  t1 <- sapply(x, function(x) c(Obs = length(x), NAs = sum(is.na(x)), Mean = mean(x, na.rm = TRUE), StDev = sd(x, na.rm = TRUE), IQR = IQR(x, na.rm = TRUE), Min = min(x, na.rm = TRUE), quantile( x, probs = probs, na.rm = TRUE ), Max = max(x, na.rm = TRUE)))
  # print(t1)
  # print(mynames)
  # print(class(t1))
  # print(is.integer(t1))
  # print( formatC(t1, format = "f") == formatC(t1, format = "fg") )
  # print(dim(t1))
  # if(xvec2 & !xvec1) colnames(t1) <- mynames
  colnames(t1) <- mynames
  # if(print){
  #  # t2 <- rbind(
  #  #   formatC(t1[c(1:2),,drop = FALSE], digits = digits, format = "fg"),
  #  #   formatC(t1[-c(1:2),,drop = FALSE], digits = digits, format = "f")
  #  # )
  #  print(data.frame(.my.prettyNum(t1)))
  # }
  if(print){
   t2 <- data.frame(.my.prettyNum(t1))
   colnames(t2) <- mynames
   print(t2)
  } 
 } else if(!is.vector(x) & !is.matrix(x) & !is.data.frame(x)){
  stop("Provide list, vector, matrix, or data.frame")
 } else {
  t1 <- apply(x, 2, function(x) c(Obs = length(x), NAs = sum(is.na(x)), Mean = mean(x, na.rm = TRUE), StDev = sd(x, na.rm = TRUE), IQR = IQR(x, na.rm = TRUE), Min = min(x, na.rm = TRUE), quantile( x, probs = probs, na.rm = TRUE ), Max = max(x, na.rm = TRUE)))
  # print(t1)
  # print(mynames)
  print(is.integer(t1))
  print(class(t1))
  # print(dim(t1))
  if(xvec2 & !xvec1) colnames(t1) <- mynames
  if(print) print(data.frame(.my.prettyNum(t1)))#t1 <- data.frame(formatC(t1, digits = digits, ...)); print(t1)
  if(print){
   t2 <- data.frame(.my.prettyNum(t1))
   colnames(t2) <- mynames
   print(t2)
  } 
 }
 tymch <- t(t1)
 class(tymch) <- "sfgtrehet"
 return(tymch)
}

.rownames.Intercept.change <- function(x){
 x <- gsub("(Intercept)","intercept",x, fixed = TRUE)
 x
}

.teRad <- function(Y,X,M,N,K,
                   Yr,Xr,Kref,
                   rts,base,ifqh,
                   print.level=0){
 .C("radial",
    as.double(Y),
    as.double(X),
    as.integer(M),
    as.integer(N),
    as.integer(K),
    as.double(Yr),
    as.double(Xr),
    as.integer(Kref),
    as.integer(rts),
    as.integer(base),
    as.integer(ifqh),
    as.integer(print.level),
    te = double(K) )$te
}

.teRad1 <- function(Y,X,M,N,K,
                   Yr,Xr,Kref,
                   rts,base,ifqh,
                   print.level=0){
 t1 <- .C("radial",
    Y = as.double(Y),
    X = as.double(X),
    as.integer(M),
    as.integer(N),
    as.integer(K),
    Yr = as.double(Yr),
    Xr = as.double(Xr),
    as.integer(Kref),
    as.integer(rts),
    as.integer(base),
    as.integer(ifqh),
    as.integer(print.level),
    te = double(K) )
 return(list(Y = t1$Y,X = t1$X,Yr = t1$Yr,Xr = t1$Xr))
}

.teNonrad <- function(Y,X,M,N,K,
                      Yr,Xr,Kref,
                      rts,base,ifqh,
                      print.level=0){
 .C("nonradial",
    as.double(Y),
    as.double(X),
    as.integer(M),
    as.integer(N),
    as.integer(K),
    as.double(Yr),
    as.double(Xr),
    as.integer(Kref),
    as.integer(rts),
    as.integer(base),
    as.integer(ifqh),
    as.integer(print.level),
    te = double(K) )$te
}

.dots <- function(nrep, message = NULL, width = 50, character = "."){
 if (!is.numeric(width)) {
  stop("'width' should be numeric")
 }
 if (width != 50 & width != 60 & width != 70 & width != 80 & width != 90  & width != 100){
  stop("'width' should be 50, 60, 70, 80, 90, or 100")
 }
 if (nrep == 0){
  if (!is.null(message)){
   cat("",message,"\n", sep = "")
  }
  cat("____|___ 1 ___|___ 2 ___|___ 3 ___|___ 4 ___|___ 5", sep = "")
  if (width == 60) {
   cat(" ___|___ 6", sep = "")
  }
  if (width == 70) {
   cat(" ___|___ 6 ___|___ 7", sep = "")
  }
  if (width == 80) {
   cat(" ___|___ 6 ___|___ 7 ___|___ 8", sep = "")
  }
  if (width == 90) {
   cat(" ___|___ 6 ___|___ 7 ___|___ 8 ___|___ 9", sep = "")
  }
  if (width == 100) {
   cat(" ___|___ 6 ___|___ 7 ___|___ 8 ___|___ 9 ___|___ 10", sep = "")
  }
  cat("\n")
 }
 else {
  if (nrep/width != floor(nrep/width)){
   cat("",character,"", sep = "")
  }
  else {
   cat("",character," ",nrep,"\n", sep = "")
  }
 }
}

.smplHom <- function(teRef,terfl,Kr,mybw,scVarHom,YorX,ba){
 # sample with replacement from Farrell measures
 # (including reflected ones)
 newSample <- sample( seq_len(2*Kr), Kr, replace = TRUE)
 bStar     <- terfl[newSample]
 # generate the sequence randomly
 epsStar   <- rnorm(Kr)
 tStar     <- bStar + mybw * epsStar
 tStar     <- ifelse(tStar >= 1, tStar, 2-tStar)
 # correct the sequence
 tStar     <- scVarHom * (tStar - mean(bStar)) + mean(bStar)
 # get new xobs if input-based, new yobs if output-based
 if (ba == 1){
  MatStar  <- YorX * ( teRef * tStar )
 }
 else {
  MatStar  <- YorX * ( teRef / tStar )
 }
 return(MatStar)
}

.smplHomTE <- function(terfl,Kr,mybw,scVarHom,ba){
 # sample with replacement from Farrell measures
 # (including reflected ones)
 newSample <- sample( seq_len(2*Kr), Kr, replace = TRUE)
 bStar     <- terfl[newSample]
 # generate the sequence randomly
 epsStar   <- rnorm(Kr)
 tStar     <- bStar + mybw * epsStar
 tStar     <- ifelse(tStar >= 1, tStar, 2-tStar)
 # correct the sequence
 tStar     <- scVarHom * (tStar - mean(bStar)) + mean(bStar)
 # inverse if input-based
 if (ba == 1){
  tStar    <- 1/tStar
 }
 return(tStar)
}

.smplHet1 <- function(Zt,L1,L2,M1,cons1,cons2,onesN,Kr,nZ,ba,M){
 ZStar1 <- NULL
 # print(class(Zt))
 if (ba == 1){
  ymax <- apply(Zt[, 1:M, drop = FALSE], MARGIN = 2, FUN = max)
  # print(ymax)
 } 	else {
  xmin <- apply(Zt[, M:(nZ-1), drop = FALSE], MARGIN = 2, FUN = min)
  # print(xmin)
 }
 cons3 <- cons1[1,] * cons2[1,]
 for(ww in seq_len(Kr)){
  flag <- TRUE
  while( flag ){
   # step 5
   newSample <- sample( seq_len(2*Kr), 1, replace = TRUE)
   ZStar <- Zt[newSample,]
   # step 6
   epsStar <- rnorm(nZ)
   if(newSample <= Kr){
    epsStar <- epsStar %*% L1
   } else {
    epsStar <- epsStar %*% L2
   }
   # step 7
   ZStar <- ZStar + as.vector(epsStar) * cons3
   if(ba == 1){
    flag0 <- max( ZStar[1:M] - ymax ) > 0
   }
   else {
    flag0 <- min( ZStar[M:(nZ-1)] - xmin ) < 0
   }
   flag <- min(ZStar) < 0 || flag0
  }
  ZStar1 <- rbind(ZStar1, ZStar)
 }

 ZStar <- ZStar1

 #  # toss possible negative values
 #  nonNegV <- apply(ZStar, MARGIN = 1, FUN = function(x) ifelse(min(x) < 0, FALSE, TRUE))
 #  # nonNegV <- rep(TRUE, Kr)
 #  ZStar <- ZStar[nonNegV,]
 nonNegV <- rep(TRUE, Kr)

 # step 8
 ZStar[,nZ] <- ifelse(ZStar[,nZ] > 1, ZStar[,nZ], 2 - ZStar[,nZ])

 # step 9
 if (ba == 1){
  teb <- 1 / ZStar[,nZ]
  yrb <- ZStar[, 1:M]
  xrb <- cbind( 1, tan(ZStar[,(M+1):(nZ-1)]) )
 } 	else {
  teb <- ZStar[,nZ]
  yrb <- cbind( 1, tan(ZStar[,1:(M-1)]) )
  xrb <- ZStar[, M:(nZ-1)]
 }
 return( list(Yrb = yrb, Xrb = xrb, teb = teb, Krb = sum(nonNegV)) )
}

.smplHet <- function(Zt,L1,L2,M1,cons1,cons2,onesN,Kr,nZ,ba,M,N){
 MinZ <- ifelse(M == 1, M, M - 1)
 if (ba == 1){
  ymax <- apply(Zt[, 1:M, drop = FALSE], MARGIN = 2, FUN = max)
 } 	else {
  xmin <- apply(Zt[, (MinZ+1):(nZ-1), drop = FALSE], MARGIN = 2, FUN = min)
 }
 # first try
 # step 5
 newSample <- sample( seq_len(2*Kr), Kr, replace = TRUE)
 ZStar <- Zt[newSample,]
 zBarStar <- matrix( colMeans( ZStar ), nrow  = 1)

 # step 6
 epsStar <- matrix(rnorm(Kr * nZ), nrow = Kr, ncol = nZ)
 for(ww in seq_len(Kr)){
  if(newSample[ww] <= Kr){
   epsStar[ww,] <- epsStar[ww, , drop = F] %*% L1
  } else {
   epsStar[ww,] <- epsStar[ww, , drop = F] %*% L2
  }
 }
 # step 7
 ZStar <- kronecker (onesN, zBarStar) + cons1 * ( M1 %*% ZStar + cons2 * epsStar)

 # toss values outside the frontier (including possible negative values)
 nonNeg <- apply(ZStar, MARGIN = 1, FUN = function(x) ifelse(min(x) < 0, FALSE, TRUE))
 if(ba == 1){
  withinFr <- apply(ZStar[,1:M, drop = FALSE], MARGIN = 1, FUN = function(z) min(ymax-z) > 0)
 }
 else {
  withinFr <- apply(ZStar[,(MinZ+1):(nZ-1), drop = FALSE], MARGIN = 1, FUN = function(z) max(xmin-z) < 0)
 }
 Zgood <- withinFr & nonNeg
 ZStar0 <- ZStar[Zgood,]

 # complete Zstar if not full
 if( nrow(ZStar0) < Kr){
  while (nrow(ZStar0) < Kr){
   # step 5
   newSample <- sample( seq_len(2*Kr), Kr, replace = TRUE)
   ZStar <- Zt[newSample,]
   zBarStar <- matrix( colMeans( ZStar ), nrow  = 1)
   # step 6
   epsStar <- matrix(rnorm(Kr * nZ), nrow = Kr, ncol = nZ)
   for(ww in seq_len(Kr)){
    if(newSample[ww] <= Kr){
     epsStar[ww,] <- epsStar[ww, , drop = F] %*% L1
    } else {
     epsStar[ww,] <- epsStar[ww, , drop = F] %*% L2
    }
   }
   # step 7
   ZStar <- kronecker (onesN, zBarStar) + cons1 * ( M1 %*% ZStar + cons2 * epsStar)
   # toss values outside the frontier (including possible negative values)
   nonNeg <- apply(ZStar, MARGIN = 1, FUN = function(x) ifelse(min(x) < 0, FALSE, TRUE))
   if(ba == 1){
    withinFr <- apply(ZStar[,1:M, drop = FALSE], MARGIN = 1, FUN = function(z) min(ymax-z) > 0)
   }
   else {
    withinFr <- apply(ZStar[,(MinZ+1):(nZ-1), drop = FALSE], MARGIN = 1, FUN = function(z) max(xmin-z) < 0)
   }
   Zgood <- withinFr & nonNeg
   ZStar0 <- rbind(ZStar0, ZStar[Zgood,])
  }
 }
 ZStar <- ZStar0[seq_len(Kr),]
 # print(nrow(ZStar))

 # step 8
 ZStar[,nZ] <- ifelse(ZStar[,nZ] > 1, ZStar[,nZ], 2 - ZStar[,nZ])

 # step 9
 if (ba == 1){
  teb <- 1 / ZStar[,nZ]
  if ( N == 1){
   xrb <- ZStar[, (M+1), drop = FALSE]
  }
  else {
   xrb <- cbind( 1, tan(ZStar[,(M+1):(nZ-1)]) )
  }
  yrb <- ZStar[, 1:M, drop = FALSE]
 }
 else {
  teb <- ZStar[,nZ]
  if ( M == 1){
   yrb <- ZStar[, 1, drop = FALSE]
  }
  else {
   yrb <- cbind( 1, tan(ZStar[,1:MinZ]) )
  }
  xrb <- ZStar[, (MinZ+1):(nZ-1), drop = FALSE]
 }
 return( list(Yrb = yrb, Xrb = xrb, teb = teb, Krb = nrow(ZStar)) )
}

.biasAndCI <- function(te, teboot, msub  = 0, K, Kr, M, N,
                       level, smoothed, forceLargerOne){

 if ( forceLargerOne ){
  teff   <- 1/te
  teboot <- 1/teboot
 }
 else {
  teff <- te
 }

 if (smoothed) {
  con1 <- 1
  con2 <- 1
 }
 else {
  con1 <- msub ^ ( 2 / (M + N + 1) )
  con2 <- Kr ^ (-2 / (M + N + 1) )
 }
 con3 <- con1 * con2

 # 	# bias-correction
 # 	tebm <- colMeans(teboot, na.rm = TRUE)
 # 	bias <- con3 * (tebm - teff)
 # 	tebc <- teff - bias
 # 	vari <- apply(teboot, 2, var, na.rm = TRUE)
 # 	BovV <- bias^2 / vari * 3
 #
 # 	# CI
 # 	teOverTEhat <- sweep(teboot, 2, FUN = "/", teff)
 # 	teOverTEhat <- (teOverTEhat - 1) * con1
 # 	quans <- apply(teOverTEhat, 2, quantile, probs = (100 + c(-level, level))/200, na.rm = T)
 # 	LL <- teff / (1 + quans[2] * con2)
 # 	UU <- teff / (1 + quans[1] * con2)
 #
 # 	# drop if inference is based on less than 100
 # 	reaB <- apply(teboot, 2, function(x) sum( !is.na(x) ) )
 #
 # 	bias <- ifelse(reaB < 100, NA, bias)
 # 	vari <- ifelse(reaB < 100, NA, vari)
 # 	tebc <- ifelse(reaB < 100, NA, tebc)
 # 	LL <- ifelse(reaB < 100, NA, LL)
 # 	UU <- ifelse(reaB < 100, NA, UU)

 bias <- vari <- reaB <- BovV <- tebc <- LL <- UU <- numeric(K)

 for (i in seq_len(K) ) {
  TEi <- teboot[,i]
  TEi <- TEi[!is.na(TEi)]
  reaB[i] <- length(TEi)
  if( reaB[i] < 100 ){
   bias[i] = NA
   vari[i] = NA
   BovV[i] = NA
   tebc[i] = NA
   LL[i]   = NA
   UU[i]   = NA
  }
  else {
   bias[i] = con3 * ( mean(TEi) - teff[i] )
   vari[i] = var(TEi)
   BovV[i] = (bias[i])^2 / vari[i] * 3
   tebc[i] = teff[i] - bias[i]
   teOverTEhat = (TEi / teff[i] - 1) * con1
   quans = quantile(teOverTEhat, probs = (100 + c(-level, level))/200, na.rm = TRUE)
   LL[i] = teff[i] / (1 + quans[2] * con2)
   UU[i] = teff[i] / (1 + quans[1] * con2)
  }
 }

 return(list(reaB=reaB, bias=bias, vari=vari, BovV=BovV, tebc=tebc, LL=LL, UU=UU))

}


.pvalsTestOne <- function(seCrs, seCrsMean, te1boot, te2boot, ba){
 seCrsB <- te1boot / te2boot
 seCrsMeanB <- rowMeans(te1boot, na.rm = TRUE) / rowMeans(te2boot, na.rm = TRUE)
 seCrsMeanB <- na.omit( seCrsMeanB )
 if(ba ==1){
  # global CRS
  pgCRS <- mean(seCrsMeanB <= seCrsMean)
  # local CRS
  plCRS <- rowMeans( apply(seCrsB, MARGIN = 1, FUN = function(z) z <= seCrs), na.rm = TRUE )
 }
 else {
  pgCRS <- mean(seCrsMeanB >= seCrsMean)
  # local CRS
  plCRS <- rowMeans( apply(seCrsB, MARGIN = 1, FUN = function(z) z >= seCrs), na.rm = TRUE )
 }

 # count those that are not NA
 nonNa <- apply(seCrsB, MARGIN = 2, FUN = function(z) length(na.omit(z)))
 plCRS <- ifelse(nonNa >= 100, plCRS, NA)

 return(list(pgCRS = pgCRS, plCRS = plCRS, nonNa = nonNa, reaB = length(seCrsMeanB)))
}

.pvalsTestTwo <- function(seNrs, seNrsMean, te1boot, te2boot, ba, performGlobal, s.inefficient){
 seNrsB <- te1boot / te2boot
 # cat("ncol(te1boot) = ",ncol(te1boot),"\n")
 # if(ncol(te1boot)==1) seNrsB <- matrix(seNrsB, ncol = 1)
 if ( performGlobal ){
  seNrsMeanB <- rowMeans(te1boot, na.rm = TRUE) / rowMeans(te2boot, na.rm = TRUE)
  seNrsMeanB <- na.omit( seNrsMeanB )

  if(ba ==1){
   # global NRS
   pgNRS <- mean(seNrsMeanB <= seNrsMean)
   # local NRS
   plNRS <- rowMeans( apply(seNrsB, MARGIN = 1, FUN = function(z) z <= seNrs), na.rm = TRUE )
  }
  else {
   # global NRS
   pgNRS <- mean(seNrsMeanB >= seNrsMean)
   # local NRS
   plNRS <- rowMeans( apply(seNrsB, MARGIN = 1, FUN = function(z) z >= seNrs), na.rm = TRUE )
  }
 }
 else {
  seNrs1 <- seNrs[s.inefficient]
  if(ba ==1){
   # local NRS
   if(ncol(te1boot) == 1){
    plNRS <- mean( apply(seNrsB, MARGIN = 1, FUN = function(z) z <= seNrs1), na.rm = TRUE )
   }
   else {
    plNRS <- rowMeans( apply(seNrsB, MARGIN = 1, FUN = function(z) z <= seNrs1), na.rm = TRUE )
   }
  }
  else {
   # local NRS
   if(ncol(te1boot) == 1){
    plNRS <- mean( apply(seNrsB, MARGIN = 1, FUN = function(z) z >= seNrs1), na.rm = TRUE )
   }
   else {
    plNRS <- rowMeans( apply(seNrsB, MARGIN = 1, FUN = function(z) z >= seNrs1), na.rm = TRUE )
   }
  }
 }
 # count those that are not NA
 nonNa <- apply(seNrsB, MARGIN = 2, FUN = function(z) length(na.omit(z)))
 plNRS <- ifelse(nonNa >= 100, plNRS, NA)

 # write into appropriate cells
 pineffdrs <- rep(NA, length(s.inefficient))

 if (performGlobal){
  pineffdrs <- ifelse(s.inefficient, plNRS, pineffdrs)
 }
 else {
  pineffdrs[s.inefficient] <- plNRS
 }
 if ( performGlobal ){
  tymch <- list(pgNRS = ifelse(performGlobal, pgNRS, NA), pineffdrs = pineffdrs, nonNa = nonNa, reaB = length(seNrsMeanB))
 }
 else {
  tymch <- list(pgNRS = ifelse(performGlobal, pgNRS, NA), pineffdrs = pineffdrs, nonNa = nonNa)
 }
 return(tymch)
}

# empirical distribution function
.edf <- function(Y,y){
 if(length(y) == 1){
  f1 <- Y<=y
  # f1 <- y-Y > 1e-6
 }
 else {
  f3 <- cbind(y, t(Y))
  f2 <- apply(f3, MARGIN = 1, FUN = function(z) z[-1] <= z[1])
  # f2 <- apply(f3, MARGIN = 1, FUN = function(z)  z[1]-z[-1]> 1e-6)
  f1 <- rowMeans(f2)==1
 }
 return(mean(f1))
}

# empirical joint distribution function
.ejdf <- function(Y, X, y, x){
 # 1
 if(length(y) == 1){
  f1 <- Y<=y
 }
 else {
  f3 <- cbind(y, t(Y))
  f2 <- apply(f3, MARGIN = 1, FUN = function(z) z[-1] <= z[1])
  f1 <- rowMeans(f2)==1
 }
 # 2
 if(length(x) == 1){
  g1 <- X<=x
 }
 else {
  g3 <- cbind(x, t(X))
  g2 <- apply(g3, MARGIN = 1, FUN = function(z) z[-1] <= z[1])
  g1 <- rowMeans(g2)==1
 }
 sum(f1 & g1) / length(f1)
}

# empirical joint distribution function
.ejdfedf1 <- function(Y, X, y, x){
 # 1
 if(length(y) == 1){
  f1 <- Y<=y
 }
 else {
  f3 <- cbind(y, t(Y))
  f2 <- apply(f3, MARGIN = 1, FUN = function(z) z[-1] <= z[1])
  f1 <- rowMeans(f2)==1
 }
 # 2
 if(length(x) == 1){
  g1 <- X<=x
 }
 else {
  g3 <- cbind(x, t(X))
  g2 <- apply(g3, MARGIN = 1, FUN = function(z) z[-1] <= z[1])
  g1 <- rowMeans(g2)==1
 }
 # mean(f1 & g1) - mean(f1)*mean(g1)
 mean(f1)
}

.ejdfedf <- function(Y, X, y, x){
 # 1
 f1 <- Y <= y
 # 2
 g2 <- X
 for(i in seq_len(ncol(X))){
  g2[,i] <- X[,i] <= x[i]
 }
 g1 <- ( rowMeans(g2) ==1 )
 mean(f1 & g1) - mean(f1)*mean(g1)
 # mean(f1)
}

.t4n <- function(w,d, print=FALSE){
 n <- length(d)
 tt <- numeric(n)
 for(i in 1:n){
  tt[i] <- .ejdfedf(Y = d, y = d[i], X = w, x = w[i,]) # - .edf(Y = d, y = d[i]) * .edf(Y = w, y = w[i,]) #
 }
 if(print) print(cbind(tt,d))
 return(sum(tt^2))
}

# begin parallel computing
.run.boots.hom.rts <- function(inp, ...){
 # Read the arguments
 args <- list(...)[[1]]
 Y <- args$Y
 X <- args$X
 ba <- args$ba
 teRef <- args$teRef
 terfl <- args$terfl
 mybw <- args$mybw
 scVarHom <- args$scVarHom

 M  <- ncol(Y)
 N  <- ncol(X)
 K  <- nrow(Y)
 # Prepare output structure
 te1boot <- rep(NA, K)
 te2boot <- rep(NA, K)

 for(ii in seq_len(K)){
  # begin homogeneous
  # print(1)
  if (ba == 2) {
   # step 1: sampling
   Yb <- .smplHom(teRef,terfl,Kr=K,mybw,scVarHom,Y,ba)
   # step 2: applying DEA
   # CRS or NiRS
   te1boot[ii] <- .teRad(Y[ii,],X[ii,],M,N,1,t(Yb),t(X),K,args$rtsHo,ba,
                         0,print.level=0)
   # VRS
   te2boot[ii] <- .teRad(Y[ii,],X[ii,],M,N,1,t(Yb),t(X),K,1,ba,
                         0,print.level=0)
  }
  else {
   # step 1: sampling
   Xb <- .smplHom(teRef,terfl,Kr=K,mybw,scVarHom,X,ba)
   # step 2: applying DEA
   # CRS or NiRS
   te1boot[ii] <- .teRad(Y[ii,],X[ii,],M,N,1,t(Y),t(Xb),K,args$rtsHo,ba,
                         0,print.level=0)
   # VRS
   te2boot[ii] <- .teRad(Y[ii,],X[ii,],M,N,1,t(Y),t(Xb),K,1,ba,
                         0,print.level=0)
  }
  # end homogeneous
  # end for ii
 }
 return (c(te1boot, te2boot, K, 0))
}

.run.boots.hom.bc <- function(inp, ...){
 # Read the arguments
 args <- list(...)[[1]]
 Y <- args$Y
 X <- args$X
 Yr <- args$Yr
 Xr <- args$Xr
 rt <- args$rt
 ba <- args$ba
 teRef <- args$teRef
 terfl <- args$terfl
 mybw <- args$mybw
 scVarHom <- args$scVarHom

 M  <- ncol(Y)
 N  <- ncol(X)
 K  <- nrow(Y)
 Kr <- nrow(Yr)

 # begin homogeneous
 # print(1)
 # teB <- numeric(K)
 if (ba == 2) {
  # step 1: sampling
  Yrb <- .smplHom(teRef,terfl,Kr,mybw,scVarHom,Yr,ba)
  # step 2: applying DEA
  teB <- .teRad(t(Y),t(X),M,N,K,t(Yrb),t(Xr),Kr,rt,ba,
                0,print.level=0)
 }
 else {
  # step 1: sampling
  Xrb <- .smplHom(teRef,terfl,Kr,mybw,scVarHom,Xr,ba)
  # step 2: applying DEA
  teB <- .teRad(t(Y),t(X),M,N,K,t(Yr),t(Xrb),Kr,rt,ba,
                0,print.level=0)
 }
 return (c(teB, K, 1))
}

.run.boots.het.rts <- function(inp, ...){
 # Read the arguments
 args <- list(...)[[1]]
 Y <- args$Y
 X <- args$X
 Yr <- args$Yr
 Xr <- args$Xr
 Zt <- args$Zt
 L1 <- args$L1
 L2 <- args$L2
 M1 <- args$M1
 cons1 <- args$cons1
 cons2 <- args$cons2
 onesN <- args$onesN
 ba <- args$ba

 M  <- ncol(Y)
 N  <- ncol(X)
 K  <- nrow(Y)
 Kr <- nrow(Yr)
 nZ <- ncol(Zt)
 # Prepare output structure
 te1boot <- rep(NA, K)
 te2boot <- rep(NA, K)

 for(ii in seq_len(K)){
  # begin heterogeneous
  # print(2)
  # step 1: sampling
  smplHet <- .smplHet(Zt,L1,L2,M1,cons1,cons2,onesN,Kr,nZ,ba,M,N)
  # step 2: get efficiency under H0
  teRefB  <- .teRad(t(smplHet$Yrb),t(smplHet$Xrb),M,N,smplHet$Krb,
                    t(Yr),t(Xr),Kr,rts=3,ba,0,print.level=0)
  # step 3: get efficient xRef or yRef
  if (ba == 1) {
   Yrb <- smplHet$Yrb
   Xrb <- smplHet$Xrb * teRefB / smplHet$teb
  }
  else {
   Yrb <- smplHet$Yrb * teRefB / smplHet$teb
   Xrb <- smplHet$Xrb
  }
  # handle infeasible
  teRefB.good <- !is.na(teRefB)
  # modify the existing matrices by tossing the missing values
  Yrb <- Yrb[teRefB.good, , drop = FALSE]
  Xrb <- Xrb[teRefB.good, , drop = FALSE]
  # new number of obs in reference
  KrB <- sum(teRefB.good)
  cat("my sample is ",KrB,"\n", sep = "")
  # step 4: get the bootstrapped TE
  # CRS or NiRS
  te1boot[ii] <- .teRad(Y[ii,],X[ii,],M,N,1,t(Yrb),t(Xrb),KrB,args$rtsHo,ba,
                        0,print.level=0)
  # VRS
  te2boot[ii] <- .teRad(Y[ii,],X[ii,],M,N,1,t(Yrb),t(Xrb),KrB,1,ba,
                        0,print.level=0)
  # end heterogeneous
  # end for ii
 }
 return (c(te1boot, te2boot, K, 0))
}

.run.boots.het.bc <- function(inp, ...){
 # Read the arguments
 args <- list(...)[[1]]
 Y <- args$Y
 X <- args$X
 Yr <- args$Yr
 Xr <- args$Xr
 rt <- args$rt
 Zt <- args$Zt
 L1 <- args$L1
 L2 <- args$L2
 M1 <- args$M1
 cons1 <- args$cons1
 cons2 <- args$cons2
 onesN <- args$onesN
 ba <- args$ba

 M  <- ncol(Y)
 N  <- ncol(X)
 K  <- nrow(Y)
 Kr <- nrow(Yr)
 nZ <- ncol(Zt)

 # begin heterogeneous
 # print(2)
 # step 1: sampling
 smplHet <- .smplHet(Zt,L1,L2,M1,cons1,cons2,onesN,Kr,nZ,ba,M,N)
 # step 2: get efficiency under H0
 teRefB  <- .teRad(t(smplHet$Yrb),t(smplHet$Xrb),M,N,smplHet$Krb,
                   t(Yr),t(Xr),Kr,rt,ba, 0,print.level=0)
 # step 3: get efficient xRef or yRef
 if (ba == 1) {
  Yrb <- smplHet$Yrb
  Xrb <- smplHet$Xrb * teRefB / smplHet$teb
 }
 else {
  Yrb <- smplHet$Yrb * teRefB / smplHet$teb
  Xrb <- smplHet$Xrb
 }
 # handle infeasible
 teRefB.good <- !is.na(teRefB)
 # modify the existing matrices by tossing the missing values
 Yrb <- Yrb[teRefB.good, , drop = FALSE]
 Xrb <- Xrb[teRefB.good, , drop = FALSE]
 # new number of obs in reference
 KrB <- sum(teRefB.good)
 # step 4: get the bootstrapped TE
 teB <- .teRad(t(Y),t(X),M,N,K,t(Yrb),t(Xrb),KrB,rt,ba,
               0,print.level=0)
 # end heterogeneous
 return (c(teB, KrB, 1))
}

.run.dots <- function(xs, nrep, args){
 width <- args$width
 character <- "."
 #  mylen <- length(xs[[nrep]])
 #  if (mylen > 0){
 #   if(xs[[nrep]][ mylen ] == 1){
 #    # boot for BC
 #    K <- mylen - 2
 #   }
 #   else {
 #    # boot for RTS test
 #    K <- (mylen - 2 )/2
 #   }
 #   # Pre-Last value in each output is 'KrB'
 #   KrB <- xs[[nrep]][ mylen-1 ]
 #   # cat("krB = ",KrB," K = ", K, " \n", sep = "")
 #   if (KrB/K < 0.80){
 #    character <- "?"
 #   }
 #   else if (KrB/K < 0.95){
 #    character <- "@"
 #   }
 #   else if (KrB/K < 1){
 #    character <- "x"
 #   }
 #  }

 #  if (!is.numeric(width)) {
 #   stop("'width' should be numeric")
 #  }
 #  if (width != 50 & width != 60 & width != 70 &
 #      width != 80 & width != 90  & width != 100){
 #   stop("'width' should be 50, 60, 70, 80, 90, or 100")
 #  }
 #  if (nrep == 0){
 #   if (!is.null(message)){
 #    cat("",message,"\n", sep = "")
 #   }
 #   cat("____|___ 1 ___|___ 2 ___|___ 3 ___|___ 4 ___|___ 5", sep = "")
 #   if (width == 60) {
 #    cat(" ___|___ 6", sep = "")
 #   }
 #   if (width == 70) {
 #    cat(" ___|___ 6 ___|___ 7", sep = "")
 #   }
 #   if (width == 80) {
 #    cat(" ___|___ 6 ___|___ 7 ___|___ 8", sep = "")
 #   }
 #   if (width == 90) {
 #    cat(" ___|___ 6 ___|___ 7 ___|___ 8 ___|___ 9", sep = "")
 #   }
 #   if (width == 100) {
 #    cat(" ___|___ 6 ___|___ 7 ___|___ 8 ___|___ 9 ___|___ 10", sep = "")
 #   }
 #   cat("\n")
 #  }
 #  else {
 if (nrep/width != floor(nrep/width)){
  cat("",character,"", sep = "")
 }
 else {
  cat("",character," ",nrep,"\n", sep = "")
 }
 # }
}

.run.boots.subs.bc <- function(inp, ...){
 # Read the arguments
 args <- list(...)[[1]]
 Y <- args$Y
 X <- args$X
 Yr <- args$Yr
 Xr <- args$Xr
 rt <- args$rt
 ba <- args$ba
 msub <- args$msub

 M  <- ncol(Y)
 N  <- ncol(X)
 K  <- nrow(Y)
 Kr <- nrow(Yr)

 # begin subsampling
 # print(3)
 # step 1: sub-sampling
 newSample <- sample( seq_len(Kr), msub, replace = TRUE)
 ystar <- Yr[newSample, , drop = FALSE]
 xstar <- Xr[newSample, , drop = FALSE]
 # step 2: applying DEA
 teB <- .teRad(t(Y),t(X),M,N,K,t(ystar),t(xstar),msub,rt,ba,
               0,print.level=0)
 mychar <- "."
 # end subsampling
 return (c(teB, K, 1))
}

# end parallel computing

.prepareYXZ.cs <- function(formula, ln.var.u.0i = NULL, ln.var.v.0i = NULL, mean.u.0i = NULL, data, subset, sysnframe = 1, ...) {
 # needed.frame <- sys.nframe() - 1
 needed.frame <- sys.nframe() - sysnframe
 # cat("sys.nframe() = ",sys.nframe(),"\n")
 # cat("needed.frame = ",needed.frame,"\n")
 mf0 <- match.call(expand.dots = FALSE, call = sys.call(sys.parent(n = needed.frame)))
 if ( is.null(ln.var.u.0i) ){
  ln.var.u.0i <- ~ 1
  ku <- 1
 } else {
  ku <- 17
 }
 if ( is.null(ln.var.v.0i) ){
  ln.var.v.0i <- ~ 1
  kv <- 1
 } else {
  kv <- 17
 }
 if ( is.null(mean.u.0i) ){
  mean.u.0i <- ~ 1
  kdel <- 1
 } else {
  kdel <- 17
 }
 
 form1 <- Formula(as.formula(paste("",deparse(formula, width.cutoff = 500L)," | ",ln.var.u.0i[2]," | ",ln.var.v.0i[2]," | ",mean.u.0i[2],"", sep = "")))
 
 form <- Formula(as.formula(paste("",deparse(formula, width.cutoff = 500L)," + ",ln.var.u.0i[2]," + ",ln.var.v.0i[2]," + ",mean.u.0i[2],"", sep = "")))
 
 # check if it is a matrix
 datasupplied <- !(match("data", names(mf0), 0) == 0)
 
 
 if(datasupplied){
  # begin get a logical vector equal TRUE if !missing
  # first using data and subset to get x without NA
  mf <- mf0
  mf$formula <- formula( form )
  m <- match(c("formula", "data", "subset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, sys.frame(sys.parent(n = needed.frame)))
  mt <- attr(mf, "terms")
  X <- as.matrix(model.matrix(mt, mf))
  # now get the names in the entire data
  # esample <- seq_len( nrow(data) ) %in% as.numeric(rownames(X))
  esample <- rownames(data) %in% rownames(X)
  if(sum(esample)==0){
   stop("no valid data points")
  }
  # print(table(esample))
  # end get a logical vector equal TRUE if !missing
  
  # get the data
  # print(form1)
  dataesample <- model.frame(form1, data = data[esample,])
  # print(dataesample)
  Y <- as.matrix(model.part(form1, data = dataesample, lhs = 1, drop = FALSE))
  # Y <- as.matrix( model.matrix(formula(form1, lhs = 1, rhs = 0), data = data[esample,]))
  # print(Y)
  X <- as.matrix( model.matrix(formula(form1, lhs = 0, rhs = 1), data = dataesample))
  # print(X)
  n <- nrow(Y)
  k <- ncol(X)
  if(ku == 1){
   Zu <- matrix(1, nrow = sum(esample), ncol = 1)
   colnames(Zu) <- "(Intercept)"
  }
  else {
   Zu <- as.matrix( model.matrix(formula(form1, lhs = 0, rhs = 2), data = dataesample))
   ku <- ncol(Zu)
  }
  if(kv == 1){
   Zv <- matrix(1, nrow = sum(esample), ncol = 1)
   colnames(Zv) <- "(Intercept)"
  }
  else {
   Zv <- as.matrix( model.matrix(formula(form1, lhs = 0, rhs = 3), data = dataesample))
   kv <- ncol(Zv)
  }
  if(kdel == 1){
   Zdel <- matrix(1, nrow = sum(esample), ncol = 1)
   colnames(Zdel) <- "(Intercept)"
  }
  else {
   Zdel <- as.matrix( model.matrix(formula(form1, lhs = 0, rhs = 4), data = dataesample))
   kdel <- ncol(Zdel)
  }
  #   print(head(Y))
  #   print(head(X))
  #   print(head(Zu))
  #   print(head(Zv))
  #   print(head(Zdel))
 }
 # if data are not supplied
 else {
  # begin get a logical vector equal TRUE if !missing
  
  # first using data and subset to get XZ without NA
  mf <- mf0
  mf$formula <- formula( form )
  subsetsupplied <- !(match("subset", names(mf), 0) == 0)
  if(subsetsupplied) stop("Subset with matrices is not allowed")
  m <- match(c("formula"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- as.name("model.frame")
  # mf <- eval(mf, parent.frame())
  mf <- eval(mf, sys.frame(sys.parent(n = needed.frame)))
  mt <- attr(mf, "terms")
  Y <- as.matrix(model.response(mf))
  n <- nrow(Y)
  XZ <- as.matrix(model.matrix(mt, mf))
  # get a logical vector equal TRUE if !missing
  with.na <- model.frame(mt, na.action = na.pass)
  esample <- rownames(with.na) %in% rownames(XZ)
  if(sum(esample)==0){
   stop("no valid data points")
  }
  # print(table(esample))
  # end get a logical vector equal TRUE if !missing
  
  # get the data
  #   print(head(Y))
  #   print(head(XZ))
  
  # get X
  mf <- mf0
  m <- match(c("formula"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, sys.frame(sys.parent(n = needed.frame)))
  mt <- attr(mf, "terms")
  X <- as.matrix(model.matrix(mt, mf))
  # print(head(X))
  k <- ncol(X)
  
  # get Zu
  if(ku == 1){
   Zu <- matrix(1, nrow = sum(esample), ncol = 1)
   colnames(Zu) <- "(Intercept)"
  }
  else {
   mf <- mf0
   mf$formula <- formula( Formula(as.formula(paste("",deparse(formula)," + ",ln.var.u.0i[2],"", sep = ""))))
   m <- match(c("formula"), names(mf), 0L)
   mf <- mf[c(1L, m)]
   mf[[1L]] <- as.name("model.frame")
   mf <- eval(mf, sys.frame(sys.parent(n = needed.frame)))
   mt <- attr(mf, "terms")
   XZu <- as.matrix(model.matrix(mt, mf))
   Zu <- XZu[, -(2:k), drop = FALSE]
   ku <- ncol(Zu)
   # print(head(Zu))
  }
  
  # get Zv
  if(kv == 1){
   Zv <- matrix(1, nrow = sum(esample), ncol = 1)
   colnames(Zv) <- "(Intercept)"
  }
  else {
   mf <- mf0
   mf$formula <- formula( Formula(as.formula(paste("",deparse(formula)," + ",ln.var.u.0i[2]," + ",ln.var.v.0i[2],"", sep = ""))))
   m <- match(c("formula"), names(mf), 0L)
   mf <- mf[c(1L, m)]
   mf[[1L]] <- as.name("model.frame")
   mf <- eval(mf, sys.frame(sys.parent(n = needed.frame)))
   mt <- attr(mf, "terms")
   XZuZv <- as.matrix(model.matrix(mt, mf))
   Zv <- XZuZv[, -(2:(k+ku-1)), drop = FALSE]
   kv <- ncol(Zv)
   # print(head(Zv))
  }
  
  # get Zdel
  if(kdel == 1){
   Zdel <- matrix(1, nrow = sum(esample), ncol = 1)
   colnames(Zdel) <- "(Intercept)"
  }
  else {
   Zdel <- XZ[, -(2:(k+ku-1+kv-1)), drop = FALSE]
   # print(head(Zdel))
  }
  
  # print(head(Zu))
  # print(head(Zv))
  # print(head(Zdel))
 }
 # colnames(X)[1] <- colnames(Zu)[1] <- colnames(Zv)[1] <- colnames(Zdel)[1] <- "Intercept"
 
 # if(sum(colnames(X) == "(Intercept)") > 0){
 #  # print(which(colnames(X) == "(Intercept)"))
 #  colnames(X)[which(colnames(X) == "(Intercept)")] <- "Intercept"
 # }
 # if(sum(colnames(Zu) == "(Intercept)") > 0){
 #  # print(which(colnames(Zu) == "(Intercept)"))
 #  colnames(Zu)[which(colnames(Zu) == "(Intercept)")] <- "Intercept"
 # }
 # if(sum(colnames(Zv) == "(Intercept)") > 0){
 #  # print(which(colnames(Zv) == "(Intercept)"))
 #  colnames(Zv)[which(colnames(Zv) == "(Intercept)")] <- "Intercept"
 # }
 # if(sum(colnames(Zdel) == "(Intercept)") > 0){
 #  # print(which(colnames(Zv) == "(Intercept)"))
 #  colnames(Zdel)[which(colnames(Zdel) == "(Intercept)")] <- "Intercept"
 # }
 
 colnames(X) <- .rownames.Intercept.change(colnames(X))
 colnames(Zu) <- .rownames.Intercept.change(colnames(Zu))
 colnames(Zv) <- .rownames.Intercept.change(colnames(Zv))
 colnames(Zdel) <- .rownames.Intercept.change(colnames(Zdel))
 
 # print(colnames(Zdel))
 
 tymch <- list(Y = Y, X = X, Zu = Zu, Zv = Zv, Zdel = Zdel, n = n, k = k, ku = ku, kv = kv, kdel = kdel, esample = esample)
 class(tymch) <- "npsf"
 return(tymch)
 
}

# Half-normal model

# Log-likelihood
.ll.hn <- function(theta, prod, k, kv, ku, kdel = NULL, y, Zv, Zu, Zdel = NULL, X) {
 s <- ifelse(prod, -1, 1)
 beta <- theta[1:k]
 gv <- theta[(k+1):(k+kv)]
 gu <- theta[-c(1:k, (k+1):(k+kv))]
 e <- y-X%*%beta
 sig <- sqrt(exp(Zv%*%gv) + exp(Zu%*%gu))
 lmd <- sqrt(exp(Zu%*%gu)/exp(Zv%*%gv))

 # Log-likelihood
 llf <- sum(log(2) - log(sqrt(2*pi)) - log(sig) + pnorm(s*e*lmd/sig, log.p = TRUE) - 0.5*(e/sig)^2)
 return(llf)
}

# Gradient
.g.hn <- function(theta, prod, k, kv, ku, kdel = NULL, y, Zv, Zu, Zdel = NULL, X) {
 if(prod == TRUE){ s = -1 } else {s = 1}
 beta <- theta[1:k]
 gv <- theta[(k+1):(k+kv)]
 gu <- theta[-c(1:k, (k+1):(k+kv))]
 e <- y - X%*%beta
 expu <- exp(Zu%*%gu); expv <- exp(Zv%*%gv)
 sig <- sqrt(expv + expu); lmd <- sqrt(expu/expv)
 ls <- lmd/sig; g <- dnorm(s*e*ls)/pnorm(s*e*ls)
 gels <- g*e*ls;

 # Gradient
 gb <- t(X)%*%(-s*g*ls + e/sig^2)
 ggv <- 0.5*t(Zv)%*%(expv/sig^2 * ((e/sig)^2 - 1) - s*gels*(1 + expv/sig^2))
 ggu <- 0.5*t(Zu)%*%(expu/sig^2 * ((e/sig)^2 - 1) - s*gels*(expu/sig^2 - 1))
 grad <- rbind(gb, ggv, ggu)
 return(grad)
}

# Hessian
.hess.hn <- function(theta, prod, k, kv, ku, kdel = NULL, y, Zv, Zu, Zdel = NULL, X) {
 s <- ifelse(prod, -1, 1)
 beta <- theta[1:k]
 gv <- theta[(k+1):(k+kv)]
 gu <- theta[-c(1:k, (k+1):(k+kv))]
 e <- y - X%*%beta
 expu <- exp(Zu%*%gu); expv <- exp(Zv%*%gv)
 sig <- sqrt(expv + expu); lmd <- sqrt(expu/expv)
 ls <- lmd/sig; els <- e*ls;
 g <- dnorm(s*e*ls)/pnorm(s*e*ls); gels <- g*els

 # Hessian
 Hb <- t(as.numeric(-s*g*(els + s*g)*ls^2 - sig^(-2))*X)%*%X
 Hbgu <- t(as.numeric(0.5*g*ls*(1 - expu/sig^2)*((g + s*els)*els - s*1) - e*expu/sig^4)*X)%*%Zu
 Hbgv <- t(as.numeric(-s*0.5*g*ls*(1 + expv/sig^2)*((els + s*g)*els - 1) - e*expv/sig^4)*X)%*%Zv
 Hgvu <- t(0.5*as.numeric(-s*0.5*gels*(1 - expu/sig^2)*(1 + expv/sig^2)*(1 - s*(g + s*els)*els) - expv*expu*(2*(e/sig)^2 - s*gels - 1)/sig^4)*Zv)%*%Zu
 Hgu <- t(0.5*as.numeric((1 - expu/sig^2)*(expu/sig^2 * ((e/sig)^2 - s*gels - 1) + s*0.5*gels*(1 - expu/sig^2)*(1 - s*(g + s*els)*els)) - (expu*e/sig^3)^2)*Zu)%*%Zu
 Hgv <- t(0.5*as.numeric(expv/sig^2 * ((1 - expv/sig^2)*((e/sig)^2 - s*gels - 1) - e^2*expv/sig^4) - s*0.5*gels*(1 + expv/sig^2)^2 * ((els + s*g)*els - 1))*Zv)%*%Zv
 H <- cbind(rbind(Hb, t(Hbgv), t(Hbgu)), rbind(Hbgv, Hgv, t(Hgvu)), rbind(Hbgu, Hgvu, Hgu))
 return(H)
}

# Truncated-normal model

#Log-likelihood
.ll.tn <- function(theta, prod, k, kv, ku, kdel, y, Zv, Zu, X, Zdel) {
 s <- ifelse(prod, -1, 1)
 beta <- theta[1:k]
 gv <- theta[(k+1):(k+kv)]
 gu <- theta[(k+kv+1):(k+kv+ku)]
 delta <- theta[(k+kv+ku+1):(k+kv+ku+kdel)]
 e <- y-X%*%beta
 sig <- sqrt(exp(Zv%*%gv) + exp(Zu%*%gu))
 lmd <- sqrt(exp(Zu%*%gu)/exp(Zv%*%gv))
 mu <- Zdel%*%delta

 # Log-likelihood
 llf <- sum(-log(sqrt(2*pi)) - log(sig) - pnorm(mu/sqrt(exp(Zu%*%gu)), log.p = TRUE) + pnorm(mu/(sig*lmd) + s*e*lmd/sig, log.p = TRUE) - 0.5*(((e - s*mu)^2)/sig^2))
 return(llf)
}

#Gradient

.g.tn <- function(theta, prod, k, kv, ku, kdel, y, Zv, Zu, X, Zdel) {
 s <- ifelse(prod, -1, 1)
 beta <- theta[1:k]
 gv <- theta[(k+1):(k+kv)]
 gu <- theta[(k+kv+1):(k+kv+ku)]
 delta <- theta[(k+kv+ku+1):(k+kv+ku+kdel)]
 e <- y-X%*%beta
 expu <- exp(Zu%*%gu); expv <- exp(Zv%*%gv)
 sig <- sqrt(exp(Zv%*%gv) + exp(Zu%*%gu))
 lmd <- sqrt(exp(Zu%*%gu)/exp(Zv%*%gv))
 mu <- Zdel%*%delta; ls <- lmd/sig
 g1 <- dnorm(ls*(mu/lmd^2 + s*e))/pnorm(ls*(mu/lmd^2 + s*e))
 g2 <- dnorm(mu/sqrt(expu))/pnorm(mu/sqrt(expu)); ml = mu/lmd^2
 d1 <- 0.5*mu/(sig*lmd)*(1 - expv/sig^2)
 d2 <- -0.5*s*e*ls * (1 + expv/sig^2)
 d3 = -0.5*mu/(sig*lmd)*(1 + expu/sig^2)
 d4 = 0.5*s*e*ls * (1 - expu/sig^2)

 # Gradient
 gb <- t(X)%*%((e - s*mu)/sig^2 - s*g1*ls)
 gdel <- t(Zdel)%*%(-g2/sqrt(expu) + g1/(sig*lmd) + s*(e - s*mu)/sig^2)
 gv <- t(Zv)%*%(-0.5*expv/sig^2 + g1*(d1 + d2) + 0.5*expv*(e - s*mu)^2/sig^4)
 gu  <- t(Zu)%*%(-0.5*expu/sig^2 + 0.5*g2*mu/sqrt(expu) + g1*(d3 + d4) + 0.5*expu*(e - s*mu)^2/sig^4)
 grad <- rbind(gb, gv, gu, gdel)
 return(grad)
}

#Hessian
.hess.tn <- function(theta, prod, k, kv, ku, kdel, y, Zv, Zu, X, Zdel) {
 s <- ifelse(prod, -1, 1)
 beta <- theta[1:k]
 gv <- theta[(k+1):(k+kv)]
 gu <- theta[(k+kv+1):(k+kv+ku)]
 delta <- theta[(k+kv+ku+1):(k+kv+ku+kdel)]
 e <- y-X%*%beta
 expu <- exp(Zu%*%gu); expv <- exp(Zv%*%gv)
 sig <- sqrt(exp(Zv%*%gv) + exp(Zu%*%gu))
 lmd <- sqrt(exp(Zu%*%gu)/exp(Zv%*%gv))
 mu <- Zdel%*%delta; ls <- lmd/sig
 g1 <- dnorm(ls*(mu/lmd^2 + s*e))/pnorm(ls*(mu/lmd^2 + s*e))
 gr1u <- ls*(0.5*(1 - expu/sig^2)*(mu/lmd^2 + s*e) - mu/lmd^2)*g1*(-s*ls*(e + s*mu/lmd^2) - g1)
 gr1v <- ls*(mu/lmd^2 - 0.5*(1 + expv/sig^2)*(mu/lmd^2 + s*e))*g1*(-s*ls*(e + s*mu/lmd^2) - g1)
 g2 <- dnorm(mu/sqrt(expu))/pnorm(mu/sqrt(expu)); ml = mu/lmd^2

 # Hessian
 Hb <- t(as.numeric(-1/sig^2 - s*g1*ls^2*(ls*(e + s*ml) + s*g1))*X)%*%X
 Hbdel <- t(as.numeric(-s/sig^2 - s*g1*(ls/lmd)^2*(-s*ls*(e + s*ml) - g1))*X)%*%Zdel
 Hbgv <- t(as.numeric(-(e - s*mu)*expv/sig^4 + s*0.5*ls*(1 + expv/sig^2)*g1 - s*ls*gr1v)*X)%*%Zv
 Hbgu <- t(as.numeric(-(e - s*mu)*expu/sig^4 - s*0.5*ls*(1 - expu/sig^2)*g1 - s*ls*gr1u)*X)%*%Zu
 Hdel <- t(as.numeric(g2/expu*(mu/sqrt(expu) + g2) - g1*ls/(lmd^3*sig) * (ls*(ml + s*e) + g1) - sig^(-2))*Zdel)%*%Zdel
 Hdelgv <- t(as.numeric(1/(lmd*sig) * (0.5*(1 - expv/sig^2)*g1 + gr1v) - s*(e - s*mu)*expv/sig^4)*Zdel)%*%Zv
 Hdelgu <- t(as.numeric(-s*(e - s*mu)*expu/sig^4 + 1/(sig*lmd)*(gr1u - 0.5*g1*(1 + expu/sig^2)) - 0.5*g2/sqrt(expu)*(-1 + mu/sqrt(expu)*(mu/sqrt(expu) + g2)))*Zdel)%*%Zu
 Hgv <- t(as.numeric(0.5*expv/sig^2 * ((1 - expv/sig^2)*((e - s*mu)^2/sig^2 - 1) - (e - s*mu)^2*expv/sig^4) + (ml - 0.5*(1 + expv/sig^2)*(ml + s*e))*(gr1v*ls - 0.5*g1*ls*(1 + expv/sig^2)) + g1*ls*(ml - 0.5*(expv/sig^2 * (1 - expv/sig^2)*(ml + s*e) + (1 + expv/sig^2)*ml)))*Zv)%*%Zv
 Hgu <- t(as.numeric((0.5*(1 - expu/sig^2)*(ml + s*e) - ml)*ls*(gr1u + 0.5*g1*(1 - expu/sig^2)) + g1*ls*(0.5*(expu/sig^2 - 1)*(expu/sig^2*(ml + s*e) + ml) + ml) + 0.5*(expu/sig^2 * ((1 - expu/sig^2)*((e - s*mu)^2/sig^2 - 1) - expu*(e - s*mu)^2/sig^4) + 0.5*g2*mu/sqrt(expu)*(-1 + mu/sqrt(expu)*(mu/sqrt(expu) + g2))))*Zu)%*%Zu
 Hgvgu <- t(as.numeric(0.5*expv*expu/sig^4 * (1 - 2*((e - s*mu)/sig)^2) + (ml - 0.5*(1 + expv/sig^2)*(ml + s*e))*ls*(gr1u + 0.5*(1 - expu/sig^2)*g1) + g1*ls*(-ml + 0.5*(expv*expu/sig^4*(ml + s*e) + (1 + expv/sig^2)*ml)))*Zv)%*%Zu

 H <- cbind(rbind(Hb, t(Hbgv),  t(Hbgu),t(Hbdel)), rbind(Hbgv, Hgv, t(Hgvgu), Hdelgv),rbind(Hbgu, Hgvgu, Hgu, Hdelgu), rbind(Hbdel, t(Hdelgv), t(Hdelgu), Hdel ))
 return(H)
}



# Technical efficiencies and prediction intervals
.u2efftnm <- function( e, su, sv, mu, alpha = 0.05, prod) {
 # if(prod){sn = -1} else {sn = 1}
 sn <- ifelse(prod, -1, 1)
 s  <- sqrt(su^2 + sv^2);  m1 <- (sn*su^2 * e + mu*sv^2)/s^2
 s1 <- su * sv / s;  z  <- m1 / s1
 point.est.mean <- m1 + s1 * dnorm(-z) / pnorm(z)
 point.est.mode <- ifelse( m1 >= 0, m1, 0 )
 te_jlms_mean <- exp( sn*point.est.mean)
 te_jlms_mode <- exp( sn*point.est.mode)
 zl    <- qnorm( 1 - alpha / 2 * pnorm(z) )
 zu    <- qnorm( 1 - ( 1 - alpha/2 ) * pnorm(z) )
 te_l  <- exp( sn*m1 - zl * s1 )
 te_u  <- exp( sn*m1 - zu * s1 )
 te_bc <- exp( sn*m1 + .5 * s1^2) * pnorm( sn*s1 + z ) / pnorm(z)
 tymch <- data.frame(te_l, te_jlms_mean, te_jlms_mode,te_bc,te_u)
 colnames(tymch) <- c("Lower bound","JLMS", "Mode", "BC","Upper bound" )
 return(tymch)
}



# Marginal effects
.me = function(theta, Zu, Zdel, ku, kdel, n, dist = c("h", "t")){

 mat.equal <- function(x, y) is.matrix(x) && is.matrix(y) && ncol(x) == ncol(y) && all(colnames(x) == colnames(y))

 gu = theta[1:ku,1,drop = F]
 expu = exp(Zu%*%gu)
 if(dist == "t"){
  delta = theta[-c(1:ku),1,drop = F]
  mu <- Zdel%*%delta
  if(length(delta) == 1) delta = rep(0, max(ku, kdel))
 }
 if(length(gu) == 1) gu = rep(0, max(ku, kdel))
 meff = matrix(NA, ncol = max(ku, kdel) - 1, nrow = n)

 if(dist == "h"){
  if(ncol(Zu) == 1) {
   warning("Marginal effects are not returned: scale of half-normal distribution of inefficiency term is not expressed as function of any exogenous variables", call. = FALSE)
   meff = NULL
  } else {
   arg = sqrt(1/(2*pi)) * sqrt(expu)
   for(i in 2:ku){
    meff[,i-1] = arg*gu[i]
    meff = round(meff, digits = 5)
   }
   colnames(meff) = rownames(gu)[-1]
  }} else if(ncol(Zu) == 1 & ncol(Zdel) == 1){
   warning("Marginal effects are not returned: mean or variance of (pre-)truncated normal distribution of inefficiency term is not expressed as function of any exogenous variables", call. = FALSE)
   meff = NULL
  }else  if(dist == "t" & (mat.equal(Zu, Zdel)| all(delta == 0) | all(gu == 0))){
   arg = mu/sqrt(expu)
   g = dnorm(arg)/pnorm(arg)
   arg1 = (1 - arg*g - g^2); arg2 = sqrt(expu)/2 * ((1 + arg^2)*g + arg*g^2)
   for(i in 2:max(ku, kdel)){
    meff[,i-1] = arg1*delta[i] + arg2*gu[i]
    meff = round(meff, digits = 5)
   }
   if(all(gu == 0)){colnames(meff) = rownames(delta)[-1]} else {colnames(meff) = rownames(gu)[-1]}
  }  else {
   warning("Marginal effects are not returned: mean and variance of (pre-)truncated normal distribution of inefficiency term are expressed as functions not of the same exogenous variables", call. = FALSE)
   meff = NULL
  }
 return( meff)
}


.skewness <- function (x, na.rm = FALSE, type = 3)
{
 if (any(ina <- is.na(x))) {
  if (na.rm)
   x <- x[!ina]
  else return(NA)
 }
 if (!(type %in% (1:3)))
  stop("Invalid 'type' argument.")
 n <- length(x)
 # x <- x - mean(x)
 y <- sqrt(n) * sum(x^3)/(sum(x^2)^(3/2))
 if (type == 2) {
  if (n < 3)
   stop("Need at least 3 complete observations.")
  y <- y * sqrt(n * (n - 1))/(n - 2)
 }
 else if (type == 3)
  y <- y * ((1 - 1/n))^(3/2)
 y
}

# Print the estimation results

.printoutcs = function(x, digits, k, kv, ku, kdel, na.print = "NA", dist, max.name.length, mycutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), mysymbols = c("***", "**", "*", ".", " ")){

 Cf = cbind(ifelse(x[,1, drop = FALSE]> 999, formatC(x[,1, drop = FALSE], digits = 1, format = "e",width = 10), formatC(x[,1, drop = FALSE], digits = digits, format = "f", width = 10)),
            ifelse(x[,2, drop = FALSE]>999, formatC(x[,2, drop = FALSE], digits = 1, format = "e", width = 10), formatC(x[,2, drop = FALSE], digits = digits, format = "f", width = 10)),
            ifelse(x[,3, drop = FALSE]>999, formatC(x[,3, drop = FALSE], digits = 1, format = "e", width = 7), formatC(x[,3, drop = FALSE], digits = 2, format = "f", width = 7)),
            ifelse(x[,4, drop = FALSE]>999, formatC(x[,4, drop = FALSE], digits = 1, format = "e", width = 10), formatC(x[,4, drop = FALSE], digits = digits, format = "f", width = 10)),
            formatC(mysymbols[findInterval(x = x[,4], vec = mycutpoints)], flag = "-"))

 # mycutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1)
 # mysymbols = c("***", "**", "*", ".", " ")
 #
 # pvals <- m0$table[,4]
 #
 # findInterval(x = pvals, vec = mycutpoints)
 #
 # cbind(pvals,mysymbols[findInterval(x = pvals, vec = cutpoints)])
 #
 # pval_sym <- mysymbols[findInterval(x = x[,4], vec = cutpoints)]


 # cat("               Coef.        SE       z       P>|z|\n", sep = "")
 row.names(Cf) <- formatC(row.names(Cf), width = max(nchar(row.names(Cf))), flag = "-")
 cat("",rep(" ", max.name.length+6),"Coef.        SE       z       P>|z|\n", sep = "")
 dimnames(Cf)[[2]] <- rep("", dim(Cf)[[2]])
 cat("",rep("_", max.name.length+42-1),"", "\n", "Frontier", "\n", sep = "")
 print.default(Cf[1:k,,drop=FALSE], quote = FALSE, right = TRUE, na.print = na.print)
 cat("",rep("-", max.name.length+42-1),"", "\n", "Random noise component: log(sigma_v^2)", "\n", sep = "")
 # dimnames(Cf)[[2]] <- rep("", dim(Cf)[[2]])
 print.default(Cf[(k+1):(k+kv),,drop=FALSE], quote = FALSE, right = TRUE, na.print = na.print)
 cat("",rep("-", max.name.length+42-1),"", "\n", "Inefficiency component: log(sigma_u^2)", "\n", sep = "")
 print.default(Cf[(k+kv+1):(k+kv+ku),,drop=FALSE], quote = FALSE, right = TRUE, na.print = na.print)
 if(dist == "t"){
  cat("",rep("-", max.name.length+42-1),"", "\n", "Mu (location parameter)", "\n", sep = "")
  print.default(Cf[(k+kv+ku+1):(k+kv+ku+kdel),,drop=F], quote = FALSE, right = TRUE, na.print = na.print)
 }
 if(nrow(Cf[-c(1:(k+kv+ku+kdel)),,drop=FALSE]) >= 1){
  cat("",rep("-", max.name.length+42-1),"", "\n", "Parameters of compound error distribution", "\n", sep = "")
  print.default(Cf[-c(1:(k+kv+ku+kdel)),,drop=FALSE], quote = FALSE, right = TRUE, na.print = na.print)
 }
 cat("",rep("_", max.name.length+42-1),"", "\n", sep = "")
 cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
 invisible(x)
}



.timing <- function(x, wording = "Time elapsed is")
{
 if(x < 60){
  # cat("\n")
  cat("",wording,"",x," seconds\n", sep = "")
  # cat("\n")
 } else {
  if(x >= 60 & x < 60*60){
   minutes <- floor(x/60)
   seconds <- round(x - minutes * 60,1)
   # cat("\n")
   cat("",wording,"",minutes," minute(s) and ",seconds," second(s)\n", sep = "")
   # cat("\n")
  } else {
   if(x >= 60*60 & x < 60*60*24){
    hours   <- floor(x / 60 / 60)
    minutes <- round( (x - hours * 60 *60) / 60, 1)
    seconds <- floor(x - hours * 60 *60 - minutes * 60)
    # cat("\n")
    cat("",wording,"",hours," hour(s) and ",minutes," minute(s) \n", sep = "")
    # cat("\n")
   } else {
    if(x >= 60*60*24){
     days    <- floor(x / 60 / 60 / 24)
     hours   <- round( (x - days * 60 * 60 * 24) / 60 /60 ,1)
     minutes <- floor( (x - days * 60 * 60 * 24 - hours * 60 *60) / 60)
     seconds <- floor(x - days * 60 * 60 * 24 - hours * 60 *60 - minutes * 60)
     # cat("\n")
     cat("",wording,"",days," day(s) and ",hours," hour(s)\n", sep = "")
     # cat("\n")
    }
   }
  }
 }
}

# library(matrixcalc)


.mlmaximize <- function(theta0, ll, gr = NULL, hess = NULL, alternate = NULL, BHHH = F, level = 0.99, step.back = .Machine$double.eps, reltol = .Machine$double.eps, lmtol = sqrt(.Machine$double.eps), steptol = sqrt(.Machine$double.eps), digits = 4, when.backedup = sqrt(.Machine$double.eps), max.backedup = 17, print.level = 6, only.maximize = FALSE, maxit = 150, n = 100, ...){

 theta00 <- theta0

 k4 <- length(theta0)

 if(print.level >= 6){
  cat("\n=================")
  cat(" Initial values:\n\n", sep = "")
  print(theta0)
 }

#  if(print.level >= 2){
#   cat("\n=================")
#   cat(" Maximization:\n\n", sep = "")
#  }

 # step.back = 2^-217

 ll0 <- ll(theta0, ...)
 ltol <- reltol * (abs(ll0) + reltol)
 typf <- ll0
 theta1 <- theta0

 iter <- iter.total <- backedup <- backedups <- wasconcave <- wasconcaves <- 0

 if( is.na(ll0) | ll0 == -Inf ){
  if(print.level >= 2){
   cat("Could not compute ll at starting values: trying something else\n")
  }
  iter1 <- backedups
  repeat{
   iter1 <- iter1 + 1
   theta0 <- theta0 * runif(length(theta0), 0.96, 1.05) # not sure what to do here
   ll0 <- ll( theta0, ... )
   if(!is.na(ll0)) break
   if(iter1 == 55){
    stop("it's not gonna happen...")
   }
  }
  # backedups <- iter1
 }

 delta1 <- gHg <- s1 <- 1
 h1 <- tryCatch( 2, error = function(e) e )
 cant.invert.hess <- FALSE

 if(print.level >= 2){
  cat(paste("Iteration ",formatC(iter, width = 3)," (at starting values):       log likelihood = ",format(ll0, digits = 13),"\n\n", sep = ""), sep = "")
 }

 repeat{
  iter.total <- iter.total + 1
  # cat("backedup = ",backedup," backedups = ",backedups,"\n", sep = "")
  if(print.level >= 6){
   print(theta0)
  }

  # cumulate how many times did it backed-up in a row
  if(s1 < when.backedup){
   backedup <- backedup + 1
  } else {
   backedup <- 0
  }
  # print(s1)
  # cumulate how many times was concave
  if( inherits(h1, "error") ){
   wasconcave <- wasconcave + 1
  } else {
   wasconcave <- 0
  }

  # try different values if was concave more than @@@ times
  if(wasconcave == max.backedup){
   # start over
   if(print.level >= 2){
    cat("Not concave ",max.backedup," times in a row: trying something else (not concave ",wasconcaves+1," times in total)\n")
   }
   iter <- wasconcave <- backedup <- 0
   wasconcaves <- wasconcaves + 1
   theta0 <- theta0*runif(length(theta0), 0.96, 1.05) # not sure what to do here
   ll0 <- ll( theta0, ... )
   # if(print.theta) print(theta0)
   if(print.level >= 2){
    cat(paste("Iteration ",formatC(iter, width = 3),"  (at slightly perturbed starting values):       log likelihood = ",format(ll0, digits = 13),"\n\n", sep = ""), sep = "")
   }
  }

  # try different values if backed-up more than @@@ times
  if(backedup == max.backedup){
   # start over
   if(print.level >= 2){
    cat("Backed-up ",max.backedup," times in a row: trying something else (backup-up ",backedups+1," times in total)\n", sep = "")
   }
   iter <- backedup <- wasconcave <- 0
   backedups <- backedups + 1
   wasconcaves <- wasconcaves + 1
   theta0 <- theta0*runif(length(theta0), 0.96, 1.04) # not sure what to do here
   ll0 <- ll( theta0, ... )
   if(print.level >= 6){
    print(theta0)
   }
   if(print.level >= 2){
    cat(paste("Iteration ",formatC(iter, width = 3)," (at slightly perturbed starting values):       log likelihood = ",format(ll0, digits = 13),"\n\n", sep = ""), sep = "")
   }
  }

  # see if it calculated ll
  if( is.na(ll0) | ll0 == -Inf | ll0 == Inf | ll0 == 0 ){
   if(print.level >= 2){
    cat("Could not compute ll: trying something else\n")
   }
   iter1 <- backedups
   repeat{
    iter1 <- iter1 + 1
    # theta0 <- c( cons0, beta0, mu = 0, eta = 0, lnsv2 = -1*iter1/2, lnsu2 = -1*iter1/2)
    theta0 <- theta00*runif(length(theta0), 0.96, 1.04) # not sure what to do here
    ll0 <- ll( theta0, ... )
    if(!is.na(ll0) & ll0 != 0) break
    if(iter1 == 15){
     stop("it's not gonna happen... could not compute at initial and find feasible values")
    }
   }
   iter <- backedup <- 0
   backedups <- iter1
   if(print.level >= 2){
    cat(paste("Iteration ",formatC(iter, width = 3)," (at slightly perturbed starting values):       log likelihood = ",format(ll0, digits = 13),"\n\n", sep = ""), sep = "")
   }
  }
  if(backedups == 5){
   stop("it's not gonna happen... backuped 5 times")
  }
  if(wasconcaves == 5){
   stop("it's not gonna happen... not concave 5 times")
  }
  iter <- iter + 1
  delta3 <- 1
  # step 2: direction vector

  # previous theta

  if(iter.total > 1) theta1 <- theta0 - s1 * d0

  # BHHH (faster, but different):
  # The Hessian is approximated as the negative
  # of the sum of the outer products of the gradients
  # of individual observations,
  # -t(gradient) %*% gradient = - crossprod( gradient )
  g1 <- gr(theta0, ...)
  # print(g1)
  if(!is.null(alternate)) BHHH <- floor(iter/alternate) %% 2 == 1
  h0 <- hess(theta0,  ...)
  # print(h0)
  # check if negative definite
  eigen1 <- eigen( h0 )
  # eigen.tol <- k4 * max(abs(eigen1$values)) * .Machine$double.eps # this is for positive definiteness
  eigen.val <- ifelse(eigen1$values < .Machine$double.eps^.1, 0, eigen1$values)
  # hess.pos.def <- sum(eigen1$values > eigen.tol) == k4
  hess.neg.def <- !any(eigen.val >= 0)
  # 1. replace negative with small ones
  # eigen.val <- ifelse(eigen1$values < 0, .0001, eigen.val)
  # 2. replace negative with absolut values
  eigen.val <- abs(eigen1$values)
  # print(hess.neg.def)
  # make it negative definite if it is not already
  if(!hess.neg.def){
   h0_ <- matrix(0, k4, k4)
   # eigen1 <- eigen( h0 )
   for( i in seq_len( k4 ) ){
    # h0_ <- h0_ - abs(eigen1$values[i]) * eigen1$vectors[,i] %*% t(eigen1$vectors[,i])
    h0_ <- h0_ - eigen.val[i] * eigen1$vectors[,i] %*% t(eigen1$vectors[,i])
   }
   h0.old <- h0
   h0 <- h0_
  }
  # print( is.negative.definite(h0) )
  # remember hessian and negative of its inverse from previous iter that could have been inverted
  if( !cant.invert.hess ){
   h0_previous <- h0
   h1_previous <- h1
  }
  # easier to invert positive definite matrix
  h1 <- tryCatch( qr.solve(-h0), error = function(e) e )
  # check if it can be inverted
  cant.invert.hess <- FALSE
  cant.invert.hess <- inherits(h1, "error")
  if( cant.invert.hess ){
   # print(h1)
   if(print.level >= 2){
    cat(paste("cannot invert Hessian, using eigenvalues\n", sep = ""), sep = "")
   }
   # this was just to get the uninvertable hessian
   # return(list(hess = h0, grad = g1))
   # stop("this")
   # @14@ this
   eig1 <- eigen( -h0_previous )
   d0 <- rep(0, length(theta0))
   # eig2 <- ifelse(eig1$values < eps1, 1, eig1$values)
   for (i in 1:length(eig1$values)){
    d0 <- d0 + (g1%*%eig1$vectors[,i])%*%eig1$vectors[,i] / eig1$values[i]
   }
   # @14@ could be done easier
   # d0 <- qr.solve(-h0, g1, tol = 1e-10)
   gHg <- sum( g1 * d0)
   # in the part of the ortogonal subspace where the eigenvalues
   # are negative or small positive numbers, use steepest ascent
   # in other subspace use NR step
   # d0 <- ifelse(eigen(-h0, only.values = TRUE)$values < reltol, g1, d0)
   gg <- sqrt( crossprod(g1) )
   gHg <- gg
   # d0 <- g1
   # d0
  } else {
   d0 <- as.vector( h1 %*% g1 )
   gg <- sqrt( crossprod(g1) )
   # h1.old <- solve(-h0.old)
   gHg <- as.vector( t(g1) %*% h1 %*% g1 )
  }
  # gg_scaled <- gg * max( crossprod(theta0), crossprod(theta1) ) / max( abs(ll0), abs(typf))
  # theta_rel <- max( abs(theta0 - theta1) / apply( cbind( abs(theta0),abs(theta1) ), 1, max) )
  theta_rel <- max( abs(theta0 - theta1) / (abs(theta1)+1) )


  # begin stopping criteria calculated using new values of g1 and h1
  if(s1 > when.backedup*10^-100 & delta1 != 17.17){ # if(s1 > when.backedup*10^-100 & !cant.invert.hess){
   if(abs(gHg) < lmtol & iter.total > 1){
    conv.crite <- 1
    if(print.level >= 2){
     cat("\nConvergence given g inv(H) g' = ",abs(gHg)," < lmtol\n", sep = "")
    }
    break
   }
   if(theta_rel < steptol & iter.total > 2){
    conv.crite <- 3
    # print(theta_rel)
    if(print.level >= 2){
     cat("\nConvergence given relative change in parameters = ",theta_rel," < steptol\n", sep = "")
    }
    break
   }
  }
  # end stopping criteria
  # use steepest ascent when backed-up
  if(s1 < when.backedup*10^-3){
   # eig1 <- eigen( -h0 )
   d0 <- g1
   # d0 <- ifelse(eig1$values < reltol, g1, d0)
   # theta0 <- theta0 - 1 * d0
  }
  # print(d0)
  # step 3: new guess
  # a: s = 1
  # b: funct(theta0 + d0) > funct(theta0)
  s1 <- 1
  theta1 <- theta0 + s1 * d0
  # print(12)
  # print(theta1)
  ll1 <- ll( theta1, ... )
  # print(13)
  delta2 <- ll1 - ll0
  flag <- (!is.na(delta2) & delta2 != -Inf & delta2 > 0)
  # begin Cases
  if( flag ){
   # begin Case 1: f(theta1) > f(theta0)
   ll.temp <- ll0
   # check if s1 = 2, 3, ... increases f even more
   while( flag ){
    if(print.level >= 6){
     cat(paste("\t\tCase 1: s = ",s1,", delta = ",delta2,"\n", sep = ""), sep = "")
    }
    ll0 <- ll1
    s1 <- s1 + 1
    theta1 <- theta0 + s1 * d0
    ll1 <- ll( theta1, ... )
    delta2 <- ll1 - ll0
    flag <- (!is.na(delta2) & delta2 != -Inf & delta2 > 0)
   }
   # overall delta
   delta1 <- ll0 - ll.temp
   delta_rel <- abs(delta1 / ll.temp)
   # print(delta_rel)
   s1 <- s1 - 1
   # overwrite the values
   theta0 <- theta0 + s1 * d0
   # end Case 1: f(theta1) > f(theta0)
  } else {
   # begin Case 2: f(theta1) < f(theta0)
   # check only if s1=1/2 increases f
   s1 <- 0.5
   theta1 <- theta0 + s1 * d0
   ll1 <- ll( theta1, ... )
   # cat(" ll1 = ",ll1,", ll0 = ",ll0,", s = ",s1,"\n", sep = "")
   delta2 <- ll1 - ll0
   if(print.level >= 6){
    cat(paste("\t\tCase 2: s = ",s1,", delta = ",delta2,"\n", sep = ""), sep = "")
   }
   flag2 <- (!is.na(delta2) & delta2 != -Inf & delta2 > 0)
   # end Case 2: f(theta1) < f(theta0)
   if( flag2 ){
    # begin Case 2a: f(theta1) > f(theta0)
    ll.temp <- ll0
    # check if s1=1/2^2,1/2^3,... increases f even more
    while( flag2 ){
     if(print.level >= 6){
      cat(paste("\t\t\tCase 2a: s = ",s1,", delta = ",delta2,"\n", sep = ""), sep = "")
     }
     ll0 <- ll1
     s1 <- 0.5 * s1
     theta1 <- theta0 + s1 * d0
     ll1 <- ll( theta1, ... )
     # cat(" ll1 = ",ll1,", ll0 = ",ll0,", s = ",s1,"\n", sep = "")
     delta2 <- ll1 - ll0
     flag2 <- (!is.na(delta2) & delta2 != -Inf & delta2 > ltol)
    }
    # overall delta
    delta1 <- ll0 - ll.temp
    delta_rel <- abs(delta1 / ll.temp)
    # print(delta_rel)
    s1 <- 2 * s1
    # overwrite the values
    theta0 <- theta0 + s1 * d0
    # end Case 2a: f(theta1) > f(theta0)
   } else {
    # begin Case 2b: f(theta1) < f(theta0)
    ll.temp <- ll0
    # try s1=1/2^2,1/2^3,... so that f(theta1) > f(theta0)
    while ( !flag2 & s1 > step.back ){
     s1 <- 0.5 * s1
     theta1 <- theta0 + s1 * d0
     ll1 <- ll( theta1, ... )
     delta2 <- ll1 - ll0
     if(print.level >= 6){
      cat(paste("\t\t\tCase 2b: s = ",s1,", delta = ",delta2,"\n", sep = ""), sep = "")
     }
     flag2 <- (!is.na(delta2) & delta2 != -Inf & delta2 > 0)
    }
    if( !flag2 | s1 < step.back ){
     # stop("provide different starting values")
     delta1 <- 17.17
    } else {
     # overwrite the values
     delta1 <- delta2
     delta_rel <- abs(delta1 / ll.temp)
     ll0 <- ll1
     theta0 <- theta0 + s1 * d0
    }
    # end Case 2b: f(theta1) < f(theta0)
   }
  }


  if(print.level >= 2){
   if( cant.invert.hess ){
    cat(paste("Iteration ",formatC(iter, width = 3)," (hessian is ",ifelse(BHHH, "BHHH", "analytical"),", ",formatC(iter.total, width = 3)," in total):   log likelihood = ",format(ll0, digits = 13)," (not concave)\n", sep = ""), sep = "")
   } else if (s1 <= when.backedup) {
    cat(paste("Iteration ",formatC(iter, width = 3)," (hessian is ",ifelse(BHHH, "BHHH", "analytical"),", ",formatC(iter.total, width = 3)," in total):   log likelihood = ",format(ll0, digits = 13)," (backed up)\n", sep = ""), sep = "")
   } else {
    cat(paste("Iteration ",formatC(iter, width = 3)," (hessian is ",ifelse(BHHH, "BHHH", "analytical"),", ",formatC(iter.total, width = 3)," in total):   log likelihood = ",format(ll0, digits = 13),"\n", sep = ""), sep = "")
   }
  }

  # printing criteria
  if(print.level >= 5){
   if( cant.invert.hess ){
    cat(paste(" (in iter ",formatC(iter, width = 3),": delta = ",format(delta1, digits = 6),"; s = ",format(s1, digits = 6),"; quasi-gHg = ",format(gHg, digits = 6),"; theta_rel_change = ",format(theta_rel, digits = 6),")\n\n", sep = ""), sep = "")
    if(print.level >= 5.5){
     print(theta0)
     cat("\n")
    }
   } else {
    cat(paste(" (in iter ",formatC(iter, width = 3),": delta = ",format(delta1, digits = 6),"; s = ",format(s1, digits = 6),"; gHg = ",format(gHg, digits = 6),"; theta_rel_change = ",format(theta_rel, digits = 6),")\n\n", sep = ""), sep = "")
    if(print.level >= 5.5){
     print(theta0)
     cat("\n")
    }
   }
  }
  # print(s1)
  if(s1 > when.backedup^2 & !cant.invert.hess){ # if(s1 > when.backedup^2 & !cant.invert.hess){
   # ltol <- reltol * (abs(ll0) + reltol)
   # print(cant.invert.hess)
   if(delta1 > 0 & !is.na(delta_rel) & delta_rel < ltol & iter.total > 1){
    conv.crite <- 2
    if(print.level >= 2){
     cat("\nConvergence given relative change in log likelihood = ",delta_rel," < ltol\n", sep = "")
    }
    break
   }
  }
  if(iter.total > maxit){
   cat("\n Maximum number of iterations (",maxit,") reached without convergence\n", sep = "")
   break
  }
 } # end repeat

 if( !only.maximize & !cant.invert.hess){
  names(ll0) <- NULL
  colnames(h1) <- rownames(h1) <- names(g1) <- names(theta0)

  # sqrt(crossprod(g1))

  b0 <- theta0
  sd0 <- sqrt( diag( h1 ) )
  t0 <- b0 / sd0
  p0 <- pt(abs(t0), n-length(b0), lower.tail = FALSE) * 2
  t10 <- qt((1-0.01*level)/2, n-length(b0), lower.tail = FALSE)
  t17 <- cbind( b0, sd0, t0, p0, b0 - t10*sd0, b0 + t10*sd0)
  # t17 <- cbind( b0, sd0, t0, p0)
  colnames(t17) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)", paste("",level,"_CI_LB", sep = ""), paste("",level,"_CI_UB", sep = ""))
  # colnames(t17) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  # t17

  if(print.level >= 2){
   cat(paste("\nFinal log likelihood = ",formatC(ll0,digits=7,format="f"),"\n\n", sep = ""), sep = "")
  }
  # cat(paste("Stoc. frontier normal/",distribution,"\n", sep = ""), sep = "")
  if(print.level >= 5.5){
   cat("\nCoefficients:\n\n", sep = "")
   printCoefmat(t17[,1:4], digits = digits)
  }

  return(list(par = theta0, table = t17, gradient = g1, vcov = h1, ll = ll0, gg = gg, gHg = gHg, delta_rel = delta_rel, ltol = ltol, theta_rel_ch = theta_rel, conv.crite = conv.crite))
 } else {
  return(list(par = theta0, gradient = g1, vcov = h1, ll = ll0, gg = gg, gHg = gHg, delta_rel = delta_rel, ltol = ltol, theta_rel_ch = theta_rel, conv.crite = conv.crite))
 }


}
.su1 <- function(x, mat.var.in.col = TRUE, digits = 4, probs = c(0.1, 0.25, 0.5, 0.75, 0.9), print = TRUE, transpose = FALSE){

 xvec2 <- xvec1 <- FALSE

 if(is.matrix(x)){
  if(min(dim(x)) == 1){
   xvec1 <- TRUE
   # mynames <- deparse(substitute(x))
   x <- as.vector(x)
  } else {
   if(!mat.var.in.col){
    x <- t(x)
   }
   mynames <- paste("Var", seq_len(ncol(x)), sep = "")
  }
  # print(x)
  # mynames <- colnames(x)
 } # end if matrix

 if(is.vector(x)){
  xvec2 <- TRUE
  mynames <- deparse(substitute(x))
  x <- data.frame(Var1 = x)
 } # end if vector

 # cat("nymanes", sep ="")
 # print(mynames)

 if(!is.vector(x) & !is.matrix(x) & !is.data.frame(x)){
  stop("Provide vector, matrix, or data.frame")
 } else {
  t1 <- apply(x, 2, function(x) c(Obs = length(x), NAs = sum(is.na(x)), Mean = mean(x, na.rm = TRUE), StDev = sd(x, na.rm = TRUE), IQR = IQR(x, na.rm = TRUE), Min = min(x, na.rm = TRUE), quantile( x, probs = probs, na.rm = TRUE ), Max = max(x, na.rm = TRUE)))
  # print(class(t1))
  # print(dim(t1))
  if(xvec2 & !xvec1) colnames(t1) <- mynames

  if(print){
   if(transpose){
    tymch <- formatC(t1, digits = 4, format = "f", width = 4+1)
   } else {
    tymch <- formatC(t(t1), digits = 4, format = "f", width = 4+1)
   }
   tymch <- gsub(".0000","",tymch, fixed = TRUE)
   print(tymch, quote = FALSE)
  }
 }
 return(t1)
}

.prepareYXZ.panel.simple <- function(formula, data, it, subset,
                                     ln.var.v.it = ln.var.v.it, 
                                     ln.var.u.0i = ln.var.u.0i, 
                                     mean.u.0i = mean.u.0i,
                                     print.level = 1, ...) {
 needed.frame <- sys.nframe() - 1
 mf0 <- match.call(expand.dots = FALSE, call = sys.call(sys.parent(n = 1)))
 # print(mf0)
 
 # check if it is a matrix
 datasupplied <- !(match("data", names(mf0), 0) == 0)
 subssupplied <- !(match("subset", names(mf0), 0) == 0)
 
 if(subssupplied & !datasupplied){
  stop("Cannot specify 'subset' without specifying 'data'\n")
 }
 
 itsupplied <- !(match("it", names(mf0), 0) == 0)
 if(!itsupplied){
  stop("Panel structure must be specified by using 'it' argument\n")
 }
 
 if(length(it) != 2){
  stop("Invalid panel structure in 'it' argument\n")
 }
 
 if ( is.null(ln.var.v.it) ){
  ln.var.v.it <- ~ 1
  kvi <- 1
 } else {
  kvi <- 17
 }
 if ( is.null(ln.var.u.0i) ){
  ln.var.u.0i <- ~ 1
  ku0 <- 1
 } else {
  ku0 <- 17
 }
 if ( is.null(mean.u.0i) ){
  mean.u.0i <- ~ 1
  kdeli <- 1
 } else {
  kdeli <- 17
 }
 
 form3 <- as.formula(paste0("~", paste0(it, collapse = "+")))
 form1 <- as.Formula(formula,ln.var.v.it,ln.var.u.0i,mean.u.0i,form3)
 
 # print(form2)
 
 if(datasupplied){
  data.order <- as.numeric(rownames(data)) # seq_len(nrow(data))
  
  # print(rownames(data)[1:100])
  # print(as.numeric(rownames(data))[1:100])
  
  # begin get a logical vector equal TRUE if !missing
  # first using data and subset to get x without NA
  mf <- mf0
  mf$formula <- form1 #formula( form )
  m <- match(c("formula", "data", "subset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, sys.frame(sys.parent(n = needed.frame)))
  mt <- attr(mf, "terms")
  #   cat("Print mt\n")
  #   print(mt)
  #   cat("Print mf\n")
  #   print(mf)
  X <- model.matrix(mt, mf)
  # Y <- as.matrix(model.response(mf, "numeric"))
  # print(dim(X))
  # X <- model.frame(mt, mf)
  #   cat("Print X\n")
  #   print(X)
  #   cat("End print X\n")
  # now get the names in the entire data
  esample.nu <- as.numeric(rownames(X))
  # print(esample.nu[1:100])
  # print(old.order)
  esample <- data.order %in% esample.nu
  # print(length(esample))
  # print(table(esample))
  if(sum(esample)==0){
   stop("no valid data points")
  }
  # print(table(esample))
  # end get a logical vector equal TRUE if !missing
  
  # get the data
  # print(form1)
  #   dataesample <- model.frame(form1, data = data)
  #   dataesample <- dataesample[esample.nu,]
  #   print(as.numeric(rownames(dataesample)))
  # esample1 <- esample
  #   esample <- as.numeric(rownames(dataesample)) %in% esample.nu
  #   # print(table(esample))
  #   print(c(length(esample), nrow(dataesample)))
  # dataesample <- dataesample[as.numeric(rownames(dataesample)) %in% esample.nu,]
  # print(dim(dataesample))
  # print(dataesample)
  Y <- as.matrix(model.part(form1, data = mf, lhs = 1, drop = FALSE))
  # print(length(Y))
  # print(length(model.response(dataesample)))
  # Y <- as.matrix( model.matrix(formula(form1, lhs = 1, rhs = 0), data = data[esample,]))
  # print(Y)
  X <- as.matrix( model.matrix(formula(form1, lhs = 0, rhs = 1), data = mf))
  # print(X)
  nt <- nrow(Y)
  k <- ncol(X)
  # Zvi
  if(is.null(ln.var.v.it)){
   Zvi <- matrix(1, nrow = sum(esample), ncol = 1)
  }
  else {
   Zvi <- as.matrix( model.matrix(formula(form1, lhs = 0, rhs = 2), data = mf))
   kvi <- ncol(Zvi)
  }
  # Zu0
  if(is.null(ln.var.u.0i)){
   Zu0 <- matrix(1, nrow = sum(esample), ncol = 1)
  }
  else {
   Zu0 <- as.matrix( model.matrix(formula(form1, lhs = 0, rhs = 3), data = mf))
   ku0 <- ncol(Zu0)
  }
  # kdeli
  if(is.null(mean.u.0i)){
   zdeli <- matrix(1, nrow = sum(esample), ncol = 1)
  }
  else {
   zdeli <- as.matrix( model.matrix(formula(form1, lhs = 0, rhs = 4), data = mf))
   kdeli <- ncol(zdeli)
  }
  # IT
  itvar <- as.matrix( model.matrix(formula(form1, lhs = 0, rhs = 5), data = mf))[,-1]
  sorted <- order(itvar[,1],itvar[,2])
  old.order <- seq(nrow(itvar))
  old.sorted <- order(old.order[sorted])
  itvar <- itvar[sorted,]
  idvar <- itvar[,1]
  timevar <- itvar[,2]
  ii <- unique(idvar)
  n <- length(ii)
  t.years <- length(unique(itvar[,2]))
  
  # The size of each panel
  
  t0 <- as.vector( by(data = itvar[,1], INDICES = itvar[,1], function(x) sum(!is.na(x))) )
  t1 <- summary( t0 )
  
  # summary of the panel data
  
  # if(print.level >= 1){
  #
  #  cat("\n=================")
  #  cat(" Summary of the panel data: \n\n", sep = "")
  #  cat("   Number of obs       (NT) =",nt,"", "\n")
  #  cat("   Number of groups     (N) =",n,"", "\n")
  #  cat("   Obs per group: (T_i) min =",t1["Min."],"", "\n")
  #  cat("                        avg =",t1["Mean"],"", "\n")
  #  cat("                        max =",t1["Max."],"", "\n\n")
  # }
  
  ids <- sort( unique(idvar) )
  idlenmax <- t1["Max."]
  dat.descr <- c(NT = nt, N = n, t1["Min."], t1["Mean"], t1["Max."])
  #   print(head(Y))
  #   print(head(X))
  #   print(head(Zu))
  #   print(head(Zv))
  #   print(head(Zdel))
 }
 # if data are not supplied
 else {
  # begin get a logical vector equal TRUE if !missing
  
  # first using data and subset to get XZ without NA
  mf <- mf0
  mf$formula <- formula( formula )
  subsetsupplied <- !(match("subset", names(mf), 0) == 0)
  if(subsetsupplied) stop("Subset with matrices is not allowed")
  m <- match(c("formula"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- as.name("model.frame")
  # mf <- eval(mf, parent.frame())
  mf <- eval(mf, sys.frame(sys.parent(n = needed.frame)))
  mt <- attr(mf, "terms")
  Y <- as.matrix(model.response(mf))
  n <- nrow(Y)
  XZ <- as.matrix(model.matrix(mt, mf))
  # get a logical vector equal TRUE if !missing
  with.na <- model.frame(mt, na.action = na.pass)
  esample.nu <- as.numeric(rownames(XZ))
  esample <- rownames(with.na) %in% esample.nu
  if(sum(esample)==0){
   stop("no valid data points")
  }
  # print(table(esample))
  # end get a logical vector equal TRUE if !missing
  
  # get the data
  #   print(head(Y))
  #   print(head(XZ))
  
  # get X
  mf <- mf0
  m <- match(c("formula"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, sys.frame(sys.parent(n = needed.frame)))
  mt <- attr(mf, "terms")
  X <- as.matrix(model.matrix(mt, mf))
  # print(head(X))
  k <- ncol(X)
  
  # get Zvi
  if(is.null(ln.var.v.it)){
   Zvi <- matrix(1, nrow = sum(esample), ncol = 1)
  }
  else {
   mf <- mf0
   mf$formula <- formula( Formula(as.formula(paste("",deparse(formula)," + ",ln.var.v.it[2],"", sep = ""))))
   m <- match(c("formula"), names(mf), 0L)
   mf <- mf[c(1L, m)]
   mf[[1L]] <- as.name("model.frame")
   mf <- eval(mf, sys.frame(sys.parent(n = needed.frame)))
   mt <- attr(mf, "terms")
   XZs <- as.matrix(model.matrix(mt, mf))
   Zvi <- XZs[, -(2:k), drop = FALSE]
   kvi <- ncol(Zvi)
   # print(head(Zu))
  }
  # get Zu0
  if(is.null(ln.var.u.0i)){
   Zu0 <- matrix(1, nrow = sum(esample), ncol = 1)
  }
  else {
   mf <- mf0
   mf$formula <- formula( Formula(as.formula(paste("",deparse(formula)," + ",ln.var.v.it[2]," + ",ln.var.u.0i[2],"", sep = ""))))
   m <- match(c("formula"), names(mf), 0L)
   mf <- mf[c(1L, m)]
   mf[[1L]] <- as.name("model.frame")
   mf <- eval(mf, sys.frame(sys.parent(n = needed.frame)))
   mt <- attr(mf, "terms")
   XZs <- as.matrix(model.matrix(mt, mf))
   Zu0 <- XZs[, -(2:(k+kvi-1)), drop = FALSE]
   ku0 <- ncol(Zu0)
   # print(head(Zv))
  }
  # get zdeli
  if(is.null(mean.u.0i)){
   zdeli <- matrix(1, nrow = sum(esample), ncol = 1)
  }
  else {
   mf <- mf0
   mf$formula <- formula( Formula(as.formula(paste("",deparse(formula)," + ",ln.var.v.it[2]," + ",ln.var.u.0i[2]," + ",mean.u.0i[2],"", sep = ""))))
   m <- match(c("formula"), names(mf), 0L)
   mf <- mf[c(1L, m)]
   mf[[1L]] <- as.name("model.frame")
   mf <- eval(mf, sys.frame(sys.parent(n = needed.frame)))
   mt <- attr(mf, "terms")
   XZs <- as.matrix(model.matrix(mt, mf))
   zdeli <- XZs[, -(2:(k+kvi-1+ku0-1)), drop = FALSE]
   kdeli <- ncol(zdeli)
   # print(head(Zv))
  }
  # IT
  itvar <- XZ[, -(2:(k+kvi-1+ku0-1+kvi-1+kdeli-1)), drop = FALSE]
  sorted <- order(itvar[,1],itvar[,2])
  itvar <- itvar[sorted,]
  idvar <- itvar[,1]
  timevar <- itvar[,2]
  ii <- unique(idvar)
  n <- length(ii)
  t.years <- length(unique(itvar[,2]))
  
  # The size of each panel
  
  t0 <- as.vector( by(data = itvar[,1], INDICES = itvar[,1], function(x) sum(!is.na(x))) )
  t1 <- summary( t0 )
  
  # summary of the panel data
  
  # if(print.level >= 1){
  #  cat("\n=================")
  #  cat(" Summary of the panel data: \n\n", sep = "")
  #  cat("   Number of obs       (NT) =",nt,"", "\n")
  #  cat("   Number of groups     (N) =",n,"", "\n")
  #  cat("   Obs per group: (T_i) min =",t1["Min."],"", "\n")
  #  cat("                        avg =",t1["Mean"],"", "\n")
  #  cat("                        max =",t1["Max."],"", "\n")
  # }
  
  ids <- sort( unique(idvar) )
  idlenmax <- t1["Max."]
  dat.descr <- c(NT = nt, N = n, t1["Min."], t1["Mean"], t1["Max."])
  
  #   print(head(Zu))
  #   print(head(Zv))
  #   print(head(Zdel))
 }
 
 # sort all
 Y <- Y[sorted,, drop = FALSE]
 X <- X[sorted,, drop = FALSE]
 
 Zu0 <- Zu0[sorted,, drop = FALSE]
 Zvi <- Zvi[sorted,, drop = FALSE]
 zdeli <- zdeli[sorted,, drop = FALSE]
 # print(head(zdeli))
 
 abbr.length <- 121
 names_x  <- abbreviate(colnames(X), abbr.length, strict = TRUE, dot = FALSE, method = "both.sides")
 names_u0 <- abbreviate(colnames(Zu0), abbr.length, strict = TRUE, dot = FALSE, method = "both.sides")
 names_vi <- abbreviate(colnames(Zvi), abbr.length, strict = TRUE, dot = FALSE, method = "both.sides")
 names_udeli <- abbreviate(colnames(zdeli), abbr.length, strict = TRUE, dot = FALSE, method = "both.sides")
 # print(names_udeli)
 max.name.length <- max(nchar(c(names_x,names_u0,names_vi,names_udeli)))
 
 names_u0 <- colnames(Zu0)
 names_udeli <- colnames(zdeli)
 
 # print(names_udeli)
 
 # print(c(dim(Zv0), length(idvar)))
 
 # print(as.vector(unlist( by(data = Zv0[,-1,drop = F], INDICES = idvar, function(x) apply(x, 2, sd)) )))
 
 if(!is.null(mean.u.0i)){
  # print(as.vector(unlist(by(data = zdeli[,-1, drop = F], INDICES = idvar, FUN = function(x) apply(x, 2, sd)))))
  if(sum(as.vector(unlist(by(data = zdeli[,-1, drop = F], INDICES = idvar, FUN = function(x) apply(x, 2, sd)))), na.rm = TRUE) > 0){
   stop("Time-varying variable(s) in 'mean.u.0i'", call. = FALSE)
  } 
 }
 zdeli <- matrix(unlist(by(zdeli, idvar, colMeans)), nrow = n, byrow = TRUE)
 colnames(zdeli) <- names_udeli
 # print(head(zdeli))
 
 
 
 if(!is.null(ln.var.u.0i)){
  # print(as.vector(unlist(by(data = Zu0[,-1, drop = F], INDICES = idvar, FUN = function(x) apply(x, 2, sd)))))
  if(sum(as.vector(unlist(by(data = Zu0[,-1, drop = F], INDICES = idvar, FUN = function(x) apply(x, 2, sd)))), na.rm = TRUE) > 0){
   stop("Time-varying variable(s) in 'ln.var.u.0i'", call. = FALSE)
  }
 }
 Zu0 <- matrix(unlist(by(Zu0, idvar, colMeans)), nrow = n, byrow = TRUE)
 colnames(Zu0) <- names_u0
 # print(head(Zu0))
 
 
 # esample <- esample[sorted,, drop = FALSE]
 # esample.nu <- esample.nu[sorted]
 
 # colnames(X)[1] <- colnames(Zvi)[1] <- colnames(Zu0)[1] <- colnames(zdeli)[1] <- "(Intercept)"
 
 
 
 tymch <- list(Formula = form1, Y = Y, X = X,
               nt = nt, n = n, Kb = k,
               idvar = idvar, timevar = timevar, ii = ii, 
               ids = ids, idlenmax = idlenmax,
               t0 = t0, itvar = itvar,
               dat.descr = dat.descr, t_i = t1, old.sorted = old.sorted,
               Zvi = Zvi, Zu0 = Zu0, Zdeli = zdeli,
               Kvi = kvi, Ku0 = ku0, Kdeli = kdeli,
               esample = esample, esample.nu = esample.nu)
 
 # cat("3\n")
 # print(head(zdeli))
 
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
 # cat("0")
 class(tymch) <- "npsf"
 return(tymch)
 
}

# Log-likelihood

.ll.panel <- function(theta, prod, coef.fixed, my.n, Ktheta, k, kv, ku, kdel, yit, zvit, zui, xit, zdeli, model, timevar, maxT, ids, idvar, t0, eff.time.invariant = FALSE, mean.u.0i.zero = FALSE) {
 
 # print(theta)

 # theta <- theta[!coef.fixed]
 # print(coef.fixed)
 s <- ifelse(prod, -1, 1)
 
 beta <- theta[1:k]
 gv <- theta[(k+1):(k+kv)]
 gu <- theta[(k+kv+1):(k+kv+ku)]
 
 eit <- as.vector( yit - xit%*%beta )
 # print(summary(as.vector(eit)))
 sigu2i <- as.vector( exp(zui%*%gu) )
 sigv2it <- as.vector( exp(zvit%*%gv) )
 
 if(!mean.u.0i.zero){
  delta <- theta[(k+kv+ku+1):(k+kv+ku+kdel)]
  mui <- as.vector( zdeli%*%delta )
 } else {
  mui <- rep(0, my.n)
 }
 # print(k+kv+ku+kdel)
 # print(kdel)
 # print(Ktheta)
 
 maxT_all <- maxT
 
 if(eff.time.invariant){
  Git <- rep(1, length(timevar))
 } else {
  if(model == "BC1992"){
   # print(theta)
   # print(c( Ktheta) )
   # print(length(theta))
   # print(theta[c( Ktheta) ])
   # eta <- theta[12 ]
   # cat("eta = ", theta[11 ],"\n")
   # print(eta)
   eta <- theta[c( Ktheta )]
   Git = exp(-eta*(timevar - maxT_all))
  } else if(model == "K1990modified"){
   eta <- theta[c( (Ktheta-1) : (Ktheta) )]
   Git = 1 + eta[1]*(timevar - maxT_all) + eta[2]*(timevar - maxT_all)^2
  } else if(model == "K1990"){
   eta <- theta[c( (Ktheta-1) : (Ktheta) )]
   Git = (1 + exp(eta[1]*timevar + eta[2]*timevar^2))^(-1)
  } else {
   stop("Unknown model\n")
  }
  # cat("eta\n")
  # print(eta)
 }
 
 # define indices
 m.end <- cumsum(t0)
 m.begin <- c(1,(m.end+1)[-my.n])

 # Log-likelihood
 llf = 0
 for(i in 1:length(ids)){
  # cat("i = ", i, "ll = ",llf,"\n")
  sample.i <- (m.begin[i]):(m.end[i])
  ei = eit[sample.i]
  sigv2i = sigv2it[sample.i]
  Gi = Git[sample.i]
  Ti = length(ei)
  lmdi = mui[i]/sigu2i[i]
  sstart = (1/sigu2i[i] + sum(Gi^2/sigv2i))
  
  llf = llf - 1/2*log(sstart) - 1/2*(mui[i]*lmdi + sum(ei^2/sigv2i)) + 1/2*((lmdi + s*sum(ei*Gi/sigv2i))^2)/sstart - Ti/2 * log(2*pi) - 1/2*sum(log(sigv2i)) - 1/2*log(sigu2i[i]) - pnorm(mui[i]/sqrt(sigu2i[i]), log.p = TRUE) + pnorm((lmdi + s*sum(ei*Gi/sigv2i))/sqrt(sstart), log.p = TRUE)
  
  # if (i == 2){
   # cat("\n")
   # print(Gi)
   # cat("sigu2i[i]\n")
   # print(sigu2i[i] )
   # print(sigv2i)
   # cat("1/2*log(sstart)\n")
   # print(1/2*log(sstart))
   # cat("1/2*(mui[i]*lmdi + sum(ei^2/sigv2i))\n")
   # print(1/2*(mui[i]*lmdi + sum(ei^2/sigv2i)))
   # cat("1/2*((lmdi + s*sum(ei*Gi/sigv2i))^2)/sstart\n")
   # print(1/2*((lmdi + s*sum(ei*Gi/sigv2i))^2)/sstart)
   # cat("1/2*sum(log(sigv2i))\n")
   # print(1/2*sum(log(sigv2i)))
   # cat("1/2*log(sigu2i[i])\n")
   # print(1/2*log(sigu2i[i]))
   # cat("pnorm(mui[i]/sqrt(sigu2i[i]), log.p = TRUE)\n")
   # print(pnorm(mui[i]/sqrt(sigu2i[i]), log.p = TRUE))
   # cat("pnorm((lmdi + s*sum(ei*Gi/sigv2i))/sqrt(sstart), log.p = TRUE)\n")
   # print(pnorm((lmdi + s*sum(ei*Gi/sigv2i))/sqrt(sstart), log.p = TRUE))
   # cat("\n")
  # }
 }
 
 return(llf)
}

.gr.hess.panel <- function(theta, prod, coef.fixed, my.n, Ktheta, k, kv, ku, kdel, yit, zvit, zui, xit, zdeli, model, timevar, maxT, ids, idvar, t0, eff.time.invariant = FALSE, mean.u.0i.zero = FALSE) {
 
 s <- ifelse(prod, -1, 1)
 beta <- theta[1:k]
 gv <- theta[(k+1):(k+kv)]
 gu <- theta[(k+kv+1):(k+kv+ku)]
 
 eit <- as.vector( yit - xit%*%beta )
 sigu2i <- as.vector( exp(zui%*%gu) )
 sigv2it <- as.vector( exp(zvit%*%gv) )
 
 if(!mean.u.0i.zero){
  delta <- theta[(k+kv+ku+1):(k+kv+ku+kdel)]
  mui <- as.vector( zdeli%*%delta )
  lmdi = mui/sigu2i
 } else {
  mui <- lmdi <- rep(0, my.n)
 }
 
 
 maxT_all <- maxT
 
 if(eff.time.invariant){
  Git <- rep(1, length(timevar))
 } else {
  if(model == "BC1992"){
   # print(theta)
   # print(c( Ktheta) )
   # print(length(theta))
   # print(theta[c( Ktheta) ])
   # eta <- theta[12 ]
   # cat("eta = ", theta[11 ],"\n")
   # print(eta)
   eta <- theta[c( Ktheta )]
   Git = exp(-eta*(timevar - maxT_all))
  } else if(model == "K1990modified"){
   eta <- theta[c( (Ktheta-1) : (Ktheta) )]
   Git = 1 + eta[1]*(timevar - maxT_all) + eta[2]*(timevar - maxT_all)^2
  } else if(model == "K1990"){
   eta <- theta[c( (Ktheta-1) : (Ktheta) )]
   Git = (1 + exp(eta[1]*timevar + eta[2]*timevar^2))^(-1)
  } else {
   stop("Unknown model\n")
  }
  # cat("eta\n")
  # print(eta)
 }
 # Keta <- length(coef.fixed) - Ktheta
 # cat.print(Keta)
 # cat.print(length(coef.fixed))
 # cat.print(Ktheta)
 # # cat("Git\n")
 # print(Git)
 
 # define indices
 m.end <- cumsum(t0)
 m.begin <- c(1,(m.end+1)[-my.n])
 
 gb = ggu = ggv = gdel = geta = 0
 
 Hb = Hbgv = Hbgu = Hbdel = Hbeta = Hgvgu = Hgvdel = Hgv = Hdel = Hdelgu = Hdeleta = Hetagu = Hgu = Hetagv = Heta = 0;
 # Hb = Hbgv = Hbgu = Hbdel = Hbeta = Hgvgu = Hgvdel = Hgv = Hdel = Hdelgu = Hdeleta = Hetagu = Hgu = Hetagv = Heta = 0;
 # if(eff.time.invariant) Hbeta <- matrix(0, nrow = Keta, ncol = k)
 
 for(i in 1:length(ids)){

  sample.i <- (m.begin[i]):(m.end[i])
  ei = eit[sample.i]
  sigv2i = sigv2it[sample.i]
  Gi = Git[sample.i]
  
  sstart = (1/sigu2i[i] + sum(Gi^2/sigv2i))
  
  xi = xit[sample.i, , drop = FALSE]; 
  zvi = zvit[sample.i, , drop = FALSE];
  Ti = length(ei)
  timevari = timevar[sample.i]
  c1 = (1/sigu2i[i] + sum(Gi^2/sigv2i))^(-1);
  c2 = ei*Gi/sigv2i; 
  c2s = sum(c2);
  c3 = t(xi)%*%(Gi/sigv2i)
  a1i = mui[i]/sqrt(sigu2i[i]); 
  a2i = (lmdi[i] + s*sum(ei*Gi/sigv2i))*sqrt(c1)
  d1 = dnorm(a1i)/pnorm(a1i); 
  d2 = dnorm(a2i)/pnorm(a2i)
  d3 = (1 - d2*(a2i + d2)) # this is A
  d4 = d1*(a1i + d1)
  
  if(!eff.time.invariant){
   if(model == "BC1992"){
    tBC1 = (timevari - maxT_all[sample.i]);
   } else if (model == "K1990modified"){
    tBC2 = cbind((timevari - maxT_all[sample.i]), (timevari - maxT_all[sample.i])^2);
   } else if (model == "K1990"){
    Gk = (Gi^(-1)-1)
    tSK = cbind(timevari, timevari^2)
   } else {
    stop("Unknown model\n")
   }
  }
  
  # if (i == 2){
  #  cat("\n")
  #  cat("Gi\n")
  #  print(Gi)
  #  cat(" c1*c3%*%t(c3)*d3\n")
  #  print( c1*c3%*%t(c3)*d3)
  #  cat.print(c1)
  #  cat.print(c3)
  #  cat.print(d3)
  #  cat.print(a2i)
  #  cat.print(d2)
  #  cat("head(xi)\n")
  #  print(head(xi))
  #  cat("sigu2i[i]\n")
  #  print(sigu2i[i] )
  #  print(sigv2i)
  #  cat("1/2*log(sstart)\n")
  #  print(1/2*log(sstart))
  #  cat("1/2*(mui[i]*lmdi + sum(ei^2/sigv2i))\n")
  #  print(1/2*(mui[i]*lmdi + sum(ei^2/sigv2i)))
  #  cat("1/2*((lmdi + s*sum(ei*Gi/sigv2i))^2)/sstart\n")
  #  print(1/2*((lmdi + s*sum(ei*Gi/sigv2i))^2)/sstart)
  #  cat("1/2*sum(log(sigv2i))\n")
  #  print(1/2*sum(log(sigv2i)))
  #  cat("1/2*log(sigu2i[i])\n")
  #  print(1/2*log(sigu2i[i]))
  #  cat("pnorm(mui[i]/sqrt(sigu2i[i]), log.p = TRUE)\n")
  #  print(pnorm(mui[i]/sqrt(sigu2i[i]), log.p = TRUE))
  #  cat("pnorm((lmdi + s*sum(ei*Gi/sigv2i))/sqrt(sstart), log.p = TRUE)\n")
  #  print(pnorm((lmdi + s*sum(ei*Gi/sigv2i))/sqrt(sstart), log.p = TRUE))
  #  cat("\n")
  # }
  
  # Gradient
  
  gb = gb + t(xi)%*%((ei/sigv2i)- s*(Gi/sigv2i)*(lmdi[i]*c1 + s*c2s*c1 + sqrt(c1)*d2))
  ggu = ggu + zui[i,,drop=FALSE]/sigu2i[i]*(0.5*c1 + mui[i]*(0.5*mui[i] + 0.5*d1*sqrt(sigu2i[i]) - sqrt(c1)*d2) + a2i*(c1*a2i/2 - sqrt(c1)*mui[i] + c1*d2/2)) - zui[i,,drop=FALSE]*0.5
  ggv = ggv + t(zvi)%*%(0.5*c1*(Gi^2/sigv2i) + 0.5*(ei^2/sigv2i) - 0.5*rep.int(1, Ti) + a2i*sqrt(c1)*((a2i + d2)*0.5*sqrt(c1)*(Gi^2/sigv2i) - s*(ei*Gi/sigv2i)) - s*sqrt(c1)*d2*(ei*Gi/sigv2i))
  gdel = gdel + zdeli[i,,drop=FALSE]/sigu2i[i] *(mui[i]*(c1/sigu2i[i] - 1) + s*c2s*c1 - d1 *sqrt(sigu2i[i]) + d2*sqrt(c1))
  
  # Hessian
  
  Hb = Hb - t(xi)%*%sweep(xi, 1, sigv2i, "/") + c1*c3%*%t(c3)*d3
  
  Hbgu = Hbgu + s*c3%*%zui[i,,drop=FALSE]*c1*(lmdi[i]*d3 - c1/sigu2i[i]*(lmdi[i] + s*sum(c2)) + sqrt(c1)/(2*sigu2i[i])*d2*(a2i*(a2i + d2) - 1))
  Hbdel = Hbdel - s*c3%*%(zdeli[i,]/sigu2i[i]*d3*c1)
  Hbgv = Hbgv + t(xi)%*%(-sweep(zvi, 1,ei/sigv2i, "*") + s*sweep(zvi, 1, Gi/sigv2i, "*")*(sqrt(c1)*(a2i + d2))) - s*(t(xi)%*%(Gi/sigv2i))%*%t(t(zvi)%*%(-s*c1*c2*d3 + (c1^1.5*(Gi^2/sigv2i))*(a2i + 0.5*d2*(1 - a2i*(a2i + d2)))))
  Hgvgu = Hgvgu + t(zvi)%*%(c1*s*c2*(mui[i]*d3 - sqrt(c1)*a2i/2*(1+d3)) + c1^1.5*a2i*0.5*(Gi^2/sigv2i)*(sqrt(c1)*a2i*0.5*(3+d3) - mui[i]*(1+d3)) + 0.5*c1^1.5*(Gi^2/sigv2i)*(sqrt(c1) + d2*(1.5*sqrt(c1)*a2i - mui[i])) - 0.5*c1^1.5*s*d2*c2)%*%zui[i,,drop=FALSE]/sigu2i[i]
  Hgvdel = Hgvdel + t(zvi)%*%(-s*c1*c2*d3 + (Gi^2/sigv2i)*(0.5*c1^1.5*(a2i*(1+d3) + d2)))%*%zdeli[i,]/sigu2i[i]
  Hgv = Hgv +  t(zvi)%*%sweep(zvi,1,Gi^2/sigv2i, "*")*0.5*c1*(-d2*a2i - a2i^2 -1) + t(zvi)%*%(Gi^2/sigv2i)%*%t(t(zvi)%*%(Gi^2/sigv2i))*0.5*c1^2*(1 + 0.5*a2i^2*(3+d3) + d2*1.5*a2i) - 0.5*t(zvi)%*%sweep(zvi,1,ei^2/sigv2i, "*") + t(zvi)%*%sweep(zvi,1,c2, "*")*s*sqrt(c1)*(d2 + a2i) - (d2+a2i*(1+d3))*0.5*c1^1.5*s*(t(zvi)%*%(c2)%*%t(t(zvi)%*%(Gi^2/sigv2i)) + t(zvi)%*%(Gi^2/sigv2i)%*%t(t(zvi)%*%(c2))) + t(zvi)%*%(c2)%*%t(t(zvi)%*%(c2))*c1*d3
  Hdel = Hdel + t(zdeli[i,,drop=FALSE]/sigu2i[i] *( - 1 + d3*c1/sigu2i[i] + d4))%*%zdeli[i,,drop=FALSE]
  Hdelgu =  Hdelgu + t(zdeli[i,,drop=FALSE])%*%zui[i,,drop=FALSE]/sigu2i[i]*(c1^1.5/sigu2i[i]*0.5*a2i*(1+d3) - mui[i]*c1/sigu2i[i]*(1+d3) + d2*sqrt(c1)*(0.5*c1/sigu2i[i] - 1) - 0.5*d1*(mui[i]*(a1i + d1) - sqrt(sigu2i[i])) + mui[i] - s*c1*sum(c2))
  Hgu = Hgu + (c1/sigu2i[i]*(0.5*c1 + (sqrt(c1)*a2i - mui[i])^2) - 0.5*(c1 + (sqrt(c1)*a2i - mui[i])^2) + d1*mui[i]/4*(mui[i]*(a1i + d1) - sqrt(sigu2i[i])) + d2*(-d2*2*c1/sigu2i[i]*(mui[i]/sqrt(2) - sqrt(c1/8)*a2i)^2 - c1*a2i/sigu2i[i]*(mui[i] - sqrt(c1)*a2i/2)^2 + sqrt(c1)*mui[i]*(1 - c1/sigu2i[i]) - c1*a2i/2*(1 - 1.5*c1/sigu2i[i])))*t(zui[i,,drop=FALSE])%*%zui[i,,drop=FALSE]/sigu2i[i]
  
  if(!eff.time.invariant){
   if(model == "BC1992"){
    geta = geta + c1*sum(Gi^2*tBC1/sigv2i) - s*sqrt(c1)*d2*sum(ei*Gi*tBC1/sigv2i) + sqrt(c1)*a2i*(sqrt(c1)*(a2i + d2)*sum(Gi^2*tBC1/sigv2i) - s*sum(ei*Gi*tBC1/sigv2i))
    Hbeta = Hbeta + s*t(xi)%*%(Gi*tBC1/sigv2i)*sqrt(c1)*(a2i + d2) - s*t(xi)%*%(Gi/sigv2i)*c1*(-s*sum(c2*tBC1)*d3 + sqrt(c1)*a2i*(1 + d3)*sum(Gi^2*tBC1/sigv2i) + sqrt(c1)*d2*sum(Gi^2*tBC1/sigv2i))
    Hdeleta = Hdeleta + t((sum(Gi^2*tBC1/sigv2i)*c1^1.5*(a2i*(1+d3) + d2)-s*c1*sum(c2*tBC1)*d3)*zdeli[i,,drop=FALSE]/sigu2i[i])
    Hetagu = Hetagu + t((sum(Gi^2*tBC1/sigv2i)*c1^1.5/sigu2i[i]*(sqrt(c1) -a2i*mui[i]*(1+d3) + sqrt(c1)*a2i^2*0.5*(3+d3) - d2*(mui[i] - 1.5*a2i*sqrt(c1))) + s*c1/sigu2i[i]*sum(c2*tBC1)*(mui[i]*d3 - d2*0.5*sqrt(c1) - 0.5*sqrt(c1)*a2i*(1+d3)))*zui[i,,drop=FALSE])
    Hetagv = Hetagv - t(zvi)%*%(Gi^2*tBC1/sigv2i)*c1*(1 + a2i^2 + a2i*d2) + sum(Gi^2*tBC1/sigv2i)*t(zvi)%*%(Gi^2/sigv2i)*c1^2*(1 + 0.5*a2i^2*(3 + d3) + 1.5*a2i*d2) + t(zvi)%*%(c2*tBC1)*s*sqrt(c1)*(d2 + a2i) + sum(c2*tBC1)*t(zvi)%*%(c2)*c1*d3 - (a2i*(1+d3) + d2)*s*c1^1.5*(0.5*sum(c2*tBC1)*t(zvi)%*%(Gi^2/sigv2i) + sum(Gi^2*tBC1/sigv2i)*t(zvi)%*%(c2))
    Heta = Heta - 2*c1*sum(Gi^2*tBC1^2/sigv2i)*(1 + a2i^2 + d2*a2i) + c1^2*sum(Gi^2*tBC1/sigv2i)^2*(2 + a2i^2*(3 + d3) + 3*d2*a2i) + sum(c2*tBC1^2)*s*sqrt(c1)*(d2 + a2i) - 2*s*c1^1.5*sum(c2*tBC1)*sum(Gi^2*tBC1/sigv2i)*(d2 + a2i*(1 + d3)) + sum(c2*tBC1)^2*c1*d3
    
   } else if(model == "K1990modified"){
    geta =  geta - c1*t(Gi/sigv2i)%*%tBC2 + 
     s*sqrt(c1)*d2*t(ei/sigv2i)%*%tBC2 - 
     sqrt(c1)*a2i*(sqrt(c1)*(a2i + d2)*t(Gi/sigv2i)%*%tBC2 - s*t(ei/sigv2i)%*%tBC2)
    Hbeta = Hbeta - s*t(xi)%*%sweep(tBC2, 1, sigv2i, "/")*sqrt(c1)*(a2i + d2) - s*(t(xi)%*%(Gi/sigv2i))%*%t(t(tBC2)%*%(c1*(s*(ei/sigv2i)*d3 - sqrt(c1)*a2i*(1 + d3)*(Gi/sigv2i) - sqrt(c1)*d2*(Gi/sigv2i))))
    Hdeleta = Hdeleta + t(zdeli[i,,drop=FALSE])%*%(t(t(tBC2)%*%(s*c1*(ei/sigv2i)*d3 - c1^1.5*(a2i*(1+d3) + d2)*(Gi/sigv2i)))/sigu2i[i])
    Hetagu = Hetagu + t(zui[i,,drop=FALSE])%*%(- t(Gi/sigv2i)%*%tBC2*c1^1.5/sigu2i[i]*(sqrt(c1) -a2i*mui[i]*(1+d3) + sqrt(c1)*a2i^2*0.5*(3+d3) - d2*(mui[i] - 1.5*a2i*sqrt(c1))) - t(ei/sigv2i)%*%tBC2*s*c1/sigu2i[i]*(mui[i]*d3 - d2*0.5*sqrt(c1) - 0.5*sqrt(c1)*a2i*(1+d3)))
    Hetagv = Hetagv + t(zvi)%*%sweep(tBC2, 1, Gi/sigv2i, "*")*c1*(1 + a2i^2 + a2i*d2) - (t(zvi)%*%(Gi^2/sigv2i))%*%(t(Gi/sigv2i)%*%tBC2)*c1^2*(1 + 0.5*a2i^2*(3 + d3) + 1.5*a2i*d2) - t(zvi)%*%sweep(tBC2, 1, ei/sigv2i, "*")*s*sqrt(c1)*(d2 + a2i) - (t(zvi)%*%(c2))%*%(t(ei/sigv2i)%*%tBC2)*c1*d3 + (a2i*(1+d3) + d2)*s*c1^1.5*(0.5*(t(zvi)%*%(Gi^2/sigv2i))%*%(t(ei/sigv2i)%*%tBC2) + (t(zvi)%*%(c2))%*%(t(Gi/sigv2i)%*%tBC2))
    Heta = Heta - c1*(1 + a2i^2 + a2i*d2)*t(tBC2)%*%sweep(tBC2, 1, sigv2i, "/") + c1^2*(2 + a2i^2*(3 + d3) + d2*a2i*3)*t(t(Gi/sigv2i)%*%tBC2)%*%(t(Gi/sigv2i)%*%tBC2) - c1^1.5*s*(d2+a2i*(1+d3))*t(t(ei/sigv2i)%*%tBC2)%*%(t(Gi/sigv2i)%*%tBC2) + c1*d3*t(t(ei/sigv2i)%*%tBC2)%*%(t(ei/sigv2i)%*%tBC2) - c1^1.5*s*(d2 + a2i*(1+d3))*t(t(Gi/sigv2i)%*%tBC2)%*%(t(ei/sigv2i)%*%tBC2)
    
   } else if(model == "K1990"){
    geta =  geta + c1*t(Gi^3*(Gi^(-1)-1)/sigv2i)%*%tSK -  s*sqrt(c1)*d2*t(ei*Gi^2*(Gi^(-1)-1)/sigv2i)%*%tSK + sqrt(c1)*a2i*(sqrt(c1)*(a2i + d2)*t(Gi^3*(Gi^(-1)-1)/sigv2i)%*%tSK - s*t(ei*Gi^2*(Gi^(-1)-1)/sigv2i)%*%tSK)
    Hbeta = Hbeta +  s*t(xi)%*%sweep(tSK, 1, Gi^2*Gk/sigv2i, "*")*sqrt(c1)*(a2i + d2) - s*(t(xi)%*%(Gi/sigv2i))%*%t(t(tSK)%*%(c1*(-s*(ei*Gi^2*Gk/sigv2i)*d3 + sqrt(c1)*a2i*(1 + d3)*(Gi^3*Gk/sigv2i) + sqrt(c1)*d2*(Gi^3*Gk/sigv2i))))
    
    Hdeleta = Hdeleta + t(zdeli[i,,drop=FALSE])%*%(s*t(t(tSK)%*%(-s*c1*(c2*Gi*Gk)*d3 + c1^1.5*(a2i*(1+d3) + d2)*(Gi^3*Gk/sigv2i)))/sigu2i[i])
    
    Hetagu = Hetagu + t(zui[i,,drop=FALSE])%*%(t(Gi^3*Gk/sigv2i)%*%tSK*c1^1.5/sigu2i[i]*(sqrt(c1) -a2i*mui[i]*(1+d3) + sqrt(c1)*a2i^2*0.5*(3+d3) - d2*(mui[i] - 1.5*a2i*sqrt(c1))) + t(c2*Gi*Gk)%*%tSK*s*c1/sigu2i[i]*(mui[i]*d3 - d2*0.5*sqrt(c1) - 0.5*sqrt(c1)*a2i*(1+d3)))
    
    Hetagv = Hetagv - t(zvi)%*%sweep(tSK, 1, Gi^3*Gk/sigv2i, "*")*c1*(1 + a2i^2 + a2i*d2) + (t(zvi)%*%(Gi^2/sigv2i))%*%(t(Gi^3*Gk/sigv2i)%*%tSK)*c1^2*(1 + 0.5*a2i^2*(3 + d3) + 1.5*a2i*d2) + t(zvi)%*%sweep(tSK, 1, c2*Gi*Gk, "*")*s*sqrt(c1)*(d2 + a2i) + (t(zvi)%*%(c2))%*%(t(c2*Gi*Gk)%*%tSK)*c1*d3 - (a2i*(1+d3) + d2)*s*c1^1.5*(0.5*(t(zvi)%*%(Gi^2/sigv2i))%*%(t(c2*Gi*Gk)%*%tSK) + (t(zvi)%*%(c2))%*%(t(Gi^3*Gk/sigv2i)%*%tSK))
    
    Heta = Heta + c1*(1 + a2i^2 + a2i*d2)*t(tSK)%*%sweep(tSK, 1, Gi^3*Gk*(1- 3*Gi*Gk)/sigv2i, "*") + c1^2*(2 + a2i^2*(3 + d3) + d2*a2i*3)*t(t(Gi^3*Gk/sigv2i)%*%tSK)%*%(t(Gi^3*Gk/sigv2i)%*%tSK) - s*sqrt(c1)*(a2i + d2)*t(tSK)%*%sweep(tSK, 1, ei*Gi^2*Gk*(1 - 2*Gi*Gk)/sigv2i, "*") - c1^1.5*s*(d2 + a2i*(1 + d3))*(t(t(ei*Gi^2*Gk/sigv2i)%*%tSK)%*%(t(Gi^3*Gk/sigv2i)%*%tSK) + t(t(Gi^3*Gk/sigv2i)%*%tSK)%*%(t(ei*Gi^2*Gk/sigv2i)%*%tSK)) + c1*d3*t(t(ei*Gi^2*Gk/sigv2i)%*%tSK)%*%(t(ei*Gi^2*Gk/sigv2i)%*%tSK)
    # print(geta)
   }
   else {
    stop("Unknown model\n")
   }
  }
  
  # if(i == 2){
  #  cat("Hb\n")
  #  print(Hb)
  # }

 }
 
 
 
 if(mean.u.0i.zero){
  
  if(eff.time.invariant){
   H <- cbind(rbind(Hb,	t(Hbgv),	t(Hbgu)),
              rbind(Hbgv, 	Hgv,		t(Hgvgu)),
              rbind(Hbgu, 	Hgvgu,		Hgu))
  } else {
   H <- cbind(rbind(Hb,	t(Hbgv),	t(Hbgu),	t(Hbeta)),
              rbind(Hbgv, 	Hgv,		t(Hgvgu),	t(Hetagv)),
              rbind(Hbgu, 	Hgvgu,		Hgu,		 	t(Hetagu)),
              rbind(Hbeta,	Hetagv,		Hetagu,		Heta ))
  }
 } else {
  
  if(eff.time.invariant){
   H <- cbind(rbind(Hb,	t(Hbgv),	t(Hbgu),		t(Hbdel)),
              rbind(Hbgv, 	Hgv,		t(Hgvgu),		t(Hgvdel)),
              rbind(Hbgu, 	Hgvgu,		Hgu,			Hdelgu),
              rbind(Hbdel,	Hgvdel,		t(Hdelgu),		Hdel))
   
  } else {
   H <- cbind(rbind(Hb,	t(Hbgv),	t(Hbgu),		t(Hbdel),	t(Hbeta)),
              rbind(Hbgv, 	Hgv,		t(Hgvgu),		t(Hgvdel),	t(Hetagv)),
              rbind(Hbgu, 	Hgvgu,		Hgu,			Hdelgu, 	t(Hetagu)),
              rbind(Hbdel,	Hgvdel,		t(Hdelgu),		Hdel,		t(Hdeleta) ),
              rbind(Hbeta,	Hetagv,		Hetagu,			Hdeleta,	Heta ))
  }
  
 }
 
 # if(eff.time.invariant){
 #  grad <- rbind(gb, as.matrix(ggv),t(as.matrix(ggu)), t(as.matrix(gdel)))
 #  H <- cbind(rbind(Hb,	t(Hbgv),	t(Hbgu),		t(Hbdel)),
 #             rbind(Hbgv, 	Hgv,		t(Hgvgu),		t(Hgvdel)),
 #             rbind(Hbgu, 	Hgvgu,		Hgu,			Hdelgu),
 #             rbind(Hbdel,	Hgvdel,		t(Hdelgu),		Hdel))
 # } else {
 #  grad <- rbind(gb, as.matrix(ggv),t(as.matrix(ggu)), t(as.matrix(gdel)), t(geta))
 #  H <- cbind(rbind(Hb,	t(Hbgv),	t(Hbgu),		t(Hbdel),	t(Hbeta)),
 #             rbind(Hbgv, 	Hgv,		t(Hgvgu),		t(Hgvdel),	t(Hetagv)),
 #             rbind(Hbgu, 	Hgvgu,		Hgu,			Hdelgu, 	t(Hetagu)),
 #             rbind(Hbdel,	Hgvdel,		t(Hdelgu),		Hdel,		t(Hdeleta) ),
 #             rbind(Hbeta,	Hetagv,		Hetagu,			Hdeleta,	Heta ))
 # }
 
 grad <- rbind(gb, as.matrix(ggv),t(as.matrix(ggu)), t(as.matrix(gdel)), t(geta))
 grad <- as.vector(grad)
 
 # cat.print(Hbeta)
 # 
 # H <- cbind(rbind(Hb,	t(Hbgv),	t(Hbgu),		t(Hbdel),	t(Hbeta)), 
 #            rbind(Hbgv, 	Hgv,		t(Hgvgu),		t(Hgvdel),	t(Hetagv)),
 #            rbind(Hbgu, 	Hgvgu,		Hgu,			Hdelgu, 	t(Hetagu)), 
 #            rbind(Hbdel,	Hgvdel,		t(Hdelgu),		Hdel,		t(Hdeleta) ),
 #            rbind(Hbeta,	Hetagv,		Hetagu,			Hdeleta,	Heta ))
 
 
 grad <- grad[!coef.fixed]
 # print(H)
 # H <- H[!coef.fixed,!coef.fixed]
 # print(H)
 
 # if(mean.u.0i.zero){
 #  grad <- grad[-(k+kv+ku+kdel)]
 #  H <- H[-(k+kv+ku+kdel),-(k+kv+ku+kdel)]
 # } 
 
 # print(length(grad))
 # print(dim(H))
 # print(length(theta))
 
 
 names(grad) <- colnames(H) <- rownames(H) <- names(theta)
 return(list(grad = grad, hessian1 = H))
 
}

# Gradient

.gr.panel <- function(theta, prod, coef.fixed, my.n, Ktheta, k, kv, ku, kdel, yit, zvit, zui, xit, zdeli, model, timevar, maxT, ids, idvar, t0, eff.time.invariant = FALSE, mean.u.0i.zero = FALSE) {
 s <- ifelse(prod, -1, 1)
 beta <- theta[1:k]
 gv <- theta[(k+1):(k+kv)]
 gu <- theta[(k+kv+1):(k+kv+ku)]
 
 eit <- as.vector( yit - xit%*%beta )
 sigu2i <- as.vector( exp(zui%*%gu) )
 sigv2it <- as.vector( exp(zvit%*%gv) )
 
 if(!mean.u.0i.zero){
  delta <- theta[(k+kv+ku+1):(k+kv+ku+kdel)]
  mui <- as.vector( zdeli%*%delta )
  lmdi = mui/sigu2i
 } else {
  mui <- lmdi <- rep(0, my.n)
 }
 
 maxT_all <- maxT
 
 if(eff.time.invariant){
  Git <- rep(1, length(timevar))
 } else {
  if(model == "BC1992"){
   # print(theta)
   # print(c( Ktheta) )
   # print(length(theta))
   # print(theta[c( Ktheta) ])
   # eta <- theta[12 ]
   # cat("eta = ", theta[11 ],"\n")
   # print(eta)
   eta <- theta[c( Ktheta )]
   Git = exp(-eta*(timevar - maxT_all))
  } else if(model == "K1990modified"){
   eta <- theta[c( (Ktheta-1) : (Ktheta) )]
   Git = 1 + eta[1]*(timevar - maxT_all) + eta[2]*(timevar - maxT_all)^2
  } else if(model == "K1990"){
   eta <- theta[c( (Ktheta-1) : (Ktheta) )]
   Git = (1 + exp(eta[1]*timevar + eta[2]*timevar^2))^(-1)
  } else {
   stop("Unknown model\n")
  }
  # cat("eta\n")
  # print(eta)
 }
 
 # define indices
 m.end <- cumsum(t0)
 m.begin <- c(1,(m.end+1)[-my.n])
 
 # Gradient
 gb = 0; ggu = 0; ggv = 0; gdel = 0; geta = 0
 for(i in 1:length(ids)){
  # maxT <- t0[i]
  sample.i <- (m.begin[i]):(m.end[i])
  ei = eit[sample.i]
  sigv2i = sigv2it[sample.i]
  Gi = Git[sample.i]
  
  xi = xit[sample.i, , drop = FALSE]; 
  zvi = zvit[sample.i, , drop = FALSE];
  Ti = length(ei)
  timevari = timevar[sample.i]
  c1 = (1/sigu2i[i] + sum(Gi^2/sigv2i))^(-1); 
  c2 = sum(ei*Gi/sigv2i);
  a1i = mui[i]/sqrt(sigu2i[i]); 
  a2i = (lmdi[i] + s*sum(ei*Gi/sigv2i))*sqrt(c1)
  d1 = dnorm(a1i)/pnorm(a1i); 
  d2 = dnorm(a2i)/pnorm(a2i)
  
  if(!eff.time.invariant){
   if(model == "BC1992"){
    tBC1 = (timevari - maxT_all[sample.i]);
   } else if (model == "K1990modified"){
    tBC2 = cbind((timevari - maxT_all[sample.i]), (timevari - maxT_all[sample.i])^2);
   } else if (model == "K1990"){
    tSK = cbind(timevari, timevari^2)
   } else {
    stop("Unknown model\n")
   }
  }
  
  gb = gb + t(xi)%*%((ei/sigv2i)- s*(Gi/sigv2i)*(lmdi[i]*c1 + s*c2*c1 + sqrt(c1)*d2))
  ggu = ggu + zui[i,,drop=FALSE]/sigu2i[i]*(0.5*c1 + mui[i]*(0.5*mui[i] + 0.5*d1*sqrt(sigu2i[i]) - sqrt(c1)*d2) + a2i*(c1*a2i/2 - sqrt(c1)*mui[i] + c1*d2/2)) - zui[i,,drop=FALSE]*0.5
  ggv = ggv + t(zvi)%*%(0.5*c1*(Gi^2/sigv2i) + 0.5*(ei^2/sigv2i) - 0.5*rep.int(1, Ti) + a2i*sqrt(c1)*((a2i + d2)*0.5*sqrt(c1)*(Gi^2/sigv2i) - s*(ei*Gi/sigv2i)) - s*sqrt(c1)*d2*(ei*Gi/sigv2i))
  gdel = gdel + zdeli[i,,drop=FALSE]/sigu2i[i] *(mui[i]*(c1/sigu2i[i] - 1) + s*c2*c1 - d1 *sqrt(sigu2i[i]) + d2*sqrt(c1))
  
  if(!eff.time.invariant){
   if(model == "BC1992"){
    geta = geta + c1*sum(Gi^2*tBC1/sigv2i) - s*sqrt(c1)*d2*sum(ei*Gi*tBC1/sigv2i) + sqrt(c1)*a2i*(sqrt(c1)*(a2i + d2)*sum(Gi^2*tBC1/sigv2i) - s*sum(ei*Gi*tBC1/sigv2i))
   } else if(model == "K1990modified"){
    geta =  geta - c1*t(Gi/sigv2i)%*%tBC2 + 
     s*sqrt(c1)*d2*t(ei/sigv2i)%*%tBC2 - 
     sqrt(c1)*a2i*(sqrt(c1)*(a2i + d2)*t(Gi/sigv2i)%*%tBC2 - s*t(ei/sigv2i)%*%tBC2)
   } else if(model == "K1990"){
    geta =  geta + c1*t(Gi^3*(Gi^(-1)-1)/sigv2i)%*%tSK -  s*sqrt(c1)*d2*t(ei*Gi^2*(Gi^(-1)-1)/sigv2i)%*%tSK + sqrt(c1)*a2i*(sqrt(c1)*(a2i + d2)*t(Gi^3*(Gi^(-1)-1)/sigv2i)%*%tSK - s*t(ei*Gi^2*(Gi^(-1)-1)/sigv2i)%*%tSK)
    # print(geta)
   }
   else {
    stop("Unknown model\n")
   }
  }
 }
 grad <- rbind(gb, as.matrix(ggv),t(as.matrix(ggu)), t(as.matrix(gdel)), t(geta))
 grad <- grad[!coef.fixed]
 grad <- as.vector(grad)
 names(grad) <- names(theta)
 return(grad)
}

# .gr.panel <- function(theta, prod, coef.fixed, my.n, Ktheta, k, kv, ku, kdel, yit, zvit, zui, xit, zdeli, model, timevar, maxT, ids, idvar, t0, eff.time.invariant = FALSE, mean.u.0i.zero = FALSE){
#  numDeriv::grad(func = .ll.panel, x = theta, prod = prod, my.n = my.n, k = k, kv = kv, ku = ku, kdel = kdel, yit = yit, zvit = zvit, zui = zui, xit = xit, zdeli = zdeli, eff.time.invariant = eff.time.invariant, model = model, timevar = timevar, maxT = maxT, ids = ids, idvar = idvar, t0 = t0)
# }

# Hessian

.hess.panel <- function(theta, prod, coef.fixed, my.n, Ktheta, k, kv, ku, kdel, yit, zvit, zui, xit, zdeli, model, timevar, maxT, ids, idvar, t0, eff.time.invariant = FALSE, mean.u.0i.zero = FALSE) {
 s <- ifelse(prod, -1, 1)
 beta <- theta[1:k]
 gv <- theta[(k+1):(k+kv)]
 gu <- theta[(k+kv+1):(k+kv+ku)]
 
 eit <- as.vector( yit - xit%*%beta )
 sigu2i <- as.vector( exp(zui%*%gu) )
 sigv2it <- as.vector( exp(zvit%*%gv) )
 
 if(!mean.u.0i.zero){
  delta <- theta[(k+kv+ku+1):(k+kv+ku+kdel)]
  mui <- as.vector( zdeli%*%delta )
  lmdi = mui/sigu2i
 } else {
  mui <- lmdi <- rep(0, my.n)
 }

 maxT_all <- maxT
 
 if(eff.time.invariant){
  Git <- rep(1, length(timevar))
 } else {
  if(model == "BC1992"){
   # print(theta)
   # print(c( Ktheta) )
   # print(length(theta))
   # print(theta[c( Ktheta) ])
   # eta <- theta[12 ]
   # cat("eta = ", theta[11 ],"\n")
   # print(eta)
   eta <- theta[c( Ktheta )]
   Git = exp(-eta*(timevar - maxT_all))
  } else if(model == "K1990modified"){
   eta <- theta[c( (Ktheta-1) : (Ktheta) )]
   Git = 1 + eta[1]*(timevar - maxT_all) + eta[2]*(timevar - maxT_all)^2
  } else if(model == "K1990"){
   eta <- theta[c( (Ktheta-1) : (Ktheta) )]
   Git = (1 + exp(eta[1]*timevar + eta[2]*timevar^2))^(-1)
  } else {
   stop("Unknown model\n")
  }
  # cat("eta\n")
  # print(eta)
 }
 
 # define indices
 m.end <- cumsum(t0)
 m.begin <- c(1,(m.end+1)[-my.n])
 
 # Hessian
 Hb = Hbgv = Hbgu = Hbdel = Hbeta = Hgvgu = Hgvdel = Hgv = Hdel = Hdelgu = Hdeleta = Hetagu = Hgu = Hetagv = Heta = 0;
 
 for(i in 1:length(ids)){
  sample.i <- (m.begin[i]):(m.end[i])
  # maxT <- maxT_all[sample.i]
  xi = xit[sample.i, , drop=FALSE]; 
  zvi = zvit[sample.i, , drop=FALSE];
  ei = eit[sample.i]; 
  sigv2i = sigv2it[sample.i]
  Gi = Git[sample.i]; 
  Ti = length(ei); 
  timevari = timevar[sample.i];
  
  if(!eff.time.invariant){
   if(model == "BC1992"){
    tBC1 = (timevari - maxT_all[sample.i]);
   } else if (model == "K1990modified"){
    tBC2 = cbind((timevari - maxT_all[sample.i]), (timevari - maxT_all[sample.i])^2);
   } else if (model == "K1990"){
    Gk = (Gi^(-1)-1)
    tSK = cbind(timevari, timevari^2)
   } else {
    stop("Unknown model\n")
   }
  }

  c1 = (1/sigu2i[i] + sum(Gi^2/sigv2i))^(-1); 
  c2 = ei*Gi/sigv2i; 
  c3 = t(xi)%*%(Gi/sigv2i)
  a1i = mui[i]/sqrt(sigu2i[i]); 
  a2i = (lmdi[i] + s*sum(ei*Gi/sigv2i))*sqrt(c1)
  d1 = dnorm(a1i)/pnorm(a1i); 
  d2 = dnorm(a2i)/pnorm(a2i); 
  d3 = (1 - d2*(a2i + d2)); # this is A
  d4 = d1*(a1i + d1)
  
  Hb = Hb - t(xi)%*%sweep(xi, 1, sigv2i, "/") + c1*c3%*%t(c3)*d3
  
  Hbgu = Hbgu + s*c3%*%zui[i,,drop=FALSE]*c1*(lmdi[i]*d3 - c1/sigu2i[i]*(lmdi[i] + s*sum(c2)) + sqrt(c1)/(2*sigu2i[i])*d2*(a2i*(a2i + d2) - 1))
  Hbdel = Hbdel - s*c3%*%(zdeli[i,]/sigu2i[i]*d3*c1)
  Hbgv = Hbgv + t(xi)%*%(-sweep(zvi, 1,ei/sigv2i, "*") + s*sweep(zvi, 1, Gi/sigv2i, "*")*(sqrt(c1)*(a2i + d2))) - s*(t(xi)%*%(Gi/sigv2i))%*%t(t(zvi)%*%(-s*c1*c2*d3 + (c1^1.5*(Gi^2/sigv2i))*(a2i + 0.5*d2*(1 - a2i*(a2i + d2)))))
  Hgvgu = Hgvgu + t(zvi)%*%(c1*s*c2*(mui[i]*d3 - sqrt(c1)*a2i/2*(1+d3)) + c1^1.5*a2i*0.5*(Gi^2/sigv2i)*(sqrt(c1)*a2i*0.5*(3+d3) - mui[i]*(1+d3)) + 0.5*c1^1.5*(Gi^2/sigv2i)*(sqrt(c1) + d2*(1.5*sqrt(c1)*a2i - mui[i])) - 0.5*c1^1.5*s*d2*c2)%*%zui[i,,drop=FALSE]/sigu2i[i]
  Hgvdel = Hgvdel + t(zvi)%*%(-s*c1*c2*d3 + (Gi^2/sigv2i)*(0.5*c1^1.5*(a2i*(1+d3) + d2)))%*%zdeli[i,]/sigu2i[i]
  Hgv = Hgv +  t(zvi)%*%sweep(zvi,1,Gi^2/sigv2i, "*")*0.5*c1*(-d2*a2i - a2i^2 -1) + t(zvi)%*%(Gi^2/sigv2i)%*%t(t(zvi)%*%(Gi^2/sigv2i))*0.5*c1^2*(1 + 0.5*a2i^2*(3+d3) + d2*1.5*a2i) - 0.5*t(zvi)%*%sweep(zvi,1,ei^2/sigv2i, "*") + t(zvi)%*%sweep(zvi,1,c2, "*")*s*sqrt(c1)*(d2 + a2i) - (d2+a2i*(1+d3))*0.5*c1^1.5*s*(t(zvi)%*%(c2)%*%t(t(zvi)%*%(Gi^2/sigv2i)) + t(zvi)%*%(Gi^2/sigv2i)%*%t(t(zvi)%*%(c2))) + t(zvi)%*%(c2)%*%t(t(zvi)%*%(c2))*c1*d3
  Hdel = Hdel + t(zdeli[i,,drop=FALSE]/sigu2i[i] *( - 1 + d3*c1/sigu2i[i] + d4))%*%zdeli[i,,drop=FALSE]
  Hdelgu =  Hdelgu + t(zdeli[i,,drop=FALSE])%*%zui[i,,drop=FALSE]/sigu2i[i]*(c1^1.5/sigu2i[i]*0.5*a2i*(1+d3) - mui[i]*c1/sigu2i[i]*(1+d3) + d2*sqrt(c1)*(0.5*c1/sigu2i[i] - 1) - 0.5*d1*(mui[i]*(a1i + d1) - sqrt(sigu2i[i])) + mui[i] - s*c1*sum(c2))
  Hgu = Hgu + (c1/sigu2i[i]*(0.5*c1 + (sqrt(c1)*a2i - mui[i])^2) - 0.5*(c1 + (sqrt(c1)*a2i - mui[i])^2) + d1*mui[i]/4*(mui[i]*(a1i + d1) - sqrt(sigu2i[i])) + d2*(-d2*2*c1/sigu2i[i]*(mui[i]/sqrt(2) - sqrt(c1/8)*a2i)^2 - c1*a2i/sigu2i[i]*(mui[i] - sqrt(c1)*a2i/2)^2 + sqrt(c1)*mui[i]*(1 - c1/sigu2i[i]) - c1*a2i/2*(1 - 1.5*c1/sigu2i[i])))*t(zui[i,,drop=FALSE])%*%zui[i,,drop=FALSE]/sigu2i[i]
  
  if(!eff.time.invariant){
   if(model == "BC1992"){
    Hbeta = Hbeta + s*t(xi)%*%(Gi*tBC1/sigv2i)*sqrt(c1)*(a2i + d2) - s*t(xi)%*%(Gi/sigv2i)*c1*(-s*sum(c2*tBC1)*d3 + sqrt(c1)*a2i*(1 + d3)*sum(Gi^2*tBC1/sigv2i) + sqrt(c1)*d2*sum(Gi^2*tBC1/sigv2i))
    Hdeleta = Hdeleta + t((sum(Gi^2*tBC1/sigv2i)*c1^1.5*(a2i*(1+d3) + d2)-s*c1*sum(c2*tBC1)*d3)*zdeli[i,,drop=FALSE]/sigu2i[i])
    Hetagu = Hetagu + t((sum(Gi^2*tBC1/sigv2i)*c1^1.5/sigu2i[i]*(sqrt(c1) -a2i*mui[i]*(1+d3) + sqrt(c1)*a2i^2*0.5*(3+d3) - d2*(mui[i] - 1.5*a2i*sqrt(c1))) + s*c1/sigu2i[i]*sum(c2*tBC1)*(mui[i]*d3 - d2*0.5*sqrt(c1) - 0.5*sqrt(c1)*a2i*(1+d3)))*zui[i,,drop=FALSE])
    Hetagv = Hetagv - t(zvi)%*%(Gi^2*tBC1/sigv2i)*c1*(1 + a2i^2 + a2i*d2) + sum(Gi^2*tBC1/sigv2i)*t(zvi)%*%(Gi^2/sigv2i)*c1^2*(1 + 0.5*a2i^2*(3 + d3) + 1.5*a2i*d2) + t(zvi)%*%(c2*tBC1)*s*sqrt(c1)*(d2 + a2i) + sum(c2*tBC1)*t(zvi)%*%(c2)*c1*d3 - (a2i*(1+d3) + d2)*s*c1^1.5*(0.5*sum(c2*tBC1)*t(zvi)%*%(Gi^2/sigv2i) + sum(Gi^2*tBC1/sigv2i)*t(zvi)%*%(c2))
    Heta = Heta - 2*c1*sum(Gi^2*tBC1^2/sigv2i)*(1 + a2i^2 + d2*a2i) + c1^2*sum(Gi^2*tBC1/sigv2i)^2*(2 + a2i^2*(3 + d3) + 3*d2*a2i) + sum(c2*tBC1^2)*s*sqrt(c1)*(d2 + a2i) - 2*s*c1^1.5*sum(c2*tBC1)*sum(Gi^2*tBC1/sigv2i)*(d2 + a2i*(1 + d3)) + sum(c2*tBC1)^2*c1*d3
    
   } else if (model == "K1990modified"){
    Hbeta = Hbeta - s*t(xi)%*%sweep(tBC2, 1, sigv2i, "/")*sqrt(c1)*(a2i + d2) - s*(t(xi)%*%(Gi/sigv2i))%*%t(t(tBC2)%*%(c1*(s*(ei/sigv2i)*d3 - sqrt(c1)*a2i*(1 + d3)*(Gi/sigv2i) - sqrt(c1)*d2*(Gi/sigv2i))))
    Hdeleta = Hdeleta + t(zdeli[i,,drop=FALSE])%*%(t(t(tBC2)%*%(s*c1*(ei/sigv2i)*d3 - c1^1.5*(a2i*(1+d3) + d2)*(Gi/sigv2i)))/sigu2i[i])
    Hetagu = Hetagu + t(zui[i,,drop=FALSE])%*%(- t(Gi/sigv2i)%*%tBC2*c1^1.5/sigu2i[i]*(sqrt(c1) -a2i*mui[i]*(1+d3) + sqrt(c1)*a2i^2*0.5*(3+d3) - d2*(mui[i] - 1.5*a2i*sqrt(c1))) - t(ei/sigv2i)%*%tBC2*s*c1/sigu2i[i]*(mui[i]*d3 - d2*0.5*sqrt(c1) - 0.5*sqrt(c1)*a2i*(1+d3)))
    Hetagv = Hetagv + t(zvi)%*%sweep(tBC2, 1, Gi/sigv2i, "*")*c1*(1 + a2i^2 + a2i*d2) - (t(zvi)%*%(Gi^2/sigv2i))%*%(t(Gi/sigv2i)%*%tBC2)*c1^2*(1 + 0.5*a2i^2*(3 + d3) + 1.5*a2i*d2) - t(zvi)%*%sweep(tBC2, 1, ei/sigv2i, "*")*s*sqrt(c1)*(d2 + a2i) - (t(zvi)%*%(c2))%*%(t(ei/sigv2i)%*%tBC2)*c1*d3 + (a2i*(1+d3) + d2)*s*c1^1.5*(0.5*(t(zvi)%*%(Gi^2/sigv2i))%*%(t(ei/sigv2i)%*%tBC2) + (t(zvi)%*%(c2))%*%(t(Gi/sigv2i)%*%tBC2))
    Heta = Heta - c1*(1 + a2i^2 + a2i*d2)*t(tBC2)%*%sweep(tBC2, 1, sigv2i, "/") + c1^2*(2 + a2i^2*(3 + d3) + d2*a2i*3)*t(t(Gi/sigv2i)%*%tBC2)%*%(t(Gi/sigv2i)%*%tBC2) - c1^1.5*s*(d2+a2i*(1+d3))*t(t(ei/sigv2i)%*%tBC2)%*%(t(Gi/sigv2i)%*%tBC2) + c1*d3*t(t(ei/sigv2i)%*%tBC2)%*%(t(ei/sigv2i)%*%tBC2) - c1^1.5*s*(d2 + a2i*(1+d3))*t(t(Gi/sigv2i)%*%tBC2)%*%(t(ei/sigv2i)%*%tBC2)
    
   } else if (model == "K1990"){
    Hbeta = Hbeta +  s*t(xi)%*%sweep(tSK, 1, Gi^2*Gk/sigv2i, "*")*sqrt(c1)*(a2i + d2) - s*(t(xi)%*%(Gi/sigv2i))%*%t(t(tSK)%*%(c1*(-s*(ei*Gi^2*Gk/sigv2i)*d3 + sqrt(c1)*a2i*(1 + d3)*(Gi^3*Gk/sigv2i) + sqrt(c1)*d2*(Gi^3*Gk/sigv2i))))
    
    Hdeleta = Hdeleta + t(zdeli[i,,drop=FALSE])%*%(s*t(t(tSK)%*%(-s*c1*(c2*Gi*Gk)*d3 + c1^1.5*(a2i*(1+d3) + d2)*(Gi^3*Gk/sigv2i)))/sigu2i[i])
    
    Hetagu = Hetagu + t(zui[i,,drop=FALSE])%*%(t(Gi^3*Gk/sigv2i)%*%tSK*c1^1.5/sigu2i[i]*(sqrt(c1) -a2i*mui[i]*(1+d3) + sqrt(c1)*a2i^2*0.5*(3+d3) - d2*(mui[i] - 1.5*a2i*sqrt(c1))) + t(c2*Gi*Gk)%*%tSK*s*c1/sigu2i[i]*(mui[i]*d3 - d2*0.5*sqrt(c1) - 0.5*sqrt(c1)*a2i*(1+d3)))
    
    Hetagv = Hetagv - t(zvi)%*%sweep(tSK, 1, Gi^3*Gk/sigv2i, "*")*c1*(1 + a2i^2 + a2i*d2) + (t(zvi)%*%(Gi^2/sigv2i))%*%(t(Gi^3*Gk/sigv2i)%*%tSK)*c1^2*(1 + 0.5*a2i^2*(3 + d3) + 1.5*a2i*d2) + t(zvi)%*%sweep(tSK, 1, c2*Gi*Gk, "*")*s*sqrt(c1)*(d2 + a2i) + (t(zvi)%*%(c2))%*%(t(c2*Gi*Gk)%*%tSK)*c1*d3 - (a2i*(1+d3) + d2)*s*c1^1.5*(0.5*(t(zvi)%*%(Gi^2/sigv2i))%*%(t(c2*Gi*Gk)%*%tSK) + (t(zvi)%*%(c2))%*%(t(Gi^3*Gk/sigv2i)%*%tSK))
    
    Heta = Heta + c1*(1 + a2i^2 + a2i*d2)*t(tSK)%*%sweep(tSK, 1, Gi^3*Gk*(1- 3*Gi*Gk)/sigv2i, "*") + c1^2*(2 + a2i^2*(3 + d3) + d2*a2i*3)*t(t(Gi^3*Gk/sigv2i)%*%tSK)%*%(t(Gi^3*Gk/sigv2i)%*%tSK) - s*sqrt(c1)*(a2i + d2)*t(tSK)%*%sweep(tSK, 1, ei*Gi^2*Gk*(1 - 2*Gi*Gk)/sigv2i, "*") - c1^1.5*s*(d2 + a2i*(1 + d3))*(t(t(ei*Gi^2*Gk/sigv2i)%*%tSK)%*%(t(Gi^3*Gk/sigv2i)%*%tSK) + t(t(Gi^3*Gk/sigv2i)%*%tSK)%*%(t(ei*Gi^2*Gk/sigv2i)%*%tSK)) + c1*d3*t(t(ei*Gi^2*Gk/sigv2i)%*%tSK)%*%(t(ei*Gi^2*Gk/sigv2i)%*%tSK)
   } else {
    stop("Unknown model\n")
   }
  }
 }
 
 # if(eff.time.invariant){
 #  H <- cbind(rbind(Hb,	t(Hbgv),	t(Hbgu),		t(Hbdel)), 
 #             rbind(Hbgv, 	Hgv,		t(Hgvgu),		t(Hgvdel)),
 #             rbind(Hbgu, 	Hgvgu,		Hgu,			Hdelgu), 
 #             rbind(Hbdel,	Hgvdel,		t(Hdelgu),		Hdel))
 # } else {
 #  H <- cbind(rbind(Hb,	t(Hbgv),	t(Hbgu),		t(Hbdel),	t(Hbeta)), 
 #             rbind(Hbgv, 	Hgv,		t(Hgvgu),		t(Hgvdel),	t(Hetagv)),
 #             rbind(Hbgu, 	Hgvgu,		Hgu,			Hdelgu, 	t(Hetagu)), 
 #             rbind(Hbdel,	Hgvdel,		t(Hdelgu),		Hdel,		t(Hdeleta) ),
 #             rbind(Hbeta,	Hetagv,		Hetagu,			Hdeleta,	Heta ))
 # }

  H <- H[!coef.fixed,!coef.fixed]
 # print(H)
 
 # if(mean.u.0i.zero){
 #  grad <- grad[-(k+kv+ku+kdel)]
 #  H <- H[-(k+kv+ku+kdel),-(k+kv+ku+kdel)]
 # } 
 
 # print(length(grad))
 # print(dim(H))
 # print(length(theta))
 
 
 colnames(H) <- rownames(H) <- names(theta)
 return(H)
}

# .hess.panel <- function(theta, prod, my.n, k, kv, ku, kdel, yit, zvit, zui, xit, zdeli, model, timevar, maxT, ids, idvar, t0, eff.time.invariant = FALSE){
#  numDeriv::hessian(func = .ll.panel, x = theta, prod = prod, my.n = my.n, k = k, kv = kv, ku = ku, kdel = kdel, yit = yit, zvit = zvit, zui = zui, xit = xit, zdeli = zdeli, eff.time.invariant = eff.time.invariant, model = model, timevar = timevar, maxT = maxT, ids = ids, idvar = idvar, t0 = t0)
# }
# 
# .hess.panel <- function(theta, prod, my.n, k, kv, ku, kdel, yit, zvit, zui, xit, zdeli, model, timevar, maxT, ids, idvar, t0, eff.time.invariant = FALSE){
#  hessian1(func = .ll.panel, at = theta, prod = prod, my.n = my.n, k = k, kv = kv, ku = ku, kdel = kdel, yit = yit, zvit = zvit, zui = zui, xit = xit, zdeli = zdeli, eff.time.invariant = eff.time.invariant, model = model, timevar = timevar, maxT = maxT, ids = ids, idvar = idvar, t0 = t0)
# }

.is.negative.definite <- function (x, tol = 1e-16)
{
 # if (!is.square.matrix(x))
 # stop("argument x is not a square matrix")
 # if (!is.symmetric.matrix(x))
 # stop("argument x is not a symmetric matrix")
 # if (!is.numeric(x))
 # stop("argument x is not a numeric matrix")
 eigenvalues <- eigen(x, only.values = TRUE)$values
 # cat("Eigenvalues\n")
 # print(eigenvalues)
 n <- nrow(x)
 for (i in 1:n) {
  if (abs(eigenvalues[i]) < tol) {
   eigenvalues[i] <- 0
  }
 }
 if (any(eigenvalues >= 0)) {
  return(FALSE)
 }
 return(TRUE)
}

.make.neg.def <- function(hess){
 # K4 <- ncol(hess)
 # h0_ <- matrix(0, K4, K4)
 eigen1 <- eigen( hess )
 # eigen.val <- abs(eigen1$values)
 #   for( i in seq_len( K4 ) ){
 #     h0_ <- h0_ - abs(eigen1$values[i]) * eigen1$vectors[,i] %*% t(eigen1$vectors[,i])
 #     # h0_ <- h0_ - eigen.val[i] * eigen1$vectors[,i] %*% t(eigen1$vectors[,i])
 #   }
 return(-crossprod(t(eigen1$vectors)*abs(eigen1$values), t(eigen1$vectors)))
}

.make.invertible <- function(hess){
 K <- ncol(hess)
 ok <- FALSE
 adj <- sqrt(.Machine$double.eps)
 # tr0 <- sum( 1/eigen(H,only.values = T)$values )
 i <- 0
 while(!ok & i < 1000){
  i <- 1 + i
  hess <- hess + adj*diag(rep(1, K))
  mytry <- tryCatch( solve(-hess), error = function(e) e )
  ok <- !inherits( mytry, "error")
  # print(mytry)
  # ok <- all(diag(H) != 0)
  adj <- adj * 2
  # print(adj)
 }
 # print(c(i, adj))
 return(hess)
}

.make.zero.one <- function(qq){
 for (i in 1:length(qq)) {
  qqq <- qq[i]
  if(abs(qqq) > 1){
   while(abs(qqq) > 1){
    qqq <- qqq / 10
   }
   qq[i] <- qqq
  }
 }
 return(qq)
}

.mlmaximize.panel <- function(theta0, ll, gr = NULL, hess = NULL, gr.hess = NULL, alternate = NULL, BHHH = F, level = 0.99, step.back = .Machine$double.eps^.5, reltol =  sqrt(.Machine$double.eps) , lmtol =  sqrt(.Machine$double.eps) , steptol =  .Machine$double.eps, digits = 4, when.backedup = sqrt(.Machine$double.eps), max.backedup = 17, print.level = 6, only.maximize = FALSE, maxit = 150, n = 100, ...){
 
 theta00 <- theta0
 
 k4 <- length(theta0)
 
 # if(print.level >= 6){
 #  cat("\n=================")
 #  cat(" Initial values:\n\n", sep = "")
 #  print(theta0)
 # }
 # 
 # if(print.level >= 2){
 #  cat("\n=================")
 #  cat(" Maximization:\n\n", sep = "")
 # }
 
 # step.back = 2^-217
 
 ll0 <- ll(theta0, ...)
 ltol <- reltol * (abs(ll0) + reltol)
 # print(ltol)
 typf <- ll0
 theta1 <- theta0
 
 iter <- iter.total <- backedup <- backedups <- wasconcave <- wasconcaves <- 0
 
 if( is.na(ll0) | ll0 == -Inf ){
  if(print.level >= 2){
   cat("Could not compute ll at starting values: trying something else\n")
  }
  iter1 <- backedups
  repeat{
   iter1 <- iter1 + 1
   theta0 <- theta0 * runif(length(theta0), 0.98, 1.02) # not sure what to do here
   ll0 <- ll( theta0, ... )
   if(!is.na(ll0)) break
   if(iter1 == 55){
    stop("it's not gonna happen...")
   }
  }
  # backedups <- iter1
 }
 
 delta1 <- gHg <- s1 <- 1
 h1 <- tryCatch( 2, error = function(e) e )
 cant.invert.hess <- FALSE
 
 if(print.level >= 2){
  cat(paste("Iteration ",formatC(iter, width = 3)," (at starting values):       log likelihood = ",format(ll0, digits = 13),"\n\n", sep = ""), sep = "")
 }
 
 repeat{
  iter.total <- iter.total + 1
  # cat("backedup = ",backedup," backedups = ",backedups,"\n", sep = "")
  if(print.level >= 6){
   print(theta0)
  }
  
  # cumulate how many times did it backed-up in a row
  if(s1 < when.backedup){
   backedup <- backedup + 1
  } else {
   backedup <- 0
  }
  # print(s1)
  # cumulate how many times was concave
  if( inherits(h1, "error") ){
   wasconcave <- wasconcave + 1
  } else {
   wasconcave <- 0
  }
  
  # try different values if was concave more than @@@ times
  if(wasconcave == max.backedup){
   # start over
   if(print.level >= 2){
    cat("Not concave ",max.backedup," times in a row: trying something else (not concave ",wasconcaves+1," times in total)\n")
   }
   iter <- wasconcave <- backedup <- 0
   wasconcaves <- wasconcaves + 1
   theta0 <- theta0*runif(length(theta0), 0.9999, 1.0001) # not sure what to do here
   ll0 <- ll( theta0, ... )
   # if(print.theta) print(theta0)
   if(print.level >= 2){
    cat(paste("Iteration ",formatC(iter, width = 3),"  (at slightly perturbed starting values):       log likelihood = ",format(ll0, digits = 13),"\n\n", sep = ""), sep = "")
   }
  }
  
  # try different values if backed-up more than @@@ times
  if(backedup == max.backedup){
   # start over
   if(print.level >= 2){
    cat("Backed-up ",max.backedup," times in a row: trying something else (backup-up ",backedups+1," times in total)\n", sep = "")
   }
   iter <- backedup <- wasconcave <- 0
   backedups <- backedups + 1
   wasconcaves <- wasconcaves + 1
   theta0 <- theta0*runif(length(theta0), 0.9999, 1.0001) # not sure what to do here
   ll0 <- ll( theta0, ... )
   if(print.level >= 6){
    print(theta0)
   }
   if(print.level >= 2){
    cat(paste("Iteration ",formatC(iter, width = 3)," (at slightly perturbed starting values):       log likelihood = ",format(ll0, digits = 13),"\n\n", sep = ""), sep = "")
   }
  }
  
  # see if it calculated ll
  if( is.na(ll0) | ll0 == -Inf | ll0 == Inf | ll0 == 0 ){
   if(print.level >= 2){
    cat("Could not compute ll: trying something else\n")
   }
   iter1 <- backedups
   repeat{
    iter1 <- iter1 + 1
    # theta0 <- c( cons0, beta0, mu = 0, eta = 0, lnsv2 = -1*iter1/2, lnsu2 = -1*iter1/2)
    theta0 <- theta00*runif(length(theta0), 0.9999, 1.0001) # not sure what to do here
    ll0 <- ll( theta0, ... )
    if(!is.na(ll0) & ll0 != 0) break
    if(iter1 == 15){
     stop("it's not gonna happen... could not compute at initial and find feasible values")
    }
   }
   iter <- backedup <- 0
   backedups <- iter1
   if(print.level >= 2){
    cat(paste("Iteration ",formatC(iter, width = 3)," (at slightly perturbed starting values):       log likelihood = ",format(ll0, digits = 13),"\n\n", sep = ""), sep = "")
   }
  }
  if(backedups == 5){
   stop("it's not gonna happen... backuped 5 times")
  }
  if(wasconcaves == 5){
   stop("it's not gonna happen... not concave 5 times")
  }
  iter <- iter + 1
  delta3 <- 1
  # step 2: direction vector
  
  # previous theta
  
  if(iter.total > 1) theta1 <- theta0 - s1 * d0
  
  # BHHH (faster, but different):
  # The Hessian is approximated as the negative
  # of the sum of the outer products of the gradients
  # of individual observations,
  # -t(gradient) %*% gradient = - crossprod( gradient )

  # gradient and hessian together
  g1_h0 <- gr.hess(theta0, ...)
  g1 <- g1_h0$grad
  h0 <- g1_h0$hessian
  
  
  # gradient and hessian separately
  # g1 <- gr(theta0, ...)
  # # print(g1)
  # if(!is.null(alternate)) BHHH <- floor(iter/alternate) %% 2 == 1
  # h0 <- hess(theta0,  ...)
  
  # print(h0)
  
  
  # check if negative definite

  
  # eigen1 <- eigen( h0 )
  # # eigen.tol <- k4 * max(abs(eigen1$values)) * .Machine$double.eps # this is for positive definiteness
  # eigen.val <- ifelse(eigen1$values < .Machine$double.eps^.1, 0, eigen1$values)
  # eigen.val <- eigen1$values
  # # hess.pos.def <- sum(eigen1$values > eigen.tol) == k4
  # hess.neg.def <- !any(eigen.val >= 0)
  # 1. replace negative with small ones
  # eigen.val <- ifelse(eigen1$values < 0, .0001, eigen.val)
  # 2. replace negative with absolut values
  # eigen.val <- abs(eigen1$values)

  
  hess.neg.def <- .is.negative.definite(h0)
  h00 <- h0
  if(!hess.neg.def) h0 <- .make.neg.def(h0)
  # print(hess.neg.def)
  # make it negative definite if it is not already
  # if(!hess.neg.def){
  #  # cat(" !hess.neg.def; k4 =",k4,"length(eigen1$vectors[,i]) = ",length(eigen1$vectors[,1]),"\n")
  #  h0_ <- matrix(0, k4, k4)
  #  # eigen1 <- eigen( h0 )
  #  for( i in seq_len( k4 ) ){
  #   # h0_ <- h0_ - abs(eigen1$values[i]) * eigen1$vectors[,i] %*% t(eigen1$vectors[,i])
  #   h0_ <- h0_ - eigen.val[i] * crossprod(t(eigen1$vectors[,i]))
  #  }
  #  h0.old <- h0
  #  h0 <- h0_
  # }
  # print( is.negative.definite(h0) )
  # remember hessian and negative of its inverse from previous iter that could have been inverted
  # if( !cant.invert.hess ){
  #  h0_previous <- h0
  #  h1_previous <- h1
  # }
  # easier to invert positive definite matrix
  # h1 <- tryCatch( qr.solve(-h0, tol = 1e-16), error = function(e) e )
  h1 <- tryCatch( solve(-h0), error = function(e) e )
  if(inherits(h1, "error")) h0 <- .make.invertible(h0)
  h1 <- solve(-h0)
  # print(diag(h1))
  # check if it can be inverted
  cant.invert.hess <- FALSE
  cant.invert.hess <- inherits(h1, "error")
  # print(cant.invert.hess)
  if( cant.invert.hess ){
   # # print(h1)
   # if(print.level >= 2){
   #  cat(paste("cannot invert Hessian, using eigenvalues\n", sep = ""), sep = "")
   # }
   # # this was just to get the uninvertable hessian
   # # return(list(hess = h0, grad = g1))
   # # stop("this")
   # # @14@ this
   # eig1 <- eigen( -h0_previous )
   # d0 <- rep(0, length(theta0))
   # # eig2 <- ifelse(eig1$values < eps1, 1, eig1$values)
   # for (i in 1:length(eig1$values)){
   #  d0 <- d0 + (g1%*%eig1$vectors[,i])%*%eig1$vectors[,i] / eig1$values[i]
   # }
   # # @14@ could be done easier
   # # d0 <- qr.solve(-h0, g1, tol = 1e-134)
   # # print(g1)
   # # print(dim(g1))
   # # print(d0)
   # # print(dim(d0))
   # gHg <- sum( as.vector(g1) * d0)
   # # in the part of the ortogonal subspace where the eigenvalues
   # # are negative or small positive numbers, use steepest ascent
   # # in other subspace use NR step
   # # d0 <- ifelse(eigen(-h0, only.values = TRUE)$values < reltol, g1, d0)
   # gg <- sqrt( crossprod(g1) )
   # gHg <- gg
   # # d0 <- g1
   # # d0
  } else {
   # print(h1)
   # cat("dim H_inv",dim(h1)," length g1 ", length(g1), "\n")
   d0 <- as.vector( h1 %*% g1 )
   # d0 <- g1 # steepest
   gg <- sqrt( crossprod(g1) )
   # h1.old <- solve(-h0.old)
   gHg <- as.vector( t(g1) %*% h1 %*% g1 )
  }
  # gg_scaled <- gg * max( crossprod(theta0), crossprod(theta1) ) / max( abs(ll0), abs(typf))
  # theta_rel <- max( abs(theta0 - theta1) / apply( cbind( abs(theta0),abs(theta1) ), 1, max) )
  theta_rel <- max( abs(theta0 - theta1) / (abs(theta1)+1) )
  
  
  # begin stopping criteria calculated using new values of g1 and h1
  # if(s1 > when.backedup*10^-100 & delta1 != 17.17){ # if(s1 > when.backedup*10^-100 & !cant.invert.hess){
  if(s1 > when.backedup*1e-10 & delta1 != 17.17){ # if(s1 > when.backedup*10^-100 & !cant.invert.hess){
   if(abs(gHg) < lmtol & iter.total > 1){
    if(print.level >= 2){
     cat("\nConvergence given abs(gHg) = ",abs(gHg)," < lmtol\n", sep = "")
    }
    break
   }
   if(theta_rel < steptol*1e-100 & iter.total > 2){
    # print(theta_rel)
    if(print.level >= 2){
     cat("\nConvergence given relative change in theta = ",theta_rel," < steptol\n", sep = "")
    }
    break
   }
  }
  # end stopping criteria
  
  
  # use steepest ascent when backed-up
  if(s1 < when.backedup | delta1 == 17.17){ # was 1e-3
   # print("try this\n\n")
   ll0 <- ll( theta1, ... )
   eig1 <- eigen( -h0 )
   # cat("Gradient\n")
   # print(g1)
   # cat("Eigenvalues of the Hessian\n")
   # print(eig1$values)
   # cat("step before transformation 0\n")
   # print(d0)
   d0 <- as.vector( solve(-.make.neg.def(h0)) %*% g1 )
   # cat("step before transformation 1\n")
   # print(d0)
   # d0 <- ifelse(eig1$values < reltol, 0, d0)
   d0 <- ifelse(eig1$values < -1e-2, 0, d0)
   # d0 <- ifelse(abs(d0) > 1, .make.zero.one(d0), d0)
   # cat("step after transformation\n")
   # print(d0)
   # cat("\n\n")
   # theta0 <- theta0 - 1 * d0
  }
  
  
  
  # print(d0)
  # step 3: new guess
  # a: s = 1
  # b: funct(theta0 + d0) > funct(theta0)
  s1 <- 1
  # s1_ <- Inf
  # try.back.up.n <- 30
  # s1s <- 2^(3:-11)#seq(-5-1e-4, 5-1e-4, length.out = try.back.up.n)
  # my.loc.ll.max <- -1e100
  # for(qq in 1:15){
  #  theta1_ <- theta0 + s1s[qq] * d0
  #  ll1_ <- ll( theta1_, ... )
  #  if(is.finite(ll1_) & ll1_> my.loc.ll.max & ll1_ - my.loc.ll.max < 1e6){
  #   my.loc.ll.max <- ll1_
  #   s1_ <- s1s[qq]
  #  }
  # }
  # if(is.finite(s1_)){
  #  s1 <- s1_
  #  cat("opt s1 =", s1_,"opt ll ==",my.loc.ll.max, "ll before =" ,ll0,"\n")
  # } 
  # print(s1)
  # print(d0)
  # print(dim(d0))
  # print(theta0)
  # print(dim(theta0))
  theta1 <- theta0 + s1 * d0
  # print(12)
  # print(theta1)
  ll1 <- ll( theta1, ... )
  # print(13)
  delta2 <- ll1 - ll0
  flag <- (!is.na(delta2) & is.finite(delta2) & delta2 > 0 & delta2 < 1e+7)
  # begin Cases
  if( flag ){
   # begin Case 1: f(theta1) > f(theta0)
   ll.temp <- ll0
   # check if s1 = 2, 3, ... increases f even more
   while( flag ){
    if(print.level >= 6){
     cat(paste("\t\tCase 1: s = ",s1,", delta = ",delta2,"\n", sep = ""), sep = "")
    }
    # cat("ll1 from before = ",ll1,"\n")
    ll0 <- ll1
    s1 <- s1 + 1
    theta1 <- theta0 + s1 * d0
    ll1 <- ll( theta1, ... )
    delta2 <- ll1 - ll0
    flag <- (!is.na(delta2) & is.finite(delta2) & delta2 > 0)
    if(is.infinite(ll1)){
     cat("s1 = ",s1,"\n")
     cat("ll1 = ",ll1,"\n")
     cat("ll0 = ",ll0,"\n")
     cat("ll.temp = ",ll.temp,"\n")
     cat("ll( theta0, ... ) = ",ll( theta0, ... ),"\n")
     cat("ll( theta1, ... ) = ",ll( theta1, ... ),"\n")
     # s1 <- s1 - 1
    }
   }
   # overall delta
   delta1 <- ll0 - ll.temp
   delta_rel <- abs(delta1 / ll.temp)
   # print(delta_rel)
   s1 <- s1 - 1
   # overwrite the values
   theta0 <- theta0 + s1 * d0
   # end Case 1: f(theta1) > f(theta0)
  } else {
   # begin Case 2: f(theta1) < f(theta0)
   # check only if s1=1/2 increases f
   s1 <- 0.5
   # s1_ <- Inf
   # try.back.up.n <- 30
   # s1s <- s1s <- 2^(3:-26)#seq(1e-10, 1-1e-4, length.out = try.back.up.n)
   # my.loc.ll.max <- -1e100
   # for(qq in 1:try.back.up.n){
   #  theta1_ <- theta0 + s1s[qq] * d0
   #  ll1_ <- ll( theta1_, ... )
   #  if(is.finite(ll1_) & ll1_> my.loc.ll.max & ll1_ - my.loc.ll.max  < 1e7){
   #   my.loc.ll.max <- ll1_
   #   s1_ <- s1s[qq]
   #  }
   # }
   # if(is.finite(s1_)){
   #  s1 <- s1_
   #  cat("opt s1 =", s1_,"opt ll ==",my.loc.ll.max, "ll before =" ,ll0,"\n")
   # }
   theta1 <- theta0 + s1 * d0
   ll1 <- ll( theta1, ... )
   # cat(" ll1 = ",ll1,", ll0 = ",ll0,", s = ",s1,"\n", sep = "")
   delta2 <- ll1 - ll0
   if(print.level >= 6){
    cat(paste("\t\tCase 2: s = ",s1,", delta = ",delta2,"\n", sep = ""), sep = "")
   }
   flag2 <- (!is.na(delta2) & is.finite(delta2) & delta2 > 0)
   # end Case 2: f(theta1) < f(theta0)
   if( flag2 ){
    # begin Case 2a: f(theta1) > f(theta0)
    ll.temp <- ll0
    # check if s1=1/2^2,1/2^3,... increases f even more
    while( flag2 ){
     if(print.level >= 6){
      cat(paste("\t\t\tCase 2a: s = ",s1,", delta = ",delta2,"\n", sep = ""), sep = "")
     }
     ll0 <- ll1
     s1 <- 0.5 * s1
     theta1 <- theta0 + s1 * d0
     ll1 <- ll( theta1, ... )
     # cat(" ll1 = ",ll1,", ll0 = ",ll0,", s = ",s1,"\n", sep = "")
     delta2 <- ll1 - ll0
     flag2 <- (!is.na(delta2) & is.finite(delta2) & delta2 > ltol)
    }
    # overall delta
    delta1 <- ll0 - ll.temp
    delta_rel <- abs(delta1 / ll.temp)
    # print(delta_rel)
    s1 <- 2 * s1
    # overwrite the values
    theta0 <- theta0 + s1 * d0
    # end Case 2a: f(theta1) > f(theta0)
   } else {
    # begin Case 2b: f(theta1) < f(theta0)
    ll.temp <- ll0
    # try s1=1/2^2,1/2^3,... so that f(theta1) > f(theta0)
    while ( !flag2 & s1 > step.back^2 ){
     s1 <- 0.5 * s1
     theta1 <- theta0 + s1 * d0
     ll1 <- ll( theta1, ... )
     delta2 <- ll1 - ll0
     if(print.level >= 6){
      cat(paste("\t\t\tCase 2b: s = ",s1,", delta = ",delta2,"\n", sep = ""), sep = "")
     }
     flag2 <- (!is.na(delta2) & is.finite(delta2) & delta2 > 0)
    }
    if( !flag2 | s1 < step.back^2 ){
     # stop("provide different starting values")
     delta1 <- 17.17
    } else {
     # overwrite the values
     delta1 <- delta2
     delta_rel <- abs(delta1 / ll.temp)
     ll0 <- ll1
     theta0 <- theta0 + s1 * d0
    }
    # end Case 2b: f(theta1) < f(theta0)
   }
  }
  
  
  if(print.level >= 2){
   if( cant.invert.hess ){
    cat(paste("Iteration ",formatC(iter, width = 3)," (hessian is ",ifelse(BHHH, "BHHH", "analytical"),", ",formatC(iter.total, width = 3)," in total):   log likelihood = ",format(ll0, digits = 13)," (not concave)\n", sep = ""), sep = "")
   } else if (s1 <= when.backedup) {
    cat(paste("Iteration ",formatC(iter, width = 3)," (hessian is ",ifelse(BHHH, "BHHH", "analytical"),", ",formatC(iter.total, width = 3)," in total):   log likelihood = ",format(ll0, digits = 13)," (backed up)\n", sep = ""), sep = "")
   } else {
    cat(paste("Iteration ",formatC(iter, width = 3)," (hessian is ",ifelse(BHHH, "BHHH", "analytical"),", ",formatC(iter.total, width = 3)," in total):   log likelihood = ",format(ll0, digits = 13),"\n", sep = ""), sep = "")
   }
  }
  
  # printing criteria
  if(print.level >= 5){
   if( cant.invert.hess ){
    cat(paste(" (in iter ",formatC(iter, width = 3),": delta = ",format(delta1, digits = 6),"; s = ",format(s1, digits = 6),"; quasi-gHg = ",format(gHg, digits = 6),"; theta_rel_change = ",format(theta_rel, digits = 6),")\n\n", sep = ""), sep = "")
    if(print.level >= 6.5){
     print(theta0)
    }
   } else {
    cat(paste(" (in iter ",formatC(iter, width = 3),": delta = ",format(delta1, digits = 6),"; s = ",format(s1, digits = 6),"; gHg = ",format(abs(gHg), digits = 6),"; theta_rel_change = ",format(theta_rel, digits = 6),")\n\n", sep = ""), sep = "")
    if(print.level >= 5.5){
     print(theta0)
    }
   }
   cat("Par\n")
   print(theta0)
   cat("Gradient\n")
   print(g1)
   cat("\n")
  }
  # print(s1)
  if(s1 > when.backedup & !cant.invert.hess){ # if(s1 > when.backedup^2 & !cant.invert.hess){
   # ltol <- reltol * (abs(ll0) + reltol)
   # print(cant.invert.hess)
   if(delta1 > 0 & !is.na(delta_rel) & delta_rel < ltol*1e-16 & iter.total > 1){
    if(print.level >= 2){
     cat("\nConvergence given delta = ",delta_rel," < dtol\n", sep = "")
    }
    g1_h0 <- gr.hess(theta0, ...)
    g1 <- g1_h0$grad
    h0 <- g1_h0$hessian
    
    # g1 <- gr(theta0, ...)
    # h0 <- hess(theta0,  ...)
    # h1 <- tryCatch( qr.solve(-h0), error = function(e) e )
    break
   }
  }
  if(iter.total > maxit){
   cat("\n Maximum number of iterations (",iter.total,") reached without convergence\n", sep = "")
   cat(" Only 'theta' from iteration ",iter," will be returned\n\n", sep = "")
   print(theta0)
   return(list(par = theta0, converged = 0))
   break
   cat(" Only 'theta' from iteration ",iter," will be returned\n", sep = "")
  }
 } # end repeat
 
 if( !only.maximize & !cant.invert.hess){
  names(ll0) <- NULL
  colnames(h1) <- rownames(h1) <- names(g1) <- names(theta0)
  
  # sqrt(crossprod(g1))
  
  b0 <- theta0
  suppressWarnings( sd0 <- sqrt( diag( h1 ) ) )
  sd.nas <- is.na(sd0)
  # cat.print(sum(sd.nas))
  if(sum(sd.nas) > 0){
   sd1 <- sqrt(diag( solve(-.make.neg.def(h00), tol = 1e-268 ) ))
   sd0[sd.nas] <- sd1[sd.nas]
  }
  # cat.print(diag( solve(-.make.neg.def(h00) ) ))
  # cat("1\n")
  # cat.print(diag( solve(-.make.invertible(h00) ) ))
  # cat.print(diag( solve(-.make.neg.def(h00), tol =1e-68 ) ))
  t0 <- b0 / sd0
  p0 <- pt(abs(t0), n-length(b0), lower.tail = FALSE) * 2
  t10 <- qt((1-0.01*level)/2, n-length(b0), lower.tail = FALSE)
  t17 <- cbind( b0, sd0, t0, p0, b0 - t10*sd0, b0 + t10*sd0)
  # t17 <- cbind( b0, sd0, t0, p0)
  colnames(t17) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)", paste("",level,"_CI_LB", sep = ""), paste("",level,"_CI_UB", sep = ""))
  # colnames(t17) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  # t17
  
  if(print.level >= 2){
   cat(paste("\nFinal log likelihood = ",format(ll0, digits = 13),"\n\n", sep = ""), sep = "")
  }
  # cat(paste("Stoc. frontier normal/",distribution,"\n", sep = ""), sep = "")
  if(print.level >= 7){
   cat("\nCoefficients:\n\n", sep = "")
   printCoefmat(t17[,1:4], digits = digits)
  }
  
  return(list(par = theta0, table = t17, gradient = g1, vcov = h1, ll = ll0, gg = gg, gHg = gHg, theta_rel_ch = theta_rel, converged = 1))
 } else {
  return(list(par = theta0, gradient = g1, vcov = h1, ll = ll0, gg = gg, gHg = gHg, theta_rel_ch = theta_rel, converged = 1))
 }
 
}

# Print the estimation results
.printpanel2nd = function(x, digits, kb, kvi, ku0, kdeli, model, eff.time.invariant, mean.u.0i.zero, na.print = "NA", max.name.length, mycutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), mysymbols = c("***", "**", "*", ".", " ")){
 
 Cf = cbind(ifelse(x[,1, drop = F]> 999, formatC(x[,1, drop = F], digits = 1, format = "e",width = 10), formatC(x[,1, drop = F], digits = digits, format = "f", width = 10)),
            ifelse(x[,2, drop = F]>999, formatC(x[,2, drop = F], digits = 1, format = "e", width = 10), formatC(x[,2, drop = F], digits = digits, format = "f", width = 10)),
            ifelse(x[,3, drop = F]>999, formatC(x[,3, drop = F], digits = 1, format = "e", width = 7), formatC(x[,3, drop = F], digits = 2, format = "f", width = 7)),
            ifelse(x[,4, drop = F]>999, formatC(x[,4, drop = F], digits = 1, format = "e", width = 10), formatC(x[,4, drop = FALSE], digits = digits, format = "f", width = 10)),
            formatC(mysymbols[findInterval(x = x[,4], vec = mycutpoints)], flag = "-"))
 
 row.names(Cf) <- formatC(row.names(Cf), width = max(nchar(row.names(Cf))), flag = "-")
 cat("",rep(" ", max.name.length+6),"Coef.        SE       z       P>|z|\n", sep = "")
 dimnames(Cf)[[2]] <- rep("", dim(Cf)[[2]])
 cat("",rep("_", max.name.length+42-1),"", "\n", "Frontier", "\n", sep = "")
 print.default(Cf[1:kb,,drop=FALSE], quote = FALSE, right = TRUE, na.print = na.print)
 
 
 cat("",rep("-", max.name.length+42-1),"", "\n", "Random noise component: log(sigma_vit^2)", "\n", sep = "")
 print.default(Cf[(kb+1):(kb+kvi),,drop=FALSE], quote = FALSE, right = TRUE, na.print = na.print)
 
 cat("",rep("-", max.name.length+42-1),"", "\n", "Inefficiency component: log(sigma_u0i^2)", "\n", sep = "")
 
 cat("(id-specific time-invariant)", "\n", sep = "")
 print.default(Cf[(kb+kvi+1):(kb+kvi+ku0),,drop=FALSE], quote = FALSE, right = TRUE, na.print = na.print)
 
 if(!mean.u.0i.zero){
  cat("",rep("-", max.name.length+42-1),"", "\n", "Inefficiency component: conditional mean of the ", "\n", sep = "")
  cat("truncated-normal distribution", "\n", sep = "")
  cat("(id-specific time-invariant)", "\n", sep = "")
  print.default(Cf[(kb+kvi+ku0+1):(kb+kvi+ku0+kdeli),,drop=FALSE], quote = FALSE, right = TRUE, na.print = na.print)
 } else {
  kdeli <- 0
 }
 

 if(!eff.time.invariant){
  if(model == "BC1992"){
   Keta <- 1
   components <- "component"
  } else {
   Keta <- 2
   components <- "components"
  }
  cat(rep("-", max.name.length+42-1),"\n Decay ",components,"\n", sep = "")
  cat("(function of inefficiency time variation)", "\n", sep = "")
  print.default(Cf[(kb+kvi+ku0+kdeli+1):(kb+kvi+ku0+kdeli+Keta),,drop=FALSE], quote = FALSE, right = TRUE, na.print = na.print)
 } else {
  Keta <- 0
 }
 
 if(nrow(x) > kb+kvi+ku0+kdeli+Keta){
  cat("",rep("-", max.name.length+42-1),"", "\n", "Auxiliary parameters", "\n", sep = "")
  print.default(Cf[(kb+kvi+ku0+kdeli+Keta+1):nrow(x),,drop=FALSE], quote = FALSE, right = TRUE, na.print = na.print)
  
 }
 
 cat("",rep("_", max.name.length+42-1),"", "\n", sep = "")
 invisible(x)
}

# Technical efficiencies and prediction intervals
.u2efftnm.panel <- function( eit, sigv2it, sigu2i, mui, Git, alpha = alpha, prod = prod, cost.eff.less.one, ids, idvar, timevar, t0, it.names) {
 # if(prod){sn = -1} else {sn = 1}
 sn <- sn1 <- ifelse(prod, -1, 1)
 if(!prod & cost.eff.less.one) sn <- -1
 # sn <- -1
 
 te_jlms_mean <- te_jlms_mode <- te_bc <- te_l <- te_u <- c()
 
 # define indices
 m.end <- cumsum(t0)
 m.begin <- c(1,(m.end+1)[-length(t0)])
 
 for(i in 1:length(ids)){
  
  
 # for(i in ids){
  # sample.i <- idvar == i
  sample.i <- (m.begin[i]):(m.end[i])
  ei = eit[sample.i]; 
  sigv2i = sigv2it[sample.i]
  Gi = Git[sample.i];
  
  s2i = (1/sigu2i[i] + sum(Gi^2/sigv2i))^(-1)
  mi = (mui[i]/sigu2i[i] + sn1*sum(ei*Gi/sigv2i))*s2i
  zi <- mi / sqrt(s2i)
  
  point.est.mean <- mi + sqrt(s2i) * dnorm(zi) / pnorm(zi)
  point.est.mode <- ifelse( mi >= 0, mi, 0 )
  
  # if(i < 2){
  #  cat("\n current i = ",i ,"\n")
  #  cat.print(sample.i)
  #  cat.print(ei)
  #  cat.print(mui[i])
  #  cat.print(sigv2i)
  #  cat.print(Gi)
  #  cat.print(s2i)
  #  cat.print(mi)
  #  cat.print(zi)
  #  cat("\n")
  # }
  
  te_jlms_mean <- append(te_jlms_mean, exp(sn*Gi*point.est.mean))
  te_jlms_mode <- append(te_jlms_mode, exp(sn*Gi*point.est.mode))
  te_bc <- append(te_bc, exp(sn*mi*Gi + 0.5*s2i*Gi^2) * pnorm(zi + sn * sqrt(s2i)*Gi)/pnorm(zi))
  zl    <- qnorm( 1 - alpha / 2 * pnorm(zi) )
  zu    <- qnorm( 1 - ( 1 - alpha/2 ) * pnorm(zi) )
  te_l  <- append(te_l, exp(Gi*(sn*mi - zl*sqrt(s2i))))
  te_u  <- append(te_u, exp(Gi*(sn*mi - zu*sqrt(s2i))))
 }
 tymch <- data.frame(idvar, timevar, te_l, te_jlms_mean, te_jlms_mode, te_bc, te_u)
 colnames(tymch) <- c(it.names, "Lower bound","JLMS", "Mode", "BC","Upper bound" )
 row.names(tymch) <- NULL
 return(tymch)
}

# truncreg ---------------------------------------------------------------------

# a and b are lower and upper truncations
.h.trunc <- function(al, be, a1 = 1, b1 = 1, a = -Inf, b = Inf, print = FALSE)
 (ifelse( b == Inf, 0, b1*dnorm(be) ) - ifelse( a == -Inf, 0, a1*dnorm(al) )  ) /
 ( pnorm(be) - pnorm(al) )

# log likelihood
.ll.trunc <- function(y, x, beta, sigma, a, b) - 0.5*length(y)*log(2*pi) - length(y)*log(sigma) - 0.5*sigma^(-2)*sum((y - x %*% beta)^2) - sum( log( pnorm((b - x %*% beta)/sigma, log.p = FALSE) - pnorm((a - x %*% beta)/sigma, log.p = FALSE) ) )

# function with formula or matrices ---------------------------------------
.prepareYXllul <- function(formula, ll = -Inf, ul = Inf, data, subset, sysnframe = 1, ...) {
 # needed.frame <- sys.nframe() - 1
 needed.frame <- sys.nframe() - sysnframe
 # cat("sys.nframe() = ",sys.nframe(),"\n")
 # cat("needed.frame = ",needed.frame,"\n")
 mf0 <- match.call(expand.dots = FALSE, call = sys.call(sys.parent(n = needed.frame)))

 ll.clas <- class(ll)
 if (ll.clas == "formula"){
  if (length(all.vars(ll)) == 1){
   ll.form <- TRUE
  } else {
   stop("invalid formula for lower limit for left-truncation; 'll' must be '~ variableName'")
  }
 } else if (ll.clas == "numeric") {
  if (length(ll) == 1){
   ll.form <- FALSE
  } else {
   stop("invalid lower limit for left-truncation; 'll' must be a scalar")
  }
 } else {
  stop("invalid lower limit for left-truncation; 'll' must be a formula or a scalar")
 }

 ul.clas <- class(ul)
 if (ul.clas == "formula"){
  if (length(all.vars(ul)) == 1){
   ul.form <- TRUE
  } else {
   stop("invalid formula for upper limit for right-truncation; 'ul' must be '~ variableName'")
  }
 } else if (ul.clas == "numeric") {
  if (length(ul) == 1){
   ul.form <- FALSE
  } else {
   stop("invalid upper limit for right-truncation; 'ul' must be a scalar")
  }
 } else {
  stop("invalid upper limit for right-truncation; 'ul' must be a formula or a scalar")
 }

 form1 <- Formula::Formula(formula)

 # cat(" print(form1)\n", sep = "")
 # print(form1)

 if (ll.form) {
  form1 <- Formula(as.formula(paste("",deparse(form1, width.cutoff = 500L)," | ",ll[2]," ", sep = "")))
 }
 if (ul.form) {
  form1 <- Formula(as.formula(paste("",deparse(form1, width.cutoff = 500L)," | ",ul[2]," ", sep = "")))
 }

 # cat(" print(form1)\n", sep = "")
 # print(form1)

 # check if it is a matrix
 datasupplied <- !(match("data", names(mf0), 0) == 0)
 subssupplied <- !(match("subset", names(mf0), 0) == 0)

 if(subssupplied & !datasupplied){
  stop("Cannot specify 'subset' without specifying 'data'\n")
 }

 if(datasupplied){
  # begin get a logical vector equal TRUE if !missing
  # first using data and subset to get x without NA
  mf <- mf0
  mf$formula <- form1 #formula( form )
  m <- match(c("formula", "data", "subset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, sys.frame(sys.parent(n = needed.frame)))
  mt <- attr(mf, "terms")
  X <- as.matrix(model.matrix(mt, mf))
  # now get the names in the entire data
  # esample <- seq_len( nrow(data) ) %in% as.numeric(rownames(X))
  esample <- rownames(data) %in% rownames(X)
  if(sum(esample)==0){
   stop("no valid data points")
  }
  # print(table(esample))
  # end get a logical vector equal TRUE if !missing

  # get the data
  # print(form1)
  dataesample <- model.frame(form1, data = data[esample,])
  # cat(" dim(dataesample)\n", sep = "")
  # print(dim(dataesample))
  # print(head(dataesample))

  # this is my full LHS
  Y <- model.part(form1, data = dataesample, lhs = 1, drop = TRUE)
  # cat(" print(Y)\n", sep = "")
  # print(Y)
  n.full <- length(Y)
  # Y <- as.matrix( model.matrix(formula(form1, lhs = 1, rhs = 0), data = data[esample,]))
  # print(Y)
  X <- as.matrix( model.matrix(formula(form1, lhs = 0, rhs = 1), data = dataesample))
  # print(X)

  if (ll.form) {
   LL <- model.matrix(formula(form1, lhs = 0, rhs = 2), data = dataesample)[,2]
  } else {
   LL <- rep(ll, n.full)
  }

  if (ul.form) {
   if (ll.form) {
    UL <- model.matrix(formula(form1, lhs = 0, rhs = 3), data = dataesample, drop = TRUE)[,2]
   } else {
    UL <- model.matrix(formula(form1, lhs = 0, rhs = 2), data = dataesample, drop = TRUE)[,2]
   }
  } else {
   UL <- rep(ul, n.full)
  }

  # cat(" LL\n", sep = "")
  # print(LL)
  # cat(" UL\n", sep = "")
  # print(UL)


  # flag those that are required for regression
  flag <- (Y > LL & Y < UL)
  # cat(" dim(X)\n", sep = "")
  # print(dim(X))
  # cat(" flag\n", sep = "")
  # print(flag)


  # get subsets
  # cat(" Y\n", sep = "")
  Y  <- Y[flag]
  n <- length(Y)
  # cat(" X\n", sep = "")
  X  <- X[flag, , drop = FALSE]
  # print(dim(X))
  # cat(" LL\n", sep = "")
  LL <- LL[flag]
  # print(LL)
  # cat(" UL\n", sep = "")
  UL <- UL[flag]
  # print(UL)
 }
 # if data are not supplied
 else {
  # begin get a logical vector equal TRUE if !missing

  # first using data and subset to get XZ without NA
  mf <- mf0
  mf$formula <- formula( form1 )
  subsetsupplied <- !(match("subset", names(mf), 0) == 0)
  if(subsetsupplied) stop("Subset with matrices is not allowed")
  m <- match(c("formula"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- as.name("model.frame")
  # mf <- eval(mf, parent.frame())
  mf <- eval(mf, sys.frame(sys.parent(n = needed.frame)))
  mt <- attr(mf, "terms")
  Y <- as.matrix(model.response(mf))
  n <- nrow(Y)
  XZ <- as.matrix(model.matrix(mt, mf))
  # print(head(XZ))
  # get a logical vector equal TRUE if !missing
  with.na <- model.frame(mt, na.action = na.pass)
  esample <- rownames(with.na) %in% rownames(XZ)
  if(sum(esample)==0){
   stop("no valid data points")
  }
  # print(table(esample))
  # end get a logical vector equal TRUE if !missing

  # get the data

  # get number of vars in X
  mf <- mf0
  m <- match(c("formula"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, sys.frame(sys.parent(n = needed.frame)))
  mt <- attr(mf, "terms")
  X <- as.matrix(model.matrix(mt, mf))
  # print(head(X))
  k <- ncol(X)

  Y <- Y[esample,]
  n.full <- length(Y)
  X <- XZ[esample,]

  if (ll.form) {
   LL <- XZ[esample,(k+1)]
  } else {
   LL <- rep(ll, n.full)
  }

  if (ul.form) {
   if (ll.form) {
    UL <- XZ[esample,(k+2)]
   } else {
    UL <- XZ[esample,(k+1)]
   }
  } else {
   UL <- rep(ul, n.full)
  }

  # cat(" LL\n", sep = "")
  # print(LL)
  # cat(" UL\n", sep = "")
  # print(UL)


  # flag those that are required for regression
  flag <- (Y > LL & Y < UL)
  # cat(" dim(X)\n", sep = "")
  # print(dim(X))
  # cat(" flag\n", sep = "")
  # print(flag)

  # get subsets
  # cat(" Y\n", sep = "")
  Y  <- Y[flag]
  n <- length(Y)
  # cat(" X\n", sep = "")
  X  <- X[flag, , drop = FALSE]
  # print(dim(X))
  # cat(" LL\n", sep = "")
  LL <- LL[flag]
  # print(LL)
  # cat(" UL\n", sep = "")
  UL <- UL[flag]
  # print(UL)

 }

 tymch <- list(Y = Y, X = X, LL = LL, UL = UL, n = n, n.full = n.full, ll.form = ll.form, ul.form = ul.form, esample = esample, nontruncsample = flag)
 class(tymch) <- "npsf"
 return(tymch)
}

