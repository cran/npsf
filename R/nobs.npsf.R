nobs.npsf <- function( object, ... ) {
 if(is.null(object)){
  stop( gettextf("There are no coefficients in '%s'.", deparse(substitute(object))) )
 }
 mymodel <- object$model
 valid.model <- mymodel %in% c("teradial", "tenonradial", "teradialbc", "nptestrts", "nptestind", "sfsc") # possible.models
 if(!valid.model){
  stop( gettextf("Not valid model in '%s'.", deparse(substitute(object))) )
 }
 # select the model
 if(mymodel == "tenonradial" | mymodel == "teradial" | mymodel == "teradialbc" | mymodel == "nptestrts" | mymodel == "nptestind"){
  n <- object$K
 }
 if(mymodel == "sfsc"){
  n <- nrow(object$eff)
 }
 return( n )
}