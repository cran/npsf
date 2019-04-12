nobs.npsf <- function( object, ... ) {
 if(is.null(object)){
  stop( gettextf("There are no coefficients in '%s'.", deparse(substitute(object))) )
 }
 mymodel <- object$model
 valid.model <- mymodel %in% c("teradial", "tenonradial", "teradialbc", "nptestrts", "nptestind", "sf_sc", "truncreg.cs", "sf_p_1", "sf_p_2_K1990", "sf_p_2_K1990modified", "sf_p_2_BC1992","sf_p_4comp_homosk") # possible.models
 if(!valid.model){
  stop( gettextf("Not valid model in '%s'.", deparse(substitute(object))) )
 }
 # select the model
 if(mymodel == "tenonradial" | mymodel == "teradial" | mymodel == "teradialbc" | mymodel == "nptestrts" | mymodel == "nptestind"){
  n <- object$K
 } else {
  n <- object$n
 }
 return( n )
}