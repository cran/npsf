coef.npsf <- function( object, ... ) {
  if(is.null(object$coef) | is.null(object)){
    stop( gettextf("There are no coefficients in '%s'.", deparse(substitute(object))) )
  }
  return( object$coef )
}