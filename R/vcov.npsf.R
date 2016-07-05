vcov.npsf <- function( object, ... ) {
 if(is.null(object$vcov) | is.null(object)){
  stop( gettextf("There is no variance-covariance matrix in '%s'.", deparse(substitute(object))) )
 }
 return( object$vcov )
}