print.npsf <- function( x, digits = NULL, ... ) {

   if( is.null( digits ) ) {
      digits <- max( 3, getOption( "digits" ) - 2 )
   }
   cat( "\nCall:\n" )
   cat( deparse( x$call ) )
   cat( "\n\n" )
   if(x$model == "sfsc"){
    if(!is.null(x$coef)){
     # printCoefmat(x$table, digits = digits)
     print.default(coef(x), digits = digits)
    }
   }
   invisible( x )
}