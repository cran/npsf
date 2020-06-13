# is

primes <- function(n = NULL, which = NULL, random.primes = FALSE, seed = 7){
  # https://primes.utm.edu/lists/small/100000.txt
  # https://primes.utm.edu
  
  # none supplied
  if(is.null(n) & is.null(which)){
    stop("Either 'n' or 'which' must be supplied")
  }
  # both supplied
  if(!is.null(n) & !is.null(which)){
    stop("Both 'n' and 'which' cannot be supplied")
  }
  # n supplied
  if(!is.null(n)){
    if(!is.numeric(n)){
      stop("'n' must be an numeric")
    } else {
      # cat.print(1)
      n <- as.integer(n)
      if(length(n) == 1){
        # cat.print(2)
        if(n > 100008){
          stop("'n' should be smaller than 100,008")
        } else {
          # cat.print(3)
          if(random.primes){
            # cat.print(4)
            set.seed(seed)
            which2choose <- sample.int(100008, size = n)
          } else {
            # cat.print(5)
            which2choose <- 1:n
          }
        }
      } else {
        stop("'n' must be one integer")
      }
    }
  } else {
    # which supplied
    if(!is.numeric(which)){
      stop("'which' must be an numeric")
    } else {
      # cat.print(6)
      which2choose <- as.integer(which)
    }
  }
  

  # if(!is.numeric(which)){
  #   stop("'which' must be an numeric")
  # } else {
  #   which <- as.integer(which)
  #   if(length(which) == 1){
  #     if(which > 100008){
  #       stop("'which' should be smaller than 100,008")
  #     } else {
  #       if(random.primes){
  #         set.seed(seed)
  #         which2choose <- sample.int(100008, size = which)
  #       } else {
  #         which2choose <- 1:which
  #       }
  #     }
  #   } else {
  #     which2choose <- which
  #   }
  # }
  # cat.print(which2choose)
  n <- length(which2choose)
  tymch <- .C("Primes", as.integer(which2choose-1), as.integer(n), double(n))[[3]]
  # class(tymch) <- "npsf"
  return(tymch)
}