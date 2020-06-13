# bases should be a vector with integers up to 100,008
# bases contains order numbers of prime numbers
# if bases is     NULL, the result is n x dim, where dim = n.bases
# if bases is not NULL, the result is n x dim, where dim = length(bases)
# if bases is not NULL, argument n.bases is ignored

halton <- function(n = 1, n.bases = 1, bases = NULL, 
                   start = 0, random.primes = FALSE, seed = 7, 
                   scale.coverage = FALSE, shuffle = FALSE){
  set.seed(seed)
  if(is.null(bases)){
    # cat.print(1)
    if(!is.numeric(n.bases)){
      stop("'n.bases' must be numeric")
    }
    if(length(n.bases) != 1)(
      stop("Length of 'n.bases' must be 1")
    )
    dim <- n.bases
    bases <- 1:n.bases
  } else {
    # cat.print(2)
    if(!is.numeric(bases)){
      stop("'bases' must be numeric")
    }
    bases <- as.integer(bases)
    dim <- length(bases)
  }
  # cat.print(bases)
  tymch <- matrix(nrow = n, ncol = dim)
  # cat.print(dim(tymch))
  my100008Primes <- npsf::primes(which = bases, random.primes = random.primes, seed = seed)
  if(length(n.bases) == 1 & random.primes){
    my100008Primes <- npsf::primes(n = n.bases, random.primes = random.primes, seed = seed)
  }
  # cat.print(my100008Primes)
  for(i in 1:dim){
    # cat.print(as.integer(start:(n-1+start)))
    # cat.print(as.integer(n))
    # cat.print(as.integer(my100008Primes[i]))
    # cat.print(double(n))
    # cat.print(bases[i])
    tymch0 <- .C("HaltonSeq", 
                 indices = as.integer(start:(n-1+start)), 
                 howmany = as.integer(n), 
                 thisbase = as.integer(my100008Primes[i]), 
                 Hdr = double(n))[[4]]
    # cat.print(tymch0)
    tymch[,i] <- tymch0
    if(shuffle){
      
      tymch[,i] <- tymch[,i][sample.int(n)]
    }
    if(scale.coverage){
      mymin <- min(c(tymch[,i], 1/n, 1-max(tymch[,i])))
      tymch[,i] <- npsf::rescale(tymch[,i], mymin, 1-mymin)
    }
  }
  # class(tymch) <- "npsf"
  return(tymch)
}