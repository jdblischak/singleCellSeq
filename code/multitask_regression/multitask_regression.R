#script.dir <- dirname(sys.frame(1)$ofile)

require(rstan)
require(doMC)

sm=stan_model("../code/multitask_regression/multitask_regression.stan")

cholsolve <- function(a,b) {
  ch <- chol(a)
  backsolve( ch, backsolve( t(ch), b, upper.tri=F ), upper.tri = T)
}

eigeninv <- function(a) {
  s <- svd(a)
  s$d <- pmax(s$d,1e-4)
  s$v %*% diag(1/s$d) %*% t(s$u)
}

multitask_regression <- function(x,y) {
  K <- length(x)
  N <- length(y[[1]])
  P <- ncol(x[[1]])
  
  dat=list(N=N,P=P,K=K,x=x,y=y)
  
  o <- optimizing(sm, dat, as_vector=F)
  cat("Done\n")
  precSigma=eigeninv(o$par$Sigma)
  o$par$coeffs <- do.call(cbind,foreach(k=1:K) %dopar% 
    #solve( precSigma + t(x[[k]]) %*% x[[k]] / o$par$noiseVar, precSigma %*% o$par$mu + t(x[[k]]) %*% y[[k]] / o$par$noiseVar)
    #cholsolve( precSigma + t(x[[k]]) %*% x[[k]] / o$par$noiseVar, precSigma %*% o$par$mu + t(x[[k]]) %*% y[[k]] / o$par$noiseVar)
      { eigeninv( precSigma + t(x[[k]]) %*% x[[k]] / o$par$noiseVar ) %*% ( precSigma %*% o$par$mu + t(x[[k]]) %*% y[[k]] / o$par$noiseVar) } 
    )
  o
}