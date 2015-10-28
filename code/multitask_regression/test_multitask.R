require(mvtnorm)

source("multitask_regression.R")

N=50
P=2
K=300

mu=c(-1,1)
sigma=diag(P)
noiseSD=.3

x=list()
y=list()
b=list()
for (k in 1:K){
  x[[k]]=cbind(runif(N),1)
  b[[k]]=rmvnorm(1,mu,sigma)
  y[[k]]=as.numeric(x[[k]] %*% t(b[[k]]) + noiseSD * rnorm(N))
}

fitted <- multitask_regression(x,y)

fitted$par$coeffs

b=t(do.call(rbind,b))

qplot( as.numeric(b), as.numeric(fitted$par$coeffs) ) + theme_bw(base_size = 16) +xlab("true coefficient") + ylab("inferred coefficient")
