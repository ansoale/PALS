#' compute directions using principal asymmetric least squares
#'
#' @param  X A design or model matrix which cannot include categorical variables
#' @param Y A continuous response variable
#' @param lambda A tuning parameter which must be greater than 0
#' @param tau The quantile levels which must be at least two
#' @return Returns the directions and their corresponding eigenvalues, and the coefficients of the expectile estimates
#' @export

LPALS <- function(X,Y,lambda=1,tau=seq(0.1,0.9,by=0.1)) {

  if(!require(KERE)) install.packages("KERE")
  library(KERE)


  n <- dim(X)[1]
  p <- dim(X)[2]
  h <- length(tau)

  Z <- stdzjoint(X)
  stdev <- matpower(cov(X),-1/2)

  kern <- vanilladot()
  coefs <- matrix(0,p,h)
  for(i in 1:h){
    model <- KERE(Z,Y,kern,lambda,omega=tau[i])
    eta <- apply(model$alpha[-1]%*%Z, 2, sum)
    coefs[,i] <- eta
  }
  beta <- solve(t(stdev))%*%coefs
  M <- beta%*%t(beta)
  eig <- eigen(M)

  return(list(dir = eig$vectors,eigval = eig$values,coefs=coefs))

}
