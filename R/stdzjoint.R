#' Standardize the input matrix
#'
#' @param  X A matrix of any dimension
#' @return Returns a joint standardized matrix
#' @export

stdzjoint  <- function(X){
  Xbar <- apply(X, 2,mean)
  invSd <- matpower(cov(X),-1/2)
  Z <-  t(invSd %*% (t(X) - Xbar))
  return(Z)
}
