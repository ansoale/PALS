#' computes the projection matrix
#'
#' @param  X A matrix of any dimension
#' @return Returns a square matrix
#' @export


proj <- function(X){
  if(!class(X)=="matrix") X <- as.matrix(X)

  return(X%*%matpower(crossprod(X),-1)%*%t(X))
}
