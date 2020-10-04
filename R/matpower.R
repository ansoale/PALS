#' computes the matrix power
#'
#' @param  A A square matrix
#' @param n The power of the matrix we want to compute
#' @return Returns a square matrix
#' @export

matpower <- function(A,n){
  A = round((A+t(A))/2,7)
  eig = eigen(A)
  return(eig$vectors%*%diag((eig$values)^n)%*%t(eig$vectors))
}
