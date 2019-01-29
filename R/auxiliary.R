# auxiliary functions -----------------------------------------------------
# 1. check_datamatrix  : check whether input data is a valid matrix.
# 2. check_sqmat       : check if square matrix
# 3. check_pd          : (at least) positive semidefinite
# 4. get_invroot       : compute an appropriate scaler matrix.
# 5. get_mean2_sigmasq : weird computation
# 6. get_covv2_tau     : weird computation, again



# 1. check_datamatrix -----------------------------------------------------
#' @keywords internal
#' @noRd
check_datamatrix <- function(A){
  cond1 = (is.matrix(A)||(inherits(A, "Matrix")))
  cond2 = (!any(is.infinite(A)))
  cond3 = (!any(is.na(A)))
  if (cond1&&cond2&&cond3){
    return(TRUE)
  } else {
    return(FALSE)
  }
}
# 2. check_sqmat ----------------------------------------------------------
#' @keywords internal
#' @noRd
check_sqmat <- function(A){
  cond1 = (is.matrix(A)||(inherits(A, "Matrix")))
  cond2 = ((nrow(A)==ncol(A)))
  if (cond1&&cond2){
    return(TRUE)
  } else{
    return(FALSE)
  }
}
# 3. check_pd -------------------------------------------------------------
#' @keywords internal
#' @noRd
check_pd <- function(A){
  eigens = eigen(A, symmetric=TRUE, only.values=TRUE)
  if (any(eigens$values<0)){
    return(FALSE)
  } else {
    return(TRUE)
  }
}

# 4. get_invroot ----------------------------------------------------------
get_invroot <- function(X){
  eigs = eigen(X)
  if (any(eigs$values < .Machine$double.eps*10)){
    stop("** The desired covariance 'Sigma0' is invalid.")
  }
  out = eigs$vectors %*% diag((eigs$values)^(-0.5)) %*% t(eigs$vectors)
  return(out)
}


# 5. get_mean2_sigmasq ----------------------------------------------------
# Given a vector, compute (1/nx)*x'(I-(1/(1+gamma))*H_1)*x
#' @keywords internal
#' @noRd
get_mean2_sigmasq <- function(x, gamma){
  nx  = length(x)

  sx2 = sum(x*x)         # <x,x> = |x|^2
  vv  = rep(1,nx)
  Hx  = outer(vv,vv)/sum(vv*vv)   # projection matrix

  return((sx2 - sum(as.vector(Hx%*%x)*x)/(1+gamma))/nx)
}

# 6. get_covv2_tau --------------------------------------------------------
#' @keywords internal
#' @noRd
get_covv2_tau <- function(xi, xj, gamma){
  nx = length(xi)
  Hj = outer(xj,xj)/sum(xj*xj)

  return((sum(xi*xi) - sum(as.vector(Hj%*%xi)*xi)/(1+gamma))/nx)
}
