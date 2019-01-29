#' One-Sample Covariance Test
#'
#' @param data an \eqn{(n\times p)} data matrix where each row is an observation.
#' @param Sigma0 a \eqn{(p\times p)} given covariance matrix.
#' @param a0 shape parameter for inverse-gamma prior.
#' @param b0 scale parameter for inverse-gamma prior.
#' @param gamma non-negative variance scaling parameter.
#'
#' @return a named list containing: \describe{
#' \item{log.BF.mat}{a \eqn{(p\times p)} matrix of pairwise log Bayes factors.}
#' }
#'
#' @examples
#' \dontrun{
#' ## generate data from multivariate normal with trivial covariance.
#' data = matrix(rnorm(100*5), nrow=100)
#'
#' ## run test with different parameters
#' out1 = testcov1(data)
#' out2 = testcov1(data, a0=5.0, b0=5.0) # change some params
#'
#' ## visualize two Bayes Factor matrices
#' par(mfrow=c(1,2), pty="s")
#' image(exp(out1$log.BF.mat)[,5:1], main="default")
#' image(exp(out2$log.BF.mat)[,5:1], main="a0=b0=5.0")
#' }
#'
#'@export
testcov1 <- function(data, Sigma0=diag(ncol(data)), a0=2.0, b0=2.0, gamma=1.0){
  ###########################################################################
  # Preprocessing : Inputs
  # 1. data
  if (!check_datamatrix(data)){
    stop("* testcov1 : an input matrix is invalid.")
  }
  n = nrow(data)
  p = ncol(data)
  if ((nrow(data)==1)||(ncol(data)==1)){
    stop("* testcov1 : invalid input matrix X.")
  }
  # 2. valid Sigma0
  if ((!check_sqmat(Sigma0))||(!check_pd(Sigma0))||(!isSymmetric(Sigma0,tol=sqrt(.Machine$double.eps)))||(nrow(Sigma0)!=p)){
    stop("* testcov1 : a given matrix for null hypothess 'Sigma0' is invalid.")
  }
  # 3. extra arguments
  extra.args = list(a0=a0, b0=b0, gamma=gamma)

  ###########################################################################
  # Preprocessing : Adjust the data for testing I
  scaler = get_invroot(Sigma0)
  X.centered = scale(data, center=TRUE, scale=FALSE)
  X.adjusted = (matrix(X.centered,nrow=n) %*% scaler)

  ###########################################################################
  # Main Computation and Return
  output = testcov1.engine(X.adjusted, extra.args)
  return(output)
}




# main computation engine -------------------------------------------------
#' @keywords internal
#' @noRd
testcov1.engine <- function(X, params){
  #################################################################
  ## Params
  parnames = names(params)
  if ("a0" %in% parnames){    a0 = params$a0  } else {    a0 = 2.0} ## a0
  if ("b0" %in% parnames){    b0 = params$b0  } else {    b0 = 2.0} ## a0
  if ("gamma" %in% parnames){gamma= params$gamma}else{ gamma = 1.0}
  ## Parameter Value Warning
  if ((length(a0)!=1)||(a0<=0)){      stop("* testcov1 : 'a0' should be nonnegative number.")}
  if ((length(b0)!=1)||(b0<=0)){      stop("* testcov1 : 'b0' should be nonnegative number.")}
  if ((length(gamma)!=1)||(gamma<=0)){stop("* testcov1 : 'gamma' should be nonnegative number.")}

  #################################################################
  ## MAIN RUN BY JAY
  p = ncol(X)
  n = nrow(X)
  log.BF.mat = matrix(0, ncol=p,nrow=p) # log Bayes factors

  for(i in 1:p){
    Xi   = matrix(X[,i], ncol=1)
    sXi2 = sum((Xi)^2)
    for(j in (1:p)[-i]){
      Xj   = matrix(X[,j], ncol=1)
      sXj2 = sum((Xj)^2)
      log.BF.mat[i,j] = a0*log(b0) - lgamma(a0) +
        1/2 * log(gamma/(1+gamma)) + lgamma(n/2 + a0) +
        (1/2 * sXi2) - (n/2 + a0) * log(1/2*(sXi2-sum(Xi*Xj)^2/(sXj2*(1+gamma))) + b0)
    }
  }
  diag(log.BF.mat) = -Inf # just to fill out the diagonal parts
  output = list()
  output$log.BF.mat = log.BF.mat
  return(output)
}
