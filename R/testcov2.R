#' Two-Sample Covariance Test
#'
#' @param X an \eqn{(n_1\times p)} data matrix where each row is an observation.
#' @param Y an \eqn{(n_2\times p)} data matrix where each row is an observation.
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
#' ## generate X from standard, isotropic normal distribution
#' ## generate Y where the last column is only differed.
#' X = matrix(rnorm(200*5), nrow=200)
#' Y = cbind(matrix(rnorm(50*4),nrow=50),rnorm(50,sd=5))
#'
#' ## run test with different parameters
#' testcov2(X, Y)
#' testcov2(X, Y, a0=5.0, b0=5.0) # change some params
#' }
#'
#' @export
testcov2 <- function(X, Y, a0=2.0, b0=2.0, gamma=1.0){
  ###########################################################################
  # Preprocessing : Inputs
  # 1. data
  if (!check_datamatrix(X)){
    stop("* testcov2 : an input matrix X is invalid.")
  }
  if (!check_datamatrix(Y)){
    stop("* testcov2 : an input matrix Y is invalid.")
  }
  p  = ncol(X)
  if (ncol(Y)!=p){
    stop("* testcov2 : two inputs X and Y must be of same dimension.")
  }
  if ((nrow(X)==1)||(n2 = nrow(Y)==1)||(p==1)){
    stop("* testcov2 : invalid inputs.")
  }
  # 2. parameters
  extra.args = list(a0=a0, b0=b0, gamma=gamma)

  ###########################################################################
  # Main Computation
  new.X  = as.matrix(scale(X, center=TRUE, scale=FALSE))
  new.Y  = as.matrix(scale(Y, center=TRUE, scale=FALSE))
  output = testcov2.engine(new.X, Y, extra.args)
  return(output)
}


# main computation engine -------------------------------------------------
#' @keywords internal
#' @noRd
testcov2.engine <- function(X, Y, params){
  #################################################################
  ## Params
  parnames = names(params)
  if ("a0" %in% parnames){    a0 = params$a0  } else {    a0 = 2.0} ## a0
  if ("b0" %in% parnames){    b0 = params$b0  } else {    b0 = 2.0} ## a0
  if ("gamma" %in% parnames){gamma= params$gamma}else{ gamma = 1.0}
  ## Parameter Value Warning
  if ((length(a0)!=1)||(a0<=0)){      stop("* testcov2 : 'a0' should be nonnegative number.")}
  if ((length(b0)!=1)||(b0<=0)){      stop("* testcov2 : 'b0' should be nonnegative number.")}
  if ((length(gamma)!=1)||(gamma<=0)){stop("* testcov2 : 'gamma' should be nonnegative number.")}

  #################################################################
  ## Main Iteration
  n1 = nrow(X)
  n2 = nrow(Y)
  p  = ncol(X)
  n  = n1+n2

  log.BF.mat = matrix(rep(0,p*p),nrow=p)
  term.const = 0.5*log(gamma/(1+gamma)) + lgamma(n1/2 + a0) + lgamma(n2/2 + a0) - lgamma(n/2 + a0) + a0*log(b0) - lgamma(a0)
  for(i in 1:p){
    Xi = as.vector(X[,i])
    Yi = as.vector(Y[,i])
    Zi = as.vector(c(Xi,Yi))
    for(j in (1:p)[-i]){
      Xj = as.vector(X[,j])
      Yj = as.vector(Y[,j])
      Zj = as.vector(c(Xj,Yj))

      term1 = (n1/2 + a0)*log(b0 + (n1/2)*get_covv2_tau(Xi,Xj,gamma))
      term2 = (n2/2 + a0)*log(b0 + (n2/2)*get_covv2_tau(Yi,Yj,gamma))
      term3 = (n/2  + a0)*log(b0 + (n/2)*get_covv2_tau(Zi,Zj,gamma))

      log.BF.mat[i,j] = term.const-term1-term2+term3
    }
  }
  diag(log.BF.mat) = -Inf

  output = list()
  output$log.BF.mat = log.BF.mat
  return(output)
}
