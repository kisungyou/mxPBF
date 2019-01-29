#' Two-Sample Mean Test
#'
#'
#' @param X an \eqn{(n_1\times p)} data matrix where each row is an observation.
#' @param Y an \eqn{(n_2\times p)} data matrix where each row is an observation.
#' @param a0 shape parameter for inverse-gamma prior.
#' @param b0 scale parameter for inverse-gamma prior.
#' @param gamma non-negative variance scaling parameter.
#'
#' @return a named list containing: \describe{
#' \item{log.BF.vec}{a length-\eqn{p} vector of log Bayes factors.}
#' }
#'
#' @examples
#' \dontrun{
#' ## generate X from standard, isotropic normal distribution
#' ## generate Y where the last column is only differed.
#' X = matrix(rnorm(200*5), nrow=200)
#' Y = cbind(matrix(rnorm(50*4),nrow=50),rnorm(50,mean=5))
#'
#' ## run test with different parameters
#' testmean2(X, Y)
#' testmean2(X, Y, a0=5.0, b0=5.0) # change some params
#' }
#'
#' @export
testmean2 <- function(X, Y, a0=2.0, b0=2.0, gamma=1.0){
  ###########################################################################
  # Preprocessing : Inputs
  # 1. data
  if (!check_datamatrix(X)){
    stop("* testmean2 : an input matrix X is invalid.")
  }
  if (!check_datamatrix(Y)){
    stop("* testmean2 : an input matrix Y is invalid.")
  }
  p  = ncol(X)
  if (ncol(Y)!=p){
    stop("* testmean2 : two inputs X and Y must be of same dimension.")
  }
  if ((nrow(X)==1)||(n2 = nrow(Y)==1)||(p==1)){
    stop("* testmean2 : invalid inputs.")
  }
  # 2. parameters
  extra.args = list(a0=a0, b0=b0, gamma=gamma)

  ###########################################################################
  # Main Computation
  output = testmean2.engine(X, Y, extra.args)
  return(output)
}


# main computation engine -------------------------------------------------
#' @keywords internal
#' @noRd
testmean2.engine <- function(X, Y, params){
  #################################################################
  ## Params
  parnames = names(params)
  if ("a0" %in% parnames){    a0 = params$a0  } else {    a0 = 2.0} ## a0
  if ("b0" %in% parnames){    b0 = params$b0  } else {    b0 = 2.0} ## a0
  if ("gamma" %in% parnames){gamma= params$gamma}else{ gamma = 1.0}
  ## Parameter Value Warning
  if ((length(a0)!=1)||(a0<=0)){      stop("* testmean2 : 'a0' should be nonnegative number.")}
  if ((length(b0)!=1)||(b0<=0)){      stop("* testmean2 : 'b0' should be nonnegative number.")}
  if ((length(gamma)!=1)||(gamma<=0)){stop("* testmean2 : 'gamma' should be nonnegative number.")}

  #################################################################
  ## Main Iteration
  n1 = nrow(X)
  n2 = nrow(Y)
  p  = ncol(X)
  n  = n1+n2

  log.BF.vec = rep(0,p)
  log.gammas = log(gamma/(1+gamma))
  for (j in 1:p){
    Xj = as.vector(X[,j])
    Yj = as.vector(Y[,j])
    Zj = as.vector(c(Xj,Yj))
    term1 = 2*b0 + n*get_mean2_sigmasq(Zj,gamma)
    term2 = 2*b0 + n1*get_mean2_sigmasq(Xj,gamma) + n2*get_mean2_sigmasq(Yj,gamma)
    log.BF.vec[j] = 0.5*log.gammas + ((n/2) + a0)*log(term1/term2)
  }

  output = list()
  output$log.BF.vec = log.BF.vec
  return(output)
}
