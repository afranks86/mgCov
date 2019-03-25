##############################################################
### Functions for Sampling
#############################################################

#' @title Sample the variance
#'
#' @description Sample isotropic variance from complete conditional
#'
#' @param Y n x p data matrix
#' @param YV  Y x V where V is a basis a p x s basis for shared subspace
#' @param O s x s matrix of eigenvectors
#' @param omega
#' @param n
#' @param nu0 prior params for inv-chisq
#' @param s20 prior params for inv-chisq
#'
#' @return single sample of varaince
#'
#' @examples
#'
#'
#' @export
sampleSigma2 <- function(Y, YV, O, omega, n, nu0=1, s20=1) {

  YU <- YV %*% O

  p <- ncol(Y)
  a <- ( nu0 + n*p)/2

  b <- ( nu0*s20 + sum(Y^2) - tr(t(YU) %*% YU %*% diag(omega, nrow=length(omega))))/2

  1/rgamma(1, a, b)
}

#' @title Sample from a truncated gamma
#'
#' @description Sample from a gamma(a, b) trucated at 1
#'
#' @param a gamma shape param
#' @param b gamma rate param
#' @param scale
#'
#' @return
#' @export
#'
#' @examples
tgamma <- function(a, b, scale) {

  logM <- -1 * pgamma(scale, a, b, log.p=TRUE)
  qgamma(exp(log(runif(length(a))) - logM), a, b)

}

#' @title Sample omega, lambda/(lambda-1)
#'
#' @param YV  Y x V where V is a basis a p x s basis for shared subspace
#' @param O s x s matrix of eigenvectors
#' @param s2 isotropic variance
#' @param n number of samples
#'
#' @return omega, a sample omega
#'
#' @examples
#'
#' @export
sampleOmega <- function(YV, O, s2, n) {

  R <- ncol(O)

  Rtilde <- min(n, R)
  O <- O[, 1:Rtilde]

  YVO <- YV %*% O

  cvec <- diag(t(YVO) %*% YVO) /(2*s2)

  ## add divide by n for stability?
  g <- tgamma(rep(n/2 + 1, length(cvec)), cvec/n, n)
  g <- g/n

  omega <- 1 - g

  if(n < R)
    omega <- c(omega, rep(0, R - n))

  omega

}

#' @title Sample the eigenvectors
#'
#' @description
#'
#' @param YV  Y x V where V is a basis a p x s basis for shared subspace
#' @param O s x s matrix of eigenvectors
#' @param s2 isotropic variance
#' @param omega omega
#'
#' @return
#'
#' @examples
#'
#' @export
sampleO <- function(YV, O, s2, omega) {

  A <- t(YV) %*% YV / (2*s2)
  B <- diag(omega, nrow=length(omega))

  ord <- order(omega, decreasing=TRUE)
  revOrd <- order(ord)

  Btilde <- B[ord, ord]
  Otilde <- O[, ord]

  O <- rbing.matrix.gibbs(A, Btilde, Otilde)

  O[, revOrd]
}
