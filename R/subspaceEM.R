#' @title EM inference for a shared subspace
#'
#' @description
#' EM algorithm for shared subspace estimation.
#'
#' @param Ylist a list of n_k x p data matrices for k = 1... K groups
#' @param P number of features
#' @param S dimension of shared subspace of variation
#' @param R
#' @param Q
#' @param nvec (n_1, ... n_K)
#' @param Vstart Initialization
#' @param lambda
#' @param EM_iters Total number of EM iterations before terminating
#' @param M_iters  Total number of iterations in the numerical M-step
#' @param verbose
#'
#' @return
#'
#' @examples
#' 
#' @importFrom Matrix bdiag
#'
#' @export
subspaceEM <- function(Ylist, S, R=S, Q=S-R, Vstart=NULL, stiefelAlgo=2, 
                       lambda=0,
                       EM_iters=10,
                       M_iters=100,
                       verbose=FALSE) {


  P <- ncol(Ylist[[1]])
  
  if(missing(S))
    stop("Must provide a rank for the shared subspace estimation")
  
  nvec <- numeric(length(Ylist))
  for(k in 1:length(Ylist)) {
    
    if(ncol(Ylist[[k]]) != P)
      stop("Each data matrix must have the same number of columns")

    nvec[k] <- nrow(Ylist[[k]])
  }
  

  
  ## Function to optimize
  F <- function(V, PhiList, PrecVec) {
    obj <- 0
    for(k in 1:length(PhiList)) {
      VY <- t(V) %*% t(Ylist[[k]])
      obj <- obj +
        1/2 * tr(  VY %*% t(VY) %*% PhiList[[k]] ) -
        1/2 * PrecVec[k] * tr(VY %*% t(VY))
    }
    obj + lambda*sum(abs(V))
  }

  ## dF(X)/dX
  dF <- function(V, PhiList, PrecVec) {
    G <- 0
    for(k in 1:length(PhiList)) {

      G <- G + t(Ylist[[k]]) %*% ((Ylist[[k]] %*% V) %*% PhiList[[k]]) -
        (PrecVec[k] * t(Ylist[[k]])) %*%  (Ylist[[k]] %*% V)

    }
    G + lambda*sign(V)
  }


  if(is.null(Vstart)) {
    Vstart = rustiefel(P, S)
  }

  Vnew <-  Vstart

  convCheck <- Inf
  iter <- 0
  while(convCheck > 1e-6 & iter < EM_iters ) {
    ## ---------- E-step -----------

    ## E[ 1/sigma^2 * (psi+I)^(-1) | V]
    PrecVec <- rep(1, length(Ylist))
    PhiList <- list()
    V1 <- Vnew[, 1:(S-Q)]
    if(Q > 0) {
      V2 <- Vnew[, (S-Q+1):S]
      V2Yk <- lapply(1:length(Ylist), function(k) t(V2) %*% t(Ylist[[k]])*sqrt(PrecVec[k]))
      V2Ssum <- Reduce('+', lapply(V2Yk, function(x) x %*% t(x)))

      ## Diag Q  to ensure inversion is stable
      PhiShared <- solve(V2Ssum + 1e-5*diag(Q)) * (sum(nvec) + Q + 1)

    } else{
      V2 <- matrix(nrow=nrow(V1), ncol=0)
      PhiShared <- matrix(nrow=0, ncol=0)
    }

    for(k in 1:length(Ylist)) {

      ## Diag (S-Q)  to ensure inversion is stable
      V1Y <- t(V1) %*% t(Ylist[[k]])
      PsiK <- solve(V1Y %*% t(V1Y) + 1e-5*diag(S-Q)) * (nvec[k] + (S-Q) +1)
      PhiList[[k]] <- as.matrix(bdiag(PsiK, 1/PrecVec[k]*PhiShared))

    }

    ## E[ 1/sigma^2  | V]
    for(k in 1:length(Ylist)) {
      PrecVec[k] <- (nvec[k] * (P - S) + 2) /
        (sum(Ylist[[k]]^2) - tr((t(Vnew) %*% t(Ylist[[k]])) %*% (Ylist[[k]] %*% Vnew )))
    }

    ## Objective function to optimize in M-step
    F_t <- function(V) F(V, PhiList, PrecVec)
    dF_t <- function(V) dF(V, PhiList, PrecVec)

    ## ------- M-step -----------
    
    Vnew <- optStiefel(F_t, dF_t, Vinit=Vstart, method="bb",
                       maxIters= M_iters,
                       searchParams = NULL, verbose=verbose)

    ## ---- Check for convergence ------

    convCheck <- 1 - (norm(t(Vstart) %*% Vnew, type="F")/sqrt(S))
    Vstart <- Vnew
    iter <- iter + 1
    if(verbose) {
      print(PrecVec)
      print(convCheck)
    }

  }

  list(V=Vnew, PhiList=PhiList, PrecVec=PrecVec)
}


#' Title
#'
#' @return
#' @export
#'
#' @examples
subspaceOptim <- function(Ylist, S, R=S, Q=S-R, Vstart=NULL, stiefelAlgo=2, 
                          lambda=0,
                          maxIters=10000,
                          verbose=FALSE) {
  
  
  P <- ncol(Ylist[[1]])
  
  ###### Marginal Likelihood
  
  dF <- function(V) {
    
    YVlist <- lapply(Ylist, function(x) x %*% V)
    
    G <- 0
    for(k in 1:length(Ylist)) {
      
      VSV <- t(YVlist[[k]]) %*% YVlist[[k]] + 1e-5 * diag(S)
      VSVinv <- solve(VSV)
      
      G <- G +
        nvec[k]/2 * (VSVinv) %*% t(YVlist[[k]]) %*% Ylist[[k]]

    }
    t(G)
  }
  
  F<- function(V) {

    YVlist <- lapply(Ylist, function(x) x %*% V)
    obj <- 0
    for(k in 1:length(Ylist)) {
      VSV <- t(YVlist[[k]]) %*% YVlist[[k]] + 1e-5*diag(S)
      obj <- obj + nvec[k]/2 * log(det(VSV))  +
        nvec[k]*(P-S)/2 * log(sum(Ylist[[k]]^2) - tr(VSV))
    }
    obj
  }
  
  Vnew <- optStiefel(F, dF, Vinit=Vstart, 
                     method="bb",
                     maxIters= maxIters,
                     searchParams = NULL, 
                     verbose=verbose)
  
  Vnew
  
}



#' @title Subspace initialization
#'
#' @param Ylist 
#' @param S 
#'
#' @return
#' @export
#'
#' @examples
subspaceInit <- function(Ylist, S) {

  P <- ncol(Ylist[[1]])
  ngroups <- length(Ylist)
  nvec <- sapply(Ylist, function(Yk) nrow(Yk))

  isoVar <- sapply(1:ngroups, function(i) median(apply(Ylist[[i]], 2, var)))
  weightsList <- sapply(1:ngroups, function(i) {
    evals <- svd(Ylist[[i]]/sqrt(nvec[i]))$d^2
    dp <- suppressWarnings(
      ifelse(evals/isoVar[i] <= (1+sqrt(P/nvec[i])), 
             0, 
             sqrt((1-P/nvec[i]/(evals/isoVar[i]-1)^2) / (1+P/nvec[i]/(evals/isoVar[i]-1))))
    )
    weights <- 1/(1-dp) - 1
    weights[1:min(nvec[i], S)]
  })
  
  ## Weighted initialization
  Vinit <- svd(do.call(cbind, lapply(1:ngroups, function(k) {
    svdk <- svd(t(Ylist[[k]]))$u[, 1:min(nvec[k], S)] %*%  diag(weightsList[[k]][1:min(nvec[k], S)])
  })))$u[, 1:S]
  
  Vinit <- Vinit %*% rustiefel(S, S)
  Vinit
  
}
