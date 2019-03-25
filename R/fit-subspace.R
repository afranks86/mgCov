## P=number of features, S=subspace dimension
## psi_k has dimension R + Q <= S, where Q eigenvectors are common
## and R eigenvectors vary across groups.
##
## Ylist: a list of n x p data matrices
## init: a list containing initial V, U's, Omega's,
## nvec: vector of number of samples per group

#' @title Bayesian inference for the spiked covariance model
#'
#' @description Run the Bayesian inference fo
#'
#' @param Ylist a list of n_k x p data matrices for k = 1... K groups
#' @param S dimension of the shared subspace (default)
#' @param R
#' @param Q
#' @param niters number of iterations to
#' @param nskip
#' @param init
#' @param verbose
#' @param sigmaTruthList
#' @param draw
#' @param printLoss
#'
#' @return
#' @export
#'
#' @examples
#'
#' @export
fitBayesianSpike <- function(V, Ylist,
                            R=S, Q=S-R,
                            niters=100, nskip=1, init=NULL,
                            verbose=TRUE, sigmaTruthList=NULL, draw=c(),
                            printLoss=FALSE) {

    ngroups = length(Ylist)
    nvec = sapply(Ylist, nrow)

    niters <- niters
    nskip <- nskip
    nsamps <- floor(niters/nskip)

    S <- ncol(V)
    P <- ncol(Ylist[[1]])

    for(k in 1:length(Ylist)) {
        if(ncol(Ylist[[k]]) != P)
          stop("All matrices in Ylist must have the same number of columns")
    }

    if(R + Q > S)
        stop("R + Q must be less than S")

    Osamps <- array(dim=c(S, R + Q, ngroups, nsamps))
    omegaSamps <- array(dim=c(R + Q, ngroups, ncol=nsamps))
    s2samps <- matrix(nrow=ngroups, ncol=nsamps)

    initItems <- c("V", "Ulist", "OmegaList")
    if(is.null(V)) {
      stop("Need to specify V")
    }

    if( is.null(init) ) {

        OmegaList <- Olist <- Ulist <- list()

        SigInvInitList <- list()
        s2vec <- rexp(ngroups)
        for(k in 1:ngroups) {
            OmegaList[[k]] <- rep(1/2, (R+Q))

            Ok <- matrix(0, nrow=S, ncol=(R+Q))
            Ok[1:R, 1:R] <- diag(R)
            if(Q > 0)
                Ok[(S-Q+1):S, (R+1):(R+Q)] <- diag(Q)

            Olist[[k]] <- Ok
            Ulist[[k]] <- V %*% Ok

            SigInvInitList[[k]] <-
                getSigmaInv(P, Ulist[[k]], OmegaList[[k]], s2vec[k])
        }


        init <- list(V=Vinit, Ulist=Ulist, OmegaList=OmegaList, s2vec=s2vec)

    } else if (any(!(initItems %in% names(init)))) {
        stop("Some parameters unintialized!")
    }

    Ulist <- init$Ulist
    Olist <- lapply(Ulist, function(u) t(V) %*% u)
    OmegaList <- init$OmegaList
    s2vec <- init$s2vec

    YVlist <- list()
    for(k  in 1:ngroups) {
        omega <- OmegaList[[k]]
        Lambda <- diag((omega/(1-omega)))

        YVlist[[k]] <- Ylist[[k]] %*% V
    }


    draw <- c(draw, O=TRUE, s2=TRUE, omega=TRUE)
    draw <- draw[unique(names(draw))]

    for ( i in 1:niters ) {

        if(Q > 0) {

            YVpooled <- lapply(1:ngroups, function(k) YVlist[[k]]/sqrt(s2vec[k])) %>%
                Reduce('rbind', .)

            ## common omega
            Omega2 <- sampleOmega(YVpooled[, (R+1):(R+Q)], Olist[[k]][(R+1):(R+Q), (R+1):(R+Q)],
                                  1, sum(nvec))

            ## Sample common O
            O2 <- sampleO(YVpooled[, (R+1):(R+Q)], Olist[[k]][(R+1):(R+Q), (R+1):(R+Q)], 1,
                          Omega2)

            ## order eigenvalues from largest to smallest
            ord_omega <- order(Omega2, decreasing=TRUE)
            Omega2 <- Omega2[ord_omega]
            O2 <- O2[, ord_omega]

        } else {
            Omega2 <- c()
            O2 <- matrix(nrow=0, ncol=0)
        }

        for ( k in 1:ngroups ) {

            ## Sample sigma^2_k
            if(draw["s2"])
                s2vec[k] <- sampleSigma2(Ylist[[k]], YVlist[[k]], Olist[[k]], OmegaList[[k]], nvec[k])


            ## Sample omegas_k, don't requrie ordered
            if(draw["omega"]) {

                Omega1 <- sampleOmega(YVlist[[k]], Olist[[k]][, 1:R],
                                         s2vec[k], nvec[k])



            }
                ## Sample O_k
            if(draw["O"]) {

                O1 <- sampleO(YVlist[[k]][, 1:R], Olist[[k]][1:R, 1:R], s2vec[k],
                                  (OmegaList[[k]])[1:R])

                ## order eigenvalues from largest to smallest
                ord_omega <- order(Omega1, decreasing=TRUE)
                Omega1 <- Omega1[ord_omega]
                O1 <- O1[, ord_omega]

                Ok <- as.matrix(bdiag(O1, O2))
            }

            OmegaList[[k]] <- c(Omega1, Omega2)

            Olist[[k]] <- Ok

            ## save samples
            if(i %% nskip==0) {
                Osamps[, , k, i/nskip] <- Olist[[k]]
                omegaSamps[, k, i/nskip] <- OmegaList[[k]]
                s2samps[k, i/nskip] <- s2vec[k]
            }

        }

        if (i %% nskip==0) {
            if(verbose & !printLoss)
                print(sprintf("Iteration %i", i))
        }

        if(i %% nskip==0 & verbose & printLoss) {
            stop("Needs updating")
            sl <- 0
            if(is.null(sigmaTruthInvList)) {
                ## Print loss relative to starting point
                for(k in 1:ngroups) {

                    omegak <- OmegaList[[k]]
                    SigHat <- s2vec[k]*(Ulist[[k]] %*%
                                        diag(omegak/(1-omegak)) %*%
                                        t(Ulist[[k]])+diag(P))

                    SigmaStartInvK <- with(initSS,
                                           1/s2vec[k]*(diag(P) - Ulist[[k]] %*% diag(OmegaList[[k]]) %*% t(Ulist[[k]])))

                    sl <- sl + steinsLoss(SigHat, SigmaStartInvK)
                }
            } else {
                for(k in 1:ngroups) {
                    omegak <- OmegaList[[k]]
                    SigHat <- s2vec[k]*(Ulist[[k]] %*%
                                        diag(omegak/(1-omegak)) %*%
                                        t(Ulist[[k]])+diag(P))

                    sl <- sl + steinsLoss(SigHat, sigmaTruthInvList[[k]])

                }
            }
            print(sprintf("Iteration %i, Loss %f", i, sl/ngroups))
        }
    }

    list(S=S, R=R, Q=Q, V=V, init=init, ngroups=ngroups, Osamps=Osamps,
         omegaSamps=omegaSamps, s2samps=s2samps)
}

