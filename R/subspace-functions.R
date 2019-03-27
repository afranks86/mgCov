library(mvtnorm)
library(magrittr)
library(scales)
library(car)

## Compute frobenius norms
## H = Hierachical eigen pooling
## e.g. compare trace (U_k^T U_J)^2
## SS = Shares space
## e.g. compare || U_k^T U_J || (frobenius)

createNormMatrix <- function(eigenlist, vindices, type="H") {
    ngroups <- length(eigenlist)
    normsMat <- matrix(NA, nrow=ngroups, ncol=ngroups)

    for( i in 1:ngroups ) {
        for( j in i:ngroups ){
            Ui <- eigenlist[[i]]$vectors[, vindices]
            Uj <- eigenlist[[j]]$vectors[, vindices]
            if(type == "H") {
                normsMat[i, j] <- sum((t(Ui)%*%Uj)^2)/sqrt(ncol(Ui)) ## these are the same
            } else if (type == "SS") {

                normsMat[i, j] <- norm(t(Ui)%*%Uj, type="F")/sqrt(ncol(Ui))
            }
        }
    }
    diag(normsMat) <- NA
    normsMat[lower.tri(normsMat)] <- t(normsMat)[lower.tri(normsMat)]

    return(normsMat)
}

#' Title TODO
#'
#' @param Osamps
#' @param OmegaSamps
#' @param s2samps
#' @param nsamps
#' @param groups
#' @param probRegion
#' @param hline
#' @param col
#' @param pch
#' @param lty
#' @param ymax
#' @param type
#' @param plotPoints
#' @param polar
#' @param cex.axis
#' @param cex.pts
#' @param splitGroups
#' @param splitPts
#'
#' @return
#' @export
#'
#' @examples
posteriorPlot <- function(Osamps, OmegaSamps, s2samps, nsamps, groups_to_plot,
                          probRegion=0.95, hline=NULL,  ymax=NULL, type = "mag",
                          plotPoints=TRUE, polar=FALSE) {

  ngroups <- length(groups_to_plot)

  group_names <- names(groups_to_plot)
  if(is.null(group_names))
    group_names = factor(1:ngroups)

  if(type=="mag") {
    ylab <- expression(lambda[1])
  }
  else if(type=="logmag") {
    ylab <- expression("log"[2]~"(" ~ lambda[1] ~ ")")
  } else if (type =="logratio") {
    ylab <- expression("log"[2]~"(" ~ lambda[1]/lambda[2] ~ ")")
  } else {
    ylab <- expression("(" ~ lambda[1]/lambda[2] ~ ")")
  }

  if(is.null(ymax))
    ymax <- 1.1*max(OmegaSamps/(1-OmegaSamps))

  ##plot(0, 0, xlim=c(-pi/2, pi/2),
  ##     ylim=c(0, ymax), cex=0, xlab=, ylab=ylab, xaxt="n", cex.axis=cex.axis, cex.lab=1.5)
  ## axis(1, at=seq(-pi/2, pi/2, by=pi/4), labels=expression(-pi/2, -pi/4, 0, pi/4, pi/2), cex.axis=cex.axis, cex.lab=1.5)

  group_pts_angle <- group_pts_eval <- group_type <- c()
  count <- 1
  for(g in groups_to_plot) {

    pmPsi <- getPostMeanPsi(Osamps[, , g , ], OmegaSamps[, g, ],
                            s2samps[g, ], nsamps)

    eigPsi <- eigen(pmPsi)
    pmValues <- eigPsi$values
    pmVectors <- eigPsi$vectors
    maxIndex <- which.max(pmValues)

    hp <- getHullPoints(nsamps, OmegaSamps[, g, ], Osamps[, , g, ],
                        type=type, probRegion=probRegion)
    pts <- hp$pts
    hullPoints <- hp$hullPoints

    group_pts_angle <- c(group_pts_angle, pts[1, ])
    group_pts_eval <- c(group_pts_eval, pts[2, ])
    group_type <- c(group_type, rep(group_names[count], length(pts[2, ])))
    count <- count + 1

  }

  posterior_summaries <- tibble(angle = group_pts_angle, eval = group_pts_eval, Group = factor(group_type))

  p <- ggplot(posterior_summaries) +
    geom_point(aes(x=angle, y=eval, col=Group)) +
    theme_bw(base_size=20) + xlim(c(-pi/2, pi/2)) + ylim(c(0, ymax))

  if(!is.null(hline))
    p <- p + geom_hline(yintercept=hline, lty=2)


  p + ylab(ylab) + xlab(expression("angle, acos("~U[1]^T*V[1]~")"))

}


posteriorPlotPolar <- function(Osamps, OmegaSamps, s2samps, nsamps, groups,
                          probRegion=0.95, hline=NULL, col=NULL,
                          pch=NULL, lty=NULL, ymax=NULL, logRatio=FALSE,
                          plotPoints=TRUE, polar=FALSE, cex.axis=1.5) {

    ngroups <- length(groups)

    if(is.null(col)){
        col <- 1:ngroups
    }
    if(is.null(pch)){
        pch=rep(19, ngroups)
    }
    if(is.null(lty)){
        lty=rep(1, ngroups)
    }

    radius_vec <- angle_vec <- type_vec <- c()

    for(g in groups) {

        pmPsi <- getPostMeanPsi(Osamps[, , g , ], OmegaSamps[, g, ],
                                s2samps[g, ], nsamps)

        eigPsi <- eigen(pmPsi)
        pmValues <- eigPsi$values
        pmVectors <- eigPsi$vectors
        maxIndex <- which.max(pmValues)

        pmPoint <- pmValues[maxIndex]*pmVectors[, maxIndex]
        if(pmPoint[1] < 0)
            pmPoint <- pmPoint * c(-1, -1)

        hp <- getHullPoints(nsamps, OmegaSamps[, g, ], Osamps[, , g, ],
                            logRatio=logRatio)
        pts <- hp$pts

        radius_vec <- c(radius_vec, pts[2, ], pts[2, ])
        angle_vec <- c(angle_vec, pts[1, ], pts[1, ]-pi)
        type_vec <- c(type_vec, rep(g, 2*length(pts[1,])))

    }

    tibble(radius=radius_vec, angle=angle_vec, type=type_vec) %>% ggplot(aes(y=radius, x=angle, col=as.factor(type))) +
        geom_point() + coord_polar(theta="x", clip="off") + ylim(0, 2) + xlim(c(-pi, pi)) + theme_bw() + scale_color_discrete_qualitative(alpha=0.3)

}


## type is "mag", "ratio",
getHullPoints <- function(nsamps, OmegaSamps, Osamps, type="mag",
                          probRegion=0.95) {

    PointsList <- lapply(1:nsamps, function(i) {
        LambdaSamp <- OmegaSamps[, i]/(1-OmegaSamps[, i])
        maxIndex <- which.max(LambdaSamp)

        if(type == "mag")
            yval <- LambdaSamp[maxIndex]
        else{
            yval <- LambdaSamp[maxIndex]/LambdaSamp[-maxIndex]
        }

        O1 <- Osamps[, maxIndex, i]
        angle <- atan(O1[2]/O1[1])

        c(angle, yval)

    })

    pts <- simplify2array(PointsList)
    allPts <- pts
    if(type == "logratio") {
        allPts[2, ] <- log2(allPts[2, ])
        pts[2, ] <- log2(pts[2, ])
    }

    numPtsToRemove <- round(nsamps*(1-probRegion))
    while(numPtsToRemove > 0) {
        hullPoints <- chull(pts[1, ], pts[2, ])
        if(length(hullPoints) > numPtsToRemove) {
            hullPoints <- sample(hullPoints, numPtsToRemove)
            pts <- pts[, -hullPoints]
            numPtsToRemove <- 0
        } else{
            pts <- pts[, -hullPoints]
            numPtsToRemove <- numPtsToRemove - length(hullPoints)
        }
    }



    hullPoints <- chull(pts[1, ], pts[2, ])

    list(allPts=allPts, pts=pts, hullPoints=hullPoints)

}

#' @title Plot covariance contours and component loadings
#'
#' @param Vsub
#' @param Osamps
#' @param OmegaSamps
#' @param s2samps
#' @param ngroups
#' @param hline
#' @param col
#' @param pch
#' @param lty
#' @param ymax
#' @param type
#' @param plotPoints
#' @param polar
#' @param cex.axis
#' @param cex.pts
#'
#' @return
#' @importFrom car ellipse
#' @importFrom gridExtra tableGrob grid.arrange
#' @importFrom grid grid.rect gpar
#' @importFrom ggrepel geom_label_repel
#' @importFrom ggforce geom_ellipse
#'
#' @examples
#'
#' @export
covarianceBiplot <- function(Vsub, Osamps, omegaSamps, s2samps, groups_to_plot=1:nrow(s2samps), nlabeled=20) {

  if(ncol(Vsub) != 2) {
    stop("Please provide 2-dimensional subspace")
  }

  npos <- round(nlabeled/4)
  nneg <- round(nlabeled/4)


  xlimits <- c(-1.1, 1.1)*max(abs(Vsub[, 1:2]))
  ylimits <- c(-1.1, 1.1)*max(abs(Vsub[, 1:2]))

  p <- ggplot(as_tibble(Vsub[, 1:2]), colnames=c("V1", "V2")) +
    geom_point(aes(x=V1, y=V2), size=0.5, col="light grey") +
    theme_bw(base_size=20) + xlim(xlimits) + ylim(ylimits)

  pos_x_indices <- order(Vsub[, 1], decreasing=TRUE)[1:npos]
  neg_x_indices <- order(Vsub[, 1], decreasing=FALSE)[1:nneg]
  pos_y_indices <- setdiff(order(Vsub[, 2], decreasing=TRUE), c(pos_x_indices, neg_x_indices))[1:npos]
  neg_y_indices <- setdiff(order(Vsub[, 2], decreasing=FALSE), c(pos_x_indices, pos_y_indices))[1:nneg]

  lambda_max <- lambda_min <- angle <- c()

  ngroups <- length(groups_to_plot)
  for(k in groups_to_plot) {

    pmPsi <- getPostMeanPsi(Osamps[, , k , ],
                            omegaSamps[, k, ],
                            s2samps[k, ], 100)

    eigK <- eigen(pmPsi)
    lambda <- eigK$values
    evecs <- eigK$vectors

    maxIndex <- which.max(lambda)
    lamRatio <- lambda[maxIndex]/lambda[-maxIndex]

    lambda_max <- c(lambda_max, lambda[maxIndex])
    lambda_min <- c(lambda_min, lambda[-maxIndex])
    angle = c(angle, atan(evecs[2, maxIndex]/evecs[1, maxIndex]))

  }

  group_names <- names(groups_to_plot)
  if(is.null(group_names))
    group_names = factor(1:ngroups)

  ellipse_tibble <- tibble(lambda_max=lambda_max,
                           lambda_min=lambda_min,
                           angle=angle,
                           Group=group_names)

  ellipse_tibble <- ellipse_tibble %>%
    mutate(sd_max = sqrt(lambda_max / max(lambda_max))) %>%
    mutate(sd_min = sqrt(lambda_min / max(lambda_max)))


  ## Plot Ellipses
  ## 95% contour
  p <- p + geom_ellipse(data=ellipse_tibble,
                        aes(x0 = 0, y0 = 0, a = sd_max*2*max(xlimits)*0.25, b = sd_min*2*max(xlimits)*0.25, angle = angle, col=Group), size=1.1)

  ## 50% contour
  p <- p + geom_ellipse(data=ellipse_tibble,
                        aes(x0 = 0, y0 = 0, a = sd_max*0.67*max(xlimits)*0.25, b = sd_min*0.67*max(xlimits)*.25, angle = angle, col=Group), size=1.1)


  ## Add labels
  all_indices <- c(pos_x_indices, neg_x_indices, pos_y_indices, neg_y_indices)
  label_data <- tibble(x=Vsub[all_indices, 1], y=Vsub[all_indices, 2],
                       label=rownames(Vsub)[all_indices], type=rep(c("pos_x", "neg_x", "pos_y", "neg_y"), each=nlabeled/4))

  p <- p + geom_point(data=label_data, aes(x=x, y=y), col="red", size=1.5)

  p <- p + geom_label_repel(
    data  = subset(label_data, type=="pos_x"),
    aes(x=x, y=y, label=label),
    nudge_x = 0.5,
    force=1,
    segment.size  = 0.5,
    segment.color = "grey50",
    direction     = "y",
    hjust         = 0,
    size = 2
  ) +
    geom_label_repel(
      data = subset(label_data, type=="neg_x"),
      aes(x=x, y=y, label=label),
      force=1,
      nudge_x = -0.5,
      segment.size  = 0.5,
      segment.color = "grey50",
      direction     = "y",
      hjust         = 1,
      size= 2
    ) +
    geom_label_repel(
      data = subset(label_data, type=="pos_y"),
      aes(x=x, y=y, label=label),
      force=1,
      nudge_y = 0.05,
      segment.size  = 0.5,
      segment.color = "grey50",
      direction = "both",
      size = 2
    ) +
    geom_label_repel(
      data = subset(label_data, type=="neg_y"),
      aes(x=x, y=y, label=label),
      force=1,
      nudge_y = -0.05,
      segment.size  = 0.5,
      segment.color = "grey50",
      direction     = "both",
      size=2
    )

  p

}

eigenvalueDists <- function(OmegaSamps, nsamps, groups,
                            probRegion=0.95, hline=NULL, col=NULL,
                            pch=NULL, lty=NULL, ymax=30, logRatio=FALSE,
                            plotPoints=TRUE, polar=FALSE, cex.axis=1.5) {

  ngroups <- length(groups)

  if(is.null(col)){
    col <- 1:ngroups
  }
  if(is.null(pch)){
    pch=rep(19, ngroups)
  }
  if(is.null(lty)){
    lty=rep(1, ngroups)
  }
  par(mar=c(5.1, 5.1, 4.1, 2.1))

  plot(0, 0, xlim=c(-pi/2, pi/2),
       ylim=c(0, ymax), cex=0, xlab=expression("angle, acos("~U[1]^T*V[1]~")"), ylab=ylab, xaxt="n", cex.axis=cex.axis, cex.lab=1.5)
  axis(1, at=seq(-pi/2, pi/2, by=pi/4), labels=expression(-pi/2, -pi/4, 0, pi/4, pi/2), cex.axis=cex.axis, cex.lab=1.5)

  LambdaSamps <- OmegaSamps[1, , ] / (1 - OmegaSamps[1, , ])

  tib <- as_tibble(LambdaSamps)
  tib
  tib$type <- types
  tib %>% gather(key=Sample, value=Lambda, -type) %>% mutate(Lambda = as.numeric(Lambda)) %>% ggplot(aes(type, Lambda)) + geom_violin(aes(fill=type)) + ylim(c(0, 300))

}





#' @title Plot biplot and posterior plots side by side
#'
#' @param Vsub
#' @param Osamps
#' @param OmegaSamps
#' @param s2samps
#' @param ngroups
#' @param hline
#' @param col
#' @param pch
#' @param lty
#' @param ymax
#' @param type
#' @param plotPoints
#' @param polar
#' @param cex.axis
#' @param cex.pts
#'
#' @return
#' @importFrom car ellipse
#' @importFrom gridExtra tableGrob grid.arrange
#' @importFrom grid grid.rect gpar
#' @importFrom ggrepel geom_label_repel
#' @importFrom ggforce geom_ellipse
#'
#' @examples
#'
#' @export
create_plots <- function(V, Osamps, omegaSamps, s2samps, group1=1, group2=2, view=c(1,2), to_plot=1:dim(Osamps)[3], nlabeled=20, group_names=NULL, R=ncol(V)) {


  g1 <- group1
  g2 <- group2

  pmPsi1 <- getPostMeanPsi(Osamps[, , g1 , ],
                           omegaSamps[, g1, ],
                           s2samps[g1, ], 100)
  pmPsi2 <- getPostMeanPsi(Osamps[, , g2 , ],
                           omegaSamps[, g2, ],
                           s2samps[g2, ], 100)

  O <- svd(pmPsi1[1:R, 1:R] - pmPsi2[1:R, 1:R])$u[, view]

  ## rotate_basis()

  ngroups <- dim(Osamps)[3]

  Vstar <- V[, 1:R] %*% O
  Osamps_proj <- array(dim = c(2, 2, ngroups, 1000))
  omegaSamps_proj <- array(dim = c(2, ngroups, 1000))

  for(i in 1:100) {
    for(k in 1:ngroups) {

      tmp <- t(O) %*% Osamps[1:R, 1:R, k, i]
      omega <- omegaSamps[1:R, k, i]
      eig <- eigen(tmp %*% diag(omega / (1-omega)) %*% t(tmp))
      Osamps_proj[, , k, i] <- eig$vectors
      lambda <- eig$values
      omegaSamps_proj[, k, i] <- lambda/(lambda+1)
    }
  }

  posterior_plot <- posteriorPlot(Osamps_proj, omegaSamps_proj,
                                  samples$s2samps, nsamps=50,
                                  groups_to_plot=to_plot,
                                  probRegion=0.95)

  biplot <- covarianceBiplot(Vstar, Osamps_proj, omegaSamps_proj, samples$s2samps,
                             groups_to_plot=to_plot, nlabeled=40)

  if(is.null(dev.list()))
    dev.new(width=14, height=7)

  posterior_plot + biplot

}

#' @export
rotate_basis <- function(V, Osamps, omegaSamps, s2samps, group1=1, group2=2) {

  g1 <- 1
  g2 <- 2

  pmPsi1 <- getPostMeanPsi(Osamps[, , g1 , ],
                           omegaSamps[, g1, ],
                           s2samps[g1, ], 100)
  pmPsi2 <- getPostMeanPsi(Osamps[, , g2 , ],
                           omegaSamps[, g2, ],
                           s2samps[g2, ], 100)

  R <- s
  O <- svd(pmPsi1[1:R, 1:R] - pmPsi2[1:R, 1:R])

  Vstar <- V[, 1:R] %*% O

  list(rotV=Vstar, rotMat=O)

}




#' @export
compute_variance_explained <- function(V, Ylist, nvec, s2vec) {

    P <- nrow(V)
    S <- ncol(V)
    ngroups <- length(Ylist)

    ## Use qudratic fomr to compute goodness of fit
    evalRatiosQuad <- sapply(1:S, function(M) {
        sapply(1:ngroups, function(k) {

            YVt <- t(V[, 1:M]) %*% t(Ylist[[k]])
            numer <- tr(YVt %*% t(YVt))/nvec[k]

            evals <- svd(Ylist[[k]])$d[1:min(M, nvec[k])]^2/nvec[k]
            b <- (s2vec[k]*P/nvec[k] - evals - 1)

            quadSol <- ifelse(b^2 - 4*evals < 0, (evals - (s2vec[k]*P/nvec[k])), suppressWarnings((-b + sqrt(b^2 - 4*evals))/2))
            denom <- sum(quadSol)

            (numer/denom)
        })
    })

    evalRatiosQuad
}







steinsLoss <- function(C1, C2inv) {

    sum(diag(C1 %*% C2inv)) - log(det(C1 %*% C2inv)) - nrow(C1)

}




getSigmaInv <- function(P, U, Omega, s2) {

    1/s2*(diag(P) - U %*% diag(Omega) %*% t(U))

}

getPostMeanSigmaInv <- function(P, USamps, OmegaSamps, s2vec, nsamps) {
    SigmaInvList <- lapply(1:nsamps, function(i) {
        1/s2vec[i]*(diag(P) - USamps[, , i] %*% diag(OmegaSamps[, i]) %*% t(USamps[, , i]))
    })
    apply(simplify2array(SigmaInvList), c(1, 2), mean)
}

## posterior mean of (t(V) x Sig x V)^(-1)
getPostMeanSigmaProjInv <- function(S, V, USamps, OmegaSamps, s2vec, nsamps) {
    SigmaProjInvList <- lapply(1:nsamps, function(i) {

        1/s2vec[i]*(diag(S) - (t(V) %*% USamps[, , i]) %*%
                    diag(OmegaSamps[, i]) %*%
                    t(t(V) %*% USamps[, , i]))
    })
    apply(simplify2array(SigmaProjInvList), c(1, 2), mean)
}

#' @export
getPostMeanPsi <- function(Osamps, OmegaSamps, s2vec, nsamps) {
    PsiList <- lapply(1:nsamps, function(i) {
        s2vec[i]*(Osamps[, , i] %*% diag(OmegaSamps[, i]/(1-OmegaSamps[, i])) %*% t(Osamps[, , i]))
    })

    apply(simplify2array(PsiList), c(1, 2), function(x) mean(x, na.rm=TRUE))
}

getPostMeanSigma <- function(P, Ulist, OmegaList, s2vec) {

    SigmaList <- sapply(1:nsamps, function(i) {
        s2vec[i]*(Ulist[, , i] %*% diag((1-omega)/omega) %*% t(Ulist[, , i]) + diag(P))
    })
    apply(simplify2array(SigmaList), c(1, 2), mean)
}

## Get mean predictions for each group relative to truth

#' Title TODO
#'
#' @param Osamps
#' @param OmegaSamps
#' @param s2samps
#' @param nsamps
#' @param groups
#' @param probRegion
#' @param hline
#' @param col
#' @param pch
#' @param lty
#' @param ymax
#' @param type
#' @param plotPoints
#' @param polar
#' @param cex.axis
#' @param cex.pts
#' @param splitGroups
#' @param splitPts
#'
#' @return
#'
#'
#' @export
makePrediction <- function(Y, V, Osamps, omegaSamps, s2samps,
                           ngroups = dim(Osamps)[1],
                           nsamps=dim(Osamps)[4],
                           numToAvg=nsamps/2) {

    YV <- Y %*% V
    SigmaHat <- list()

    ll_array <- array(0, dim=c(nrow(Y), numToAvg, ngroups))
    for(i in (nsamps-numToAvg+1):nsamps) {

        for(k  in 1:ngroups) {
            Ok <- Osamps[, , k, i]
            om <- omegaSamps[, k, i]
            Lam <- diag(om/(1-om))
            s2 <- s2samps[k, i]
            SigmaHat[[k]] <-  s2*(Ok %*% Lam %*% t(Ok) + diag(nrow(Ok)))
        }

        ll <- computeLL(YV, SigmaHat)
        ll_array[, nsamps - i, ] <- ll
    }
    median_ll <- apply(ll_array, c(1, 3), median)
    dim(median_ll)

    priorWeights <- rep(1/length(SigmaHat), length(SigmaHat))
    probs <- apply(median_ll, 1, function(probVec) {
        normProbs <- probVec - max(probVec)
        normProbs <- exp(normProbs) * priorWeights /
            sum(exp(normProbs) * priorWeights)

        normProbs
    })

    t(probs)

}

computeLL <- function(Y, SigmaList, priorWeights=rep(1/length(SigmaList),
                                                     length(SigmaList))) {

    priorWeights <- priorWeights / sum(priorWeights)
    unnormalizedProbs <- sapply(SigmaList, function(Sigma) {
        dmvnorm(as.matrix(Y), sigma=Sigma, log=TRUE)
    })

    unnormalizedProbs
}


computeMembershipProbabilities <-  function(Y, SigmaList,
                                            priorWeights=rep(1/length(SigmaList),
                                                             length(SigmaList))) {

    priorWeights <- priorWeights / sum(priorWeights)
    unnormalizedProbs <- sapply(SigmaList, function(Sigma) {
        dmvnorm(Y, sigma=Sigma, log=TRUE)
    })

    if( !class(unnormalizedProbs) == "matrix")
        unnormalizedProbs <- matrix(unnormalizedProbs, nrow=1, ncol=length(SigmaList))

    probs <- apply(unnormalizedProbs, 1, function(probVec) {
        normProbs <- probVec - max(probVec)
        normProbs <- exp(normProbs) * priorWeights /
            sum(exp(normProbs) * priorWeights)

        normProbs
    })

    t(probs)
}


## Optimal threshold from Gavish, Donoho 2014
#' Title
#'
#' @param Y
#'
#' @return
#' @export
#'
#' @examples
getRank <- function(Y) {

  svals <- svd(Y)$d

  m <- max(nrow(Y), ncol(Y))
  n <- min(nrow(Y), ncol(Y))

  if(m==n) {
    rank <- sum(svals > 2.858*median(svals))
  } else {
    beta <- n/m
    omeg <- 0.56*beta^3 - 0.95*beta^2 + 1.82*beta + 1.43
    rank <- sum(svals > omeg*median(svals))
  }

  rank
}


sampleV2 <- function(Slist, Ulist, s2vec, OmegaList, V, method="gibbs") {

    K <- length(Ulist)
    S <- ncol(V)
    P <- nrow(V)

    ## First save the Bk matrices for each k
    BkList <- OkList <- list()
    for( k in 1:K ){
        OkList[[k]] <- t(V) %*% Ulist[[k]]
         omegaK <- OmegaList[[k]]
         BkList[[k]] <- OkList[[k]] %*% (diag(omegaK, nrow=length(omegaK)) /
                                (2*s2vec[k])) %*% t(OkList[[k]])
    }

    ## Randomly sample a column,  i
    for( i in sample(1:S) ) {

      N <- NullC(V[, -i])

        ## Get the linear term
        ## for j neq i
        C <- rep(0, P)
        for( j in setdiff(1:S, i)) {
            Ak_bij <- matrix(0, P, P)
            for(k in 1:K) {
                Ak <- Slist[[k]]
                bk_ij <- BkList[[k]][i, j]
                Ak_bij <- Ak_bij + Ak*bk_ij
            }
            C <- C+t(V[, j]) %*% Ak_bij
        }
        C <- 2*C

        ## Get the quadratic term
        Ak_bii <- matrix(0, P, P)
        for(k in 1:K) {
            Ak <- Slist[[k]]
            bk_ii <- BkList[[k]][i, i]
            Ak_bii <- Ak_bii+Ak*bk_ii
        }
        A <- Ak_bii

        Ctilde <- as.vector(C %*% N)
        Atilde <- t(N) %*% A %*% N

        NV <- as.vector(t(N) %*% V[, i])

        V[, i] <- N %*% R.rbmf.vector.mises(Atilde, Ctilde, NV)

    }

    ## Propose swaps to handle multi-modality
    if(ncol(V) > 1) {

        if(ncol(V)==2) {
            Vswap <- V[, 2:1]
        } else {
            perm <- sample(2:ncol(V))
            onePos <- sample(2:ncol(V), size=1)
            perm <- c(perm[1:(onePos-1)], onePos, perm[onePos:length(perm)])
            Vswap <- V[, sample(1:ncol(V))]
        }

        logdens <- function(Vcur) {
            sum(
                sapply(1:length(Slist), function(k) {
                    Ok <- OkList[[k]]
                    omegaK <- OmegaList[[k]]
                    B <- Ok %*% (diag(omegaK, nrow=length(omegaK)) /
                                 (2*s2vec[k])) %*% t(Ok)

                    tr(Slist[[k]] %*% Vcur %*%
                       B %*% t(Vcur))
                })
            )
        }

        logMFswap <- logdens(Vswap)
        logMFcur <- logdens(V)
        log.mh <- logMFswap - logMFcur
        u <- runif(1,0,1)
        if ( log(u) < log.mh ) {
            V <- Vswap
            print("SWAP")
        }
    }


    V
}

