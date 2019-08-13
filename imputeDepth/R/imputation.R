################################################################################
## File:             imputation.R                                             ##
## Created by:       Pavlo Mozharovskyi                                       ##
## Last revised:     24.06.2018                                               ##
##                                                                            ##
## Contains functions for depth-based imputation.                             ##
##                                                                            ##
################################################################################

kth <- function(x, k){
  a <- .C("KthElement",
          as.double(x),
          as.integer(length(x)),
          as.integer(k),
          kthel = double(1))
  return (a$kthel)
}

imputeEllP <- function(point, Sigma.inv){
  point.new <- point
  d <- length(point)
  index <- which(is.na(point))
  if (length(index) == 1){
    point.new[index] <- -point.new[-index] %*% Sigma.inv[index,-index] /
      Sigma.inv[index,index]
  }else{
    index <- which(is.na(point))
    A <- Sigma.inv[index,index]
    b <- -Sigma.inv[index,(1:d)[-index], drop = FALSE] %*% point[-index]
    point.new[index] <- solve(A) %*% b
  }
  return (point.new)
}

#' Single imputation by Mahalanobis depth
#'
#' Imputes missing values using Mahalanobis depth function.
#'
#' @param X Data set with missing entries.
#'
#' @param num.iter Number of imputation iteraions, 100 by default.
#'
#' @param alpha Parameter of the Minimum Covariance Determinan estimator
#' used for calculating covariance matrix required for imputation with the
#' Mahalanobis.
#'
#' @param X.start Starting data set; if \code{NULL}, slightly perturbed mean
#' imputation based on existing entries serves as a starting point.
#'
#' @param mu Distributinal mean; if specified together with \code{Sigma}, one
#' imputing iteration is performed only.
#'
#' @param Sigma Distributinal covariance matrix; if specified together with
#' \code{mu}, one imputing iteratin is performed only.
#'
#' @return Imputed data set.
impute.depth.Mahalanobis <- function(X, num.iter = 100, X.start = NULL,
                                     alpha = 1, mu = NULL, Sigma = NULL){
  if (!is.null(mu) && !is.null(Sigma)){
    miss <- which(is.na(X), arr.ind=T)
    X.prep <- X
    X.prep <- t(t(X.prep) - mu)
    Inv.Sigma <- solve(Sigma)
    miss.rowi <- unique(miss[,1])
    X.new <- X.prep
    X.new <- X.new[miss.rowi,,drop=FALSE]
    X.prep[miss.rowi,] <- t(apply(X.new, 1, imputeEllP, Inv.Sigma))
    X.prep <- t(t(X.prep) + mu)
    return(X.prep)
  }
  if (is.null(X.start)){
    X.prep <- X
    miss.label <- is.na(X.prep)
    # TODO: The following line is dangerous when each observation has missings
    X.prep[miss.label] <- matrix(rep(colMeans(X.prep, na.rm = TRUE),
                                     nrow(X.prep)), nrow = nrow(X.prep),
                                 byrow = TRUE)[miss.label] +
      rnorm(n = sum(miss.label), mean = 0, sd = 0.001)
  }else{
    X.prep <- X.start
  }
  miss <- which(is.na(X), arr.ind=T)
  for (i in 1:num.iter){
    if (alpha > 0.99){
      mu.tmp <- colMeans(X.prep)
      Sigma.tmp <- cov(X.prep)
    }else{
      mcd.est <- covMcd(X.prep, alpha = alpha)
      mu.tmp <- mcd.est$center
      Sigma.tmp <- mcd.est$cov
    }
    X.prep <- t(t(X.prep) - mu.tmp)
    Inv.Sigma.tmp <- solve(Sigma.tmp)
    miss.rowi <- unique(miss[,1])
    X.new <- X.prep
    X.new[miss] <- NA
    X.new <- X.new[miss.rowi,,drop=FALSE]
    X.prep[miss.rowi,] <- t(apply(X.new, 1, imputeEllP, Inv.Sigma.tmp))
    X.prep <- t(t(X.prep) + mu.tmp)
  }
  return (X.prep)
}

imp.depth.spatial <- function(X, rob = 1, max.iter = 25, hist = FALSE){
  Mt <- t(X)
  misst.indices <- which(is.na(Mt))
  Mt[misst.indices] <- t(matrix(rep(colMeans(X, na.rm = TRUE), nrow(X)),
                                nrow = nrow(X),
                                byrow = TRUE))[misst.indices]
  cat("Spatial depth imputation:")
  if (hist){
    ress <- list("")
    ress[[1]] <- t(Mt)
  }
  for (i in 1:max.iter){
    MtOld <- Mt
    if (rob >= 0.99){
      w.est <- 0
    }else{
      mcd.est <- covMcd(t(Mt), alpha = rob)
      mcd.eig <- eigen(mcd.est$cov)
      w.est <- as.vector(mcd.eig$vectors %*% diag(1/sqrt(mcd.eig$values)))
    }
    a <- .C("ImputeSpatialOneIter",
            as.double(as.vector(Mt)),
            as.integer(ncol(Mt)),
            as.integer(nrow(Mt)),
            as.double(w.est),
            as.integer(misst.indices - 1),
            as.integer(length(misst.indices)),
            imputed = double(length(misst.indices)))
    Mt[misst.indices] <- a$imputed
    if (hist){
      ress[[i + 1]] <- t(Mt)
    }
    dif <- max(abs(MtOld - Mt))
    cat(" ", i, "(", round(dif, 6), ")", sep = "")
    if (dif < 0.001){break}
  }
  cat(".\n")
  if (hist){
    return (list(X = t(Mt), hist = ress))
  }else{
    return (t(Mt))
  }
}

imp.depth.projection <- function(X, r = 10000){
  Mt <- t(X)
  misst.indices <- which(is.na(Mt))
  Mt[misst.indices] <- t(matrix(rep(colMeans(X, na.rm = TRUE), nrow(X)),
                                nrow = nrow(X),
                                byrow = TRUE))[misst.indices]
  cat("Projection depth (", r, " dirs) imputation:", sep = "")
  for (i in 1:30){
    MtOld <- Mt
    if (i == 1 || i == 6){
      a <- .C("ImputeProjectionOneIter",
              as.double(as.vector(Mt)),
              as.integer(ncol(Mt)),
              as.integer(nrow(Mt)),
              as.integer(misst.indices - 1),
              as.integer(length(misst.indices)),
              as.integer(0),
              as.integer(r),
              dirs = double(r * nrow(Mt)),
              imputed = double(length(misst.indices)))
      dirs <- a$dirs
    }else{
      a <- .C("ImputeProjectionOneIter",
              as.double(as.vector(Mt)),
              as.integer(ncol(Mt)),
              as.integer(nrow(Mt)),
              as.integer(misst.indices - 1),
              as.integer(length(misst.indices)),
              as.integer(1),
              as.integer(r),
              as.double(dirs),
              imputed = double(length(misst.indices)))
    }
    Mt[misst.indices] <- a$imputed
    dif <- max(abs(MtOld - Mt))
    cat(" ", i, "(", round(dif, 6), ")", sep = "")
    if (dif < 0.001){break}
  }
  cat(".\n")
  return (t(Mt))
}

imp.depth.zonoid.o <- function(X, rob = 1, hist = FALSE){
  Mt <- t(X)
  misst.indices <- which(is.na(Mt))
  miss.rows <- rowSums(is.na(X)) > 0.5
  Mt[misst.indices] <- t(matrix(rep(colMeans(X, na.rm = TRUE), nrow(X)),
                                nrow = nrow(X),
                                byrow = TRUE))[misst.indices]
  #cat("Zonoid depth imputation:")
  if (hist){
    ress <- list("")
    ress[[1]] <- t(Mt)
  }
  for (i in 1:25){
    # Save the previous-iteration data
    MtOld <- Mt
    # Separate points on the convex hull
    tMt <- t(Mt)
    rows.qh <- sort(unique(as.vector(convhulln(tMt))))
    miss.rows.qh <- which(miss.rows & (1:nrow(X) %in% rows.qh))
    misst.indices.qh <- which(t(matrix(rep(miss.rows & (1:nrow(X) %in% rows.qh),
                                           ncol(X)),
                                       ncol = ncol(X)) & is.na(X)))
    misst.indices <- which(t(!matrix(rep(miss.rows & (1:nrow(X) %in% rows.qh),
                                         ncol(X)), ncol = ncol(X)) & is.na(X)))
    # Impute the points on the convex hull if such exist
    qh.run <- FALSE
    if (length(misst.indices.qh) > 0 && i < 6){
      qh.run <- TRUE
      if (rob >= 0.99){
        w.est <- 0
      }else{
        mcd.est <- covMcd(tMt, alpha = rob)
        mcd.eig <- eigen(mcd.est$cov)
        w.est <- as.vector(mcd.eig$vectors %*% diag(1/sqrt(mcd.eig$values)))
      }
      a <- .C("ImputeSpatialOneIter",
              as.double(as.vector(Mt)),
              as.integer(ncol(Mt)),
              as.integer(nrow(Mt)),
              as.double(w.est),
              as.integer(misst.indices.qh - 1),
              as.integer(length(misst.indices.qh)),
              imputed = double(length(misst.indices.qh)))
      Mt[misst.indices.qh] <- a$imputed
    }else{
      misst.indices <- sort(c(misst.indices, misst.indices.qh))
    }
    # Input the points inside of the convex hull
    a <- .C("ImputeZonoidOneIter",
            as.double(as.vector(Mt)),
            as.integer(ncol(Mt)),
            as.integer(nrow(Mt)),
            as.integer(misst.indices - 1),
            as.integer(length(misst.indices)),
            imputed = double(length(misst.indices)))
    Mt[misst.indices] <- a$imputed
    if (hist){
      ress[[i + 1]] <- t(Mt)
    }
    # Check whether to continue
    dif <- max(abs(MtOld - Mt))
    cat(" ", i, "(", round(dif, 6), ifelse(qh.run, "qh", ""), ")", sep = "")
    if (dif < 0.00001){break}
  }
  cat(".\n")
  if (hist){
    return (list(X = t(Mt), hist = ress))
  }else{
    return (t(Mt))
  }
}

imp.depth.halfspace.o <- function(X, r = 1000, rob = 0.5){
  Mt <- t(X)
  misst.indices <- which(is.na(Mt))
  miss.rows <- rowSums(is.na(X)) > 0.5
  Mt[misst.indices] <- t(matrix(rep(colMeans(X, na.rm = TRUE), nrow(X)),
                                nrow = nrow(X),
                                byrow = TRUE))[misst.indices] +
    rnorm(n = length(misst.indices), mean = 0, sd = 0.001)
  cat("Tukey depth (", r, " dirs) imputation:", sep = "")
  for (i in 1:25){
    # Save the previous-iteration data
    MtOld <- Mt
    # Separate points on the convex hull
    tMt <- t(Mt)
    rows.qh <- sort(unique(as.vector(convhulln(tMt))))
    miss.rows.qh <- which(miss.rows & (1:nrow(X) %in% rows.qh))
    misst.indices.qh <- which(t(matrix(rep(miss.rows & (1:nrow(X) %in% rows.qh),
                                           ncol(X)),
                                       ncol = ncol(X)) & is.na(X)))
    misst.indices <- which(t(!matrix(rep(miss.rows & (1:nrow(X) %in% rows.qh),
                                         ncol(X)), ncol = ncol(X)) & is.na(X)))
    # Impute the points on the convex hull if such exist
    qh.run <- FALSE
    if (length(misst.indices.qh) > 0 && i < 6){
      qh.run <- TRUE
      if (rob >= 0.99){
        w.est <- 0
      }else{
        mcd.est <- covMcd(tMt, alpha = rob)
        mcd.eig <- eigen(mcd.est$cov)
        w.est <- as.vector(mcd.eig$vectors %*% diag(1/sqrt(mcd.eig$values)))
      }
      a <- .C("ImputeSpatialOneIter",
              as.double(as.vector(Mt)),
              as.integer(ncol(Mt)),
              as.integer(nrow(Mt)),
              as.double(w.est),
              as.integer(misst.indices.qh - 1),
              as.integer(length(misst.indices.qh)),
              imputed = double(length(misst.indices.qh)))
      Mt[misst.indices.qh] <- a$imputed
    }else{
      misst.indices <- sort(c(misst.indices, misst.indices.qh))
    }
    # Input the points inside of the convex hull
    if (i == 1 || i == 6){
      a <- .C("ImputeHalfspaceOneIter",
              as.double(as.vector(Mt)),
              as.integer(ncol(Mt)),
              as.integer(nrow(Mt)),
              as.integer(misst.indices - 1),
              as.integer(length(misst.indices)),
              as.integer(0),
              as.integer(r),
              dirs = double(r * nrow(Mt)),
              imputed = double(length(misst.indices)))
      dirs <- a$dirs
    }else{
      a <- .C("ImputeHalfspaceOneIter",
              as.double(as.vector(Mt)),
              as.integer(ncol(Mt)),
              as.integer(nrow(Mt)),
              as.integer(misst.indices - 1),
              as.integer(length(misst.indices)),
              as.integer(1),
              as.integer(r),
              as.double(dirs),
              imputed = double(length(misst.indices)))
    }
    Mt[misst.indices] <- a$imputed
    # Check whether to continue
    dif <- max(abs(MtOld - Mt))
    cat(" ", i, "(", round(dif, 6), ifelse(qh.run, "qh", ""), ")", sep = "")
    if (dif < 0.001){break}
  }
  cat(".\n")
  return (t(Mt))
}

imp.depth.halfspace.ex.o <- function(X, rob = 0.5){
  Mt <- t(X)
  misst.indices <- which(is.na(Mt))
  miss.rows <- rowSums(is.na(X)) > 0.5
  Mt[misst.indices] <- t(matrix(rep(colMeans(X, na.rm = TRUE), nrow(X)),
                                nrow = nrow(X),
                                byrow = TRUE))[misst.indices]
  cat("Tukey depth imputation:", sep = "")
  for (i in 1:19){
    # Save the previous-iteration data
    MtOld <- Mt
    # Separate points on the convex hull
    tMt <- t(Mt)
    rows.qh <- sort(unique(as.vector(convhulln(tMt))))
    miss.rows.qh <- which(miss.rows & (1:nrow(X) %in% rows.qh))
    misst.indices.qh <- which(t(matrix(rep(miss.rows & (1:nrow(X) %in% rows.qh),
                                           ncol(X)),
                                       ncol = ncol(X)) & is.na(X)))
    misst.indices <- which(t(!matrix(rep(miss.rows & (1:nrow(X) %in% rows.qh),
                                         ncol(X)), ncol = ncol(X)) & is.na(X)))
    # Impute the points on the convex hull if such exist
    qh.run <- FALSE
    if (length(misst.indices.qh) > 0 && i < 6){
      qh.run <- TRUE
      if (rob >= 0.99){
        w.est <- 0
      }else{
        mcd.est <- covMcd(tMt, alpha = rob)
        mcd.eig <- eigen(mcd.est$cov)
        w.est <- as.vector(mcd.eig$vectors %*% diag(1/sqrt(mcd.eig$values)))
      }
      a <- .C("ImputeSpatialOneIter",
              as.double(as.vector(Mt)),
              as.integer(ncol(Mt)),
              as.integer(nrow(Mt)),
              as.double(w.est),
              as.integer(misst.indices.qh - 1),
              as.integer(length(misst.indices.qh)),
              imputed = double(length(misst.indices.qh)))
      Mt[misst.indices.qh] <- a$imputed
    }else{
      misst.indices <- sort(c(misst.indices, misst.indices.qh))
    }
    # Input the points inside of the convex hull
    a <- .C("ImputeHalfspaceExOneIter",
            as.double(as.vector(Mt)),
            as.integer(ncol(Mt)),
            as.integer(nrow(Mt)),
            as.integer(misst.indices - 1),
            as.integer(length(misst.indices)),
            imputed = double(length(misst.indices)))
    Mt[misst.indices] <- a$imputed
    # Check whether to continue
    dif <- max(abs(MtOld - Mt))
    cat(" ", i, "(", round(dif, 6), ifelse(qh.run, "qh", ""), ")", sep = "")
    if (dif < 0.001){break}
  }
  cat(".\n")
  return (t(Mt))
}

#' Single local-depth imputation
#'
#' Imputes missing values using local depth with a pluged-in specified depth
#' function.
#'
#' @param X Data set with missing entries.
#'
#' @param par.loc Localisation parameter for the local depth.
#'
#' @param depth A notion of depth to be pluged-in in the local depth chosen
#' among "Mahalanobis", "zonoid" (default), "Tukey".
#'
#' @param parMcd.impute Parameter of the Minimum Covariance Determinan estimator
#' used for calculating covariance matrix required for imputation with the
#' Mahalanobis or spatial depth.
#'
#' @param max.iter Maximum number of imputation iteraions, 25 by default.
#'
#' @param eps.precision Numerical precision, i.e. maximum deviation from
#' previous iteration to continue, \code{0.001} by default.
#'
#' @param core.count Number of cores to be used, if \code{1} (default) then no
#' parallelisation involved.
#'
#' @param cl.type Type of the cluster, \code{"SOCK"} by default.
#'
#' @param verbosity Verbosity parameter: \code{0} for no messages at
#' all, \code{1} (default) for minimal messages, \code{2} for a diagnostic mode
#' detailing each imputation iteration.
#'
#' @param hist Logical indicating whether the history of each iteration should
#' be returned, \code{FALSE} by default.
#'
#' @return Imputed data set.
impute.depth.local <- function(X, par.loc = 0.5, depth = "Mahalanobis",
                               parMcd.impute = 0.5, max.iter = 25,
                               eps.precision = 0.001, core.count = 1,
                               cl.type = "SOCK", verbosity = 1, hist = FALSE){
  Mt <- t(X)
  misst.indices <- which(is.na(Mt))
  miss.rows <- rowSums(is.na(X)) > 0.5
  Mt[misst.indices] <- t(matrix(rep(colMeans(X, na.rm = TRUE), nrow(X)),
                                nrow = nrow(X),
                                byrow = TRUE))[misst.indices]

  Mt <- t(imputeKnn(X))

  cat("Local (", depth, ") depth imputation:", sep = "")
  if (hist){
    ress <- list("")
    ress[[1]] <- t(Mt)
  }
  for (i in 1:max.iter){
    # Save the previous-iteration data
    MtOld <- Mt
    # Separate points on the convex hull
    tMt <- t(Mt)
    rows.qh <- sort(unique(as.vector(convhulln(tMt))))
    miss.rows.qh <- which(miss.rows & (1:nrow(X) %in% rows.qh))
    misst.indices.qh <- which(t(matrix(rep(miss.rows & (1:nrow(X) %in% rows.qh),
                                           ncol(X)),
                                       ncol = ncol(X)) & is.na(X)))
    misst.indices <- which(t(!matrix(rep(miss.rows & (1:nrow(X) %in% rows.qh),
                                         ncol(X)), ncol = ncol(X)) & is.na(X)))
    # Impute the points on the convex hull if such exist
    qh.run <- FALSE
    if (length(misst.indices.qh) > 0 && i < -1){
      qh.run <- TRUE
      if (parMcd.impute >= 0.99){
        w.est <- 0
      }else{
        mcd.est <- covMcd(tMt, alpha = parMcd.impute)
        mcd.eig <- eigen(mcd.est$cov)
        w.est <- as.vector(mcd.eig$vectors %*% diag(1/sqrt(mcd.eig$values)))
      }
      a <- .C("ImputeSpatialOneIter",
              as.double(as.vector(Mt)),
              as.integer(ncol(Mt)),
              as.integer(nrow(Mt)),
              as.double(w.est),
              as.integer(misst.indices.qh - 1),
              as.integer(length(misst.indices.qh)),
              imputed = double(length(misst.indices.qh)))
      Mt[misst.indices.qh] <- a$imputed
    }else{
      misst.indices <- sort(c(misst.indices, misst.indices.qh))
    }
    # Input the points inside of the convex hull
    depth.int = 0
    if (toupper(depth) == "HALFSPACE" ||
        toupper(depth) == "TUKEY" ||
        toupper(depth) == "LOCATION"){
      depth.int = 1
    }
    if (toupper(depth) == "ZONOID"){
      depth.int = 2
    }
    if (toupper(depth) == "MAHALANOBIS"){
      depth.int = 3
    }
    a <- .C("ImputeLocalOneIter",
            as.double(as.vector(Mt)),
            as.integer(ncol(Mt)),
            as.integer(nrow(Mt)),
            as.double(par.loc),
            as.integer(depth.int),
            as.integer(misst.indices - 1),
            as.integer(length(misst.indices)),
            imputed = double(length(misst.indices)))
    Mt[misst.indices] <- a$imputed
    # Check whether to continue
    if (hist){
      ress[[i + 1]] <- t(Mt)
    }
    dif <- max(abs(MtOld - Mt))
    cat(" ", i, "(", round(dif, 6), ifelse(qh.run, "qh", ""), ")", sep = "")
    if (dif < eps.precision){break}
  }
  cat(".\n")
  if (hist){
    return (list(X = t(Mt), hist = ress))
  }else{
    return (t(Mt))
  }
}

depth.local <- function(objects, data, par.loc = 0.5, depth = "Mahalanobis"){
  objects <- as.matrix(objects, ncol = ncol(data))
  #print(objects)
  #print(typeof(objects))
  depth.int = 0
  if (toupper(depth) == "HALFSPACE" ||
      toupper(depth) == "TUKEY" ||
      toupper(depth) == "LOCATION"){
    depth.int = 1
  }
  if (toupper(depth) == "ZONOID"){
    depth.int = 2
  }
  if (toupper(depth) == "MAHALANOBIS"){
    depth.int = 3
  }
  depth <- .C("LocalDepth", as.double(t(objects)),
              as.double(t(data)),
              as.integer(nrow(data)),
              as.integer(ncol(data)),
              as.integer(nrow(objects)),
              as.double(par.loc),
              as.integer(depth.int),
              res = double(nrow(objects)))
  return (depth$res)
}

tune.impute.depth.local <- function(X, pars.loc = 1:10/10, l = 10, p = 0.03,
                                    depth = "Mahalanobis", rob = 1){
  # Initial imputation
  miss.label <- is.na(X)
  miss.n <- sum(miss.label)
  X.prep <- X
  # Cross-validation
  errors <- rep(0, length(pars.loc))
  for (i in 1:l){
    cat("CV for", pars.loc[j], ":\n")
    X.cvi <- as.matrix(X.prep)
    # Check on p portion of newly introduced missing values
    miss.i.added <- sample(which(!is.na(X)), ceiling(length(X) * p))
    X.cvi[miss.i.added] <- NA
    for (j in 1:length(pars.loc)){
      X.impi <- depth.impute.local(X.cvi, par.loc = pars.loc[j],
                                   depth = depth, rob = rob)
      errors[j] <- errors[j] + sum((X.impi[miss.i.added] -
                                      X.prep[miss.i.added])^2)
    }
  }
  print(errors)
  par.best.j <- which.min(errors)
  #  return (kNN(data.frame(X), k = k.best)[,1:d])
  return (pars.loc[par.best.j])
}

tune.impute.depth.local_old <- function(X, pars.loc = 1:10/10, l = 10,
                                    depth = "Mahalanobis", rob = 1){
  # Initial imputation
  miss.label <- is.na(X)
  miss.n <- sum(miss.label)
  X.prep <- X
  X.prep[miss.label] <- matrix(rep(colMeans(X.prep, na.rm = TRUE),
                                   nrow(X.prep)), nrow = nrow(X.prep),
                               byrow = TRUE)[miss.label]
  # Cross-validation
  errors <- rep(0, length(pars.loc))
  for (i in 1:l){
    X.cvi <- as.matrix(X.prep)
    X.cvi[sample(1:(nrow(X) * ncol(X)), miss.n)] <- NA
    for (j in 1:length(pars.loc)){
      X.impi <- depth.impute.local(X.cvi, par.loc = pars.loc[j],
                                   depth = depth, rob = rob)
#      X.impi <- kNN(data.frame(X.cvi), k = ks[j])[,1:d]
      errors[j] <- errors[j] + sum((X.impi - X.prep)^2)
    }
  }
  print(errors)
  par.best.j <- which.min(errors)
#  return (kNN(data.frame(X), k = k.best)[,1:d])
  return (pars.loc[par.best.j])
}

impute.depth.Tukey.worker <- function(misst.indices, ...){
  Mt <- list(...)[[1]]$Mt
  res <- .C("ImputeHalfspaceExOneIter",
            as.double(as.vector(Mt)),
            as.integer(ncol(Mt)),
            as.integer(nrow(Mt)),
            as.integer(misst.indices - 1),
            as.integer(length(misst.indices)),
            imputed = double(length(misst.indices)))
  return (res$imputed)
}

impute.depth.rTukey.worker <- function(misst.indices, ...){
  args <- list(...)[[1]]
  res <- .C("ImputeHalfspaceOneIter",
            as.double(as.vector(args$Mt)),
            as.integer(ncol(args$Mt)),
            as.integer(nrow(args$Mt)),
            as.integer(misst.indices - 1),
            as.integer(length(misst.indices)),
            as.integer(0),
            as.integer(args$nProj),
            dirs = double(args$nProj * nrow(args$Mt)),
            imputed = double(length(misst.indices)))
  return (res$imputed)
}

impute.depth.extrTukey.worker <- function(misst.indices, ...){
  args <- list(...)[[1]]
  res <- .C("ImputeHalfExtrOneIter",
            as.double(as.vector(args$Mt)),
            as.integer(ncol(args$Mt)),
            as.integer(nrow(args$Mt)),
            as.integer(misst.indices - 1),
            as.integer(length(misst.indices)),
            as.integer(0),
            as.integer(args$nProj),
            as.integer(args$exactInner),
            as.integer(args$kPar),
            dirs = double(args$nProj * nrow(args$Mt)),
            imputed = double(length(misst.indices)))
  return (res$imputed)
}

print.fun <-  function(outputs, B, args){
  pb <- args$mypb
  len.one.run <- args$len.one.run
  outputs <- unlist(outputs)
  outputs <- outputs[!is.null(outputs)]
  setTxtProgressBar(pb, length(outputs)/len.one.run)
}

#' Single imputation by data depth
#'
#' Imputes missing values using specified depth function.
#'
#' @param X Data set with missing entries.
#'
#' @param depth A notion of depth to be used for imputation chosen among
#' "Mahalanobis", "zonoid" (default), "Tukey", "spatial", "projection",
#' "randomTukey", "simlicial".
#'
#' @param parMcd.impute Parameter of the Minimum Covariance Determinan estimator
#' used for calculating covariance matrix required for imputation with the
#' Mahalanobis or spatial depth.
#'
#' @param depth.outsiders A notion of depth to be used for imputation of
#' outsiders, i.e. (here) points lying on the convex hull of the data set,
#' chosen among "extremity" (default), "Mahalanobis", "spatial", "projection".
#'
#' @param parMcd.outsiders Parameter of the Minimum Covariance Determinan
#' estimator used for calculating covariance matrix required for imputation
#' of the outsiders with the Mahalanobis or spatial depth.
#'
#' @param max.iter Maximum number of imputation iteraions, 50 by default.
#'
#' @param max.iter.outsiders Maximum number of outsider-imputing iterations,
#' 5 by default.
#'
#' @param n.proj Number of projections used to calculate (approximate)
#' projection and random Tukey depths if chosen for \code{depth} and (or)
#' \code{depth.outsiders}, \code{1000} by default.
#'
#' @param p.extreme Portion of points in the tail to consider for the depth
#' extension by the extreme value theory, \code{0.1} by default.
#'
#' @param eps.precision Numerical precision, i.e. maximum deviation from
#' previous iteration to continue, \code{0.001} by default.
#'
#' @param core.count Number of cores to be used, if \code{1} (default) then no
#' parallelisation involved.
#'
#' @param cl.type Type of the cluster, \code{"SOCK"} by default.
#'
#' @param verbosity Verbosity parameter: \code{0} for no messages at
#' all, \code{1} (default) for minimal messages, \code{2} for a diagnostic mode
#' detailing each imputation iteration.
#'
#' @param hist Logical indicating whether the history of each iteration should
#' be returned, \code{FALSE} by default.
#'
#' @return Imputed data set.
impute.depth <- function(X, depth = "zonoid", parMcd.impute = 1,
                         depth.outsiders = "spatial", parMcd.outsiders = 1,
                         max.iter = 50, max.iter.outsiders = 5,
                         n.proj = 1000, p.extreme = 0.1, eps.precision = 0.001,
                         core.count = 1, cl.type = c("SOCK", "MPI"),
                         verbosity = 1, hist = FALSE){
  # Detect missing values
  Mt <- t(X)
  misst.indices <- which(is.na(Mt))
  miss.rows <- rowSums(is.na(X)) > 0.5
  if (verbosity >= 2){
    cat("Imputing", sum(miss.rows), "points having in total",
        length(misst.indices), "missing entries.\n")
  }
  # Initial imputation
  if (verbosity >= 2){
    cat("Initial imputation by mean...")
  }
  Mt[misst.indices] <- t(matrix(rep(colMeans(X, na.rm = TRUE), nrow(X)),
                                nrow = nrow(X),
                                byrow = TRUE))[misst.indices]
  if (verbosity >= 2){
    cat("done.\n")
  }
  if (verbosity >= 1){
    cat("Single imputation using", depth, "depth")
    if (toupper(depth) != "MAHALANOBIS" &&
        toupper(depth) != "PROJECTION" &&
        toupper(depth) != "SPATIAL"){
      cat(" (and", depth.outsiders, "depth for outsiders)")
    }
    cat(":")
  }
  if (verbosity >= 2){
    cat("\n")
  }
  for (i in 1:max.iter){
    if (verbosity >= 2){
      cat("Starting iteration ", i, ":\n", sep = "")
    }
    # Save the previous-iteration data
    MtOld <- Mt
    qh.run <- FALSE
    misst.indices.qh <- NULL
    # Impute outrisers if required by the depth notion
    if (i <= max.iter.outsiders &&
        (toupper(depth) != "MAHALANOBIS" &&
        toupper(depth) != "PROJECTION" &&
        toupper(depth) != "SPATIAL" &&
        toupper(depth.outsiders) != "EXTREMEVAL") &&
        toupper(depth.outsiders) != "NO"){
      # Calculate auxilary statistics where and if needed for imputing outsiders
      if (toupper(depth.outsiders) == "MAHALANOBIS" ||
          toupper(depth.outsiders) == "SPATIAL"){
        if (verbosity >= 2){
          cat("Calculating covariance estimates for outsiders...")
        }
        if (parMcd.outsiders < 0.99){
          if (verbosity >= 2){
            cat("by MCD...")
          }
          est.cov <- NULL
          try(est.cov <- covMcd(t(Mt), alpha = parMcd.outsiders)$cov)
          if (is.null(est.cov)){
            est.cov <- cov(t(Mt))
            warning("MCD covariance estimate for outsiders not computed due to a precision problem, replaced by the moment estimate.")
          }
          if (toupper(depth.outsiders) == "SPATIAL"){
            mcd.eig <- eigen(est.cov)
            w.est <- as.vector(mcd.eig$vectors %*% diag(1/sqrt(mcd.eig$values)))
          }
        }else{
          if (verbosity >= 2){
            cat("by moment...")
          }
          if (toupper(depth.outsiders) == "MAHALANOBIS"){
            est.cov <- cov(t(Mt))
          }
          if (toupper(depth.outsiders) == "SPATIAL"){
            w.est <- 0
          }
        }
        if (verbosity >= 2){
          cat("done.\n")
        }
      }
      # Separate points on the convex hull
      rows.qh <- sort(unique(as.vector(convhulln(t(Mt)))))
      miss.rowi.qh <- which(miss.rows & (1:nrow(X) %in% rows.qh))
      misst.indices.qh <- which(t(matrix(
        rep(miss.rows & (1:nrow(X) %in% rows.qh), ncol(X)),
        ncol = ncol(X)) & is.na(X)))
      misst.indices <- which(t(!matrix(
        rep(miss.rows & (1:nrow(X) %in% rows.qh), ncol(X)),
        ncol = ncol(X)) & is.na(X)))
      # Impute the points on the convex hull if such exist
      qh.run <- FALSE
      if (length(misst.indices.qh) > 0 && i <= max.iter.outsiders){
        qh.run <- TRUE
        # Call the outsider imputation procedure
        if (toupper(depth.outsiders) == "MAHALANOBIS"){
          mu.tmp <- rowMeans(Mt)
          X.prep <- Mt - mu.tmp
          Inv.Sigma.tmp <- solve(est.cov)
          #Inv.Sigma.tmp <- NULL
          #try(Inv.Sigma.tmp <- solve(est.cov), silent = TRUE)
          #if (is.null(Inv.Sigma.tmp)){
          #  warning("Error in solve.default(cov, ...) :
          #          system is computationally singular.")
          #}else{
            miss.rowi <- miss.rowi.qh
            X.new <- X.prep
            X.new[misst.indices.qh] <- NA
            X.new <- t(X.new)
            X.new <- X.new[miss.rowi,,drop=FALSE]
            X.prep <- t(X.prep)
            X.prep[miss.rowi,] <- t(apply(X.new, 1, imputeEllP, Inv.Sigma.tmp))
            X.prep <- t(t(X.prep) + mu.tmp)
            Mt[misst.indices.qh] <- t(X.prep)[misst.indices.qh]
          #}
        }
        if (toupper(depth.outsiders) == "SPATIAL"){
          res <- .C("ImputeSpatialOneIter",
                    as.double(as.vector(Mt)),
                    as.integer(ncol(Mt)),
                    as.integer(nrow(Mt)),
                    as.double(w.est),
                    as.integer(misst.indices.qh - 1),
                    as.integer(length(misst.indices.qh)),
                    imputed = double(length(misst.indices.qh)))
          Mt[misst.indices.qh] <- res$imputed
        }
        if (toupper(depth.outsiders) == "PROJECTION"){
          res <- .C("ImputeProjectionOneIter",
                    as.double(as.vector(Mt)),
                    as.integer(ncol(Mt)),
                    as.integer(nrow(Mt)),
                    as.integer(misst.indices.qh - 1),
                    as.integer(length(misst.indices.qh)),
                    as.integer(0),
                    as.integer(n.proj),
                    dirs = double(n.proj * nrow(Mt)),
                    imputed = double(length(misst.indices.qh)))
          Mt[misst.indices.qh] <- res$imputed
        }
      }
    }else{
      rows.qh <- NULL
      misst.indices <- sort(c(misst.indices, misst.indices.qh))
    }
    # Calculate auxilary statistics where and if needed
    if (toupper(depth) == "MAHALANOBIS" ||
        toupper(depth) == "SPATIAL"){
      if (verbosity >= 2){
        cat("Calculating covariance estimates...")
      }
      if (parMcd.impute < 0.99){
        if (verbosity >= 2){
          cat("by MCD...")
        }
        est.cov <- NULL
        try(est.cov <- covMcd(t(Mt), alpha = parMcd.outsiders)$cov)
        if (is.null(est.cov)){
          est.cov <- cov(t(Mt))
          warning("MCD covariance estimate for imputation not computed due to a precision problem, replaced by the moment estimate.")
        }
        if (toupper(depth) == "SPATIAL"){
          mcd.eig <- eigen(est.cov)
          w.est <- as.vector(mcd.eig$vectors %*% diag(1/sqrt(mcd.eig$values)))
        }
      }else{
        if (verbosity >= 2){
          cat("by moment...")
        }
        if (toupper(depth) == "MAHALANOBIS"){
          est.cov <- cov(t(Mt))
        }
        if (toupper(depth) == "SPATIAL"){
          w.est <- 0
        }
      }
      if (verbosity >= 2){
        cat("done.\n")
      }
    }
    # Input the points inside of the convex hull
    if (toupper(depth) == "MAHALANOBIS"){
      mu.tmp <- rowMeans(Mt)
      X.prep <- Mt - mu.tmp
      Inv.Sigma.tmp <- solve(est.cov)
      miss.rowi <- which(miss.rows & !(1:nrow(X) %in% rows.qh))
      X.new <- X.prep
      X.new[misst.indices] <- NA
      X.new <- t(X.new)
      X.new <- X.new[miss.rowi,,drop=FALSE]
      X.prep <- t(X.prep)
      X.prep[miss.rowi,] <- t(apply(X.new, 1, imputeEllP, Inv.Sigma.tmp))
      X.prep <- t(t(X.prep) + mu.tmp)
      Mt[misst.indices] <- t(X.prep)[misst.indices]
    }
    if (toupper(depth) == "SPATIAL"){
      res <- .C("ImputeSpatialOneIter",
                as.double(as.vector(Mt)),
                as.integer(ncol(Mt)),
                as.integer(nrow(Mt)),
                as.double(w.est),
                as.integer(misst.indices - 1),
                as.integer(length(misst.indices)),
                imputed = double(length(misst.indices)))
      Mt[misst.indices] <- res$imputed
    }
    if (toupper(depth) == "PROJECTION"){
      res <- .C("ImputeProjectionOneIter",
                as.double(as.vector(Mt)),
                as.integer(ncol(Mt)),
                as.integer(nrow(Mt)),
                as.integer(misst.indices - 1),
                as.integer(length(misst.indices)),
                as.integer(0),
                as.integer(n.proj),
                dirs = double(n.proj * nrow(Mt)),
                imputed = double(length(misst.indices)))
      Mt[misst.indices] <- res$imputed
    }
    if (toupper(depth) == "ZONOID"){
      res <- .C("ImputeZonoidOneIter",
                as.double(as.vector(Mt)),
                as.integer(ncol(Mt)),
                as.integer(nrow(Mt)),
                as.integer(misst.indices - 1),
                as.integer(length(misst.indices)),
                imputed = double(length(misst.indices)))
      Mt[misst.indices] <- res$imputed
    }
    if (toupper(depth) == "SIMPLICIAL"){
      res <- .C("ImputeSimplicialExOneIter",
                as.double(as.vector(Mt)),
                as.integer(ncol(Mt)),
                as.integer(nrow(Mt)),
                as.integer(misst.indices - 1),
                as.integer(length(misst.indices)),
                imputed = double(length(misst.indices)))
      Mt[misst.indices] <- res$imputed
    }
    if (toupper(depth) == "TUKEY" ||
        toupper(depth) == "HALFSPACE" ||
        toupper(depth) == "LOCATION"){
      if (toupper(depth.outsiders) == "EXTREMEVAL"){
        kPar <- ceiling(p.extreme * ncol(Mt))
        if (core.count > 1){
          miss.rowi <- which(miss.rows & !(1:nrow(X) %in% rows.qh))
          inputs <- list()
          for (j in 1:length(miss.rowi)){
            # Choose missing coordinates of the current point
            inputs[[j]] <- misst.indices[floor(
              (misst.indices - 1) / miss.rowi[j]) + 1 == miss.rowi[j]]
          }
          # Parallel call
          res <- performParallel(count = core.count, x = inputs,
                                 fun = impute.depth.extrTukey.worker,
                                 printfun = print.fun,
                                 printargs = list(
                                   mypb = txtProgressBar(min = 0,
                                                         max =
                                                           length(miss.rowi),
                                                         style = 3),
                                   len.one.run = 1), printrepl = 1,
                                 cltype = cl.type, ... = list(Mt = Mt,
                                                              nProj = n.proj,
                                                              kPar = kPar,
                                                              exactInner = 1))
          # Place the imputed values into the data set
          for (j in 1:length(miss.rowi)){
            Mt[input[[j]]] <- res[[j]]
          }
        }else{
          res <- .C("ImputeHalfExtrOneIter",
                    as.double(as.vector(Mt)),
                    as.integer(ncol(Mt)),
                    as.integer(nrow(Mt)),
                    as.integer(misst.indices - 1),
                    as.integer(length(misst.indices)),
                    as.integer(0),
                    as.integer(n.proj),
                    as.integer(1),
                    as.integer(kPar),
                    dirs = double(n.proj * nrow(Mt)),
                    imputed = double(length(misst.indices)))
          Mt[misst.indices] <- res$imputed
        }
      }else{
        if (core.count > 1){
          miss.rowi <- which(miss.rows & !(1:nrow(X) %in% rows.qh))
          inputs <- list()
          for (j in 1:length(miss.rowi)){
            # Choose missing coordinates of the current point
            inputs[[j]] <- misst.indices[floor(
              (misst.indices - 1) / miss.rowi[j]) + 1 == miss.rowi[j]]
          }
          # Parallel call
          res <- performParallel(count = core.count, x = inputs,
                                 fun = impute.depth.Tukey.worker,
                                 printfun = print.fun,
                                 printargs = list(
                                   mypb = txtProgressBar(min = 0,
                                                         max =
                                                           length(miss.rowi),
                                                         style = 3),
                                   len.one.run = 1), printrepl = 1,
                                 cltype = cl.type, ... = list(Mt = Mt))
          # Place the imputed values into the data set
          for (j in 1:length(miss.rowi)){
            Mt[input[[j]]] <- res[[j]]
          }
        }else{
          res <- .C("ImputeHalfspaceExOneIter",
                    as.double(as.vector(Mt)),
                    as.integer(ncol(Mt)),
                    as.integer(nrow(Mt)),
                    as.integer(misst.indices - 1),
                    as.integer(length(misst.indices)),
                    imputed = double(length(misst.indices)))
          Mt[misst.indices] <- res$imputed
        }
      }
    }
    if (toupper(depth) == "RANDOMTUKEY"){
      if (toupper(depth.outsiders) == "EXTREMEVAL"){
        kPar <- ceiling(p.extreme * ncol(Mt))
        if (core.count > 1){
          miss.rowi <- which(miss.rows & !(1:nrow(X) %in% rows.qh))
          inputs <- list()
          for (j in 1:length(miss.rowi)){
            # Choose missing coordinates of the current point
            inputs[[j]] <- misst.indices[floor(
              (misst.indices - 1) / miss.rowi[j]) + 1 == miss.rowi[j]]
          }
          # Parallel call
          res <- performParallel(count = core.count, x = inputs,
                                 fun = impute.depth.extrTukey.worker,
                                 printfun = print.fun,
                                 printargs = list(
                                   mypb = txtProgressBar(min = 0,
                                                         max =
                                                           length(miss.rowi),
                                                         style = 3),
                                   len.one.run = 1), printrepl = 1,
                                 cltype = cl.type, ... = list(Mt = Mt,
                                                              nProj = n.proj,
                                                              kPar = kPar,
                                                              exactInner = 0))
          # Place the imputed values into the data set
          for (j in 1:length(miss.rowi)){
            Mt[input[[j]]] <- res[[j]]
          }
        }else{
          res <- .C("ImputeHalfExtrOneIter",
                    as.double(as.vector(Mt)),
                    as.integer(ncol(Mt)),
                    as.integer(nrow(Mt)),
                    as.integer(misst.indices - 1),
                    as.integer(length(misst.indices)),
                    as.integer(0),
                    as.integer(n.proj),
                    as.integer(0),
                    as.integer(kPar),
                    dirs = double(n.proj * nrow(Mt)),
                    imputed = double(length(misst.indices)))
          Mt[misst.indices] <- res$imputed
        }
      }else{
        if (core.count > 1){
          miss.rowi <- which(miss.rows & !(1:nrow(X) %in% rows.qh))
          inputs <- list()
          for (j in 1:length(miss.rowi)){
            # Choose missing coordinates of the current point
            inputs[[j]] <- misst.indices[floor(
              (misst.indices - 1) / miss.rowi[j]) + 1 == miss.rowi[j]]
          }
          # Parallel call
          res <- performParallel(count = core.count, x = inputs,
                                 fun = impute.depth.rTukey.worker,
                                 printfun = print.fun,
                                 printargs = list(
                                   mypb = txtProgressBar(min = 0,
                                                         max =
                                                           length(miss.rowi),
                                                         style = 3),
                                   len.one.run = 1), printrepl = 1,
                                 cltype = cl.type, ... = list(Mt = Mt,
                                                              nProj = n.proj))
          # Place the imputed values into the data set
          for (j in 1:length(miss.rowi)){
            Mt[input[[j]]] <- res[[j]]
          }
        }else{
          res <- .C("ImputeHalfspaceOneIter",
                    as.double(as.vector(Mt)),
                    as.integer(ncol(Mt)),
                    as.integer(nrow(Mt)),
                    as.integer(misst.indices - 1),
                    as.integer(length(misst.indices)),
                    as.integer(0),
                    as.integer(n.proj),
                    dirs = double(n.proj * nrow(Mt)),
                    imputed = double(length(misst.indices)))
          Mt[misst.indices] <- res$imputed
        }
      }
    }
    # Check whether to continue
    dif <- max(abs(MtOld - Mt))
    if (verbosity == 1){
      cat(" ", i, "(", round(dif, 6), ifelse(qh.run, "qh", ""), ")", sep = "")
    }
    if (verbosity == 2){
      cat("Qhull ", ifelse(qh.run, "used", "not used"),
          ", convergence precision: ", round(dif, 6), ".\n", sep = "")
    }
    if (dif < eps.precision){break}
  }
  if (verbosity == 1){
    cat(".\n")
  }
  if (i == max.iter){
    warning(paste("Algorithm did not converge after max.iter = ", i,
                  " iterations.\n", sep = ""))
  }
  return (t(Mt))
}

# Function for improper imputation of a single point
imputeEllPRnd <- function(point, Sigma, depths){
  # Preliminary computations
  point.new <- point
  d <- length(point)
  n <- length(depths)
  index <- which(is.na(point))
  Sigma.inv <- solve(Sigma)
  # Conditional imputation
  if (length(index) == 1){
    point.new[index] <- -point.new[-index] %*% Sigma.inv[index,-index] /
      Sigma.inv[index,index]
  }else{
    index <- which(is.na(point))
    A <- Sigma.inv[index,index]
    b <- -Sigma.inv[index,(1:d)[-index], drop = FALSE] %*% point[-index]
    point.new[index] <- solve(A) %*% b
  }
  # Adding noise: draw from the conditional distribution
  depth <- as.numeric(1/(1 + point.new %*% Sigma.inv %*% point.new))
  if (depth <= depths[1]){
    rnd.depth <- depths[1]
  }else{
    if(depth >= depths[n]){
      rnd.depth <- depths[n]
    }else{
      n.new <- min(which(depth < unique(depths))) - 1
      if (n.new <= 1){
        rnd.depth <- unique(depths)[1]
      }else{
        # Create depth's conditional c.d.f.
        probs <- table(depths)[1:n.new]/sum(table(depths)[1:n.new])
        depths.new <- unique(depths)[1:n.new]
        probs <- probs / sqrt(1/depths.new - 1)^(d - 1) *
          sqrt((1/depths.new - 1) - (1/depth - 1))^(length(index) - 1) *
          (sqrt(1/depths.new - 1) / sqrt((1/depths.new - 1) - (1/depth - 1)))
        probs[n.new] <- probs[n.new - 1]
        probs <- cumsum(probs)
        # Draw the depth of the point
        rnd.prob <- runif(1, probs[1], probs[n.new])
        upperIndex <- min(which(probs > rnd.prob))
        rnd.depth <- (rnd.prob - probs[upperIndex - 1]) /
          (probs[upperIndex] - probs[upperIndex - 1]) *
          (depths[upperIndex] - depths[upperIndex - 1]) + depths[upperIndex - 1]
      }
    }
  }
  # Project the point into the sample's space
  d.na <- sum(is.na(point))
  dir <- rnorm(d.na)
  dir <- dir / sqrt(sum(dir^2))
  if (length(index) > 1){
    Sigma.cond <- Sigma[index,index,drop = FALSE] -
      Sigma[index,-index,drop = FALSE] %*%
      solve(Sigma[-index,-index,drop = FALSE]) %*%
      Sigma[-index,index,drop = FALSE]
    e <- eigen(Sigma.cond)
    L <- e$vectors %*% diag(sqrt(e$values)) %*% t(e$vectors)
    dir <- as.numeric(L %*% dir)
  }
  # Solve quadratic equation
  d.tmp <- 1/rnd.depth - 1
  r <- rep(0, d)
  r[is.na(point)] <- dir
  a.tmp <- r %*% Sigma.inv %*% r
  b.tmp <- 2 * point.new %*% Sigma.inv %*% r
  c.tmp <- point.new %*% Sigma.inv %*% point.new - d.tmp
  Dis <- b.tmp^2 - 4 * a.tmp * c.tmp
  if (Dis <= 0){
    alpha <- 0
  }else{
    x1 <- (-b.tmp + sqrt(Dis)) / (2 * a.tmp)
    x2 <- (-b.tmp - sqrt(Dis)) / (2 * a.tmp)
    if (x1 < 0 && x2 < 0){cat("Negative discriminant error!\n")}
    if (x1 > x2){alpha <- x1}else{alpha <- x2}
  }
  point.cur <- point.new + r * as.vector(alpha)
  # Optimize instead of solving quadratic equation
  #   min.cur <- 0
  #   max.cur <- 1000
  #   point.cur <- point.new
  #   while(abs(max.cur - min.cur) > sqrt(.Machine$double.eps)){
  #     mid.cur <- (min.cur + max.cur) / 2
  #     point.cur[is.na(point)] <- (point.new[is.na(point)] + dir * mid.cur)
  #     depth.cur <- as.numeric(1/(1 + point.cur %*% Sigma.inv %*% point.cur))
  #     if (depth.cur < rnd.depth){
  #       max.cur <- mid.cur
  #     }else{
  #       min.cur <- mid.cur
  #     }
  #   }
  return (point.cur)
}

#' Improper imputation for elliptical family
#'
#' Depth-based improper imputation for an elliptically symmetric distribution.
#'
#' @param m Number of improper-imputed data sets.
#'
#' @param iter.burnin Number of burn-in iterations.
#'
#' @param iter.skip Sampling period.
#'
#' @param verbosity Verbosity parameter: \code{0} for no messages at
#' all, \code{1} (default) for minimal messages, \code{2} for a diagnostic mode
#' detailing each imputation iteration, \code{3} for the graphical mode.
#'
#' @return A list with each element being an improper-imputed data set.
impute.ell.improper <- function(X, num = 5, iter.burnin = 10, iter.skip = 10,
                                verbosity = 1){
  if (verbosity >= 2){
    cat("Improper imputation for elliptical family:\n")
  }
  # Perform initial imputation with the single imputation
  if (verbosity >= 2){
    cat("Initial (single) imputation with Mahalanobis depth...\n")
  }
  X.prep <- X
  n <- nrow(X.prep)
  d <- ncol(X.prep)
  miss.label <- is.na(X.prep)
  if (any(rowSums(miss.label) > d - 0.5)){
    stop("There are rows with all missing entries, cannot proceed.")
  }
  miss = which(is.na(X), arr.ind=T)
  X.prep <- impute.depth(X, depth = "Mahalanobis", verbosity = verbosity - 1)
  #if (verbosity >= 2){
  #  cat("done.\n")
  #}
  # Diagnostics (1)
  # means <- NULL
  # cors <- NULL
  Xs <- list("") # structure for imputed data sets
  if (verbosity >= 1){
    cat("Burn-in ")
  }
  for (iter in 1:(iter.burnin + (iter.skip + 1) * (num - 1) + 1)){
    # Calculate points' depths
    mu.tmp <- colMeans(X.prep)
    X.prep <- t(t(X.prep) - mu.tmp)
    Sigma <- cov(X.prep)
    Inv.Sigma <- solve(Sigma)
    depths <- sort(1/(1 + (X.prep %*% Inv.Sigma * X.prep) %*% rep(1, d)))
    # Prepare structures for imputation of points at once
    miss.rowi <- unique(miss[,1])
    X.new <- X.prep
    X.new[miss] <- NA
    X.new <- X.new[miss.rowi,]
    X.prep[miss.rowi,] <- t(apply(X.new, 1, imputeEllPRnd, Sigma, depths))
    X.prep <- t(t(X.prep) + mu.tmp)
    # Diagnostics (2)
    # means <- rbind(means, colMeans(X.prep))
    # cors <- rbind(cors, as.vector(cor(X.prep)))
    if (verbosity >= 2){
      cat(".")
    }
    if (iter == 1){
      if (verbosity >= 3){
        plot(cbind(depths, 1:length(depths)/length(depths)), type = "l",
             lwd = 2, col = "red", xlim = c(0, 1), ylim = c(0, 1),
             main = "Depth empirical distribution",
             xlab = "Depth", ylab = "Cumlative probability")
        grid()
        #lines(cbind(depths, 1:length(depths)/length(depths)), col = "red")
      }
    }else{
      if (iter == iter.burnin){
        if (verbosity >= 1){
          cat(" passed")
        }
        if (verbosity >= 2){
          cat("\n")
        }
      }
      if ((iter > iter.burnin) &&
          ((iter - iter.burnin) %% (iter.skip + 1) == 1)){
        iter.sample <- floor((iter - iter.burnin) / (iter.skip + 1)) + 1
        Xs[[iter.sample]] <- X.prep
        if (verbosity == 1){
          cat(paste(" (smpl ", iter.sample, ")", sep = ""))
          if (iter.sample %% 25 == 0){
            cat("\n")
          }
        }
        if (verbosity >= 2){
          cat(paste(" sample ", iter.sample, " taken.\n", sep = ""))
        }
        if (verbosity >= 3){
          lines(cbind(depths, 1:length(depths)/length(depths)), col = "orange",
                lwd = 2)
        }
      }else{
        if (verbosity >= 3){
          lines(cbind(depths, 1:length(depths)/length(depths)), col = "green")
        }
      }
    }
  }
  # Diagnostics (3)
  # res <- list(Xs = Xs, means = means, cors = cors)
  if (verbosity == 1){
    cat(".\n")
  }
  # Diagnostics (4)
  # return (res)
  return (Xs)
}

#' Multiple imputation for elliptical family
#'
#' Depth-based multiple imputation for an elliptically symmetric distribution.
#'
#' @param X Data set with missing entries.
#'
#' @param m Number of multiply-imputed data sets.
#'
#' @param iter.burnin Number of burn-in iterations.
#'
#' @param verbosity Verbosity parameter: \code{0} for no messages at
#' all, \code{1} (default) for minimal messages, \code{2} for a diagnostic mode
#' detailing each imputation iteration, \code{3} for the graphical mode.
#'
#' @return A list with each element being a multiply imputed data set.
impute.ell.multiple <- function(X, m = 5, iter.burnin = 10, verbosity = 1){
  if (verbosity >= 2){
    cat("Multiple imputation (m = ", m, ") for elliptical family:\n", sep = "")
  }
  n <- nrow(X)
  d <- ncol(X)
  miss = which(is.na(X),arr.ind=T)
  miss.label <- is.na(X)
  if (any(rowSums(miss.label) > d - 0.5)){
    stop("There are rows with all missing entries, cannot proceed.")
  }
  # Perform initial imputation with the single imputation
  if (verbosity >= 2){
    cat("Initial (single) imputation with Mahalanobis depth...\n")
  }
  X.prep.tmp <- impute.depth(X, depth = "Mahalanobis",
                             verbosity = verbosity - 1)
  #if (verbosity >= 2){
  #  cat("done.\n")
  #}
  Xs <- list("") # structure for imputed data sets
  if (verbosity == 1){
    cat("Multiple imputation", sep = "")
  }
  for (iter in 1:m){
    X.prep <- X.prep.tmp
    if (verbosity >= 2){
      cat("Drawing multiple sample ", iter, ":\n", sep = "")
    }
    # Draw bootstrap sample
    if (verbosity >= 2){
      cat("Drawing bootstrap sample...")
    }
    sampleBi <- sample(x = 1:n, size = n, replace = TRUE)
    if (verbosity >= 2){
      cat("done.\n")
    }
    if (verbosity >= 2){
      cat("Burn-in ")
    }
    for (i in 1:iter.burnin){
      mu.tmp <- colMeans(X.prep[sampleBi,])
      Sigma <- cov(X.prep[sampleBi,])
      Inv.Sigma <- solve(Sigma)
      X.prep <- t(t(X.prep) - mu.tmp)
      depths <- sort(1/(1 + (X.prep %*% Inv.Sigma * X.prep) %*% rep(1, d)))
      miss.rowi <- unique(miss[,1])
      X.new <- X.prep
      X.new[miss] <- NA
      X.new <- X.new[miss.rowi,]
      X.prep[miss.rowi,] <- t(apply(X.new, 1, imputeEllPRnd, Sigma, depths))
      X.prep <- t(t(X.prep) + mu.tmp)
      # Diagnostics
      if (i == 1){
        if (verbosity >= 3){
          plot(cbind(depths, 1:length(depths)/length(depths)), type = "l",
               lwd = 2, col = "red", xlim = c(0, 1), ylim = c(0, 1),
               main = "Depth empirical distribution",
               xlab = "Depth", ylab = "Cumlative probability")
          grid()
          #lines(cbind(depths, 1:length(depths)/length(depths)), col = "red")
        }
        # lines(cbind(depths, 1:length(depths)/length(depths)), col = "red")
      }else{
        if (i == iter.burnin){
          if (verbosity >= 2){
            cat(" passed.\n")
          }
          if (verbosity >= 3){
            lines(cbind(depths, 1:length(depths)/length(depths)),
                  col = "orange", lwd = 2)
          }
          # lines(cbind(depths, 1:length(depths)/length(depths)),
          #       col = "orange")
        }else{
          if (verbosity >= 2){
            cat(".")
          }
          if (verbosity >= 3){
            lines(cbind(depths, 1:length(depths)/length(depths)), col = "green")
          }
          # lines(cbind(depths, 1:length(depths)/length(depths)),
          #       col = "green")
        }
      }
    }
    Xs[[iter]] <- X.prep
    if (verbosity >= 2){
      cat(paste("Multiple sample ", iter, " taken.\n", sep = ""))
    }
    if (verbosity == 1){
      cat(" ", iter, sep = "")
    }
  }
  if (verbosity == 1){
    cat(".\n", sep = "")
  }
  return (Xs)
}
