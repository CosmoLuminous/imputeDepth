################################################################################
## File:             imputation.R                                             ##
## Created by:       Pavlo Mozharovskyi                                       ##
## Last revised:     23.12.2016                                               ##
##                                                                            ##
## Contains functions for the depth-based singular value decomposition of the ##
## shape matrix.                                                              ##
##                                                                            ##
################################################################################

median.hfsp <- function(X, r = 1000){
  a <- .C("MedianHfsp",
          as.double(t(X)),
          as.integer(nrow(X)),
          as.integer(ncol(X)),
          as.integer(r),
          med = double(ncol(X)))
  return (a$med)
}

median.proj <- function(X, r = 1000){
  a <- .C("MedianProj",
          as.double(t(X)),
          as.integer(nrow(X)),
          as.integer(ncol(X)),
          as.integer(r),
          med = double(ncol(X)))
  return (a$med)
}

bin.searchl <- function(x, val){
  a <- .C("BinSearchlRoutine",
          as.double(x),
          as.integer(length(x)),
          as.double(val),
          pos = integer(1))$pos
  return (a)
}

bin.searchr <- function(x, val){
  a <- .C("BinSearchrRoutine",
          as.double(x),
          as.integer(length(x)),
          as.double(val),
          pos = integer(1))$pos
  return (a)
}

median.qs <- function(x){
  a <- .C("MedianQs", as.double(x), as.integer(length(x)), res = double(1))$res
  return (a)
}

depths.projection <- function(X, r = 1000){
  a <- .C("DepthsProj",
          as.double(t(X)),
          as.integer(nrow(X)),
          as.integer(ncol(X)),
          as.integer(r),
          res = double(nrow(X)))$res
  return (a)
}

depths.halfspace <- function(X, r = 1000){
  a <- .C("DepthsHfsp",
          as.double(t(X)),
          as.integer(nrow(X)),
          as.integer(ncol(X)),
          as.integer(r),
          res = double(nrow(X)))$res
  return (a)
}

depths.spatial <- function(X, sd = TRUE){
  a <- .C("DepthsSptl",
          as.double(t(X)),
          as.integer(nrow(X)),
          as.integer(ncol(X)),
          as.integer(sd),
          res = double(nrow(X)))$res
  return (a)
}

an.svd <- function(X){
  a <- .C("AnSvd",
          as.double(X),
          as.integer(nrow(X)),
          as.integer(ncol(X)),
          s = double(ncol(X)),
          u = double(nrow(X) * ncol(X)),
          v = double(ncol(X) * ncol(X)))
  return (a)
}

get.w <- function(X){
  a <- .C("GetW",
          as.double(t(X)),
          as.integer(nrow(X)),
          as.integer(ncol(X)),
          mu = double(ncol(X)),
          w = double(ncol(X) * ncol(X)))
  return (list(mu = a$mu,
               w = matrix(a$w, nrow = ncol(X), byrow = TRUE)))
}

svd.depth.proj <- function(X, r = 10000000 / nrow(X)){
  a <- .C("SvdDepthProj",
          as.double(t(X)),
          as.integer(nrow(X)),
          as.integer(ncol(X)),
          as.integer(r),
          as.integer(0),
          s = double(ncol(X)),
          v = double(ncol(X) * ncol(X)),
          e = double(1),
          c = double(ncol(X)))
  return (list(s = a$s,
               v = matrix(a$v, nrow = ncol(X), byrow = TRUE),
               e = a$e,
               c = a$c))
}

cov.depth.proj <- function(X, r = 10000000 / nrow(X)){
  svdPrj <- svd.depth.proj(X, r)
  return (svdPrj$v %*%
            diag(svdPrj$s^2 * svdPrj$e^2 / qnorm(.75)^2) %*%
            t(svdPrj$v))
}

svd.depth.hfsp <- function(X, r = 10000000 / nrow(X) / log(nrow(X))){
  a <- .C("SvdDepthHfsp",
          as.double(t(X)),
          as.integer(nrow(X)),
          as.integer(ncol(X)),
          as.integer(r),
          as.integer(0),
          as.integer(nrow(X) / 2),
          s = double(ncol(X)),
          v = double(ncol(X) * ncol(X)),
          e = double(1),
          c = double(ncol(X)))
  return (list(s = a$s,
               v = matrix(a$v, nrow = ncol(X), byrow = TRUE),
               e = a$e,
               c = a$c))
}

cov.depth.hfsp <- function(X, r = 10000000 / nrow(X) / log(nrow(X))){
  svdHsp <- svd.depth.hfsp(X, r)
  return (svdHsp$v %*%
            diag(svdHsp$s^2 * svdHsp$e^2 / qnorm(.75)^2) %*%
            t(svdHsp$v))
}
