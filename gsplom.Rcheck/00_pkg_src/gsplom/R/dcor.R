##########################################################################
# Copyright (c) 2014 Andrew Yates
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
##########################################################################

# TODO: consider usage of foreach package for parallelization
dcorMatrix <- function(M, verbose=TRUE) {
  if (sum(is.na(M))>0) {
    stop("dcorMatrix does not support missing values. Use dcorMatrixNA instead.")
  }
  # Compute all pairs-rows distance correlation
  # ========================================
  # U_ij = mean of abs diff of vector x_i from element x_ij (x_i distance matrix row/col means, entry j)
  U <- apply(M, 2, function(v) rowMeans(abs(M-v))) # each execution on column v in M becomes column in U.
  # UU_i = mean of all entries of U_i* (mean of distance matrix for x_i)
  UU <- rowMeans(U) # this doesn't account for NAs
  DCOV <- matrix(0,nrow(M),nrow(M))
  n <- ncol(M)
  
  # upper triangle of DCOV**2 matrix
  for (j in 1:(n-1)) {
    if(verbose && j%%100==0)
      message(sprintf("Computing distance correlation up through sample %d of %d...", j, n))
    A <- abs(M[,(j+1):n]-M[,j]) - U[,(j+1):n] - U[,j] + UU
    DCOV <- DCOV + A %*% t(A)
  }
  # diagonal of DCOV**2 matrix (=distance variance, DVAR**2)
  DCOV <- DCOV*2
  for (j in 1:n) {
    A <- UU - 2*U[,j]
    DCOV <- DCOV + A %*% t(A)
  }
  
  DCOV <- sqrt(DCOV)
  DSTD <- sqrt(diag(DCOV))
  # divide DSTD from both columns and rows.
  DCOR <- sweep((DCOV / DSTD),MARGIN=2,DSTD,FUN="/")
  rownames(DCOR) <- colnames(DCOR) <- rownames(M)
  DCOR
}

dcorSingle <- function(x,y) {
  # Compute euclidean distance correlation between two 1-dimension vectors.
  if (is.null(x)) {
    stopifnot(x==y)
    return (NA)
  }
  stopifnot(length(x)==length(y))
  stopifnot(sum(is.na(x))==0)
  stopifnot(sum(is.na(y))==0)
  n <- length(x)
  if (n <= 1) return (NA)
  dx <- sapply(x, function(a) abs(x-a))
  dy <- sapply(y, function(a) abs(y-a))
  dx.row.means <- apply(dx,1,mean)
  dx.mean <- mean(dx)
  dy.row.means <- apply(dy,1,mean)
  dy.mean <- mean(dy)
  dxc <- t(t(dx-dx.row.means)-dx.row.means) + dx.mean
  dyc <- t(t(dy-dy.row.means)-dy.row.means) + dy.mean
  dcov <- sqrt(sum(dxc * dyc))
  xdvar <- sqrt(sum(dxc * dxc))
  ydvar <- sqrt(sum(dyc * dyc))
  dcor <- dcov / sqrt(xdvar*ydvar)
  dcor
}

dcorMatrixNA <- function(M, do.rank=FALSE, verbose=TRUE) {
  # Compute all-pairs-rows distance correlation in a loop, accounting for missing values.
  # WARNING: this is much slower than dcorMatrix because it is not vectorized.
  DCOR <- matrix(NA, nrow(M), nrow(M))
  SIZE <- matrix(NA, nrow(M), nrow(M))
  for (i in 1:(nrow(M)-1)) {
    if(verbose && (i%%100==0||i==1))
      message(sprintf("Computing distance correlation up through variable %d of %d...", i, nrow(M)))
    for (j in (i+1):nrow(M)) {
      y <- M[i,]
      x <- M[j,]
      mask <- !(is.na(x)|is.na(y))
      x <- x[mask]
      y <- y[mask]
      # warn if masked sample size is less than 20.
      stopifnot(length(x)==length(y))
      if (length(x)<20)
        warning(sprintf("Masked sample size for row %d, column %d is %d < 20.", i, j, length(x)))
      if (do.rank) {
        x <- rank(x)
        y <- rank(y)
      }
      d <- dcorSingle(x,y)
      DCOR[i,j] <- DCOR[j,i] <- d
      SIZE[i,j] <- SIZE[j,i] <- length(x)
    }
  }
  diag(DCOR) <- 1
  diag(SIZE) <- apply(M,1,function(x) sum(!is.na(x)))
  rownames(DCOR) <- colnames(DCOR) <- rownames(M)
  rownames(SIZE) <- colnames(SIZE) <- rownames(M)
  list(DCOR=DCOR, SIZE=SIZE)
}
