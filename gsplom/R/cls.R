##########################################################################
# Copyright (c) 2014 Andrew Yates
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
##########################################################################

percentileStd <- function(M, percentile=0.03) {
  # All-equal uncertainty interval from 3% percentile row sd
  rep(quantile(apply(M,1,sd), percentile), nrow(M))
}

# Partition ranks by minimum regression error
mse.f <- function(s,ss,n) ss - s**2/n
adaptiveRegress <- function(v) {
  x <- sort(v)
  x2 <- x**2
  k <- 0; mse <- Inf
  sum.low <- 0; sqsum.low <- 0
  sum.high <- sum(x); sqsum.high <- sum(x2)
  n <- length(v)
  stopifnot(n >= 2)

  for (i in 1:(n-1)) {
    sum.low <- sum.low + x[i]
    sqsum.low <- sqsum.low + x2[i]
    sum.high <- sum.high - x[i]
    sqsum.high <- sqsum.high - x2[i]
    mse.i <- mse.f(sum.low, sqsum.low, i) + mse.f(sum.high, sqsum.high, n-i)
    if (mse.i < mse) {
      k <- i; mse <- mse.i
    }
  }
  thresh <- (x[k]+x[k+1])/2
  list("thresh"=thresh, "k"=k)
}

logicClassMatrix <- function(M, thresholds=NULL, intervals=NULL, z=3, min.cnt=0, frac.conf=0.2) {
  # All pairs rows.
  if (is.null(thresholds)) {
    thresholds <- apply(M, 1, function(v) adaptiveRegress(v)[[1]])
  }
  if (is.null(intervals)) {
    intervals <- percentileStd(M)
  }
  # Quantize each row by low, uncertain, high
  LUH <- matrix(0, nrow(M), ncol(M))
  for (i in 1:nrow(LUH)) {
    LUH[i,M[i,] < thresholds[i]-intervals[i]] <- -1
    LUH[i,M[i,] > thresholds[i]+intervals[i]] <- 1
  }
  # Counts per quadrant: Q_{row}_{col}
  QLL <- (LUH==-1) %*% t(LUH==-1) # row low, col low
  QLH <- (LUH==-1) %*% t(LUH==1)  # row low, col high
  QHL <- (LUH==1) %*% t(LUH==-1)  # row high, col low
  QHH <- (LUH==1) %*% t(LUH==1)   # row high, col high

  RL <- QLL + QLH               # row variable (y-axis) is low
  RH <- QHL + QHH               # row variable (y-axis) is high
  CL <- QLL + QHL               # col variable (x-axis) is low
  CH <- QLH + QHH               # row variable (x-axis) is high
  ALL <- QLL + QLH + QHL + QHH

  # to avoid division by zero, set each margin to at least 1
  ZERO.MARGIN = (RL==0)|(RH==0)|(CL==0)|(CH==0)
  RL[RL==0] <- 1
  RH[RH==0] <- 1
  CL[CL==0] <- 1
  CH[CH==0] <- 1
  ALL[ALL==0] <- 1
  TOO.FEW <- ALL < frac.conf*ncol(M)

  # Dense/Sparse tests
  SLL <- getSparse(QLL, RL, CL, ALL, z=z, min.cnt=min.cnt)
  SLH <- getSparse(QLH, RL, CH, ALL, z=z, min.cnt=min.cnt)
  SHL <- getSparse(QHL, RH, CL, ALL, z=z, min.cnt=min.cnt)
  SHH <- getSparse(QHH, RH, CH, ALL, z=z, min.cnt=min.cnt)

  # Assign class enumerations.
  CLS <- matrix(0, nrow(M), nrow(M))
  CLS[!SLL & SLH & !SHL & !SHH]  <- 1
  CLS[!SLL & SLH & SHL & !SHH]   <- 2
  CLS[!SLL & !SLH & SHL & !SHH]  <- 3
  CLS[!SLL & !SLH & !SHL & !SHH] <- 4
  CLS[!SLL & !SLH & !SHL & SHH]  <- 5
  CLS[SLL & !SLH & !SHL & SHH]   <- 6
  CLS[SLL & !SLH & !SHL & !SHH]  <- 7
  CLS[ZERO.MARGIN] <- 0
  CLS[TOO.FEW] <- 0
  rownames(CLS) <- colnames(CLS) <- rownames(M)
  CLS
}

getSparse <- function(Q, Row, Col, All, z=3, min.cnt=0) {
  Exp <- Row*Col/All
  Z <- (Exp-Q) / sqrt(Exp)  # marginal probability test
  D <- Q <= min.cnt         # number of points test
  (Z > z) | D               # sparse matrix for single quadrant for all pairs
}
