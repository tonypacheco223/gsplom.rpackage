##########################################################################
# Copyright (c) 2014 Andrew Yates
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
##########################################################################

# Glyph Distance to NA class.
NAdist <- 2 # round 2.05 to 2 to keep integer values
# Convert class enumerations to binary representations
#   Bit enumerations:
#     2 | 1
#     -----
#     3 | 4
enumToBits <- function(CLS) {
  BitCLS <- mat.or.vec(nrow(CLS), ncol(CLS))
  BitCLS[CLS==1] <- 14  # 1110
  BitCLS[CLS==2] <- 10  # 1010
  BitCLS[CLS==3] <- 11  # 1011
  BitCLS[CLS==4] <- 15  # 1111
  BitCLS[CLS==5] <- 7   # 0111
  BitCLS[CLS==6] <- 5   # 0101
  BitCLS[CLS==7] <- 13  # 1101
  BitCLS
}
# xor of row to all other rows in M
rowXor <- function(row, M) matrix(bitwXor(row,t(M)),nrow(M),nrow(M),byrow=TRUE)
# Number of bits set in 4bit ming distance
Xor2Hamm <- function(XOR) {
  sigma <- (bitwAnd(XOR,1)>0)+(bitwAnd(XOR,2)>0)+(bitwAnd(XOR,4)>0)+(bitwAnd(XOR,8)>0)
  matrix(sigma,nrow(XOR),nrow(XOR))
}
# Force NA hamming distance to always be NAdist except NA<->NA is 0.
naHamm <- function(row, M, HAMM) {
  rowNA <- row==0
  clsNA <- M==0
  # select columns where row==NA, then select rows where exisiting class is not NA
  HAMM[,rowNA][!clsNA[,rowNA]] <- NAdist
  # select columns where row!=NA, then select rows where existing class is NA
  HAMM[,!rowNA][clsNA[,!rowNA]] <- NAdist
  HAMM
}
# Compute row-to-matrix by row Hamming distance and account for NA classes
rowXorNA <- function(row, M) naHamm(row, M, Xor2Hamm(rowXor(row, M)))
# Sum hamming distances
sumRowHamm <- function(row, M) apply(rowXorNA(row, M), 1, sum)

logicClassDist <- function(CLS) {
  BitCLS <- enumToBits(CLS)
  CLS.Dist <- apply(BitCLS, 1, function(row) sumRowHamm(row, BitCLS))
  rownames(CLS.Dist) <- colnames(CLS.Dist) <- rownames(CLS)
  CLS.Dist
}

# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# TODO: find how to run this in package unit tests without exposing this function
# ------------------------------
# TEST.CLS
# ----------
##   A B C D E
## A 2 2 1 4 4
## B 2 2 3 4 4
## C 3 1 2 7 4
## D 4 4 7 2 0
## E 4 4 4 0 2

# Expected TEST.CLS Distance Matrix
# ----------
##    A  B  C  D E
## A  0  2  4 10 9
## B  2  0  4 10 9
## C  4  4  0 10 8
## D 10 10 10  0 5
## E  9  9  8  5 0

test_logicClassDist <- function() {
  TEST.CLS <- matrix(c(2,2,1,4,4,2,2,3,4,4,3,1,2,7,4,4,4,7,2,0,4,4,4,0,2),5,5,byrow=TRUE)
  rownames(TEST.CLS) <- colnames(TEST.CLS) <- c("A","B","C","D","E")
  D <- logicClassDist(TEST.CLS)
  Dt <- logicClassDist(t(TEST.CLS))
  # Hand compute expected result.
  Expected <- matrix(c(0,2,4,10,9,2,0,4,10,9,4,4,0,10,8,10,10,10,0,5,9,9,8,5,0),5,5,byrow=TRUE)
  # checkTrue(all(D==Expected))
  # checkTrue(all(Dt==Expected))
  stopifnot(all(D==Expected))
  stopifnot(all(Dt==Expected))
}
