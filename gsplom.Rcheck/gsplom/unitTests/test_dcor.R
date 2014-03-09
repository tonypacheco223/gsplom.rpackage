library(energy)

test_common_use <- function() {
  M <- matrix(4,20,20)
  D <- dcorMatrix(M)
}

test_dcor_single <- function() {
  for (i in 1:100) {
    x <- rnorm(30)
    y <- rnorm(30)
    checkEqualsNumeric(dcor(x,y), dcorSingle(x,y), tolerance=1.0e-6)
  }
}
test_dcor_single_NA <- function() {
  checkTrue(is.na(dcorSingle(c(1),c(1))))
  checkTrue(is.na(dcorSingle(c(),c())))
}

test_all_dcor_matrix <- function() {
  M <- matrix(rnorm(1600),20,80)
  D1 <- dcorMatrix(M)
  R2 <- dcorMatrixNA(M)
  D3 <- matrix(0, 20, 20)
  for (i in 1:20) {
    for (j in 1:20) {
      D3[i,j] <- dcor(M[i,],M[j,])
    }
  }
  checkTrue(all(abs(D1-R2$DCOR)<0.000001))
  checkTrue(all(abs(D3-D1)<0.000001))
}

test_dcor_missing <- function() {
  M <- matrix(rnorm(400),20,20)
  M[1,1] <- NA
  R <- dcorMatrixNA(M)
  checkEquals(R$DCOR[1,1],1)
  checkEquals(R$DCOR[2,2],1)
  checkEquals(R$DCOR[1,2],R$DCOR[2,1])
  checkEquals(R$SIZE[1,1]+1,R$SIZE[2,2])
  d <- dcor(M[1,2:20], M[2,2:20])
  checkEqualsNumeric(d, R$DCOR[1,2], tolerance=1.0e-6)
}
