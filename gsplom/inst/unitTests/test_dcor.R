test_common_use <- function() {
  M <- matrix(4,20,20)
  D <- dcorMatrix(M)
}

test_dcor_single_NA <- function() {
  checkTrue(is.na(dcorSingle(c(1),c(1))))
  checkTrue(is.na(dcorSingle(c(),c())))
}

test_all_dcor_matrix <- function() {
  M <- matrix(rnorm(1600),20,80)
  D1 <- dcorMatrix(M)
  R2 <- dcorMatrixNA(M)
  checkTrue(all(abs(D1-R2$DCOR)<0.000001))
}

test_dcor_missing <- function() {
  M <- matrix(rnorm(400),20,20)
  M[1,1] <- NA
  R <- dcorMatrixNA(M)
  checkEquals(R$DCOR[1,1],1)
  checkEquals(R$DCOR[2,2],1)
  checkEquals(R$DCOR[1,2],R$DCOR[2,1])
  checkEquals(R$SIZE[1,1]+1,R$SIZE[2,2])
}
