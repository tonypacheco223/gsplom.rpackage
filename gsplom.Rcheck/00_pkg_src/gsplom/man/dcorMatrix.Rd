\name{dcorMatrix}
\alias{dcorMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Compute Distance Correlation Matrix.
}
\description{
  Efficiently compute an all-pairs-rows distance correlation matrix from a
  data matrix. 
}
\usage{
dcorMatrix(M, verbose = TRUE)
}
\arguments{
  \item{M}{
     A numeric matrix where rows are variables and columns
     are samples with no missing values (\code{NA}). 
   }
  \item{verbose}{
     Whether to output status information as the result is computed. Default is TRUE.
   }
}
%\details{
%}
\value{
  A numeric matrix where the entry at row-i, column-j is the distance
  correlation between row-i and row-j in \code{M}. Rows and columns are
  labeled by the row labels of \code{M}.
}
\author{
  Andrew D. Yates
}
\note{
  Uses absolute difference (Euclidan Distance in one dimension) as a
  distance function in the distance correlation computation. Does not
  support missing values (\code{NA}); see \code{dcorMatrixNA}.
  }
\examples{
## Load data so that samples (countries) are columns and econometric
##   variables are rows.
M <- t(LifeCycleSavings)
## Compute all pairs distance correlation between econometric variables.
DCOR <- dcorMatrix(M)
}
