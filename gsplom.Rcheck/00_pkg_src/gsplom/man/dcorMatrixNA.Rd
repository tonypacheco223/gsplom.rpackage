\name{dcorMatrixNA}
\alias{dcorMatrixNA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Compute Distance Correlation Matrix with Missing Values in Data.
}
\description{
  A slower method of computing all-pairs-rows distance correlation
  that excludes samples with at least one missing value in each computation.
}
\usage{
dcorMatrixNA(M, do.rank = FALSE, verbose = TRUE)
}
\arguments{
  \item{M}{
     A numeric matrix where rows are variables and columns
     are samples. May have missing values (\code{NA}). 
   }
  \item{do.rank}{
     Whether to convert values to ranks after removing samples with
     missing values. Default is FALSE.
   }
  \item{verbose}{
     Whether to output status information as the result is computed. Default is TRUE.
   }
}
%\details{
%}
\value{
  A list containing two matrices:
  
  \code{DCOR}: A numeric matrix where the entry at row-i, column-j is the distance
  correlation between row-i and row-j in \code{M} after removing samples
  where there is at least one missing value in either row-i or
  row-j. If only one or zero samples remain after removing missing
  values, the entry value is \code{NA}. Rows and columns are labeled by
  the row labels of \code{M}.
  
  \code{SIZE}: An integer matrix where the entry at row-i, column-j is
  the number of samples used to compute DCOR[i,j] after handling missing
  values.
}
\author{
  Andrew D. Yates
}
\note{
  Uses absolute difference (Euclidan Distance in one dimension) as a
  distance function in the distance correlation computation. This
  function is much less efficient than \code{dcorMatrix} and should
  only be used when missing values are in the data.
  }
\examples{
## Load data so that samples (countries) are columns and econometric
##   variables are rows.
M <- t(LifeCycleSavings)
## Insert a missing value.
M[1,1] <- NA
## Compute all-pairs distance correlation between econometric variables.
Results <- dcorMatrixNA(M)
}
