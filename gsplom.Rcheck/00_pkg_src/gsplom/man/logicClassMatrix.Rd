\name{logicClassMatrix}
\alias{logicClassMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Compute Logical Dependency Class Matrix.
}
\description{
  Efficiently compute an all-pairs-rows logical dependency class
  matrix from a data matrix.
}
\usage{
logicClassMatrix(M, thresholds=NULL, intervals=NULL, z=3, min.cnt=0, frac.conf=0.2)
}
\arguments{
  \item{M}{
     A numeric matrix where rows are variables and columns
     are samples.
   }
   \item{thresholds}{
     A vector of high/low thresholds corresponding to each row of \code{M}.
     Default is NULL. If NULL, then thresholds are computed automatically
     using \code{adaptiveRegress}.
   }
   \item{intervals}{
     A vector of uncertainty intervals corresponding to each row of
     \code{M}. Default is NULL. If NULL, then intervals are uniformly
     set to the 3rd percentile row standard deviation.
   }
   \item{z}{
     A parameter of how sparse a quadrant must be
     relative to the margins to be called as empty. Higher z
     means that a quadrant must be more relatively sparse to be
     called as sparse. Default is 3.
   }
   \item{min.cnt}{
     A parameter to force a quadrant to be called as sparse if a fewer or
     equal number of samples are in it.
   }
   \item{frac.conf}{
     A parameter of what fraction of points must be in a quadrant and
     not in the uncertainty region to not classifiy as logical
     dependency class \emph{NA}.
   }
}
%\details{
%}
\value{
  A numeric matrix where the entry at row-i, column-j is an integer
  enumeration of a logical dependency class between row-i and row-j
  in \code{M}. Rows and columns are labeled by the row labels of \code{M}.
}
\note{
  Supports missing values. Missing values in either variable per pair
  are called to be in the uncertainty region and ignored during classification.
  }
\author{
  Andrew D. Yates
}
\examples{
M <- t(LifeCycleSavings)
CLS <- logicClassMatrix(M, intervals=apply(M,1,sd)/2, frac.conf=0.2)
M[1,1] <- NA
CLS.NA <- logicClassMatrix(M, intervals=apply(M,1,function(x)sd(x,na.rm=TRUE))/2, frac.conf=0.2)
}
