%% ==================================================
\name{dcorSingle}
\alias{dcorSingle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Single Pair Distance Correlation.
}
\description{
  Compute Distance Correlation between two 1-dimensional vectors.
}
\usage{
dcorSingle(x, y)
}
\arguments{
  \item{x}{
     A numeric vector.
   }
  \item{y}{
     A numeric vector of the same length as \code{x}.
   }
}
%\details{
%}
\value{
  Numeric value between 0 and 1 of the distance correlation between \code{x}
  and \code{y}. Returns \code{NA} if distance correlation is undefined.
}
\author{
  Andrew D. Yates
}
\note{
  Uses absolute difference (Euclidan Distance in one dimension) as a
  distance function in the distance correlation computation. Does not
  support missing values; remove samples with at least one missing value
  from \code{x} and \code{y} prior to using \code{dcorSingle}.
  }
\examples{
## Load data so that samples (countries) are columns and econometric
##   variables are rows.
M <- t(LifeCycleSavings)
## Compute distance correlation between two econometric variables.
d <- dcorSingle(M[1,], M[2,])
}
