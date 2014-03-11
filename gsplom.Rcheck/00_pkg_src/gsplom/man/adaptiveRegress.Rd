\name{adaptiveRegress}
\alias{adaptiveRegress}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  High/Low binary partition by adaptive regression.
}
\description{
  Sort values from low to high, fit a step-up function, and return a
  partition. Removes \code{NA} values before computation.
}
\usage{
adaptiveRegress(v)
}
\arguments{
  \item{v}{
    A numeric vector to quantize.
   }
}
%\details{
%}
\value{
  A list containing two elements:
  
  \code{thresh}: high/low threshold value
  
  \code{k}: highest rank of low partition
}
\author{
  Andrew D. Yates
}
  \examples{
v <- sample(1:100)
R <- adaptiveRegress(v)
}
