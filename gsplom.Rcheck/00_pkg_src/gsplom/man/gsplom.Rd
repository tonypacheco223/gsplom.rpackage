\name{gsplom}
\alias{gsplom}
\alias{gsplomCore}
\alias{gsplomPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Glyph SPLOM Workflow Wrapper
}
\description{
  Generate a Glyph SPLOM from a data matrix.
}
\usage{
gsplom(M, \dots)
}
\arguments{
  \item{M}{
     A numeric matrix where rows are variables and columns
     are samples.
   }
   \item{...}{Special plotting options.}
}
%\details{
%}
\value{
  A list of intermediate results and the gsplom plot itself sent to the
  current open plotting device.
}
\author{
  Andrew D. Yates
}
\note{
  This is a high level function that completes the entire Glyph SPLOM
  workflow on a data matrix using the recommended default values. See
  \code{gsplomCore} to plot a Glyph SPLOM from a pre-computed all-pairs
  distance correlation matrix and logical dependency class matrix. Other
  parameters are related to special Glyph SPLOM options.
  }
  \examples{
\dontrun{
## Load data so that samples (countries) are columns and econometric
##   variables are rows.
M <- t(LifeCycleSavings)
pdf("/mypath/gsplom.pdf")
R <- gsplom(M)
dev.off()
}
}
