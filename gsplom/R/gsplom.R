##########################################################################
# Copyright (c) 2014 Andrew Yates
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
##########################################################################
source("R/clsDist.R")
library("fastcluster")

# Generated using make.color.bins()
#   core colors: GLYPH.COLS <- c("#ffffff", "#00b271", "#0089d9", "#3424b3", "#000000", "#a40e78", "#d82a36", "#eb5918", "#ffffff")
#   N: 15
GLYPH.COLORS <- c("#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#EEF9F5","#DDF4EC","#CCEFE2","#BBEAD9","#AAE5CF","#99E0C6","#88DBBC","#77D5B3","#66D0A9","#55CBA0","#44C696","#32C18D","#21BC83","#10B77A","#00B271","#FFFFFF","#EEF7FC","#DDEFF9","#CCE7F7","#BBDFF4","#AAD7F2","#99CFEF","#88C7ED","#77C0EA","#66B8E8","#55B0E5","#44A8E3","#32A0E0","#2198DE","#1090DB","#0089D9","#FFFFFF","#F1F0F9","#E3E1F4","#D6D3EF","#C8C4EA","#BBB6E5","#ADA7E0","#A098DB","#928AD6","#857BD1","#776DCC","#6A5EC7","#5C4FC2","#4F41BD","#4132B8","#3424B3","#FFFFFF","#EEEEEE","#DDDDDD","#CCCCCC","#BBBBBB","#AAAAAA","#999999","#888888","#777777","#666666","#555555","#444444","#323232","#212121","#101010","#000000","#FFFFFF","#F8EEF6","#F2DEED","#ECCEE3","#E6BEDB","#E0AED2","#DA9EC9","#D48EC0","#CE7EB7","#C86EAE","#C25EA5","#BC4E9C","#B63E93","#B02E8A","#AA1E81","#A40E78","#FFFFFF","#FCF0F1","#F9E2E4","#F7D4D6","#F4C6C9","#F2B8BC","#EFA9AE","#EC9BA1","#EA8D93","#E77F86","#E57179","#E2626B","#DF545E","#DD4650","#DA3843","#D82A36","#FFFFFF","#FDF3EF","#FCE8E0","#FBDDD0","#F9D2C1","#F8C7B2","#F7BCA2","#F5B193","#F4A683","#F39B74","#F19065","#F08555","#EF7A46","#ED6F36","#EC6427","#EB5918","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF")

gsplom <- function(M, ...) {
  message("Computing all-pairs-rows Distance Correlation...")
  if (any(is.na(M))) {
    message("Input M has missing values. Using dcorMatrixNA.")
    DCOR <- dcorMatrixNA(M)
  } else {
    DCOR <- dcorMatrix(M)
  }
  message("Computing all-pairs-rows Logical Dependency Class...")
  CLS  <- logicClassMatrix(M)
  R <- gsplomCore(DCOR, CLS, ...)
  R$CLS <- CLS
  R$DCOR <- DCOR
  R
}

gsplomCore <- function(DCOR, CLS, doLabels=TRUE, ...) {
  R <- list()
  message("Computing dCor distance matrix...")
  DCOR.Dist <- as.matrix(dist(DCOR))
  message("Computing Logical Dependency Class distance matrix...")
  CLS.Dist <- logicClassDist(CLS)
  
  DIST <- DCOR.Dist + sqrt(CLS.Dist)/2
  dCorRowMeans <- rowMeans(DCOR, na.rm=TRUE)
  R$Rhclust <- as.dendrogram(hclust(as.dist(DIST), method="average"))
  R$Rhclust <- reorder(R$Rhclust, dCorRowMeans)
  R$order <- order.dendrogram(R$Rhclust)
  if (doLabels) {
    labels <- rownames(DCOR)
  } else {
    labels <- NA
  }
  gsplomPlot(DCOR, CLS, R$order, labels=labels, ...)
  R
}

gsplomPlot <- function(DCOR, CLS, order, labels, MIN=0.2, MAX=1, asGlyphs=TRUE, doDrawLabs=TRUE, cex=0.7, ...) {
  breaks <- makeBreaks(MAX, MIN)
  offsets <- makeOffsets(breaks)
  ## Select the greatest offset less than x.
  choose.offset <- function(x) offsets[tail(which(breaks<=x),1)]
  N <- 15
  N.BREAKS <- c(sapply(0:8, function(x) rep(x,N+1)+breaks[1:(N+1)]), 9)
  ## Make image heatmap from matrices
  OFFSET <- apply(DCOR, c(1,2), choose.offset)
  G <- (CLS+OFFSET)[order,order]
  if (asGlyphs) G <- expandCls(G)
  w <- ncol(G); h <- nrow(G)
  Img <- t(G)[,seq(h,1,-1)]

  ## Get labels
  if(!is.null(labels)) {
    labels <- labels[order]
    labCol <- labels
    labRow <- rev(labels)
    if (asGlyphs) {
      labCol <- sapply(1:(length(labels)*2), function(i) expandNames(i,labCol))
      labRow <- sapply(1:(length(labels)*2), function(i) expandNames(i,labRow))
    } 
  }

  image(1:w, 1:h, Img, xlab="", ylab="", col=GLYPH.COLORS, breaks=N.BREAKS, axes=FALSE, ...)
  if(!is.null(labels)) {
    axis(1, 1:w, labels = labCol, las = 2, line = -0.5, tick = 0, cex.axis = cex)
    axis(4, 1:h, labels = labRow, las = 2, line = -0.5, tick = 0, cex.axis = cex)
  }
}

makeBreaks <- function(MAX=1, MIN=0, N=15, MOST=1, LEAST=0) {
  th <- (MAX-MIN)/N   ## bin for above and below bin threshold
  c(LEAST, sapply(0:(N-1), function(i) MIN+i*th), MOST)
}
makeOffsets <- function(breaks, N=15) {
  offsets <- sapply(1:(N+1), function(i) (breaks[i]+breaks[i+1])/2)
  offsets <- c(offsets, tail(offsets,1)) ## offset beyond max value is still max value
  offsets
}
expandCls <- function(CLS, pad=FALSE, bg=NA) {
  G <- matrix(bg, nrow(CLS)*2, ncol(CLS)*2)
  for(i in 0:(nrow(CLS)-1)) {
    for(j in 0:(ncol(CLS)-1)) {
      gly <- toGlyph(CLS[i+1,j+1], bg)
      G[(i*2+1):(i*2+2),(j*2+1):(j*2+2)] <- gly
    }
  }
  G
}

### Expand CLS matrix into 2x2 glyphs.
## --------------------
# 0,  1,2,3,  4,  5,6,7
toGlyph <- function(z, bg=NA) {
  r <- NaN
  if(z >= 0 && z < 1)  # NA    (no significant dependency)
    r <- matrix(c(bg,bg,bg,bg), nrow=2)
  if(z >= 1 && z < 2)  # HIH   (high x implies high y)
    r <- matrix(c(z,z,z,bg), nrow=2)
  if(z >= 2 && z < 3)  # PC    (positive correlation)
    r <- matrix(c(bg,z,z,bg), nrow=2)
  if(z >= 3 && z < 4)  # LIL   (low x implies low y)
    r <- matrix(c(bg,z,z,z), nrow=2)
  if(z >= 4 && z < 5)  # UNL   (unspecified non-linear)
    r <- matrix(c(z,z,z,z), nrow=2)
  if(z >= 5 && z < 6)  # HIL   (high x implies low y)
    r <- matrix(c(z,z,bg,z), nrow=2)
  if(z >= 6 && z < 7)  # NC    (negative correlation)   
    r <- matrix(c(z,bg,bg,z), nrow=2)
  if(z >= 7 && z < 8)  # LIH   (low x implies high y) 
    r <- matrix(c(z,bg,z,z), nrow=2)
  r
}

expandNames <- function(i, name.list) { 
  if (i%%2==1) {
    name.list[(i+1)/2]
  } else {
    ""
  }
}
