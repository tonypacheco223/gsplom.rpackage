library(gsplom)
M <- t(LifeCycleSavings)
CLS <- logicClassMatrix(M)
DCOR <- dcorMatrix(M)
source("gsplom/R/clsDist.R")
CLS.Dist <- logicClassDist(CLS)
DCOR.Dist <- as.matrix(dist(DCOR))
DIST <- DCOR.Dist + sqrt(CLS.Dist)/2

# vector with labels
pdf("~/Desktop/pix_glyph.pdf", width=10, height=10)
gsplom(M)
dev.off()

# rasterized, pixel perfect images
png("~/Desktop/pix_glyph.png", nrow(M), nrow(M))
par(mar = c(0,0,0,0))
gsplom(M, doLabels=FALSE, useRaster=TRUE)
dev.off()

png("~/Desktop/pix_single.png", nrow(M), nrow(M))
par(mar = c(0,0,0,0))
gsplom(M, doLabels=FALSE, useRaster=TRUE, asGlyphs=FALSE)
dev.off()
