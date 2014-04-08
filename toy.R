library(gsplom)
# Load a sample dataset built into R
M <- t(LifeCycleSavings)

# Vector image with labels
# WARNING: vector images with labels can be unwieldily to render to screen for large (10k+) GSPLOMs!
#  consider rendering single pixel variants at scale.
pdf("~/Desktop/pix_glyph.pdf", width=10, height=10) # dimensions in inches
R <- gsplom(M) # R contains 
dev.off()


# Rasterized, pixel perfect images
# Single pixel image, colors only, no glyphs or labels
png("~/Desktop/pix_glyph.png", width=nrow(M), height=nrow(M)) # dimensions in pixels
par(mar = c(0,0,0,0))
gsplom(M, doLabels=FALSE, useRaster=TRUE)
dev.off()

# 2x2 pixel image with glyphs
png("~/Desktop/pix_single.png", width=nrow(M)*2, height=nrow(M)*2) # dimensions in pixels
par(mar = c(0,0,0,0))
gsplom(M, doLabels=FALSE, useRaster=TRUE, asGlyphs=FALSE)
dev.off()

CLS <- logicClassMatrix(M)
DCOR <- dcorMatrix(M)
