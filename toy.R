library(gsplom)
# Load a sample dataset built into R
M <- t(LifeCycleSavings)

# Vector image with labels
# WARNING: vector images with labels can be unwieldily to render to screen for large (10k+) GSPLOMs!
#  consider rendering single pixel variants at scale.
pdf("~/Desktop/gsplom_vector_glyphs.pdf", width=10, height=10) # dimensions in inches
gsplomResults <- gsplom(M) 
dev.off()

# gsplomResults is a list containing information about how the glyph splom was generated.
# CONTENTS:
#   gsplom$Rhclust: a dendrogram object from hierarchical clustering.
#   gsplom$order: the reordering of rows and columns from the original data matrix M
#   gsplom$CLS: all-pairs-rows logical dependency class enumerations in order of the original data matrix M
#   gsplom$DCOR: all-pairs-rows distance correlations in order of the original data matrix M

# Class enumerations (see Eurovis publication)
#   1: YiX, 2: PC, 3: XiY, 4: UNL, 5: MX, 6: NC, 7: OR, 0: NA

# Some useful things to do with gsplomResults:
# ------------------------------
# Save entire result object to file as an R binary object
save(gsplomResults, file="myResults.Rdata")
# Reload those saved results from a file into memory.
myResultsFromDisk <- load("myResults.Rdata")
# Save DCOR and CLS matrices as csv files that can be opened in Excel.
write.csv(gsplomResults$DCOR, file="~/Desktop/DCOR.csv")
write.csv(gsplomResults$CLS, file="~/Desktop/CLS_ENUM.csv")
# Save DCOR and CLS matrices so that they align with the GSPLOM image.
write.csv(gsplomResults$DCOR[gsplomResults$order, gsplomResults$order], file="~/Desktop/DCOR_ordered.csv")
write.csv(gsplomResults$CLS[gsplomResults$order, gsplomResults$order], file="~/Desktop/CLS_ENUM_ordered.csv")
# Plot hierarchical clustering dendrogram
pdf("~/Desktop/myDendrogram.pdf", width=8, height=4) # units in inches
plot(gsplomResults$Rhclust)
dev.off()

# Additional plotting examples using bitmap images
# ------------------------------

# Rasterized, pixel perfect images
# Single pixel image, colors only, no glyphs or labels
png("~/Desktop/pix_glyph.png", width=nrow(M), height=nrow(M)) # dimensions in pixels
par(mar = c(0,0,0,0))
R <- gsplom(M, doLabels=FALSE, useRaster=TRUE)
dev.off()

# 2x2 pixel image with glyphs
png("~/Desktop/pix_single.png", width=nrow(M)*2, height=nrow(M)*2) # dimensions in pixels
par(mar = c(0,0,0,0))
R <- gsplom(M, doLabels=FALSE, useRaster=TRUE, asGlyphs=FALSE)
dev.off()

# Don't need the gsplom plot, only logical classes and distance correlations?
# ------------------------------
CLS <- logicClassMatrix(M)
DCOR <- dcorMatrix(M)
