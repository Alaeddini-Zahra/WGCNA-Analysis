setwd("C:/Users/Atefeh/Desktop/ovarian")

# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
enableWGCNAThreads()


# Load the expression and trait data saved in the first part
lnames <- load(file = "Tumor.RData");
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames <- load(file = "Tumormodules.RData");
lnames

data.Expression.Tumor <- as.data.frame(data.Expression.Tumor)

###############################################
###########Exporting to Cytoscape##############
###############################################
# Recalculate topological overlap if needed
TOM <- TOMsimilarityFromExpr(data.Expression.Tumor, power = 12)
colnames(TOM) <- colnames(data.Expression.Tumor)
rownames(TOM) <- colnames(data.Expression.Tumor)
TOM1 <- as.data.frame(TOM)
TOM <- as.matrix(TOM1)
# Select modulesTOM# Select modules

#modules = c("darkmagenta","darkorange","floralwhite","greenyellow","lightcyan","magenta", "orange","plum1","sienna3","yellowgreen","gold","paleturquoise","skyblue3");

modules = c('black')
# Select module probes
probes <- names(data.Expression.Tumor)
inModule <- is.finite(match(moduleColors1, modules));
modProbes <- probes[inModule];

                               
# Select the corresponding Topological Overlap
modTOM <- TOM[inModule, inModule];



dimnames(modTOM) <- list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt <- exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-Tumor-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-Tumor-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeAttr = moduleColors1[inModule]);

