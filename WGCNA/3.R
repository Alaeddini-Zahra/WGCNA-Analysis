setwd("C:/Users/Atefeh/Desktop/data")
library(WGCNA)
library(flashClust)

##### change the global default setting. So every data frame you create after executing
##### that line will not auto-convert to factors unless explicitly told to do so.
##### not change string as factors
options(stringsAsFactors = FALSE)

enableWGCNAThreads()
# Load the data saved in the first part
lnames <- load(file = "Tumor2.RData")
#The variable lnames contains the names of loaded variables..
lnames

###############################################
#####Choosing the soft-thresholding power######
###############################################
# Choose a set of soft-thresholding powers
powers <- c(c(1:10), seq(from = 12, to=50, by=2))
# Call the network topology analysis function to detect the best soft thresholding powers (beta)
sft <- pickSoftThreshold(data.Expression.Tumor, powerVector = powers, verbose = 5, networkType = "signed hybrid")
#sft <- as.data.frame(sft)
#sft <- sft[order(sft$fitIndices.SFT.R.sq),]
#write.csv(sft, file = "Tumor_power.csv", quote = F, row.names = F)
#pwr <- read.csv('Tumor_power.csv')
#View(pwr)

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));

text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.9,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

###############################################
#####Co-expression similarity and adjacency####
###############################################
softPower = 12
adjacency = adjacency(data.Expression.Tumor, type = "signed hybrid", power = softPower);


###############################################
############Topological Overlap Matrix#########
###############################################
# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM

###############################################
############Clustering using TOM###############
###############################################
# Call the hierarchical clustering function
geneTree1 = flashClust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree1, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree1, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)
# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree1, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

###############################################
##############Merging of modules###############
###############################################
# Calculate eigengenes
MEList = moduleEigengenes(data.Expression.Tumor, colors = dynamicColors)
MEs1 = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs1);
# Cluster module eigengenes
METree = flashClust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.8 ###########################################
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(data.Expression.Tumor, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

sizeGrWindow(12, 9)
#pdf(file = "geneDendro-7.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree1, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#dev.off()
table(mergedColors)
# Rename to moduleColors
moduleColors1 = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels1 = match(moduleColors1, colorOrder)-1;
MEs1 = mergedMEs;
# Save module colors and labels for use in subsequent parts
save(MEs1, moduleLabels1, moduleColors1, geneTree1, file = "Tumormodules.RData")


