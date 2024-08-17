setwd("C:/Users/Atefeh/Desktop/data")
library(WGCNA)
library(flashClust)



##### change the global default setting. So every data frame you create after executing
##### that line will not auto-convert to factors unless explicitly told to do so.
##### not change string as factors

options(stringsAsFactors = FALSE)
################################################## Normal tissue
##### Recall data and save Normal data samples' name into an array
Normal_data <- read.csv("Normal_samples.csv")
array_name <- names(Normal_data)

##### transpose the Normal samples data frame and remove the Gene.symbols row
data <- as.data.frame(t(Normal_data[, -1]))
data_ <- data.matrix(data)

##### assign column names with gene.symbols and row with sample ID
colnames(data_) <- Normal_data$Gene.Symbol
rownames(data_) <- names(Normal_data)[-1]
View(data_)

##### diagnosis the good sample genes with verbose 3
gsg <- goodSamplesGenes(data_, verbose = 3)
gsg$allOK

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(data_)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(data_)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  data_ <- data_[gsg$goodSamples, gsg$goodGenes]
}


##### clustering Normal genes
sampleTree = flashClust(dist(data_), method = "average");
sizeGrWindow(12,9)
pdf(file = "Normal_sample_Clustering.pdf", width = 12, height = 9)
par(cex = 0.9)
par(mar = c(0,6,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

dev.off()
data.Expression.Normal=data_
save(data.Expression.Normal, file = "Normal.RData")

##### removing outliers
Normal.omit.outlier <- Normal_data[ , -which(names(Normal_data) %in% c("GSM5319934", 'GSM5319899', 'GSM5319916'))]
write.csv(Normal.omit.outlier, file = "Normal.omit.outlier.csv", quote = F, row.names = F)

################################################Tumor tissue
##### Recall data and save Normal data samples' name into an array
Tumor_data <- read.csv("Tumor_samples.csv")
array.names <- names(Tumor_data)

##### transpose the Tumor samples data frame and remove the Gene.symbols row
data1 <- as.data.frame(t(Tumor_data[, -1]))
data1_ <- data.matrix(data1)

##### assign column names with gene.symbols and row with sample ID
colnames(data1_) = Tumor_data$Gene.Symbol
rownames(data1_) = names(Tumor_data)[-1]
View(data1_)


##### diagnosis the good sample genes with verbose 3
gsg = goodSamplesGenes(data1_, verbose = 3)
gsg$allOK

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(data1_)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(data1_)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  data1_ = data1_[gsg$goodSamples, gsg$goodGenes]
}


##### clustering Tumor genes

sampleTree = flashClust(dist(data1_), method = "average");
sizeGrWindow(12,9)
pdf(file = "Tumor.sample.Clustering.pdf", width = 12, height = 9)
par(cex = 0.9)
par(mar = c(0,6,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

dev.off()
data.Expression.Tumor=data1_
save(data.Expression.Tumor, file = "Tumor.RData")

##### removing outliers
Tumor.omit.outlier <- Tumor_data[ , -which(names(Tumor_data) %in% c("GSM5319853", "GSM5319839","GSM5319886"))]
write.csv(Tumor.omit.outlier, file = "Tumor.omit.outlier.csv", quote = F, row.names = F) 


