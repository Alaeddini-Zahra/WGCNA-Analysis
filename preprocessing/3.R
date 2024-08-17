setwd("C:/Users/Atefeh/Desktop/data")


series_matrix <- read.delim("GSE174570_series_matrix.txt", comment.char = "!")
View(series_matrix)
write.csv(series_matrix, 'GSE117290_series_matrix.csv', quote = F, row.names = F)

##### recall driver deg genes
DDG <- read.csv("Selected Driver-Deg Genes.csv")
View(DDG)

##### ID matching 
ID_matching <- match(DDG$ID, series_matrix$ID_REF)
common_ID <- series_matrix[ID_matching,]
View(common_ID)

##### append a column for genes name
common_ID[, 116] <- DDG$Gene.Symbol
View(common_ID)

##### replace the ID_column with new genes' column
common_ID[, 1] <- common_ID[, 116]
View(common_ID)

##### remove the appended column and assign new name to the first column
common_ID <- common_ID[, -116]
colnames(common_ID)[1] <- "Gene.Symbol"
View(common_ID)

##### save file as Deg
write.csv(common_ID, "Deg.csv" , quote = F, row.names = F)
