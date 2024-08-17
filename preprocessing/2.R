setwd("C:/Users/Atefeh/Desktop/data/")

##### Recall the DEG genes
df <- read.csv("main_filter.csv")
class(df)
View(df)

##### Recall the driver genes
df1 <- read.csv("cancer driver genes.csv")
class(df1)
View(df1)

##### matching Gene.symbol column in df and df1
match_genes <- match(df1$Gene.symbol, df$Gene.Symbol)
match_genes <- df[match_genes,]
View(match_genes)

match_genes <- na.omit(match_genes)
View(match_genes)

##### Save selected driver-DEG Genes file
write.table(match_genes, "Selected Driver-Deg Genes.txt", sep="\t", quote = F, row.names = F)

