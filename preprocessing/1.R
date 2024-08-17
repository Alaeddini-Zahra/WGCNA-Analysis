library(dplyr)
library(plyr)
result <- read.delim("C:/Users/Atefeh/Desktop/ovarian/GSE115635_table.txt")

##### select a subset based on the adj.P.Val
new_table <- subset(result, logFC < 1)

##### removing missing values in Gene.symbol column

df <- new_table[!(new_table$Gene.symbol == "" ) ,] 

##### removing duplicates based on the Gene.symbol column

t1 <- df %>% group_by(Gene.symbol) %>% filter (! duplicated(Gene.symbol))

##### save filtered data as main_filter

write.table(t1, "C:/Users/Atefeh/Desktop/data/main_filter.txt", quote = F, sep = '\t', row.names = F)

