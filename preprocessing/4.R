##### divide Deg data.frame into Normal and Tumor sample data.frame

library(dplyr)
setwd('C:/Users/Atefeh/Desktop/data')

data <- read.csv('Deg.csv')

Tumor_samples <- select(data, 1:58)
View(Tumor_samples)

write.csv(Tumor_samples, 'Tumor_samples.csv', quote = F, row.names = F)

Normal_samples <- select(data, 1, 59:115)
View(Normal_samples)

write.csv(Normal_samples, 'Normal_samples.csv', quote = F, row.names = F)
