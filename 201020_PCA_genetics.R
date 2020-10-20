#The current goal is to take this PCA, which includes an analysis for 6 populations at 6 time points, and plot each population individually. 
library(data.table)
library(reshape2)
library(ggfortify)



all <- read.csv("/Users/katrina/Desktop/PALTE_final/draft_revision_work/PCA_data.csv", stringsAsFactors = F)

cast_genetics <- (dcast(all, Sample~Gene, mean, value.var ="Frequency", fill = 0))
cast_genetics <- as.data.frame(cast_genetics, header = T)

cast_genetics <- as.data.frame(apply(cast_genetics[,2:580], 2, function(x) as.numeric(as.character(x))))

#renaming the samples so that all days are the same population name
cast_genetics$Sample <- c("B1","B1","B1","B1","B1","B1","B2","B2","B2","B2","B2","B2","B3","B3","B3","B3","B3", "B3","P1","P1","P1","P1","P2","P2","P2","P2","P3","P3","P3","P3")
cast_genetics$Time <- c("17","25","44","66","75","90","17","25","44","66","75","90","17","25","44","66","75","90","17","44","66","90","17","44","66","90","17","44","66","90")

#write.csv(cast_genetics, "/Users/katrina/Desktop/PALTE_final/draft_revision_work/genetics_PCA.csv")

#import data for PCA
#cast_genetics <- read.csv("/Users/katrina/Desktop/PALTE_final/draft_revision_work/genetics_PCA.csv", stringsAsFactors = F)

#apply(cast_genetics[,2:581], function(x) as.numeric(as.character(x)))
#now do PCA analysis on this data set.
genetics_pca <- prcomp(cast_genetics[,2:579], center = T, scale = T)


heatmap(genetics_pca)
summary(genetics_pca)
str(genetics_pca)

library(ggfortify)

theme_set(theme_bw())
autoplot(genetics_pca, data= cast_genetics, colour = "Sample", shape = "Time", size = 8)+ scale_color_manual(values=c("blue4","royalblue3","lightskyblue1","red4","orangered1","sandybrown","Black")) + 
  theme(legend.text=element_text(size=30)) + 
  theme(axis.text = element_text(size = 28), axis.title=element_text(size=30, face="bold")) + 
  theme(legend.title = element_text(size = 30))



all <- read.csv("/Users/katrina/Desktop/PALTE_final/draft_revision_work/all_201004.csv", stringsAsFactors = F)
heatmap(all[,4:7])



