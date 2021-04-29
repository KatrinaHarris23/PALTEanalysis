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
                                     
                                     
#trying to see if days 44-90 taken out from B1 and B2 help differentiate between environment at all. 
#need to get rid of the rows that contain time points 44-90 for the B1 and B2 populations.
cast_genetics_p3 <-cast_genetics[-c(3,4,5,6,9,10,11,12),]
#now need to get rid of columns that have no variation these are all of the columns containing mutations that were only in the population time points that I just deleted. PCA analysis will fail if I don't do this.
cast_genetics_p3_2 <- cast_genetics_p3[,-c(6,11,12,13,16,19,20,25,27,29,31,42,44,47,49,55,56,57,58,59,61,63,64,70,76,79,80,88,91,96,97,98,101,103,107,109,111,114,123,124,126,133,139,142,143,144,148,155,156,159,160,163,164,166,168,170,177,179,182,185)]
cast_genetics_p3_3 <- cast_genetics_p3_2[,-c(126,130,134,136,139,142,143,146,149,150,152,154,155,156,157,158,161,170,175,176,177,180,184,187,188,192,195,197,198,199,200,201,205,207,208,215,216,218,220,221,225,226,227,231,234,236,240,243,244,245,246,247,252,255,256,260,267,268,269,271,272,278,280,283,293,294,296,311,312,317,321,322,323,324,326,328,329,331,332,335,336,338,352,359,364,378,379,381,383,386,388,391,393,398,401,402,404,420,422,424,429,434,435,436,437,438,439,443,444,448,453,457,458,459,462,463,464,466,467,470,478,480,484,491,492,496,497,501,503,509,511,514,515,517)]


genetics_pca_p3 <- prcomp(cast_genetics_p3_3[,1:385], center = T, scale = T)
which(apply(cast_genetics_p3_3, 2, var)==0) #this line tests to make sure that there are no columns with a variance of 0. 


theme_set(theme_bw())
autoplot(genetics_pca_p3, data= cast_genetics_p3_3, colour = "Sample", shape = "Time", size = 8)+ scale_color_manual(values=c("blue4","royalblue3","lightskyblue1","red4","orangered1","sandybrown","Black")) + 
  theme(legend.text=element_text(size=30)) + 
  theme(axis.text = element_text(size = 28), axis.title=element_text(size=30, face="bold")) + 
  theme(legend.title = element_text(size = 30))                                     
                                     
                                     



