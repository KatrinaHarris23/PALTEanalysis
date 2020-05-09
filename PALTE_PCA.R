#going to analyze the PALTE phenotypic traits that Imeasured in a PCA plot. I will be using the biofilm production, maximum growth rate (Vmax), and the swimming motility for this analysis. 

#this is a file of the evolved data, using the average values for each population statistic.
phenotypes <- read.csv("/Users/katrina/Desktop/PALTE_final/PALTE_PCA2.csv", stringsAsFactors = F)
phenotypes <- phenotypes[1:28,1:4] #now 7 observations (rows) of 4 variables (columns)
#phenotypes_numeric <- as.data.frame(apply(phenotypes[,2:4],2,function(x) as.numeric(as.character(x))))


phenotypes.pca <- prcomp(phenotypes[,2:4], center = T, scale = T)
summary(phenotypes.pca)
str(phenotypes.pca)

library(ggfortify)

autoplot(phenotypes.pca, data = phenotypes, colour = "Population", shape = ("Population"), size = 8)  + scale_color_manual(values=c("blue4","royalblue3","lightskyblue1","red4","orangered1","tomato1","Black")) +scale_shape_manual(values=c(15,15,15,16,16,16,18)) + theme(legend.text=element_text(size=20)) + theme(axis.text = element_text(size = 18), axis.title=element_text(size=20,face="bold"))



+ axis.text=element_text(size=18)
