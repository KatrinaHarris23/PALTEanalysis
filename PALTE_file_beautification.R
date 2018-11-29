# This script is me just taking the output of filtered mutations for each PALTE population (output from the filtering scripts I made), and changing the formatting so that they are much more human readable. These should at some point be added to the end of each one of those populations, but for the moment this is how it is working. 

#input is a csv file

#output is a csv file -- more human readable


library("vegan")
library("plyr")
library("RColorBrewer")
library("ggplot2")
library("data.table")
library("dplyr")
library("reshape2")
#library("xlsx")
library("scales")

theme_set(theme_bw())

setwd("/Users/katrina/Desktop/working")


#make the B1 population pretty -> splitting columns and renaming columns to make the files more human readable - also makes it so that I can use the defaults
B1 <- read.csv("/Users/katrina/Desktop/working/B1/B1_allfilters.csv")
View(B1)
splitB1_1 <- colsplit(B1$X, ";;", names = c("desc_gene_annot","Position","Mutation"))
splitB1_2 <- colsplit(splitB1_1$desc_gene_annot, "::", names = c("Description", "Gene","Annotation"))

B1 <- cbind(splitB1_2,splitB1_1[,2:3],B1)
View(B1)
ncol(B1)
B1 <- B1[,1:13]
B1 <- B1[-6]
colnames(B1)
B1names <- c("Description", "Gene", "Annotation", "Position", "Mutation", "0", "17", "25", "44", "66", "75", "90")
colnames(B1) <- B1names

View(B1) #I want description and mutation columns
colnames(B1)
write.csv(B1,file="/Users/katrina/Desktop/working/B1_pretty.csv")



######################### B2

B2 <- read.csv("/Users/katrina/Desktop/working/B2/B2_allfilters.csv")
View(B2)
splitB2_1 <- colsplit(B2$X, ";;", names = c("desc_gene_annot","Position","Mutation"))
splitB2_2 <- colsplit(splitB2_1$desc_gene_annot, "::", names = c("Description", "Gene","Annotation"))

B2 <- cbind(splitB2_2,splitB2_1[,2:3],B2)
View(B2)
ncol(B2)
B2 <- B2[,1:13]
B2 <- B2[-6]
colnames(B2)
B2names <- c("Description", "Gene", "Annotation", "Position", "Mutation", "0", "17", "25", "44", "66", "75", "90")
colnames(B2) <- B2names

View(B2) #I want description and mutation columns
write.csv(B2,file="/Users/katrina/Desktop/working/B2_pretty.csv")


######################## B3 population

B3 <- read.csv("/Users/katrina/Desktop/working/B3/B3_allfilters.csv")
View(B3)
splitB3_1 <- colsplit(B3$X, ";;", names = c("desc_gene_annot","Position","Mutation"))
splitB3_2 <- colsplit(splitB3_1$desc_gene_annot, "::", names = c("Description", "Gene","Annotation"))

B3 <- cbind(splitB3_2,splitB3_1[,2:3],B3)
ncol(B3)
B3 <- B3[,1:13]
B3 <- B3[-6]
colnames(B3)
B3names <- c("Description", "Gene", "Annotation", "Position", "Mutation", "0", "17", "25", "44", "66", "75", "90")
colnames(B3) <- B3names

View(B3) #I want description and mutation columns
write.csv(B3,file="/Users/katrina/Desktop/working/B3_pretty.csv")





############################ P1
P1 <- read.csv("/Users/katrina/Desktop/working/P1/P1_allfilters.csv")
View(P1)
splitP1_1 <- colsplit(P1$X, ";;", names = c("desc_gene_annot","Position","Mutation"))
splitP1_2 <- colsplit(splitP1_1$desc_gene_annot, "::", names = c("Description", "Gene","Annotation"))

P1 <- cbind(splitP1_2,splitP1_1[,2:3],P1)
View(P1)
ncol(P1)
P1 <- P1[,1:11]
P1 <- P1[-6]
colnames(P1)
P1names <- c("Description", "Gene", "Annotation", "Position", "Mutation", "0", "17", "44", "75", "90")
colnames(P1) <- P1names
nrow(P1)
View(P1) #I want description and mutation columns
write.csv(P1,file="/Users/katrina/Desktop/working/P1_pretty.csv")




############################ P2
P2 <- read.csv("/Users/katrina/Desktop/working/P2/P2_allfilters.csv")
splitP2_1 <- colsplit(P2$X, ";;", names = c("desc_gene_annot","Position","Mutation"))
splitP2_2 <- colsplit(splitP2_1$desc_gene_annot, "::", names = c("Description", "Gene","Annotation"))

P2 <- cbind(splitP2_2,splitP2_1[,2:3],P2)
P2 <- P2[,1:11]
P2 <- P2[-6]
colnames(P2)
P2names <- c("Description", "Gene", "Annotation", "Position", "Mutation", "0", "17", "44", "75", "90")
colnames(P2) <- P2names

View(P2) #I want description and mutation columns
write.csv(P2,file="/Users/katrina/Desktop/working/P2_pretty.csv")




############################ P3
P3 <- read.csv("/Users/katrina/Desktop/working/P3/P3_allfilters.csv")
splitP3_1 <- colsplit(P3$X, ";;", names = c("desc_gene_annot","Position","Mutation"))
splitP3_2 <- colsplit(splitP3_1$desc_gene_annot, "::", names = c("Description", "Gene","Annotation"))

P3 <- cbind(splitP3_2,splitP3_1[,2:3],P3)
P3 <- P3[,1:11]
P3 <- P3[-6]
colnames(P3)
P3names <- c("Description", "Gene", "Annotation", "Position", "Mutation", "0", "17", "44", "75", "90")
colnames(P3) <- P3names

View(P3) #I want description and mutation columns
write.csv(P3, file="/Users/katrina/Desktop/working/P3_pretty.csv")



