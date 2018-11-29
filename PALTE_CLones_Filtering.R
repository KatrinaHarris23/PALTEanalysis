#analysis of isolates from the PALTE that Kenny published in his Jbact paper in 2016. I got his fastq files from UNH and I resequenced the 7 lettered isolates (A,D,H,M,O,P,V) so there will be 3 files to import (ancester, UNH sequencing, and pitt sequencing)
library("vegan")
library("plyr")
library("RColorBrewer")
library("ggplot2")
library("data.table")
library("dplyr")
library("reshape2")
library("xlsx")
library("scales")

theme_set(theme_bw())
setwd("/Users/katrina/Desktop/working")


##get rid of spaces in all columns except description.
#convert arrows ← and → o < and > get rid of all commas, remove Â, ‑,–

#import all three files and see how long they are
ancestor <- read.csv("KBH5_WT_Breseq_Output.csv", header = TRUE)
View(ancestor) #going to want the SeqID column 
nrow(ancestor)#434

#7 clones sequenced at pitt
PittClones <- read.csv("Clones_Breseq_Output.csv", header = TRUE)
View(PittClones) #going to want the position column
nrow(PittClones)#2374

#the 20 clones kenny published
UNHClones <- read.csv("UNHClones_Breseq_Output.csv", header = TRUE)
View(UNHClones) #going to want the position column
nrow(UNHClones)#6262


#take out the mutations from each clonal data set that are found in the ancestor
PittClonesfiltered <- PittClones[!(PittClones$Position %in% ancestor$SeqID),]
View(PittClonesfiltered)
nrow(PittClonesfiltered)#813
write.csv(PittClonesfiltered, file = "PittSequencedClonesFiltered.csv")

UNHfiltered <- UNHClones[!(UNHClones$Position %in% ancestor$SeqID),]
View(UNHfiltered)
nrow(UNHfiltered)#1944
write.csv(UNHfiltered, file = "UNHSequencedClonesFiltered.csv")
