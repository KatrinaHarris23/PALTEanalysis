#this analysis is done off of the P2 population resurected from freezer stocks at days 17, 44, 66, and 90. They are isolated both in planktonic conditions and after 1 additional day of biofilm selection. 

#samples were sequenced on the MiGS 190828 run. they were trimmed so that at each time point there were the same number of reads for the two samples. 

#→ to >
#← to <
#+ to plus
#Δ to delta
#remove spaces from everywhere except description
#remove , 
#remve % symbol from population data 


#Average coverage from Breseq: 
#P2_17_B -
#P2_44_B -
#P2_66_B -
#P2_90_B -
#P2_17_P -
#P2_44_P -
#P2_66_P -
#P2_90_P -

library(ggplot2)
theme_set(theme_bw())

#read in the ancestral mutations
ancestor_snps <- read.csv("/Users/katrina/Desktop/working/KBH5_WT_Breseq_Output.csv", header=TRUE, stringsAsFactors = F)
#head(ancestor_snps) #the SeqID column is the one I want 
nrow(ancestor_snps)#434

#read in the mutations from the sequencing run
P2_comparison <- read.csv("/Users/katrina/Desktop/scripts/P2_standard_Breseq_Output.csv", header=TRUE, stringsAsFactors = F)
#head(P2_comparison) #SeqID column is the one that I want
nrow(P2_comparison) #5191

#create a data frame with only those mutaitons found in the seqeunced samples but not found in the ancestral strain (at the position level)
P2_comparison_filtered <- P2_comparison[!(P2_comparison$SeqID %in% ancestor_snps$SeqID),]
nrow(P2_comparison_filtered)#3272


#write the file to a folder to have for later use
write.csv(P2_comparison_filtered, "/Users/katrina/Desktop/P2_comparison/P2_comparison_filtered.csv")


#Going to filter the popuations according to mutational frequencies through time to see if I can get rid of a lot of false positives.

#inport the data
filtering <- read.csv("/Users/katrina/Desktop/P2_comparison/P2_comparison_filtered.csv", stringsAsFactors = F)

#combine all of the extra data into one column
filtering$info <- paste(filtering$SeqID, filtering$Position, filtering$Annotation, filtering$Gene, filtering$Description, sep = ":::")


#now split into individual populations, based on the name of the population
B_17 <- filtering[(filtering$Sample == "17B"),]
B_44 <- filtering[(filtering$Sample == "44B"),]
B_66 <- filtering[(filtering$Sample == "66B"),]
B_90 <- filtering[(filtering$Sample == "P2_90B_standard"),]

P_17 <- filtering[(filtering$Sample == "P2_17P_standard"),]
P_44 <- filtering[(filtering$Sample == "P2_44P_standard"),]
P_66 <- filtering[(filtering$Sample == "66P"),]
P_90 <- filtering[(filtering$Sample == "90P"),]

#create one data frame for the planktonic data and one for the biofilm data

planktonic <- rbind(P_17,P_44,P_66,P_90)
biofilm <- rbind(B_17,B_44,B_66,B_90)

#View(planktonic)

#needto put these in a time sequence, so I need to melt and cast each of the data frames.
#first planktonic
library(reshape2)
m_planktonic <- melt(planktonic, id=c("X","Sample","SeqID","Position","Annotation","Gene","Description","info"), measure.vars = c("Mutation"))

#now biofilm
m_biofilm <- melt(biofilm, id=c("X","Sample","SeqID","Position","Annotation","Gene","Description","info"), measure.vars = c("Mutation"))

#and I need to cast the data frames so that each mutation is looked at through time.
#first planktonic
cast_planktonic <- t(dcast(m_planktonic, Sample~info, mean, value.var = "value", fill = 0)) #cast the data in the order that I want it
cast_planktonic <- as.data.frame(cast_planktonic, header=T) #make sure it is in a data frame type
colnames(cast_planktonic) <- as.character(unlist(cast_planktonic[1,])) #make sure the column names are what I want them to be
cast_planktonic$"0" <- 0.0 #add in a time 0 
planktonic_timeseries <- cast_planktonic[-1,] # the row names are the same as the first row, so I et rid of the first row
#reorder the columns
plank_col_order <- c("0","P2_17P_standard","P2_44P_standard","66P","90P") #the order I want the columns to be

library(data.table)
setcolorder(planktonic_timeseries, plank_col_order) #set the column order

#then biofilm
cast_biofilm <- t(dcast(m_biofilm, Sample~info, mean, value.var = "value", fill = 0)) #cast the data in the order that I want it
cast_biofilm <- as.data.frame(cast_biofilm, header=T) #make sure it is in a data frame type
colnames(cast_biofilm) <- as.character(unlist(cast_biofilm[1,])) #make sure the column names are what I want them to be
cast_biofilm$"0" <- 0.0 #add in a time 0 
biofilm_timeseries <- cast_biofilm[-1,] # the row names are the same as the first row, so I et rid of the first row
#reorder the columns
biof_col_order <- c("0","17B","44B","66B","P2_90B_standard") #the order I want the columns to be
setcolorder(biofilm_timeseries, biof_col_order) #set the column order

#convert all of the frequencies into type numeric, becuase they are currently factors
# first do planktonic 
planktonic_frequencies <- (apply(planktonic_timeseries[,1:5],2,function(x) as.numeric(as.character(x))))
rownames(planktonic_frequencies) <- rownames(planktonic_timeseries)

#then do biofilm
biofilm_frequencies <- (apply(biofilm_timeseries[,1:5],2,function(x) as.numeric(as.character(x))))
rownames(biofilm_frequencies) <- rownames(biofilm_timeseries)


#I don't know what is wrong with the above that it won't allow me to just use it with rowsums without coercing into a list, but since I can't figure it out I am just going to save it to a file and then load that file. 
write.csv(planktonic_frequencies, "/Users/katrina/Desktop/P2_comparison/planktonic_frequencies.csv")
write.csv(biofilm_frequencies, "/Users/katrina/Desktop/P2_comparison/biofilm_frequencies.csv")

plank_freqs <- read.csv("/Users/katrina/Desktop/P2_comparison/planktonic_frequencies.csv", stringsAsFactors = F)
biof_freqs <- read.csv("/Users/katrina/Desktop/P2_comparison/biofilm_frequencies.csv", stringsAsFactors = F)

#change the time 0 to be values of 0.0
plank_freqs$X0 <- 0.0
biof_freqs$X0 <- 0.0

#first looking at the planktonic file
#a column with the number of time points that are above 0 frequncy
plank_freqs$count <- rowSums(plank_freqs!=0.0)
#a column with the total frequency of each mutation
plank_freqs$sums <- rowSums(plank_freqs[,2:6])

#now I am going to implement the exact same filters that I used in the metagenomic sequencing.

#filter 1!!! select rows with at least 10% total frequncy
filter1_p <- subset(plank_freqs, plank_freqs$sums >= 10)
#filter 2!!! select only rows that appear in more than 1 day
filter2_p <- subset(filter1_p, filter1_p$count > 1)
#filter 3!!! remove if at 95% frequency at firt sequenced time point.
filter3_p <- subset(filter2_p, filter2_p$P2_17P_standard < 95)
#filter 4!!! remove if the mutation doesn't change in frequency by at least 10% 
filter3_p$min <- min(filter3_p[,2:5]) #find the smallest frequency
filter3_p$max <- max(filter3_p[,2:5]) #find the largest frequency

filter4_p <- subset(filter3_p, filter3_p$max-filter3_p$min >= 10)

#just seeing where we lose mutations
nrow(plank_freqs)#693
nrow(filter1_p)#553
nrow(filter2_p)#553
nrow(filter3_p)#537
nrow(filter4_p)#537


#now filter the biofilm ones
#View(biof_freqs)
biof_freqs$X0 <- 0.0
#a column with the number of time points that are above 0 frequncy
biof_freqs$count <- rowSums(biof_freqs[,3:6]!=0.0)
#a column with the total frequency of each mutation
biof_freqs$sums <- rowSums(biof_freqs[,3:6])

#now I am going to implement the exact same filters that I used in the metagenomic sequencing.
#filter 1!!! select rows with at least 10% total frequncy
filter1_b <- subset(biof_freqs, biof_freqs$sums >= 10)
#filter 2!!! select only rows that appear in more than 1 day

filter2_b <- subset(filter1_b, filter1_b$count > 1)
#filter 3!!! remove if at 95% frequency at firt sequenced time point

filter3_b <- subset(filter2_b, filter2_b$X17B < 95)
#filter 4!!! remove if the mutation doesn't change in frequency by at least 10% 
filter3_b$min <- min(filter3_b[,3:6]) #find the smallest frequency
filter3_b$max <- max(filter3_b[,3:6]) #find the largest frequency

filter4_b <- subset(filter3_b, filter3_b$max-filter3_b$min >= 10)

#just seeing where we lose mutations
nrow(biof_freqs)#778
nrow(filter1_b)#637
nrow(filter2_b)#322 I have no idea why there are so many mutations seen at only one time point.
nrow(filter3_b)#306
nrow(filter4_b)#306

#now save these filtered files: 
write.csv(filter4_p, "/Users/katrina/Desktop/P2_comparison/plankotnic_filtered.csv")
write.csv(filter4_b, "/Users/katrina/Desktop/P2_comparison/biofilm_filtered.csv")

#add a biofilm/planktonic designator column
filter4_b$sample <- "biofilm"
filter4_p$sample <- "planktonic"

head(filter4_b)
#now split into the different time points
filtered_b17 <- filter4_b[,c(1,3,11)]
filtered_b44 <- filter4_b[,c(1,4,11)]
filtered_b66 <- filter4_b[,c(1,5,11)]
filtered_b90 <- filter4_b[,c(1,6,11)]

filtered_p17 <- filter4_p[,c(1,3,11)]
filtered_p44 <- filter4_p[,c(1,4,11)]
filtered_p66 <- filter4_p[,c(1,5,11)]
filtered_p90 <- filter4_p[,c(1,6,11)]

colnames(filtered_b17) <- c("info","X17","Sample")
colnames(filtered_b44) <- c("info","X44","Sample")
colnames(filtered_b66) <- c("info","X66","Sample")
colnames(filtered_b90) <- c("info","X90","Sample")
colnames(filtered_p17) <- c("info","X17","Sample")
colnames(filtered_p44) <- c("info","X44","Sample")
colnames(filtered_p66) <- c("info","X66","Sample")
colnames(filtered_p90) <- c("info","X90","Sample")

filtered_day17 <- rbind(filtered_b17, filtered_p17)
filtered_day44 <- rbind(filtered_b44, filtered_p44)
filtered_day66 <- rbind(filtered_b66, filtered_p66)
filtered_day90 <- rbind(filtered_b90, filtered_p90)


#now make the tables compare the frequencies of the same mutations between the two samples.
#day 17
m_filtered_17 <- melt(filtered_day17, id=c("info","Sample"), measure.vars="X17")
#cast the data
cast_filtered_17 <- t(dcast(m_filtered_17, Sample~info, mean, value.var="value", fill = 0))
cast_filtered_17 <- as.data.frame(cast_filtered_17, header=T)
colnames(cast_filtered_17)<- c("Biofilm","Planktonic")
cast_filtered_17 <- cast_filtered_17[-1,]

split_17_filtered <- colsplit(rownames(cast_filtered_17), ":::", names = c("Position","Mutation", "Annotation", "Gene", "Description")) 

cast_filtered_17$Position <- split_17_filtered$Position
cast_filtered_17$Mutation <- split_17_filtered$Mutation
cast_filtered_17$Annotation <- split_17_filtered$Annotation
cast_filtered_17$Gene <- split_17_filtered$Gene
cast_filtered_17$Description <- split_17_filtered$Description

#and write the file
write.csv(cast_filtered_17, "/Users/katrina/Desktop/P2_comparison/day17_filtered.csv")

#day 44
m_filtered_44 <- melt(filtered_day44, id=c("info","Sample"), measure.vars="X44")
#cast the data
cast_filtered_44 <- t(dcast(m_filtered_44, Sample~info, mean, value.var="value", fill = 0))
cast_filtered_44 <- as.data.frame(cast_filtered_44, header=T)
colnames(cast_filtered_44) <- c("Biofilm", "Planktonic")
cast_filtered_44 <- cast_filtered_44[-1,]

split_44_filtered <- colsplit(rownames(cast_filtered_44), ":::", names = c("Position","Mutation", "Annotation", "Gene", "Description")) 

cast_filtered_44$Position <- split_44_filtered$Position
cast_filtered_44$Mutation <- split_44_filtered$Mutation
cast_filtered_44$Annotation <- split_44_filtered$Annotation
cast_filtered_44$Gene <- split_44_filtered$Gene
cast_filtered_44$Description <- split_44_filtered$Description

#and write the file
write.csv(cast_filtered_44, "/Users/katrina/Desktop/P2_comparison/day44_filtered.csv")

#day 66
m_filtered_66 <- melt(filtered_day66, id=c("info","Sample"), measure.vars="X66")
#cast the data
cast_filtered_66 <- t(dcast(m_filtered_66, Sample~info, mean, value.var="value", fill = 0))
cast_filtered_66 <- as.data.frame(cast_filtered_66, header=T)
colnames(cast_filtered_66) <- c("Biofilm", "Planktonic")
cast_filtered_66 <- cast_filtered_66[-1,]

split_66_filtered <- colsplit(rownames(cast_filtered_66), ":::", names = c("Position","Mutation", "Annotation", "Gene", "Description")) 

cast_filtered_66$Position <- split_66_filtered$Position
cast_filtered_66$Mutation <- split_66_filtered$Mutation
cast_filtered_66$Annotation <- split_66_filtered$Annotation
cast_filtered_66$Gene <- split_66_filtered$Gene
cast_filtered_66$Description <- split_66_filtered$Description

#and write the file
write.csv(cast_filtered_66, "/Users/katrina/Desktop/P2_comparison/day66_filtered.csv")


#day 90
m_filtered_90 <- melt(filtered_day90, id=c("info","Sample"), measure.vars="X90")
#cast the data
cast_filtered_90 <- t(dcast(m_filtered_90, Sample~info, mean, value.var="value", fill = 0))
cast_filtered_90 <- as.data.frame(cast_filtered_90, header=T)
colnames(cast_filtered_90) <- c("Biofilm", "Planktonic")
cast_filtered_90 <- cast_filtered_90[-1,]

split_90_filtered <- colsplit(rownames(cast_filtered_90), ":::", names = c("Position","Mutation", "Annotation", "Gene", "Description")) 
#head(split_90_filtered)
cast_filtered_90$Position <- split_90_filtered$Position
cast_filtered_90$Mutation <- split_90_filtered$Mutation
cast_filtered_90$Annotation <- split_90_filtered$Annotation
cast_filtered_90$Gene <- split_90_filtered$Gene
cast_filtered_90$Description <- split_90_filtered$Description

#and write the file
write.csv(cast_filtered_90, "/Users/katrina/Desktop/P2_comparison/day90_filtered.csv")








#############
#plot according to enrichment

#####
#this is the same for all of them It just depends on which file you import for if it is the unfiltered or the filtered data.

#these import the unfiltered data
#d17 <- read.csv("/Users/katrina/Desktop/B1_comparison/B1_17_comparison.csv", stringsAsFactors = F) 
#d44 <- read.csv("/Users/katrina/Desktop/B1_comparison/B1_44_comparison.csv", stringsAsFactors = F) 
#d66 <- read.csv("/Users/katrina/Desktop/B1_comparison/B1_66_comparison.csv", stringsAsFactors = F) 
#d90 <- read.csv("/Users/katrina/Desktop/B1_comparison/B1_90_comparison.csv", stringsAsFactors = F) 

#these import the filtered data
d17 <- read.csv("/Users/katrina/Desktop/P2_comparison/day17_filtered.csv", stringsAsFactors = F)
d44 <- read.csv("/Users/katrina/Desktop/P2_comparison/day44_filtered.csv", stringsAsFactors = F)
d66 <- read.csv("/Users/katrina/Desktop/P2_comparison/day66_filtered.csv", stringsAsFactors = F)
d90 <- read.csv("/Users/katrina/Desktop/P2_comparison/day90_filtered.csv", stringsAsFactors = F)

d17 <- d17[!(d17$Biofilm == 0),]
d17 <- d17[!(d17$Planktonic == 0),]
d44 <- d44[!(d44$Biofilm == 0),]
d44 <- d44[!(d44$Planktonic == 0),]
d66 <- d66[!(d66$Biofilm == 0),]
d66 <- d66[!(d66$Planktonic == 0),]
d90 <- d90[!(d90$Biofilm == 0),]
d90 <- d90[!(d90$Planktonic == 0),]




#create a color scheme, where they are colored based on their relative frequencies. These color schemes can be changed later.
d17$color <- 1
for (i in 1:nrow(d17)) {
  if (d17[i,2] == 0) 
    d17[i,9] <- "Biofilm enriched"
  else if (d17[i,2] < d17[i,3])
    d17[i,9] <- "Planktonic enriched"
  else if (d17[i,2] == d17[i,3])
    d17[i,9] <- "equal"
  else
    d17[i,9] <- "other"
}
d17$color <- as.factor(d17$color)
#create a color scheme, where they are colored based on their relative frequencies. These color schemes can be changed later.
d44$color <- 1
for (i in 1:nrow(d44)) {
  if (d44[i,2] == 0) 
    d44[i,9] <- "Biofilm enriched"
  else if (d44[i,2] < d44[i,3])
    d44[i,9] <- "Planktonic enriched"
  else if (d44[i,2] == d44[i,3])
    d44[i,9] <- "equal"
  else
    d44[i,9] <- "other"
}
d44$color <- as.factor(d44$color)
#create a color scheme, where they are colored based on their relative frequencies. These color schemes can be changed later.
d66$color <- 1
for (i in 1:nrow(d66)) {
  if (d66[i,2] == 0) 
    d66[i,9] <- "Biofilm enriched"
  else if (d66[i,2] < d66[i,3])
    d66[i,9] <- "Planktonic enriched"
  else if (d66[i,2] == d66[i,3])
    d66[i,9] <- "equal"
  else
    d66[i,9] <- "other"
}
d66$color <- as.factor(d66$color)
#create a color scheme, where they are colored based on their relative frequencies. These color schemes can be changed later.
d90$color <- ""
for (i in 1:nrow(d90)) {
  if (d90[i,2] == 0) 
    d90[i,9] <- "Biofilm enriched"
  else if (d90[i,2] < d90[i,3])
    d90[i,9] <- "Planktonic enriched"
  else if (d90[i,2] == d90[i,3])
    d90[i,9] <- "equal"
  else
    d90[i,9] <- "other"
}
d90$color <- as.factor(d90$color)

#colnames(d90)

#print the files with the color column so that I can manually change them in excel if I want to. 

write.csv(d17,"/Users/katrina/Desktop/200117/d17_colors2.csv")
write.csv(d44,"/Users/katrina/Desktop/200117/d44_colors2.csv")
write.csv(d66,"/Users/katrina/Desktop/200117/d66_colors2.csv")
write.csv(d90,"/Users/katrina/Desktop/200117/d90_colors2.csv")



View(d17)
#############



#calculating correlation statistics





######
#calculate the spearman's coefficient to see if there is any statistical significance to the trends seen 
#first need to make sure that there is a value in both of the frequency columns 

#day 17 statistics
d17_spearman <- d17[,2:3]
d17_spearman$total <- d17_spearman[,1] + d17_spearman[,2]
d17_spearman <- d17_spearman[(d17_spearman$total != 0.0),]

cor.test(~ Planktonic + Biofilm,
         data = d17_spearman, 
         method = "pearson", 
         continuity = F, 
         conf.level = 0.95)
#Pearson for Day 17
#Pearson's product-moment correlation: t = 23.416, df = 108, p-value < 2.2e-16
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
  0.8768542 0.9403362
sample estimates:
  cor 
0.9140275



#day 44 statistics
d44_spearman <- d44[,2:3]
d44_spearman$total <- d44_spearman[,1] + d44_spearman[,2]
d44_spearman <- d44_spearman[(d44_spearman$total != 0.0),]

cor.test(~ Planktonic + Biofilm,
         data = d44_spearman, 
         method = "pearson", 
         continuity = F, 
         conf.level = 0.95)
#Day 44
#Pearson's product-moment correlation: t = 24.944, df = 127, p-value < 2.2e-16
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
  0.8765344 0.9366260
sample estimates:
  cor 
0.9113131



#day 66 statistics
d66_spearman <- d66[,2:3]
d66_spearman$total <- d66_spearman[,1] + d66_spearman[,2]
d66_spearman <- d66_spearman[(d66_spearman$total != 0.0),]

cor.test(~ Planktonic + Biofilm,
         data = d66_spearman, 
         method = "pearson", 
         continuity = F, 
         conf.level = 0.95)
#Day 66 
#Pearson's product-moment correlation: t = 25.703, df = 151, p-value < 2.2e-16
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
  0.8677429 0.9280199
sample estimates:
  cor 
0.902196



#day 90 statistics
d90_spearman <- d90[,2:3]
d90_spearman$total <- d90_spearman[,1] + d90_spearman[,2]
d90_spearman <- d90_spearman[(d90_spearman$total != 0.0),]

cor.test(~ Planktonic + Biofilm,
         data = d90_spearman, 
         method = "pearson", 
         continuity = F, 
         conf.level = 0.95)
#data:  Planktonic and Biofilm t = 15.193, df = 123, p-value < 2.2e-16
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
  0.7365510 0.8611536
sample estimates:
  cor 
0.8076889 



#plot each time point in ggplot with coloring according to the categories we put them in above.
#need to set the color so that I don't get the grey background in the ggplots.

ggplot(d17, aes(x =Biofilm, y=Planktonic, size = 3)) + 
  geom_point() + 
  xlab("Biofilm selection") +
  ylab("Planktonic selection") +
  ggtitle("A. Day 17") +
  theme(aspect.ratio=1)+
  geom_abline(intercept=0)+
  ylim(0,100)+
  xlim(0,100)+
  theme(legend.text=element_text(size=18))

ggplot(d44, aes(x =Biofilm, y=Planktonic, size = 3)) + 
  geom_point() + 
  xlab("Biofilm selection") +
  ylab("Planktonic selection") +
  ggtitle("A. Day 44") +
  theme(aspect.ratio=1)+
  geom_abline(intercept=0)+
  ylim(0,100)+
  xlim(0,100)+
  theme(legend.text=element_text(size=18))

ggplot(d66, aes(x =Biofilm, y=Planktonic, size = 3)) + 
  geom_point() + 
  xlab("Biofilm selection") +
  ylab("Planktonic selection") +
  ggtitle("C. Day 66") +
  theme(aspect.ratio=1)+
  geom_abline(intercept=0)+
  ylim(0,100)+
  xlim(0,100)+
  theme(legend.text=element_text(size=18))

ggplot(d90, aes(x =Biofilm, y=Planktonic, size = 3)) + 
  geom_point() + 
  xlab("Biofilm selection") +
  ylab("Planktonic selection") +
  ggtitle("D. Day 90") +
  theme(aspect.ratio=1)+
  geom_abline(intercept=0)+
  ylim(0,100)+
  xlim(0,100)+
  theme(legend.text=element_text(size=18))

