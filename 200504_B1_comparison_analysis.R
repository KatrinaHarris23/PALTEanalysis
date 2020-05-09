#this analysis is done off of the B1 population resurected from freezer stocks at days 17, 44, 66, and 90. They are isolated both in planktonic conditions and after 1 additional day of biofilm selection. 

#samples were sequenced on runs on 160611 and 160628. this analysis is using reads from both runs. 

#Coverage: 
#17_B - 190.1
#44_B - 149.3
#66_B - 165.8
#90_B - 130.4
#17_P - 272.0
#44_P - 144.4
#66_P - 188.5
#90_P - 161.3
#####
library(ggplot2)
theme_set(theme_bw())

#read in the ancestral mutations
ancestor_snps <- read.csv("/Users/katrina/Desktop/working/KBH5_WT_Breseq_Output.csv", header=TRUE, stringsAsFactors = F)
#head(ancestor_snps) #the SeqID column is the one I want 
nrow(ancestor_snps)#434

#read in the mutations from the sequencing run
B1_comparison <- read.csv("/Users/katrina/Desktop/200117/Breseq_Output.csv", header=TRUE, stringsAsFactors = F)
#head(B1_comparison) #Seq.ID column is the one that I want
nrow(B1_comparison) #5191

#create a data frame with only those mutaitons found in the seqeunced samples but not found in the ancestral strain (at the position level)
B1_comparison_filtered <- B1_comparison[!(B1_comparison$SeqID %in% ancestor_snps$SeqID),]
nrow(B1_comparison_filtered)#2903

#write the file to a folder to have for later use
write.csv(B1_comparison_filtered, "/Users/katrina/Desktop/200117/B1_comparison_filtered2.csv")


#Going to filter the popuations according to mutational frequencies through time to see if I can get rid of a lot of false positives.

#inport the data
filtering <- read.csv("/Users/katrina/Desktop/200117/B1_comparison_filtered2.csv", stringsAsFactors = F)


#combine all of the extra data into one column
filtering$info <- paste(filtering$SeqID, filtering$Position, filtering$Annotation, filtering$Gene, filtering$Description, sep = ":::")

#now split into individual populations, based on the name of the population
B_17 <- filtering[(filtering$Sample == "B1_17B_standard"),]
B_44 <- filtering[(filtering$Sample == "B1_44B_standard"),]
B_66 <- filtering[(filtering$Sample == "B1_66B_standard"),]
B_90 <- filtering[(filtering$Sample == "B1_90B_standardized"),]

P_17 <- filtering[(filtering$Sample == "B1_17P_standard"),]
P_44 <- filtering[(filtering$Sample == "B1_44P_standard"),]
P_66 <- filtering[(filtering$Sample == "B1_66P_standard"),]
P_90 <- filtering[(filtering$Sample == "B1_90P_standardized"),]

#create one data frame for the planktonic data and one for the biofilm data

planktonic <- rbind(P_17,P_44,P_66,P_90)
biofilm <- rbind(B_17,B_44,B_66,B_90)

#View(planktonic)

#needto put these in a time sequence, so I need to melt and cast each of the data frames.
#first planktonic
library(data.table)
m_planktonic <- melt(planktonic, id=c("X","Sample","Evidence","SeqID","Position","Annotation","Gene","Description","info"), measure.vars = c("Mutation"))

#now biofilm
m_biofilm <- melt(biofilm, id=c("X","Sample","Evidence","SeqID","Position","Annotation","Gene","Description","info"), measure.vars = c("Mutation"))

#and I need to cast the data frames so that each mutation is looked at through time.
#first planktonic
cast_planktonic <- t(dcast(m_planktonic, Sample~info, mean, value.var = "value", fill = 0)) #cast the data in the order that I want it
cast_planktonic <- as.data.frame(cast_planktonic, header=T) #make sure it is in a data frame type
colnames(cast_planktonic) <- as.character(unlist(cast_planktonic[1,])) #make sure the column names are what I want them to be
cast_planktonic$"0" <- 0.0 #add in a time 0 
planktonic_timeseries <- cast_planktonic[-1,] # the row names are the same as the first row, so I et rid of the first row
#reorder the columns
plank_col_order <- c("0","B1_17P_standard","B1_44P_standard","B1_66P_standard", "B1_90P_standardized")#,"90P") #the order I want the columns to be
setcolorder(planktonic_timeseries, plank_col_order) #set the column order

#then biofilm
cast_biofilm <- t(dcast(m_biofilm, Sample~info, mean, value.var = "value", fill = 0)) #cast the data in the order that I want it
cast_biofilm <- as.data.frame(cast_biofilm, header=T) #make sure it is in a data frame type
colnames(cast_biofilm) <- as.character(unlist(cast_biofilm[1,])) #make sure the column names are what I want them to be
cast_biofilm$"0" <- 0.0 #add in a time 0 
biofilm_timeseries <- cast_biofilm[-1,] # the row names are the same as the first row, so I et rid of the first row
#reorder the columns
biof_col_order <- c("0","B1_17B_standard","B1_44B_standard","B1_66B_standard","B1_90B_standardized") #the order I want the columns to be
setcolorder(biofilm_timeseries, biof_col_order) #set the column order

#convert all of the frequencies into type numeric, becuase they are currently factors
# first do planktonic 
planktonic_frequencies <- (apply(planktonic_timeseries[,1:5],2,function(x) as.numeric(as.character(x))))
rownames(planktonic_frequencies) <- rownames(planktonic_timeseries)

#then do biofilm
biofilm_frequencies <- (apply(biofilm_timeseries[,1:5],2,function(x) as.numeric(as.character(x))))
rownames(biofilm_frequencies) <- rownames(biofilm_timeseries)


#I don't know what is wrong with the above that it won't allow me to just use it with rowsums without coercing into a list, but since I can't figure it out I am just going to save it to a file and then load that file. 
write.csv(planktonic_frequencies, "/Users/katrina/Desktop/200117/planktonic_frequencies2.csv")
write.csv(biofilm_frequencies, "/Users/katrina/Desktop/200117/biofilm_frequencies2.csv")





plank_freqs <- read.csv("/Users/katrina/Desktop/200117/planktonic_frequencies2.csv", stringsAsFactors = F)
biof_freqs <- read.csv("/Users/katrina/Desktop/200117/biofilm_frequencies2.csv", stringsAsFactors = F)

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
filter2_p <- subset(filter1_p, filter1_p$count > 2)
#filter 3!!! remove if at 95% frequency at firt sequenced time point.
filter3_p <- subset(filter2_p, filter2_p$B1_17P_standard < 95)
#filter 4!!! remove if the mutation doesn't change in frequency by at least 10% 
filter3_p$min <- min(filter3_p[,2:6]) #find the smallest frequency
filter3_p$max <- max(filter3_p[,2:6]) #find the largest frequency

filter4_p <- subset(filter3_p, filter3_p$max-filter3_p$min >= 10)

#just seeing where we lose mutations
nrow(plank_freqs)#773
nrow(filter1_p)#573
nrow(filter2_p)#216
nrow(filter3_p)#199
nrow(filter4_p)#199


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
#filter 3!!! remove if at 95% frequency at firt sequenced time point.

filter3_b <- subset(filter2_b, filter2_b$B1_17B_standard < 95)
#filter 4!!! remove if the mutation doesn't change in frequency by at least 10% 
filter3_b$min <- min(filter3_b[,3:6]) #find the smallest frequency
filter3_b$max <- max(filter3_b[,3:6]) #find the largest frequency

filter4_b <- subset(filter3_b, filter3_b$max-filter3_b$min >= 10)

#just seeing where we lose mutations
nrow(biof_freqs)#827
nrow(filter1_b)#605
nrow(filter2_b)#212
nrow(filter3_b)#196
nrow(filter4_b)#196

#now save these filtered files: 
write.csv(filter4_p, "/Users/katrina/Desktop/200117/plankotnic_filtered2.csv")
write.csv(filter4_b, "/Users/katrina/Desktop/200117/biofilm_filtered2.csv")

#add a biofilm/planktonic designator column
filter4_b$sample <- "biofilm"
filter4_p$sample <- "planktonic"


#now split into the different time points
filtered_b17 <- filter4_b[,c(1,3,11)]
filtered_b44 <- filter4_b[,c(1,4,11)]
filtered_b66 <- filter4_b[,c(1,5,11)]
filtered_b90 <- filter4_b[,c(1,6,11)]

filtered_p17 <- filter4_p[,c(1,3,10)]
filtered_p44 <- filter4_p[,c(1,4,10)]
filtered_p66 <- filter4_p[,c(1,5,10)]
filtered_p90 <- filter4_p[,c(1,6,11)] #where did it go? 

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
write.csv(cast_filtered_17, "/Users/katrina/Desktop/200117/day17_filtered2.csv")

#day 44
m_filtered_44 <- melt(filtered_day44, id=c("info","Sample"), measure.vars="X44")
#cast the data
cast_filtered_44 <- t(dcast(m_filtered_44, Sample~info, mean, value.var="value", fill = 0))
cast_filtered_44 <- as.data.frame(cast_filtered_44, header=T)
colnames(cast_filtered_44) <- c("Biofilm", "Planktonic")
cast_filtered_44 <- cast_filtered_44[-1,]

library(reshape)
split_44_filtered <- colsplit(rownames(cast_filtered_44), ":::", names = c("Position","Mutation", "Annotation", "Gene", "Description")) 

cast_filtered_44$Position <- split_44_filtered$Position
cast_filtered_44$Mutation <- split_44_filtered$Mutation
cast_filtered_44$Annotation <- split_44_filtered$Annotation
cast_filtered_44$Gene <- split_44_filtered$Gene
cast_filtered_44$Description <- split_44_filtered$Description

#and write the file
write.csv(cast_filtered_44, "/Users/katrina/Desktop/200117/day44_filtered2.csv")

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
write.csv(cast_filtered_66, "/Users/katrina/Desktop/200117/day66_filtered2.csv")


#day 90
m_filtered_90 <- melt(filtered_day90, id=c("info","Sample"), measure.vars="X90")
#cast the data
cast_filtered_90 <- t(dcast(m_filtered_90, Sample~info, mean, value.var="value", fill = 0))
cast_filtered_90 <- as.data.frame(cast_filtered_90, header=T)
colnames(cast_filtered_90) <- c("Biofilm", "Planktonic")
cast_filtered_90 <- cast_filtered_90[-1,]

split_90_filtered <- colsplit(rownames(cast_filtered_90), ":::", names = c("Position","Mutation", "Annotation", "Gene", "Description")) 

cast_filtered_90$Position <- split_90_filtered$Position
cast_filtered_90$Mutation <- split_90_filtered$Mutation
cast_filtered_90$Annotation <- split_90_filtered$Annotation
cast_filtered_90$Gene <- split_90_filtered$Gene
cast_filtered_90$Description <- split_90_filtered$Description

#and write the file
write.csv(cast_filtered_90, "/Users/katrina/Desktop/200117/day90_filtered2.csv")




#############
#this is the same for all of them It just depends on which file you import for if it is the unfiltered or the filtered data.

#these import the unfiltered data
#d17 <- read.csv("/Users/katrina/Desktop/B1_comparison/B1_17_comparison.csv", stringsAsFactors = F) 
#d44 <- read.csv("/Users/katrina/Desktop/B1_comparison/B1_44_comparison.csv", stringsAsFactors = F) 
#d66 <- read.csv("/Users/katrina/Desktop/B1_comparison/B1_66_comparison.csv", stringsAsFactors = F) 
#d90 <- read.csv("/Users/katrina/Desktop/B1_comparison/B1_90_comparison.csv", stringsAsFactors = F) 

#these import the filtered data
d17 <- read.csv("/Users/katrina/Desktop/200117/day17_filtered2.csv", stringsAsFactors = F)
d44 <- read.csv("/Users/katrina/Desktop/200117/day44_filtered2.csv", stringsAsFactors = F)
d66 <- read.csv("/Users/katrina/Desktop/200117/day66_filtered2.csv", stringsAsFactors = F)
d90 <- read.csv("/Users/katrina/Desktop/200117/day90_filtered2.csv", stringsAsFactors = F)

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
#Pearson's product-moment correlation: t = 9.3883, df = 36, p-value = 3.263e-11
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
#  0.7157341 0.9156521
#sample estimates:
#  cor 
#0.8426192 



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
#Pearson's product-moment correlation: t = 8.6438, df = 165, p-value = 4.546e-15
#50 percent confidence interval: 0.4440735 0.6546474
#sample estimates: cor - 0.558287



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
#Pearson's product-moment correlation: t = 21.539, df = 86, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
#  0.8779308 0.9459549
#sample estimates:
#  cor 
#0.9184857 



#day 90 statistics
d90_spearman <- d90[,2:3]
d90_spearman$total <- d90_spearman[,1] + d90_spearman[,2]
d90_spearman <- d90_spearman[(d90_spearman$total != 0.0),]

cor.test(~ Planktonic + Biofilm,
         data = d90_spearman, 
         method = "pearson", 
         continuity = F, 
         conf.level = 0.95)
#data:  Planktonic and Biofilm t = 24.57, df = 80, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
#  0.9077790 0.9607665
#sample estimates:
#  cor 
#0.9396752



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

#obs1 is the random binomial distribution of a first observation of the sequencing data. this can be used to see which planktonic frequencies




#this is a binomial distribution. this is what the sequencing data on two runs should look like.

#find the binomial distribution confidence intervals
library(binom)
x <- (binom.confint(obs1, 1000, conf.int = 0.95, methods= "exact")$lower)
y <- (binom.confint(obs2, 1000, conf.int = 0.95, methods= "exact")$lower)

plot(obs1, obs2)





