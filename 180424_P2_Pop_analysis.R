#this is the P2 population of the LTE that I am trying to automate the analysis of

#The input for this script is a csv file of a breseq output. Can use the breseq_cat script that nate/Chris Dietrick made to get an excel file. You do have to do a little bit of manipulation on this file in excel first, described in comments below.

#this script will take the list of mutations and subtract mutations:1. found in the ancestor, 2. that don't reach a cumulative frequency of 10%, 3. that only appear at one time point, 4. that are fixed from the first measured time point, 5. that are almost fixed from the first time point, and 6. that do not change in frequency by at least 10% over the course of the experiment.

#the output of this script is an excel sheet that is formatted to go directly into Katya's matlab scripts to define mutational cohorts, infer ancestry, and make muller plots. there are many places that data frames can be printed for other purposes.

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

#first thing is downloading time course breseq data from beagle
#/home/kah231/scripts/BreseqCat.py -d /home/kah231/PA/analysis/LTEpitt_seq/P2
#and the ancester
#/home/kah231/scripts/SingleBreseqCat.py -f /home/kah231/PA/analysis/LTEpitt_seq/KBH5_WT/output/index.html
#
#
##get rid of spaces in all columns except description.
#convert arrows ← and → o < and > get rid of all commas, remove Â, ‑,–
#need to make sure days are just days not "P217" etc.
##saved the snp tab as a csv format


#FIle locations
#/Users/katrina/Desktop/working
#called:
#  P2_Breseq_Output.csv
#  KBH5_WT_Breseq_Output.csv

setwd("/Users/katrina/Desktop/working")

#ancestral/background SNPs
ancestor_snps <- read.csv("KBH5_WT_Breseq_Output.csv", header=TRUE)
head(ancestor_snps)
View(ancestor_snps) #want SeqID column because that is where the positions are

#P2 population SNPs
P2_snps <- read.csv("P2_Breseq_Output.csv",header=TRUE)
View(P2_snps) #again SeqID is what we want to match
nrow(P2_snps) #3941 total

#get the SNPs that are found in both files according to the position numbers that are
#actually stored in the SeqID column
#filter1 = no mutations found in the ancestral clone
P2_filter1 <- P2_snps[ !(P2_snps$SeqID %in% ancestor_snps$SeqID), ]
#see how many rows are in the data frame now
nrow(P2_filter1) #2981
#create a data frame of the shared mutations (the ones taken out with filter 1)
P2_filter1_refmatch <- P2_snps[ (P2_snps$SeqID %in% ancestor_snps$SeqID), ]
nrow(P2_filter1_refmatch)#960
write.csv(P2_filter1_refmatch, file = "P2_ancestral.csv")


#will create a df with taking out based on the gene name. I don't think this will be what I want, I think it will be too strict but I am doing it anyway
#I don't like using this because some genes really can get more than one mutation, especially if it is a big gene
P2_genelevelfilter <- P2_snps[ !(P2_snps$Gene %in% ancestor_snps$Gene), ]
nrow(P2_genelevelfilter)#2757
#now see how many it took out
P2_genelevelfilter_refmatch <- P2_snps[ (P2_snps$Gene %in% ancestor_snps$Gene), ]
nrow(P2_genelevelfilter_refmatch)#1184

#this is done to the original data set through filter1
#remove '%' symbol
P2_filter1$Mutation <- gsub( "%", "", as.character(P2_filter1$Mutation), n)
P2_filter1$Mutation <- as.numeric(as.character(P2_filter1$Mutation))

#combine annotation::gene::description in one column
P2_filter1$desc_gene_annot <-  paste(P2_filter1$Annotation,P2_filter1$Gene,P2_filter1$Description, sep="::")
P2_filter1$desc_gene_annot <- as.factor(P2_filter1$desc_gene_annot)
#combine desc_gene_annot and position (SeqID) in one column
P2_filter1$details <- paste(P2_filter1$desc_gene_annot, P2_filter1$SeqID, sep=";;")
P2_filter1$details <- as.factor(P2_filter1$details)
#combine details with the actual mutation
P2_filter1$info <- paste(P2_filter1$details, P2_filter1$Position, sep=";;")
P2_filter1$info <- as.factor(P2_filter1$info)
View(P2_filter1)

#melt data frame for casting
m_P2_filter1 <- melt(P2_filter1, id=c("Sample","Evidence","SeqID","Position","Annotation","Gene","Description","desc_gene_annot", "details", "info"),measure.vars = c("Mutation"))
head(m_P2_filter1)
View(m_P2_filter1)

#cast data frame - organizing with each mutation as the rows and the frequency of that mutation on a given day as the columns
P2_cast1 <- t(dcast(m_P2_filter1,Sample~info,mean, value.var = "value",fill=0))
P2_cast1 <- as.data.frame(P2_cast1,header=TRUE)
colnames(P2_cast1) <- as.character(unlist(P2_cast1[1,]))
colnames(P2_cast1)
View(P2_cast1)
P2_cast1$"0" <- 0.0
View(P2_cast1)
P2_cast1 <- P2_cast1[-1,]
View(P2_cast1)

#need to reorder the columns in ascending order
P2_column_order <- c("0","17", "44","66","90")
colnames(P2_cast1)
setcolorder(P2_cast1, P2_column_order)
View(P2_cast1)

nrow(P2_cast1)#1632

#transpose the matrix
t_P2 <- as.data.frame(t(P2_cast1))
View(t_P2)
#figure out what class the frequency values are in the matrix - they need to be numeric
class(t_P2[2,2]) #factor
ncol(t_P2)#1632 <- for sanity this check should match up with what you found previously for a nrow count

#convert frequency values to numeric class - start as "character"
t_P2[,2:1632] <-(apply(t_P2[,2:1632], 2, function(x) as.numeric(as.character(x))))
P2 <- transpose(t_P2[,2:1632])
colnames(P2) <- rownames(t_P2)
rownames(P2) <- colnames(t_P2[,2:1632])
View(P2)

class(P2[2,2]) #yay, it's numeric now!!

#adds a count number on each row (mutation) that tells how many columns (days) have a frequency above 0
P2$count <- rowSums(P2!=0.0)
#sums up the rows to tell you the total % that you get - flaw is that it also adds the count, so need to subtract that value
P2$Sums <- rowSums(P2)-P2[,6]
nrow(P2)#1631
View(P2) #these last two are mostly just sanity checks

#filter 2!!!!! select only rows with greater than 10% total frequency
P2_filter2 <- (subset(P2, P2$Sums >= 10)) #greater than 10%
nrow(P2_filter2)#898

P2_filter2_out <- (subset(P2, P2$Sums < 10)) #put all of the filtered out mutations in one place

#filter 3!! select only rows that appear in more than 1 day
P2_filter3 <- (subset(P2_filter2, P2_filter2$count > 1))
nrow(P2_filter3)#678
P2_filter3_out <- (subset(P2_filter2, P2_filter2$count <= 1))
View(P2_filter3)

#filter 4 -- remove all mutations at 100% across all measured time points
#4 time points so 400 value --> problem with this is mutations can start at 100 and dip slightly
P2_filter4 <- (subset(P2_filter3, P2_filter3$Sums < 400))
nrow(P2_filter4) #667
P2_filter4_out <- (subset(P2_filter3, P2_filter3$Sums >= 400))

#filter out if the first time points that are 95% or above
P2_filter5 <- (subset(P2_filter4, P2_filter4$"17" < 95))
nrow(P2_filter5) #666
P2_filter5_out <- (subset(P2_filter4, P2_filter4$"17" >= 95))

#filter out if the HIGHEST frequency isn't 10, not if the combined total frequency doesn't get to 10

#filter out if the change in frequency from the first time point to the last time point does not change by at least 10%, having the additive value be above 10 isn't stringent enough.
############need to change this because it does not do what I thought it did
#P2_filter6 <- (subset(P2_filter5, (P2_filter5$"17"+P2_filter5$"90")/2 >= 10))
#nrow(P2_filter6)#248
#filter6_removed <- (subset(P2_filter5, (P2_filter5$"17"+P2_filter5$"90")/2 < 10))
#View(filter6_removed)
#I really should be checking what mutations are being taken out with each filter...
#View(P2_filter6)
###############
#what about subtracting 90 from 17 and thenn if the absolute value isn't above 10 filter it out
P2_filter62 <- (subset(P2_filter5, abs(P2_filter5$"17"-P2_filter5$"90") >= 10))
nrow(P2_filter62) #92
P2_filter62_out <- (subset(P2_filter5, abs(P2_filter5$"17"-P2_filter5$"90") < 10))
View(P2_filter62)

not_real_mutations <- rbind(P2_filter2_out, P2_filter3_out, P2_filter4_out, P2_filter5_out, P2_filter62_out)
write.csv(not_real_mutations, file = "P2_Filtered_out.csv")


ncol(P2_filter62)#7
P2_filter7 <- P2_filter62[,-c(6,7)] #remove columns with count and sums
View(P2_filter7)

#split the row names into columns again
P2_split <- P2_filter7
P2_split$info <- rownames(P2_filter7)
View(P2_split)

P2_split4 = transform(P2_split, info =colsplit(P2_split$info,';;', names = c('desc_gene_annot','position', 'Mutation')))
#View(P2_split4)
P2_split4 <- as.data.frame(P2_split4)
colnames(P2_split4)#6
info <- (P2_split4$info)
head(info)
P2_split4 <- P2_split4[,c(1:5)]
View(P2_split4)

#P2_split4$desc_gene_annot <- rownames(P2_split4)
#P2_split4$desc_gene_annot <- gsub(";;.*","",P2_split4$desc_gene_annot)
#does_this_work <- merge(P2_split4, info, by="desc_gene_annot")
#this does work but I like my way better because I know what it is doing and I still have the rownames
#View(does_this_work )

View(P2_split4)
View(info)
P2_split5 <- P2_split4
P2_split5$Position <- info$position
P2_split5$Mutation <- info$Mutation
View(P2_split5)
#P2_split6 <- P2_split5[,-6]
#View(P2_split6)

#rename columns after splitting
#colnames(P2_split6)
colnames(P2_split5) <- c("0","17", "44","66","90","Position", "Mutation")
View(P2_split5)

#write this to a file that I can find
write.csv(P2_split5,file="P2_allfilters.csv")

#oh plotting....
#have to melt it first...
#melt - with whatever i want to keep
#transform
View(P2_split5)
nrow(P2_split5)#92

t_P2_plot <- t(P2_split5)
nrow(t_P2_plot)#7
t_P2_plot_2 <- t_P2_plot[-c(6,7),]
m_P2_plot <- melt(t_P2_plot_2, id=c("desc_gene_annot"),Value.name = "frequency")
head(m_P2_plot)

#plot
#plot(NA, xlim=c(0,100), ylim=c(0,100))
#lines(m_P2_plot$X1, m_P2_plot$value)
#lines(c(0,17,44,66,90), t_final_P2_filter3[,1], add=TRUE)
#lines(c(0,17,44,66,90), t_final_P2_filter3[,2], add=TRUE)
#this is how I could add things one at a time

#but I want to know how to use ggplot

colnames(m_P2_plot) <- c("day", "mutation","value")
m_P2_plot$value <- as.numeric(as.character(m_P2_plot$value)) # you have to change from a factor to a number this way. You have to co to a character first always. if I had just gone to a number it would have given me the level of factor that the previous data point was.
View(m_P2_plot)

ggplot(m_P2_plot,aes(x=day,y=value,color=mutation)) +theme(text = element_text(size=20),legend.text=element_text(size=10),legend.position="none") +geom_point(size=4) +geom_line()


#this data set is going to be formatted to go into katya's matlab scripts. The data frame needs specific columns with particular data. they will have to be in a specific order and named correctly, but first I just need to make them holding the correct information.
View(P2_split5)
P2_Muller <- P2_split5[,c(1:5)]
P2_Muller$Population <- "P2"
P2_Muller$Population2 <- 1 #make sure this is a number
P2_Muller$Chromosome <- 1 #make sure this is a number
P2_Muller$Position <- P2_split5$Position
P2_Muller$Class <- "SNP"
P2_Muller$Mutation <- P2_split5$Mutation
P2_Muller$Gene <- ""
P2_Muller$AminoAcid <- ""
P2_Muller$Class2 <- ""
P2_Muller$Amino <- ""
P2_Muller$NearestDownstreamGene <- ""
P2_Muller$Distance <- ""
P2_Muller$Trajectory <- 1:nrow(P2_Muller)

colnames(P2_Muller)
#now put the columns in the correct order
Muller_col_order <- c("Population", "Population2","Trajectory","Chromosome","Position","Class","Mutation","Gene","AminoAcid","Class2","Amino","NearestDownstreamGene","Distance","0","17","44","66","90")
setcolorder(P2_Muller,Muller_col_order)

#now I need to name them what they are actually supposed to be named
colnames(P2_Muller)

colnames(P2_Muller) <-c("Population", "Population number","Trajectory","Chromosome","Position","Class","Mutation","Gene","Amino Acid","Class","Amino","Nearest Downstream Gene","Distance","0","17","44","66","90")

View(P2_Muller)
#need to remove the rownames
#rownames(P2_Muller) <- c()
#View(P2_Muller)
#decided not to do this because I can just print without including the row names. this will also just print out the row names, they will just be the numbers instead of the descriptions


#latest problem is that the frequencies need to be percentages and the column names for the frequencies need to be numbers.
#first solve the frequencies to percentages problem - should be able to do with scales package
#would like to keep 2 decimal points if possible

P2_Muller_try <- P2_Muller
P2_Muller_try$`90` <- as.numeric(as.character(P2_Muller_try$`90`))
P2_Muller_try$`66` <- as.numeric(as.character(P2_Muller_try$`66`))
P2_Muller_try$`44` <- as.numeric(as.character(P2_Muller_try$`44`))
P2_Muller_try$`17` <- as.numeric(as.character(P2_Muller_try$`17`))
P2_Muller_try$`0` <- as.numeric(as.character(P2_Muller_try$`0`))

#the following does work. I divide everything by 100 so that when in excel I can change to "percent" type and it will be the correct value. Remember to keep 1 decimal when changing the type in excel or it will round everything.
P2_Muller_try$`90` <- (P2_Muller_try$`90`/100)
P2_Muller_try$`66` <- (P2_Muller_try$`66`/100)
P2_Muller_try$`44` <- (P2_Muller_try$`44`/100)
P2_Muller_try$`17` <- (P2_Muller_try$`17`/100)
P2_Muller_try$`0` <- (P2_Muller_try$`0`/100)


View(P2_Muller_try)

#now to write this file so that I can use it as an input to matlab. The matlab file requires it to be a .xlsx file so I can just write to that type of file. need to make sure that I don't print out the row names or they will be the first column. I do NEED the column names though.
write.xlsx(P2_Muller_try, file="P2_Muller.xlsx", row.names = FALSE)
#this file needs to be loaded into matlab for Katya's scripts.






#see what the last filter took out to make sure the mutations are what I want them to be. This check needs to be done on all mutations that are taken out from step 1 on.
filter62_removed <- (subset(P2_filter5, abs(P2_filter5$"17"-P2_filter5$"90") < 10))
View(filter62_removed)#yay it removed what I wanted it to!!
nrow(filter62_removed)#587
#I need to look at the pileups for these mutations to se why they are coming up in my populations.
#there are way more of these "nonsense" mutations than there dhould be if it were just normal sequencing errors. so the question is where is the source of the variation. I will look at the pileup at the areas of these mutations in IGV to try to determine the answer to this.

split_removed <- filter62_removed
split_removed$info <- rownames(filter62_removed)
split_removed2 = transform(split_removed, info =colsplit(split_removed$info,';;', names = c('desc_gene_annot','position', 'Mutation')))
split_removed2 <- as.data.frame(split_removed2)
colnames(split_removed2)
View(split_removed2)
info2 <-(split_removed2$info)
View(info2)
ncol(split_removed2)
split_removed2 <- split_removed2[,-c(8)]
View(split_removed2)
colnames(info2)
split_removed2$position <- info2$position
split_removed2$Mutation <- info2$Mutation
View(split_removed2)
write.csv(split_removed2, file = "P2_removed.csv")


#######################################
#the following is a section of MatLab code (so won't work to run in R)

#problem I ran into is that there is apparently something wrong with the excel file... it won't read in
#manipulations done to excel file once exported
#change column names of frequencies to "number" type
#change the frequencies to "percentage" type
#
names = ["P2_muller.xlsx"]; %take in the excel file you got from R
sheets = ["Sheet1"]; %it is an excel workbook so you have to tell it which sheet its on
tNum = 5; %the number of timepoints
timepoints = [0,17,44,66,90]; %the names of the timepoints
xUnits = ["Days"]; %units of the timepoints
time_series_import_KMF %this reads in the data
get_genotypes_KMF_KBH %this will determine the genotypes --> need to learn how it is defining a genotype
%genotype_plots_KMF %this plots the genotypes. It will give 3 plots. the first is unclustered, there is an intermediate, and the third is the final fully clustered data set.


order_clusters_KMF % this will determine the order in which the genotypes showed up

%ordered_cluster_plots_KMF %visualize the final clusters with the ancestry incorporated --> should save this figure.

frequencies = squeeze(genneststotal(1, any(squeeze(genneststotal(1, :, :)), 2), 2:trajSize)); % define the frequencies variable
nests = squeeze(genneststotal(1, any(squeeze(genneststotal(1, :, :)), 2), nameSize:end)); % define the nests variable
nests = nests(:, any(nests, 1)); #still defninng nests
#muller_plots_KMF(frequencies, nests, timepoints) % problem with this is that it always outputs with the x axis saying "time in generations" I want to be able to change this to whatever I want to

csvwrite("/Users/katrina/Desktop/working/P2/P2timepoints.csv", timepoints)
csvwrite("/Users/katrina/Desktop/working/P2/P2frequencies.csv",frequencies)
csvwrite("/Users/katrina/Desktop/working/P2/P2nests.csv",nests)

