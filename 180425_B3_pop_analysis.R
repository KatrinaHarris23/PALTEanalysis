#this is the B3 population of the LTE that I am trying to automate the analysis of

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

setwd("/Users/katrina/Desktop/working/B3")

#first thing is downloading time course breseq data from beagle
#/home/kah231/scripts/BreseqCat.py -d /home/kah231/PA/analysis/LTEpitt_seq/B3
#and the ancester
#/home/kah231/scripts/SingleBreseqCat.py -f /home/kah231/PA/analysis/LTEpitt_seq/KBH5_WT/output/index.html
#
#
##get rid of spaces in all columns except description.
#convert arrows ← and → o < and > get rid of all commas, remove Â, ‑,–
#need to make sure days are just days not "B317" etc.
##saved the snp tab as a csv format


#FIle locations
#/Users/katrina/Desktop/working
#called:
#  B3_Breseq_Output.csv
#  KBH5_WT_Breseq_Output.csv

#ancestral/background SNPs
ancestor_snps <- read.csv("/Users/katrina/Desktop/working/KBH5_WT_Breseq_Output.csv", header=TRUE)
head(ancestor_snps)
View(ancestor_snps) #want SeqID column because that is where the positions are

#B3 population SNPs
B3_snps <- read.csv("B3_Breseq_Output.csv",header=TRUE)
View(B3_snps) #again SeqID is what we want to match
nrow(B3_snps) #6771

#get the SNPs that are found in both files according to the position numbers that are
#actually stored in the SeqID column
#filter1 = no mutations found in the ancestral clone
B3_filter1 <- B3_snps[ !(B3_snps$SeqID %in% ancestor_snps$SeqID), ]
#see how many rows are in the data frame now
nrow(B3_filter1) #5333
#create a data frame of the shared mutations (the ones taken out with filter 1)
B3_filter1_refmatch <- B3_snps[ (B3_snps$SeqID %in% ancestor_snps$SeqID), ]
nrow(B3_filter1_refmatch)#1438
write.csv(B3_filter1_refmatch, file = "B3_ancestral.csv")


#will create a df with taking out based on the gene name. I don't think this will be what I want, I think it will be too strict but I am doing it anyway
#I don't like using this because some genes really can get more than one mutation, especially if it is a big gene
B3_genelevelfilter <- B3_snps[ !(B3_snps$Gene %in% ancestor_snps$Gene), ]
nrow(B3_genelevelfilter)#4956
#now see how many it took out
B3_genelevelfilter_refmatch <- B3_snps[ (B3_snps$Gene %in% ancestor_snps$Gene), ]
nrow(B3_genelevelfilter_refmatch)#1815

#this is done to the original data set through filter1
#remove '%' symbol
B3_filter1$Mutation <- gsub( "%", "", as.character(B3_filter1$Mutation), n)
B3_filter1$Mutation <- as.numeric(as.character(B3_filter1$Mutation))

#combine annotation::gene::description in one column
B3_filter1$desc_gene_annot <-  paste(B3_filter1$Annotation,B3_filter1$Gene,B3_filter1$Description, sep="::")
B3_filter1$desc_gene_annot <- as.factor(B3_filter1$desc_gene_annot)
#combine desc_gene_annot and position (SeqID) in one column
B3_filter1$details <- paste(B3_filter1$desc_gene_annot, B3_filter1$SeqID, sep=";;")
B3_filter1$details <- as.factor(B3_filter1$details)
#combine details with the actual mutation
B3_filter1$info <- paste(B3_filter1$details, B3_filter1$Position, sep=";;")
B3_filter1$info <- as.factor(B3_filter1$info)
View(B3_filter1)

#melt data frame for casting
m_B3_filter1 <- melt(B3_filter1, id=c("Sample","Evidence","SeqID","Position","Annotation","Gene","Description","desc_gene_annot", "details", "info"),measure.vars = c("Mutation"))
head(m_B3_filter1)
View(m_B3_filter1)

#cast data frame - organizing with each mutation as the rows and the frequency of that mutation on a given day as the columns
B3_cast1 <- t(dcast(m_B3_filter1,Sample~info,mean, value.var = "value",fill=0))
B3_cast1 <- as.data.frame(B3_cast1,header=TRUE)
colnames(B3_cast1) <- as.character(unlist(B3_cast1[1,]))
colnames(B3_cast1)
View(B3_cast1)
B3_cast1$"0" <- 0.0
View(B3_cast1)
B3_cast1 <- B3_cast1[-1,]
View(B3_cast1)

#need to reorder the columns in ascending order
B3_column_order <- c("0","17", "25", "44","66","75", "90")
colnames(B3_cast1)
setcolorder(B3_cast1, B3_column_order)
View(B3_cast1)

nrow(B3_cast1)#2349

#transpose the matrix
t_B3 <- as.data.frame(t(B3_cast1))
View(t_B3)
#figure out what class the frequency values are in the matrix - they need to be numeric
class(t_B3[2,2]) #factor
ncol(t_B3)#2349 <- for sanity this check should match up with what you found previously for a nrow count

#convert frequency values to numeric class - start as "character"
t_B3[,2:2349] <-(apply(t_B3[,2:2349], 2, function(x) as.numeric(as.character(x))))
B3 <- transpose(t_B3[,2:2349])
colnames(B3) <- rownames(t_B3)
rownames(B3) <- colnames(t_B3[,2:2349])
View(B3)

class(B3[2,2]) #yay, it's numeric now!!

#adds a count number on each row (mutation) that tells how many columns (days) have a frequency above 0
B3$count <- rowSums(B3!=0.0)
#sums up the rows to tell you the total % that you get - flaw is that it also adds the count, so need to subtract that value
B3$Sums <- rowSums(B3)-B3[,8]
nrow(B3)#2348
View(B3) #these last two are mostly just sanity checks

#filter 2!!!!! select only rows with greater than 10% total frequency
B3_filter2 <- (subset(B3, B3$Sums >= 10)) #greater than 10%
nrow(B3_filter2)#1268

B3_filter2_out <- (subset(B3, B3$Sums < 10)) #put all of the filtered out mutations in one place

#filter 3!! select only rows that appear in more than 1 day
B3_filter3 <- (subset(B3_filter2, B3_filter2$count > 1))
nrow(B3_filter3)#1011
B3_filter3_out <- (subset(B3_filter2, B3_filter2$count <= 1))

#filter 4 -- remove all mutations at 100% across all measured time points
#4 time points so 400 value --> problem with this is mutations can start at 100 and dip slightly
B3_filter4 <- (subset(B3_filter3, B3_filter3$Sums < 600))
nrow(B3_filter4) #998
B3_filter4_out <- (subset(B3_filter3, B3_filter3$Sums >= 600))

#filter out if the first time points that are 95% or above
B3_filter5 <- (subset(B3_filter4, B3_filter4$"17" < 95))
nrow(B3_filter5) #997
B3_filter5_out <- (subset(B3_filter4, B3_filter4$"17" >= 95))

#filter out if the HIGHEST frequency isn't 10, not if the combined total frequency doesn't get to 10

#filter out if the change in frequency from the first time point to the last time point does not change by at least 10%, having the additive value be above 10 isn't stringent enough.
############need to change this because it does not do what I thought it did
#B3_filter6 <- (subset(B3_filter5, (B3_filter5$"17"+B3_filter5$"90")/2 >= 10))
#nrow(B3_filter6)#248
#filter6_removed <- (subset(B3_filter5, (B3_filter5$"17"+B3_filter5$"90")/2 < 10))
#View(filter6_removed)
#I really should be checking what mutations are being taken out with each filter...
#View(B3_filter6)
###############
#what about subtracting 90 from 17 and thenn if the absolute value isn't above 10 filter it out
B3_filter62 <- (subset(B3_filter5, abs(B3_filter5$"17"-B3_filter5$"90") >= 10))
nrow(B3_filter62) #141
B3_filter62_out <- (subset(B3_filter5, abs(B3_filter5$"17"-B3_filter5$"90") < 10))
View(B3_filter62)

not_real_mutations <- rbind(B3_filter2_out, B3_filter3_out, B3_filter4_out, B3_filter5_out, B3_filter62_out)
write.csv(not_real_mutations, file = "B3_Filtered_out.csv")


ncol(B3_filter62)#9
B3_filter7 <- B3_filter62[,-c(8,9)] #remove columns with count and sums
View(B3_filter7)

#split the row names into columns again
B3_split <- B3_filter7
B3_split$info <- rownames(B3_filter7)
View(B3_split)

B3_split4 = transform(B3_split, info =colsplit(B3_split$info,';;', names = c('desc_gene_annot','position', 'Mutation')))
#View(B3_split4)
B3_split4 <- as.data.frame(B3_split4)
colnames(B3_split4)#8
info <- (B3_split4$info)
head(info)
B3_split4 <- B3_split4[,c(1:7)]
View(B3_split4)

#B3_split4$desc_gene_annot <- rownames(B3_split4)
#B3_split4$desc_gene_annot <- gsub(";;.*","",B3_split4$desc_gene_annot)
#does_this_work <- merge(B3_split4, info, by="desc_gene_annot")
#this does work but I like my way better because I know what it is doing and I still have the rownames
#View(does_this_work )

View(B3_split4)
View(info)
B3_split5 <- B3_split4
B3_split5$Position <- info$position
B3_split5$Mutation <- info$Mutation
View(B3_split5)
#B3_split6 <- B3_split5[,-6]
#View(B3_split6)

#rename columns after splitting
#colnames(B3_split6)
colnames(B3_split5) <- c("0","17","25", "44","66","75","90","Position", "Mutation")
View(B3_split5)

#write this to a file that I can find
write.csv(B3_split5,file="B3_allfilters.csv")

#oh plotting....
#have to melt it first...
#melt - with whatever i want to keep
#transform
View(B3_split5)
nrow(B3_split5)#141

t_B3_plot <- t(B3_split5)
nrow(t_B3_plot)#9
t_B3_plot_2 <- t_B3_plot[-c(8,9),]
m_B3_plot <- melt(t_B3_plot_2, id=c("desc_gene_annot"),Value.name = "frequency")
head(m_B3_plot)

#plot
#plot(NA, xlim=c(0,100), ylim=c(0,100))
#lines(m_B3_plot$X1, m_B3_plot$value)
#lines(c(0,17,44,66,90), t_final_B3_filter3[,1], add=TRUE)
#lines(c(0,17,44,66,90), t_final_B3_filter3[,2], add=TRUE)
#this is how I could add things one at a time

#but I want to know how to use ggplot

colnames(m_B3_plot) <- c("day", "mutation","value")
m_B3_plot$value <- as.numeric(as.character(m_B3_plot$value)) # you have to change from a factor to a number this way. You have to co to a character first always. if I had just gone to a number it would have given me the level of factor that the previous data point was.
View(m_B3_plot)

ggplot(m_B3_plot,aes(x=day,y=value,color=mutation)) +theme(text = element_text(size=20),legend.text=element_text(size=10),legend.position="none") +geom_point(size=4) +geom_line()


#####now to make the file pretty 


B3 <- read.csv("/Users/katrina/Desktop/working/B3/B3_allfilters.csv")


splitB3_1 <- colsplit(B3$X, ";;", names = c("desc_gene_annot","Position","Mutation"))
splitB3_2 <- colsplit(splitB3_1$desc_gene_annot, "::", names = c("Description", "Gene","Annotation"))
View(splitB3_1)

nrow(splitB3_1)
nrow(splitB3_2)

B3 <- cbind(splitB3_2,splitB3_1[,2:3],B3)
ncol(B3)
B3 <- B3[,1:13]
B3 <- B3[-6]
colnames(B3)
B3names <- c("Description", "Gene", "Annotation", "Position", "Mutation", "0", "17", "25", "44", "66", "75", "90")
colnames(B3) <- B3names

#View(B1) #I want description and mutation columns
colnames(B3)
write.csv(B3,file="/Users/katrina/Desktop/working/B3/B3_pretty.csv")




#this data set is going to be formatted to go into katya's matlab scripts. The data frame needs specific columns with particular data. they will have to be in a specific order and named correctly, but first I just need to make them holding the correct information.
View(B3)

B3_Muller <- B3[,6:12]
B3_Muller$Population <- "B1"
B3_Muller$Population2 <- 1L #make sure this is a number
B3_Muller$Chromosome <- 1L #make sure this is a number
B3_Muller$Position <- B3$Position
B3_Muller$Class <- "SNP"
B3_Muller$Mutation <- B3$Mutation
B3_Muller$Gene <- B3$Gene
B3_Muller$AminoAcid <- B3$Description
B3_Muller$Class2 <- ""
B3_Muller$Amino <- ""
B3_Muller$NearestDownstreamGene <- ""
B3_Muller$Distance <- ""
B3_Muller$Trajectory <- 1:nrow(B3_Muller)

colnames(B3_Muller)
#now put the columns in the correct order
Muller_col_order <- c("Population", "Population2","Trajectory","Chromosome","Position","Class","Mutation","Gene","AminoAcid","Class2","Amino","NearestDownstreamGene","Distance","0","17","25", "44","66","75", "90")
setcolorder(B3_Muller,Muller_col_order)

#now I need to name them what they are actually supposed to be named
colnames(B3_Muller)

colnames(B3_Muller) <-c("Population", "Population number","Trajectory","Chromosome","Position","Class","Mutation","Gene","Amino Acid","Class","Amino","Nearest Downstream Gene","Distance","0","17","25","44","66","75", "90")

View(B3_Muller)
#need to remove the rownames
#rownames(B3_Muller) <- c()
#View(B3_Muller)
#decided not to do this because I can just print without including the row names. this will also just print out the row names, they will just be the numbers instead of the descriptions


#latest problem is that the frequencies need to be percentages and the column names for the frequencies need to be numbers.
#first solve the frequencies to percentages problem - should be able to do with scales package
#would like to keep 2 decimal points if possible

B3_Muller_try <- B3_Muller
B3_Muller_try$`90` <- as.numeric(as.character(B3_Muller_try$`90`))
B3_Muller_try$`75` <- as.numeric(as.character(B3_Muller_try$`75`))
B3_Muller_try$`66` <- as.numeric(as.character(B3_Muller_try$`66`))
B3_Muller_try$`44` <- as.numeric(as.character(B3_Muller_try$`44`))
B3_Muller_try$`25` <- as.numeric(as.character(B3_Muller_try$`25`))
B3_Muller_try$`17` <- as.numeric(as.character(B3_Muller_try$`17`))
B3_Muller_try$`0` <- as.numeric(as.character(B3_Muller_try$`0`))

#the following does work. I divide everything by 100 so that when in excel I can change to "percent" type and it will be the correct value. Remember to keep 1 decimal when changing the type in excel or it will round everything.
B3_Muller_try$`90` <- (B3_Muller_try$`90`/100)
B3_Muller_try$`75` <- (B3_Muller_try$`75`/100)
B3_Muller_try$`66` <- (B3_Muller_try$`66`/100)
B3_Muller_try$`44` <- (B3_Muller_try$`44`/100)
B3_Muller_try$`25` <- (B3_Muller_try$`25`/100)
B3_Muller_try$`17` <- (B3_Muller_try$`17`/100)
B3_Muller_try$`0` <- (B3_Muller_try$`0`/100)


View(B3_Muller_try)

#now to write this file so that I can use it as an input to matlab. The matlab file requires it to be a .xlsx file so I can just write to that type of file. need to make sure that I don't print out the row names or they will be the first column. I do NEED the column names though.
write.csv(B3_Muller_try, file="B3_Muller.csv", row.names = FALSE)
#this file needs to be loaded into matlab for Katya's scripts.



#now to analyze the mutations... 

Mutations_analysis <- function(Mutation_Data, AminoAcid="description", Bases = "Mutation") {
  #Mutation_Data is the input CSV file that is a breseq output with rows containing different mutations and columns containing the various characteristics for those mutations.
  #AminoAcid is the column name that holds the breseq information on amino acid changes. this column will look like "*342C(TGA>TGC)" or say coding, pseudogene, etc. The default that will be looked for is "Description". This is case sensitive!
  #Bases is the column name that holds the breseq information for the nucleotide mutations. These are things like A>C, or T>G in the breseq output. This is case sensitive!!!
  
  ##############
  #first i am going to deal with the nucleotide level - the mutation column. This uses the Mutation_Data that you put in and grabs the information in the column that you specified under Bases. It looks for the > symbol, because that is how breseq separates the two bases, and splits the data. It then creates a new data set called Nucleotides that has 2 columns containing the original base (from the ancestor) and the mutated base.
  Nucleotides <- colsplit(Mutation_Data[,Bases], ">", names = c("original", "mutant"))
  #View(Nucleotides)
  #I want to calculate the total number of mutations present in the sample.
  allmutations <- nrow(Nucleotides)
  #I want to determine the number of mutations that are not just substituting one base from another. These are indels, because this is how breseq represents them.
  indel <- sum(grepl("[^ACGT]", Nucleotides$original))
  
  
  #I need to find all of the different combinations for base replacements. To do this I am going to find the index for each base, and then I will look at that index in the second column and see what the new base is.
  C <- grep("C", Nucleotides$original) #find placeswhere the original was C
  CT <- sum(ifelse(Nucleotides$mutant[C]=="T",1,0)) # find when there was a C to T transition
  CA <- sum(ifelse(Nucleotides$mutant[C]=="A",1,0)) # find when there was a C to A transition
  CG <- sum(ifelse(Nucleotides$mutant[C]=="G",1,0)) # find when there was a C to G transition
  
  Ts <- grep("T", Nucleotides$original) #find when the original was a T
  TC <- sum(ifelse(Nucleotides$mutant[Ts]=="C",1,0)) #find when there were T to C transitions
  TG <- sum(ifelse(Nucleotides$mutant[Ts]=="G",1,0)) #find when there were T to G transitions
  TA <- sum(ifelse(Nucleotides$mutant[Ts]=="A",1,0)) #find when there were T to A transitions
  
  G <- grep("G", Nucleotides$original) #find placeswhere the original was G
  GA <- sum(ifelse(Nucleotides$mutant[G]=="A",1,0)) # find when there was a G to A transition
  GT <- sum(ifelse(Nucleotides$mutant[G]=="T",1,0)) # find when there was a G to T transition
  GC <- sum(ifelse(Nucleotides$mutant[G]=="C",1,0)) # find when there was a G to C transition
  
  A <- grep("A", Nucleotides$original) #find placeswhere the original was A
  AG <- sum(ifelse(Nucleotides$mutant[A]=="G",1,0)) # find when there was a A to G transition
  AC <- sum(ifelse(Nucleotides$mutant[A]=="C",1,0)) # find when there was a A to C transition
  AT <- sum(ifelse(Nucleotides$mutant[A]=="T",1,0)) # find when there was a A to T transition
  
  # Now that I have the numbers of all of the possible base changes, I can look for the
  transitions <- sum(c(CT,TC,GA,AG)) #there are 4 options for transitions. C>T, T>C, G>A, A>G. this adds up all of those changes
  transversions <- AT+AC+GC+GT+CA+CG+TG+TA # need to do to check that the sums of the transition categories actually do add up to the number of transitions that there should be (assuming transitions and indel numbers are correct) when I turn this into a function I need to stop if transversions != trans -- should be fine but just an extra error checking step.
  
  ###############
  
  ### now at the Amino acid level
  #have to get the amino acid column that the user specifies out of the input data
  Protein <- colsplit(Mutation_Data[,AminoAcid], "\\(", names = c("AA", "DNA"))
  
  #there are a few options that you can get for this one and I can't just split the column. I need to look for all of them in the strings FOr this I need to use regular expressions. I will have to use gregexpr which returns a list of positions and then I will have to find the length of that to determine the number of them. The regular expressios that I will use for each option are as follows.
  
  # reg expressions to use "coding", "intergenic", "pseudogene",
  #"[A-Z][0-9]+[A-Z]" #this looks for a base, followed by a number that can be any size 1 or more, followed by a base.
  #"\\*[0-9]+[A-Z]" #this looks for an asterisk, followed by at least 1 number, followed by a base
  #"[A-Z][0-9]*\\*" #this looks for a base, followed by at least 1 number, followed by an asterisk
  
  coding = sum(grepl("coding",Protein$AA)) #Breseq's coding region
  intergenic = sum(grepl("intergenic",Protein$AA)) #intergenic breseq designation
  pseudogene = sum(grepl("pseudogene",Protein$AA)) #pseudogene breseq designation
  prematurestop = sum(lengths(regmatches(Protein$AA,gregexpr("[A-Z][0-9]*\\*", Protein$AA)))) #these are when you have a coding amino acid that gets mutated into a stop codon
  elongating = sum(lengths(regmatches(Protein$AA,gregexpr("\\*[0-9]*[A-Z]", Protein$AA)))) #these are stop codons that get mutated into coding amino acids that elongate your protein.
  aamutation = sum(lengths(regmatches(Protein$AA,gregexpr("[A-Z][0-9]+[A-Z]", Protein$AA)))) # these are all of the mutations that dont fit other categories. so these mutations change one amino acid to another with no other breseq designation.
  
  #I now need to determine if the amino acid mutations category are synonymous or nonsynonymous. The above just determines the number of leter to leter strings exist. Now I need to actually look at that subset of mutations and determine if
  aas <- lengths(regmatches(Protein$AA,gregexpr("[A-Z][0-9]+[A-Z]", Protein$AA))) #this returns a list of logicals that tell you if there is an aamutation at that spot (1) or not (0)
  #aas
  aminos <- as.matrix(Protein$AA[aas==1]) #this find the amino acid changes at the previously found indexes, so these are the actual identities of the aamutation mutations.
  aminos2 <- colsplit(aminos, "[0-9]+", names = c("first", "last")) # splitting the breseq annotation. I am taking out the number, because the position is irrelevent, and separating the two letters into two columns containing the first, original aa, and the last, or mutated, aa.
  
  synonymous <- sum(ifelse(aminos2$first ==aminos2$last,1,0)) #if the letters are the same before and after the mutation it is synonymous
  nonsynonymous <- sum(ifelse(aminos2$first == aminos2$last,0,1)) #if the letters are different then it is nonsynonymous
  dnds <- nonsynonymous/synonymous
  
  # I am now making a table of all of the mutational types that I can print out later. The other thing that I would like to do is give it a name specific to the data set that you put in, but I don't know how to do that. For the moment it will just always return the same named file each time, so you have to change the name before you use the code again. 
  table<- matrix(c("Mutations: ", allmutations,
                   "Nucleotide level mutations", "",
                   "Indels: ", indel,
                   "Transitions: ",transitions,
                   "C>T: ", CT,
                   "T>C: ", TC,
                   "A>G: ", AG,
                   "G>A: ", GA,
                   "Transversions: ", transversions,
                   "A>T: ", AT,
                   "A>C: ", AC,
                   "G>C: ", GC,
                   "G>T: ", GT,
                   "C>A: ", CA,
                   "C>G: ", CG,
                   "T>G: ", TG,
                   "T>A: ", TA,
                   "Amino acid level mutations", "",
                   "Coding: ", coding,
                   "Intergenic: ", intergenic,
                   "Pseudogene: ", pseudogene,
                   "Premature stop: ", prematurestop,
                   "Elongating: ", elongating,
                   "Synonymous: ", synonymous,
                   "Non Synonymous: ", nonsynonymous,
                   "dN/dS: ", dnds), ncol = 2, byrow=T)
  
  #write out the file. Adding col.names won't actually give the columns names but it will keep things in two different columns instead of compiling it all.
  write.csv(table, file = "Mutations_table.csv", col.names = T)
  
}

setwd("/Users/katrina/Desktop/working/B3/")
Mutations_analysis(B3, "Description", "Mutation")


#for running Deitrick's scripts
#python muller.py -i /Users/katrina/Desktop/working/B2/B2_muller.xlsx -o /Users/katrina/Desktop/working/B2/genotypes
















#see what the last filter took out to make sure the mutations are what I want them to be. This check needs to be done on all mutations that are taken out from step 1 on.
filter62_removed <- (subset(B3_filter5, abs(B3_filter5$"17"-B3_filter5$"90") < 10))
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
write.csv(split_removed2, file = "B3_removed.csv")


#######################################
#the following is a section of MatLab code (so won't work to run in R)

#problem I ran into is that there is apparently something wrong with the excel file... it won't read in
#manipulations done to excel file once exported
#change column names of frequencies to "number" type
#change the frequencies to "percentage" type
#
names = ["B3_muller.xlsx"]; %take in the excel file you got from R
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



csvwrite("/Users/katrina/Desktop/working/B3/B3timepoints.csv", timepoints)
csvwrite("/Users/katrina/Desktop/working/B3/B3frequencies.csv",frequencies)
csvwrite("/Users/katrina/Desktop/working/B3/B3nests.csv",nests)
