#looking at the amount of paralellism seen in the six evolved populations. 



#parallelism studies

######
#operon level parallelism within populations
######
B1_operon_occur <- data.frame(table(B1$operon))
B1_operon_parallelism <- B1_operon_occur[B1_operon_occur$Freq>1,]

B2_operon_occur <- data.frame(table(B2$operon))
B2_operon_parallelism <- B2_operon_occur[B2_operon_occur$Freq>1,]

B3_operon_occur <- data.frame(table(B3$operon))
B3_operon_parallelism <- B3_operon_occur[B3_operon_occur$Freq>1,]

P1_operon_occur <- data.frame(table(P1$operon))
P1_operon_parallelism <- P1_operon_occur[P1_operon_occur$Freq>1,]

P2_operon_occur <- data.frame(table(P2$operon))
P2_operon_parallelism <- P2_operon_occur[P2_operon_occur$Freq>1,]

P3_operon_occur <- data.frame(table(P3$operon))
P3_operon_parallelism <- P3_operon_occur[P3_operon_occur$Freq>1,]

#write out files.
write.csv(B1_operon_parallelism, "/Users/katrina/Desktop/muller_v_0.3/parallelism/operon/B1_operon_parallelism.csv")
write.csv(B2_operon_parallelism, "/Users/katrina/Desktop/muller_v_0.3/parallelism/operon/B2_operon_parallelism.csv")
write.csv(B3_operon_parallelism, "/Users/katrina/Desktop/muller_v_0.3/parallelism/operon/B3_operon_parallelism.csv")
write.csv(P1_operon_parallelism, "/Users/katrina/Desktop/muller_v_0.3/parallelism/operon/P1_operon_parallelism.csv")
write.csv(P2_operon_parallelism, "/Users/katrina/Desktop/muller_v_0.3/parallelism/operon/P2_operon_parallelism.csv")
write.csv(P3_operon_parallelism, "/Users/katrina/Desktop/muller_v_0.3/parallelism/operon/P3_operon_parallelism.csv")


###### 
#operon level parallelism between populations 
#######
all_operon_occur <- data.frame(table(all_pops$operon))
all_operon_parallelism <- all_operon_occur[all_operon_occur$Freq>1,]
write.csv(all_operon_parallelism, "/Users/katrina/Desktop/muller_v_0.3/parallelism/operon/all_operon.csv")

#####
#table of operon level parallelism within populations
#######
library("vegan")
library("plyr")
library("RColorBrewer")
library("ggplot2")
library("data.table")
library("dplyr")
library("reshape2")
library("scales")

all_operon_parallelism$Pop <- "all"
B1_operon_parallelism$Pop <- "B1"
B2_operon_parallelism$Pop <- "B2"
B3_operon_parallelism$Pop <- "B3"
P1_operon_parallelism$Pop <- "P1"
P2_operon_parallelism$Pop <- "P2"
P3_operon_parallelism$Pop <- "P3"

operon_parallelism_within <- rbind(all_operon_parallelism, B1_operon_parallelism, B2_operon_parallelism, B3_operon_parallelism, P1_operon_parallelism, P2_operon_parallelism, P3_operon_parallelism)

m_operon_paralellism_within <- melt(operon_parallelism_within,id=c("Pop","Var1"),measure.vars="Freq")
operon_within_cast <- t(dcast(m_operon_paralellism_within,Pop~Var1,mean,value.var="value",fill=0))
operon_within_final <- as.data.frame(operon_within_cast, header=T)
colnames(operon_within_final) <- as.character(unlist(operon_within_cast[1,]))
operon_within_final <- operon_within_final[-1,]

nrow(operon_within_final)#114

#print whole data frame:remember this is within population
write.csv(operon_within_final, "/Users/katrina/Desktop/muller_v_0.3/parallelism/operon/operon_level_within_pops.csv")

######
#print out a table of all cases of operon level parallelism across populations.
#####
operon_of_interest <- rownames(operon_within_final)

B1_operon_interest <- as.data.frame(B1_operon_occur[(B1_operon_occur$Var1 %in% operon_of_interest),])
B1_operon_interest$Pop <- "B1"
B2_operon_interest <- as.data.frame(B2_operon_occur[(B2_operon_occur$Var1 %in% operon_of_interest),])
B2_operon_interest$Pop <- "B2"
B3_operon_interest <- as.data.frame(B3_operon_occur[(B3_operon_occur$Var1 %in% operon_of_interest),])
B3_operon_interest$Pop <- "B3"

P1_operon_interest <- as.data.frame(P1_operon_occur[(P1_operon_occur$Var1 %in% operon_of_interest),])
P1_operon_interest$Pop <- "P1"
P2_operon_interest <- as.data.frame(P2_operon_occur[(P2_operon_occur$Var1 %in% operon_of_interest),])
P2_operon_interest$Pop <- "P2"
P3_operon_interest <- as.data.frame(P3_operon_occur[(P3_operon_occur$Var1 %in% operon_of_interest),])
P3_operon_interest$Pop <- "P3"

#combine into one big table - remember this is between pops
operon_between_pops <- rbind(all_operon_parallelism, B1_operon_interest, B2_operon_interest, B3_operon_interest, P1_operon_interest, P2_operon_interest, P3_operon_interest)

m_operon_between <- melt(operon_between_pops,id=c("Pop","Var1"),measure.vars="Freq")
operon_between_cast <- t(dcast(m_operon_between,Pop~Var1,mean,value.var="value",fill=0))
operon_between_final <- as.data.frame(operon_between_cast, header = T)
colnames(operon_between_final) <- as.character(unlist(operon_between_cast[1,]))
operon_between_final <- operon_between_final[-1,]

write.csv(operon_between_final, "/Users/katrina/Desktop/muller_v_0.3/parallelism/operon/operon_between_pops.csv")

#find functional identifiers for the different operons with parallelism
targetoperons <- rownames(operon_between_final)
operon_identity <- all_pops[(all_pops$operon %in% targetoperons),]
operon_identity <- operon_identity[,c(17,20)] #only keep operon and pseudocap function
colnames(operon_identity)
write.csv(operon_identity, "/Users/katrina/Desktop/muller_v_0.3/parallelism/operon/operon_identity.csv")
#this doesn't actually work, because there are genes with different predicted functions within operons. So with 174 mutations I need plus an extra ~20 functions.


#####
#gene level paralellism within populations (need to change paths so they are also in gene folder within parallelism)
########### 
#find genes that are found more than once
B1_gene_occur <- data.frame(table(B1$Gene)) #Find all of the unique instances and the number of times they are found. 
B1_gene_paralellism <- B1_gene_occur[B1_gene_occur$Freq>1,] #only take the ones found more than once

B2_gene_occur <- data.frame(table(B2$Gene))
B2_gene_paralellism <- B2_gene_occur[B2_gene_occur$Freq>1,]

B3_gene_occur <- data.frame(table(B3$Gene))
B3_gene_paralellism <- B3_gene_occur[B3_gene_occur$Freq>1,]

P1_gene_occur <- data.frame(table(P1$Gene))
P1_gene_paralellism <- P1_gene_occur[P1_gene_occur$Freq>1,]

P2_gene_occur <- data.frame(table(P2$Gene))
P2_gene_paralellism <- P2_gene_occur[P2_gene_occur$Freq>1,]

P3_gene_occur <- data.frame(table(P3$Gene))
P3_gene_paralellism <- P3_gene_occur[P3_gene_occur$Freq>1,]


#write out files for the gene level paralellism found in all populations
write.csv(B1_gene_paralellism, "/Users/katrina/Desktop/muller_v_0.3/parallelism/gene/B1_gene_paralellism.csv")
write.csv(B2_gene_paralellism, "/Users/katrina/Desktop/muller_v_0.3/parallelism/gene/B2_gene_paralellism.csv")
write.csv(B3_gene_paralellism, "/Users/katrina/Desktop/muller_v_0.3/parallelism/gene/B3_gene_paralellism.csv")
write.csv(P1_gene_paralellism, "/Users/katrina/Desktop/muller_v_0.3/parallelism/gene/P1_gene_paralellism.csv")
write.csv(P2_gene_paralellism, "/Users/katrina/Desktop/muller_v_0.3/parallelism/gene/P2_gene_paralellism.csv")
write.csv(P3_gene_paralellism, "/Users/katrina/Desktop/muller_v_0.3/parallelism/gene/P3_gene_paralellism.csv")


#######
#gene level paralellism between populations
######
#now looking at all of the mutations as one batch

all_pops_occur <- data.frame(table(all_pops$Gene))
all_gene_paralellism <- all_pops_occur[all_pops_occur$Freq>1,]

#View(all_paralellism)
write.csv(all_gene_paralellism, "/Users/katrina/Desktop/muller_v_0.3/parallelism/gene/all_pops_gene.csv")

#######
#find the numbers of all functional categories in all populations
#######
all_pops_function <- data.frame(table(P3$pseudocap.function))
write.csv(all_pops_function, "/Users/katrina/Desktop/muller_v_0.3/parallelism/function/all_pops_functions.csv")

#########
#Create a table of all gene level paralellism data showing the number in all populations, and the number in the subsequent popuations. 
##########
all_paralell <- as.data.frame(read.csv("/Users/katrina/Desktop/muller_v_0.3/parallelism/gene/all_pops_gene.csv", stringsAsFactors = FALSE, header=TRUE))
B1_paralell <- as.data.frame(read.csv("/Users/katrina/Desktop/muller_v_0.3/parallelism/gene/B1_gene_paralellism.csv", stringsAsFactors = FALSE, header=TRUE))
B2_paralell <- as.data.frame(read.csv("/Users/katrina/Desktop/muller_v_0.3/parallelism/gene/B2_gene_paralellism.csv", stringsAsFactors = FALSE, header=TRUE))
B3_paralell <- as.data.frame(read.csv("/Users/katrina/Desktop/muller_v_0.3/parallelism/gene/B3_gene_paralellism.csv", stringsAsFactors = FALSE, header=TRUE))
P1_paralell <- as.data.frame(read.csv("/Users/katrina/Desktop/muller_v_0.3/parallelism/gene/P1_gene_paralellism.csv", stringsAsFactors = FALSE, header=TRUE))
P2_paralell <- as.data.frame(read.csv("/Users/katrina/Desktop/muller_v_0.3/parallelism/gene/P2_gene_paralellism.csv", stringsAsFactors = FALSE, header=TRUE))
P3_paralell <- as.data.frame(read.csv("/Users/katrina/Desktop/muller_v_0.3/parallelism/gene/P3_gene_paralellism.csv", stringsAsFactors = FALSE, header=TRUE))



#View(paralellism)
library("vegan")
library("plyr")
library("RColorBrewer")
library("ggplot2")
library("data.table")
library("dplyr")
library("reshape2")
library("scales")


all_paralell$Pop <- "all"
B1_paralell$Pop <- "B1"
B2_paralell$Pop <- "B2"
B3_paralell$Pop <- "B3"
P1_paralell$Pop <- "P1"
P2_paralell$Pop <- "P2"
P3_paralell$Pop <- "P3"

#head(B1_paralell)
all_paralell <- all_paralell[,-1]
B1_paralell <- B1_paralell[,-1]
B2_paralell <- B2_paralell[,-1]
B3_paralell <- B3_paralell[,-1]
P1_paralell <- P1_paralell[,-1]
P2_paralell <- P2_paralell[,-1]
P3_paralell <- P3_paralell[,-1]

paralellism <- rbind(all_paralell, B1_paralell,B2_paralell,B3_paralell,P1_paralell,P2_paralell,P3_paralell)
#View(paralellism)

#melt data frame for casting
m_paralellism <- melt(paralellism,id=c("Pop","Var1"),measure.vars="Freq")

paralellism_cast <- t(dcast(m_paralellism,Pop~Var1,mean,value.var="value",fill=0))
paralellism_final <- as.data.frame(paralellism_cast, header=TRUE)
colnames(paralellism_final) <- as.character(unlist(paralellism_cast[1,]))
paralellism_final <- paralellism_final[-1,]

nrow(paralellism_final) #165
#print whole data frame: this is within population
write.csv(paralellism_final, "/Users/katrina/Desktop/muller_v_0.3/parallelism/gene/within_population_paralellism_table.csv", col.names = TRUE)


#####
#now to find all of the populations that have between population gene paralellism
#######
#I can use the tables of gene level paralellism from the previous code here
genes_of_interest <- rownames(paralellism_final)
#View(genes_of_interest)


B1_of_interest <- B1_gene_occur[(B1_gene_occur$Var1 %in% genes_of_interest),]
B1_of_interest$Pop <- "B1"
B2_of_interest <- B2_gene_occur[(B2_gene_occur$Var1 %in% genes_of_interest),]
B2_of_interest$Pop <- "B2"
B3_of_interest <- B3_gene_occur[(B3_gene_occur$Var1 %in% genes_of_interest),]
B3_of_interest$Pop <- "B3"
P1_of_interest <- P1_gene_occur[(P1_gene_occur$Var1 %in% genes_of_interest),]
P1_of_interest$Pop <- "P1"
P2_of_interest <- P2_gene_occur[(P2_gene_occur$Var1 %in% genes_of_interest),]
P2_of_interest$Pop <- "P2"
P3_of_interest <- P3_gene_occur[(P3_gene_occur$Var1 %in% genes_of_interest),]
P3_of_interest$Pop <- "P3"

#combine to be able to make the big table
#head(all_paralellism)

paralellism_all <- rbind(all_paralell, B1_of_interest,B2_of_interest,B3_of_interest,P1_of_interest, P2_of_interest, P3_of_interest)

m_paralellism_all <- melt(paralellism_all,id=c("Pop","Var1"),measure.vars="Freq")
paralellism_all_cast <- t(dcast(m_paralellism_all,Pop~Var1,mean,value.var="value",fill=0))
paralellism_all_final <- as.data.frame(paralellism_all_cast, header=TRUE)
colnames(paralellism_all_final) <- as.character(unlist(paralellism_all_cast[1,]))
paralellism_all_final <- paralellism_all_final[-1,]
#View(paralellism_all_final)  

write.csv(paralellism_all_final, "/Users/katrina/Desktop/muller_v_0.3/parallelism/gene/Gene_level_Paralellism_all_populations.csv", col.names = TRUE)


#find gene and function identifiers for the different genes with paralellism
targetgenes <- rownames(paralellism_all_final)
gene_identity <- all_pops[(all_pops$Gene %in% targetgenes),]
gene_identity <- gene_identity[,c(19,17,20)] #only want the operon and pseudocap function information
colnames(gene_identity)
write.csv(gene_identity, "/Users/katrina/Desktop/muller_v_0.3/parallelism/gene/gene_level_identity.csv")

#will need to add these to the printed out table separately.




#position level paralellism

######

#position level paralellism - across populations only
################
#position level paralellism

#now to look at the position level paralellism. For this all I am going to do is look at the all_pops data frame

#find the unique mutations in all of the populations together
mutation_occur <- data.frame(table(all_pops$Position))
mutation_paralellism <- mutation_occur[mutation_occur$Freq>1,]
#View(mutation_paralellism)
mutation_paralellism$Pop <- "all"
#find the unique mutations in each population
B1_position_occur <- data.frame(table(B1$Position)) 
B2_position_occur <- data.frame(table(B2$Position)) 
B3_position_occur <- data.frame(table(B3$Position)) 
P1_position_occur <- data.frame(table(P1$Position)) 
P2_position_occur <- data.frame(table(P2$Position)) 
P3_position_occur <- data.frame(table(P3$Position)) 

#FIND ALL OF THE INSTANCES OF THE MUTATION PARALELLISM IN THE INDIVIDUAL POPULATIONS
B1_mut_want <- B1_position_occur[(B1_position_occur$Var1 %in% mutation_paralellism$Var1),]
B1_mut_want$Pop <- "B1"
B2_mut_want <- B2_position_occur[(B2_position_occur$Var1 %in% mutation_paralellism$Var1),]
B2_mut_want$Pop <- "B2"
B3_mut_want <- B3_position_occur[(B3_position_occur$Var1 %in% mutation_paralellism$Var1),]
B3_mut_want$Pop <- "B3"
P1_mut_want <- P1_position_occur[(P1_position_occur$Var1 %in% mutation_paralellism$Var1),]
P1_mut_want$Pop <- "P1"
P2_mut_want <- P2_position_occur[(P2_position_occur$Var1 %in% mutation_paralellism$Var1),]
P2_mut_want$Pop <- "P2"
P3_mut_want <- P3_position_occur[(P3_position_occur$Var1 %in% mutation_paralellism$Var1),]
P3_mut_want$Pop <- "P3"

#combine to make a big table
mut_paralell <- rbind(mutation_paralellism, B1_mut_want, B2_mut_want, B3_mut_want, P1_mut_want, P2_mut_want, P3_mut_want)

#head(mut_paralell)

m_mut_paralell <- melt(mut_paralell,id=c("Pop","Var1"),measure.vars="Freq" )
mut_paralellism_cast <-t(dcast(m_mut_paralell,Pop~Var1,mean,value.var="value",fill=0))
mut_paralell_final <- as.data.frame(mut_paralellism_cast, header=T)
colnames(mut_paralell_final) <- as.character(unlist(mut_paralellism_cast[1,]))
mut_paralell_final <- mut_paralell_final[-1,]

#View(mut_paralell_final)

mut_paralell_final$Gene <- "" #wish I knew how to assign a gene name for the position, will do manually for the sake of time


write.csv(mut_paralell_final, "/Users/katrina/Desktop/muller_v_0.3/parallelism/position/Mutation_level_Paralellism_all_populations.csv")
#View(B3)


#find gene and function identifiers for the different positions with paralellism 
identity <- all_pops[(all_pops$Position %in% mutation_paralellism$Var1),]
colnames(identity)
identity <- identity[,c(11,12,16,17,19,20)]

#head(identity)
write.csv(identity, "/Users/katrina/Desktop/muller_v_0.3/parallelism/position/mut_level_identity.csv")

#####
######


#doing within environment analysis


#operon level parallelism within environments
######
#make dfs for both environemnts
biofilm <- rbind(B1,B2,B3)
plank <- rbind(P1,P2,P3)

#find parallelism for biofilm environemnt
biofilm_operon_occur <- data.frame(table(biofilm$operon)) #will find all unique operons and give a count of how many times its seen
biofilm_operon_parallelism <- biofilm_operon_occur[(biofilm_operon_occur$Freq>1),] #gets all of the operons that are seen more than once, so all of the cases where there is parallelism

#find parallelism for planktonic enivroment
plank_operon_occur <- data.frame(table(plank$operon))
plank_operon_parallelism <- plank_operon_occur[(plank_operon_occur$Freq>1),]

#find parallelism for all of the populations together
all_operon_occur <- data.frame(table(all_pops$operon))
all_operon_parallel <- all_operon_occur[(all_operon_occur$Freq > 1),]

#write out all files. 
write.csv(biofilm_operon_parallelism, "/Users/katrina/Desktop/muller_v_0.3/parallelism/by_environment/operon/biofilm_operon_parallelism.csv")
write.csv(plank_operon_parallelism, "/Users/katrina/Desktop/muller_v_0.3/parallelism/by_environment/operon/plank_operon_parallelism.csv")
write.csv(all_operon_parallel, "/Users/katrina/Desktop/muller_v_0.3/parallelism/by_environment/operon/all_operon_parallelism.csv")



#table of within environment operon level parallelism 
all_operon_parallel <- as.data.frame(all_operon_parallel)
all_operon_parallel$Pop <- "all"

plank_operon_parallelism <- as.data.frame(plank_operon_parallelism)
plank_operon_parallelism$Pop <- "Planktonic"

biofilm_operon_parallelism <- as.data.frame(biofilm_operon_parallelism)
biofilm_operon_parallelism$Pop <- "Biofilm"

library("reshape2")
by_env <- rbind(all_operon_parallel, plank_operon_parallelism, biofilm_operon_parallelism)

m_by_env <- melt(by_env, id=c("Pop", "Var1"), measure.vars="Freq")
by_env_cast <- t(dcast(m_by_env, Pop~Var1, mean, value.var="value",fill=0))
by_env_final <- as.data.frame(by_env_cast, header = T)
colnames(by_env_final) <- as.character(unlist(by_env_cast[1,]))
by_env_final <- by_env_final[-1,]

#nrow(by_env_final) #171

#print whithin population data frame. 
write.csv(by_env_final, "/Users/katrina/Desktop/muller_v_0.3/parallelism/by_environment/operon/operon_parallel_within_environment.csv")


#######
#find cases of operon parallelism across environments. 
#####

operon_of_interest <- rownames(by_env_final)
biofilm_interest <- as.data.frame(biofilm_operon_occur[(biofilm_operon_occur$Var1 %in% operon_of_interest),])
biofilm_interest$Pop <- "Biofilm"
plank_interest <- as.data.frame(plank_operon_occur[(plank_operon_occur$Var1 %in% operon_of_interest),])
plank_interest$Pop <- "Planktonic"

between_pops <- rbind(all_operon_parallel, biofilm_interest, plank_interest)
m_between_pops <- melt(between_pops, id=c("Pop", "Var1"), measure.vars = "Freq")
between_pops_cast <- t(dcast(m_between_pops,Pop~Var1, mean, value.var="value",fill=0))
between_pops_final <- as.data.frame(between_pops_cast, header=T)
colnames(between_pops_final) <- as.character((unlist(between_pops_cast[1,])))
between_pops_final <- between_pops_final[-1,]

write.csv(between_pops_final, "/Users/katrina/Desktop/muller_v_0.3/parallelism/by_environment/operon/operon_parallel_across_environments.csv")












#####

######
#gene level parallelism within environments
######
#find parallelism for biofilm environemnt

biofilm_gene_occur <- data.frame(table(biofilm$Gene)) #will find all unique genes and give a count of how many times its seen
biofilm_gene_parallelism <- biofilm_gene_occur[(biofilm_gene_occur$Freq>1),] #gets all of the genes that are seen more than once, so all of the cases where there is parallelism

#find parallelism for planktonic enivroment
plank_gene_occur <- data.frame(table(plank$Gene))
plank_gene_parallelism <- plank_gene_occur[(plank_gene_occur$Freq>1),]

#find parallelism for all of the populations together
all_gene_occur <- data.frame(table(all_pops$Gene))
all_gene_parallel <- all_gene_occur[(all_gene_occur$Freq > 1),]

#write out all files. 
write.csv(biofilm_gene_parallelism, "/Users/katrina/Desktop/muller_v_0.3/parallelism/by_environment/gene/biofilm_gene_parallelism.csv")
write.csv(plank_gene_parallelism, "/Users/katrina/Desktop/muller_v_0.3/parallelism/by_environment/gene/plank_gene_parallelism.csv")
write.csv(all_gene_parallel, "/Users/katrina/Desktop/muller_v_0.3/parallelism/by_environment/gene/all_gene_parallelism.csv")



#table of within environment gene level parallelism 
all_gene_parallel <- as.data.frame(all_gene_parallel)
all_gene_parallel$Pop <- "all"

plank_gene_parallelism <- as.data.frame(plank_gene_parallelism)
plank_gene_parallelism$Pop <- "Planktonic"

biofilm_gene_parallelism <- as.data.frame(biofilm_gene_parallelism)
biofilm_gene_parallelism$Pop <- "Biofilm"

library("reshape2")
by_env <- rbind(all_gene_parallel, plank_gene_parallelism, biofilm_gene_parallelism)

m_by_env <- melt(by_env, id=c("Pop", "Var1"), measure.vars="Freq")
by_env_cast <- t(dcast(m_by_env, Pop~Var1, mean, value.var="value",fill=0))
by_env_final <- as.data.frame(by_env_cast, header = T)
colnames(by_env_final) <- as.character(unlist(by_env_cast[1,]))
by_env_final <- by_env_final[-1,]

#nrow(by_env_final) #165

#print whithin population data frame. 
write.csv(by_env_final, "/Users/katrina/Desktop/muller_v_0.3/parallelism/by_environment/gene/gene_parallel_within_environment.csv")


#######
#find cases of gene parallelism across environments. 
#####

gene_of_interest <- rownames(by_env_final)
biofilm_interest <- as.data.frame(biofilm_gene_occur[(biofilm_gene_occur$Var1 %in% gene_of_interest),])
biofilm_interest$Pop <- "Biofilm"
plank_interest <- as.data.frame(plank_gene_occur[(plank_gene_occur$Var1 %in% gene_of_interest),])
plank_interest$Pop <- "Planktonic"

between_pops <- rbind(all_gene_parallel, biofilm_interest, plank_interest)
m_between_pops <- melt(between_pops, id=c("Pop", "Var1"), measure.vars = "Freq")
between_pops_cast <- t(dcast(m_between_pops,Pop~Var1, mean, value.var="value",fill=0))
between_pops_final <- as.data.frame(between_pops_cast, header=T)
colnames(between_pops_final) <- as.character((unlist(between_pops_cast[1,])))
between_pops_final <- between_pops_final[-1,]

#nrow(between_pops_final)#165 <- sanity check
write.csv(between_pops_final, "/Users/katrina/Desktop/muller_v_0.3/parallelism/by_environment/gene/gene_parallel_across_environments.csv")

####find gene and function identifiers for the different genes with parallelism
targetgenes <- rownames(between_pops_final)
gene_identity <- all_pops[(all_pops$Gene %in% targetgenes),]
gene_identity <- gene_identity[,c(10,11,8)]
#head(gene_identity)
write.csv(gene_identity, "/Users/katrina/Desktop/muller_v_0.3/parallelism/by_environment/gene/gene_identity.csv")






######
#find numbers of pseudocap functions and put into a table
########
#now finding the pseudocap functional class numbers
B1_function_occur <- data.frame(table(B1$pseudocap.function))
B2_function_occur <- data.frame(table(B2$pseudocap.function))
B3_function_occur <- data.frame(table(B3$pseudocap.function))
P1_function_occur <- data.frame(table(P1$pseudocap.function))
P2_function_occur <- data.frame(table(P2$pseudocap.function))
P3_function_occur <- data.frame(table(P3$pseudocap.function))

#write the functional analysis out to files
write.csv(B1_function_occur, "/Users/katrina/Desktop/muller_v_0.3/parallelism/function/B1_pseudocap_functions.csv")
write.csv(B2_function_occur, "/Users/katrina/Desktop/muller_v_0.3/parallelism/function/B2_pseudocap_functions.csv")
write.csv(B3_function_occur, "/Users/katrina/Desktop/muller_v_0.3/parallelism/function/B3_pseudocap_functions.csv")
write.csv(P1_function_occur, "/Users/katrina/Desktop/muller_v_0.3/parallelism/function/P1_pseudocap_functions.csv")
write.csv(P2_function_occur, "/Users/katrina/Desktop/muller_v_0.3/parallelism/function/P2_pseudocap_functions.csv")
write.csv(P3_function_occur, "/Users/katrina/Desktop/muller_v_0.3/parallelism/function/P3_pseudocap_functions.csv")

#figure out how many there are for the whole batch
all_function_occur <- data.frame(table(all_pops$pseudocap.function))
all_function_occur$Pop <- "all"
B1_function_occur$Pop <- "B1"
B2_function_occur$Pop <- "B2"
B3_function_occur$Pop <- "B3"
P1_function_occur$Pop <- "P1"
P2_function_occur$Pop <- "P2"
P3_function_occur$Pop <- "P3"


function_classes <- rbind(all_function_occur, B1_function_occur, B2_function_occur, B3_function_occur, P1_function_occur, P2_function_occur, P3_function_occur)

#head(function_classes)

m_function_classes <- melt(function_classes, id=c("Pop","Var1"), measure.vars="Freq")
function_classes_cast <- t(dcast(m_function_classes, Pop~Var1,mean,Value.var="value",fill=0))
function_classes_final <- as.data.frame(function_classes_cast, header=T)
colnames(function_classes_final) <- as.character(unlist(function_classes_cast[1,]))
function_classes_final <- function_classes_final[-1,]

#head(function_classes_final)

write.csv(function_classes_final, "/Users/katrina/Desktop/muller_v_0.3/parallelism/function/pseudocap_function_numbers.csv")


#I added to this file in excel. I added the total number of genes found in the UCBPP PA14 genome in each one of these functional categories. This will allow me to determine if I get certain categories called more often than would be expected purely by chance.
########

#trying to determine if you see any functional categories hit more often than would be expected by chance. Incomplete methods here - don't know how to do necessarily
#######
#--> this is only really an issue in functional cateogries. In genes/positions and even operons, given the mutation rates and the size of the genome, if mutations were completely random, you would never expect to see one gene mutated more than once.


#pseudocap_all <- read.csv("/Users/katrina/Desktop/muller_v_0.3/parallelism/function/complete_pseudocap_function_numbers.csv", stringsAsFactors = F)

#pseudocap_all[2,2]

# for anosim fucntion the data needs rows beig the samples and columns being the different variables. so the read in table needs to be transposed
t_pseudocap_all <- read.csv("/Users/katrina/Desktop/muller_v_0.3/parallelism/function/t_complete_pseudocap.csv", stringsAsFactors = F)

#need to take out the sample identifier row
t_pseudocap_all <- t_pseudocap_all[,-1]
View(t_pseudocap_all)

#it still feels like I am missing something. There has to be a way to incorporate the nmber of mutations seen total into account or the number of genes in the genome. 

#made data frames for each population individually 
B1_pseudocap <- t_pseudocap_all[,c(1,3,9)]
B2_pseudocap <- pseudocap_all[,c(1,4,9)]
B3_pseudocap <- pseudocap_all[,c(1,5,9)]
P1_pseudocap <- pseudocap_all[,c(1,6,9)]
P2_pseudocap <- pseudocap_all[,c(1,7,9)]
P3_pseudocap <- pseudocap_all[,c(1,8,9)]

#generate the euclidean distance for each populaiton against the total found in the genome. 
B1_euclidean_dist <- dist(B1_pseudocap, method = "euclidean")
B2_euclidean_dist <- dist(B2_pseudocap, method = "euclidean")
B3_euclidean_dist <- dist(B3_pseudocap, method = "euclidean")
P1_euclidean_dist <- dist(P1_pseudocap, method = "euclidean")
P2_euclidean_dist <- dist(P2_pseudocap, method = "euclidean")
P3_euclidean_dist <- dist(P3_pseudocap, method = "euclidean")


all_euclidean <- dist(t_pseudocap_all, method = "euclidean")
#the distance matricies now have to be visualized somehow. Kenny used wards minimum variance method. I am using ward.D2 because it seems like that is the more accepted ward equation. 

B1_clustering <- hclust(B1_euclidean_dist, method="ward.D2", members = NULL)
head(B1_clustering)

plot(B1_clustering)

all_clustering <- hclust(all_euclidean, method = "ward.D2", members = NULL)
plot(all_clustering)

#determine the euclidean distance matrix for all
anosim(t_pseudocap_all,grouping=c("all","B","B","B","P","P","P","PA14genome"), permutations = 999, distance = "euclidean")

anosim(t_pseudocap_all[c(2,3,4,8),], grouping= c("b","b","b", "PA14"), permutations = 999, distance = "euclidean")
anosim(t_pseudocap_all[c(5,6,7,8),], grouping= c("p","p","p", "PA14"), permutations = 999, distance = "euclidean")

anosim(t_pseudocap_all[c(2,3,4,5,6,7),], grouping= c("b","b","b", "p","p","p"), permutations = 999, distance = "euclidean")

#to use adonis, you take the output of the dist function 
adonis()



#######
#want to find cases where there is parallelism and they are found in the comparison analysis. 

#####
#read in the B1 comparison mutations with the metagenomic info
B1_comparison <- read.csv("/Users/katrina/Desktop/B1_comparison/filtered/d90_metainfo.csv", stringsAsFactors = F)
B2_comparison <- read.csv("/Users/katrina/Desktop/B2_comparison/d90_metainfo.csv", stringsAsFactors = F)

head(B1_comparison)
#select the ones that are only in the metagenome 
B1_comparison_metagenome <- B1_comparison[(B1_comparison$metagenome == "yes"),]
B2_comparison_metagenome <- B2_comparison[(B2_comparison$metagenome == "yes"),]

#read in the position level parallelism table
position_parallelism <- read.csv("/Users/katrina/Google Drive/CooperLab/Katrina/PA/PA_LTE/final_files (as of 6_23_19)/Parallelism/position/Mutation_level_Paralellism_all_populations.csv", stringsAsFactors = F)

#find the cases where there is a position level match for the B1 comparison and the parallelism 

position_targets <- position_parallelism[(position_parallelism$Position %in% B2_comparison_metagenome$Position %in% B1_comparison_metagenome),]
#there are none... now we will try the gene level


gene_parallelism <- read.csv("/Users/katrina/Google Drive/CooperLab/Katrina/PA/PA_LTE/final_files (as of 6_23_19)/Parallelism/gene/Gene_level_Paralellism_all_populations.csv", stringsAsFactors = F)

#issue is that the gene names are different in this file from what I am using in the comparison analysis. 
#import the edited data set that has both the default gene names and the new ones. 
all_edited <- read.csv("/Users/katrina/Desktop/muller_v_0.3/paralellism_all/all_operon_intergenic_info.csv", stringsAsFactors = F)
head(all_edited)
gene_names <- all_edited[(all_edited$Gene %in% gene_parallelism$X),]

#find cases where there is gene level parallelism between all three data sets. 
gene_targets_B1 <- gene_names[(gene_names$Gene.via.gbk.file %in% B1_comparison_metagenome$Gene),]

gene_targets_B1_B2 <- gene_targets_B1[(gene_targets_B1$Gene.via.gbk.file %in% B2_comparison_metagenome$Gene),]

#save the tables to files.
write.csv(gene_targets_B1_B2, "/Users/katrina/Desktop/B2_comparison/gene_parallelism_B1_B2.csv")

#read in the paired down file 
final_gene_parallelism <- read.csv("/Users/katrina/Desktop/B2_comparison/gene_paralellism_meta_comparison.csv", stringsAsFactors = F)

#now get the origninal gene parallelism data to know which populations they are in. 

gene_parallelism_pops <- gene_parallelism[(gene_parallelism$X %in% final_gene_parallelism$Gene),]

write.csv(gene_parallelism_pops, "/Users/katrina/Desktop/B2_comparison/gene_paralellism_meta_comparison_data.csv")




