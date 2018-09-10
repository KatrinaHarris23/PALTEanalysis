#PALTE mutational cohort plots

B1_pretty <- read.csv("/Users/katrina/Desktop/latest_PALTE/B1_pretty.csv")
B2_pretty <- read.csv("/Users/katrina/Desktop/latest_PALTE//B2_pretty.csv")
B3_pretty <- read.csv("/Users/katrina/Desktop/latest_PALTE//B3_pretty.csv")
P1_pretty <- read.csv("/Users/katrina/Desktop/latest_PALTE//P1_pretty.csv")
P2_pretty <- read.csv("/Users/katrina/Desktop/latest_PALTE//P2_pretty.csv")
P3_pretty <- read.csv("/Users/katrina/Desktop/latest_PALTE//P3_pretty.csv")



##mutational plots. first grab all of the frequencies and then plot them. 
#View(B1_pretty)
#ncol(B1_pretty) #13

B1_NT <- B1_pretty[,7:13]
B2_NT <- B2_pretty[,7:13]
B3_NT <- B3_pretty[,7:13]


#ncol(P1_pretty) #11
P1_NT <- P1_pretty[,7:11]
P2_NT <- P2_pretty[,7:11]
P3_NT <- P3_pretty[,7:11]


#change the column names so that we are working in generations
colnames(B1_NT) <- c("0","113","167","293","440","500","600")
colnames(B2_NT) <- c("0","113","167","293","440","500","600")
colnames(B3_NT) <- c("0","113","167","293","440","500","600")

colnames(P1_NT) <- c("0","113","293","500","600")
colnames(P2_NT) <- c("0","113","293","500","600")
colnames(P3_NT) <- c("0","113","293","500","600")



b_generations <- c(0,113,167,293,440,500,600)
p_generations <- c(0,113,293,500,600)


View(P1_NT)

layout(matrix(c(1,2,3,4,5,6),2))


plot(NA, xlim=c(0,600), ylim=c(0,100), ylab = "Frequency", xlab = "Generations", main = "B1") #, main = "Mutation frequencies over time") # 
for(i in seq_len(nrow(B1_NT))){
  lines(b_generations,B1_NT[i,], type="l", col=i)
} #need to get them colored differently, but will work for now. 


plot(NA, xlim=c(0,600), ylim=c(0,100), ylab = "", xlab = "", main = "P1") #, main = "Mutation frequencies over time") # 
for(i in seq_len(nrow(P1_NT))){
  lines(p_generations,P1_NT[i,], type="l", col=i)
} #need to get them colored differently, but will work for now. 


plot(NA, xlim=c(0,600), ylim=c(0,100), ylab = "", xlab = "", main = "B2") #, main = "Mutation frequencies over time") # 
for(i in seq_len(nrow(B2_NT))){
  lines(b_generations,B2_NT[i,], type="l", col=i)
} #need to get them colored differently, but will work for now. 

plot(NA, xlim=c(0,600), ylim=c(0,100), ylab = "", xlab = "", main = "P2") #, main = "Mutation frequencies over time") # 
for(i in seq_len(nrow(P2_NT))){
  lines(p_generations,P2_NT[i,], type="l", col=i)
} #need to get them colored differently, but will work for now.

plot(NA, xlim=c(0,600), ylim=c(0,100), ylab = "", xlab = "", main = "B3") #, main = "Mutation frequencies over time") # 
for(i in seq_len(nrow(B3_NT))){
  lines(b_generations,B3_NT[i,], type="l", col=i)
} #need to get them colored differently, but will work for now. 

plot(NA, xlim=c(0,600), ylim=c(0,100), ylab = "", xlab = "", main = "P3") #, main = "Mutation frequencies over time") # 
for(i in seq_len(nrow(P3_NT))){
  lines(p_generations,P3_NT[i,], type="l", col=i)
} #need to get them colored differently, but will work for now. 





#####Now plotting the cohort trajectories

B1_cohort_trajectories <- read.csv("/Users/katrina/Desktop/working/B1/genotypes/B1_muller.trajectories.csv", sep = "\t")
B2_cohort_trajectories <- read.csv("/Users/katrina/Desktop/working/B2/genotypes/B2_muller.trajectories.csv", sep = "\t")
B3_cohort_trajectories <- read.csv("/Users/katrina/Desktop/working/B3/genotypes/B3_muller.trajectories.csv", sep = "\t")
P1_cohort_trajectories <- read.csv("/Users/katrina/Desktop/working/P1/genotypes/P1_muller.trajectories.csv", sep = "\t")
P2_cohort_trajectories <- read.csv("/Users/katrina/Desktop/working/P2/genotypes/P2_muller.trajectories.csv", sep = "\t")
P3_cohort_trajectories <- read.csv("/Users/katrina/Desktop/working/P3/genotypes/P3_muller.trajectories.csv", sep = "\t")

B1_cohort_trajectories <- B1_cohort_trajectories[,4:10]
B2_cohort_trajectories <- B2_cohort_trajectories[,4:10]
B3_cohort_trajectories <- B3_cohort_trajectories[,4:10]

P1_cohort_trajectories <- P1_cohort_trajectories[,4:8]
P2_cohort_trajectories <- P2_cohort_trajectories[,4:8]
P3_cohort_trajectories <- P3_cohort_trajectories[,4:8]





#B1_cohorts_edited <- B1_cohorts_edited*100
layout(matrix(c(1,2,3,4,5,6),2))

plot(NA, xlim=c(0,600), ylim=c(0,1), ylab = "Frequency", xlab = "Generations", main = "B1", cex.axis=2, cex.lab=2, cex.main=3) #, main = "Mutation frequencies over time") # 
for(i in seq_len(nrow(B1_cohort_trajectories))){
  lines(b_generations,B1_cohort_trajectories[i,], type="l", col=1)
}
#lines(b_generations, B1_cohorts[4,], type="l", col = "red", lwd=3)

#need to get them colored differently, but will work for now. 

plot(NA, xlim=c(0,600), ylim=c(0,1), ylab = "", xlab = "", main = "P1", cex.axis=2, cex.main=3) #, main = "Mutation frequencies over time") # 
for(i in seq_len(nrow(P1_cohort_trajectories))){
  lines(p_generations,P1_cohort_trajectories[i,], type="l", col=1)
} #need to get them colored differently, but will work for now. 

plot(NA, xlim=c(0,600), ylim=c(0,1), ylab = "", xlab = "", main = "B2", cex.axis=2, cex.main=3) #, main = "Mutation frequencies over time") # 
for(i in seq_len(nrow(B2_cohort_trajectories))){
  lines(b_generations,B2_cohort_trajectories[i,], type="l", col=1)
} #need to get them colored differently, but will work for now. 
#lines(b_generations,b2_cohort_lineage[8,], type="l", col="red", lwd=3)

plot(NA, xlim=c(0,600), ylim=c(0,1), ylab = "", xlab = "", main = "P2", cex.axis=2, cex.main=3) #, main = "Mutation frequencies over time") # 
for(i in seq_len(nrow(P2_cohort_trajectories))){
  lines(p_generations,P2_cohort_trajectories[i,], type="l", col=1)
} #need to get them colored differently, but will work for now. 

plot(NA, xlim=c(0,600), ylim=c(0,1), ylab = "", xlab = "", main = "B3", cex.axis=2, cex.main=3) #, main = "Mutation frequencies over time") # 
for(i in seq_len(nrow(B3_cohort_trajectories))){
  lines(b_generations,B3_cohort_trajectories[i,], type="l", col=1)
} #need to get them colored differently, but will work for now. 
#lines(b_generations, b3_cohort_lineage[5,], type="l", col="red", lwd=3)

plot(NA, xlim=c(0,600), ylim=c(0,1), ylab = "", xlab = "", main = "P3", cex.axis=2, cex.main=3) #, main = "Mutation frequencies over time") # 
for(i in seq_len(nrow(P3_cohort_trajectories))){
  lines(p_generations,P3_cohort_trajectories[i,], type="l", col=1)
} #need to get them colored differently, but will work for now. 
#lines(p_generations, p3_cohort_lineage[4,], type="l", col="red", lwd=3)





#### now taking all of the mutation breakdown tables and putting them in one table

B1_mut_breakdown <- read.csv("/Users/katrina/Desktop/working/B1/B1_Mutations_table.csv")
B2_mut_breakdown <- read.csv("/Users/katrina/Desktop/working/B2/B2_Mutations_table.csv")
B3_mut_breakdown <- read.csv("/Users/katrina/Desktop/working/B3/B3_Mutations_table.csv")
P1_mut_breakdown <- read.csv("/Users/katrina/Desktop/working/P1/P1_Mutations_table.csv")
P2_mut_breakdown <- read.csv("/Users/katrina/Desktop/working/P2/Mutations_table.csv")
P3_mut_breakdown <- read.csv("/Users/katrina/Desktop/working/P3/P3_Mutations_table.csv")


muttable <- B1_mut_breakdown[,2:3]
names = c("Category", "B1")
colnames(muttable) <- names

muttable$B2 <- as.numeric(B2_mut_breakdown$V2, round = 2)
muttable$B3 <- B3_mut_breakdown$V2
muttable$P1 <- P1_mut_breakdown$V2
muttable$P2 <- P2_mut_breakdown$V2
muttable$P3 <- P3_mut_breakdown$V2
View(muttable)

write.csv(muttable, file = "/Users/katrina/Desktop/latest_PALTE/PALTE_mutational_breakdown.csv")



###now to see if I can extract the cohort members from the rows of frequency data 
######I need to figure out how to code this 


B1_cohorts <- read.csv("/Users/katrina/Desktop/latest_PALTE/B1_muller.genotypemembers.csv")
B2_cohorts <- read.csv("/Users/katrina/Desktop/latest_PALTE/B2_muller.genotypemembers.csv")
B3_cohorts <- read.csv("/Users/katrina/Desktop/latest_PALTE/B3_muller.genotypemembers.csv")
P1_cohorts <- read.csv("/Users/katrina/Desktop/latest_PALTE/P1_muller.genotypemembers.csv")
P2_cohorts <- read.csv("/Users/katrina/Desktop/latest_PALTE/P2_muller.genotypemembers.csv")
P3_cohorts <- read.csv("/Users/katrina/Desktop/latest_PALTE/P3_muller.genotypemembers.csv")






#B1_indexing <- sapply(B1_cohorts, function(x) length(x)) #get the number of members in each cohort

#i = 1 #set i to 1
#for(i in seq_len(length(B1_indexing))){ #iterate through the indexing list
#  new <- list(rep(i, B1_indexing[i])) #say that you want the number of the #index printed the same number of times as their are items in that train car #of the list.
#  B1_index <- c(B1_index,new) #update the index to include the new numbers
#}


#B1_indecies <- unlist(B1_index)


#B1_lineage_cohort <- cbind(as.numeric(B1_cohortmembers),as.numeric(B1_indecies))


#B1_sorted_lineage_cohort <- B1_lineage_cohort[order(B1_lineage_cohort[,1]),] # this sorts based on the trajectory



#B1_pretty_cohorts <- B1_pretty #creating a new data set instead of just editing the last one
#B1_pretty_cohorts$Cohorts <- B1_sorted_lineage_cohort[,2] #add the cohort number to the data set.
#write.csv(B1_pretty_cohorts, file = "/Users/katrina/Desktop/working/PALTE/B1_with_cohorts.csv")



#B1_freqs <- B1_pretty_cohorts[,7:13]
#B1_zero <- colSums(B1_freqs!=0.0) #adds up all of the non zero numbers in the columns. gives a named list


#Biofilm_freqs <- cbind(B1_zero,B2_zero,B3_zero)
#Planktonic_freqs <- cbind(P1_zero,P2_zero,P3_zero)



########
#trying to get all of the populations together in one big matrix, well, one for biofilm and one for planktonic. 

B1 <- read.csv("/Users/katrina/Desktop/latest_PALTE/B1_cohorts.csv", sep = ",", fill = T)
B2 <- read.csv("/Users/katrina/Desktop/latest_PALTE/B2_cohorts.csv", sep = ",", fill = T)
B3 <- read.csv("/Users/katrina/Desktop/latest_PALTE/B3_cohorts.csv", sep = ",", fill = T)

P1 <- read.csv("/Users/katrina/Desktop/latest_PALTE/P1_cohorts.csv", sep = ",", fill = T)
P2 <- read.csv("/Users/katrina/Desktop/latest_PALTE/P2_cohorts.csv", sep = ",", fill = T)
P3 <- read.csv("/Users/katrina/Desktop/latest_PALTE/P3_cohorts.csv", sep = ",", fill = T)

View(B1)


B1 <- B1[-1]
B2 <- B2[-1]
B3 <- B3[-1]
P1 <- P1[-1]
P2 <- P2[-1]
P3 <- P3[-1]

B1$Population <- "B1"
B2$Population <- "B2"
B3$Population <- "B3"
P1$Population <- "P1"
P2$Population <- "P2"
P3$Population <- "P3"

ncol(P1)
ncol(P2)
ncol(P3)


PALTE_allbiofilms <- rbind(B1,B2,B3)
PALTE_allplanktonic <- rbind(P1,P2,P3)

write.csv(PALTE_allbiofilms, file = "/Users/katrina/Desktop/latest_PALTE/PALTE_allBiofilms.csv", row.names = F)
write.csv(PALTE_allplanktonic, file = "/Users/katrina/Desktop/latest_PALTE/PALTE_allPlanktonic.csv", row.names = F)

#‰ÛÔ and â€‘ were replaced with - in excel in these master files that have my notes in them


#I don't think I actually plotted the cohort trajectories....

B1_cohortsagain <- read.csv("/Users/katrina/Desktop/working/B1/genotypes/B1_muller.genotypes.csv", sep="")
B2_cohortsagain <- read.csv("/Users/katrina/Desktop/working/B2/genotypes/B2_muller.genotypes.csv", sep="")
B3_cohortsagain <- read.csv("/Users/katrina/Desktop/working/B3/genotypes/B3_muller.genotypes.csv", sep="")
P1_cohortsagain <- read.csv("/Users/katrina/Desktop/working/P1/genotypes/P1_muller.genotypes.csv", sep="")
P2_cohortsagain <- read.csv("/Users/katrina/Desktop/working/P2/genotypes/P2_muller.genotypes.csv", sep="")
P3_cohortsagain <- read.csv("/Users/katrina/Desktop/working/P3/genotypes/P3_muller.genotypes.csv", sep="")

#remove all of the members
B1_cohortsagain <- B1_cohortsagain[-8]
B2_cohortsagain <- B2_cohortsagain[-8]
B3_cohortsagain <- B3_cohortsagain[-8]
P1_cohortsagain <- P1_cohortsagain[-6]
P2_cohortsagain <- P2_cohortsagain[-6]
P3_cohortsagain <- P3_cohortsagain[-6]


#Now plot the cohorts

layout(matrix(c(1,2,3,4,5,6),2))

plot(NA, xlim=c(0,600), ylim=c(0,1), ylab = "Frequency", xlab = "Generations", main = "B1", cex.axis=2, cex.lab=2, cex.main=3) #, main = "Mutation frequencies over time") # 
for(i in seq_len(nrow(B1_cohortsagain))){
  lines(b_generations,B1_cohortsagain[i,], type="l", col=1)
}

plot(NA, xlim=c(0,600), ylim=c(0,1), ylab = "Frequency", xlab = "Generations", main = "P1", cex.axis=2, cex.lab=2, cex.main=3) #, main = "Mutation frequencies over time") # 
for(i in seq_len(nrow(P1_cohortsagain))){
  lines(p_generations,P1_cohortsagain[i,], type="l", col=1)
}

plot(NA, xlim=c(0,600), ylim=c(0,1), ylab = "Frequency", xlab = "Generations", main = "B2", cex.axis=2, cex.lab=2, cex.main=3) #, main = "Mutation frequencies over time") # 
for(i in seq_len(nrow(B2_cohortsagain))){
  lines(b_generations,B2_cohortsagain[i,], type="l", col=1)
}

plot(NA, xlim=c(0,600), ylim=c(0,1), ylab = "Frequency", xlab = "Generations", main = "P2", cex.axis=2, cex.lab=2, cex.main=3) #, main = "Mutation frequencies over time") # 
for(i in seq_len(nrow(P2_cohortsagain))){
  lines(p_generations,P2_cohortsagain[i,], type="l", col=1)
}

plot(NA, xlim=c(0,600), ylim=c(0,1), ylab = "Frequency", xlab = "Generations", main = "B3", cex.axis=2, cex.lab=2, cex.main=3) #, main = "Mutation frequencies over time") # 
for(i in seq_len(nrow(B3_cohortsagain))){
  lines(b_generations,B3_cohortsagain[i,], type="l", col=1)
}

plot(NA, xlim=c(0,600), ylim=c(0,1), ylab = "Frequency", xlab = "Generations", main = "P3", cex.axis=2, cex.lab=2, cex.main=3) #, main = "Mutation frequencies over time") # 
for(i in seq_len(nrow(P3_cohortsagain))){
  lines(p_generations,P3_cohortsagain[i,], type="l", col=1)
}


write.csv(B1_cohortsagain, file = "/Users/katrina/Desktop/working/B1/genotypes/B1_cohortfreqs.csv")
write.csv(B2_cohortsagain, file = "/Users/katrina/Desktop/working/B2/genotypes/B2_cohortfreqs.csv")
write.csv(B3_cohortsagain, file = "/Users/katrina/Desktop/working/B3/genotypes/B3_cohortfreqs.csv")
write.csv(P1_cohortsagain, file = "/Users/katrina/Desktop/working/P1/genotypes/P1_cohortfreqs.csv")
write.csv(P2_cohortsagain, file = "/Users/katrina/Desktop/working/P2/genotypes/P2_cohortfreqs.csv")
write.csv(P3_cohortsagain, file = "/Users/katrina/Desktop/working/P3/genotypes/P3_cohortfreqs.csv")





#############plotting the number of cohorts that reach frequencies higher than 50%

#these values were derived from the files from the cohort data sets printed out above.
B1 <- c(0,3,4,11,11,17,17)
B2 <- c(0,2,6,8,11,16,16)
B3 <- c(0,0,0,1,3,3,3)
P1 <- c(0,0,1,1,1)
P2 <- c(0,0,1,2,2)
P3 <- c(0,0,1,2,2)

b_generations <- c(0,113,167,293,440,500,600)
p_generations <- c(0,113,293,500,600)

plot(b_generations, B1, type = "l", col = "red", xlab = "Time (Generations)", ylab = "Genotypes", main = "Number of high frequency genotypes", lwd = 2)
lines(b_generations, B2, type = "l", col = "orange", lwd = 2)
lines(b_generations, B3, type = "l", col = "purple", lwd = 2)
lines(p_generations, P1, type = "l", col = "black", lwd =  2)
lines(p_generations, P2, type = "l", col = "blue", lwd = 2)
lines(p_generations, P3, type = "l", col = "green", lwd = 2)
legend("topleft", legend = c("B1","B2","B3","P1","P2","P3"), lty = 1, col = c("red", "orange","purple", "black", "blue","green"), lwd = 2)








############### now to plot the number of high frequency mutations present. these were derived from the data sets that have cohorts associated with mutations above. 





####### determine how many mutations overlap between populations 

B1 <- read.csv("/Users/katrina/Desktop/latest_PALTE/B1_cohorts.csv", sep = ",", fill = T)
B2 <- read.csv("/Users/katrina/Desktop/latest_PALTE/B2_cohorts.csv", sep = ",", fill = T)
B3 <- read.csv("/Users/katrina/Desktop/latest_PALTE/B3_cohorts.csv", sep = ",", fill = T)

P1 <- read.csv("/Users/katrina/Desktop/latest_PALTE/P1_cohorts.csv", sep = ",", fill = T)
P2 <- read.csv("/Users/katrina/Desktop/latest_PALTE/P2_cohorts.csv", sep = ",", fill = T)
P3 <- read.csv("/Users/katrina/Desktop/latest_PALTE/P3_cohorts.csv", sep = ",", fill = T)

View(B1)

#position level parallelism

B1_B2_mut <- B1[ (B1$Position %in% B2$Position), ]
B2_B3_mut <- B2[ (B2$Position %in% B3$Position), ]
biofilm_sharedNT <- B1_B2_mut[(B1_B2_mut$Position %in% B2_B3_mut$Position),]
nrow(biofilm_sharedNT)#12
View(biofilm_sharedNT)
write.csv(biofilm_sharedNT, file = "/Users/katrina/Desktop/latest_PALTE/Biofilm_shared_nucleotide.csv")


P1_P2_mut <- P1[ (P1$Position %in% P2$Position), ]
P2_P3_mut <- P2[ (P2$Position %in% P3$Position), ]
plank_sharedNT <- P1_P2_mut[(P1_P2_mut$Position %in% P2_P3_mut$Position),]
write.csv(plank_sharedNT, file = "/Users/katrina/Desktop/latest_PALTE/Planktonic_shared_nucleotide.csv")

#these mutations now need to be taken out of the populations tables of mutations. 
nrow(B1)#228
B1_Keep <- B1[!(B1$Position %in% biofilm_sharedNT$Position),]
nrow(B1_Keep)#216 < yes, this is right because there are 12 shared mutations anmongst all of the populations.
B2_keep <- B2[!(B2$Position %in% biofilm_sharedNT$Position),]
B3_keep <- B3[!(B3$Position %in% biofilm_sharedNT$Position),]

P1_keep <- P1[!(P1$Position %in% plank_sharedNT$Position),]
P2_keep <- P2[!(P2$Position %in% plank_sharedNT$Position),]
P3_keep <- P3[!(P3$Position %in% plank_sharedNT$Position),]

#now write the files with all of these filters. 
write.csv(B1_Keep, file = "/Users/katrina/Desktop/latest_PALTE/B1_updated.csv")
write.csv(B2_keep, file = "/Users/katrina/Desktop/latest_PALTE/B2_updated.csv")
write.csv(B3_keep, file = "/Users/katrina/Desktop/latest_PALTE/B3_updated.csv")
write.csv(P1_keep, file = "/Users/katrina/Desktop/latest_PALTE/P1_updated.csv")
write.csv(P2_keep, file = "/Users/katrina/Desktop/latest_PALTE/P2_updated.csv")
write.csv(P3_keep, file = "/Users/katrina/Desktop/latest_PALTE/P3_updated.csv")





#see how many mutations are shared at the position level

b1_2 <- B1_Keep[(B1_Keep$Position %in% B2_keep$Position),]
b1_3 <- B1_Keep[(B1_Keep$Position %in% B3_keep$Position),]
biofilm <- b1_2[(b1_2$Position %in% b1_3$Position),] #no longer have overlap between all of the biofilm populations, which is good.


b1_p1 <- B1_Keep[(B1_Keep$Position %in% P1_keep$Position),]
b1_p2 <- B1_Keep[(B1_Keep$Position %in% P2_keep$Position),]
b1_p3 <- B1_Keep[(B1_Keep$Position %in% P3_keep$Position),]

b2_3 <- B2_keep[(B2_keep$Position %in% B3_keep$Position),]
b2_p1 <- B2_keep[(B2_keep$Position %in% P1_keep$Position),]
b2_p2 <- B2_keep[(B2_keep$Position %in% P2_keep$Position),]
b2_p3 <- B2_keep[(B2_keep$Position %in% P3_keep$Position),]

b3_p1 <- B3_keep[(B3_keep$Position %in% P1_keep$Position),]
b3_p2 <- B3_keep[(B3_keep$Position %in% P2_keep$Position),]
b3_p3 <- B3_keep[(B3_keep$Position %in% P3_keep$Position),]

p1_p2 <- P1_keep[(P1_keep$Position %in% P2_keep$Position),]
p1_p3 <- P1_keep[(P1_keep$Position %in% P3_keep$Position),]

p2_p3 <- P2_keep[(P2_keep$Position %in% P3_keep$Position),]






Table <- matrix(c("B1 and B2", nrow(b1_2), 
                  "B1 and B3", nrow(b1_3),
                  "B1 and P1", nrow(b1_p1),
                  "B1 and P2", nrow(b1_p2),
                  "B1 and P3", nrow(b1_p3),
                  "B2 and B3", nrow(b2_3),
                  "B2 and P1", nrow(b2_p1),
                  "B2 and P2", nrow(b2_p2),
                  "B2 and P3", nrow(b2_p3),
                  "B3 and P1", nrow(b3_p1),
                  "B3 and P2", nrow(b3_p2),
                  "B3 and P3", nrow(b3_p3),
                  "P1 and P2", nrow(p1_p2),
                  "P1 and P3", nrow(p1_p3),
                  "P2 and P3", nrow(p2_p3)), ncol = 2, byrow=T)
write.csv(Table, file = "/Users/katrina/Desktop/latest_PALTE/NT-levelparallelism.csv", col.names = T)

#make data frames that include all of these pairs

write.csv(b1_2, file = "/Users/katrina/Desktop/latest_PALTE/b1_2.csv")
write.csv(b1_3, file = "/Users/katrina/Desktop/latest_PALTE/b1_3.csv")
write.csv(b1_p1, file = "/Users/katrina/Desktop/latest_PALTE/b1_p1.csv")
write.csv(b1_p2, file = "/Users/katrina/Desktop/latest_PALTE/b1_p2.csv")
write.csv(b1_p3, file = "/Users/katrina/Desktop/latest_PALTE/b1_p3.csv")

write.csv(b2_3, file = "/Users/katrina/Desktop/latest_PALTE/b2_3.csv")
write.csv(b2_p1, file = "/Users/katrina/Desktop/latest_PALTE/b2_p1.csv")
write.csv(b2_p2, file = "/Users/katrina/Desktop/latest_PALTE/b2_p2.csv")
write.csv(b2_p3, file = "/Users/katrina/Desktop/latest_PALTE/b2_p3.csv")

write.csv(b3_p1, file = "/Users/katrina/Desktop/latest_PALTE/b3_p1.csv")
write.csv(b3_p2, file = "/Users/katrina/Desktop/latest_PALTE/b3_p2.csv")
write.csv(b3_p3, file = "/Users/katrina/Desktop/latest_PALTE/b3_p3.csv")

write.csv(p1_p2, file = "/Users/katrina/Desktop/latest_PALTE/p1_p2.csv")
write.csv(p1_p3, file = "/Users/katrina/Desktop/latest_PALTE/p1_p3.csv")

write.csv(p2_p3, file = "/Users/katrina/Desktop/latest_PALTE/p2_p3.csv")




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

setwd("/Users/katrina/Desktop/latest_PALTE/parallelism")

Mutations_analysis(b1_2, "Description", "Mutation")
Mutations_analysis(b1_3, "Description", "Mutation")
Mutations_analysis(b1_p1, "Description", "Mutation")

Mutations_analysis(b1_p2, "Description", "Mutation")
Mutations_analysis(b1_p3, "Description", "Mutation")

Mutations_analysis(b2_3, "Description", "Mutation")
Mutations_analysis(b2_p1, "Description", "Mutation")
Mutations_analysis(b2_p2, "Description", "Mutation")
Mutations_analysis(b2_p3, "Description", "Mutation")

Mutations_analysis(b3_p1, "Description", "Mutation")
Mutations_analysis(b3_p2, "Description", "Mutation")
Mutations_analysis(b3_p3, "Description", "Mutation")

Mutations_analysis(p1_p2, "Description", "Mutation")
Mutations_analysis(p1_p3, "Description", "Mutation")

Mutations_analysis(p2_p3, "Description", "Mutation")


b12 <- read.csv("b1_2Mutations_table.csv", stringsAsFactors = F)
b13 <- read.csv("b1_3Mutations_table.csv", stringsAsFactors = F)
b1p1 <- read.csv("b1_p1Mutations_table.csv", stringsAsFactors = F)
b1p2 <- read.csv("b1_p2Mutations_table.csv", stringsAsFactors = F)
b1p3 <- read.csv("b1_p3Mutations_table.csv", stringsAsFactors = F)
b23 <- read.csv("b2_3Mutations_table.csv", stringsAsFactors = F)
b2p1 <- read.csv("b2_p1Mutations_table.csv", stringsAsFactors = F)
b2p2 <- read.csv("b2_p2Mutations_table.csv", stringsAsFactors = F)
b2p3 <- read.csv("b2_p3Mutations_table.csv", stringsAsFactors = F)
p12 <- read.csv("p1_p2Mutations_table.csv", stringsAsFactors = F)
p13 <- read.csv("p1_p3Mutations_table.csv", stringsAsFactors = F)
p23 <- read.csv("p2_p3Mutations_table.csv", stringsAsFactors = F)

View(b12)

sumtable <- b12[,2:3]
colnames(sumtable) <- c("category", "b1b2")
sumtable$b1b3 <- b13$V2
sumtable$b1p1 <- b1p1$V2
sumtable$b1p2 <- b1p2$V2
sumtable$b1p3 <- b1p3$V2
sumtable$b23 <- b23$V2
sumtable$b2p1 <- b2p1$V2
sumtable$b2p2 <- b2p2$V2
sumtable$b2p3 <- b2p3$V2
sumtable$p12 <- p12$V2
sumtable$p13 <- p13$V2
sumtable$p23 <- p23$V2

View(sumtable)

write.csv(sumtable, file="sharedmutations.csv")


####before I move on I want to create one data frame that contains info for all of the parallelism between populations 


b1_2 <- b1_2[1:7]
b1_3 <- b1_3[1:7]
b1_p1 <- b1_p1[1:7]
b1_p2 <- b1_p2[1:7]
b1_p3 <- b1_p3[1:7]
b2_3 <- b2_3[1:7]
b2_p1 <- b2_p1[1:7]
b2_p2 <- b2_p2[1:7]
b2_p3 <- b2_p3[1:7]
b3_p1 <- b3_p1[1:7]
b3_p2 <- b3_p2[1:7]
b3_p3 <- b3_p3[1:7]
p1_p2 <- p1_p2[1:7]
p1_p3 <- p1_p3[1:7]
p2_p3 <- p2_p3[1:7]


b1_2$comparison <- "b1_2"
b1_3$comparison <- "b1_3"
b1_p1$comparison <- "b1_p1"
b1_p2$comparison <- "b1_p2"
b1_p3$comparison <- "b1_p3"
b2_3$comparison <- "b2_3"
b2_p1$comparison <- "b2_p1"
b2_p2$comparison <- "b2_p2"
b2_p3$comparison <- "b2_p3"
b3_p1$comparison <- "b3_p1"
b3_p2$comparison <- "b3_p2"
b3_p3$comparison <- "b3_p3"
p1_p2$comparison <- "p1_p2"
p1_p3$comparison <- "p1_p3"
p2_p3$comparison <- "p2_p3"

ncol(b1_2)
ncol(b1_p2)

x <- rbind(b1_2,b1_3,b1_p1,b1_p2,b1_p3)
z <- rbind(b2_3,b2_p1,b2_p2,b2_p3)
y <- rbind(b3_p1,b3_p2,b3_p3) 
q <- rbind(p1_p2,p1_p3,p2_p3)

dim(x)
dim(y)
dim(q)
dim(z)
colnames(x) <- c("X","Cohort","Description","Gene","Annotation","Position","Mutation","comparison")

all <- rbind(x,y,q,z)
View(all)
dim(all)

write.csv(all, file = "/Users/katrina/Desktop/latest_PALTE/parallelism/NTlevel_parallelism_calls.csv")




#### now to look at the gene level

b1_2 <- B1_Keep[(B1_Keep$Gene %in% B2_keep$Gene),]
b1_3 <- B1_Keep[(B1_Keep$Gene %in% B3_keep$Gene),]


b1_p1 <- B1_Keep[(B1_Keep$Gene %in% P1_keep$Gene),]
b1_p2 <- B1_Keep[(B1_Keep$Gene %in% P2_keep$Gene),]
b1_p3 <- B1_Keep[(B1_Keep$Gene %in% P3_keep$Gene),]

b2_3 <- B2_keep[(B2_keep$Gene %in% B3_keep$Gene),]
b2_p1 <- B2_keep[(B2_keep$Gene %in% P1_keep$Gene),]
b2_p2 <- B2_keep[(B2_keep$Gene %in% P2_keep$Gene),]
b2_p3 <- B2_keep[(B2_keep$Gene %in% P3_keep$Gene),]

b3_p1 <- B3_keep[(B3_keep$Gene %in% P1_keep$Gene),]
b3_p2 <- B3_keep[(B3_keep$Gene %in% P2_keep$Gene),]
b3_p3 <- B3_keep[(B3_keep$Gene %in% P3_keep$Gene),]

p1_p2 <- P1_keep[(P1_keep$Gene %in% P2_keep$Gene),]
p1_p3 <- P1_keep[(P1_keep$Gene %in% P3_keep$Gene),]

p2_p3 <- P2_keep[(P2_keep$Gene %in% P3_keep$Gene),]


Table <- matrix(c("B1 and B2", nrow(b1_2), 
                  "B1 and B3", nrow(b1_3),
                  "B1 and P1", nrow(b1_p1),
                  "B1 and P2", nrow(b1_p2),
                  "B1 and P3", nrow(b1_p3),
                  "B2 and B3", nrow(b2_3),
                  "B2 and P1", nrow(b2_p1),
                  "B2 and P2", nrow(b2_p2),
                  "B2 and P3", nrow(b2_p3),
                  "B3 and P1", nrow(b3_p1),
                  "B3 and P2", nrow(b3_p2),
                  "B3 and P3", nrow(b3_p3),
                  "P1 and P2", nrow(p1_p2),
                  "P1 and P3", nrow(p1_p3),
                  "P2 and P3", nrow(p2_p3)), ncol = 2, byrow=T)
write.csv(Table, file = "/Users/katrina/Desktop/latest_PALTE/gene-levelparallelism.csv", col.names = T)

#make data frames that include all of these pairs

write.csv(b1_2, file = "/Users/katrina/Desktop/latest_PALTE/b1_2.csv")
write.csv(b1_3, file = "/Users/katrina/Desktop/latest_PALTE/b1_3.csv")
write.csv(b1_p1, file = "/Users/katrina/Desktop/latest_PALTE/b1_p1.csv")
write.csv(b1_p2, file = "/Users/katrina/Desktop/latest_PALTE/b1_p2.csv")
write.csv(b1_p3, file = "/Users/katrina/Desktop/latest_PALTE/b1_p3.csv")

write.csv(b2_3, file = "/Users/katrina/Desktop/latest_PALTE/b2_3.csv")
write.csv(b2_p1, file = "/Users/katrina/Desktop/latest_PALTE/b2_p1.csv")
write.csv(b2_p2, file = "/Users/katrina/Desktop/latest_PALTE/b2_p2.csv")
write.csv(b2_p3, file = "/Users/katrina/Desktop/latest_PALTE/b2_p3.csv")

write.csv(b3_p1, file = "/Users/katrina/Desktop/latest_PALTE/b3_p1.csv")
write.csv(b3_p2, file = "/Users/katrina/Desktop/latest_PALTE/b3_p2.csv")
write.csv(b3_p3, file = "/Users/katrina/Desktop/latest_PALTE/b3_p3.csv")

write.csv(p1_p2, file = "/Users/katrina/Desktop/latest_PALTE/p1_p2.csv")
write.csv(p1_p3, file = "/Users/katrina/Desktop/latest_PALTE/p1_p3.csv")

write.csv(p2_p3, file = "/Users/katrina/Desktop/latest_PALTE/p2_p3.csv")


setwd("/Users/katrina/Desktop/latest_PALTE/parallelism")

Mutations_analysis(b1_2, "Description", "Mutation")
Mutations_analysis(b1_3, "Description", "Mutation")
Mutations_analysis(b1_p1, "Description", "Mutation")

Mutations_analysis(b1_p2, "Description", "Mutation")
Mutations_analysis(b1_p3, "Description", "Mutation")

Mutations_analysis(b2_3, "Description", "Mutation")
Mutations_analysis(b2_p1, "Description", "Mutation")
Mutations_analysis(b2_p2, "Description", "Mutation")
Mutations_analysis(b2_p3, "Description", "Mutation")

Mutations_analysis(b3_p1, "Description", "Mutation")
Mutations_analysis(b3_p2, "Description", "Mutation")
Mutations_analysis(b3_p3, "Description", "Mutation")

Mutations_analysis(p1_p2, "Description", "Mutation")
Mutations_analysis(p1_p3, "Description", "Mutation")

Mutations_analysis(p2_p3, "Description", "Mutation")


b12 <- read.csv("b1_2Mutations_table.csv", stringsAsFactors = F)
b13 <- read.csv("b1_3Mutations_table.csv", stringsAsFactors = F)
b1p1 <- read.csv("b1_p1Mutations_table.csv", stringsAsFactors = F)
b1p2 <- read.csv("b1_p2Mutations_table.csv", stringsAsFactors = F)
b1p3 <- read.csv("b1_p3Mutations_table.csv", stringsAsFactors = F)
b23 <- read.csv("b2_3Mutations_table.csv", stringsAsFactors = F)
b2p1 <- read.csv("b2_p1Mutations_table.csv", stringsAsFactors = F)
b2p2 <- read.csv("b2_p2Mutations_table.csv", stringsAsFactors = F)
b2p3 <- read.csv("b2_p3Mutations_table.csv", stringsAsFactors = F)
p12 <- read.csv("p1_p2Mutations_table.csv", stringsAsFactors = F)
p13 <- read.csv("p1_p3Mutations_table.csv", stringsAsFactors = F)
p23 <- read.csv("p2_p3Mutations_table.csv", stringsAsFactors = F)

View(b12)

sumtable <- b12[,2:3]
colnames(sumtable) <- c("category", "b1b2")
sumtable$b1b3 <- b13$V2
sumtable$b1p1 <- b1p1$V2
sumtable$b1p2 <- b1p2$V2
sumtable$b1p3 <- b1p3$V2
sumtable$b23 <- b23$V2
sumtable$b2p1 <- b2p1$V2
sumtable$b2p2 <- b2p2$V2
sumtable$b2p3 <- b2p3$V2
sumtable$p12 <- p12$V2
sumtable$p13 <- p13$V2
sumtable$p23 <- p23$V2

View(sumtable)

write.csv(sumtable, file="sharedGeneMutations.csv")



####before I move on I want to create one data frame that contains info for all of the parallelism between populations 


b1_2 <- b1_2[1:7]
b1_3 <- b1_3[1:7]
b1_p1 <- b1_p1[1:7]
b1_p2 <- b1_p2[1:7]
b1_p3 <- b1_p3[1:7]
b2_3 <- b2_3[1:7]
b2_p1 <- b2_p1[1:7]
b2_p2 <- b2_p2[1:7]
b2_p3 <- b2_p3[1:7]
b3_p1 <- b3_p1[1:7]
b3_p2 <- b3_p2[1:7]
b3_p3 <- b3_p3[1:7]
p1_p2 <- p1_p2[1:7]
p1_p3 <- p1_p3[1:7]
p2_p3 <- p2_p3[1:7]


b1_2$comparison <- "b1_2"
b1_3$comparison <- "b1_3"
b1_p1$comparison <- "b1_p1"
b1_p2$comparison <- "b1_p2"
b1_p3$comparison <- "b1_p3"
b2_3$comparison <- "b2_3"
b2_p1$comparison <- "b2_p1"
b2_p2$comparison <- "b2_p2"
b2_p3$comparison <- "b2_p3"
b3_p1$comparison <- "b3_p1"
b3_p2$comparison <- "b3_p2"
b3_p3$comparison <- "b3_p3"
p1_p2$comparison <- "p1_p2"
p1_p3$comparison <- "p1_p3"
p2_p3$comparison <- "p2_p3"

ncol(b1_2)
ncol(b1_p2)

x <- rbind(b1_2,b1_3,b1_p1,b1_p2,b1_p3)
z <- rbind(b2_3,b2_p1,b2_p2,b2_p3)
y <- rbind(b3_p1,b3_p2,b3_p3) 
q <- rbind(p1_p2,p1_p3,p2_p3)

dim(x)
dim(y)
dim(q)
dim(z)
colnames(x) <- c("X","Cohort","Description","Gene","Annotation","Position","Mutation","comparison")

all <- rbind(x,y,q,z)
View(all)
dim(all)

write.csv(all, file = "/Users/katrina/Desktop/latest_PALTE/parallelism/Genelevel_parallelism_calls.csv")






######################## determine the diversity of mutation types seen in the final files (so these are the _updated files)

B1 <- read.csv("/Users/katrina/Desktop/latest_PALTE/B1_updated.csv")
B2 <- read.csv("/Users/katrina/Desktop/latest_PALTE/B2_updated.csv")
B3 <- read.csv("/Users/katrina/Desktop/latest_PALTE/B3_updated.csv")
P1 <- read.csv("/Users/katrina/Desktop/latest_PALTE/P1_updated.csv")
P2 <- read.csv("/Users/katrina/Desktop/latest_PALTE/P2_updated.csv")
P3 <- read.csv("/Users/katrina/Desktop/latest_PALTE/P3_updated.csv")

setwd("/Users/katrina/Desktop/latest_PALTE/mutsummary")
Mutations_analysis(B1, "Description", "Mutation")
Mutations_analysis(B2, "Description", "Mutation")
Mutations_analysis(B3, "Description", "Mutation")
Mutations_analysis(P1, "Description", "Mutation")
Mutations_analysis(P2, "Description", "Mutation")
Mutations_analysis(P3, "Description", "Mutation")


setwd("/Users/katrina/Desktop/latest_PALTE/mutsummary")

B1_mutation_breakdown <- read.csv("B1_Mutations_table.csv", stringsAsFactors = F)
B2_mutation_breakdown <- read.csv("B2_Mutations_table.csv", stringsAsFactors = F)
B3_mutation_breakdown <- read.csv("B3_Mutations_table.csv", stringsAsFactors = F)
P1_mutation_breakdown <- read.csv("P1_Mutations_table.csv", stringsAsFactors = F)
P2_mutation_breakdown <- read.csv("P2_Mutations_table.csv", stringsAsFactors = F)
P3_mutation_breakdown <- read.csv("P3_Mutations_table.csv", stringsAsFactors = F)

muttable <- B1_mutation_breakdown[,2:3]
muttable$B2 <- B2_mutation_breakdown$V2
muttable$B3 <- B3_mutation_breakdown$V2
muttable$P1 <- P1_mutation_breakdown$V2
muttable$P2 <- P2_mutation_breakdown$V2
muttable$P3 <- P3_mutation_breakdown$V2

View(muttable)

columnnames <- c("Mutation Class", "B1","B2","B3","P1","P2","P3")
colnames(muttable)  <- columnnames 

View(muttable)
write.csv(muttable, file = "PALTE_mutational_breakdown.csv")













nrow(B1_final)
uniqueN(B1_final$Position)
which(duplicated(B1_final$Position)) #the rows at which the first duplicate happens.







############################################################
## I am going to try to figure out diversity first. I will look at between different time points in the B1 population. 

B1 <- read.csv("/Users/katrina/Desktop/latest_PALTE/B1_updated.csv", header = T, stringsAsFactors = F)

#trying to get the differences i diversity over time in the B1 population. 
B1_diversity <- B1[-c(1,2,3,4,5,6,7,8)]
t_B1_diversity <- t(B1_diversity)
B1_shannon <- diversity(t_B1_diversity, index="shannon", MARGIN = 1)
B1_simpson <- diversity(t_B1_diversity, index="simpson", MARGIN = 1)
B1_invsimpson <- diversity(t_B1_diversity, index="invsimpson", MARGIN = 1)

B2 <- read.csv("/Users/katrina/Desktop/latest_PALTE/B2_updated.csv", header = T, stringsAsFactors = F)

#trying to get the differences i diversity over time in the B1 population. 
t_B2_diversity <- t(B2[-c(1,2,3,4,5,6,7,8)])
B2_shannon <- diversity(t_B2_diversity, index="shannon", MARGIN = 1)
B2_simpson <- diversity(t_B2_diversity, index="simpson", MARGIN = 1)
B2_invsimpson <- diversity(t_B2_diversity, index="invsimpson", MARGIN = 1)

B3 <- read.csv("/Users/katrina/Desktop/latest_PALTE/B3_updated.csv", header = T, stringsAsFactors = F)

#trying to get the differences i diversity over time in the B1 population. 
t_B3_diversity <- t(B3[-c(1,2,3,4,5,6,7,8)])
B3_shannon <- diversity(t_B3_diversity, index="shannon", MARGIN = 1)
B3_simpson <- diversity(t_B3_diversity, index="simpson", MARGIN = 1)
B3_invsimpson <- diversity(t_B3_diversity, index="invsimpson", MARGIN = 1)


P1 <- read.csv("/Users/katrina/Desktop/latest_PALTE/P1_updated.csv", header = T, stringsAsFactors = F)

#trying to get the differences i diversity over time in the B1 population. 
t_P1_diversity <- t(P1[-c(1,2,3,4,5,6,7,8)])
P1_shannon <- diversity(t_P1_diversity, index="shannon", MARGIN = 1)
P1_simpson <- diversity(t_P1_diversity, index="simpson", MARGIN = 1)
P1_invsimpson <- diversity(t_P1_diversity, index="invsimpson", MARGIN = 1)

P2 <- read.csv("/Users/katrina/Desktop/latest_PALTE/P2_updated.csv", header = T, stringsAsFactors = F)

#trying to get the differences i diversity over time in the B1 population. 
t_P2_diversity <- t(P2[-c(1,2,3,4,5,6,7,8)])
P2_shannon <- diversity(t_P2_diversity, index="shannon", MARGIN = 1)
P2_simpson <- diversity(t_P2_diversity, index="simpson", MARGIN = 1)
P2_invsimpson <- diversity(t_P2_diversity, index="invsimpson", MARGIN = 1)


P3 <- read.csv("/Users/katrina/Desktop/latest_PALTE/P3_updated.csv", header = T, stringsAsFactors = F)

#trying to get the differences i diversity over time in the B1 population. 
t_P3_diversity <- t(P3[-c(1,2,3,4,5,6,7,8)])
P3_shannon <- diversity(t_P3_diversity, index="shannon", MARGIN = 1)
P3_simpson <- diversity(t_P3_diversity, index="simpson", MARGIN = 1)
P3_invsimpson <- diversity(t_P3_diversity, index="invsimpson", MARGIN = 1)

alphadiversitybiofilm <- matrix(c("Diversity Matrix", "0", "17", "25","44","66","75","90",
                           "B1 Shannons", t(B1_shannon),
                           "B1 Simpsons", t(B1_simpson),
                           "B1 invsimpson", t(B1_invsimpson),
                           "B2 Shannons", t(B2_shannon),
                           "B2 Simpsons", t(B2_simpson),
                           "B2 invsimpson", t(B2_invsimpson),
                           "B3 Shannons", t(B3_shannon),
                           "B3 Simpsons", t(B3_simpson),
                           "B3 invsimpson", t(B3_invsimpson)), ncol=8, byrow=T)
                           
alphadiversityplanktonic <- matrix(c("Diversity Matrix", "0", "17","44","66","90",
                           "P1 Shannons", t(P1_shannon),
                           "P1 Simpsons", t(P1_simpson),
                           "P1 invsimpson", t(P1_invsimpson),
                           "P2 Shannons", t(P2_shannon),
                           "P2 Simpsons", t(P2_simpson),
                           "P2 invsimpson", t(P2_invsimpson),
                           "P3 Shannons", t(P3_shannon),
                           "P3 Simpsons", t(P3_simpson),
                           "P3 invsimpson", t(P3_invsimpson)), ncol=6, byrow=T)

View(alphadiversitybiofilm)
View(alphadiversityplanktonic)
write.csv(alphadiversitybiofilm, file = "/Users/katrina/Desktop/latest_PALTE/AlphaDiversityBiofilm.csv")
write.csv(alphadiversityplanktonic, file = "/Users/katrina/Desktop/latest_PALTE/AlphaDiversityPlanktonic.csv")

View(P3)



################################################## not working
#Now to compare the 90 day time points for all of the data
day90 <- read.csv("/Users/katrina/Desktop/latest_PALTE/day_90.csv", header = T, stringsAsFactors = F)
dim(day90)
colnames(day90) <- c("Position", "B1","B2","B3","P1","P2","P3")
day90 <- day90[3:740,2:7]



day90 <- t(day90)
View(day90)

write.csv(day90, file = "/Users/katrina/Desktop/latest_PALTE/t_day90.csv")

table = read.csv("/Users/katrina/Desktop/latest_PALTE/t_day90.csv")
rownames(table) <- table$X
table <- subset(table, select = -c(X))
shannonDay90 <- diversity(table, index = "shannon", MARGIN = 1)
simpsonDay90 <- diversity(table, index = "simpson", MARGIN = 1)
invsimpsonDay90 <- diversity(table, index = "invsimpson", MARGIN = 1)



View(shannonDay90)

expshannon <- matrix(c(exp(shannonDay90[1]), exp(shannonDay90[2]), exp(shannonDay90[3]), exp(shannonDay90[4]), exp(shannonDay90[5]), exp(shannonDay90[6])), ncol = 6, byrow = T) #this is actually called the effective number of specieseffective number of
colnames(expshannon) <- c("B1","B2","B3","P1","P2","P3")


Day90diversity <- matrix(c("Diversity type", "B1","B2","B3","P1","P2","P3",
                           "Shannon",shannonDay90,
                           "expshannon", expshannon,
                           "Simpson", simpsonDay90,
                           "invsimpson", invsimpsonDay90), ncol =7, byrow=T )
View(Day90diversity)

write.csv(Day90diversity, file = "/Users/katrina/Desktop/latest_PALTE/parallelism/AlphaDay90Diversity.csv")




typeof(day90[3,3])







#figure out how many mutations at each time point for all of the populations. 
View(B1)

B10 <- length(which(B1$X0 != 0))
B117 <- length(which(B1$X17 != 0))
B125 <- length(which(B1$X25 != 0))
B144 <- length(which(B1$X44 != 0))
B166 <- length(which(B1$X66 != 0))
B175 <- length(which(B1$X75 != 0))
B190 <- length(which(B1$X90 != 0))

B20 <- length(which(B2$X0 != 0))
B217 <- length(which(B2$X17 != 0))
B225 <- length(which(B2$X25 != 0))
B244 <- length(which(B2$X44 != 0))
B266 <- length(which(B2$X66 != 0))
B275 <- length(which(B2$X75 != 0))
B290 <- length(which(B2$X90 != 0))

B30 <- length(which(B3$X0 != 0))
B317 <- length(which(B3$X17 != 0))
B325 <- length(which(B3$X25 != 0))
B344 <- length(which(B3$X44 != 0))
B366 <- length(which(B3$X66 != 0))
B375 <- length(which(B3$X75 != 0))
B390 <- length(which(B3$X90 != 0))

P10 <- length(which(P1$X0 != 0))
P117 <- length(which(P1$X17 != 0))
P144 <- length(which(P1$X44 != 0))
P166 <- length(which(P1$X66 != 0))
P190 <- length(which(P1$X90 != 0))

P20 <- length(which(P2$X0 != 0))
P217 <- length(which(P2$X17 != 0))
P244 <- length(which(P2$X44 != 0))
P266 <- length(which(P2$X66 != 0))
P290 <- length(which(P2$X90 != 0))

P30 <- length(which(P3$X0 != 0))
P317 <- length(which(P3$X17 != 0))
P344 <- length(which(P3$X44 != 0))
P366 <- length(which(P3$X66 != 0))
P390 <- length(which(P3$X90 != 0))


b_generations <- c(0,113,167,293,440,500,600)
p_generations <- c(0,113,293,440,600)

plot(b_generations, c(B10, B117, B125, B144,B166,B175,B190), type = "l", col = "red", xlab = "Time (Generations)", ylab = "number of mutations", main = "Number of mutations over time", lwd = 2)
lines(b_generations, c(B20, B217, B225, B244,B266,B275,B290), type = "l", col = "orange", lwd = 2)
lines(b_generations, c(B30, B317, B325, B344,B366,B375,B390), type = "l", col = "purple", lwd = 2)
lines(p_generations, c(P10, P117, P144,P166,P190), type = "l", col = "black", lwd = 2)
lines(p_generations, c(P20, P217, P244,P266,P290), type = "l", col = "blue", lwd = 2)
lines(p_generations, c(P30, P317, P344,P366,P390), type = "l", col = "green", lwd = 2)
par(xpd=TRUE)
legend("topleft", legend = c("B1","B2","B3","P1","P2","P3"), lty = 1, col = c("red", "orange","purple", "black", "blue","green"), lwd = 2)





