#this script analyzes the mutational dynamics seen in the evolved pseudomonas aeruginosa populations for the PALTE project. 



######
### inport csv files for the different populations and all populations in one sheet, and print out files for each individual population.
######
#this file has all of the mutations, but I had not put all of the intergenic/operon information in it.
#all_pops <- read.csv("/Users/katrina/Desktop/muller_v_0.3/all_edited.csv", stringsAsFactors = FALSE)

#this file has all of the mutations after r filtering, and it has pa14 locus tags, operon labels and intergenic region information. This is in the format of the output of version 0.5.1 of the muller scripts
all_pops <- read.csv("/Users/katrina/Desktop/muller_v_0.3/paralellism_all/all_operon_intergenic_info.csv", stringsAsFactors = FALSE)

View(all_pops)


#these files are the same as the files 190424__.csv on the drive in CooperLab > Katrina > PA> PALTE > working files
B1 <- all_pops[(all_pops$Population == "B1"),]
write.csv(B1, "/Users/katrina/Desktop/muller_v_0.3/Population_files_gene_names/B1_final.csv")
B2 <- all_pops[(all_pops$Population == "B2"),]
write.csv(B2, "/Users/katrina/Desktop/muller_v_0.3/Population_files_gene_names/B2_final.csv")
B3 <- all_pops[(all_pops$Population == "B3"),]
write.csv(B3, "/Users/katrina/Desktop/muller_v_0.3/Population_files_gene_names/B3_final.csv")





#planktonic populations have 2 less time points measured, so I have to take out the empty columns in the data before printing or the muller plot scripts will incorporate these. It is easiest to take them out here. They are the 25 and 75 day time points
P1 <- all_pops[(all_pops$Population == "P1"),]
P1 <- P1[,-c(4,7)]
#colnames(P1)
#View(P1)
write.csv(P1, "/Users/katrina/Desktop/muller_v_0.3/Population_files_gene_names/P1_final.csv")
P2 <- all_pops[(all_pops$Population == "P2"),]
P2 <- P2[,-c(4,7)]
write.csv(P2, "/Users/katrina/Desktop/muller_v_0.3/Population_files_gene_names/P2_final.csv")
P3 <- all_pops[(all_pops$Population == "P3"),]
P3 <- P3[,-c(4,7)]
write.csv(P3, "/Users/katrina/Desktop/muller_v_0.3/Population_files_gene_names/P3_final.csv")

#######
#read in files from the 6 evolved populations
#######

B1 <- read.csv("/Users/katrina/Desktop/muller_v_0.3/Population_files_gene_names/B1_final.csv", stringsAsFactors = FALSE)
B2 <- read.csv("/Users/katrina/Desktop/muller_v_0.3/Population_files_gene_names/B2_final.csv", stringsAsFactors = FALSE)
B3 <- read.csv("/Users/katrina/Desktop/muller_v_0.3/Population_files_gene_names/B3_final.csv", stringsAsFactors = FALSE)
P1 <- read.csv("/Users/katrina/Desktop/muller_v_0.3/Population_files_gene_names/P1_final.csv", stringsAsFactors = FALSE)
P2 <- read.csv("/Users/katrina/Desktop/muller_v_0.3/Population_files_gene_names/P2_final.csv", stringsAsFactors = FALSE)
P3 <- read.csv("/Users/katrina/Desktop/muller_v_0.3/Population_files_gene_names/P3_final.csv", stringsAsFactors = FALSE)



######



#Analysis to look at myutational dynamics through time in the populations




#determine the number of mutations seen in each population at each time point
#this is the number of mutations seen above a frequency of 0 at each given time point.
########

#B1
B1_17 <- B1[!(B1$X17 == 0),]
B1_25 <- B1[!(B1$X25 == 0),]
B1_44 <- B1[!(B1$X44 == 0),]
B1_66 <- B1[!(B1$X66 == 0),]
B1_75 <- B1[!(B1$X75 == 0),]
B1_90 <- B1[!(B1$X90 == 0),]

B1_numbers <- matrix(ncol = 2, data = c("Time point", "Number of mutations","0",0L, "17", nrow(B1_17), "25", nrow(B1_25), "44", nrow(B1_44),"66",nrow(B1_66), "75",nrow(B1_75),"90",nrow(B1_90)), byrow=T)
#View(B1_numbers)
#B2
B2_17 <- B2[!(B2$X17 == 0),]
B2_25 <- B2[!(B2$X25 == 0),]
B2_44 <- B2[!(B2$X44 == 0),]
B2_66 <- B2[!(B2$X66 == 0),]
B2_75 <- B2[!(B2$X75 == 0),]
B2_90 <- B2[!(B2$X90 == 0),]

B2_numbers <- matrix(ncol = 2, data = c("Time point", "Number of mutations","0",0L, "17", nrow(B2_17), "25", nrow(B2_25), "44", nrow(B2_44),"66",nrow(B2_66), "75",nrow(B2_75),"90",nrow(B2_90)), byrow=T)


#B3
B3_17 <- B3[!(B3$X17 == 0),]
B3_25 <- B3[!(B3$X25 == 0),]
B3_44 <- B3[!(B3$X44 == 0),]
B3_66 <- B3[!(B3$X66 == 0),]
B3_75 <- B3[!(B3$X75 == 0),]
B3_90 <- B3[!(B3$X90 == 0),]

B3_numbers <- matrix(ncol = 2, data = c("Time point", "Number of mutations","0",0L, "17", nrow(B3_17), "25", nrow(B3_25), "44", nrow(B3_44),"66",nrow(B3_66), "75",nrow(B3_75),"90",nrow(B3_90)), byrow=T)



#P1
P1_17 <- P1[!(P1$X17 == 0),]
P1_44 <- P1[!(P1$X44 == 0),]
P1_66 <- P1[!(P1$X66 == 0),]
P1_90 <- P1[!(P1$X90 == 0),]

P1_numbers <- matrix(ncol = 2, data = c("Time point", "Number of mutations","0",0L, "17", nrow(P1_17), "44", nrow(P1_44),"66",nrow(P1_66), "90",nrow(P1_90)), byrow=T)

#P2
P2_17 <- P2[!(P2$X17 == 0),]
P2_44 <- P2[!(P2$X44 == 0),]
P2_66 <- P2[!(P2$X66 == 0),]
P2_90 <- P2[!(P2$X90 == 0),]

P2_numbers <- matrix(ncol = 2, data = c("Time point", "Number of mutations","0",0L, "17", nrow(P2_17), "44", nrow(P2_44),"66",nrow(P2_66), "90",nrow(P2_90)), byrow=T)

#P3
P3_17 <- P3[!(P3$X17 == 0),]
P3_44 <- P3[!(P3$X44 == 0),]
P3_66 <- P3[!(P3$X66 == 0),]
P3_90 <- P3[!(P3$X90 == 0),]

P3_numbers <- matrix(ncol = 2, data = c("Time point", "Number of mutations","0",0L, "17", nrow(P3_17), "44", nrow(P3_44),"66",nrow(P3_66), "90",nrow(P3_90)), byrow=T)



#plot all population trajectories together. different color for all pops and different line type for different environments. 
plot(B1_numbers, type = "l", xlab="Days", ylab = "Number of mutations")
lines(B2_numbers, col = 2)
lines(B3_numbers, col = 3)
lines(P1_numbers, col = 4, lty = 2)
lines(P2_numbers, col = 5, lty = 2)
lines(P3_numbers, col = 6, lty = 2)
legend(x = "topleft", legend=c("B1", "B2","B3","P1","P2","P3"), lty = c(1,1,1,2,2,2), col=c(1,2,3,4,5,6))


#make large table combining all data points for the different populations. 

biofilm_numbers <- cbind(B1_numbers, B2_numbers[,2], B3_numbers[,2])
colnames(biofilm_numbers) <-c("Time point", "B1", "B2","B3")
biofilm_numbers <-biofilm_numbers[-1,]

planktonic_numbers <- cbind(P1_numbers, P2_numbers[,2], P3_numbers[,2])
colnames(planktonic_numbers) <- c("Time point","P1", "P2", "P3" )
planktonic_numbers <- planktonic_numbers[-1,]

#to combine all populations together, I need to introduce rows for the time points that were measured in biofilm but not in planktonic. If I don't do this, I won't be able to merge all of them into one table.
empty <- matrix(nrow=1, ncol = 4)
planktonic_numbers <- rbind(planktonic_numbers[1:2,], empty, planktonic_numbers[3:4,], empty, planktonic_numbers[5,])



all <- cbind(biofilm_numbers, planktonic_numbers[,2:4])

#write out the table to a file so that I can use it later. 
write.csv(all, "/Users/katrina/Desktop/muller_v_0.3/Dynamics/above_0_numbers.csv")

#what experiment is chris talking about that he did counting with alfonso's data? and what are the different methods? 

#retS gacS - have a gacA in pop data don't see others in clones or pops



#######
#now to look at mutations that FIX at the different time points
######

#B1
B1_17_fixed <- B1[(B1$X17 == 1),]
B1_25_fixed <- B1[(B1$X25 == 1),]
B1_44_fixed <- B1[(B1$X44 == 1),]
B1_66_fixed <- B1[(B1$X66 == 1),]
B1_75_fixed <- B1[(B1$X75 == 1),]
B1_90_fixed <- B1[(B1$X90 == 1),]
B1_fixed <- matrix(ncol = 2, data = c("Time point", "Number of mutations","0",0L, "17", nrow(B1_17_fixed), "25", nrow(B1_25_fixed), "44", nrow(B1_44_fixed),"66",nrow(B1_66_fixed), "75",nrow(B1_75_fixed),"90",nrow(B1_90_fixed)), byrow=T)

#B2
B2_17_fixed <- B2[(B2$X17 == 1),]
B2_25_fixed <- B2[(B2$X25 == 1),]
B2_44_fixed <- B2[(B2$X44 == 1),]
B2_66_fixed <- B2[(B2$X66 == 1),]
B2_75_fixed <- B2[(B2$X75 == 1),]
B2_90_fixed <- B2[(B2$X90 == 1),]
B2_fixed <- matrix(ncol = 2, data = c("Time point", "Number of mutations","0",0L, "17", nrow(B2_17_fixed), "25", nrow(B2_25_fixed), "44", nrow(B2_44_fixed),"66",nrow(B2_66_fixed), "75",nrow(B2_75_fixed),"90",nrow(B2_90_fixed)), byrow=T)

#B3
B3_17_fixed <- B3[(B3$X17 == 1),]
B3_25_fixed <- B3[(B3$X25 == 1),]
B3_44_fixed <- B3[(B3$X44 == 1),]
B3_66_fixed <- B3[(B3$X66 == 1),]
B3_75_fixed <- B3[(B3$X75 == 1),]
B3_90_fixed <- B3[(B3$X90 == 1),]
B3_fixed <- matrix(ncol = 2, data = c("Time point", "Number of mutations","0",0L, "17", nrow(B3_17_fixed), "25", nrow(B3_25_fixed), "44", nrow(B3_44_fixed),"66",nrow(B3_66_fixed), "75",nrow(B3_75_fixed),"90",nrow(B3_90_fixed)), byrow=T)

#P1
P1_17_fixed <- P1[(P1$X17 == 1),]
P1_44_fixed <- P1[(P1$X44 == 1),]
P1_66_fixed <- P1[(P1$X66 == 1),]
P1_90_fixed <- P1[(P1$X90 == 1),]
P1_fixed <- matrix(ncol = 2, data = c("Time point", "Number of mutations","0",0L, "17", nrow(P1_17_fixed), "44", nrow(P1_44_fixed),"66",nrow(P1_66_fixed),"90",nrow(P1_90_fixed)), byrow=T)

#P2
P2_17_fixed <- P2[(P2$X17 == 1),]
P2_44_fixed <- P2[(P2$X44 == 1),]
P2_66_fixed <- P2[(P2$X66 == 1),]
P2_90_fixed <- P2[(P2$X90 == 1),]
P2_fixed <- matrix(ncol = 2, data = c("Time point", "Number of mutations","0",0L, "17", nrow(P2_17_fixed), "44", nrow(P2_44_fixed),"66",nrow(P2_66_fixed),"90",nrow(P2_90_fixed)), byrow=T)

#P3
P3_17_fixed <- P3[(P3$X17 == 1),]
P3_44_fixed <- P3[(P3$X44 == 1),]
P3_66_fixed <- P3[(P3$X66 == 1),]
P3_90_fixed <- P3[(P3$X90 == 1),]
P3_fixed <- matrix(ncol = 2, data = c("Time point", "Number of mutations","0",0L, "17", nrow(P3_17_fixed), "44", nrow(P3_44_fixed),"66",nrow(P3_66_fixed),"90",nrow(P3_90_fixed)), byrow=T)


#plot them all
plot(B1_fixed, typ="l", xlim=c(0,90), ylim=c(0,25), xlab="Days", ylab = "number of mutations", main = "Fixed mutations")
lines(B2_fixed, typ="l", col=2)
lines(B3_fixed, typ="l", col=3)
lines(P1_fixed, typ="l", col=4, lty=2)
lines(P2_fixed, typ="l", col=5, lty=2) #is still on there but the P2 and P3 fixation rates are identical.
lines(P3_fixed, typ="l", col=6, lty=2) 
legend(x = "topleft", legend=c("B1", "B2","B3","P1","P2","P3"), lty = c(1,1,1,2,2,2), col=c(1,2,3,4,5,6))


#now to print out a table of all of the fixed mutations
biofilm_fixed <- cbind(B1_fixed, B2_fixed[,2], B3_fixed[,2])
colnames(biofilm_fixed) <- c("Time point", "B1", "B2","B3")
biofilm_fixed <- biofilm_fixed[-1,]

planktonic_fixed <- cbind(P1_fixed, P2_fixed[,2], P3_fixed[,2])
colnames(planktonic_fixed) <- c("Time point","P1", "P2", "P3" )
planktonic_fixed <- planktonic_fixed[-1,]

#to combine all populations together, I need to introduce rows for the time points that were measured in biofilm but not in planktonic. If I don't do this, I won't be able to merge all of them into one table.
empty <- matrix(nrow=1, ncol = 4)
planktonic_fixed <- rbind(planktonic_fixed[1:2,], empty, planktonic_fixed[3:4,], empty, planktonic_fixed[5,])

all_fixed <- cbind(biofilm_fixed, planktonic_fixed[,2:4])

#write out the table to a file so that I can use it later. 
write.csv(all_fixed, "/Users/katrina/Desktop/muller_v_0.3/Dynamics/all_fixed_mutations.csv")

#######
#now, to find the number of polymorphic mutations I need to subtract the number of fixed mutations, from the number of non zero mutations seen at each time point.
#####
B1_polymorphic <- matrix(ncol = 1,(as.numeric(B1_numbers[2:8,2])- as.numeric(B1_fixed[2:8,2])))
B2_polymorphic <- matrix(ncol = 1, (as.numeric(B2_numbers[2:8,2])- as.numeric(B2_fixed[2:8,2])))
B3_polymorphic <- matrix(ncol = 1, (as.numeric(B3_numbers[2:8,2])- as.numeric(B3_fixed[2:8,2])))

P1_polymorphic <- matrix(ncol = 1, (as.numeric(P1_numbers[2:6,2])- as.numeric(P1_fixed[2:6,2])))
P2_polymorphic <- matrix(ncol = 1, (as.numeric(P2_numbers[2:6,2])- as.numeric(P2_fixed[2:6,2])))
P3_polymorphic <- matrix(ncol = 1, (as.numeric(P3_numbers[2:6,2])- as.numeric(P3_fixed[2:6,2])))

plot(B1_polymorphic, x=c(0,17,25,44,66,75,90), xlab = "Day", ylab = "Number of mutations", main = "Polymorophic mutations", type = "l")
lines(B2_polymorphic, x=c(0,17,25,44,66,75,90), col = 2)
lines(B3_polymorphic, x=c(0,17,25,44,66,75,90), col = 3)
lines(P1_polymorphic, x = c(0,17,44,66,90), col = 4, lty=2)
lines(P2_polymorphic, x = c(0,17,44,66,90), col = 5, lty=2)
lines(P3_polymorphic, x = c(0,17,44,66,90), col = 6, lty=2)
legend(x = "topleft", legend=c("B1", "B2","B3","P1","P2","P3"), lty = c(1,1,1,2,2,2), col=c(1,2,3,4,5,6))

#now to print out a table 
Days <- matrix(ncol = 1,c("0","17","25","44","66","75","90") )
biofilm_polymorphic <- as.data.frame(cbind(Days,B1_polymorphic, B2_polymorphic, B3_polymorphic))
colnames(biofilm_polymorphic) <- c("Time point","B1", "B2","B3")

plank_days<- matrix(ncol = 1,c("0","17","44","66","90"))
planktonic_polymorphic <- as.data.frame(cbind(plank_days,P1_polymorphic, P2_polymorphic, P3_polymorphic))
colnames(planktonic_polymorphic) <- c("Time point","P1", "P2", "P3" )


#to combine all populations together, I need to introduce rows for the time points that were measured in biofilm but not in planktonic. If I don't do this, I won't be able to merge all of them into one table.
empty <- matrix(nrow=1, ncol = 4)
colnames(empty) <- colnames(planktonic_polymorphic)
planktonic_polymorphic <- rbind(planktonic_polymorphic[1:2,], empty, planktonic_polymorphic[3:4,], empty, planktonic_polymorphic[5,])

all_polymorphic <- cbind(biofilm_polymorphic, planktonic_polymorphic[,2:4])

#write out the table to a file so that I can use it later. 
write.csv(all_polymorphic, "/Users/katrina/Desktop/muller_v_0.3/Dynamics/all_polymorphic_mutations.csv")


#######
#now looking at the number of mutations that go extinct. This means that they start at a frequency of 0, rise to a detectable frequency, and then go back to 0. 
#######
#B1
B1_25_extinct <- B1[(B1$X17 != 0 & B1$X25 == 0 & B1$X44 == 0 & B1$X66 == 0 & B1$X75 == 0 & B1$X90 == 0),]
B1_44_extinct <- B1[( B1$X25 != 0 & B1$X44 == 0 & B1$X66 == 0 & B1$X75 == 0 & B1$X90 == 0),]
B1_66_extinct <- B1[(B1$X44 != 0 & B1$X66 == 0 & B1$X75 == 0 & B1$X90 == 0),]
B1_75_extinct <- B1[(B1$X66 != 0 & B1$X75 == 0 & B1$X90 == 0),]
B1_90_extinct <- B1[(B1$X75 != 0 & B1$X90 == 0),]

B1_extinct <- matrix(ncol = 2, data = c("Time point", "Number of mutations","0",0L, "17", 0L, "25", nrow(B1_25_extinct), "44", nrow(B1_44_extinct),"66",nrow(B1_66_extinct), "75",nrow(B1_75_extinct),"90",nrow(B1_90_extinct)), byrow=T)


#B2
B2_25_extinct <- B2[(B2$X17 != 0 & B2$X25 == 0 & B2$X44 == 0 & B2$X66 == 0 & B2$X75 == 0 & B2$X90 == 0),]
B2_44_extinct <- B2[(B2$X25 != 0 & B2$X44 == 0 & B2$X66 == 0 & B2$X75 == 0 & B2$X90 == 0),]
B2_66_extinct <- B2[(B2$X44 != 0 & B2$X66 == 0 & B2$X75 == 0 & B2$X90 == 0),]
B2_75_extinct <- B2[(B1$X66 != 0 & B2$X75 == 0 & B2$X90 == 0),]
B2_90_extinct <- B2[(B1$X75 != 0 & B2$X90 == 0),]

B2_extinct <- matrix(ncol = 2, data = c("Time point", "Number of mutations","0",0L, "17", 0L, "25", nrow(B2_25_extinct), "44", nrow(B2_44_extinct),"66",nrow(B2_66_extinct), "75",nrow(B2_75_extinct),"90",nrow(B2_90_extinct)), byrow=T)

#B3
B3_25_extinct <- B3[(B2$X17 != 0 & B3$X25 == 0 & B3$X44 == 0 & B3$X66 == 0 & B3$X75 == 0 & B3$X90 == 0),]
B3_44_extinct <- B3[(B3$X25 != 0 & B3$X44 == 0 & B3$X66 == 0 & B3$X75 == 0 & B3$X90 == 0),]
B3_66_extinct <- B3[(B3$X44 != 0 & B3$X66 == 0 & B3$X75 == 0 & B3$X90 == 0),]
B3_75_extinct <- B3[(B3$X66 != 0 & B3$X75 == 0 & B3$X90 == 0),]
B3_90_extinct <- B3[(B3$X75 != 0 & B3$X90 == 0),]

B3_extinct <- matrix(ncol = 2, data = c("Time point", "Number of mutations","0",0L, "17", 0L, "25", nrow(B3_25_extinct), "44", nrow(B3_44_extinct),"66",nrow(B3_66_extinct), "75",nrow(B3_75_extinct),"90",nrow(B3_90_extinct)), byrow=T)

#P1
P1_44_extinct <- P1[(P1$X17 != 0 & P1$X44 == 0 & P1$X66 == 0 & P1$X90 == 0),]
P1_66_extinct <- P1[(P1$X44 != 0 & P1$X66 == 0 & P1$X90 == 0),]
P1_90_extinct <- P1[(P1$X66 != 0 & P1$X90 == 0),]

P1_extinct <- matrix(ncol = 2, data = c("Time point", "Number of mutations","0",0L, "17", 0L, "44", nrow(P1_44_extinct),"66",nrow(P1_66_extinct),"90",nrow(P1_90_extinct)), byrow=T)


#P2
P2_44_extinct <- P2[(P2$X17 != 0 & P2$X44 == 0 & P2$X66 == 0 & P2$X90 == 0),]
P2_66_extinct <- P2[(P2$X44 != 0 & P2$X66 == 0 & P2$X90 == 0),]
P2_90_extinct <- P2[(P2$X66 != 0 & P2$X90 == 0),]

P2_extinct <- matrix(ncol = 2, data = c("Time point", "Number of mutations","0",0L, "17", 0L, "44", nrow(P2_44_extinct),"66",nrow(P2_66_extinct),"90",nrow(P2_90_extinct)), byrow=T)

#P3
P3_44_extinct <- P3[(P3$X17 != 0 & P3$X44 == 0 & P3$X66 == 0 & P3$X90 == 0),]
P3_66_extinct <- P3[(P3$X44 != 0 & P3$X66 == 0 & P3$X90 == 0),]
P3_90_extinct <- P3[(P3$X66 != 0 & P3$X90 == 0),]

P3_extinct <- matrix(ncol = 2, data = c("Time point", "Number of mutations","0",0L, "17", 0L, "44", nrow(P3_44_extinct),"66",nrow(P3_66_extinct),"90",nrow(P3_90_extinct)), byrow=T)


#plot all of the populations on one graph.
plot(B1_extinct, type="l", main= "Extinct mutations", xlab = "Days", ylab = "Number of mutations", ylim = c(0,50))
lines(B2_extinct, col = 2)
lines(B3_extinct, col = 3)
lines(P1_extinct, col = 4, lty = 2)
lines(P2_extinct, col = 5, lty = 2)
lines(P3_extinct, col = 6, lty = 2)

legend(x = "topleft", legend=c("B1", "B2","B3","P1","P2","P3"), lty = c(1,1,1,2,2,2), col=c(1,2,3,4,5,6))



#and put it all in a table 

biofilm_extinct<- cbind(B1_extinct, B2_extinct[,2], B3_extinct[,2])
colnames(biofilm_extinct) <- c("Time point", "B1", "B2","B3")
biofilm_extinct <- biofilm_extinct[-1,]

planktonic_extinct <- cbind(P1_extinct, P2_extinct[,2], P3_extinct[,2])
colnames(planktonic_extinct) <- c("Time point","P1", "P2", "P3" )
planktonic_extinct <- planktonic_extinct[-1,]

#to combine all populations together, I need to introduce rows for the time points that were measured in biofilm but not in planktonic.

empty <- matrix(nrow=1, ncol = 4)
planktonic_extinct <- rbind(planktonic_extinct[1:2,], empty, planktonic_extinct[3:4,], empty, planktonic_extinct[5,])

all_extinct <- cbind(biofilm_extinct, planktonic_extinct[,2:4])

#write out the table to a file so that I can use it later. 
write.csv(all_extinct, "/Users/katrina/Desktop/muller_v_0.3/Dynamics/all_extinct_mutations.csv")

####### 
#now to look at the number of new mutations that show up at each time point. 
#these are mutations that have never been seen in that population before. the sum of these should be the total number of rows in the imported data frame - will need to compare the sum with the number of rows. 
########

#B1
B1_17_new <- B1[(B1$X0 == 0 & B1$X17 != 0),]
B1_25_new <- B1[(B1$X0 == 0 & B1$X17 == 0 & B1$X25 != 0),]
B1_44_new <- B1[(B1$X0 == 0 & B1$X17 == 0 & B1$X25 == 0 & B1$X44 != 0),]
B1_66_new <- B1[(B1$X0 == 0 & B1$X17 == 0 & B1$X25 == 0 & B1$X44 == 0 & B1$X66 != 0),]
B1_75_new <- B1[(B1$X0 == 0 & B1$X17 == 0 & B1$X25 == 0 & B1$X44 == 0 & B1$X66 == 0 & B1$X75 != 0),]
B1_90_new <- B1[(B1$X0 == 0 & B1$X17 == 0 & B1$X25 == 0 & B1$X44 == 0 & B1$X66 == 0 & B1$X75 == 0 & B1$X90 != 0),]

B1_new <- matrix(ncol = 2, data = c("Time point", "Number of mutations","0",0L, "17", nrow(B1_17_new), "25", nrow(B1_25_new), "44", nrow(B1_44_new),"66",nrow(B1_66_new), "75",nrow(B1_75_new),"90",nrow(B1_90_new)), byrow=T)

#B2
B2_17_new <- B2[(B2$X0 == 0 & B2$X17 != 0),]
B2_25_new <- B2[(B2$X0 == 0 & B2$X17 == 0 & B2$X25 != 0),]
B2_44_new <- B2[(B2$X0 == 0 & B2$X17 == 0 & B2$X25 == 0 & B2$X44 != 0),]
B2_66_new <- B2[(B2$X0 == 0 & B2$X17 == 0 & B2$X25 == 0 & B2$X44 == 0 & B2$X66 != 0),]
B2_75_new <- B2[(B2$X0 == 0 & B2$X17 == 0 & B2$X25 == 0 & B2$X44 == 0 & B2$X66 == 0 & B2$X75 != 0),]
B2_90_new <- B2[(B2$X0 == 0 & B2$X17 == 0 & B2$X25 == 0 & B2$X44 == 0 & B2$X66 == 0 & B2$X75 == 0 & B2$X90 != 0),]

B2_new <- matrix(ncol = 2, data = c("Time point", "Number of mutations","0",0L, "17", nrow(B2_17_new), "25", nrow(B2_25_new), "44", nrow(B2_44_new),"66",nrow(B2_66_new), "75",nrow(B2_75_new),"90",nrow(B2_90_new)), byrow=T)


#B3
B3_17_new <- B3[(B3$X0 == 0 & B3$X17 != 0),]
B3_25_new <- B3[(B3$X0 == 0 & B3$X17 == 0 & B3$X25 != 0),]
B3_44_new <- B3[(B3$X0 == 0 & B3$X17 == 0 & B3$X25 == 0 & B3$X44 != 0),]
B3_66_new <- B3[(B3$X0 == 0 & B3$X17 == 0 & B3$X25 == 0 & B3$X44 == 0 & B3$X66 != 0),]
B3_75_new <- B3[(B3$X0 == 0 & B3$X17 == 0 & B3$X25 == 0 & B3$X44 == 0 & B3$X66 == 0 & B3$X75 != 0),]
B3_90_new <- B3[(B3$X0 == 0 & B3$X17 == 0 & B3$X25 == 0 & B3$X44 == 0 & B3$X66 == 0 & B3$X75 == 0 & B3$X90 != 0),]

B3_new <- matrix(ncol = 2, data = c("Time point", "Number of mutations","0",0L, "17", nrow(B3_17_new), "25", nrow(B3_25_new), "44", nrow(B3_44_new),"66",nrow(B3_66_new), "75",nrow(B3_75_new),"90",nrow(B3_90_new)), byrow=T)



#P1
P1_17_new <- P1[(P1$X0 == 0 & P1$X17 != 0),]
P1_44_new <- P1[(P1$X0 == 0 & P1$X17 == 0  & P1$X44 != 0),]
P1_66_new <- P1[(P1$X0 == 0 & P1$X17 == 0 & P1$X44 == 0 & P1$X66 != 0),]
P1_90_new <- P1[(P1$X0 == 0 & P1$X17 == 0 & P1$X44 == 0 & P1$X66 == 0 & P1$X90 != 0),]

P1_new <- matrix(ncol = 2, data = c("Time point", "Number of mutations","0",0L, "17", nrow(P1_17_new), "44", nrow(P1_44_new),"66",nrow(P1_66_new),"90",nrow(P1_90_new)), byrow=T)


#P2
P2_17_new <- P2[(P2$X0 == 0 & P2$X17 != 0),]
P2_44_new <- P2[(P2$X0 == 0 & P2$X17 == 0  & P2$X44 != 0),]
P2_66_new <- P2[(P2$X0 == 0 & P2$X17 == 0 & P2$X44 == 0 & P2$X66 != 0),]
P2_90_new <- P2[(P2$X0 == 0 & P2$X17 == 0 & P2$X44 == 0 & P2$X66 == 0 & P2$X90 != 0),]

P2_new <- matrix(ncol = 2, data = c("Time point", "Number of mutations","0",0L, "17", nrow(P2_17_new), "44", nrow(P2_44_new),"66",nrow(P2_66_new),"90",nrow(P2_90_new)), byrow=T)


#P3
P3_17_new <- P3[(P3$X0 == 0 & P3$X17 != 0),]
P3_44_new <- P3[(P3$X0 == 0 & P3$X17 == 0  & P3$X44 != 0),]
P3_66_new <- P3[(P3$X0 == 0 & P3$X17 == 0 & P3$X44 == 0 & P3$X66 != 0),]
P3_90_new <- P3[(P3$X0 == 0 & P3$X17 == 0 & P3$X44 == 0 & P3$X66 == 0 & P3$X90 != 0),]

P3_new <- matrix(ncol = 2, data = c("Time point", "Number of mutations","0",0L, "17", nrow(P3_17_new), "44", nrow(P3_44_new),"66",nrow(P3_66_new),"90",nrow(P3_90_new)), byrow=T)

#plot all of the data points.
plot(B1_new, xlab ="Days", ylab="Number of mutations", main = "New mutations", type ="l")
lines(B2_new, col = 2)
lines(B3_new, col = 3)
lines(P1_new, col = 4, lty = 2)
lines(P2_new, col = 5, lty = 2)
lines(P3_new, col = 6, lty = 2)
legend(x = "topright", legend=c("B1", "B2","B3","P1","P2","P3"), lty = c(1,1,1,2,2,2), col=c(1,2,3,4,5,6))

#need to print out a table of the numbers of mutations.






#and put it all in a table 

biofilm_new<- cbind(B1_new, B2_new[,2], B3_new[,2])
colnames(biofilm_new) <- c("Time point", "B1", "B2","B3")
biofilm_new <- biofilm_new[-1,]

planktonic_new <- cbind(P1_new, P2_new[,2], P3_new[,2])
colnames(planktonic_new) <- c("Time point","P1", "P2", "P3" )
planktonic_new <- planktonic_new[-1,]

#to combine all populations together, I need to introduce rows for the time points that were measured in biofilm but not in planktonic.

empty <- matrix(nrow=1, ncol = 4)
planktonic_new <- rbind(planktonic_new[1:2,], empty, planktonic_new[3:4,], empty, planktonic_new[5,])

all_new <- cbind(biofilm_new, planktonic_new[,2:4])

#write out the table to a file so that I can use it later. 
write.csv(all_new, "/Users/katrina/Desktop/muller_v_0.3/Dynamics/all_new_mutations.csv")


#########
#determine the "Total" number of mutations seen through time. This is the one plot that should always increase and never decrease.
########

#B1 - should end with 147
#sanity check to make sure that I am doing this correctly: 
sum(as.numeric(B1_new[2:8,2])) == nrow(B1) #True! So I did the new calculations correctly.
B1_total <- matrix(ncol = 2, byrow = T, c("Day", "Number of mutations", 
                                          "0", sum(as.numeric(B1_new[2:2,2])), 
                                          "17", sum(as.numeric(B1_new[2:3,2])), 
                                          "25", sum(as.numeric(B1_new[2:4,2])), 
                                          "44", sum(as.numeric(B1_new[2:5,2])),
                                          "66", sum(as.numeric(B1_new[2:6,2])),
                                          "75", sum(as.numeric(B1_new[2:7,2])),
                                          "90", sum(as.numeric(B1_new[2:8,2]))))
#B2
sum(as.numeric(B2_new[2:8,2])) == nrow(B2) #True! So I did the new calculations correctly.
B2_total <- matrix(ncol = 2, byrow = T, c("Day", "Number of mutations", 
                                          "0", sum(as.numeric(B2_new[2:2,2])), 
                                          "17", sum(as.numeric(B2_new[2:3,2])), 
                                          "25", sum(as.numeric(B2_new[2:4,2])), 
                                          "44", sum(as.numeric(B2_new[2:5,2])),
                                          "66", sum(as.numeric(B2_new[2:6,2])),
                                          "75", sum(as.numeric(B2_new[2:7,2])),
                                          "90", sum(as.numeric(B2_new[2:8,2]))))
#B3
sum(as.numeric(B3_new[2:8,2])) == nrow(B3) #True! So I did the new calculations correctly.
B3_total <- matrix(ncol = 2, byrow = T, c("Day", "Number of mutations", 
                                          "0", sum(as.numeric(B3_new[2:2,2])), 
                                          "17", sum(as.numeric(B3_new[2:3,2])), 
                                          "25", sum(as.numeric(B3_new[2:4,2])), 
                                          "44", sum(as.numeric(B3_new[2:5,2])),
                                          "66", sum(as.numeric(B3_new[2:6,2])),
                                          "75", sum(as.numeric(B3_new[2:7,2])),
                                          "90", sum(as.numeric(B3_new[2:8,2]))))

#P1 -
#sanity check to make sure that I am doing this correctly: 
sum(as.numeric(P1_new[2:6,2])) == nrow(P1) #True! So I did the new calculations correctly.
P1_total <- matrix(ncol = 2, byrow = T, c("Day", "Number of mutations", 
                                          "0", sum(as.numeric(P1_new[2:2,2])), 
                                          "17", sum(as.numeric(P1_new[2:3,2])), 
                                          "44", sum(as.numeric(P1_new[2:4,2])),
                                          "66", sum(as.numeric(P1_new[2:5,2])),
                                          "90", sum(as.numeric(P1_new[2:6,2]))))

#P2 -
#sanity check to make sure that I am doing this correctly: 
sum(as.numeric(P2_new[2:6,2])) == nrow(P2) #True! So I did the new calculations correctly.
P2_total <- matrix(ncol = 2, byrow = T, c("Day", "Number of mutations", 
                                          "0", sum(as.numeric(P2_new[2:2,2])), 
                                          "17", sum(as.numeric(P2_new[2:3,2])), 
                                          "44", sum(as.numeric(P2_new[2:4,2])),
                                          "66", sum(as.numeric(P2_new[2:5,2])),
                                          "90", sum(as.numeric(P2_new[2:6,2]))))

#P3 -
#sanity check to make sure that I am doing this correctly: 
sum(as.numeric(P3_new[2:6,2])) == nrow(P3) #True! So I did the new calculations correctly.
P3_total <- matrix(ncol = 2, byrow = T, c("Day", "Number of mutations", 
                                          "0", sum(as.numeric(P3_new[2:2,2])), 
                                          "17", sum(as.numeric(P3_new[2:3,2])), 
                                          "44", sum(as.numeric(P3_new[2:4,2])),
                                          "66", sum(as.numeric(P3_new[2:5,2])),
                                          "90", sum(as.numeric(P3_new[2:6,2]))))



plot(B1_total,xlab="Days", ylab="Number of mutations", main ="Total number of mutations", type="l")
lines(B2_total, col = 2)
lines(B3_total, col = 3)
lines(P1_total, col = 4, lty = 2)
lines(P2_total, col = 5, lty = 2)
lines(P3_total, col = 6, lty = 2)
legend(x = "topleft", legend=c("B1", "B2","B3","P1","P2","P3"), lty = c(1,1,1,2,2,2), col=c(1,2,3,4,5,6))


#and put it all in a table 

biofilm_total<- cbind(B1_total, B2_total[,2], B3_total[,2])
colnames(biofilm_total) <- c("Time point", "B1", "B2","B3")
biofilm_total <- biofilm_total[-1,]

planktonic_total <- cbind(P1_total, P2_total[,2], P3_total[,2])
colnames(planktonic_total) <- c("Time point","P1", "P2", "P3" )
planktonic_total<- planktonic_total[-1,]

#to combine all populations together, I need to introduce rows for the time points that were measured in biofilm but not in planktonic.

empty <- matrix(nrow=1, ncol = 4)
planktonic_total <- rbind(planktonic_total[1:2,], empty, planktonic_total[3:4,], empty, planktonic_total[5,])

all_total <- cbind(biofilm_total, planktonic_total[,2:4])

#write out the table to a file so that I can use it later. 
write.csv(all_total, "/Users/katrina/Desktop/muller_v_0.3/Dynamics/all_total_mutations.csv")



#######
#plot them all together in one big image so that they can be compared 
#######
par(mfrow=c(2,3))

#total number of mutations
plot(B1_total,xlab="Days", ylab="Number of mutations", type="l", ylim = c(0,200), cex.lab = 1.5, lwd = 2)
lines(B2_total, col = 2, lwd = 2)
lines(B3_total, col = 3, lwd = 2)
lines(P1_total, col = 4, lty = 2, lwd = 2)
lines(P2_total, col = 5, lty = 2, lwd = 2)
lines(P3_total, col = 6, lty = 2, lwd = 2)
title("A.", adj = 0, cex.main = 2)

#all mutations above 0 at each time point 
plot(B1_numbers, type = "l", ylim = c(0,200), lwd = 2, ylab = "", xlab = "")
lines(B2_numbers, col = 2, lwd = 2)
lines(B3_numbers, col = 3, lwd = 2)
lines(P1_numbers, col = 4, lty = 2, lwd = 2)
lines(P2_numbers, col = 5, lty = 2, lwd = 2)
lines(P3_numbers, col = 6, lty = 2, lwd = 2)
title("B.", adj = 0, cex.main = 2)

#fixed mutations
plot(B1_fixed, typ="l", xlim=c(0,90), ylim = c(0,200), lwd = 2, ylab = "", xlab = "")
lines(B2_fixed, typ="l", col=2, lwd = 2)
lines(B3_fixed, typ="l", col=3, lwd = 2)
lines(P1_fixed, typ="l", col=4, lty=2, lwd = 2)
lines(P2_fixed, typ="l", col=5, lty=2, lwd = 2) #is still on there but the P2 and P3 fixation rates are identical.
lines(P3_fixed, typ="l", col=6, lty=2, lwd = 2) 
legend(x = "topright", legend=c("B1", "B2","B3","P1","P2","P3"), lty = c(1,1,1,2,2,2), col=c(1,2,3,4,5,6), lwd = 2)
title("C.", adj = 0, cex.main = 2)


#polymorphic mutations 
plot(B1_polymorphic, x=c(0,17,25,44,66,75,90), type = "l", ylim = c(0,200), lwd = 2, ylab = "", xlab = "")
lines(B2_polymorphic, x=c(0,17,25,44,66,75,90), col = 2, lwd = 2)
lines(B3_polymorphic, x=c(0,17,25,44,66,75,90), col = 3, lwd = 2)
lines(P1_polymorphic, x = c(0,17,44,66,90), col = 4, lty=2, lwd = 2)
lines(P2_polymorphic, x = c(0,17,44,66,90), col = 5, lty=2, lwd = 2)
lines(P3_polymorphic, x = c(0,17,44,66,90), col = 6, lty=2, lwd = 2)
title("D.", adj = 0, cex.main = 2)

#extinct mutations 
plot(B1_extinct, type="l", ylim = c(0,200), ylab = "", xlab = "",lwd = 2)
lines(B2_extinct, col = 2, lwd = 2)
lines(B3_extinct, col = 3, lwd = 2)
lines(P1_extinct, col = 4, lty = 2, lwd = 2)
lines(P2_extinct, col = 5, lty = 2, lwd = 2)
lines(P3_extinct, col = 6, lty = 2, lwd = 2)
title("E.", adj = 0, cex.main = 2)

#new mutations
plot(B1_new, xlab ="", ylab="", type ="l", ylim = c(0,200),lwd = 2)
lines(B2_new, col = 2,lwd = 2)
lines(B3_new, col = 3,lwd = 2)
lines(P1_new, col = 4, lty = 2,lwd = 2)
lines(P2_new, col = 5, lty = 2,lwd = 2)
lines(P3_new, col = 6, lty = 2,lwd = 2)
title("F.", adj = 0, cex.main = 2)

######
#now I need to look at these dynamics across environments
#B1, B2, B3 vs P1, P2, P3
######
#total number of mutations by environemnt - only reporting averages of the 3 replicates for each one at the present moment

#initiate the matrix
by_environment_total <- matrix(ncol = 3, nrow = 7)
#add in the days 
by_environment_total[,1] <- c("0", "17","25","44","66","75","90")
#calculate averages for biofilm replicates at each day
by_environment_total[1,2] <- mean(x = c(as.numeric(B1_total[2,2]),as.numeric(B2_total[2,2]),as.numeric(B3_total[2,2])))
by_environment_total[2,2] <- mean(x = c(as.numeric(B1_total[3,2]),as.numeric(B2_total[3,2]),as.numeric(B3_total[3,2])))
by_environment_total[3,2] <- mean(x = c(as.numeric(B1_total[4,2]),as.numeric(B2_total[4,2]),as.numeric(B3_total[4,2])))
by_environment_total[4,2] <- mean(x = c(as.numeric(B1_total[5,2]),as.numeric(B2_total[5,2]),as.numeric(B3_total[5,2])))
by_environment_total[5,2] <- mean(x = c(as.numeric(B1_total[6,2]),as.numeric(B2_total[6,2]),as.numeric(B3_total[6,2])))
by_environment_total[6,2] <- mean(x = c(as.numeric(B1_total[7,2]),as.numeric(B2_total[7,2]),as.numeric(B3_total[7,2])))
by_environment_total[7,2] <- mean(x = c(as.numeric(B1_total[8,2]),as.numeric(B2_total[8,2]),as.numeric(B3_total[8,2])))
#calculate averages of planktonic replicates at each day
by_environment_total[1,3] <- mean(x = c(as.numeric(P1_total[2,2]),as.numeric(P2_total[2,2]),as.numeric(P3_total[2,2])))
by_environment_total[2,3] <- mean(x = c(as.numeric(P1_total[3,2]),as.numeric(P2_total[3,2]),as.numeric(P3_total[3,2])))
by_environment_total[4,3] <- mean(x = c(as.numeric(P1_total[4,2]),as.numeric(P2_total[4,2]),as.numeric(P3_total[4,2])))
by_environment_total[5,3] <- mean(x = c(as.numeric(P1_total[5,2]),as.numeric(P2_total[5,2]),as.numeric(P3_total[5,2])))
by_environment_total[7,3] <- mean(x = c(as.numeric(P1_total[6,2]),as.numeric(P2_total[6,2]),as.numeric(P3_total[6,2])))

#plot total numbers of mutations by environment
plot(x = by_environment_total[,1], y = by_environment_total[,2], type = "l", xlab = "Day", ylab="Number of mutations", main = "Total number of mutations by environment", lwd = 3, ylim=c(0,200))
lines(x = by_environment_total[c(1,2,4,5,7),1], y = by_environment_total[c(1,2,4,5,7),3], col = 2, lwd = 3)
legend(x= "topleft", legend=c("Biofilm", "Planktonic"), col= c(1,2), lty=1, lwd=3)



####
#now do by environment averages of mutations above frequency of 0 at each time point.

by_environment_numbers <- matrix(ncol = 3, nrow = 7)
#add in the days 
by_environment_numbers[,1] <- c("0", "17","25","44","66","75","90")
#calculate averages for biofilm replicates at each day
by_environment_numbers[1,2] <- mean(x = c(as.numeric(B1_numbers[2,2]),as.numeric(B2_numbers[2,2]),as.numeric(B3_numbers[2,2])))
by_environment_numbers[2,2] <- mean(x = c(as.numeric(B1_numbers[3,2]),as.numeric(B2_numbers[3,2]),as.numeric(B3_numbers[3,2])))
by_environment_numbers[3,2] <- mean(x = c(as.numeric(B1_numbers[4,2]),as.numeric(B2_numbers[4,2]),as.numeric(B3_numbers[4,2])))
by_environment_numbers[4,2] <- mean(x = c(as.numeric(B1_numbers[5,2]),as.numeric(B2_numbers[5,2]),as.numeric(B3_numbers[5,2])))
by_environment_numbers[5,2] <- mean(x = c(as.numeric(B1_numbers[6,2]),as.numeric(B2_numbers[6,2]),as.numeric(B3_numbers[6,2])))
by_environment_numbers[6,2] <- mean(x = c(as.numeric(B1_numbers[7,2]),as.numeric(B2_numbers[7,2]),as.numeric(B3_numbers[7,2])))
by_environment_numbers[7,2] <- mean(x = c(as.numeric(B1_numbers[8,2]),as.numeric(B2_numbers[8,2]),as.numeric(B3_numbers[8,2])))
#calculate averages of planktonic replicates at each day
by_environment_numbers[1,3] <- mean(x = c(as.numeric(P1_numbers[2,2]),as.numeric(P2_numbers[2,2]),as.numeric(P3_numbers[2,2])))
by_environment_numbers[2,3] <- mean(x = c(as.numeric(P1_numbers[3,2]),as.numeric(P2_numbers[3,2]),as.numeric(P3_numbers[3,2])))
by_environment_numbers[4,3] <- mean(x = c(as.numeric(P1_numbers[4,2]),as.numeric(P2_numbers[4,2]),as.numeric(P3_numbers[4,2])))
by_environment_numbers[5,3] <- mean(x = c(as.numeric(P1_numbers[5,2]),as.numeric(P2_numbers[5,2]),as.numeric(P3_numbers[5,2])))
by_environment_numbers[7,3] <- mean(x = c(as.numeric(P1_numbers[6,2]),as.numeric(P2_numbers[6,2]),as.numeric(P3_numbers[6,2])))

#plot total numbers of mutations by environment
plot(x = by_environment_numbers[,1], y = by_environment_numbers[,2], type = "l", xlab = "Day", ylab="Number of mutations", main = "Detected mutations", lwd = 3, ylim = c(0,200))
lines(x = by_environment_numbers[c(1,2,4,5,7),1], y = by_environment_numbers[c(1,2,4,5,7),3], col = 2, lwd = 3)
legend(x= "topleft", legend=c("Biofilm", "Planktonic"), col= c(1,2), lty=1, lwd=3)


### 
#now do environmnet averages for fixed mutations

by_environment_fixed <- matrix(ncol = 3, nrow = 7)
#add in the days 
by_environment_fixed[,1] <- c("0", "17","25","44","66","75","90")
#calculate averages for biofilm replicates at each day
by_environment_fixed[1,2] <- mean(x = c(as.numeric(B1_fixed[2,2]),as.numeric(B2_fixed[2,2]),as.numeric(B3_fixed[2,2])))
by_environment_fixed[2,2] <- mean(x = c(as.numeric(B1_fixed[3,2]),as.numeric(B2_fixed[3,2]),as.numeric(B3_fixed[3,2])))
by_environment_fixed[3,2] <- mean(x = c(as.numeric(B1_fixed[4,2]),as.numeric(B2_fixed[4,2]),as.numeric(B3_fixed[4,2])))
by_environment_fixed[4,2] <- mean(x = c(as.numeric(B1_fixed[5,2]),as.numeric(B2_fixed[5,2]),as.numeric(B3_fixed[5,2])))
by_environment_fixed[5,2] <- mean(x = c(as.numeric(B1_fixed[6,2]),as.numeric(B2_fixed[6,2]),as.numeric(B3_fixed[6,2])))
by_environment_fixed[6,2] <- mean(x = c(as.numeric(B1_fixed[7,2]),as.numeric(B2_fixed[7,2]),as.numeric(B3_fixed[7,2])))
by_environment_fixed[7,2] <- mean(x = c(as.numeric(B1_fixed[8,2]),as.numeric(B2_fixed[8,2]),as.numeric(B3_fixed[8,2])))
#calculate averages of planktonic replicates at each day
by_environment_fixed[1,3] <- mean(x = c(as.numeric(P1_fixed[2,2]),as.numeric(P2_fixed[2,2]),as.numeric(P3_fixed[2,2])))
by_environment_fixed[2,3] <- mean(x = c(as.numeric(P1_fixed[3,2]),as.numeric(P2_fixed[3,2]),as.numeric(P3_fixed[3,2])))
by_environment_fixed[4,3] <- mean(x = c(as.numeric(P1_fixed[4,2]),as.numeric(P2_fixed[4,2]),as.numeric(P3_fixed[4,2])))
by_environment_fixed[5,3] <- mean(x = c(as.numeric(P1_fixed[5,2]),as.numeric(P2_fixed[5,2]),as.numeric(P3_fixed[5,2])))
by_environment_fixed[7,3] <- mean(x = c(as.numeric(P1_fixed[6,2]),as.numeric(P2_fixed[6,2]),as.numeric(P3_fixed[6,2])))

#plot fixed mutations by environment
plot(x = by_environment_fixed[,1], y = by_environment_fixed[,2], type = "l", xlab = "Day", ylab="Number of mutations", main = "fixed mutations", lwd = 3, ylim = c(0,200))
lines(x = by_environment_fixed[c(1,2,4,5,7),1], y = by_environment_fixed[c(1,2,4,5,7),3], col = 2, lwd = 3)
legend(x= "topleft", legend=c("Biofilm", "Planktonic"), col= c(1,2), lty=1, lwd=3)


#### 
#now to look at polymorphic mutations by environment

by_environment_polymorphic<- matrix(ncol = 3, nrow = 7)
#add in the days 
by_environment_polymorphic[,1] <- c("0", "17","25","44","66","75","90")
#calculate averages for biofilm replicates at each day
by_environment_polymorphic[1,2] <- mean(x = c(as.numeric(B1_polymorphic[1,1]),as.numeric(B2_polymorphic[1,1]),as.numeric(B3_polymorphic[1,1])))
by_environment_polymorphic[2,2] <- mean(x = c(as.numeric(B1_polymorphic[2,1]),as.numeric(B2_polymorphic[2,1]),as.numeric(B3_polymorphic[2,1])))
by_environment_polymorphic[3,2] <- mean(x = c(as.numeric(B1_polymorphic[3,1]),as.numeric(B2_polymorphic[3,1]),as.numeric(B3_polymorphic[3,1])))
by_environment_polymorphic[4,2] <- mean(x = c(as.numeric(B1_polymorphic[4,1]),as.numeric(B2_polymorphic[4,1]),as.numeric(B3_polymorphic[4,1])))
by_environment_polymorphic[5,2] <- mean(x = c(as.numeric(B1_polymorphic[5,1]),as.numeric(B2_polymorphic[5,1]),as.numeric(B3_polymorphic[5,1])))
by_environment_polymorphic[6,2] <- mean(x = c(as.numeric(B1_polymorphic[6,1]),as.numeric(B2_polymorphic[6,1]),as.numeric(B3_polymorphic[6,1])))
by_environment_polymorphic[7,2] <- mean(x = c(as.numeric(B1_polymorphic[7,1]),as.numeric(B2_polymorphic[7,1]),as.numeric(B3_polymorphic[7,1])))
#calculate averages of planktonic replicates at each day
by_environment_polymorphic[1,3] <- mean(x = c(as.numeric(P1_polymorphic[1,1]),as.numeric(P2_polymorphic[1,1]),as.numeric(P3_polymorphic[1,1])))
by_environment_polymorphic[2,3] <- mean(x = c(as.numeric(P1_polymorphic[2,1]),as.numeric(P2_polymorphic[2,1]),as.numeric(P3_polymorphic[2,1])))
by_environment_polymorphic[4,3] <- mean(x = c(as.numeric(P1_polymorphic[3,1]),as.numeric(P2_polymorphic[3,1]),as.numeric(P3_polymorphic[3,1])))
by_environment_polymorphic[5,3] <- mean(x = c(as.numeric(P1_polymorphic[4,1]),as.numeric(P2_polymorphic[4,1]),as.numeric(P3_polymorphic[4,1])))
by_environment_polymorphic[7,3] <- mean(x = c(as.numeric(P1_polymorphic[5,1]),as.numeric(P2_polymorphic[5,1]),as.numeric(P3_polymorphic[5,1])))

#plot polymorphic mutations by environment
plot(x = by_environment_polymorphic[,1], y = by_environment_polymorphic[,2], type = "l", xlab = "Day", ylab="Number of mutations", main = "polymorphic mutations", lwd = 3, ylim = c(0,200))
lines(x = by_environment_polymorphic[c(1,2,4,5,7),1], y = by_environment_polymorphic[c(1,2,4,5,7),3], col = 2, lwd = 3)
legend(x= "topleft", legend=c("Biofilm", "Planktonic"), col= c(1,2), lty=1, lwd=3)



### 
#now do environmnet averages for extinct mutations

by_environment_extinct <- matrix(ncol = 3, nrow = 7)
#add in the days 
by_environment_extinct[,1] <- c("0", "17","25","44","66","75","90")
#calculate averages for biofilm replicates at each day
by_environment_extinct[1,2] <- mean(x = c(as.numeric(B1_extinct[2,2]),as.numeric(B2_extinct[2,2]),as.numeric(B3_extinct[2,2])))
by_environment_extinct[2,2] <- mean(x = c(as.numeric(B1_extinct[3,2]),as.numeric(B2_extinct[3,2]),as.numeric(B3_extinct[3,2])))
by_environment_extinct[3,2] <- mean(x = c(as.numeric(B1_extinct[4,2]),as.numeric(B2_extinct[4,2]),as.numeric(B3_extinct[4,2])))
by_environment_extinct[4,2] <- mean(x = c(as.numeric(B1_extinct[5,2]),as.numeric(B2_extinct[5,2]),as.numeric(B3_extinct[5,2])))
by_environment_extinct[5,2] <- mean(x = c(as.numeric(B1_extinct[6,2]),as.numeric(B2_extinct[6,2]),as.numeric(B3_extinct[6,2])))
by_environment_extinct[6,2] <- mean(x = c(as.numeric(B1_extinct[7,2]),as.numeric(B2_extinct[7,2]),as.numeric(B3_extinct[7,2])))
by_environment_extinct[7,2] <- mean(x = c(as.numeric(B1_extinct[8,2]),as.numeric(B2_extinct[8,2]),as.numeric(B3_extinct[8,2])))
#calculate averages of planktonic replicates at each day
by_environment_extinct[1,3] <- mean(x = c(as.numeric(P1_extinct[2,2]),as.numeric(P2_extinct[2,2]),as.numeric(P3_extinct[2,2])))
by_environment_extinct[2,3] <- mean(x = c(as.numeric(P1_extinct[3,2]),as.numeric(P2_extinct[3,2]),as.numeric(P3_extinct[3,2])))
by_environment_extinct[4,3] <- mean(x = c(as.numeric(P1_extinct[4,2]),as.numeric(P2_extinct[4,2]),as.numeric(P3_extinct[4,2])))
by_environment_extinct[5,3] <- mean(x = c(as.numeric(P1_extinct[5,2]),as.numeric(P2_extinct[5,2]),as.numeric(P3_extinct[5,2])))
by_environment_extinct[7,3] <- mean(x = c(as.numeric(P1_extinct[6,2]),as.numeric(P2_extinct[6,2]),as.numeric(P3_extinct[6,2])))

#plot extinct mutations by environment
plot(x = by_environment_extinct[,1], y = by_environment_extinct[,2], type = "l", xlab = "Day", ylab="Number of mutations", main = "Extinct mutations", lwd = 3, ylim = c(0,200))
lines(x = by_environment_extinct[c(1,2,4,5,7),1], y = by_environment_extinct[c(1,2,4,5,7),3], col = 2, lwd = 3)
legend(x= "topleft", legend=c("Biofilm", "Planktonic"), col= c(1,2), lty=1, lwd=3)

###
#last but not least, need to do all of the new mutations by environment

by_environment_new <- matrix(ncol = 3, nrow = 7)
#add in the days 
by_environment_new[,1] <- c("0", "17","25","44","66","75","90")
#calculate averages for biofilm replicates at each day
by_environment_new[1,2] <- mean(x = c(as.numeric(B1_new[2,2]),as.numeric(B2_new[2,2]),as.numeric(B3_new[2,2])))
by_environment_new[2,2] <- mean(x = c(as.numeric(B1_new[3,2]),as.numeric(B2_new[3,2]),as.numeric(B3_new[3,2])))
by_environment_new[3,2] <- mean(x = c(as.numeric(B1_new[4,2]),as.numeric(B2_new[4,2]),as.numeric(B3_new[4,2])))
by_environment_new[4,2] <- mean(x = c(as.numeric(B1_new[5,2]),as.numeric(B2_new[5,2]),as.numeric(B3_new[5,2])))
by_environment_new[5,2] <- mean(x = c(as.numeric(B1_new[6,2]),as.numeric(B2_new[6,2]),as.numeric(B3_new[6,2])))
by_environment_new[6,2] <- mean(x = c(as.numeric(B1_new[7,2]),as.numeric(B2_new[7,2]),as.numeric(B3_new[7,2])))
by_environment_new[7,2] <- mean(x = c(as.numeric(B1_new[8,2]),as.numeric(B2_new[8,2]),as.numeric(B3_new[8,2])))
#calculate averages of planktonic replicates at each day
by_environment_new[1,3] <- mean(x = c(as.numeric(P1_new[2,2]),as.numeric(P2_new[2,2]),as.numeric(P3_new[2,2])))
by_environment_new[2,3] <- mean(x = c(as.numeric(P1_new[3,2]),as.numeric(P2_new[3,2]),as.numeric(P3_new[3,2])))
by_environment_new[4,3] <- mean(x = c(as.numeric(P1_new[4,2]),as.numeric(P2_new[4,2]),as.numeric(P3_new[4,2])))
by_environment_new[5,3] <- mean(x = c(as.numeric(P1_new[5,2]),as.numeric(P2_new[5,2]),as.numeric(P3_new[5,2])))
by_environment_new[7,3] <- mean(x = c(as.numeric(P1_new[6,2]),as.numeric(P2_new[6,2]),as.numeric(P3_new[6,2])))

#plot extinct mutations by environment
plot(x = by_environment_new[,1], y = by_environment_new[,2], type = "l", xlab = "Day", ylab="Number of mutations", main = "New mutations", lwd = 3, ylim = c(0,200))
lines(x = by_environment_new[c(1,2,4,5,7),1], y = by_environment_new[c(1,2,4,5,7),3], col = 2, lwd = 3)
legend(x= "topleft", legend=c("Biofilm", "Planktonic"), col= c(1,2), lty=1, lwd=3)


##plot all by environment in the same pane.
par(mfrow=c(2,3))

#plot total numbers of mutations by environment
plot(x = by_environment_total[,1], y = by_environment_total[,2], type = "l", xlab = "Day", ylab="Number of mutations", ylim=c(0,200), lwd=2, cex.lab = 1.5, col = 4)
lines(x = by_environment_total[c(1,2,4,5,7),1], y = by_environment_total[c(1,2,4,5,7),3], col = 2, lwd = 2)
title("A.", adj = 0, cex.main = 2)

#plot number of detected mutations by environment
plot(x = by_environment_numbers[,1], y = by_environment_numbers[,2], type = "l", xlab = "", ylab="",  lwd = 2, ylim = c(0,200), col = 4)
lines(x = by_environment_numbers[c(1,2,4,5,7),1], y = by_environment_numbers[c(1,2,4,5,7),3], col = 2, lwd = 2)
title("B.", adj = 0, cex.main = 2)

#plot fixed mutations by environment
plot(x = by_environment_fixed[,1], y = by_environment_fixed[,2], type = "l", xlab = "", ylab="", lwd = 2, ylim = c(0,200), col = 4)
lines(x = by_environment_fixed[c(1,2,4,5,7),1], y = by_environment_fixed[c(1,2,4,5,7),3], col = 2, lwd = 2)
legend(x= "topright", legend=c("Biofilm", "Planktonic"), col= c(4,2), lty=1, lwd=2)
title("C.", adj = 0, cex.main = 2)

#plot polymorphic mutations by environment
plot(x = by_environment_polymorphic[,1], y = by_environment_polymorphic[,2], type = "l", xlab = "", ylab="", lwd = 2, ylim = c(0,200), col = 4)
lines(x = by_environment_polymorphic[c(1,2,4,5,7),1], y = by_environment_polymorphic[c(1,2,4,5,7),3], col = 2, lwd = 2)
title("D.", adj = 0, cex.main = 2)

#plot extinct mutations by environment
plot(x = by_environment_extinct[,1], y = by_environment_extinct[,2], type = "l", xlab = "", ylab="" , lwd = 2, ylim = c(0,200), col = 4)
lines(x = by_environment_extinct[c(1,2,4,5,7),1], y = by_environment_extinct[c(1,2,4,5,7),3], col = 2, lwd = 2)
title("E.", adj = 0, cex.main = 2)

#plot extinct mutations by environment
plot(x = by_environment_new[,1], y = by_environment_new[,2], type = "l", xlab = "", ylab="", lwd = 2, ylim = c(0,200), col = 4)
lines(x = by_environment_new[c(1,2,4,5,7),1], y = by_environment_new[c(1,2,4,5,7),3], col = 2, lwd = 2)
title("F.", adj = 0, cex.main = 2)

######
# now to break into mutator vs non mutator.
# B1, B2, B3, P3 vs P1, P2
######

#first for total number of mutations by allele
#initiate the matrix
mut_nonmut_total <- matrix(ncol = 3, nrow = 5)
#add in the days 
mut_nonmut_total[,1] <- c("0", "17","44","66","90")
#calculate averages for mutator replicates at each point measured - can only do 4 time points since we are combinging plank and biofilm samples
mut_nonmut_total[1,2] <- mean(x = c(as.numeric(B1_total[2,2]),as.numeric(B2_total[2,2]),as.numeric(B3_total[2,2]),as.numeric(P3_total[2,2])))
mut_nonmut_total[2,2] <- mean(x = c(as.numeric(B1_total[3,2]),as.numeric(B2_total[3,2]),as.numeric(B3_total[3,2]),as.numeric(P3_total[3,2])))
mut_nonmut_total[3,2] <- mean(x = c(as.numeric(B1_total[5,2]),as.numeric(B2_total[5,2]),as.numeric(B3_total[5,2]),as.numeric(P3_total[4,2])))
mut_nonmut_total[4,2] <- mean(x = c(as.numeric(B1_total[6,2]),as.numeric(B2_total[6,2]),as.numeric(B3_total[6,2]),as.numeric(P3_total[5,2])))
mut_nonmut_total[5,2] <- mean(x = c(as.numeric(B1_total[8,2]),as.numeric(B2_total[8,2]),as.numeric(B3_total[8,2]),as.numeric(P3_total[6,2])))
#calculate averages of planktonic replicates at each day
mut_nonmut_total[1,3] <- mean(x = c(as.numeric(P1_total[2,2]),as.numeric(P2_total[2,2])))
mut_nonmut_total[2,3] <- mean(x = c(as.numeric(P1_total[3,2]),as.numeric(P2_total[3,2])))
mut_nonmut_total[3,3] <- mean(x = c(as.numeric(P1_total[4,2]),as.numeric(P2_total[4,2])))
mut_nonmut_total[4,3] <- mean(x = c(as.numeric(P1_total[5,2]),as.numeric(P2_total[5,2])))
mut_nonmut_total[5,3] <- mean(x = c(as.numeric(P1_total[6,2]),as.numeric(P2_total[6,2])))

#plot total numbers of mutations by environment
plot(x = mut_nonmut_total[,1], y = mut_nonmut_total[,2], type = "l", xlab = "Day", ylab="Number of mutations", main = "Total mutations", lwd = 3, ylim=c(0,200), col = "magenta")
lines(x = mut_nonmut_total[,1], y = mut_nonmut_total[,3], col = "purple", lwd = 3)
legend(x= "topleft", legend=c("Mutator", "Non-Mutator"), col= c("magenta","purple"), lty=1, lwd=3)


###
#now to look at all of the detected mutations at each time point 
#initiate the matrix
mut_nonmut_numbers <- matrix(ncol = 3, nrow = 5)
#add in the days 
mut_nonmut_numbers[,1] <- c("0", "17","44","66","90")
#calculate averages for mutator replicates at each point measured - can only do 4 time points since we are combinging plank and biofilm samples
mut_nonmut_numbers[1,2] <- mean(x = c(as.numeric(B1_numbers[2,2]),as.numeric(B2_numbers[2,2]),as.numeric(B3_numbers[2,2]),as.numeric(P3_numbers[2,2])))
mut_nonmut_numbers[2,2] <- mean(x = c(as.numeric(B1_numbers[3,2]),as.numeric(B2_numbers[3,2]),as.numeric(B3_numbers[3,2]),as.numeric(P3_numbers[3,2])))
mut_nonmut_numbers[3,2] <- mean(x = c(as.numeric(B1_numbers[5,2]),as.numeric(B2_numbers[5,2]),as.numeric(B3_numbers[5,2]),as.numeric(P3_numbers[4,2])))
mut_nonmut_numbers[4,2] <- mean(x = c(as.numeric(B1_numbers[6,2]),as.numeric(B2_numbers[6,2]),as.numeric(B3_numbers[6,2]),as.numeric(P3_numbers[5,2])))
mut_nonmut_numbers[5,2] <- mean(x = c(as.numeric(B1_numbers[8,2]),as.numeric(B2_numbers[8,2]),as.numeric(B3_numbers[8,2]),as.numeric(P3_numbers[6,2])))
#calculate averages of planktonic replicates at each day
mut_nonmut_numbers[1,3] <- mean(x = c(as.numeric(P1_numbers[2,2]),as.numeric(P2_numbers[2,2])))
mut_nonmut_numbers[2,3] <- mean(x = c(as.numeric(P1_numbers[3,2]),as.numeric(P2_numbers[3,2])))
mut_nonmut_numbers[3,3] <- mean(x = c(as.numeric(P1_numbers[4,2]),as.numeric(P2_numbers[4,2])))
mut_nonmut_numbers[4,3] <- mean(x = c(as.numeric(P1_numbers[5,2]),as.numeric(P2_numbers[5,2])))
mut_nonmut_numbers[5,3] <- mean(x = c(as.numeric(P1_numbers[6,2]),as.numeric(P2_numbers[6,2])))

#plot total numbers of mutations by environment
plot(x = mut_nonmut_numbers[,1], y = mut_nonmut_numbers[,2], type = "l", xlab = "Day", ylab="Number of mutations", main = "Detected mutations", lwd = 3, ylim=c(0,200), col = "magenta")
lines(x = mut_nonmut_numbers[,1], y = mut_nonmut_numbers[,3], lwd = 3, col = "purple")
legend(x= "topleft", legend=c("Mutator", "Non-mutator"), col= c("magenta","purple"), lty=1, lwd=3)


#### 
#looking at the fixed mutations by allele
#initiate the matrix
mut_nonmut_fixed <- matrix(ncol = 3, nrow = 5)
#add in the days 
mut_nonmut_fixed[,1] <- c("0", "17","44","66","90")
#calculate averages for mutator replicates at each point measured - can only do 4 time points since we are combinging plank and biofilm samples
mut_nonmut_fixed[1,2] <- mean(x = c(as.numeric(B1_fixed[2,2]),as.numeric(B2_fixed[2,2]),as.numeric(B3_fixed[2,2]),as.numeric(P3_fixed[2,2])))
mut_nonmut_fixed[2,2] <- mean(x = c(as.numeric(B1_fixed[3,2]),as.numeric(B2_fixed[3,2]),as.numeric(B3_fixed[3,2]),as.numeric(P3_fixed[3,2])))
mut_nonmut_fixed[3,2] <- mean(x = c(as.numeric(B1_fixed[5,2]),as.numeric(B2_fixed[5,2]),as.numeric(B3_fixed[5,2]),as.numeric(P3_fixed[4,2])))
mut_nonmut_fixed[4,2] <- mean(x = c(as.numeric(B1_fixed[6,2]),as.numeric(B2_fixed[6,2]),as.numeric(B3_fixed[6,2]),as.numeric(P3_fixed[5,2])))
mut_nonmut_fixed[5,2] <- mean(x = c(as.numeric(B1_fixed[8,2]),as.numeric(B2_fixed[8,2]),as.numeric(B3_fixed[8,2]),as.numeric(P3_fixed[6,2])))
#calculate averages of planktonic replicates at each day
mut_nonmut_fixed[1,3] <- mean(x = c(as.numeric(P1_fixed[2,2]),as.numeric(P2_fixed[2,2])))
mut_nonmut_fixed[2,3] <- mean(x = c(as.numeric(P1_fixed[3,2]),as.numeric(P2_fixed[3,2])))
mut_nonmut_fixed[3,3] <- mean(x = c(as.numeric(P1_fixed[4,2]),as.numeric(P2_fixed[4,2])))
mut_nonmut_fixed[4,3] <- mean(x = c(as.numeric(P1_fixed[5,2]),as.numeric(P2_fixed[5,2])))
mut_nonmut_fixed[5,3] <- mean(x = c(as.numeric(P1_fixed[6,2]),as.numeric(P2_fixed[6,2])))

#plot fixed mutations by environment
plot(x = mut_nonmut_fixed[,1], y = mut_nonmut_fixed[,2], type = "l", xlab = "Day", ylab="Number of mutations", main = "Fixed mutations", lwd = 3, ylim=c(0,200), col = "magenta")
lines(x = mut_nonmut_fixed[,1], y = mut_nonmut_fixed[,3], col = "purple", lwd = 3)
legend(x= "topleft", legend=c("Mutator", "Non-mutator"), col= c("magenta","purple"), lty=1, lwd=3)

##
#now to look at the polymorphic mutations by allele

#initiate the matrix
mut_nonmut_polymorphic <- matrix(ncol = 3, nrow = 5)
#add in the days 
mut_nonmut_polymorphic[,1] <- c("0", "17","44","66","90")
#calculate averages for mutator replicates at each point measured - can only do 4 time points since we are combinging plank and biofilm samples
mut_nonmut_polymorphic[1,2] <- mean(x = c(as.numeric(B1_polymorphic[1,]),as.numeric(B2_polymorphic[1,]),as.numeric(B3_polymorphic[1,]),as.numeric(P3_polymorphic[1,])))
mut_nonmut_polymorphic[2,2] <- mean(x = c(as.numeric(B1_polymorphic[2,]),as.numeric(B2_polymorphic[2,]),as.numeric(B3_polymorphic[2,]),as.numeric(P3_polymorphic[2,])))
mut_nonmut_polymorphic[3,2] <- mean(x = c(as.numeric(B1_polymorphic[4,]),as.numeric(B2_polymorphic[4,]),as.numeric(B3_polymorphic[4,]),as.numeric(P3_polymorphic[3,])))
mut_nonmut_polymorphic[4,2] <- mean(x = c(as.numeric(B1_polymorphic[5,]),as.numeric(B2_polymorphic[5,]),as.numeric(B3_polymorphic[5,]),as.numeric(P3_polymorphic[4,])))
mut_nonmut_polymorphic[5,2] <- mean(x = c(as.numeric(B1_polymorphic[7,]),as.numeric(B2_polymorphic[7,]),as.numeric(B3_polymorphic[7,]),as.numeric(P3_polymorphic[5,])))
#calculate averages of planktonic replicates at each day
mut_nonmut_polymorphic[1,3] <- mean(x = c(as.numeric(P1_polymorphic[1,]),as.numeric(P2_polymorphic[1,])))
mut_nonmut_polymorphic[2,3] <- mean(x = c(as.numeric(P1_polymorphic[2,]),as.numeric(P2_polymorphic[2,])))
mut_nonmut_polymorphic[3,3] <- mean(x = c(as.numeric(P1_polymorphic[3,]),as.numeric(P2_polymorphic[3,])))
mut_nonmut_polymorphic[4,3] <- mean(x = c(as.numeric(P1_polymorphic[4,]),as.numeric(P2_polymorphic[4,])))
mut_nonmut_polymorphic[5,3] <- mean(x = c(as.numeric(P1_polymorphic[5,]),as.numeric(P2_polymorphic[5,])))

#plot polymorphic mutations by allele
plot(x = mut_nonmut_polymorphic[,1], y = mut_nonmut_polymorphic[,2], type = "l", xlab = "Day", ylab="Number of mutations", main = "Polymorphic mutations", lwd = 3, ylim=c(0,200), col = "magenta")
lines(x = mut_nonmut_polymorphic[,1], y = mut_nonmut_polymorphic[,3], col = "purple", lwd = 3)
legend(x= "topleft", legend=c("Mutator", "Non-mutator"), col= c("magenta","purple"), lty=1, lwd=3)


##
#now looking at extinct mutations by allele
#initiate the matrix
mut_nonmut_extinct <- matrix(ncol = 3, nrow = 5)
#add in the days 
mut_nonmut_extinct[,1] <- c("0", "17","44","66","90")
#calculate averages for mutator replicates at each point measured - can only do 4 time points since we are combinging plank and biofilm samples
mut_nonmut_extinct[1,2] <- mean(x = c(as.numeric(B1_extinct[2,2]),as.numeric(B2_extinct[2,2]),as.numeric(B3_extinct[2,2]),as.numeric(P3_extinct[2,2])))
mut_nonmut_extinct[2,2] <- mean(x = c(as.numeric(B1_extinct[3,2]),as.numeric(B2_extinct[3,2]),as.numeric(B3_extinct[3,2]),as.numeric(P3_extinct[3,2])))
mut_nonmut_extinct[3,2] <- mean(x = c(as.numeric(B1_extinct[5,2]),as.numeric(B2_extinct[5,2]),as.numeric(B3_extinct[5,2]),as.numeric(P3_extinct[4,2])))
mut_nonmut_extinct[4,2] <- mean(x = c(as.numeric(B1_extinct[6,2]),as.numeric(B2_extinct[6,2]),as.numeric(B3_extinct[6,2]),as.numeric(P3_extinct[5,2])))
mut_nonmut_extinct[5,2] <- mean(x = c(as.numeric(B1_extinct[8,2]),as.numeric(B2_extinct[8,2]),as.numeric(B3_extinct[8,2]),as.numeric(P3_extinct[6,2])))
#calculate averages of planktonic replicates at each day
mut_nonmut_extinct[1,3] <- mean(x = c(as.numeric(P1_extinct[2,2]),as.numeric(P2_extinct[2,2])))
mut_nonmut_extinct[2,3] <- mean(x = c(as.numeric(P1_extinct[3,2]),as.numeric(P2_extinct[3,2])))
mut_nonmut_extinct[3,3] <- mean(x = c(as.numeric(P1_extinct[4,2]),as.numeric(P2_extinct[4,2])))
mut_nonmut_extinct[4,3] <- mean(x = c(as.numeric(P1_extinct[5,2]),as.numeric(P2_extinct[5,2])))
mut_nonmut_extinct[5,3] <- mean(x = c(as.numeric(P1_extinct[6,2]),as.numeric(P2_extinct[6,2])))

#plot total numbers of mutations by environment
plot(x = mut_nonmut_extinct[,1], y = mut_nonmut_extinct[,2], type = "l", xlab = "Day", ylab="Number of mutations", main = "Extinct mutations", lwd = 3, ylim=c(0,200), col="magenta")
lines(x = mut_nonmut_extinct[,1], y = mut_nonmut_extinct[,3], col = "purple", lwd = 3)
legend(x= "topleft", legend=c("Mutator", "Non-mutator"), col= c("magenta","purple"), lty=1, lwd=3)












#now looking at extinct mutations by allele
#initiate the matrix
mut_nonmut_new <- matrix(ncol = 3, nrow = 5)
#add in the days 
mut_nonmut_new[,1] <- c("0", "17","44","66","90")
#calculate averages for mutator replicates at each point measured - can only do 4 time points since we are combinging plank and biofilm samples
mut_nonmut_new[1,2] <- mean(x = c(as.numeric(B1_new[2,2]),as.numeric(B2_new[2,2]),as.numeric(B3_new[2,2]),as.numeric(P3_new[2,2])))
mut_nonmut_new[2,2] <- mean(x = c(as.numeric(B1_new[3,2]),as.numeric(B2_new[3,2]),as.numeric(B3_new[3,2]),as.numeric(P3_new[3,2])))
mut_nonmut_new[3,2] <- mean(x = c(as.numeric(B1_new[5,2]),as.numeric(B2_new[5,2]),as.numeric(B3_new[5,2]),as.numeric(P3_new[4,2])))
mut_nonmut_new[4,2] <- mean(x = c(as.numeric(B1_new[6,2]),as.numeric(B2_new[6,2]),as.numeric(B3_new[6,2]),as.numeric(P3_new[5,2])))
mut_nonmut_new[5,2] <- mean(x = c(as.numeric(B1_new[8,2]),as.numeric(B2_new[8,2]),as.numeric(B3_new[8,2]),as.numeric(P3_new[6,2])))
#calculate averages of planktonic replicates at each day
mut_nonmut_new[1,3] <- mean(x = c(as.numeric(P1_new[2,2]),as.numeric(P2_new[2,2])))
mut_nonmut_new[2,3] <- mean(x = c(as.numeric(P1_new[3,2]),as.numeric(P2_new[3,2])))
mut_nonmut_new[3,3] <- mean(x = c(as.numeric(P1_new[4,2]),as.numeric(P2_new[4,2])))
mut_nonmut_new[4,3] <- mean(x = c(as.numeric(P1_new[5,2]),as.numeric(P2_new[5,2])))
mut_nonmut_new[5,3] <- mean(x = c(as.numeric(P1_new[6,2]),as.numeric(P2_new[6,2])))

#plot total numbers of mutations by environment
plot(x = mut_nonmut_new[,1], y = mut_nonmut_new[,2], type = "l", xlab = "Day", ylab="Number of mutations", main = "New mutations", lwd = 3, ylim=c(0,200), col="magenta")
lines(x = mut_nonmut_new[,1], y = mut_nonmut_new[,3], col = "purple", lwd = 3)
legend(x= "topleft", legend=c("Mutator", "Non-mutator"), col= c("magenta","purple"), lty=1, lwd=3)


######
#plot all of the comparisons by allele
######

par(mfrow=c(2,3))

#plot total numbers of mutations by allele
plot(x = mut_nonmut_total[,1], y = mut_nonmut_total[,2], type = "l", xlab = "Day", ylab="Number of mutations", lwd = 2, ylim=c(0,200), col = "magenta", cex.lab = 1.5)
lines(x = mut_nonmut_total[,1], y = mut_nonmut_total[,3], col = "purple", lwd = 2)
title("A.", adj = 0, cex.main = 2)

#plot detected mutations by allele
plot(x = mut_nonmut_numbers[,1], y = mut_nonmut_numbers[,2], type = "l", xlab = "", ylab="", lwd = 2, ylim=c(0,200), col = "magenta")
lines(x = mut_nonmut_numbers[,1], y = mut_nonmut_numbers[,3], col = "purple", lwd = 2)
title("B.", adj = 0, cex.main = 2)

#plot fixed mutations by allele
plot(x = mut_nonmut_fixed[,1], y = mut_nonmut_fixed[,2], type = "l", xlab = "", ylab="", lwd = 2, col ="magenta", ylim=c(0,200))
lines(x = mut_nonmut_fixed[,1], y = mut_nonmut_fixed[,3], col = "purple", lwd = 2)
legend(x= "topright", legend=c("Mutator", "Non-mutator"), col= c("magenta","purple"), lty=1, lwd=2)
title("C.", adj = 0, cex.main = 2)

#plot polymorphic mutations by allele
plot(x = mut_nonmut_polymorphic[,1], y = mut_nonmut_polymorphic[,2], type = "l", xlab = "", ylab="", lwd = 2, col = "magenta", ylim=c(0,200))
lines(x = mut_nonmut_polymorphic[,1], y = mut_nonmut_polymorphic[,3], col = "purple", lwd = 2)
title("D.", adj = 0, cex.main = 2)

#plot extinct mutations by allele
plot(x = mut_nonmut_extinct[,1], y = mut_nonmut_extinct[,2], type = "l", xlab = "", ylab="", lwd = 2, col = "magenta", ylim=c(0,200))
lines(x = mut_nonmut_extinct[,1], y = mut_nonmut_extinct[,3], col = "purple", lwd = 2)
title("E.", adj = 0, cex.main = 2)

#plot new mutations by allele
plot(x = mut_nonmut_new[,1], y = mut_nonmut_new[,2], type = "l", xlab = "", ylab="", lwd = 2, col = "magenta", ylim=c(0,200))
lines(x = mut_nonmut_new[,1], y = mut_nonmut_new[,3], col = "purple", lwd = 2)
title("F.", adj = 0, cex.main = 2)


######
#print out the tables of values for both environemnt and allele
#######
#printing out tables
colnames(mut_nonmut_new) <- c("Day", "Mutator", "Non Mutator")
colnames(mut_nonmut_extinct) <- c("Day", "Mutator", "Non Mutator")
colnames(mut_nonmut_polymorphic) <- c("Day", "Mutator", "Non Mutator")
colnames(mut_nonmut_fixed) <- c("Day", "Mutator", "Non Mutator")
colnames(mut_nonmut_numbers) <- c("Day", "Mutator", "Non Mutator")
colnames(mut_nonmut_total) <- c("Day", "Mutator", "Non Mutator")

write.csv(mut_nonmut_new, "/Users/katrina/Desktop/muller_v_0.3/Dynamics/mut_nonmut_new.csv")
write.csv(mut_nonmut_extinct, "/Users/katrina/Desktop/muller_v_0.3/Dynamics/mut_nonmut_extinct.csv")
write.csv(mut_nonmut_polymorphic,"/Users/katrina/Desktop/muller_v_0.3/Dynamics/mut_nonmut_polymorphic.csv")
write.csv(mut_nonmut_fixed,"/Users/katrina/Desktop/muller_v_0.3/Dynamics/mut_nonmut_fixed.csv")
write.csv(mut_nonmut_numbers,"/Users/katrina/Desktop/muller_v_0.3/Dynamics/mut_nonmut_numbers.csv")
write.csv(mut_nonmut_total,"/Users/katrina/Desktop/muller_v_0.3/Dynamics/mut_nonmut_total.csv")



colnames(by_environment_new) <- c("Day", "Biofilm", "Planktonic")
colnames(by_environment_extinct) <- c("Day","Biofilm", "Planktonic")
colnames(by_environment_polymorphic) <- c("Day", "Biofilm", "Planktonic")
colnames(by_environment_fixed) <- c("Day", "Biofilm", "Planktonic")
colnames(by_environment_numbers) <- c("Day", "Biofilm", "Planktonic")
colnames(by_environment_total) <- c("Day", "Biofilm", "Planktonic")

write.csv(by_environment_new, "/Users/katrina/Desktop/muller_v_0.3/Dynamics/by_environment_new.csv")
write.csv(by_environment_extinct, "/Users/katrina/Desktop/muller_v_0.3/Dynamics/by_environment_extinct.csv")
write.csv(by_environment_polymorphic,"/Users/katrina/Desktop/muller_v_0.3/Dynamics/by_environment_polymorphic.csv")
write.csv(by_environment_fixed,"/Users/katrina/Desktop/muller_v_0.3/Dynamics/by_environment_fixed.csv")
write.csv(by_environment_numbers,"/Users/katrina/Desktop/muller_v_0.3/Dynamics/by_environment_numbers.csv")
write.csv(by_environment_total,"/Users/katrina/Desktop/muller_v_0.3/Dynamics/by_environment_total.csv")



######
#I want to combine some of these lines in one figure, like is seen in the Lang 2013 paper. 
######
#I want total, polymorphic, extinct, and fixed all on the same graph for the 4 different groups of populaitons 

par(mfrow=c(2,2))

#starting with biofilm
#plot biofilm total, fixed, extinct, and polymorphic
plot(x = by_environment_total[,1], y = by_environment_total[,2], type = "l", xlab = "Day", ylab="Number of mutations", lwd = 2, ylim=c(0,200), cex.lab = 1.5) #total
lines(x = by_environment_fixed[,1], y = by_environment_fixed[,2], col = 2, lwd = 2) #fixed
lines(x = by_environment_extinct[,1], y = by_environment_extinct[,2], col = 3, lwd = 2) #extinct
lines(x = by_environment_polymorphic[,1], y = by_environment_polymorphic[,2], col = 4, lwd = 2) #polymorphic
#legend(x="topleft", legend=c("Total", "Fixed","extinct", "polymorphic"), col =c(1,2,3,4), lwd = 3)
title("A.", adj = 0, cex.main = 2)

#then planktonic
#plot planktonic total, fixed, extinct, and polymorphic
plot(x = by_environment_total[c(1,2,4,5,7),1], y = by_environment_total[c(1,2,4,5,7),3], type = "l", xlab = "", ylab="", lwd = 2, ylim=c(0,200)) #total
lines(x = by_environment_fixed[c(1,2,4,5,7),1], y = by_environment_fixed[c(1,2,4,5,7),3], col = 2, lwd = 2) #fixed
lines(x = by_environment_extinct[c(1,2,4,5,7),1], y = by_environment_extinct[c(1,2,4,5,7),3], col = 3, lwd = 2) #extinct
lines(x = by_environment_polymorphic[c(1,2,4,5,7),1], y = by_environment_polymorphic[c(1,2,4,5,7),3], col = 4, lwd = 2) #polymorphic
#legend(x="topleft", legend=c("Total", "Fixed","extinct", "polymorphic"), col =c(1,2,3,4), lwd = 2)
title("B.", adj = 0, cex.main = 2)

#then mutator
#plot mutator total, fixed, extinct, and polymorphic
plot(x = mut_nonmut_total[,1], y = mut_nonmut_total[,2], type = "l", xlab = "", ylab="", lwd = 2, ylim=c(0,200)) #total
lines(x = mut_nonmut_fixed[,1], y = mut_nonmut_fixed[,2], col = 2, lwd = 2) #fixed
lines(x = mut_nonmut_extinct[,1], y = mut_nonmut_extinct[,2], col = 3, lwd = 2) #extinct
lines(x = mut_nonmut_polymorphic[,1], y = mut_nonmut_polymorphic[,2], col = 4, lwd = 2) #polymorphic
title("C.", adj = 0, cex.main = 2)

#plot NON mutator total, fixed, extinct, and polymorphic
plot(x = mut_nonmut_total[,1], y = mut_nonmut_total[,3], type = "l", xlab = "", ylab="", lwd = 2, ylim=c(0,200)) #total
lines(x = mut_nonmut_fixed[,1], y = mut_nonmut_fixed[,3], col = 2, lwd = 2) #fixed
lines(x = mut_nonmut_extinct[,1], y = mut_nonmut_extinct[,3], col = 3, lwd = 2) #extinct
lines(x = mut_nonmut_polymorphic[,1], y = mut_nonmut_polymorphic[,3], col = 4, lwd = 2) #polymorphic
#legend(x="topleft", legend=c("Total", "Fixed","extinct", "polymorphic"), col =c(1,2,3,4), lwd = 2)
title("D.", adj = 0, cex.main = 2)
