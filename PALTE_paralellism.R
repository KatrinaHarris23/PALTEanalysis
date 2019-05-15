#analysis of the PALTE populations - gene names and functional categories



######
### inport csv files for the different populations and all populations in one sheet, and print out files for each individual population.
######
all_pops <- read.csv("/Users/katrina/Desktop/muller_v_0.3/all_edited.csv", stringsAsFactors = FALSE)


#these files are the same as the files 190424__.csv on the drive in CooperLab > Katrina > PA> PALTE > working files
B1 <- all_pops[(all_pops$Population == "B1"),]
write.csv(B1, "/Users/katrina/Desktop/muller_v_0.3/Population_files_gene_names/B1_names.csv")
B2 <- all_pops[(all_pops$Population == "B2"),]
write.csv(B2, "/Users/katrina/Desktop/muller_v_0.3/Population_files_gene_names/B2_names.csv")
B3 <- all_pops[(all_pops$Population == "B3"),]
write.csv(B3, "/Users/katrina/Desktop/muller_v_0.3/Population_files_gene_names/B3_names.csv")


#planktonic populations have 2 less time points measured, so I have to take out the empty columns in the data before printing or the muller plot scripts will incorporate these. It is easiest to take them out here. They are the 25 and 75 day time points
P1 <- all_pops[(all_pops$Population == "P1"),]
P1 <- P1[,-c(18,21)]
write.csv(P1, "/Users/katrina/Desktop/muller_v_0.3/Population_files_gene_names/P1_names.csv")
P2 <- all_pops[(all_pops$Population == "P2"),]
P2 <- P2[,-c(18,21)]
write.csv(P2, "/Users/katrina/Desktop/muller_v_0.3/Population_files_gene_names/P2_names.csv")
P3 <- all_pops[(all_pops$Population == "P3"),]
P3 <- P3[,-c(18,21)]
write.csv(P3, "/Users/katrina/Desktop/muller_v_0.3/Population_files_gene_names/P3_names.csv")


#########
#comparison of my mutations and the mutations kenny found
######
#read in kenny's data
kenny_data <- read.csv("/Users/katrina/Desktop/kennys_data/kenny_mutation_data.csv", stringsAsFactors = F)

kenny_B1 <- kenny_data[(kenny_data$Population == "B1"),]
kenny_B2 <- kenny_data[(kenny_data$Population == "B2"),]
kenny_B3 <- kenny_data[(kenny_data$Population == "B3"),]


nrow(B1)#211
nrow(kenny_B1)#143

position_same_2 <- kenny_B1[(kenny_B1$Position %in% B1$Position),]
nrow(position_same_2)#72

position_same <- B1[(B1$Position %in% kenny_B1$Position),]
nrow(position_same)#71

View(position_same)
View(position_same_2)
position_only_mine <- B1[!(B1$Position %in% kenny_B1$Position),]
nrow(position_only_mine) #140
position_only_kenny <- kenny_B1[!(kenny_B1$Position %in% B1$Position),]
nrow(position_only_kenny)#71


#######
#read in files from the 6 evolved populations
#######

B1 <- read.csv("/Users/katrina/Desktop/muller_v_0.3/Population_files_gene_names/B1_names.csv", stringsAsFactors = FALSE)
B2 <- read.csv("/Users/katrina/Desktop/muller_v_0.3/Population_files_gene_names/B2_names.csv", stringsAsFactors = FALSE)
B3 <- read.csv("/Users/katrina/Desktop/muller_v_0.3/Population_files_gene_names/B3_names.csv", stringsAsFactors = FALSE)
P1 <- read.csv("/Users/katrina/Desktop/muller_v_0.3/Population_files_gene_names/P1_names.csv", stringsAsFactors = FALSE)
P2 <- read.csv("/Users/katrina/Desktop/muller_v_0.3/Population_files_gene_names/P2_names.csv", stringsAsFactors = FALSE)
P3 <- read.csv("/Users/katrina/Desktop/muller_v_0.3/Population_files_gene_names/P3_names.csv", stringsAsFactors = FALSE)



######
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
write.csv(all, "/Users/katrina/Desktop/muller_v_0.3/all_mutation_numbers.csv")

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
write.csv(all_fixed, "/Users/katrina/Desktop/muller_v_0.3/mutation_breakdown/all_fixed_mutations.csv")

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
write.csv(all_polymorphic, "/Users/katrina/Desktop/muller_v_0.3/mutation_breakdown/all_polymorphic_mutations.csv")


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
write.csv(all_extinct, "/Users/katrina/Desktop/muller_v_0.3/mutation_breakdown/all_extinct_mutations.csv")

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

#########
#determine the "Total" number of mutations seen through time. This is the one plot that should always increase and never decrease.
########

#B1 - should end with 211
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




#######
#plot them all together in one big image so that they can be compared 
#######
par(mfrow=c(2,3))

#total number of mutations
plot(B1_total,xlab="Days", ylab="Number of mutations", main ="Total number of mutations", type="l", ylim = c(0,215))
lines(B2_total, col = 2)
lines(B3_total, col = 3)
lines(P1_total, col = 4, lty = 2)
lines(P2_total, col = 5, lty = 2)
lines(P3_total, col = 6, lty = 2)
#legend(x = "topleft", legend=c("B1", "B2","B3","P1","P2","P3"), lty = c(1,1,1,2,2,2), col=c(1,2,3,4,5,6))

#all mutations above 0 at each time point 
plot(B1_numbers, type = "l", xlab="Days", ylab = "Number of mutations", main = "mutations above 0", ylim = c(0,215))
lines(B2_numbers, col = 2)
lines(B3_numbers, col = 3)
lines(P1_numbers, col = 4, lty = 2)
lines(P2_numbers, col = 5, lty = 2)
lines(P3_numbers, col = 6, lty = 2)
#legend(x = "topleft", legend=c("B1", "B2","B3","P1","P2","P3"), lty = c(1,1,1,2,2,2), col=c(1,2,3,4,5,6))

#fixed mutations
plot(B1_fixed, typ="l", xlim=c(0,90), xlab="Days", ylab = "number of mutations", main = "Fixed mutations", ylim = c(0,215))
lines(B2_fixed, typ="l", col=2)
lines(B3_fixed, typ="l", col=3)
lines(P1_fixed, typ="l", col=4, lty=2)
lines(P2_fixed, typ="l", col=5, lty=2) #is still on there but the P2 and P3 fixation rates are identical.
lines(P3_fixed, typ="l", col=6, lty=2) 
legend(x = "topright", legend=c("B1", "B2","B3","P1","P2","P3"), lty = c(1,1,1,2,2,2), col=c(1,2,3,4,5,6))


#polymorphic mutations 
plot(B1_polymorphic, x=c(0,17,25,44,66,75,90), xlab = "Day", ylab = "Number of mutations", main = "Polymorophic mutations", type = "l", ylim = c(0,215))
lines(B2_polymorphic, x=c(0,17,25,44,66,75,90), col = 2)
lines(B3_polymorphic, x=c(0,17,25,44,66,75,90), col = 3)
lines(P1_polymorphic, x = c(0,17,44,66,90), col = 4, lty=2)
lines(P2_polymorphic, x = c(0,17,44,66,90), col = 5, lty=2)
lines(P3_polymorphic, x = c(0,17,44,66,90), col = 6, lty=2)
#legend(x = "topleft", legend=c("B1", "B2","B3","P1","P2","P3"), lty = c(1,1,1,2,2,2), col=c(1,2,3,4,5,6))

#extinct mutations 
plot(B1_extinct, type="l", main= "Extinct mutations", xlab = "Days", ylab = "Number of mutations", ylim = c(0,215))
lines(B2_extinct, col = 2)
lines(B3_extinct, col = 3)
lines(P1_extinct, col = 4, lty = 2)
lines(P2_extinct, col = 5, lty = 2)
lines(P3_extinct, col = 6, lty = 2)
#legend(x = "topleft", legend=c("B1", "B2","B3","P1","P2","P3"), lty = c(1,1,1,2,2,2), col=c(1,2,3,4,5,6))

#new mutations
plot(B1_new, xlab ="Days", ylab="Number of mutations", main = "New mutations", type ="l", ylim = c(0,215))
lines(B2_new, col = 2)
lines(B3_new, col = 3)
lines(P1_new, col = 4, lty = 2)
lines(P2_new, col = 5, lty = 2)
lines(P3_new, col = 6, lty = 2)
#legend(x = "topright", legend=c("B1", "B2","B3","P1","P2","P3"), lty = c(1,1,1,2,2,2), col=c(1,2,3,4,5,6))

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
plot(x = by_environment_total[,1], y = by_environment_total[,2], type = "l", xlab = "Day", ylab="Number of mutations", main = "Total mutations", lwd = 3, ylim=c(0,200))
lines(x = by_environment_total[c(1,2,4,5,7),1], y = by_environment_total[c(1,2,4,5,7),3], col = 2, lwd = 3)

#plot number of detected mutations by environment
plot(x = by_environment_numbers[,1], y = by_environment_numbers[,2], type = "l", xlab = "Day", ylab="Number of mutations", main = "Detected mutations", lwd = 3, ylim = c(0,200))
lines(x = by_environment_numbers[c(1,2,4,5,7),1], y = by_environment_numbers[c(1,2,4,5,7),3], col = 2, lwd = 3)

#plot fixed mutations by environment
plot(x = by_environment_fixed[,1], y = by_environment_fixed[,2], type = "l", xlab = "Day", ylab="Number of mutations", main = "fixed mutations", lwd = 3, ylim = c(0,200))
lines(x = by_environment_fixed[c(1,2,4,5,7),1], y = by_environment_fixed[c(1,2,4,5,7),3], col = 2, lwd = 3)
legend(x= "topright", legend=c("Biofilm", "Planktonic"), col= c(1,2), lty=1, lwd=3)

#plot polymorphic mutations by environment
plot(x = by_environment_polymorphic[,1], y = by_environment_polymorphic[,2], type = "l", xlab = "Day", ylab="Number of mutations", main = "polymorphic mutations", lwd = 3, ylim = c(0,200))
lines(x = by_environment_polymorphic[c(1,2,4,5,7),1], y = by_environment_polymorphic[c(1,2,4,5,7),3], col = 2, lwd = 3)

#plot extinct mutations by environment
plot(x = by_environment_extinct[,1], y = by_environment_extinct[,2], type = "l", xlab = "Day", ylab="Number of mutations", main = "Extinct mutations", lwd = 3, ylim = c(0,200))
lines(x = by_environment_extinct[c(1,2,4,5,7),1], y = by_environment_extinct[c(1,2,4,5,7),3], col = 2, lwd = 3)

#plot extinct mutations by environment
plot(x = by_environment_new[,1], y = by_environment_new[,2], type = "l", xlab = "Day", ylab="Number of mutations", main = "New mutations", lwd = 3, ylim = c(0,200))
lines(x = by_environment_new[c(1,2,4,5,7),1], y = by_environment_new[c(1,2,4,5,7),3], col = 2, lwd = 3)


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
plot(x = mut_nonmut_total[,1], y = mut_nonmut_total[,2], type = "l", xlab = "Day", ylab="Number of mutations", main = "Total mutations", lwd = 3, ylim=c(0,200))
lines(x = mut_nonmut_total[,1], y = mut_nonmut_total[,3], col = 2, lwd = 3)
legend(x= "topleft", legend=c("Mutator", "Non-Mutator"), col= c(1,2), lty=1, lwd=3)


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
plot(x = mut_nonmut_numbers[,1], y = mut_nonmut_numbers[,2], type = "l", xlab = "Day", ylab="Number of mutations", main = "Detected mutations", lwd = 3, ylim=c(0,200))
lines(x = mut_nonmut_numbers[,1], y = mut_nonmut_numbers[,3], col = 2, lwd = 3)
legend(x= "topleft", legend=c("Mutator", "Non-mutator"), col= c(1,2), lty=1, lwd=3)


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
plot(x = mut_nonmut_fixed[,1], y = mut_nonmut_fixed[,2], type = "l", xlab = "Day", ylab="Number of mutations", main = "Fixed mutations", lwd = 3, ylim=c(0,200))
lines(x = mut_nonmut_fixed[,1], y = mut_nonmut_fixed[,3], col = 2, lwd = 3)
legend(x= "topleft", legend=c("Mutator", "Non-mutator"), col= c(1,2), lty=1, lwd=3)

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
plot(x = mut_nonmut_polymorphic[,1], y = mut_nonmut_polymorphic[,2], type = "l", xlab = "Day", ylab="Number of mutations", main = "Polymorphic mutations", lwd = 3, ylim=c(0,200))
lines(x = mut_nonmut_polymorphic[,1], y = mut_nonmut_polymorphic[,3], col = 2, lwd = 3)
legend(x= "topleft", legend=c("Mutator", "Non-mutator"), col= c(1,2), lty=1, lwd=3)


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
plot(x = mut_nonmut_extinct[,1], y = mut_nonmut_extinct[,2], type = "l", xlab = "Day", ylab="Number of mutations", main = "Extinct mutations", lwd = 3, ylim=c(0,200))
lines(x = mut_nonmut_extinct[,1], y = mut_nonmut_extinct[,3], col = 2, lwd = 3)
legend(x= "topleft", legend=c("Mutator", "Non-mutator"), col= c(1,2), lty=1, lwd=3)
######
#plot all of the comparisons by allele
######

par(mfrow=c(2,3))

#plot total numbers of mutations by allele
plot(x = mut_nonmut_total[,1], y = mut_nonmut_total[,2], type = "l", xlab = "Day", ylab="Number of mutations", main = "Total mutations", lwd = 3, ylim=c(0,200))
lines(x = mut_nonmut_total[,1], y = mut_nonmut_total[,3], col = 2, lwd = 3)

#plot detected mutations by allele
plot(x = mut_nonmut_numbers[,1], y = mut_nonmut_numbers[,2], type = "l", xlab = "Day", ylab="Number of mutations", main = "Detected mutations", lwd = 3, ylim=c(0,200))
lines(x = mut_nonmut_numbers[,1], y = mut_nonmut_numbers[,3], col = 2, lwd = 3)

#plot fixed mutations by allele
plot(x = mut_nonmut_fixed[,1], y = mut_nonmut_fixed[,2], type = "l", xlab = "Day", ylab="Number of mutations", main = "Fixed mutations", lwd = 3, ylim=c(0,200))
lines(x = mut_nonmut_fixed[,1], y = mut_nonmut_fixed[,3], col = 2, lwd = 3)
legend(x= "topright", legend=c("Mutator", "Non-mutator"), col= c(1,2), lty=1, lwd=3)

#plot polymorphic mutations by allele
plot(x = mut_nonmut_polymorphic[,1], y = mut_nonmut_polymorphic[,2], type = "l", xlab = "Day", ylab="Number of mutations", main = "Polymorphic mutations", lwd = 3, ylim=c(0,200))
lines(x = mut_nonmut_polymorphic[,1], y = mut_nonmut_polymorphic[,3], col = 2, lwd = 3)

#plot extinct mutations by allele
plot(x = mut_nonmut_extinct[,1], y = mut_nonmut_extinct[,2], type = "l", xlab = "Day", ylab="Number of mutations", main = "Extinct mutations", lwd = 3, ylim=c(0,200))
lines(x = mut_nonmut_extinct[,1], y = mut_nonmut_extinct[,3], col = 2, lwd = 3)

#plot new mutations by allele
plot(x = mut_nonmut_new[,1], y = mut_nonmut_new[,2], type = "l", xlab = "Day", ylab="Number of mutations", main = "New mutations", lwd = 3, ylim=c(0,200))
lines(x = mut_nonmut_new[,1], y = mut_nonmut_new[,3], col = 2, lwd = 3)



###
#and last thing is to look at the new mutations by allele 
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
plot(x = mut_nonmut_new[,1], y = mut_nonmut_new[,2], type = "l", xlab = "Day", ylab="Number of mutations", main = "New mutations", lwd = 3, ylim=c(0,200))
lines(x = mut_nonmut_new[,1], y = mut_nonmut_new[,3], col = 2, lwd = 3)
legend(x= "topleft", legend=c("Mutator", "Non-mutator"), col= c(1,2), lty=1, lwd=3)


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

write.csv(mut_nonmut_new, "/Users/katrina/Desktop/muller_v_0.3/mutation_breakdown/mut_vs_nonmut/mut_nonmut_new.csv")
write.csv(mut_nonmut_extinct, "/Users/katrina/Desktop/muller_v_0.3/mutation_breakdown/mut_vs_nonmut/mut_nonmut_extinct.csv")
write.csv(mut_nonmut_polymorphic,"/Users/katrina/Desktop/muller_v_0.3/mutation_breakdown/mut_vs_nonmut/mut_nonmut_polymorphic.csv")
write.csv(mut_nonmut_fixed,"/Users/katrina/Desktop/muller_v_0.3/mutation_breakdown/mut_vs_nonmut/mut_nonmut_fixed.csv")
write.csv(mut_nonmut_numbers,"/Users/katrina/Desktop/muller_v_0.3/mutation_breakdown/mut_vs_nonmut/mut_nonmut_numbers.csv")
write.csv(mut_nonmut_total,"/Users/katrina/Desktop/muller_v_0.3/mutation_breakdown/mut_vs_nonmut/mut_nonmut_total.csv")



colnames(by_environment_new) <- c("Day", "Biofilm", "Planktonic")
colnames(by_environment_extinct) <- c("Day","Biofilm", "Planktonic")
colnames(by_environment_polymorphic) <- c("Day", "Biofilm", "Planktonic")
colnames(by_environment_fixed) <- c("Day", "Biofilm", "Planktonic")
colnames(by_environment_numbers) <- c("Day", "Biofilm", "Planktonic")
colnames(by_environment_total) <- c("Day", "Biofilm", "Planktonic")

write.csv(by_environment_new, "/Users/katrina/Desktop/muller_v_0.3/mutation_breakdown/by_environment/by_environment_new.csv")
write.csv(by_environment_extinct, "/Users/katrina/Desktop/muller_v_0.3/mutation_breakdown/by_environment/by_environment_extinct.csv")
write.csv(by_environment_polymorphic,"/Users/katrina/Desktop/muller_v_0.3/mutation_breakdown/by_environment/by_environment_polymorphic.csv")
write.csv(by_environment_fixed,"/Users/katrina/Desktop/muller_v_0.3/mutation_breakdown/by_environment/by_environment_fixed.csv")
write.csv(by_environment_numbers,"/Users/katrina/Desktop/muller_v_0.3/mutation_breakdown/by_environment/by_environment_numbers.csv")
write.csv(by_environment_total,"/Users/katrina/Desktop/muller_v_0.3/mutation_breakdown/by_environment/by_environment_total.csv")

######
#I want to combine some of these lines in one figure, like is seen in the Lang 2013 paper. 
######
#I want total, polymorphic, extinct, and fixed all on the same graph for the 4 different categories

par(mfrow=c(2,2))

#starting with biofilm
#plot biofilm total, fixed, extinct, and polymorphic
plot(x = by_environment_total[,1], y = by_environment_total[,2], type = "l", xlab = "Day", ylab="Number of mutations", main = "Biofilm Mutations", lwd = 3, ylim=c(0,200)) #total
lines(x = by_environment_fixed[,1], y = by_environment_fixed[,2], col = 2, lwd = 3) #fixed
lines(x = by_environment_extinct[,1], y = by_environment_extinct[,2], col = 3, lwd = 3) #extinct
lines(x = by_environment_polymorphic[,1], y = by_environment_polymorphic[,2], col = 4, lwd = 3) #polymorphic
#legend(x="topleft", legend=c("Total", "Fixed","extinct", "polymorphic"), col =c(1,2,3,4), lwd = 3)

#then planktonic
#plot planktonic total, fixed, extinct, and polymorphic
plot(x = by_environment_total[c(1,2,4,5,7),1], y = by_environment_total[c(1,2,4,5,7),3], type = "l", xlab = "Day", ylab="Number of mutations", main = "Planktonic Mutations", lwd = 3, ylim=c(0,200)) #total
lines(x = by_environment_fixed[c(1,2,4,5,7),1], y = by_environment_fixed[c(1,2,4,5,7),3], col = 2, lwd = 3) #fixed
lines(x = by_environment_extinct[c(1,2,4,5,7),1], y = by_environment_extinct[c(1,2,4,5,7),3], col = 3, lwd = 3) #extinct
lines(x = by_environment_polymorphic[c(1,2,4,5,7),1], y = by_environment_polymorphic[c(1,2,4,5,7),3], col = 4, lwd = 3) #polymorphic
#legend(x="topleft", legend=c("Total", "Fixed","extinct", "polymorphic"), col =c(1,2,3,4), lwd = 3)

#then mutator
#plot mutator total, fixed, extinct, and polymorphic
plot(x = mut_nonmut_total[,1], y = mut_nonmut_total[,2], type = "l", xlab = "Day", ylab="Number of mutations", main = "Mutator Mutations", lwd = 3, ylim=c(0,200)) #total
lines(x = mut_nonmut_fixed[,1], y = mut_nonmut_fixed[,2], col = 2, lwd = 3) #fixed
lines(x = mut_nonmut_extinct[,1], y = mut_nonmut_extinct[,2], col = 3, lwd = 3) #extinct
lines(x = mut_nonmut_polymorphic[,1], y = mut_nonmut_polymorphic[,2], col = 4, lwd = 3) #polymorphic
#legend(x="topleft", legend=c("Total", "Fixed","extinct", "polymorphic"), col =c(1,2,3,4), lwd = 3)
#then non mutator


#plot NON mutator total, fixed, extinct, and polymorphic
plot(x = mut_nonmut_total[,1], y = mut_nonmut_total[,3], type = "l", xlab = "Day", ylab="Number of mutations", main = "Non-mutator Mutations", lwd = 3, ylim=c(0,200)) #total
lines(x = mut_nonmut_fixed[,1], y = mut_nonmut_fixed[,3], col = 2, lwd = 3) #fixed
lines(x = mut_nonmut_extinct[,1], y = mut_nonmut_extinct[,3], col = 3, lwd = 3) #extinct
lines(x = mut_nonmut_polymorphic[,1], y = mut_nonmut_polymorphic[,3], col = 4, lwd = 3) #polymorphic
#legend(x="topleft", legend=c("Total", "Fixed","extinct", "polymorphic"), col =c(1,2,3,4), lwd = 3)


#####
#pseudocap function analysis 
########
#now finding the pseudocap functional class numbers
B1_function_occur <- data.frame(table(B1$pseudocap.function))
B2_function_occur <- data.frame(table(B2$pseudocap.function))
B3_function_occur <- data.frame(table(B3$pseudocap.function))
P1_function_occur <- data.frame(table(P1$pseudocap.function))
P2_function_occur <- data.frame(table(P2$pseudocap.function))
P3_function_occur <- data.frame(table(P3$pseudocap.function))

#write the functional analysis out to files
write.csv(B1_function_occur, "/Users/katrina/Desktop/muller_v_0.3/paralellism_within_pop/B1_pseudocap_functions.csv")
write.csv(B2_function_occur, "/Users/katrina/Desktop/muller_v_0.3/paralellism_within_pop/B2_pseudocap_functions.csv")
write.csv(B3_function_occur, "/Users/katrina/Desktop/muller_v_0.3/paralellism_within_pop/B3_pseudocap_functions.csv")
write.csv(P1_function_occur, "/Users/katrina/Desktop/muller_v_0.3/paralellism_within_pop/P1_pseudocap_functions.csv")
write.csv(P2_function_occur, "/Users/katrina/Desktop/muller_v_0.3/paralellism_within_pop/P2_pseudocap_functions.csv")
write.csv(P3_function_occur, "/Users/katrina/Desktop/muller_v_0.3/paralellism_within_pop/P3_pseudocap_functions.csv")


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

head(function_classes)

m_function_classes <- melt(function_classes, id=c("Pop","Var1"), measure.vars="Freq")
function_classes_cast <- t(dcast(m_function_classes, Pop~Var1,mean,Value.var="value",fill=0))
function_classes_final <- as.data.frame(function_classes_cast, header=T)
colnames(function_classes_final) <- as.character(unlist(function_classes_cast[1,]))
function_classes_final <- function_classes_final[-1,]

head(function_classes_final)

write.csv(function_classes_final, "/Users/katrina/Desktop/muller_v_0.3/paralellism_within_pop/pseudocap_function_numbers.csv")


#I added to this file in excel. I added the total number of genes found in the UCBPP PA14 genome in each one of these functional categories. This will allow me to determine if I get certain categories called more often than would be expected purely by chance.

pseudocap_all <- read.csv("/Users/katrina/Desktop/muller_v_0.3/paralellism_within_pop/function/pseudocap_function_numbers_all.csv", stringsAsFactors = F)

pseudocap_all[2,2]

# for anosim fucntion the data needs rows beig the samples and columns being the different variables. so the read in table needs to be transposed
t_pseudocap_all <- read.csv("/Users/katrina/Desktop/muller_v_0.3/paralellism_within_pop/function/t_pseudocap_function_numbers.csv", stringsAsFactors = F)

#need to take out the sample identifier row
t_pseudocap_all <- t_pseudocap_all[,-1]


anosim(t_pseudocap_all,grouping=c("B1","B2","B3","P1","P2","P3","PA14genome"), permutations = 999, distance = "euclidean")









B1_pseudocap <- pseudocap_all[,c(1,3,9)]
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

all_euclidean <- dist(pseudocap_all, method = "euclidean")
#the distance matricies now have to be visualized somehow. Kenny used wards minimum variance method. I am using ward.D2 because it seems like that is the more accepted ward equation. 

B1_clustering <- hclust(B1_euclidean_dist, method="ward.D2", members = NULL)
head(B1_clustering)

plot(B1_clustering)

all_clustering <- hclust(all_euclidean, method = "ward.D2", members = NULL)
plot(all_clustering)


#I don't know what I am seeing here. this is supposed to answer "are any of them seen more often than by chance" and it gives me a tree..... 
t_all_euclidean <- dist(t_pseudocap_all, method = "euclidean")

t_all_clustering <- hclust(t_all_euclidean, method = "ward.D2", members = NULL)
plot(t_all_clustering)

#gene level paralellism and functional classes

#####
#looking at how mutational spectra change over time. For this I am going to use the 

######
#gene level paralellism
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
write.csv(B1_gene_paralellism, "/Users/katrina/Desktop/muller_v_0.3/paralellism_within_pop/B1_gene_paralellism.csv")
write.csv(B2_gene_paralellism, "/Users/katrina/Desktop/muller_v_0.3/paralellism_within_pop/B2_gene_paralellism.csv")
write.csv(B3_gene_paralellism, "/Users/katrina/Desktop/muller_v_0.3/paralellism_within_pop/B3_gene_paralellism.csv")
write.csv(P1_gene_paralellism, "/Users/katrina/Desktop/muller_v_0.3/paralellism_within_pop/P1_gene_paralellism.csv")
write.csv(P2_gene_paralellism, "/Users/katrina/Desktop/muller_v_0.3/paralellism_within_pop/P2_gene_paralellism.csv")
write.csv(P3_gene_paralellism, "/Users/katrina/Desktop/muller_v_0.3/paralellism_within_pop/P3_gene_paralellism.csv")


#now looking at all of the mutations as one batch

all_pops_occur <- data.frame(table(all_pops$Gene))
all_paralellism <- all_pops_occur[all_pops_occur$Freq>1,]

#View(all_paralellism)
write.csv(all_paralellism, "/Users/katrina/Desktop/muller_v_0.3/all_pops_paralellism.csv")

all_pops_function <- data.frame(table(P3$pseudocap.function))
write.csv(all_pops_function, "/Users/katrina/Desktop/muller_v_0.3/all_pops_paralellism.csv")


#Create a table of all paralellism data showing the number in all populations, and the number in the subsequent popuations. 

All_paralell <- read.csv("/Users/katrina/Desktop/muller_v_0.3/all_pops_paralellism.csv", stringsAsFactors = FALSE, header=TRUE)
B1_paralell <- read.csv("/Users/katrina/Desktop/muller_v_0.3/paralellism_within_pop/gene/B1_gene_paralellism.csv", stringsAsFactors = FALSE, header=TRUE)
B2_paralell <- read.csv("/Users/katrina/Desktop/muller_v_0.3/paralellism_within_pop/gene/B2_gene_paralellism.csv", stringsAsFactors = FALSE, header=TRUE)
B3_paralell <- read.csv("/Users/katrina/Desktop/muller_v_0.3/paralellism_within_pop/gene/B3_gene_paralellism.csv", stringsAsFactors = FALSE, header=TRUE)
P1_paralell <- read.csv("/Users/katrina/Desktop/muller_v_0.3/paralellism_within_pop/gene/P1_gene_paralellism.csv", stringsAsFactors = FALSE, header=TRUE)
P2_paralell <- read.csv("/Users/katrina/Desktop/muller_v_0.3/paralellism_within_pop/gene/P2_gene_paralellism.csv", stringsAsFactors = FALSE, header=TRUE)
P3_paralell <- read.csv("/Users/katrina/Desktop/muller_v_0.3/paralellism_within_pop/gene/P3_gene_paralellism.csv", stringsAsFactors = FALSE, header=TRUE)



View(paralellism)
library("vegan")
library("plyr")
library("RColorBrewer")
library("ggplot2")
library("data.table")
library("dplyr")
library("reshape2")
library("scales")


all_paralellism$Pop <- "all"

paralellism <- rbind(all_paralellism, B1_paralell,B2_paralell,B3_paralell,P1_paralell,P2_paralell,P3_paralell)
#View(paralellism)

#melt data frame for casting
m_paralellism <- melt(paralellism,id=c("Pop","Var1"),measure.vars="Freq")

paralellism_cast <- t(dcast(m_paralellism,Pop~Var1,mean,value.var="value",fill=0))
paralellism_final <- as.data.frame(paralellism_cast, header=TRUE)
colnames(paralellism_final) <- as.character(unlist(paralellism_cast[1,]))
paralellism_final <- paralellism_final[-1,]

#nrow(paralellism_final) #166
#print whole data frame: this is within population
write.csv(paralellism_final, "/Users/katrina/Desktop/muller_v_0.3/paralellism_within_pop/gene/within_population_paralellism.csv", col.names = TRUE)



#now to find all of the populations that have the paralell genes
#I can use the tables of gene level paralellism from the previous code here
genes_of_interest <- rownames(paralellism_final)
#View(genes_of_interest)

head(paralellism_final)


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
head(all_paralellism)

paralellism_all <- rbind(all_paralellism, B1_of_interest,B2_of_interest,B3_of_interest,P1_of_interest, P2_of_interest, P3_of_interest)


m_paralellism_all <- melt(paralellism_all,id=c("Pop","Var1"),measure.vars="Freq")
paralellism_all_cast <- t(dcast(m_paralellism_all,Pop~Var1,mean,value.var="value",fill=0))
paralellism_all_final <- as.data.frame(paralellism_all_cast, header=TRUE)
colnames(paralellism_all_final) <- as.character(unlist(paralellism_all_cast[1,]))
paralellism_all_final <- paralellism_all_final[-1,]
#View(paralellism_all_final)  

write.csv(paralellism_all_final, "/Users/katrina/Desktop/muller_v_0.3/paralellism_within_pop/gene/Gene_level_Paralellism_all_populations.csv", col.names = TRUE)


#find gene and function identifiers for the different genes with paralellism
targetgenes <- rownames(paralellism_all_final)
gene_identity <- all_pops[(all_pops$Gene %in% targetgenes),]
gene_identity <- gene_identity[,-c(1,2,3,4,6,8,11:23)]
gene_identity <- gene_identity[,-c(1,2)]
head(gene_identity)
write.csv(gene_identity, "/Users/katrina/Desktop/muller_v_0.3/paralellism_within_pop/gene/gene_level_identity.csv")





#position level paralellism

######

#position level paralellism
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

head(mut_paralell)

m_mut_paralell <- melt(mut_paralell,id=c("Pop","Var1"),measure.vars="Freq" )
mut_paralellism_cast <-t(dcast(m_mut_paralell,Pop~Var1,mean,value.var="value",fill=0))
mut_paralell_final <- as.data.frame(mut_paralellism_cast, header=T)
colnames(mut_paralell_final) <- as.character(unlist(mut_paralellism_cast[1,]))
mut_paralell_final <- mut_paralell_final[-1,]

#View(mut_paralell_final)

mut_paralell_final$Gene <- "" #wish I knew how to assign a gene name for the position, will do manually for the sake of time


write.csv(mut_paralell_final, "/Users/katrina/Desktop/muller_v_0.3/paralellism_within_pop/gene/Mutation_level_Paralellism_all_populations.csv", col.names = TRUE)
View(B3)


#find gene and function identifiers for the different positions with paralellism 
identity <- all_pops[(all_pops$Position %in% mutation_paralellism$Var1),]
identity <- identity[,-c(1,2,3,4,6,8,12:23)]



head(identity)
write.csv(identity, "/Users/katrina/Desktop/muller_v_0.3/paralellism_within_pop/gene/mut_level_identity.csv")





#####

#mutation analysis of all populations
######
#Mutation analysis of all populations 

#set the working directory to where you want your files to print out 
setwd("/Users/katrina/Desktop/muller_v_0.3/mutation_breakdown")

head(B1)

#load the function that I had made previously to identify different functional classes of mutations. This only takes a couple seconds to run.
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

#run the fucntion on all of the originlal populations. Make sure to rename the output file after each run, because it will just overwrite the previous.
Mutations_analysis(P3, "Amino.Acid", "Mutation")

#read in all of the output files

all_mutation_analysis <- read.csv("/Users/katrina/Desktop/muller_v_0.3/mutation_breakdown/all_pops_Mutations_table.csv", stringsAsFactors = F, header = T)
all_mutation_analysis$Pop <- "all"
all_mutation_analysis <- all_mutation_analysis[,-1]


B1_mutation_analysis <- read.csv("/Users/katrina/Desktop/muller_v_0.3/mutation_breakdown/B1_Mutations_table.csv", stringsAsFactors = F, header = T)
B2_mutation_analysis <- read.csv("/Users/katrina/Desktop/muller_v_0.3/mutation_breakdown/B2_Mutations_table.csv", stringsAsFactors = F, header = T)
B3_mutation_analysis <- read.csv("/Users/katrina/Desktop/muller_v_0.3/mutation_breakdown/B3_Mutations_table.csv", stringsAsFactors = F, header = T)
P1_mutation_analysis <- read.csv("/Users/katrina/Desktop/muller_v_0.3/mutation_breakdown/P1_Mutations_table.csv", stringsAsFactors = F, header = T)
P2_mutation_analysis <- read.csv("/Users/katrina/Desktop/muller_v_0.3/mutation_breakdown/P2_Mutations_table.csv", stringsAsFactors = F, header = T)
P3_mutation_analysis <- read.csv("/Users/katrina/Desktop/muller_v_0.3/mutation_breakdown/P3_Mutations_table.csv", stringsAsFactors = F, header = T)

B1_mutation_analysis$Pop <- "B1"
B2_mutation_analysis$Pop <- "B2"
B3_mutation_analysis$Pop <- "B3"
P1_mutation_analysis$Pop <- "P1"
P2_mutation_analysis$Pop <- "P2"
P3_mutation_analysis$Pop <- "P3"

B1_mutation_analysis <- B1_mutation_analysis[,-1]
B2_mutation_analysis <- B2_mutation_analysis[,-1]
B3_mutation_analysis <- B3_mutation_analysis[,-1]
P1_mutation_analysis <- P1_mutation_analysis[,-1]
P2_mutation_analysis <- P2_mutation_analysis[,-1]
P3_mutation_analysis <- P3_mutation_analysis[,-1]

#Combine all tables into one so that I can make a large matrix later. because i just need to merge one column from all of them, I can just use cbind
mutation_analysis_combined <- cbind(all_mutation_analysis[,c(1,2)], B1_mutation_analysis[,2], B2_mutation_analysis[,2], B3_mutation_analysis[,2], P1_mutation_analysis[,2], P2_mutation_analysis[,2], P3_mutation_analysis[,2])

colnames(mutation_analysis_combined) <- c("", "all","B1","B2","B3","P1","P2","P3")


#write out table to file to be able to use in a paper. 
write.csv(mutation_analysis_combined, file = "/Users/katrina/Desktop/muller_v_0.3/mutation_breakdown/all_mutations_analysis.csv")


######
now looking at how the mutational spectra changes over time in all of the populations. I alredy made lists of all mutations detected at the different measured time points (above), so I am going to reuse those. 

#reminder of how they were calculated above: 
#B1_17 <- B1[!(B1$X17 == 0),]

#set the working directory because the below function just prints out a file. Also need to remember to change the name of the file every time I make one.

setwd("/Users/katrina/Desktop/muller_v_0.3/mutation_breakdown")

library("reshape2") #ned to load this package for the function to work

#load the function that I made a while ago:
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

#B1 population

Mutations_analysis(P3_90, "Amino.Acid", "Mutation")


#Now combine them in big tables. 
#read in all files. 
b1_17_breakdown <- read.csv("/Users/katrina/Desktop/muller_v_0.3/mutation_breakdown/B1_17.csv")
b1_25_breakdown <- read.csv("/Users/katrina/Desktop/muller_v_0.3/mutation_breakdown/B1_25.csv")
b1_44_breakdown <- read.csv("/Users/katrina/Desktop/muller_v_0.3/mutation_breakdown/B1_44.csv")
b1_66_breakdown <- read.csv("/Users/katrina/Desktop/muller_v_0.3/mutation_breakdown/B1_66.csv")
b1_75_breakdown <- read.csv("/Users/katrina/Desktop/muller_v_0.3/mutation_breakdown/B1_75.csv")
b1_90_breakdown <- read.csv("/Users/katrina/Desktop/muller_v_0.3/mutation_breakdown/B1_90.csv")

b2_17_breakdown <- read.csv("/Users/katrina/Desktop/muller_v_0.3/mutation_breakdown/B2_17.csv")
b2_25_breakdown <- read.csv("/Users/katrina/Desktop/muller_v_0.3/mutation_breakdown/B2_25.csv")
b2_44_breakdown <- read.csv("/Users/katrina/Desktop/muller_v_0.3/mutation_breakdown/B2_44.csv")
b2_66_breakdown <- read.csv("/Users/katrina/Desktop/muller_v_0.3/mutation_breakdown/B2_66.csv")
b2_75_breakdown <- read.csv("/Users/katrina/Desktop/muller_v_0.3/mutation_breakdown/B2_75.csv")
b2_90_breakdown <- read.csv("/Users/katrina/Desktop/muller_v_0.3/mutation_breakdown/B2_90.csv")

b3_17_breakdown <- read.csv("/Users/katrina/Desktop/muller_v_0.3/mutation_breakdown/B3_17.csv")
b3_25_breakdown <- read.csv("/Users/katrina/Desktop/muller_v_0.3/mutation_breakdown/B3_25.csv")
b3_44_breakdown <- read.csv("/Users/katrina/Desktop/muller_v_0.3/mutation_breakdown/B3_44.csv")
b3_66_breakdown <- read.csv("/Users/katrina/Desktop/muller_v_0.3/mutation_breakdown/B3_66.csv")
b3_75_breakdown <- read.csv("/Users/katrina/Desktop/muller_v_0.3/mutation_breakdown/B3_75.csv")
b3_90_breakdown <- read.csv("/Users/katrina/Desktop/muller_v_0.3/mutation_breakdown/B3_90.csv")

p1_17_breakdown <- read.csv("/Users/katrina/Desktop/muller_v_0.3/mutation_breakdown/p1_17.csv")
p1_44_breakdown <- read.csv("/Users/katrina/Desktop/muller_v_0.3/mutation_breakdown/p1_44.csv")
p1_66_breakdown <- read.csv("/Users/katrina/Desktop/muller_v_0.3/mutation_breakdown/p1_66.csv")
p1_90_breakdown <- read.csv("/Users/katrina/Desktop/muller_v_0.3/mutation_breakdown/p1_90.csv")

p2_17_breakdown <- read.csv("/Users/katrina/Desktop/muller_v_0.3/mutation_breakdown/p2_17.csv")
p2_44_breakdown <- read.csv("/Users/katrina/Desktop/muller_v_0.3/mutation_breakdown/p2_44.csv")
p2_66_breakdown <- read.csv("/Users/katrina/Desktop/muller_v_0.3/mutation_breakdown/p2_66.csv")
p2_90_breakdown <- read.csv("/Users/katrina/Desktop/muller_v_0.3/mutation_breakdown/p2_90.csv")

p3_17_breakdown <- read.csv("/Users/katrina/Desktop/muller_v_0.3/mutation_breakdown/p3_17.csv")
p3_44_breakdown <- read.csv("/Users/katrina/Desktop/muller_v_0.3/mutation_breakdown/p3_44.csv")
p3_66_breakdown <- read.csv("/Users/katrina/Desktop/muller_v_0.3/mutation_breakdown/p3_66.csv")
p3_90_breakdown <- read.csv("/Users/katrina/Desktop/muller_v_0.3/mutation_breakdown/p3_90.csv")

#and combine columns to make a big table for each population

b1_breakdown <- as.data.frame(cbind(b1_17_breakdown[,2:3],b1_25_breakdown[,3],b1_44_breakdown[,3],b1_66_breakdown[,3],b1_75_breakdown[,3],b1_90_breakdown[,3]))
b2_breakdown <- as.data.frame(cbind(b2_17_breakdown[,2:3],b2_25_breakdown[,3],b2_44_breakdown[,3],b2_66_breakdown[,3],b2_75_breakdown[,3],b2_90_breakdown[,3]))
b3_breakdown <- as.data.frame(cbind(b3_17_breakdown[,2:3],b3_25_breakdown[,3],b3_44_breakdown[,3],b3_66_breakdown[,3],b3_75_breakdown[,3],b3_90_breakdown[,3]))

p1_breakdown <- as.data.frame(cbind(p1_17_breakdown[,2:3],p1_44_breakdown[,3],p1_66_breakdown[,3],p1_90_breakdown[,3]))
p2_breakdown <- as.data.frame(cbind(p2_17_breakdown[,2:3],p2_44_breakdown[,3],p2_66_breakdown[,3],p2_90_breakdown[,3]))
p3_breakdown <- as.data.frame(cbind(p3_17_breakdown[,2:3],p3_44_breakdown[,3],p3_66_breakdown[,3],p3_90_breakdown[,3]))

#change the column names to the dates.
colnames(b1_breakdown) <- c("Category", "Day 17", "Day 25", "Day 44", "Day 66", "Day 75", "Day 90")
colnames(b2_breakdown) <- c("Category", "Day 17", "Day 25", "Day 44", "Day 66", "Day 75", "Day 90")
colnames(b3_breakdown) <- c("Category", "Day 17", "Day 25", "Day 44", "Day 66", "Day 75", "Day 90")

colnames(p1_breakdown) <- c("Category", "Day 17", "Day 44", "Day 66", "Day 90")
colnames(p2_breakdown) <- c("Category", "Day 17", "Day 44", "Day 66", "Day 90")
colnames(p3_breakdown) <- c("Category", "Day 17", "Day 44", "Day 66", "Day 90")

#print out the final tables
write.csv(b1_breakdown, "/Users/katrina/Desktop/muller_v_0.3/mutation_breakdown/B1_mutations_breakdown.csv")
write.csv(b2_breakdown, "/Users/katrina/Desktop/muller_v_0.3/mutation_breakdown/B2_mutations_breakdown.csv")
write.csv(b3_breakdown, "/Users/katrina/Desktop/muller_v_0.3/mutation_breakdown/B3_mutations_breakdown.csv")
write.csv(p1_breakdown, "/Users/katrina/Desktop/muller_v_0.3/mutation_breakdown/P1_mutations_breakdown.csv")
write.csv(p2_breakdown, "/Users/katrina/Desktop/muller_v_0.3/mutation_breakdown/P2_mutations_breakdown.csv")
write.csv(p3_breakdown, "/Users/katrina/Desktop/muller_v_0.3/mutation_breakdown/P3_mutations_breakdown.csv")
######
#Now trying to get an accurate count of mutations before ANY filtering but taking out the ancestral mutations. 
######

#I am now going to look at all of the populations and measure the amount of diversity 
######
#I am going to be measuring diversity in all populations for each time point that they were measured in. 

#B1 pop first
limited_B1 <- B1[,-c(1:16)] #only selecting the frequency rows
t_limited_B1 <- t(limited_B1) #transform the matrix so that you are measuring the diversity at the various time points.
class(t_limited_B1) <- "numeric" #change all numbers to type numeric so that the diversity function will work.
B1_simpson <- diversity(t_limited_B1, index = "simpson")
B1_inv_simpson <- diversity(t_limited_B1, index = "invsimpson")
B1_shannon <- diversity(t_limited_B1, index = "shannon")


#B2 pop
limited_B2 <- B2[,-c(1:16)]
t_limited_B2 <- t(limited_B2)
class(t_limited_B2) <- "numeric"
B2_simpson <- diversity(t_limited_B2, index = "simpson")
B2_inv_simpson <- diversity(t_limited_B2, index = "invsimpson")
B2_shannon <- diversity(t_limited_B2, index = "shannon")


#B3 pop
limited_B3 <- B3[,-c(1:16)]
t_limited_B3 <- t(limited_B3)
class(t_limited_B3) <- "numeric"
B3_simpson <- diversity(t_limited_B3, index = "simpson")
B3_inv_simpson <- diversity(t_limited_B3, index = "invsimpson")
B3_shannon <- diversity(t_limited_B3, index = "shannon")


#P1 pop
limited_P1 <- P1[,-c(1:16,19,22)]
t_limited_P1 <- t(limited_P1)
class(t_limited_P1) <- "numeric"
P1_simpson <- diversity(t_limited_P1, index = "simpson")
P1_inv_simpson <- diversity(t_limited_P1, index = "invsimpson")
P1_shannon <- diversity(t_limited_P1, index = "shannon")

#P2 pop
limited_P2 <- P2[,-c(1:16,19,22)]
t_limited_P2 <- t(limited_P2)
class(t_limited_P2) <- "numeric"
P2_simpson <- diversity(t_limited_P2, index = "simpson")
P2_inv_simpson <- diversity(t_limited_P2, index = "invsimpson")
P2_shannon <- diversity(t_limited_P2, index = "shannon")


#P3 pop
limited_P3 <- P3[,-c(1:16,19,22)]
t_limited_P3 <- t(limited_P3)
class(t_limited_P3) <- "numeric"
P3_simpson <- diversity(t_limited_P3, index = "simpson")
P3_inv_simpson <- diversity(t_limited_P3, index = "invsimpson")
P3_shannon <- diversity(t_limited_P3, index = "shannon")

#now to combine into a large matrix 
biof_diversity <- cbind(B1_shannon, B2_shannon, B3_shannon, B1_simpson, B2_simpson, B3_simpson, B1_inv_simpson, B2_inv_simpson, B3_inv_simpson)

plank_diversity <- cbind(P1_shannon, P2_shannon, P3_shannon, P1_simpson, P2_simpson, P3_simpson, P1_inv_simpson, P2_inv_simpson, P3_inv_simpson)

#write both out to file so that I can do statistics in prism
write.csv (biof_diversity, file = "/Users/katrina/Desktop/muller_v_0.3/diversity/biofilm_pops_diversity.csv")
write.csv (plank_diversity, file = "/Users/katrina/Desktop/muller_v_0.3/diversity/plank_pops_diversity.csv")

#I don't know how to do ttests in R so I do stats in prism. 
#Stats to do: 
# pairwise ttests for the three biofilm populations vs the three planktonic populations. 
#pairwise ttests for the 4 mutator populations compared to the 2 non mutator populations.






