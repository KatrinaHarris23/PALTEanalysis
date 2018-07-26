#PALTE mutational cohort plots

B1_pretty <- read.csv("/Users/katrina/Desktop/working/B1_pretty.csv")
B2_pretty <- read.csv("/Users/katrina/Desktop/working/B2_pretty.csv")
B3_pretty <- read.csv("/Users/katrina/Desktop/working/B3_pretty.csv")
P1_pretty <- read.csv("/Users/katrina/Desktop/working/P1_pretty.csv")
P2_pretty <- read.csv("/Users/katrina/Desktop/working/P2_pretty.csv")
P3_pretty <- read.csv("/Users/katrina/Desktop/working/P3_pretty.csv")




View(B1_pretty)
ncol(B1_pretty) #13

B1_pretty_plot <- B1_pretty[,7:13]
B2_pretty_plot <- B2_pretty[,7:13]
B3_pretty_plot <- B3_pretty[,7:13]


ncol(P1_pretty) #11
P1_pretty_plot <- P1_pretty[,7:11]
P2_pretty_plot <- P2_pretty[,7:11]
P3_pretty_plot <- P3_pretty[,7:11]


layout(matrix(c(1,2,3,4,5,6),2))


plot(NA, xlim=c(0,90), ylim=c(0,100), ylab = "Frequency", xlab = "Time (Days)", main = "B1") #, main = "Mutation frequencies over time") # 
for(i in seq_len(nrow(B1_pretty_plot))){
  lines(c(0,17,25,44,66,75,90),B1_pretty_plot[i,], type="l", col=i)
} #need to get them colored differently, but will work for now. 


plot(NA, xlim=c(0,90), ylim=c(0,100), ylab = "Frequency", xlab = "Time (Days)", main = "P1") #, main = "Mutation frequencies over time") # 
for(i in seq_len(nrow(P1_pretty_plot))){
  lines(c(0,17,44,75,90),P1_pretty_plot[i,], type="l", col=i)
} #need to get them colored differently, but will work for now. 


plot(NA, xlim=c(0,90), ylim=c(0,100), ylab = "Frequency", xlab = "Time (Days)", main = "B2") #, main = "Mutation frequencies over time") # 
for(i in seq_len(nrow(B2_pretty_plot))){
  lines(c(0,17,25,44,66,75,90),B2_pretty_plot[i,], type="l", col=i)
} #need to get them colored differently, but will work for now. 

plot(NA, xlim=c(0,90), ylim=c(0,100), ylab = "Frequency", xlab = "Time (Days)", main = "P2") #, main = "Mutation frequencies over time") # 
for(i in seq_len(nrow(P2_pretty_plot))){
  lines(c(0,17,44,75,90),P2_pretty_plot[i,], type="l", col=i)
} #need to get them colored differently, but will work for now.

plot(NA, xlim=c(0,90), ylim=c(0,100), ylab = "Frequency", xlab = "Time (Days)", main = "B3") #, main = "Mutation frequencies over time") # 
for(i in seq_len(nrow(B3_pretty_plot))){
  lines(c(0,17,25,44,66,75,90),B3_pretty_plot[i,], type="l", col=i)
} #need to get them colored differently, but will work for now. 

plot(NA, xlim=c(0,90), ylim=c(0,100), ylab = "Frequency", xlab = "Time (Days)", main = "P3") #, main = "Mutation frequencies over time") # 
for(i in seq_len(nrow(P3_pretty_plot))){
  lines(c(0,17,44,75,90),P3_pretty_plot[i,], type="l", col=i)
} #need to get them colored differently, but will work for now. 





#now plotting all of the cohort lineages. using deitrick's script because it works better. 
B1_pythoncohorts <- read.table("/Users/katrina/Desktop/working/get_genotypes/B1_Muller.mean.tsv", sep="\t", header=T)
B2_pythoncohorts <- read.table("/Users/katrina/Desktop/working/get_genotypes/B2_Muller.mean.tsv", sep="\t", header=T)
B3_pythoncohorts <- read.table("/Users/katrina/Desktop/working/get_genotypes/B3_Muller.mean.tsv", sep="\t", header=T)
P1_pythoncohorts <- read.table("/Users/katrina/Desktop/working/get_genotypes/P1_Muller.mean.tsv", sep="\t", header=T)
P2_pythoncohorts <- read.table("/Users/katrina/Desktop/working/get_genotypes/P2_Muller.mean.tsv", sep="\t", header=T)
P3_pythoncohorts <- read.table("/Users/katrina/Desktop/working/get_genotypes/P3_Muller.mean.tsv", sep="\t", header=T)
View(B1_pythoncohorts)


B1_cohort_lineage <- B1_pythoncohorts[,-1]
B2_cohort_lineage <- B2_pythoncohorts[-1,-c(1,2)] # I don't know why the B2 cohort output looks different from everything else, it was passed through the exact same pipeline...

B3_cohort_lineage <- B3_pythoncohorts[,-1]
P1_cohort_lineage <- P1_pythoncohorts[,-1]
P2_cohort_lineage <- P2_pythoncohorts[,-1]
P3_cohort_lineage <- P3_pythoncohorts[,-1]


View(P1_cohort_lineage)
View(B2_cohort_lineage)
View(B2_pythoncohorts)


layout(matrix(c(1,2,3,4,5,6),2))

plot(NA, xlim=c(0,90), ylim=c(0,1), ylab = "Frequency", xlab = "Time (Days)", main = "B1") #, main = "Mutation frequencies over time") # 
for(i in seq_len(nrow(B1_cohort_lineage))){
  lines(c(0,17,25,44,66,75,90),B1_cohort_lineage[i,], type="l", col=i)
} #need to get them colored differently, but will work for now. 

plot(NA, xlim=c(0,90), ylim=c(0,1), ylab = "Frequency", xlab = "Time (Days)", main = "P1") #, main = "Mutation frequencies over time") # 
for(i in seq_len(nrow(P1_cohort_lineage))){
  lines(c(0,17,44,75,90),P1_cohort_lineage[i,], type="l", col=i)
} #need to get them colored differently, but will work for now. 

plot(NA, xlim=c(0,90), ylim=c(0,1), ylab = "Frequency", xlab = "Time (Days)", main = "B2") #, main = "Mutation frequencies over time") # 
for(i in seq_len(nrow(B2_cohort_lineage))){
  lines(c(0,17,25,44,66,75,90),B2_cohort_lineage[i,], type="l", col=i)
} #need to get them colored differently, but will work for now. 

plot(NA, xlim=c(0,90), ylim=c(0,1), ylab = "Frequency", xlab = "Time (Days)", main = "P2") #, main = "Mutation frequencies over time") # 
for(i in seq_len(nrow(P2_cohort_lineage))){
  lines(c(0,17,44,75,90),P2_cohort_lineage[i,], type="l", col=i)
} #need to get them colored differently, but will work for now. 

plot(NA, xlim=c(0,90), ylim=c(0,1), ylab = "Frequency", xlab = "Time (Days)", main = "B3") #, main = "Mutation frequencies over time") # 
for(i in seq_len(nrow(B3_cohort_lineage))){
  lines(c(0,17,25,44,66,75,90),B3_cohort_lineage[i,], type="l", col=i)
} #need to get them colored differently, but will work for now. 

plot(NA, xlim=c(0,90), ylim=c(0,1), ylab = "Frequency", xlab = "Time (Days)", main = "P3") #, main = "Mutation frequencies over time") # 
for(i in seq_len(nrow(P3_cohort_lineage))){
  lines(c(0,17,44,75,90),P3_cohort_lineage[i,], type="l", col=i)
} #need to get them colored differently, but will work for now. 