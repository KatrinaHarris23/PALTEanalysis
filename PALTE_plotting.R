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


#change the column names so that we are working in generations
colnames(B1_pretty_plot) <- c("0","113","167","293","440","500","600")
colnames(B2_pretty_plot) <- c("0","113","167","293","440","500","600")
colnames(B3_pretty_plot) <- c("0","113","167","293","440","500","600")

colnames(P1_pretty_plot) <- c("0","113","293","500","600")
colnames(P2_pretty_plot) <- c("0","113","293","500","600")
colnames(P3_pretty_plot) <- c("0","113","293","500","600")



b_generations <- c(0,113,167,293,440,500,600)
p_generations <- c(0,113,293,500,600)


layout(matrix(c(1,2,3,4,5,6),2))


plot(NA, xlim=c(0,600), ylim=c(0,100), ylab = "Frequency", xlab = "Generations", main = "B1") #, main = "Mutation frequencies over time") # 
for(i in seq_len(nrow(B1_pretty_plot))){
  lines(b_generations,B1_pretty_plot[i,], type="l", col=i)
} #need to get them colored differently, but will work for now. 


plot(NA, xlim=c(0,600), ylim=c(0,100), ylab = "", xlab = "", main = "P1") #, main = "Mutation frequencies over time") # 
for(i in seq_len(nrow(P1_pretty_plot))){
  lines(p_generations,P1_pretty_plot[i,], type="l", col=i)
} #need to get them colored differently, but will work for now. 


plot(NA, xlim=c(0,600), ylim=c(0,100), ylab = "", xlab = "", main = "B2") #, main = "Mutation frequencies over time") # 
for(i in seq_len(nrow(B2_pretty_plot))){
  lines(b_generations,B2_pretty_plot[i,], type="l", col=i)
} #need to get them colored differently, but will work for now. 

plot(NA, xlim=c(0,600), ylim=c(0,100), ylab = "", xlab = "", main = "P2") #, main = "Mutation frequencies over time") # 
for(i in seq_len(nrow(P2_pretty_plot))){
  lines(p_generations,P2_pretty_plot[i,], type="l", col=i)
} #need to get them colored differently, but will work for now.

plot(NA, xlim=c(0,600), ylim=c(0,100), ylab = "", xlab = "", main = "B3") #, main = "Mutation frequencies over time") # 
for(i in seq_len(nrow(B3_pretty_plot))){
  lines(b_generations,B3_pretty_plot[i,], type="l", col=i)
} #need to get them colored differently, but will work for now. 

plot(NA, xlim=c(0,600), ylim=c(0,100), ylab = "", xlab = "", main = "P3") #, main = "Mutation frequencies over time") # 
for(i in seq_len(nrow(P3_pretty_plot))){
  lines(p_generations,P3_pretty_plot[i,], type="l", col=i)
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


b1_cohort_lineage_2 <- B1_cohort_lineage *100
b2_cohort_lineage_2 <- B2_cohort_lineage *100
b3_cohort_lineage_2 <- B3_cohort_lineage *100
p1_cohort_lineage_2 <- P1_cohort_lineage *100
p2_cohort_lineage_2 <- P2_cohort_lineage *100
p3_cohort_lineage_2 <- P3_cohort_lineage *100


View(P1_cohort_lineage)
View(B2_cohort_lineage)
View(B2_pythoncohorts)


layout(matrix(c(1,2,3,4,5,6),2))

plot(NA, xlim=c(0,600), ylim=c(0,100), ylab = "Frequency", xlab = "Time (Days)", main = "B1") #, main = "Mutation frequencies over time") # 
for(i in seq_len(nrow(b1_cohort_lineage_2))){
  lines(b_generations,b1_cohort_lineage_2[i,], type="l", col=i)
} #need to get them colored differently, but will work for now. 

plot(NA, xlim=c(0,600), ylim=c(0,100), ylab = "", xlab = "", main = "P1") #, main = "Mutation frequencies over time") # 
for(i in seq_len(nrow(p1_cohort_lineage_2))){
  lines(p_generations,p1_cohort_lineage_2[i,], type="l", col=i)
} #need to get them colored differently, but will work for now. 

plot(NA, xlim=c(0,600), ylim=c(0,100), ylab = "", xlab = "", main = "B2") #, main = "Mutation frequencies over time") # 
for(i in seq_len(nrow(b2_cohort_lineage_2))){
  lines(b_generations,b2_cohort_lineage_2[i,], type="l", col=i)
} #need to get them colored differently, but will work for now. 

plot(NA, xlim=c(0,600), ylim=c(0,100), ylab = "", xlab = "", main = "P2") #, main = "Mutation frequencies over time") # 
for(i in seq_len(nrow(p2_cohort_lineage_2))){
  lines(p_generations,p2_cohort_lineage_2[i,], type="l", col=i)
} #need to get them colored differently, but will work for now. 

plot(NA, xlim=c(0,600), ylim=c(0,100), ylab = "", xlab = "", main = "B3") #, main = "Mutation frequencies over time") # 
for(i in seq_len(nrow(b3_cohort_lineage_2))){
  lines(b_generations,b3_cohort_lineage_2[i,], type="l", col=i)
} #need to get them colored differently, but will work for now. 

plot(NA, xlim=c(0,600), ylim=c(0,100), ylab = "", xlab = "", main = "P3") #, main = "Mutation frequencies over time") # 
for(i in seq_len(nrow(p3_cohort_lineage_2))){
  lines(p_generations,p3_cohort_lineage_2[i,], type="l", col=i)
} #need to get them colored differently, but will work for now. 






###now to see if I can extract the cohort members from the rows of frequency data 