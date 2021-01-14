##Have done a second experiment where I split the same population into planktonic and biofilm regimes that are propagated for a couple days. This is the new analysis. 

library(data.table)
#coverage for all the samples
#B_B1 - 328
#B_B2 - 391
#B_P1 - 375
#B_P2 - 331
#P_B1 - 410
#P_B2 - 376
#P_P1 - 353
#P_P2 - 373

#Import in the data and the ancestor
#data with all evolved population samples
populations <- read.csv("/Users/katrina/Desktop/201016_ecological_interactions/ecological_comparison_Breseq_Output.csv", stringsAsFactors = F)
ancestor <- read.csv("/Users/katrina/Desktop/PALTE_final/working/KBH5_WT_Breseq_Output.csv",stringsAsFactors = F)

#filter out mutations that were found in the ancestor. 
populations_noancestor <- populations[!(populations$position %in% ancestor$SeqID),]

#determine the number of mutations found in each sample
number_of_mutations <- table(populations_noancestor$SampleName)

#melt the data so that I can reshape it
m_pops <- melt(populations_noancestor, id=c("SampleName","evidence","position","mutation","annotation","gene","description"), measure.vars = "freq")

#reshape the data so that rows are freqencies of positions and columns are each sample

c_pops <- as.data.frame(t(dcast(m_pops, SampleName~position, mean, value.var ="value", fill = 0)))
colnames(c_pops) <- as.character(unlist(c_pops[1,]))
c_pops <- c_pops[-1,]


#now to look at the B1 population
B1 <- c_pops[,c(1,5)]
write.csv(B1, "/Users/katrina/Desktop/201016_ecological_interactions/B1.csv")
B1 <- read.csv("/Users/katrina/Desktop/201016_ecological_interactions/B1.csv", stringsAsFactors = F)

B1 <- B1[!(B1$B_B1 == 0),] #remove all cases that are not detected in the biofilm regime
B1 <- B1[!(B1$P_B1 == 0),] #remove all cases that are not detected in the planktonic regime.

#B2 population
B2 <- c_pops[,c(2,6)]
write.csv(B2, "/Users/katrina/Desktop/201016_ecological_interactions/B2.csv")
B2 <- read.csv("/Users/katrina/Desktop/201016_ecological_interactions/B2.csv", stringsAsFactors = F)

B2 <- B2[!(B2$B_B2 == 0),] #remove all cases that are not detected in the biofilm regime
B2 <- B2[!(B2$P_B2 == 0),] #remove all cases that are not detected in the planktonic regime.


#P1
P1 <- c_pops[,c(3,7)]
write.csv(P1, "/Users/katrina/Desktop/201016_ecological_interactions/P1.csv")
P1 <- read.csv("/Users/katrina/Desktop/201016_ecological_interactions/P1.csv", stringsAsFactors = F)

P1 <- P1[!(P1$B_P1 == 0),] #remove all cases that are not detected in the biofilm regime
P1 <- P1[!(P1$P_P1 == 0),] #remove all cases that are not detected in the planktonic regime.

#P2 population
P2 <- c_pops[,c(4,8)]
write.csv(P2, "/Users/katrina/Desktop/201016_ecological_interactions/P2.csv")
P2 <- read.csv("/Users/katrina/Desktop/201016_ecological_interactions/P2.csv", stringsAsFactors = F)

P2 <- P2[!(P2$B_P2 == 0),] #remove all cases that are not detected in the biofilm regime
P2 <- P2[!(P2$P_P2 == 0),] #remove all cases that are not detected in the planktonic regime.



#plot the 4 populations with biofilm selection on the x axis and planktonic selection on the y axis. Just looking at all mutations in the data set. 

#import all of the data sets: 
B1 <- read.csv("/Users/katrina/Desktop/201016_ecological_interactions/B1.csv", stringsAsFactors = F)
B1 <- B1[!(B1$B_B1 == 0),] #remove all cases that are not detected in the biofilm regime
B1 <- B1[!(B1$P_B1 == 0),] #remove all cases that are not detected in the planktonic regime.

B2 <- read.csv("/Users/katrina/Desktop/201016_ecological_interactions/B2.csv", stringsAsFactors = F)
B2 <- B2[!(B2$B_B2 == 0),] #remove all cases that are not detected in the biofilm regime
B2 <- B2[!(B2$P_B2 == 0),] #remove all cases that are not detected in the planktonic regime.


P1 <- read.csv("/Users/katrina/Desktop/201016_ecological_interactions/P1.csv", stringsAsFactors = F)
P1 <- P1[!(P1$B_P1 == 0),] #remove all cases that are not detected in the biofilm regime
P1 <- P1[!(P1$P_P1 == 0),] #remove all cases that are not detected in the planktonic regime.


P2 <- read.csv("/Users/katrina/Desktop/201016_ecological_interactions/P2.csv", stringsAsFactors = F)
P2 <- P2[!(P2$B_P2 == 0),] #remove all cases that are not detected in the biofilm regime
P2 <- P2[!(P2$P_P2 == 0),] #remove all cases that are not detected in the planktonic regime.



#trying to detect outliars using a multivariate approach
#starting with B1 population
mod <- lm(B1[,2] ~ B1[,3]) #calculating a linear regression for the data
cooksd <- cooks.distance(mod) #calculating cooks distance - generally those data points that have a cooks distance 4 or greater, are influential, or outliars. 

outliar_criteria <- 4*mean(cooksd, na.rm=T)
B1$cooksd <- cooksd
B1$significant <- cooksd>=outliar_criteria
write.csv(B1, "/Users/katrina/Desktop/201016_ecological_interactions/cooksd/B1_cooksd.csv")

plot(cooksd, pch="*", cex=2, main="Influential Obs by Cooks distance: B1")  # plot cook's distance
abline(h = 4*mean(cooksd, na.rm=T), col="red")  # add cutoff line
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4*mean(cooksd, na.rm=T),names(cooksd),""), col="red")  # add labels



#starting with B2 population
mod <- lm(B2[,2] ~ B2[,3]) #calculating a linear regression for the data
cooksd <- cooks.distance(mod) #calculating cooks distance - generally those data points that have a cooks distance 4 or greater, are influential, or outliars. 

outliar_criteria <- 4*mean(cooksd, na.rm=T)
B2$cooksd <- cooksd
B2$significant <- cooksd>=outliar_criteria
write.csv(B2, "/Users/katrina/Desktop/201016_ecological_interactions/cooksd/B2_cooksd.csv")

plot(cooksd, pch="*", cex=2, main="Influential Obs by Cooks distance: B2")  # plot cook's distance
abline(h = 4*mean(cooksd, na.rm=T), col="red")  # add cutoff line
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4*mean(cooksd, na.rm=T),names(cooksd),""), col="red")  # add labels



#starting with P1 population
mod <- lm(P1[,2] ~ P1[,3]) #calculating a linear regression for the data
cooksd <- cooks.distance(mod) #calculating cooks distance - generally those data points that have a cooks distance 4 or greater, are influential, or outliars. 

outliar_criteria <- 4*mean(cooksd, na.rm=T)
P1$cooksd <- cooksd
P1$significant <- cooksd>=outliar_criteria
write.csv(P1, "/Users/katrina/Desktop/201016_ecological_interactions/cooksd/P1_cooksd.csv")

plot(cooksd, pch="*", cex=2, main="Influential Obs by Cooks distance: P1")  # plot cook's distance
abline(h = 4*mean(cooksd, na.rm=T), col="red")  # add cutoff line
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4*mean(cooksd, na.rm=T),names(cooksd),""), col="red")  # add labels


#starting with P2 population
mod <- lm(P2[,2] ~ P2[,3]) #calculating a linear regression for the data
cooksd <- cooks.distance(mod) #calculating cooks distance - generally those data points that have a cooks distance 4 or greater, are influential, or outliars. 

outliar_criteria <- 4*mean(cooksd, na.rm=T)
P2$cooksd <- cooksd
P2$significant <- cooksd>=outliar_criteria
write.csv(P2, "/Users/katrina/Desktop/201016_ecological_interactions/cooksd/P2_cooksd.csv")

plot(cooksd, pch="*", cex=2, main="Influential Obs by Cooks distance: P2")  # plot cook's distance
abline(h = 4*mean(cooksd, na.rm=T), col="red")  # add cutoff line
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4*mean(cooksd, na.rm=T),names(cooksd),""), col="red")  # add labels




library(ggplot2)
theme_set(theme_bw())

ggplot(B1[,2:3], aes(x =B_B1, y=P_B1, size = 3)) + 
  geom_point() + 
  xlab("Biofilm selection") +
  ylab("Planktonic selection") +
  ggtitle("A. B1") +
  theme(aspect.ratio=1)+
  geom_abline(intercept=0)+
  ylim(0,100)+
  xlim(0,100)+
  theme(legend.position = "none")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  theme(plot.title = element_text(size = 40))

ggplot(B2[,2:3], aes(x =B_B2, y=P_B2, size = 3)) + 
  geom_point() + 
  xlab("Biofilm selection") +
  ylab("Planktonic selection") +
  ggtitle("B. B2") +
  theme(aspect.ratio=1)+
  geom_abline(intercept=0)+
  ylim(0,100)+
  xlim(0,100)+
  theme(legend.position = "none")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  theme(plot.title = element_text(size = 40))

ggplot(P1[,2:3], aes(x =B_P1, y=P_P1, size = 3)) + 
  geom_point() + 
  xlab("Biofilm selection") +
  ylab("Planktonic selection") +
  ggtitle("C. P1") +
  theme(aspect.ratio=1)+
  geom_abline(intercept=0)+
  ylim(0,100)+
  xlim(0,100)+
  theme(legend.position = "none")+
  theme(axis.text = element_text(size = 28), axis.title=element_text(size=30, face="bold"))+
  theme(plot.title = element_text(size = 40))

ggplot(P2[,2:3], aes(x =B_P2, y=P_P2, size = 3)) + 
  geom_point() + 
  xlab("Biofilm selection") +
  ylab("Planktonic selection") +
  ggtitle("D. P2") +
  theme(aspect.ratio=1)+
  geom_abline(intercept=0)+
  ylim(0,100)+
  xlim(0,100)+
  theme(legend.position = "none")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  theme(plot.title = element_text(size = 40))
  


#Now going to color each population based on if it is higher frequency in biofilm or planktonic conditions. 

B1$color <- 1
for (i in 1:nrow(B1)) {
  if (B1$significant = TRUE)
    B1[i,6] <- "Planktonic enriched"
  else if (B1[i,2] > B1[i,3])
    B1[i,4] <- "Biofilm enriched"
  else if (B1[i,2] == B1[i,3])
    B1[i,4] <- "equal"
  else
    B1[i,4] <- "other"
}
B1$color <- as.factor(B1$color)


B2$color <- 1
for (i in 1:nrow(B2)) {
  if (B2[i,2] < B2[i,3])
    B2[i,4] <- "Planktonic enriched"
  else if (B2[i,2] > B2[i,3])
    B2[i,4] <- "Biofilm enriched"
  else if (B2[i,2] == B2[i,3])
    B2[i,4] <- "equal"
  else
    B2[i,4] <- "other"
}
B2$color <- as.factor(B2$color)

P1$color <- 1
for (i in 1:nrow(P1)) {
  if (P1[i,2] < P1[i,3])
    P1[i,4] <- "Planktonic enriched"
  else if (P1[i,2] > P1[i,3])
    P1[i,4] <- "Biofilm enriched"
  else if (P1[i,2] == P1[i,3])
    P1[i,4] <- "equal"
  else
    P1[i,4] <- "other"
}
P1$color <- as.factor(P1$color)


P2$color <- 1
for (i in 1:nrow(P2)) {
  if (P2[i,2] < P2[i,3])
    P2[i,4] <- "Planktonic enriched"
  else if (P2[i,2] > P2[i,3])
    P2[i,4] <- "Biofilm enriched"
  else if (P2[i,2] == P2[i,3])
    P2[i,4] <- "equal"
  else
    P2[i,4] <- "other"
}
P2$color <- as.factor(P2$color)



#plot according to significance

ggplot(B1[,2:3], aes(x =B_B1, y=P_B1, size = 3, color =B1$significant)) + 
  geom_point() + 
  xlab("Biofilm selection") +
  ylab("Planktonic selection") +
  ggtitle("A. B1") +
  theme(aspect.ratio=1)+
  geom_abline(intercept=0)+
  ylim(0,100)+
  xlim(0,100)+
  theme(legend.position = "none")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  theme(plot.title = element_text(size = 40))+ 
  scale_color_manual(values=c("Black", "Red"))


ggplot(B2[,2:3], aes(x =B_B2, y=P_B2, size = 3, color =B2$significant)) + 
  geom_point() + 
  xlab("Biofilm selection") +
  ylab("Planktonic selection") +
  ggtitle("B. B2") +
  theme(aspect.ratio=1)+
  geom_abline(intercept=0)+
  ylim(0,100)+
  xlim(0,100)+
  theme(legend.position = "right")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  theme(plot.title = element_text(size = 40))+ 
  scale_color_manual(values=c("Black", "Red"))

ggplot(P1[,2:3], aes(x =B_P1, y=P_P1, size = 3, color =P1$significant)) + 
  geom_point() + 
  xlab("Biofilm selection") +
  ylab("Planktonic selection") +
  ggtitle("C. P1") +
  theme(aspect.ratio=1)+
  geom_abline(intercept=0)+
  ylim(0,100)+
  xlim(0,100)+
  theme(legend.position = "none")+
  theme(axis.text = element_text(size = 28), axis.title=element_text(size=30, face="bold"))+
  theme(plot.title = element_text(size = 40))+ 
  scale_color_manual(values=c("Black", "Red"))

ggplot(P2[,2:3], aes(x =B_P2, y=P_P2, size = 3, color =P2$significant)) + 
  geom_point() + 
  xlab("Biofilm selection") +
  ylab("Planktonic selection") +
  ggtitle("D. P2") +
  theme(aspect.ratio=1)+
  geom_abline(intercept=0)+
  ylim(0,100)+
  xlim(0,100)+
  theme(legend.position = "none")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  theme(plot.title = element_text(size = 40))+ 
  scale_color_manual(values=c("Black", "Red"))




write.csv(B2, "/Users/katrina/Desktop/201016_ecological_interactions/B2_colors.csv")




#####
#statistics

#B1
B1_spearman <- B1[,2:3]
B1_spearman$total <- B1_spearman[,1] + B1_spearman[,2]
B1_spearman <- B1_spearman[(B1_spearman$total != 0.0),]

cor.test(~ P_B1 + B_B1,
         data = B1_spearman, 
         method = "pearson", 
         continuity = F, 
         conf.level = 0.95)

#B2
B2_spearman <- B2[,2:3]
B2_spearman$total <- B2_spearman[,1] + B2_spearman[,2]
B2_spearman <- B2_spearman[(B2_spearman$total != 0.0),]

cor.test(~ P_B2 + B_B2,
         data = B2_spearman, 
         method = "pearson", 
         continuity = F, 
         conf.level = 0.95)

#P1
P1_spearman <- P1[,2:3]
P1_spearman$total <- P1_spearman[,1] + P1_spearman[,2]
P1_spearman <- P1_spearman[(P1_spearman$total != 0.0),]

cor.test(~ P_P1 + B_P1,
         data = P1_spearman, 
         method = "pearson", 
         continuity = F, 
         conf.level = 0.95)

#P2
P2_spearman <- P2[,2:3]
P2_spearman$total <- P2_spearman[,1] + P2_spearman[,2]
P2_spearman <- P2_spearman[(P2_spearman$total != 0.0),]

cor.test(~ P_P2 + B_P2,
         data = P2_spearman, 
         method = "pearson", 
         continuity = F, 
         conf.level = 0.95)
