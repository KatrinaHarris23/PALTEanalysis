library(ggmuller)
library(ggplot2)
setwd("/Users/katrina/Desktop/working/ggmuller")
theme_set(theme_bw())
#edges <- read.csv("edges_try.csv", header=TRUE)
#View(edges)

#pop <- read.csv("pop_try.csv", header=TRUE)
#View(pop)

#Muller_try <- get_Muller_df(edges, pop, threshold = .005)
#Muller_plot(Muller_try)

#the above works as a test case so now I need to use my data from matLab to create these.

#my test is population P3, just because that was the data set that matlab had up when I started.
#to build the muller diagrams in matlab you need 3 data sets: frequencies, nests, and timepoints.
#gg muller combines these into two. one they call edges and one they call pops. edges denotes the ancestry of the lineages. this is matlab's nests, but simplified. Matlab stores the whole lineage whereas gg muller only requires the most recent ancestry. I believe that I am interpreting these correctly, but we will see when the figure pops out.

#p3_edges <- read.csv("edges_P3.csv", header = TRUE)
#View(p3_edges)

#p3_pop <- read.csv("pop_p3.csv", header= TRUE)
#View(p3_pop)
#p3_muller <- get_Muller_df(p3_edges, p3_pop, threshold = .0005)
#Muller_plot(p3_muller)

###########################################################
#this section is the trials of nailing down the propper method to make these graphs with.
# p1 population was used because it appears to be the most straightforward poulation.
#p1_edges <-  read.csv("edges_p1.csv", header = TRUE)
#p1_pop <- read.csv("pop_p1.csv", header = TRUE)

#View(p1_edges)
#View(p1_pop)
#p1_muller <- get_Muller_df(p1_edges, p1_pop, threshold= .005)
#Muller_plot(p1_muller)

#p1_pop2 <- read.csv("pop_p1_2.csv", header = TRUE)
#View(p1_pop2)
#p1_muller_2 <- get_Muller_df(p1_edges, p1_pop2, threshold = .005)
#Muller_plot(p1_muller_2)

#p1_pop3 <- read.csv("pop_p1_3.csv", header = TRUE)
#View(p1_pop3)
#p1_muller_3 <- get_Muller_df(p1_edges, p1_pop3, threshold = .005)
#Muller_plot(p1_muller_3)

#p1_edges2 <- read.csv("edges_p1_2.csv", header = TRUE)
#p1_muller_4 <- get_Muller_df(p1_edges2, p1_pop2, threshold = .005)
#Muller_plot(p1_muller_4)

#p1_muller_5 <- get_Muller_df(p1_edges2, p1_pop, threshold = .005)
#Muller_plot(p1_muller_5)


#trying the p1 population without the three cohorts that I don't think can be real.

#p1_editededges <- read.csv("edges_P1_edited.csv", header = TRUE)
#p1_popedited <- read.csv("p1_edited_pop.csv", header = TRUE)
#View(p1_editededges)
#View(p1_popedited)
#Mullerp1 <- get_Muller_df(p1_editededges, p1_popedited, threshold = .005)
#Muller_plot(Mullerp1)

#p1_editededges2 <- read.csv("edges_P1_edited2.csv", header = TRUE)
#Mullerp2 <- get_Muller_df(p1_editededges2, p1_popedited, threshold = .005)
#Muller_plot(Mullerp2)
##################################################################

#P1 population final
p1_editededges2 <- read.csv("edges_P1_edited2.csv", header = TRUE)
p1_editedpop2 <- read.csv("p1_edited_pop2.csv", header = TRUE)
Mullerp3 <- get_Muller_df(p1_editededges2, p1_editedpop2, threshold = .005)
Muller_plot(Mullerp3)


#the last figure is the one that I decided is the best one for this use. It does require some manual manipulation. From the three variables produced in MatLab (Frequency, nests, and timepoints) I have manually created two data sets for ggMuller called edges and pop.

#Edges is used to determine which lineage comes from what previous lineage, with 0 being the ancester, and then the subsequent numbers correlating to the row of frequencies in the Frequency variable in Matlab. These ended up being manualy curated by me, using matLab's cohorts. This is a flaw, it meant that I can't just go from the matlab nests variable, but it is better than going from scratch. note: The columns need to be named in a specific way for both of these

#pop is a table that consists of the frequency each cohort reaches at every given time point. there are 3 columns, the first is simply the time points for all of the cohorts, the second is the identity of each cohort, and the third is the actual frequency. Note: again, ggmuller needs these columns to be named very specifically.


# now to attempt to use this method for the rest of the populations.

#This will be the B1 population
#read in the csv files for the two data frames
B1_edges <- read.csv("/Users/katrina/Desktop/working/B1/B1_edges.csv")
#View(B1_edges)
B1_pop <- read.csv("/Users/katrina/Desktop/working/B1/B1_pop.csv")
#View(B1_nests)


MullerB1 <- get_Muller_df(B1_edges, B1_pop, threshold = .005)
Muller_plot(MullerB1)

#note that the ggmuller ancestry differs from the matlab variables in that I took out clusters 2,3, and 4

#second edition just has the last 0 value switched to .0001 for each cluster

B1_pop2 <- read.csv("/Users/katrina/Desktop/working/B1/B1_pop.csv")
MullerB1_2 <- get_Muller_df(B1_edges, B1_pop2, threshold = .005)
Muller_plot(MullerB1_2)

#this is much better and is the one that is going to be used as the "final" for the moment until I can analyze it in relation to the frequency plots to make sure it makes sense. My guess is that I am going to have to change the nesting order, which has not been done here.

#Now for B2
B2_edges <- read.csv("/Users/katrina/Desktop/working/B2/B2_edges.csv")
#View(B2_edges)
#B2_pop <- read.csv("/Users/katrina/Desktop/working/B2/B2_pop.csv")
#View(B2_pop)

#B2_muller <- get_Muller_df(B2_edges, B2_pop, threshold = .005)
#Muller_plot(B2_muller)

#now to change the last zero value to .0001

B2_pop2 <- read.csv ("/Users/katrina/Desktop/working/B2/B2_pop.csv")
B2_muller2 <- get_Muller_df(B2_edges, B2_pop2, threshold = .005)
Muller_plot(B2_muller2)



#the B3 population -- Ann wants this one

B3_edges <- read.csv("/Users/katrina/Desktop/working/B3/B3_edges.csv")
#View(B3_edges)
B3_pop <- read.csv("/Users/katrina/Desktop/working/B3/B3_pop.csv")
#View(B3_pop)

B3_muller <- get_Muller_df(B3_edges, B3_pop, threshold = .005)
Muller_plot(B3_muller)



#Now the P2 population
p2_edges <- read.csv("/Users/katrina/Desktop/working/P2/P2_edges.csv")
#View(p2_edges)
p2_pop <- read.csv("/Users/katrina/Desktop/working/P2/p2_pop.csv")
#View(p2_pop)

p2_muller <- get_Muller_df(p2_edges, p2_pop, threshold = .005)
Muller_plot(p2_muller)


#p3 population
p3_edges <- read.csv("/Users/katrina/Desktop/working/P3/P3_edges.csv")
#View(p3_edges)
P3_pop <- read.csv("/Users/katrina/Desktop/working/P3/P3_pop.csv")
#View(P3_pop)

p3_muller <- get_Muller_df(p3_edges, P3_pop, threshold = .005) #this annoyingly doesn't work
