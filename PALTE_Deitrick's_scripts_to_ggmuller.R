library(ggmuller)
theme_set(theme_bw())

python muller.py -i P1_Muller.xlsx -o /Users/katrina/Desktop/working/get_genotypes/muller_diagrams-master/muller_diagrams/muller/output

python muller.py -i P2_Muller.xlsx -o /Users/katrina/Desktop/working/get_genotypes/muller_diagrams-master/muller_diagrams/muller/output

python muller.py -i P3_Muller.xlsx -o /Users/katrina/Desktop/working/get_genotypes/muller_diagrams-master/muller_diagrams/muller/output

python muller.py -i B1_Muller.xlsx -o /Users/katrina/Desktop/working/get_genotypes/muller_diagrams-master/muller_diagrams/muller/output

python muller.py -i B2_Muller.xlsx -o /Users/katrina/Desktop/working/get_genotypes/muller_diagrams-master/muller_diagrams/muller/output

python muller.py -i B3_Muller.xlsx -o /Users/katrina/Desktop/working/get_genotypes/muller_diagrams-master/muller_diagrams/muller/output

P1_pop <- read.csv("/Users/katrina/Desktop/working/get_genotypes/muller_diagrams-master/muller_diagrams/muller/output/P1_Muller.ggmuller_populations.csv")
P1_edges <- read.csv("/Users/katrina/Desktop/working/get_genotypes/muller_diagrams-master/muller_diagrams/muller/output/P1_Muller.ggmuller_edges.csv")

P2_pop <- read.csv("/Users/katrina/Desktop/working/get_genotypes/muller_diagrams-master/muller_diagrams/muller/output/P2_Muller.ggmuller_populations.csv")
P2_edges <- read.csv("/Users/katrina/Desktop/working/get_genotypes/muller_diagrams-master/muller_diagrams/muller/output/P2_Muller.ggmuller_edges.csv")

P3_pop <- read.csv("/Users/katrina/Desktop/working/get_genotypes/muller_diagrams-master/muller_diagrams/muller/output/P3_Muller.ggmuller_populations.csv")
P3_edges <- read.csv("/Users/katrina/Desktop/working/get_genotypes/muller_diagrams-master/muller_diagrams/muller/output/P3_Muller.ggmuller_edges.csv")

B1_pop <- read.csv("/Users/katrina/Desktop/working/get_genotypes/muller_diagrams-master/muller_diagrams/muller/output/B1_Muller.ggmuller_populations.csv")
B1_edges <- read.csv("/Users/katrina/Desktop/working/get_genotypes/muller_diagrams-master/muller_diagrams/muller/output/B1_Muller.ggmuller_edges.csv")

B2_pop <- read.csv("/Users/katrina/Desktop/working/get_genotypes/muller_diagrams-master/muller_diagrams/muller/output/B2_Muller.ggmuller_populations.csv")
B2_edges <- read.csv("/Users/katrina/Desktop/working/get_genotypes/muller_diagrams-master/muller_diagrams/muller/output/B2_Muller.ggmuller_edges.csv")

B3_pop <- read.csv("/Users/katrina/Desktop/working/get_genotypes/muller_diagrams-master/muller_diagrams/muller/output/B3_Muller.ggmuller_populations.csv")
B3_edges <- read.csv("/Users/katrina/Desktop/working/get_genotypes/muller_diagrams-master/muller_diagrams/muller/output/B3_Muller.ggmuller_edges.csv")


#use ggmuller to make a muller plot
MullerP1 <- get_Muller_df(P1_edges, P1_pop, threshold = .005)
MullerP2 <- get_Muller_df(P2_edges, P2_pop, threshold = .005)
MullerP3 <- get_Muller_df(P3_edges, P3_pop, threshold = .005)
MullerB1 <- get_Muller_df(B1_edges, B1_pop, threshold = .005)
MullerB2 <- get_Muller_df(B2_edges, B2_pop, threshold = .005)
MullerB3 <- get_Muller_df(B3_edges, B3_pop, threshold = .005)


Muller_plot(MullerB1, add_legend = TRUE, xlab = "Time (Days)", ylab = "Frequency") + ggtitle("B1")
Muller_plot(MullerP1, add_legend = TRUE, xlab = "Time (Days)", ylab = "Frequency")+ggtitle("P1")
Muller_plot(MullerB2, add_legend = TRUE, xlab = "Time (Days)", ylab = "Frequency")+ggtitle("B2")
Muller_plot(MullerP2, add_legend = TRUE, xlab = "Time (Days)", ylab = "Frequency")+ggtitle("P2")
Muller_plot(MullerB3, add_legend = TRUE, xlab = "Time (Days)", ylab = "Frequency")+ggtitle("B3")
Muller_plot(MullerP3, add_legend = TRUE, xlab = "Time (Days)", ylab = "Frequency")+ggtitle("P3")


