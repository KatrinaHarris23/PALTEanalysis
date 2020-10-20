
library(ggplot2)
library(ggmuller)

#B1

population <- read.table("/Users/katrina/Dropbox/My Mac (Katrinas-MBP-3)/Desktop/PALTE_final/draft_revision_work/B1.1/tables/B1_200927.populations.tsv", header=TRUE)
edges <- read.table("/Users/katrina/Dropbox/My Mac (Katrinas-MBP-3)/Desktop/PALTE_final/draft_revision_work/B1.1/tables/B1_200927.edges.tsv", header=TRUE)

Muller_df <- get_Muller_df(edges, population)
palette <- c("#FFFFFF","deeppink","#ffca00","deeppink4","#65479e","#acaad1","hotpink","hotpink4","#7e5033","#89bedc","#dbf1d6","violet","#6699ff","#c6c7e1","grey40","#cc33ff","#33ccff","#908dc2","#0b7734","grey50","#dfdfed","#ff5c00","indianred1","indianred4","#3f2819","#dbe9f6","#73c476","#539ecd","grey60","lightpink","#fc9f65","#0b559f","#790000","lightpink4","#9966ff","#ffff5a","#51228d","#bad6eb","#ea0000","#2b7bba","#f1eff6","#37a055","grey30","#796eb2","#aedea7","#bd784c","plum","plum4")

ggplot(Muller_df, aes_string(x = "Generation", y = "Frequency", group = "Group_id", fill = "Identity", colour = "Identity")) +
  geom_area() +
  theme(legend.position = "none") +
  guides(linetype = FALSE, color = FALSE) +
  scale_y_continuous(labels = 25 * (0:4), name = "Percentage", expand=c(0,0)) +
  scale_x_continuous(expand=c(0,0)) +
  scale_fill_manual(name = "Identity", values = palette) +
  scale_color_manual(values = palette)


#B2

population <- read.table("/Users/katrina/Dropbox/My Mac (Katrinas-MBP-3)/Desktop/PALTE_final/draft_revision_work/B2.2/tables/B2_200927.populations.tsv", header=TRUE)
edges <- read.table("/Users/katrina/Dropbox/My Mac (Katrinas-MBP-3)/Desktop/PALTE_final/draft_revision_work/B2.2/tables/B2_200927.edges.tsv", header=TRUE)

Muller_df <- get_Muller_df(edges, population)
palette <- c("#FFFFFF","#73c476","#fdd4c2","#b11218","#ffff07","#c7e9c0","#ff5c00","#4f3220","#8683bd","#abd0e6","#b6b6d8","#1b69af","#228a44","#3787c0","#e32f27","#fca082","#e1edf8","#686868","#c6c6c6","#82bbdb","#58a1cf","#b30000","#ccdff1","#ed9660","#9e6440","#61409b","#e2e2ef","#fb694a","#084d96")

ggplot(Muller_df, aes_string(x = "Generation", y = "Frequency", group = "Group_id", fill = "Identity", colour = "Identity")) +
  geom_area() +
  theme(legend.position = "none") +
  guides(linetype = FALSE, color = FALSE) +
  scale_y_continuous(labels = 25 * (0:4), name = "Percentage", expand=c(0,0)) +
  scale_x_continuous(expand=c(0,0)) +
  scale_fill_manual(name = "Identity", values = palette) +
  scale_color_manual(values = palette)



#B3
population <- read.table("/Users/katrina/Dropbox/My Mac (Katrinas-MBP-3)/Desktop/PALTE_final/draft_revision_work/B3.2/tables/B3_200927.populations.tsv", header=TRUE)
edges <- read.table("/Users/katrina/Dropbox/My Mac (Katrinas-MBP-3)/Desktop/PALTE_final/draft_revision_work/B3.2/tables/B3_200927.edges.tsv", header=TRUE)

Muller_df <- get_Muller_df(edges, population)
palette <- c("#FFFFFF", "grey23", "grey25", "grey27", "grey29", "grey31", "grey33", "blue", "grey37", "grey39", "grey41", "grey43", "grey45", "grey47", "grey49", "green", "grey35", "grey51", "grey53", "grey55", "grey57", "grey59", "grey61", "grey63", "grey65", "red", "grey67", "grey69", "grey71", "grey73", "grey75")

ggplot(Muller_df, aes_string(x = "Generation", y = "Frequency", group = "Group_id", fill = "Identity", colour = "Identity")) +
  geom_area() +
  theme(legend.position = "none") +
  guides(linetype = FALSE, color = FALSE) +
  scale_y_continuous(labels = 25 * (0:4), name = "Percentage", expand=c(0,0)) +
  scale_x_continuous(expand=c(0,0)) +
  scale_fill_manual(name = "Identity", values = palette) +
  scale_color_manual(values = palette)

#P1
population <- read.table("/Users/katrina/Dropbox/My Mac (Katrinas-MBP-3)/Desktop/PALTE_final/draft_revision_work/P1.1/tables/P1_200927.populations.tsv", header=TRUE)
edges <- read.table("/Users/katrina/Dropbox/My Mac (Katrinas-MBP-3)/Desktop/PALTE_final/draft_revision_work/P1.1/tables/P1_200927.edges.tsv", header=TRUE)

Muller_df <- get_Muller_df(edges, population)
palette <- c("#FFFFFF","grey20","grey23","grey25","grey27","grey29","grey31","grey33","grey35","grey37","grey39","grey41","grey43","grey45","grey47","grey49","grey51","grey53","grey55","grey57","grey59","grey61","grey63","grey65","grey67","grey69","grey71","grey73","grey75","grey77","grey79","grey81","grey83")

ggplot(Muller_df, aes_string(x = "Generation", y = "Frequency", group = "Group_id", fill = "Identity", colour = "Identity")) +
  geom_area() +
  theme(legend.position = "none") +
  guides(linetype = FALSE, color = FALSE) +
  scale_y_continuous(labels = 25 * (0:4), name = "Percentage", expand=c(0,0)) +
  scale_x_continuous(expand=c(0,0)) +
  scale_fill_manual(name = "Identity", values = palette) +
  scale_color_manual(values = palette)



#P2

population <- read.table("/Users/katrina/Dropbox/My Mac (Katrinas-MBP-3)/Desktop/PALTE_final/draft_revision_work/P2.1/tables/P2_200927.populations.tsv", header=TRUE)
edges <- read.table("/Users/katrina/Dropbox/My Mac (Katrinas-MBP-3)/Desktop/PALTE_final/draft_revision_work/P2.1/tables/P2_200927.edges.tsv", header=TRUE)

Muller_df <- get_Muller_df(edges, population)
palette <- c("#FFFFFF","grey20","grey23","grey25","grey27","grey29","grey31","grey33","grey35","grey37","grey39","grey41","grey43","grey45","grey47","grey49","grey51","grey53","grey55","grey57","grey59","grey61","grey63","grey65","grey67","grey69","grey71","grey73","grey75","grey77","grey79","grey81","grey83")

ggplot(Muller_df, aes_string(x = "Generation", y = "Frequency", group = "Group_id", fill = "Identity", colour = "Identity")) +
  geom_area() +
  theme(legend.position = "none") +
  guides(linetype = FALSE, color = FALSE) +
  scale_y_continuous(labels = 25 * (0:4), name = "Percentage", expand=c(0,0)) +
  scale_x_continuous(expand=c(0,0)) +
  scale_fill_manual(name = "Identity", values = palette) +
  scale_color_manual(values = palette)


#P3

population <- read.table("/Users/katrina/Dropbox/My Mac (Katrinas-MBP-3)/Desktop/PALTE_final/draft_revision_work/P3.1/tables/P3_200927.populations.tsv", header=TRUE)
edges <- read.table("/Users/katrina/Dropbox/My Mac (Katrinas-MBP-3)/Desktop/PALTE_final/draft_revision_work/P3.1/tables/P3_200927.edges.tsv", header=TRUE)

Muller_df <- get_Muller_df(edges, population)
palette <- c("#FFFFFF","grey15","grey17","grey19","grey21","#2070b4","grey23","grey25","grey27","green","grey29","grey31","red","grey33","grey35","grey37","grey39","grey41","grey43","grey45","grey47","grey49","grey51","grey53","grey55","grey57","grey59","grey61","grey63","grey65","grey67")

ggplot(Muller_df, aes_string(x = "Generation", y = "Frequency", group = "Group_id", fill = "Identity", colour = "Identity")) +
  geom_area() +
  theme(legend.position = "none") +
  guides(linetype = FALSE, color = FALSE) +
  scale_y_continuous(labels = 25 * (0:4), name = "Percentage", expand=c(0,0)) +
  scale_x_continuous(expand=c(0,0)) +
  scale_fill_manual(name = "Identity", values = palette) +
  scale_color_manual(values = palette)
