#trying to find mutations that are between 5 and 10% frequency

#already have tables of mutations filtered out so going to load those

B1_filteredout <- read.csv("/Users/katrina/Desktop/PALTE_final/working/B1/B1_Filtered_out.csv", stringsAsFactors = F)
B2_filteredout <- read.csv("/Users/katrina/Desktop/PALTE_final/working/B2/B2_Filtered_out.csv", stringsAsFactors = F)
B3_filteredout <- read.csv("/Users/katrina/Desktop/PALTE_final/working/B3/B3_Filtered_out.csv", stringsAsFactors = F)

P1_filteredout <- read.csv("/Users/katrina/Desktop/PALTE_final/working/P1/P1_Filtered_out.csv", stringsAsFactors = F)
P2_filteredout <- read.csv("/Users/katrina/Desktop/PALTE_final/working/P2/P2_Filtered_out.csv", stringsAsFactors = F)
P3_filteredout <- read.csv("/Users/katrina/Desktop/PALTE_final/working/P3/P3_Filtered_out.csv", stringsAsFactors = F)

#add sample names to all of them

B1_filteredout$Population <- "B1"
B2_filteredout$Population <- "B2"
B3_filteredout$Population <- "B3"
P1_filteredout$Population <- "P1"
P2_filteredout$Population <- "P2"
P3_filteredout$Population <- "P3"

P1_filteredout$X25 <- "NA"
P1_filteredout$X75 <- "NA"
P2_filteredout$X25 <- "NA"
P2_filteredout$X75 <- "NA"
P3_filteredout$X25 <- "NA"
P3_filteredout$X75 <- "NA"

order <- colnames(B1_filteredout)
setcolorder(P1_filteredout, order)
setcolorder(P2_filteredout, order)
setcolorder(P3_filteredout, order)


all_filtered_out <- rbind(B1_filteredout, B2_filteredout, B3_filteredout, P1_filteredout, P2_filteredout, P3_filteredout)

write.csv(all_filtered_out, "/Users/katrina/Desktop/PALTE_final/draft_revision_work/all_filtered_out.csv")


nrow(all_filtered_out)#11341

#filtering out all mutations that are only observed at one time point in the population
count_filter <- all_filtered_out[!(all_filtered_out$count == 1),]
nrow(count_filter) #4579 mutations left that are observed at more than one time point

#filtering out all mutations that are seen at over 90% frequency at the first sampled time point
beginning_filter <- count_filter[(count_filter$X17 <91),]
nrow(beginning_filter) #4497 mutations left

#filtering for mutations that do not get a sum of 10% frequency
frequency_filter <- beginning_filter[(beginning_filter$Sums >= 10),]
nrow(frequency_filter)

#unpack the information from the description
split <- frequency_filter$X
split = transform(frequency_filter, info =colsplit(frequency_filter$X,';;', names = c('desc_gene_annot','position', 'Mutation')))

split2 = transform(split, info = colsplit(split$info.desc_gene_annot, "::", names= c("description", "gene", "annotation")))


#figure out which positions are observed in all populations, they need to be filtered out.
six <- table(split2$info.position)
write.csv(six, "/Users/katrina/Desktop/PALTE_final/draft_revision_work/six.csv")
six <- read.csv("/Users/katrina/Desktop/PALTE_final/draft_revision_work/six.csv", stringsAsFactors = F) #i manually removed all positions that were observed less than 6 times. 

filter_parallelism <- split2[!(split2$info.position %in% six$Var1),]  #remove mutations in which the position was seen mutated six or more times. this indicates that the mutation is observed in all 6 populations, which should not be possible. 

nrow(filter_parallelism)#2956

write.csv(filter_parallelism, "/Users/katrina/Desktop/PALTE_final/draft_revision_work/5_10_cutoff.csv")

#filter out some descriptions that I know are not real from past experience
filtered <- filter_parallelism[!(filter_parallelism$info.annotation == "16S rRNA (cytosine(967)‑C(5))‑methyltransferase RsmB"),]
filtered <- filtered[!(filtered$info.annotation == "23S ribosomal RNA"),]
filtered <- filtered[!(filtered$info.gene == "PA14_RS17985</>acnA"),]
filtered <- filtered[!(filtered$info.gene == "PA14_RS21975>/>acnD"),]
filtered <- filtered[!(filtered$info.gene == "PA14_RS09995>"),]
filtered <- filtered[!(filtered$info.gene == "PA14_RS15685>/<PA14_RS15690"),]
filtered <- filtered[!(filtered$info.gene == "PA14_RS25875>"),]
filtered <- filtered[!(filtered$info.gene == "PA14_RS26645</<PA14_RS26650"),]
filtered <- filtered[!(filtered$info.gene == "PA14_RS07360<"),]
filtered <- filtered[!(filtered$info.gene == "bioF>"),]

write.csv(filtered, "/Users/katrina/Desktop/PALTE_final/draft_revision_work/5_10_cutoff.csv")


#want to filter that you need to observe it at 2 consecutive time points. did this by hand. now to import it
consecutive_filter <- read.csv("/Users/katrina/Desktop/PALTE_final/draft_revision_work/consecutive_filter.csv", stringsAsFactors = F)
table(consecutive_filter$Population)

consecutive_filter <- consecutive_filter[!(consecutive_filter$info.annotation == "type VI secretion system tip protein VgrG"),] #forgot to remove vgrG. of which there are 8, so there are repeat region

nrow(consecutive_filter) #2442




#

#import the B117 data
B117 <- read.csv("/Users/katrina/Desktop/LTE_output/B117/output_reads.csv", stringsAsFactors = F, header = F)
B117_split <- transform(B117, info = colsplit(B117$V9, "=", names = c("header", "reads")))
B117_split2 <- transform(B117_split, info = colsplit(B117_split$info.reads, "/", names = c("forward", "reverse")))
#now only keep the ones that have at least 3 reads for forward and backwards
B117_fwd <- B117_split2[(B117_split2$info.forward >2),]
B117_rev <- B117_fwd[(B117_fwd$info.reverse >2),]

#import the B125 data
B125 <- read.csv("/Users/katrina/Desktop/LTE_output/B125/output_reads.csv", stringsAsFactors = F, header = F)
B125_split <- transform(B125, info = colsplit(B125$V9, "=", names = c("header", "reads")))
B125_split2 <- transform(B125_split, info = colsplit(B125_split$info.reads, "/", names = c("forward", "reverse")))
#now only keep the ones that have at least 3 reads for forward and backwards
B125_fwd <- B125_split2[(B125_split2$info.forward >2),]
B125_rev <- B125_fwd[(B125_fwd$info.reverse >2),]

#import the B144 data
B144 <- read.csv("/Users/katrina/Desktop/LTE_output/B144/output_reads.csv", stringsAsFactors = F, header = F)
B144_split <- transform(B144, info = colsplit(B144$V9, "=", names = c("header", "reads")))
B144_split2 <- transform(B144_split, info = colsplit(B144_split$info.reads, "/", names = c("forward", "reverse")))
#now only keep the ones that have at least 3 reads for forward and backwards
B144_fwd <- B144_split2[(B144_split2$info.forward >2),]
B144_rev <- B144_fwd[(B144_fwd$info.reverse >2),]

#import the B166 data
B166 <- read.csv("/Users/katrina/Desktop/LTE_output/B166/output_reads.csv", stringsAsFactors = F, header = F)
B166_split <- transform(B166, info = colsplit(B166$V9, "=", names = c("header", "reads")))
B166_split2 <- transform(B166_split, info = colsplit(B166_split$info.reads, "/", names = c("forward", "reverse")))
#now only keep the ones that have at least 3 reads for forward and backwards
B166_fwd <- B166_split2[(B166_split2$info.forward >2),]
B166_rev <- B166_fwd[(B166_fwd$info.reverse >2),]

#import the B175 data
B175 <- read.csv("/Users/katrina/Desktop/LTE_output/B175/output_reads.csv", stringsAsFactors = F, header = F)
B175_split <- transform(B175, info = colsplit(B175$V9, "=", names = c("header", "reads")))
B175_split2 <- transform(B175_split, info = colsplit(B175_split$info.reads, "/", names = c("forward", "reverse")))
#now only keep the ones that have at least 3 reads for forward and backwards
B175_fwd <- B175_split2[(B175_split2$info.forward >2),]
B175_rev <- B175_fwd[(B175_fwd$info.reverse >2),]

#import the B190 data
B190 <- read.csv("/Users/katrina/Desktop/LTE_output/B190/output_reads.csv", stringsAsFactors = F, header = F)
B190_split <- transform(B190, info = colsplit(B190$V9, "=", names = c("header", "reads")))
B190_split2 <- transform(B190_split, info = colsplit(B190_split$info.reads, "/", names = c("forward", "reverse")))
#now only keep the ones that have at least 3 reads for forward and backwards
B190_fwd <- B190_split2[(B190_split2$info.forward >2),]
B190_rev <- B190_fwd[(B190_fwd$info.reverse >2),]


####

#import the B217 data
B217 <- read.csv("/Users/katrina/Desktop/LTE_output/B217/output_reads.csv", stringsAsFactors = F, header = F)
B217_split <- transform(B217, info = colsplit(B217$V9, "=", names = c("header", "reads")))
B217_split2 <- transform(B217_split, info = colsplit(B217_split$info.reads, "/", names = c("forward", "reverse")))
#now only keep the ones that have at least 3 reads for forward and backwards
B217_fwd <- B217_split2[(B217_split2$info.forward >2),]
B217_rev <- B217_fwd[(B217_fwd$info.reverse >2),]

#import the B225 data
B225 <- read.csv("/Users/katrina/Desktop/LTE_output/B225/output_reads.csv", stringsAsFactors = F, header = F)
B225_split <- transform(B225, info = colsplit(B225$V9, "=", names = c("header", "reads")))
B225_split2 <- transform(B225_split, info = colsplit(B225_split$info.reads, "/", names = c("forward", "reverse")))
#now only keep the ones that have at least 3 reads for forward and backwards
B225_fwd <- B225_split2[(B225_split2$info.forward >2),]
B225_rev <- B225_fwd[(B225_fwd$info.reverse >2),]

#import the B244 data
B244 <- read.csv("/Users/katrina/Desktop/LTE_output/B244/output_reads.csv", stringsAsFactors = F, header = F)
B244_split <- transform(B244, info = colsplit(B244$V9, "=", names = c("header", "reads")))
B244_split2 <- transform(B244_split, info = colsplit(B244_split$info.reads, "/", names = c("forward", "reverse")))
#now only keep the ones that have at least 3 reads for forward and backwards
B244_fwd <- B244_split2[(B244_split2$info.forward >2),]
B244_rev <- B244_fwd[(B244_fwd$info.reverse >2),]

#import the B266 data
B266 <- read.csv("/Users/katrina/Desktop/LTE_output/B266/output_reads.csv", stringsAsFactors = F, header = F)
B266_split <- transform(B266, info = colsplit(B266$V9, "=", names = c("header", "reads")))
B266_split2 <- transform(B266_split, info = colsplit(B266_split$info.reads, "/", names = c("forward", "reverse")))
#now only keep the ones that have at least 3 reads for forward and backwards
B266_fwd <- B266_split2[(B266_split2$info.forward >2),]
B266_rev <- B266_fwd[(B266_fwd$info.reverse >2),]

#import the B275 data
#B275 <- read.csv("/Users/katrina/Desktop/LTE_output/B275/output_reads.csv", stringsAsFactors = F, header = F)
#B275_split <- transform(B275, info = colsplit(B275$V9, "=", names = c("header", "reads")))
#B275_split2 <- transform(B275_split, info = colsplit(B275_split$info.reads, "/", names = c("forward", "reverse")))
#now only keep the ones that have at least 3 reads for forward and backwards
#B275_fwd <- B275_split2[(B275_split2$info.forward >2),]
#B275_rev <- B275_fwd[(B275_fwd$info.reverse >2),]

#import the B290 data
B290 <- read.csv("/Users/katrina/Desktop/LTE_output/B290/output_reads.csv", stringsAsFactors = F, header = F)
B290_split <- transform(B290, info = colsplit(B290$V9, "=", names = c("header", "reads")))
B290_split2 <- transform(B290_split, info = colsplit(B290_split$info.reads, "/", names = c("forward", "reverse")))
#now only keep the ones that have at least 3 reads for forward and backwards
B290_fwd <- B290_split2[(B290_split2$info.forward >2),]
B290_rev <- B290_fwd[(B290_fwd$info.reverse >2),]






#import the B317 data
B317 <- read.csv("/Users/katrina/Desktop/LTE_output/B317/output_reads.csv", stringsAsFactors = F, header = F)
B317_split <- transform(B317, info = colsplit(B317$V9, "=", names = c("header", "reads")))
B317_split2 <- transform(B317_split, info = colsplit(B317_split$info.reads, "/", names = c("forward", "reverse")))
#now only keep the ones that have at least 3 reads for forward and backwards
B317_fwd <- B317_split2[(B317_split2$info.forward >2),]
B317_rev <- B317_fwd[(B317_fwd$info.reverse >2),]

#import the B325 data
B325 <- read.csv("/Users/katrina/Desktop/LTE_output/B325/output_reads.csv", stringsAsFactors = F, header = F)
B325_split <- transform(B325, info = colsplit(B325$V9, "=", names = c("header", "reads")))
B325_split2 <- transform(B325_split, info = colsplit(B325_split$info.reads, "/", names = c("forward", "reverse")))
#now only keep the ones that have at least 3 reads for forward and backwards
B325_fwd <- B325_split2[(B325_split2$info.forward >2),]
B325_rev <- B325_fwd[(B325_fwd$info.reverse >2),]

#import the B344 data
B344 <- read.csv("/Users/katrina/Desktop/LTE_output/B344/output_reads.csv", stringsAsFactors = F, header = F)
B344_split <- transform(B344, info = colsplit(B344$V9, "=", names = c("header", "reads")))
B344_split2 <- transform(B344_split, info = colsplit(B344_split$info.reads, "/", names = c("forward", "reverse")))
#now only keep the ones that have at least 3 reads for forward and backwards
B344_fwd <- B344_split2[(B344_split2$info.forward >2),]
B344_rev <- B344_fwd[(B344_fwd$info.reverse >2),]

#import the B366 data
B366 <- read.csv("/Users/katrina/Desktop/LTE_output/B366/output_reads.csv", stringsAsFactors = F, header = F)
B366_split <- transform(B366, info = colsplit(B366$V9, "=", names = c("header", "reads")))
B366_split2 <- transform(B366_split, info = colsplit(B366_split$info.reads, "/", names = c("forward", "reverse")))
#now only keep the ones that have at least 3 reads for forward and backwards
B366_fwd <- B366_split2[(B366_split2$info.forward >2),]
B366_rev <- B366_fwd[(B366_fwd$info.reverse >2),]

#import the B375 data
B375 <- read.csv("/Users/katrina/Desktop/LTE_output/B375/output_reads.csv", stringsAsFactors = F, header = F)
B375_split <- transform(B375, info = colsplit(B375$V9, "=", names = c("header", "reads")))
B375_split2 <- transform(B375_split, info = colsplit(B375_split$info.reads, "/", names = c("forward", "reverse")))
#now only keep the ones that have at least 3 reads for forward and backwards
B375_fwd <- B375_split2[(B375_split2$info.forward >2),]
B375_rev <- B375_fwd[(B375_fwd$info.reverse >2),]

#import the B390 data
B390 <- read.csv("/Users/katrina/Desktop/LTE_output/B390/output_reads.csv", stringsAsFactors = F, header = F)
B390_split <- transform(B390, info = colsplit(B390$V9, "=", names = c("header", "reads")))
B390_split2 <- transform(B390_split, info = colsplit(B390_split$info.reads, "/", names = c("forward", "reverse")))
#now only keep the ones that have at least 3 reads for forward and backwards
B390_fwd <- B390_split2[(B390_split2$info.forward >2),]
B390_rev <- B390_fwd[(B390_fwd$info.reverse >2),]







#import the P117 data
P117 <- read.csv("/Users/katrina/Desktop/LTE_output/P117/output_reads.csv", stringsAsFactors = F, header = F)
P117_split <- transform(P117, info = colsplit(P117$V9, "=", names = c("header", "reads")))
P117_split2 <- transform(P117_split, info = colsplit(P117_split$info.reads, "/", names = c("forward", "reverse")))
#now only keep the ones that have at least 3 reads for forward and backwards
P117_fwd <- P117_split2[(P117_split2$info.forward >2),]
P117_rev <- P117_fwd[(P117_fwd$info.reverse >2),]

#import the P144 data
P144 <- read.csv("/Users/katrina/Desktop/LTE_output/P144/output_reads.csv", stringsAsFactors = F, header = F)
P144_split <- transform(P144, info = colsplit(P144$V9, "=", names = c("header", "reads")))
P144_split2 <- transform(P144_split, info = colsplit(P144_split$info.reads, "/", names = c("forward", "reverse")))
#now only keep the ones that have at least 3 reads for forward and backwards
P144_fwd <- P144_split2[(P144_split2$info.forward >2),]
P144_rev <- P144_fwd[(P144_fwd$info.reverse >2),]

#import the P166 data
P166 <- read.csv("/Users/katrina/Desktop/LTE_output/P166/output_reads.csv", stringsAsFactors = F, header = F)
P166_split <- transform(P166, info = colsplit(P166$V9, "=", names = c("header", "reads")))
P166_split2 <- transform(P166_split, info = colsplit(P166_split$info.reads, "/", names = c("forward", "reverse")))
#now only keep the ones that have at least 3 reads for forward and backwards
P166_fwd <- P166_split2[(P166_split2$info.forward >2),]
P166_rev <- P166_fwd[(P166_fwd$info.reverse >2),]


#import the P190 data
P190 <- read.csv("/Users/katrina/Desktop/LTE_output/P190/output_reads.csv", stringsAsFactors = F, header = F)
P190_split <- transform(P190, info = colsplit(P190$V9, "=", names = c("header", "reads")))
P190_split2 <- transform(P190_split, info = colsplit(P190_split$info.reads, "/", names = c("forward", "reverse")))
#now only keep the ones that have at least 3 reads for forward and backwards
P190_fwd <- P190_split2[(P190_split2$info.forward >2),]
P190_rev <- P190_fwd[(P190_fwd$info.reverse >2),]







#import the P217 data
P217 <- read.csv("/Users/katrina/Desktop/LTE_output/P217/output_reads.csv", stringsAsFactors = F, header = F)
P217_split <- transform(P217, info = colsplit(P217$V9, "=", names = c("header", "reads")))
P217_split2 <- transform(P217_split, info = colsplit(P217_split$info.reads, "/", names = c("forward", "reverse")))
#now only keep the ones that have at least 3 reads for forward and backwards
P217_fwd <- P217_split2[(P217_split2$info.forward >2),]
P217_rev <- P217_fwd[(P217_fwd$info.reverse >2),]

#import the P244 data
P244 <- read.csv("/Users/katrina/Desktop/LTE_output/P244/output_reads.csv", stringsAsFactors = F, header = F)
P244_split <- transform(P244, info = colsplit(P244$V9, "=", names = c("header", "reads")))
P244_split2 <- transform(P244_split, info = colsplit(P244_split$info.reads, "/", names = c("forward", "reverse")))
#now only keep the ones that have at least 3 reads for forward and backwards
P244_fwd <- P244_split2[(P244_split2$info.forward >2),]
P244_rev <- P244_fwd[(P244_fwd$info.reverse >2),]

#import the P266 data
P266 <- read.csv("/Users/katrina/Desktop/LTE_output/P266/output_reads.csv", stringsAsFactors = F, header = F)
P266_split <- transform(P266, info = colsplit(P266$V9, "=", names = c("header", "reads")))
P266_split2 <- transform(P266_split, info = colsplit(P266_split$info.reads, "/", names = c("forward", "reverse")))
#now only keep the ones that have at least 3 reads for forward and backwards
P266_fwd <- P266_split2[(P266_split2$info.forward >2),]
P266_rev <- P266_fwd[(P266_fwd$info.reverse >2),]


#import the P290 data
P290 <- read.csv("/Users/katrina/Desktop/LTE_output/P290/output_reads.csv", stringsAsFactors = F, header = F)
P290_split <- transform(P290, info = colsplit(P290$V9, "=", names = c("header", "reads")))
P290_split2 <- transform(P290_split, info = colsplit(P290_split$info.reads, "/", names = c("forward", "reverse")))
#now only keep the ones that have at least 3 reads for forward and backwards
P290_fwd <- P290_split2[(P290_split2$info.forward >2),]
P290_rev <- P290_fwd[(P290_fwd$info.reverse >2),]




#import the P317 data
P317 <- read.csv("/Users/katrina/Desktop/LTE_output/P317/output_reads.csv", stringsAsFactors = F, header = F)
P317_split <- transform(P317, info = colsplit(P317$V9, "=", names = c("header", "reads")))
P317_split2 <- transform(P317_split, info = colsplit(P317_split$info.reads, "/", names = c("forward", "reverse")))
#now only keep the ones that have at least 3 reads for forward and backwards
P317_fwd <- P317_split2[(P317_split2$info.forward >2),]
P317_rev <- P317_fwd[(P317_fwd$info.reverse >2),]

#import the P344 data
P344 <- read.csv("/Users/katrina/Desktop/LTE_output/P344/output_reads.csv", stringsAsFactors = F, header = F)
P344_split <- transform(P344, info = colsplit(P344$V9, "=", names = c("header", "reads")))
P344_split2 <- transform(P344_split, info = colsplit(P344_split$info.reads, "/", names = c("forward", "reverse")))
#now only keep the ones that have at least 3 reads for forward and backwards
P344_fwd <- P344_split2[(P344_split2$info.forward >2),]
P344_rev <- P344_fwd[(P344_fwd$info.reverse >2),]

#import the P366 data
P366 <- read.csv("/Users/katrina/Desktop/LTE_output/P366/output_reads.csv", stringsAsFactors = F, header = F)
P366_split <- transform(P366, info = colsplit(P366$V9, "=", names = c("header", "reads")))
P366_split2 <- transform(P366_split, info = colsplit(P366_split$info.reads, "/", names = c("forward", "reverse")))
#now only keep the ones that have at least 3 reads for forward and backwards
P366_fwd <- P366_split2[(P366_split2$info.forward >2),]
P366_rev <- P366_fwd[(P366_fwd$info.reverse >2),]


#import the P390 data
P390 <- read.csv("/Users/katrina/Desktop/LTE_output/P390/output_reads.csv", stringsAsFactors = F, header = F)
P390_split <- transform(P390, info = colsplit(P390$V9, "=", names = c("header", "reads")))
P390_split2 <- transform(P390_split, info = colsplit(P390_split$info.reads, "/", names = c("forward", "reverse")))
#now only keep the ones that have at least 3 reads for forward and backwards
P390_fwd <- P390_split2[(P390_split2$info.forward >2),]
P390_rev <- P390_fwd[(P390_fwd$info.reverse >2),]





#going to do this for all time points and then that is the list of all positions that you can keep

reads_keep <- rbind(B117_rev, B125_rev, B144_rev, B166_rev, B175_rev, B190_rev, B217_rev, B225_rev, B244_rev, B266_rev,  B290_rev, B317_rev, B325_rev, B344_rev, B366_rev, B375_rev, B390_rev, P117_rev, P144_rev, P166_rev, P190_rev, P217_rev, P244_rev, P266_rev, P290_rev, P317_rev, P344_rev, P366_rev, P390_rev) #positions are in V5

#only keep the positions that are found in the above reads_keep file

three_reads <- consecutive_filter[(consecutive_filter$info.position %in% reads_keep$V5),] 
nrow(three_reads) #2359

write.csv(three_reads, "/Users/katrina/Desktop/PALTE_final/draft_revision_work/three_reads.csv")

#I have filtered calls to make sure that there are no more than 2 mutations per population in a 50 bp sliding window. This got rid of a lot more mutations. 

fiftybp_filter <- read.csv("/Users/katrina/Desktop/PALTE_final/draft_revision_work/50bp_filter.csv", stringsAsFactors = F)

table(fiftybp_filter$Population)
nrow(fiftybp_filter) #1134


#I went through one by one to figure out which mutations were in regions of high variation (mapping error). this is reading in the data

final <- read.csv("/Users/katrina/Desktop/PALTE_final/draft_revision_work/final_lowfrequency_mutations.csv", stringsAsFactors = F)
final[,4:10] <- final[,4:10]/100
table(final$Population)
nrow(final)

#these are the final mutations. So I am going to split into the populations to be able to add to the muller plots. 
#read in the previous table of all mutations 10% and above
old <- read.csv("/Users/katrina/Desktop/PALTE_final/draft_revision_work/all_populations.csv", stringsAsFactors = F)

old_modified <- old[, c(2:8, 10, 11, 12, 19, 21)]
new <- final[,c(4:10, 13, 15, 16, 18, 17)]
colnames(new) <- colnames(old_modified)
all <- rbind(new, old_modified)

table(old_modified$Population)
table(all$Population)

#wirte out the large final list of mutations
write.csv(all, "/Users/katrina/Desktop/PALTE_final/draft_revision_work/all_mutations_200930.csv")

#updated the gene names to be the current locus tags so need to read in again
all <- read.csv("/Users/katrina/Desktop/PALTE_final/draft_revision_work/all_mutations_200930.csv")

#split into each population

B1 <- all[(all$Population == "B1"),]
B2 <- all[(all$Population == "B2"),]
B3 <- all[(all$Population == "B3"),]

P1 <- all[(all$Population == "P1"),]
P2 <- all[(all$Population == "P2"),]
P3 <- all[(all$Population == "P3"),]

#write out the files so that I can run them through lolipop
write.csv(B1, "/Users/katrina/Desktop/PALTE_final/draft_revision_work/B1_200927.csv")
write.csv(B2, "/Users/katrina/Desktop/PALTE_final/draft_revision_work/B2_200927.csv")
write.csv(B3, "/Users/katrina/Desktop/PALTE_final/draft_revision_work/B3_200927.csv")
write.csv(P1, "/Users/katrina/Desktop/PALTE_final/draft_revision_work/P1_200927.csv")
write.csv(P2, "/Users/katrina/Desktop/PALTE_final/draft_revision_work/P2_200927.csv")
write.csv(P3, "/Users/katrina/Desktop/PALTE_final/draft_revision_work/P3_200927.csv")

#I have added genotype information on here so will reimport the data
B1 <- read.csv( "/Users/katrina/Desktop/PALTE_final/draft_revision_work/B1_200927.csv", stringsAsFactors = F)
B2 <- read.csv( "/Users/katrina/Desktop/PALTE_final/draft_revision_work/B2_200927.csv", stringsAsFactors = F)
B3 <- read.csv( "/Users/katrina/Desktop/PALTE_final/draft_revision_work/B3_200927.csv", stringsAsFactors = F)
P1 <- read.csv( "/Users/katrina/Desktop/PALTE_final/draft_revision_work/P1_200927.csv", stringsAsFactors = F)
P2 <- read.csv( "/Users/katrina/Desktop/PALTE_final/draft_revision_work/P2_200927.csv", stringsAsFactors = F)
P3 <- read.csv( "/Users/katrina/Desktop/PALTE_final/draft_revision_work/P3_200927.csv", stringsAsFactors = F)

P1_mutations_per_genotype <- as.data.frame(table(P1$genotype))
P2_mutations_per_genotype <- as.data.frame(table(P2$genotype))
P3_mutations_per_genotype <- as.data.frame(table(P3$Genotype))
B1_mutations_per_genotype <- as.data.frame(table(B1$Genotype))
B2_mutations_per_genotype <- as.data.frame(table(B2$Genotype))
B3_mutations_per_genotype <- as.data.frame(table(B3$Genotype))

#print the tables of number of mutations per genotype
write.csv(P1_mutations_per_genotype,"/Users/katrina/Desktop/PALTE_final/draft_revision_work/P1_mutations_per_genotype.csv")
write.csv(P2_mutations_per_genotype,"/Users/katrina/Desktop/PALTE_final/draft_revision_work/P2_mutations_per_genotype.csv")
write.csv(P3_mutations_per_genotype,"/Users/katrina/Desktop/PALTE_final/draft_revision_work/P3_mutations_per_genotype.csv")
write.csv(B1_mutations_per_genotype,"/Users/katrina/Desktop/PALTE_final/draft_revision_work/B1_mutations_per_genotype.csv")
write.csv(B2_mutations_per_genotype,"/Users/katrina/Desktop/PALTE_final/draft_revision_work/B2_mutations_per_genotype.csv")
write.csv(B3_mutations_per_genotype,"/Users/katrina/Desktop/PALTE_final/draft_revision_work/B3_mutations_per_genotype.csv")




#determine mutational dynamics for all populations using the previously made function
setwd("/Users/katrina/Desktop/PALTE_final/draft_revision_work")

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
  
  synonymous <- sum(ifelse(aminos2$first == aminos2$last,1,0)) #if the letters are the same before and after the mutation it is synonymous
  nonsynonymous <- sum(ifelse(aminos2$first == aminos2$last,0,1)) #if the letters are different then it is nonsynonymous
  dnds <- (nonsynonymous/synonymous)/2.60
  
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



Mutations_analysis(all, "Amino.Acid", "Mutation")

#import all files
B1_summary <- read.csv("/Users/katrina/Desktop/PALTE_final/draft_revision_work/B1_summary.csv")
B2_summary <- read.csv("/Users/katrina/Desktop/PALTE_final/draft_revision_work/B2_summary.csv")
B3_summary <- read.csv("/Users/katrina/Desktop/PALTE_final/draft_revision_work/B3_summary.csv")
P1_summary <- read.csv("/Users/katrina/Desktop/PALTE_final/draft_revision_work/P1_summary.csv")
P2_summary <- read.csv("/Users/katrina/Desktop/PALTE_final/draft_revision_work/P2_summary.csv")
P3_summary <- read.csv("/Users/katrina/Desktop/PALTE_final/draft_revision_work/P3_summary.csv")
all_summary <- read.csv("/Users/katrina/Desktop/PALTE_final/draft_revision_work/all_summary.csv")

#combine into one big table
Mutation_breakdown <- cbind(all_summary, B1_summary$V2, B2_summary$V2, B3_summary$V2, P1_summary$V2, P2_summary$V2, P3_summary$V2 )
colnames(Mutation_breakdown) <- c("","Category", "All", "B1","B2","B3","P1","P2","P3")

#write out file
write.csv(Mutation_breakdown, "/Users/katrina/Desktop/PALTE_final/draft_revision_work/Final_breakdown.csv")

#now to calculate diversity
shannon_diversity <- as.data.frame(matrix(ncol = 8, byrow=T, c("Population", "Day 0", "Day 17", "Day 25", "Day 44", "Day 66", "Day 75", "Day 90",
                    "P1", diversity(P1$X0, index="shannon"),
                          diversity(P1$X17, index = "shannon"),
                          diversity(P1$X25, index = "shannon"),
                          diversity(P1$X44, index = "shannon"),
                          diversity(P1$X66, index = "shannon"),
                          diversity(P1$X75, index = "shannon"),
                          diversity(P1$X90, index = "shannon"),
                    "P2", diversity(P2$X0, index="shannon"),
                    diversity(P2$X17, index = "shannon"),
                    diversity(P2$X25, index = "shannon"),
                    diversity(P2$X44, index = "shannon"),
                    diversity(P2$X66, index = "shannon"),
                    diversity(P2$X75, index = "shannon"),
                    diversity(P2$X90, index = "shannon"),
                    "P3", diversity(P3$X0, index="shannon"),
                    diversity(P3$X17, index = "shannon"),
                    diversity(P3$X25, index = "shannon"),
                    diversity(P3$X44, index = "shannon"),
                    diversity(P3$X66, index = "shannon"),
                    diversity(P3$X75, index = "shannon"),
                    diversity(P3$X90, index = "shannon"),
                    "B1", diversity(B1$X0, index="shannon"),
                    diversity(B1$X17, index = "shannon"),
                    diversity(B1$X25, index = "shannon"),
                    diversity(B1$X44, index = "shannon"),
                    diversity(B1$X66, index = "shannon"),
                    diversity(B1$X75, index = "shannon"),
                    diversity(B1$X90, index = "shannon"),
                    "B2", diversity(B2$X0, index="shannon"),
                    diversity(B2$X17, index = "shannon"),
                    diversity(B2$X25, index = "shannon"),
                    diversity(B2$X44, index = "shannon"),
                    diversity(B2$X66, index = "shannon"),
                    diversity(B2$X75, index = "shannon"),
                    diversity(B2$X90, index = "shannon"),
                    "B3", diversity(B3$X0, index="shannon"),
                    diversity(B3$X17, index = "shannon"),
                    diversity(B3$X25, index = "shannon"),
                    diversity(B3$X44, index = "shannon"),
                    diversity(B3$X66, index = "shannon"),
                    diversity(B3$X75, index = "shannon"),
                    diversity(B3$X90, index = "shannon"))))
#write out the diversity table
write.csv(shannon_diversity, "/Users/katrina/Desktop/PALTE_final/draft_revision_work/Shannon_Diversity.csv")



all_17 <- all[,-c(1,3:7)]
all_25 <- all[,-c(1:2,4:7)]
all_44 <- all[,-c(1:3,5:7)]
all_66 <- all[,-c(1:4,6:7)]
all_75 <- all[,-c(1:5,7:7)]
all_90 <- all[,-c(1:6)]

all_17 <- all_17[!(all_17$X17 == 0.0),]
all_25 <- all_25[!(all_25$X25 == 0.0),]
all_44 <- all_44[!(all_44$X44 == 0),]
all_66 <- all_66[!(all_66$X66 == 0),]
all_75 <- all_75[!(all_75$X75 == 0.0),]
all_90 <- all_90[!(all_90$X90 == 0.0),]


#reformat to do a gene level analysis
m_17 <- melt(all_17, id=c("Gene", "Population"), measure.vars = "X17")
cast_17 <- as.data.frame(t(dcast(m_17, Gene~Position, mean, value.var = "value", fill = 0)))
colnames(cast_17) <- as.character(unlist(cast_17[1,]))
a_17 <- cast_17[-1,]



#determine how many mutations in each population at each time point. 
Mutation_numbers <- as.data.frame(matrix(ncol = 8, byrow=T, c("Population", "Day 0", "Day 17", "Day 25", "Day 44", "Day 66", "Day 75", "Day 90",
                     "P1", nrow(P1[!(P1$X0 == 0),]),
                     nrow(P1[!(P1$X17 == 0),]),
                     nrow(P1[!(P1$X25 == 0),]),
                     nrow(P1[!(P1$X44 == 0),]),
                     nrow(P1[!(P1$X66 == 0),]),
                     nrow(P1[!(P1$X75 == 0),]),
                     nrow(P1[!(P1$X90 == 0),]),
                     "P2", nrow(P2[!(P2$X0 == 0),]),
                     nrow(P2[!(P2$X17 == 0),]),
                     nrow(P2[!(P2$X25 == 0),]),
                     nrow(P2[!(P2$X44 == 0),]),
                     nrow(P2[!(P2$X66 == 0),]),
                     nrow(P2[!(P2$X75 == 0),]),
                     nrow(P2[!(P2$X90 == 0),]),
                     "P3", nrow(P3[!(P3$X0 == 0),]),
                     nrow(P3[!(P3$X17 == 0),]),
                     nrow(P3[!(P3$X25 == 0),]),
                     nrow(P3[!(P3$X44 == 0),]),
                     nrow(P3[!(P3$X66 == 0),]),
                     nrow(P3[!(P3$X75 == 0),]),
                     nrow(P3[!(P3$X90 == 0),]),
                     "B1", nrow(B1[!(B1$X0 == 0),]),
                     nrow(B1[!(B1$X17 == 0),]),
                     nrow(B1[!(B1$X25 == 0),]),
                     nrow(B1[!(B1$X44 == 0),]),
                     nrow(B1[!(B1$X66 == 0),]),
                     nrow(B1[!(B1$X75 == 0),]),
                     nrow(B1[!(B1$X90 == 0),]),
                     "B2", nrow(B2[!(B2$X0 == 0),]),
                     nrow(B2[!(B2$X17 == 0),]),
                     nrow(B2[!(B2$X25 == 0),]),
                     nrow(B2[!(B2$X44 == 0),]),
                     nrow(B2[!(B2$X66 == 0),]),
                     nrow(B2[!(B2$X75 == 0),]),
                     nrow(B2[!(B2$X90 == 0),]),
                     "B3", nrow(B3[!(B3$X0 == 0),]),
                     nrow(B3[!(B3$X17 == 0),]),
                     nrow(B3[!(B3$X25 == 0),]),
                     nrow(B3[!(B3$X44 == 0),]),
                     nrow(B3[!(B3$X66 == 0),]),
                     nrow(B3[!(B3$X75 == 0),]),
                     nrow(B3[!(B3$X90 == 0),]))))

#write to file
write.csv(Mutation_numbers, "/Users/katrina/Desktop/PALTE_final/draft_revision_work/Mutations_over_time.csv")
                     


#want to make tables of how many genotypes are present at each time point in each population 
B1_genotypefreqs <- read.csv("/Users/katrina/Desktop/PALTE_final/draft_revision_work/B1_genotypefreqs.csv", stringsAsFactors = F)
B2_genotypefreqs <- read.csv("/Users/katrina/Desktop/PALTE_final/draft_revision_work/B2_genotypefreqs.csv", stringsAsFactors = F)
B3_genotypefreqs <- read.csv("/Users/katrina/Desktop/PALTE_final/draft_revision_work/B3_genotypefreqs.csv", stringsAsFactors = F)
P1_genotypefreqs <- read.csv("/Users/katrina/Desktop/PALTE_final/draft_revision_work/P1_genotypefreqs.csv", stringsAsFactors = F)
P2_genotypefreqs <- read.csv("/Users/katrina/Desktop/PALTE_final/draft_revision_work/P2_genotypefreqs.csv", stringsAsFactors = F)
P3_genotypefreqs <- read.csv("/Users/katrina/Desktop/PALTE_final/draft_revision_work/P3_genotypefreqs.csv", stringsAsFactors = F)

P1_17 <- P1_genotypefreqs[(P1_genotypefreqs$Generation == 17),]
P1_44 <- P1_genotypefreqs[(P1_genotypefreqs$Generation == 44),]
P1_66 <- P1_genotypefreqs[(P1_genotypefreqs$Generation == 66),]
P1_90 <- P1_genotypefreqs[(P1_genotypefreqs$Generation == 90),]

P2_17 <- P2_genotypefreqs[(P2_genotypefreqs$Generation == 17),]
P2_44 <- P2_genotypefreqs[(P2_genotypefreqs$Generation == 44),]
P2_66 <- P2_genotypefreqs[(P2_genotypefreqs$Generation == 66),]
P2_90 <- P2_genotypefreqs[(P2_genotypefreqs$Generation == 90),]

P3_17 <- P3_genotypefreqs[(P3_genotypefreqs$Generation == 17),]
P3_44 <- P3_genotypefreqs[(P3_genotypefreqs$Generation == 44),]
P3_66 <- P3_genotypefreqs[(P3_genotypefreqs$Generation == 66),]
P3_90 <- P3_genotypefreqs[(P3_genotypefreqs$Generation == 90),]


B1_17 <- B1_genotypefreqs[(B1_genotypefreqs$Generation == 17),]
B1_25 <- B1_genotypefreqs[(B1_genotypefreqs$Generation == 25),]
B1_75 <- B1_genotypefreqs[(B1_genotypefreqs$Generation == 75),]
B1_44 <- B1_genotypefreqs[(B1_genotypefreqs$Generation == 44),]
B1_66 <- B1_genotypefreqs[(B1_genotypefreqs$Generation == 66),]
B1_90 <- B1_genotypefreqs[(B1_genotypefreqs$Generation == 90),]

B2_17 <- B2_genotypefreqs[(B2_genotypefreqs$Generation == 17),]
B2_25 <- B2_genotypefreqs[(B2_genotypefreqs$Generation == 25),]
B2_75 <- B2_genotypefreqs[(B2_genotypefreqs$Generation == 75),]
B2_44 <- B2_genotypefreqs[(B2_genotypefreqs$Generation == 44),]
B2_66 <- B2_genotypefreqs[(B2_genotypefreqs$Generation == 66),]
B2_90 <- B2_genotypefreqs[(B2_genotypefreqs$Generation == 90),]

B3_17 <- B3_genotypefreqs[(B3_genotypefreqs$Generation == 17),]
B3_25 <- B3_genotypefreqs[(B3_genotypefreqs$Generation == 25),]
B3_75 <- B3_genotypefreqs[(B3_genotypefreqs$Generation == 75),]
B3_44 <- B3_genotypefreqs[(B3_genotypefreqs$Generation == 44),]
B3_66 <- B3_genotypefreqs[(B3_genotypefreqs$Generation == 66),]
B3_90 <- B3_genotypefreqs[(B3_genotypefreqs$Generation == 90),]


P1_17_keep <- P1_17[!(P1_17$Population == 0),]
P1_44_keep <- P1_44[!(P1_44$Population == 0),]
P1_66_keep <- P1_66[!(P1_66$Population == 0),]
P1_90_keep <- P1_90[!(P1_90$Population == 0),]

P2_17_keep <- P2_17[!(P2_17$Population == 0),]
P2_44_keep <- P2_44[!(P2_44$Population == 0),]
P2_66_keep <- P2_66[!(P2_66$Population == 0),]
P2_90_keep <- P2_90[!(P2_90$Population == 0),]

P3_17_keep <- P3_17[!(P3_17$Population == 0),]
P3_44_keep <- P3_44[!(P3_44$Population == 0),]
P3_66_keep <- P3_66[!(P3_66$Population == 0),]
P3_90_keep <- P3_90[!(P3_90$Population == 0),]


B1_17_keep <- B1_17[!(B1_17$Population == 0),]
B1_25_keep <- B1_25[!(B1_25$Population == 0),]
B1_44_keep <- B1_44[!(B1_44$Population == 0),]
B1_66_keep <- B1_66[!(B1_66$Population == 0),]
B1_75_keep <- B1_75[!(B1_75$Population == 0),]
B1_90_keep <- B1_90[!(B1_90$Population == 0),]

B2_17_keep <- B2_17[!(B2_17$Population == 0),]
B2_25_keep <- B2_25[!(B2_25$Population == 0),]
B2_44_keep <- B2_44[!(B2_44$Population == 0),]
B2_66_keep <- B2_66[!(B2_66$Population == 0),]
B2_75_keep <- B2_75[!(B2_75$Population == 0),]
B2_90_keep <- B2_90[!(B2_90$Population == 0),]

B3_17_keep <- B3_17[!(B3_17$Population == 0),]
B3_25_keep <- B3_25[!(B3_25$Population == 0),]
B3_44_keep <- B3_44[!(B3_44$Population == 0),]
B3_66_keep <- B3_66[!(B3_66$Population == 0),]
B3_75_keep <- B3_75[!(B3_75$Population == 0),]
B3_90_keep <- B3_90[!(B3_90$Population == 0),]


Genotype_through_time <- as.data.frame(matrix(ncol = 8, byrow=T, c("Population", "Day 0", "Day 17", "Day 25", "Day 44", "Day 66", "Day 75", "Day 90",
                                                              "P1", "0", nrow(P1_17_keep), 
                                                              "0",nrow(P1_44_keep), nrow(P1_66_keep), 
                                                              "0", nrow(P1_90_keep), 
                                                              "P2", "0", nrow(P2_17_keep), 
                                                              "0",nrow(P2_44_keep), nrow(P2_66_keep), 
                                                              "0", nrow(P2_90_keep),
                                                              "P1", "0", nrow(P3_17_keep), 
                                                              "0",nrow(P3_44_keep), nrow(P3_66_keep), 
                                                              "0", nrow(P3_90_keep),
                                                              
                                                              "B1", "0", nrow(B1_17_keep), 
                                                              nrow(B1_25_keep),nrow(B1_44_keep), 
                                                              nrow(B1_66_keep), 
                                                              nrow(B1_75_keep), nrow(B1_90_keep),
                                                              "B2", "0", nrow(B2_17_keep), 
                                                              nrow(B2_25_keep),nrow(B2_44_keep), 
                                                              nrow(B2_66_keep), 
                                                              nrow(B2_75_keep), nrow(B2_90_keep),
                                                              "B3", "0", nrow(B3_17_keep), 
                                                              nrow(B3_25_keep),nrow(B3_44_keep), 
                                                              nrow(B3_66_keep), 
                                                              nrow(B3_75_keep), nrow(B3_90_keep))))
write.csv(Genotype_through_time, "/Users/katrina/Desktop/PALTE_final/draft_revision_work/genotypes_through_time.csv")


gene_parallelism <- as.data.frame(table(all$Gene))
gene_parallelism <- gene_parallelism[(gene_parallelism$Freq >1),]
                                    
write.csv(gene_parallelism, "/Users/katrina/Desktop/PALTE_final/draft_revision_work/gene_parallelism.csv")

#added coordinated but now need to split the column and subtract one column from the next
gene_parallelism <- read.csv("/Users/katrina/Desktop/PALTE_final/draft_revision_work/gene_parallelism.csv", stringsAsFactors = F)

Gene_parallelism <- transform(gene_parallelism, info = colsplit(gene_parallelism$X.1, "::", names = c("start","stop")))
Gene_parallelism$length = Gene_parallelism$info.start-Gene_parallelism$info.stop

write.csv(Gene_parallelism, "/Users/katrina/Desktop/PALTE_final/draft_revision_work/gene_parallelism.csv")
#last import of ths data
fishersexact <- read.csv("/Users/katrina/Desktop/PALTE_final/draft_revision_work/gene_parallelism.csv", stringsAsFactors = F)
fisher <- fishersexact[,c(4,8)]

fisher.test(matrix(c(2, 464, 874, 6537648), nrow = 2))


#grab all the extra data for this data set
fishersexact <- read.csv("/Users/katrina/Desktop/PALTE_final/draft_revision_work/gene_parallelism.csv", stringsAsFactors = F)

all <- read.csv("/Users/katrina/Desktop/PALTE_final/draft_revision_work/all_mutations_200930.csv", stringsAsFactors = F)

parallelism <- all[(all$Gene %in% fishersexact$Gene),]
write.csv(parallelism, "/Users/katrina/Desktop/PALTE_final/draft_revision_work/parallelism_info.csv")





####going to try to find where kenny's B1 clones fit
#import kenny's clones
clones <- read.csv("/Users/katrina/Desktop/PALTE_final/draft_revision_work/KennyClones.csv", stringsAsFactors = F)
B1 <- read.csv("/Users/katrina/Desktop/PALTE_final/draft_revision_work/B1_200927.csv",stringsAsFactors = F)

#only keep the clone calls that are in the all data
B1_clones_noancestor <- clones[!(clones$Position %in% ancestor$SeqID),]
B1_clones_filtered <- B1_clones_noancestor[(B1_clones_noancestor$Position %in% B1$Position),]

unfiltered_clone_table <- table(B1_clones_noancestor$Sample)
filtered_clone_table <- table(B1_clones_filtered$Sample)

write.csv(unfiltered_clone_table, "/Users/katrina/Desktop/PALTE_final/draft_revision_work/clones/B1_unfiltered_clone_table.csv")
write.csv(filtered_clone_table, "/Users/katrina/Desktop/PALTE_final/draft_revision_work/clones/B1_filtered_clone_table.csv")

write.csv(all_together, "/Users/katrina/Desktop/PALTE_final/draft_revision_work/all_201004.csv")

write.csv(B1_clones_filtered, "/Users/katrina/Desktop/PALTE_final/draft_revision_work/B1_clone_filtered.csv")

table(clones$Sample)
table(clone_filtered$Sample)
clone_filtered$frequency <- 100

m_clones <- melt(clone_filtered, id=c("Sample", "Evidence","Position","Mutation","Annotation","Gene","Description", measure.vars="frequency"))
cast_clones <- as.data.frame(t(dcast(m_clones, Position~Sample, mean, value.var = "frequency", fill = 0)))
colnames(cast_clones) <- as.character(unlist(cast_clones[1,]))
clones_table <- cast_clones[-1,]

write.csv(clones_table, "/Users/katrina/Desktop/PALTE_final/draft_revision_work/clone_table.csv")
ncol(clones_table)






#identifying mutations in the named clones from Kenny's study that are in the B1 population data.

B1 <- read.csv("/Users/katrina/Desktop/PALTE_final/draft_revision_work/B1_200927.csv",stringsAsFactors = F)
clones2 <- read.csv("/Users/katrina/Desktop/PALTE_final/draft_revision_work/clones/clones.csv",stringsAsFactors = F)

clones2_noancestor <- clones2[!(clones2$Position %in% ancestor$SeqID),]
clones2_filter <- clones2_noancestor[(clones2_noancestor$Position %in% B1$Position),]

Named_clones_filtered <- table(clones2_filter$Sample)
Named_clones_unfiltered <- table(clones2_noancestor$Sample)

write.csv(clones2_filter, "/Users/katrina/Desktop/PALTE_final/draft_revision_work/clones/clones_B1filter.csv")
write.csv(Named_clones_filtered,"/Users/katrina/Desktop/PALTE_final/draft_revision_work/clones/Named_clones_filtered_numbers.csv")
write.csv(Named_clones_unfiltered,"/Users/katrina/Desktop/PALTE_final/draft_revision_work/clones/Named_clones_unfiltered_numbers.csv")


#identifying mutations in p1 clones that are in the P1 population. 
P1_mutations <- read.csv("/Users/katrina/Desktop/PALTE_final/draft_revision_work/P1_200927.csv",stringsAsFactors = F)
ancestor <- read.csv("/Users/katrina/Desktop/PALTE_final/working/KBH5_WT_Breseq_Output.csv",stringsAsFactors = F)
P1_clones <- read.csv("/Users/katrina/Desktop/PALTE_final/draft_revision_work/P1_clones.csv",stringsAsFactors = F)

p1_clones_noancestor <- P1_clones[!(P1_clones$position %in% ancestor$SeqID),]
P1_clones_filtered <- p1_clones_noancestor[(P1_clones$position %in% P1_mutations$Position),]

unfiltered_p1_clones <- table(p1_clones_noancestor$SampleName)
filtered_p1_clones <- table(P1_clones_filtered$SampleName)

write.csv(unfiltered_p1_clones, "/Users/katrina/Desktop/PALTE_final/draft_revision_work/clones/P1_clones_unfiltered_numbers.csv")
write.csv(filtered_p1_clones, "/Users/katrina/Desktop/PALTE_final/draft_revision_work/clones/P1_clones_filtered_numbers.csv")
write.csv(p1_clones_noancestor, "/Users/katrina/Desktop/PALTE_final/draft_revision_work/clones/P1_clones_filtered.csv")
table(p1_clones_noancestor$gene)

######
#parallelism - gene level
#####
#trying to get a list of all of the mutations with high levels of parallelism only.
parallelism_genes <- c("rluB",
                       "PA14_31070",
                       "PA14_53110",
                       "pslI",
                       "PA14_22650",
                       "PA14_46110",
                       "PA14_18200",
                       "PA14_45060",
                       "PA14_20510",
                       "PA14_32830",
                       "PA14_58070",
                       "PA14_02220",
                       "PA14_58510",
                       "pntB",
                       "nosD",
                       "PA14_69010",
                       "fha1",
                       "rpoB",
                       "dsbA2",
                       "pitA",
                       "argJ",
                       "soxA",
                       "PA14_46030",
                       "cobG",
                       "mutS",
                       "napF",
                       "nqrE",
                       "nuoG",
                       "algF",
                       "PA14_13150",
                       "PA14_51840",
                       "PA14_71750",
                       "PA14_09300",
                       "PA14_47900",
                       "pchF",
                       "PA14_10770",
                       "PA14_32300",
                       "PA14_52260",
                       "PA14_01160",
                       "PA14_54810")

parallelism <- all_together[(all_together$Gene %in% parallelism_genes),]
write.csv(parallelism, "/Users/katrina/Desktop/PALTE_final/draft_revision_work/parallel_genes.csv")

parallelism2 <-all_together[(all_together$Gene == "PA14_23110"),] #for rluB
parallelism2 <-all_together[(all_together$Gene == "PA14_55770"),] #for pitA








#####
#####Now plotting the cohort trajectories

######
B1_cohort_trajectories <- read.csv("/Users/katrina/Desktop/PALTE_final/draft_revision_work/B1.1/B1_200927.genotypes.csv", stringsAsFactors = F)
B2_cohort_trajectories <- read.csv("/Users/katrina/Desktop/PALTE_final/draft_revision_work/B2.2/B2_200927.genotypes.csv", stringsAsFactors = F)
B3_cohort_trajectories <- read.csv("/Users/katrina/Desktop/PALTE_final/draft_revision_work/B3.2/B3_200927.genotypes.csv", stringsAsFactors = F)

P1_cohort_trajectories <- read.csv("/Users/katrina/Desktop/PALTE_final/draft_revision_work/P1.1/P1_200927.genotypes.csv", stringsAsFactors = F)
P2_cohort_trajectories <- read.csv("/Users/katrina/Desktop/PALTE_final/draft_revision_work/P2.1/P2_200927.genotypes.csv", stringsAsFactors = F)
P3_cohort_trajectories <- read.csv("/Users/katrina/Desktop/PALTE_final/draft_revision_work/P3.1/P3_200927.genotypes.csv", stringsAsFactors = F)

B1_cohort_trajectories <- B1_cohort_trajectories[,2:8]
B2_cohort_trajectories <- B2_cohort_trajectories[,2:8]
B3_cohort_trajectories <- B3_cohort_trajectories[,2:8]

P1_cohort_trajectories <- P1_cohort_trajectories[,2:6]
P2_cohort_trajectories <- P2_cohort_trajectories[,2:6]
P3_cohort_trajectories <- P3_cohort_trajectories[,2:6]

b_generations <- c(0,113,167,293,440,500,600)
p_generations <- c(0,113,293,500,600)


#B1_cohorts_edited <- B1_cohorts_edited*100
layout(matrix(c(1,2,3,4,5,6),2))

plot(NA, xlim=c(0,600), ylim=c(0,1), ylab = "Frequency", xlab = "Generations", main = "B1", cex.axis=2, cex.lab=2, cex.main=3) #, main = "Mutation frequencies over time") # 
for(i in seq_len(nrow(B1_cohort_trajectories))){
  lines(b_generations,B1_cohort_trajectories[i,], type="l", col=1)
}
lines(b_generations, B1_cohort_trajectories[2,], type="l", col = "red", lwd=3)


plot(NA, xlim=c(0,600), ylim=c(0,1), ylab = "", xlab = "", main = "P1", cex.axis=2, cex.main=3) #, main = "Mutation frequencies over time") # 
for(i in seq_len(nrow(P1_cohort_trajectories))){
  lines(p_generations,P1_cohort_trajectories[i,], type="l", col=1)
} #need to get them colored differently, but will work for now. 

plot(NA, xlim=c(0,600), ylim=c(0,1), ylab = "", xlab = "", main = "B2", cex.axis=2, cex.main=3) #, main = "Mutation frequencies over time") # 
for(i in seq_len(nrow(B2_cohort_trajectories))){
  lines(b_generations,B2_cohort_trajectories[i,], type="l", col=1)
} #need to get them colored differently, but will work for now. 
lines(b_generations,B2_cohort_trajectories[3,], type="l", col="red", lwd=3)

plot(NA, xlim=c(0,600), ylim=c(0,1), ylab = "", xlab = "", main = "P2", cex.axis=2, cex.main=3) #, main = "Mutation frequencies over time") # 
for(i in seq_len(nrow(P2_cohort_trajectories))){
  lines(p_generations,P2_cohort_trajectories[i,], type="l", col=1)
} #need to get them colored differently, but will work for now. 

plot(NA, xlim=c(0,600), ylim=c(0,1), ylab = "", xlab = "", main = "B3", cex.axis=2, cex.main=3) #, main = "Mutation frequencies over time") # 
for(i in seq_len(nrow(B3_cohort_trajectories))){
  lines(b_generations,B3_cohort_trajectories[i,], type="l", col=1)
} #need to get them colored differently, but will work for now. 
lines(b_generations, B3_cohort_trajectories[3,], type="l", col="red", lwd=3)

plot(NA, xlim=c(0,600), ylim=c(0,1), ylab = "", xlab = "", main = "P3", cex.axis=2, cex.main=3) #, main = "Mutation frequencies over time") # 
for(i in seq_len(nrow(P3_cohort_trajectories))){
  lines(p_generations,P3_cohort_trajectories[i,], type="l", col=1)
} #need to get them colored differently, but will work for now. 
lines(p_generations, P3_cohort_trajectories[4,], type="l", col="red", lwd=3)





######
#looking at the extent of within population gene-level parallelism
B1 <- read.csv("/Users/katrina/Desktop/PALTE_final/draft_revision_work/B1_200927.csv",stringsAsFactors = F)
B2 <- read.csv("/Users/katrina/Desktop/PALTE_final/draft_revision_work/B2_200927.csv",stringsAsFactors = F)
B3 <- read.csv("/Users/katrina/Desktop/PALTE_final/draft_revision_work/B3_200927.csv",stringsAsFactors = F)
P1 <- read.csv("/Users/katrina/Desktop/PALTE_final/draft_revision_work/P1_200927.csv",stringsAsFactors = F)
P2 <- read.csv("/Users/katrina/Desktop/PALTE_final/draft_revision_work/P2_200927.csv",stringsAsFactors = F)
P3 <- read.csv("/Users/katrina/Desktop/PALTE_final/draft_revision_work/P3_200927.csv",stringsAsFactors = F)

B1_within <- table(B1$Gene)
B2_within <- table(B2$Gene)
B3_within <- table(B3$Gene)
P1_within <- table(P1$Gene)
P2_within <- table(P2$Gene)
P3_within <- table(P3$Gene)

write.csv(B1_within, "/Users/katrina/Desktop/PALTE_final/draft_revision_work/B1_gene_level_parallelism.csv")
write.csv(B2_within, "/Users/katrina/Desktop/PALTE_final/draft_revision_work/B2_gene_level_parallelism.csv")
write.csv(B3_within, "/Users/katrina/Desktop/PALTE_final/draft_revision_work/B3_gene_level_parallelism.csv")
write.csv(P1_within, "/Users/katrina/Desktop/PALTE_final/draft_revision_work/P1_gene_level_parallelism.csv")
write.csv(P2_within, "/Users/katrina/Desktop/PALTE_final/draft_revision_work/P2_gene_level_parallelism.csv")
write.csv(P3_within, "/Users/katrina/Desktop/PALTE_final/draft_revision_work/P3_gene_level_parallelism.csv")

#I externally combined all of the genes with population identifiers that had 2 or more mutations. 
within_pop_parallelism <- read.csv("/Users/katrina/Desktop/PALTE_final/draft_revision_work/withinpop_gene_level_parallelism.csv", stringsAsFactors = F)
all_mutations <- read.csv("/Users/katrina/Desktop/PALTE_final/draft_revision_work/all_mutations_final_201031.csv", stringsAsFactors = F)

#now make a document where only the parallelism cases are saved
parallelism_within_pop <- all_mutations[(all_mutations$Gene %in% within_pop_parallelism$Gene),]
write.csv(parallelism_within_pop, "/Users/katrina/Desktop/PALTE_final/draft_revision_work/parallelism_within_pop.csv")










#####
#####Now plotting the cohort trajectories fro just genotypes with within population parallelism

######
B1_cohort_trajectories <- read.csv("/Users/katrina/Desktop/PALTE_final/draft_revision_work/withinpop_genotypes/B1_200927.genotypes.csv", stringsAsFactors = F)
B2_cohort_trajectories <- read.csv("/Users/katrina/Desktop/PALTE_final/draft_revision_work/withinpop_genotypes/B2_200927.genotypes.csv", stringsAsFactors = F)
B3_cohort_trajectories <- read.csv("/Users/katrina/Desktop/PALTE_final/draft_revision_work/withinpop_genotypes/B3_200927.genotypes.csv", stringsAsFactors = F)

P1_cohort_trajectories <- read.csv("/Users/katrina/Desktop/PALTE_final/draft_revision_work/withinpop_genotypes/P1_200927.genotypes.csv", stringsAsFactors = F)
P2_cohort_trajectories <- read.csv("/Users/katrina/Desktop/PALTE_final/draft_revision_work/withinpop_genotypes/P2_200927.genotypes.csv", stringsAsFactors = F)
P3_cohort_trajectories <- read.csv("/Users/katrina/Desktop/PALTE_final/draft_revision_work/withinpop_genotypes/P3_200927.genotypes.csv", stringsAsFactors = F)

B1_cohort_trajectories <- B1_cohort_trajectories[,2:8]
B2_cohort_trajectories <- B2_cohort_trajectories[,2:8]
B3_cohort_trajectories <- B3_cohort_trajectories[,2:8]

P1_cohort_trajectories <- P1_cohort_trajectories[,2:6]
P2_cohort_trajectories <- P2_cohort_trajectories[,2:6]
P3_cohort_trajectories <- P3_cohort_trajectories[,2:6]

b_generations <- c(0,113,167,293,440,500,600)
p_generations <- c(0,113,293,500,600)


#B1_cohorts_edited <- B1_cohorts_edited*100
layout(matrix(c(1,2,3,4,5,6),2))

plot(NA, xlim=c(0,600), ylim=c(0,1), ylab = "Frequency", xlab = "Generations", main = "B1", cex.axis=2, cex.lab=2, cex.main=3) #, main = "Mutation frequencies over time") # 
for(i in seq_len(nrow(B1_cohort_trajectories))){
  lines(b_generations,B1_cohort_trajectories[i,], type="l", col=1)
}


plot(NA, xlim=c(0,600), ylim=c(0,1), ylab = "", xlab = "", main = "P1", cex.axis=2, cex.main=3) #, main = "Mutation frequencies over time") # 
for(i in seq_len(nrow(P1_cohort_trajectories))){
  lines(p_generations,P1_cohort_trajectories[i,], type="l", col=1)
} #need to get them colored differently, but will work for now. 

plot(NA, xlim=c(0,600), ylim=c(0,1), ylab = "", xlab = "", main = "B2", cex.axis=2, cex.main=3) #, main = "Mutation frequencies over time") # 
for(i in seq_len(nrow(B2_cohort_trajectories))){
  lines(b_generations,B2_cohort_trajectories[i,], type="l", col=1)
} #need to get them colored differently, but will work for now. 

plot(NA, xlim=c(0,600), ylim=c(0,1), ylab = "", xlab = "", main = "P2", cex.axis=2, cex.main=3) #, main = "Mutation frequencies over time") # 
for(i in seq_len(nrow(P2_cohort_trajectories))){
  lines(p_generations,P2_cohort_trajectories[i,], type="l", col=1)
} #need to get them colored differently, but will work for now. 

plot(NA, xlim=c(0,600), ylim=c(0,1), ylab = "", xlab = "", main = "B3", cex.axis=2, cex.main=3) #, main = "Mutation frequencies over time") # 
for(i in seq_len(nrow(B3_cohort_trajectories))){
  lines(b_generations,B3_cohort_trajectories[i,], type="l", col=1)
} #need to get them colored differently, but will work for now. 

plot(NA, xlim=c(0,600), ylim=c(0,1), ylab = "", xlab = "", main = "P3", cex.axis=2, cex.main=3) #, main = "Mutation frequencies over time") # 
for(i in seq_len(nrow(P3_cohort_trajectories))){
  lines(p_generations,P3_cohort_trajectories[i,], type="l", col=1)
} #need to get them colored differently, but will work for now. 






#creating a collectors curve: 
#import the metagenomic data
metagenomes <- read.csv("/Users/katrina/Desktop/PALTE_final/draft_revision_work/all_mutations_final_201031.csv", stringsAsFactors = F)


metagenomes_day90 <- metagenomes[!(metagenomes$X90 == 0),]
unique_metagenome <- table(metagenomes_day90$Position)#520
nrow(metagenomes_day90)#583

#Now split into populations
B1 <- metagenomes[(metagenomes$Population == "B1"),]
B2 <- metagenomes[(metagenomes$Population == "B2"),]
B3 <- metagenomes[(metagenomes$Population == "B3"),]
P1 <- metagenomes[(metagenomes$Population == "P1"),]
P2 <- metagenomes[(metagenomes$Population == "P2"),]
P3 <- metagenomes[(metagenomes$Population == "P3"),]

#need number of populations and number of unique positions mutated.

#for 1 population - P1
one <- table(P1$Position) #109

#for 2 populations - P1 and P2
two_data <- rbind(P1, P2)
two <- table(two_data$Position)#186

#for 3 - all planktonic
three_data <- rbind(two_data, P3)
three <- table(three_data$Position) #298

#for 4 - plank and B1
four_data <- rbind(three_data, B1)
four <- table(four_data$Position) #512

#for 5 - plank, B1, and B2
five_data <- rbind(four_data, B2)
five <- table(five_data$Position) #648

#for all 6
six_data <- rbind(five_data, B3)
six <- table(six_data$Position) #724


#then do the same for clones
B1_clones <- read.csv("/Users/katrina/Desktop/PALTE_final/draft_revision_work/B1_clone_filtered.csv", stringsAsFactors = F)

named_clones <- read.csv("/Users/katrina/Desktop/PALTE_final/draft_revision_work/clones/clones_B1filter.csv", stringsAsFactors = F)

P1_clones <- read.csv("/Users/katrina/Desktop/PALTE_final/draft_revision_work/clones/P1_clones_filtered.csv", stringsAsFactors = F)

names <- table(P1_clones$Sample)
PA101 <- B1_clones[(B1_clones$Sample =="PA101"),]
PA102 <- B1_clones[(B1_clones$Sample =="PA102"),]
PA103 <- B1_clones[(B1_clones$Sample =="PA103"),]
PA104 <- B1_clones[(B1_clones$Sample =="PA104"),]
PA105 <- B1_clones[(B1_clones$Sample =="PA105"),]
PA106 <- B1_clones[(B1_clones$Sample =="PA106"),]
PA107 <- B1_clones[(B1_clones$Sample =="PA107"),]
PA108 <- B1_clones[(B1_clones$Sample =="PA108"),]
PA109 <- B1_clones[(B1_clones$Sample =="PA109"),]
PA110 <- B1_clones[(B1_clones$Sample =="PA110"),]
PA111 <- B1_clones[(B1_clones$Sample =="PA111"),]
PA112 <- B1_clones[(B1_clones$Sample =="PA112"),]
PA113 <- B1_clones[(B1_clones$Sample =="PA113"),]
PA114 <- B1_clones[(B1_clones$Sample =="PA114"),]
PA116 <- B1_clones[(B1_clones$Sample =="PA116"),]
PA117 <- B1_clones[(B1_clones$Sample =="PA117"),]
PA118 <- B1_clones[(B1_clones$Sample =="PA118"),]
PA119 <- B1_clones[(B1_clones$Sample =="PA119"),]
PA121 <- B1_clones[(B1_clones$Sample =="PA121"),]

PA20 <- named_clones[(named_clones$Sample =="20"),]
PA21 <- named_clones[(named_clones$Sample =="21"),]
PA22 <- named_clones[(named_clones$Sample =="22"),]
PA23 <- named_clones[(named_clones$Sample =="23"),]
PA24 <- named_clones[(named_clones$Sample =="24"),]
PA25 <- named_clones[(named_clones$Sample =="25"),]
PA26 <- named_clones[(named_clones$Sample =="26"),]


KBH110 <- P1_clones[(P1_clones$SampleName == "KBH_110"),]
KBH111 <- P1_clones[(P1_clones$SampleName == "KBH_111"),]
KBH112 <- P1_clones[(P1_clones$SampleName == "KBH_112"),]
KBH113 <- P1_clones[(P1_clones$SampleName == "KBH_113"),]
KBH114 <- P1_clones[(P1_clones$SampleName == "KBH_114"),]
KBH115 <- P1_clones[(P1_clones$SampleName == "KBH_115"),]
KBH116 <- P1_clones[(P1_clones$SampleName == "KBH_116"),]
KBH117 <- P1_clones[(P1_clones$SampleName == "KBH_117"),]
KBH118 <- P1_clones[(P1_clones$SampleName == "KBH_118"),]
KBH119 <- P1_clones[(P1_clones$SampleName == "KBH_119"),]


#now find the number of unique positons
#starting with planktonic
one <- table(KBH110$position)#14

two_data <- rbind(KBH110, KBH111)
two <- table(two_data$position)#28

three_data <- rbind(two_data, KBH112)
three <- table(three_data$position)#29

four_data <- rbind(three_data, KBH113)
four <- table(four_data$position)#29

five_data <- rbind(four_data, KBH114)
five <- table(five_data$position)#33

six_data <- rbind(five_data, KBH115)
six <- table(six_data$position)#35

seven_data <- rbind(six_data, KBH116)
seven <- table(seven_data$position)#36

eight_data <- rbind(seven_data, KBH117)
eight <- table(eight_data$position)#36

nine_data <- rbind(eight_data, KBH118)
nine <- table(nine_data$position)#37

ten_data <- rbind(nine_data, KBH119)
ten <- table(ten_data$position)#45
ten_data <- ten_data[,1:8] #need to change so that the columns are the same between this data set and the next data set
colnames(ten_data) <- colnames(PA20) #and changing the column names so that I can rbind

eleven_data <- rbind(ten_data, PA20)
eleven <- table(eleven_data$Position)#115

twelve_data <- rbind(eleven_data, PA21)
twelve <- table(twelve_data$Position)#115

thirteen_data <- rbind(twelve_data, PA22)
thirteen<- table(thirteen_data$Position)#128

fourteen_data <- rbind(thirteen_data, PA23)
fourteen<- table(fourteen_data$Position)#128

fifteen_data <- rbind(fourteen_data, PA24)
fifteen<- table(fifteen_data$Position)#128

sixteen_data <- rbind(fifteen_data, PA25)
sixteen<- table(sixteen_data$Position)#138

seventeen_data <- rbind(sixteen_data, PA26)
seventeen<- table(seventeen_data$Position)#138

seventeen_data <- rbind(sixteen_data, PA26)
seventeen<- table(seventeen_data$Position)#138



eightteen_data <- rbind(seventeen_data, PA101)
eightteen<- table(eightteen_data$Position)#157

nineteen_data <- rbind(eightteen_data, PA102)
nineteen<- table(nineteen_data$Position)#157

twenty_data <- rbind(nineteen_data, PA103)
twenty <- table(twenty_data$Position)#157

twentyone_data <- rbind(twenty_data, PA104)
twentyone <- table(twentyone_data$Position)#157

twentytwo_data <- rbind(twentyone_data, PA105)
twentytwo <- table(twentytwo_data$Position)#157

twentytwo_data <- rbind(twentyone_data, PA106)
twentytwo <- table(twentytwo_data$Position)#157



P1_clones_2 <- P1_clones[,1:8]
colnames(P1_clones_2) <- colnames(B1_clones)

all_clones <- rbind(P1_clones_2, B1_clones, named_clones)
all <- table(all_clones$Position)  #161 total number of unique sites mutated
nrow(all_clones) #1377 total number of mutations
all_table <- table(all_clones$Sample)






#determining the number of mutations and the number of fixation events at each time point
metagenomes <- read.csv("/Users/katrina/Desktop/PALTE_final/draft_revision_work/all_mutations_final_201031.csv", stringsAsFactors = F)

B1 <- metagenomes[(metagenomes$Population == "B1"),]
B1_17 <- B1[(B1$X17 > 0),]
B1_17_fixed <- B1_17[(B1_17$X17 == 1),]
B1_25 <- B1[(B1$X25 > 0),]
B1_25_fixed <- B1_25[(B1_25$X25 == 1),]
B1_44 <- B1[(B1$X44 > 0),]
B1_44_fixed <- B1_44[(B1_44$X44 == 1),] #11
B1_66 <- B1[(B1$X66 > 0),]
B1_66_fixed <- B1_66[(B1_66$X66 == 1),]
B1_75 <- B1[(B1$X75 > 0),]
B1_75_fixed <- B1_75[(B1_75$X75 == 1),] #19
B1_90 <- B1[(B1$X90 > 0),]
B1_90_fixed <- B1_90[(B1_90$X90 == 1),]

B2 <- metagenomes[(metagenomes$Population == "B2"),]
B2_17 <- B2[(B2$X17 > 0),]
B2_17_fixed <- B2_17[(B2_17$X17 == 1),]
B2_25 <- B2[(B2$X25 > 0),]
B2_25_fixed <- B2_25[(B2_25$X25 == 1),]
B2_44 <- B2[(B2$X44 > 0),]
B2_44_fixed <- B2_44[(B2_44$X44 == 1),]
B2_66 <- B2[(B2$X66 > 0),]
B2_66_fixed <- B2_66[(B2_66$X66 == 1),]
B2_75 <- B2[(B2$X75 > 0),]
B2_75_fixed <- B2_75[(B2_75$X75 == 1),]
B2_90 <- B2[(B2$X90 > 0),]
B2_90_fixed <- B2_90[(B2_90$X90 == 1),]


B3 <- metagenomes[(metagenomes$Population == "B3"),]
B3_17 <- B3[(B3$X17 > 0),]
B3_17_fixed <- B3_17[(B3_17$X17 == 1),]
B3_25 <- B3[(B3$X25 > 0),]
B3_25_fixed <- B3_25[(B3_25$X25 == 1),]
B3_44 <- B3[(B3$X44 > 0),]
B3_44_fixed <- B3_44[(B3_44$X44 == 1),]
B3_66 <- B3[(B3$X66 > 0),]
B3_66_fixed <- B3_66[(B3_66$X66 == 1),]
B3_75 <- B3[(B3$X75 > 0),]
B3_75_fixed <- B3_75[(B3_75$X75 == 1),]
B3_90 <- B3[(B3$X90 > 0),]
B3_90_fixed <- B3_90[(B3_90$X90 == 1),]


P1 <- metagenomes[(metagenomes$Population == "P1"),]
P1_17 <- P1[(P1$X17 > 0),]
P1_17_fixed <- P1_17[(P1_17$X17 == 1),]
P1_44 <- P1[(P1$X44 > 0),]
P1_44_fixed <- P1_44[(P1_44$X44 == 1),]
P1_66 <- P1[(P1$X66 > 0),]
P1_66_fixed <- P1_66[(P1_66$X66 == 1),]
P1_90 <- P1[(P1$X90 > 0),]
P1_90_fixed <- P1_90[(P1_90$X90 == 1),]


P2 <- metagenomes[(metagenomes$Population == "P2"),]
P2_17 <- P2[(P2$X17 > 0),]
P2_17_fixed <- P2_17[(P2_17$X17 == 1),]
P2_44 <- P2[(P2$X44 > 0),]
P2_44_fixed <- P2_44[(P2_44$X44 == 1),]
P2_66 <- P2[(P2$X66 > 0),]
P2_66_fixed <- P2_66[(P2_66$X66 == 1),]
P2_90 <- P2[(P2$X90 > 0),]
P2_90_fixed <- P2_90[(P2_90$X90 == 1),]


P3 <- metagenomes[(metagenomes$Population == "P3"),]
P3_17 <- P3[(P3$X17 > 0),]
P3_17_fixed <- P3_17[(P3_17$X17 == 1),]
P3_44 <- P3[(P3$X44 > 0),]
P3_44_fixed <- P3_44[(P3_44$X44 == 1),]
P3_66 <- P3[(P3$X66 > 0),]
P3_66_fixed <- P3_66[(P3_66$X66 == 1),]
P3_90 <- P3[(P3$X90 > 0),]
P3_90_fixed <- P3_90[(P3_90$X90 == 1),]

