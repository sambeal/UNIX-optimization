#Analysis of ASVs and depths from different UNIX parameters
#1: all default
#2: cutadapt -e0.15 
#3: cutadapt -e0.5
#4: cutadapt -e0.5 & fastx_filter 2
#5: cutadapt -e0.4 & fastx_filter 1
#6: cutadapt -e0.4 & fastx_filter 2
#7: cutadapt -e0.4, fastx_filter 2, unoise_alpha 1, minsize 8
#8: cutadapt -e0.4, fastx_filter 2, unoise_alpha 2, minsize 4
#9: cutadapt -e0.4, fastx_filter 2, unoise_alpha 1, minsize 4
#10: cutadapt -e0.4, fastx_filter 2, unoise_alpha 3, minsize 4

#set wd to "/Users/samanthabeal/Documents/MSc/Bioinformatics"
getwd()
#"/Users/samanthabeal"
setwd("Documents/MSc/Bioinformatics")

#load data
o1 <- read.csv("UNIX-optimization/organized_ASVs/output1_ASVanalysis.csv", header = TRUE)
o2 <- read.csv("UNIX-optimization/organized_ASVs/output2_ASVanalysis.csv", header = TRUE)
o3 <- read.csv("UNIX-optimization/organized_ASVs/output3_ASVanalysis.csv", header = TRUE)
o4 <- read.csv("UNIX-optimization/organized_ASVs/output4_ASVanalysis.csv", header = TRUE)
o5 <- read.csv("UNIX-optimization/organized_ASVs/output5_ASVanalysis.csv", header = TRUE)
o6 <- read.csv("UNIX-optimization/organized_ASVs/output6_ASVanalysis.csv", header = TRUE)
o7 <- read.csv("UNIX-optimization/organized_ASVs/output7_ASVanalysis.csv", header = TRUE)
o8 <- read.csv("UNIX-optimization/organized_ASVs/output8_ASVanalysis.csv", header = TRUE)
o9 <- read.csv("UNIX-optimization/organized_ASVs/output9_ASVanalysis.csv", header = TRUE)
o10 <- read.csv("UNIX-optimization/organized_ASVs/output10_ASVanalysis.csv", header = TRUE)

output <- rbind(o1,o2,o3,o4,o5,o6,o7,o8,o9,o10)

#check data
summary(output)

#change data types and col names
names(output)[8] <- "Output"
output$Fish <- as.factor(output$Fish)
output$Output <- as.factor(output$Output)
output$All.reps <- as.factor(output$All.reps)

#check data
summary(output)

#plot
ggplot(output, aes(fill=All.reps, y=Avg.depth, x=output)) + 
  geom_bar(position="dodge", stat="identity") + 
  theme_linedraw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("Parameter test") + ylab("Average read depth") + 
  labs(fill= "ASV detected in all PCR replicates?")


  geom_text(aes(label = round(avgASV$Num.ASVs, digits = 0)), position=position_dodge(width=0.9), vjust=-0.25, size=3)


#want to add a label to each bar showing how many ASVs this average was taken from
avgASV <- aggregate(Num.ASVs ~ output + All.reps, data = all_cut2, FUN = mean)

avgdepth <- aggregate(Avg.depth ~ Fish + All.reps, data = all_cut2, FUN = mean)


#plot with bars for individual fish





#true ASV comparison----
o1_true <- read.csv("UNIX-optimization/organized_ASVs/output1_TrueASVs.csv", header = TRUE)
o2_true <- read.csv("UNIX-optimization/organized_ASVs/output2_TrueASVs.csv", header = TRUE)
o3_true <- read.csv("UNIX-optimization/organized_ASVs/output3_TrueASVs.csv", header = TRUE)
o4_true <- read.csv("UNIX-optimization/organized_ASVs/output4_TrueASVs.csv", header = TRUE)
o5_true <- read.csv("UNIX-optimization/organized_ASVs/output5_TrueASVs.csv", header = TRUE)
o6_true <- read.csv("UNIX-optimization/organized_ASVs/output6_TrueASVs.csv", header = TRUE)
o7_true <- read.csv("UNIX-optimization/organized_ASVs/output7_TrueASVs.csv", header = TRUE)
o8_true <- read.csv("UNIX-optimization/organized_ASVs/output8_TrueASVs.csv", header = TRUE)
o9_true <- read.csv("UNIX-optimization/organized_ASVs/output9_TrueASVs.csv", header = TRUE)
o10_true <- read.csv("UNIX-optimization/organized_ASVs/output10_TrueASVs.csv", header = TRUE)

output_true <- rbind(o1_true,o2_true,o3_true,o4_true,o5_true,o6_true,o7_true,o8_true,o9_true,o10_true)


#Add seqs to table
#UNIX output does not have the actual sequence, only IDs (e.g. uniq1)
#need to get them from earlier in the pipeline and add them to the table
#converted denoised.fasta to .tab, imported .tab and now will reformat
#http://sequenceconversion.bugaco.com/converter/biology/sequences/fasta_to_tab.php 
#using output9 file as it had the most ASVs (see: UNIX-opt true ASV analysis.R file to see how I came to this conclusion)

#read in data
UNIXseqsS <- read.table("UNIX-optimization/organized_ASVs/SINEo9seqs.tab")

#reformat
library(tibble)
UNIXseqsS <- tibble::rownames_to_column(UNIXseqsS, "ASV")
names(UNIXseqsS)[2] <- "Sequence"

#",.*" says: replace comma (,) and every character after that (.*), with nothing "".
UNIXseqsS$ASV <- gsub(";.*", "", UNIXseqsS$ASV)


library(tidyverse)
UNIXtop_asvSeqs <- output_true %>% 
  inner_join(UNIXseqsS, by = c("ASV" = "ASV"))

UNIXtop_asvSeqs <- UNIXtop_asvSeqs %>% relocate(Sequence, .after = ASV)

#plot consistently detected ASVs
UNIXtop_asvSeqs$Sequence_length <- str_length(UNIXtop_asvSeqs$Sequence)

hist(UNIXtop_asvSeqs$Sequence_length, main="Distribution of consitently detected SINE ASVs from gDNA - UNIX outputs",
     xlab="Length of ASV seqeuence")

#align the seqs----
#make data frame of only ASV and seq then convert to .fasta
UNIXfasta <- UNIXtop_asvSeqs[c(1:2)]
names(UNIXfasta)[1] <- "ASV"
UNIXfasta$ASV <- paste0(">", UNIXfasta$ASV)
library(tidyverse)
UNIXall_print <-
  UNIXfasta %>%
  select(ASV,Sequence) %>%
  rowwise() %>%
  pivot_longer(ASV:Sequence) %>%
  select(-name)

#allseqs
write.table(UNIXall_print, 
            file = "R/UNIX v dada2/UNIXSINE2all_ASV.fasta", 
            col.names = FALSE, 
            row.names = FALSE, 
            quote = FALSE)

#this (^) made the fasta seq of all ASVs (even duplicates)
# want to align a file of only each unique seq
