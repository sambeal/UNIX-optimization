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
#11: (trying to match dada2 parameters) 
    # cutadapt -e0.1,fastq_maxee 2 --fastq_maxns 0 --fastq_truncqual 2 --fastq_minlen 50, minsize 4, alpha 2 
    #(denosing is very different b/w the two so unsure if these were the right parameters to use or not)
#12: cutadapt -e0.4, fastq_filter maxEE 2, minsize 4, alpha 5
#13: cutadapt -e0.4, fastq_filter maxEE 2, minsize 4, alpha 10


#last pushed to github: after output 13 23.03.21

#set wd to "/Users/samanthabeal/Documents/MSc/Bioinformatics"
getwd()
#"/Users/samanthabeal"
setwd("Documents/MSc/Bioinformatics")

#load data for all ASVs----
o01 <- read.csv("UNIX-optimization/organized_ASVs/output1_ASVanalysis.csv", header = TRUE)
o02 <- read.csv("UNIX-optimization/organized_ASVs/output2_ASVanalysis.csv", header = TRUE)
o03 <- read.csv("UNIX-optimization/organized_ASVs/output3_ASVanalysis.csv", header = TRUE)
o04 <- read.csv("UNIX-optimization/organized_ASVs/output4_ASVanalysis.csv", header = TRUE)
o05 <- read.csv("UNIX-optimization/organized_ASVs/output5_ASVanalysis.csv", header = TRUE)
o06 <- read.csv("UNIX-optimization/organized_ASVs/output6_ASVanalysis.csv", header = TRUE)
o07 <- read.csv("UNIX-optimization/organized_ASVs/output7_ASVanalysis.csv", header = TRUE)
o08 <- read.csv("UNIX-optimization/organized_ASVs/output8_ASVanalysis.csv", header = TRUE)
o09 <- read.csv("UNIX-optimization/organized_ASVs/output9_ASVanalysis.csv", header = TRUE)
o10 <- read.csv("UNIX-optimization/organized_ASVs/output10_ASVanalysis.csv", header = TRUE)
o11 <- read.csv("UNIX-optimization/organized_ASVs/output11_ASVanalysis.csv", header = TRUE)
o12 <- read.csv("UNIX-optimization/organized_ASVs/output12_ASVanalysis.csv", header = TRUE)
o13 <- read.csv("UNIX-optimization/organized_ASVs/output13_ASVanalysis.csv", header = TRUE)


output <- rbind(o01,o02,o03,o04,o05,o06,o07,o08,o09,o10,o11,o12,o13)

#check data
summary(output)

#change data types and col names
output$Fish <- as.factor(output$Fish)
output$Output <- as.factor(output$Output)
output$All.reps <- as.factor(output$All.reps)
output$Total.depth <- as.numeric(output$Total.depth)


#check data
summary(output)

#reorder factors (should have named outputs 01-10 instead of 1-10)
output$Output <- factor(output$Output, levels=c('output1','output2',
                                                 'output3','output4',
                                                 'output5','output6',
                                                 'output7','output8',
                                                 'output9','output10',
                                                 'output11','output12',
                                                 'output13'))


#plot total read depth across all fish 
#TOTAL READ DEPTH IS MORE INFORMATIVE AT THIS STAGE THAN AVG. DEPTH/FISH

#need to sum the total depth first, otherwise the bars will only be the highest depth/output
totaldepth <- aggregate(Total.depth ~ Output + All.reps, data = output, FUN = sum)
aggregate(Num.ASVs ~ Output + All.reps, data = output, FUN = sum)

#add labels for the number of ASVs this data came from:
#totalASVnum <- aggregate(Num.ASVs ~ Output + All.reps, data = output, FUN = sum)
#this is not helpful right now since it is all ASVs, not just unique ones
# geom_text(aes(label = round(totalASVnum$Num.ASVs, digits = 0)), position=position_dodge(width=0.9), vjust=-0.25, size=3)


#plot
ggplot(totaldepth, aes(fill=All.reps, y=Total.depth, x=Output)) + 
  geom_bar(position="dodge", stat="identity") + 
  theme_linedraw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "bottom") +
  xlab("Parameter test") + ylab("Total read depth (across all 16 fish)") + 
  ggtitle("Total depth retained across all technical replicates using different UNIX pipeline parameters") +
  labs(fill= "ASV detected in all PCR replicates?")




#look at the averages----
#avg. number of ASVs per output
avgASVnum <- aggregate(Num.ASVs ~ Output + All.reps, data = output, FUN = mean)

#avg. number of ASVs per fish per output
avgASVnumFish <- aggregate(Num.ASVs ~ Output + All.reps + Fish, data = output, FUN = mean)

#avg. depth of ASVs per output
avgdepth <- aggregate(Avg.depth ~ Output + All.reps, data = output, FUN = mean)

#avg. depth of ASVs per fish per output
avgdepthFish <- aggregate(Avg.depth ~ Output + All.reps +Fish, data = output, FUN = mean)

#plot the averages
ggplot(avgdepth, aes(fill=All.reps, y=Avg.depth, x=Output)) + 
  geom_bar(position="dodge", stat="identity") + 
  theme_linedraw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "bottom") +
  xlab("Parameter test") + ylab("Average read depth (per fish)") + 
  ggtitle("Avg. depth retained across all technical replicates using different UNIX pipeline parameters") +
labs(fill= "ASV detected in all PCR replicates?") +
  geom_text(aes(label = round(avgASVnum$Num.ASVs, digits = 0)), position=position_dodge(width=0.9), vjust=-0.25, size=3)

#text label = average number of ASVs PER FISH in the output


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
o11_true <- read.csv("UNIX-optimization/organized_ASVs/output11_TrueASVs.csv", header = TRUE)
o12_true <- read.csv("UNIX-optimization/organized_ASVs/output12_TrueASVs.csv", header = TRUE)
o13_true <- read.csv("UNIX-optimization/organized_ASVs/output13_TrueASVs.csv", header = TRUE)


output_true <- rbind(o1_true,o2_true,o3_true,o4_true,o5_true,o6_true,o7_true,o8_true,o9_true,o10_true,o11_true,o12_true,o13_true)
length(unique(output_true$ASV))
#output-1-13: 10,178 unique ASVs across all outputs

#Add seqs to table
#UNIX output does not have the actual sequence, only IDs (e.g. uniq1)
#need to get them from earlier in the pipeline and add them to the table
#converted denoised.fasta to .tab, imported .tab and now will reformat
#http://sequenceconversion.bugaco.com/converter/biology/sequences/fasta_to_tab.php 
#using output13 file as it had the most ASVs (see: UNIX-opt true ASV analysis.R file to see how I came to this conclusion)

#read in seqs in R
library (devtools)
library (tidyverse)
source_url("https://raw.githubusercontent.com/lrjoshi/FastaTabular/master/fasta_and_tabular.R")

FastaToTabular("UNIX-optimization/output13/ASV/denoised.fasta")
#The output will be stored as dna_table.csv in the current directory.
#this put me into the files enviornment
#f then q to leave, but needed to reload all previous data so need to find a better way

#read in data
UNIXseqs <- read.csv("dna_table.csv")

#reformat
UNIXseqs <- UNIXseqs[-1]
names(UNIXseqs)[1] <- "ASV"
names(UNIXseqs)[2] <- "Sequence"

#";.*" says: replace semicolon (;) and every character after that (.*), with nothing "".
UNIXseqs$ASV <- gsub(";.*", "", UNIXseqs$ASV)
#">" says: replace arrow (>) with nothing "".
UNIXseqs$ASV <- gsub(">", "", UNIXseqs$ASV)

#add data frames together
library(tidyverse)
output_true_seqs <- output_true %>% 
  inner_join(UNIXseqs, by = c("ASV" = "ASV"))

output_true_seqs <- output_true_seqs %>% relocate(Sequence, .after = ASV)

#add seq lengths
output_true_seqs$Sequence_length <- str_length(output_true_seqs$Sequence)

#check
summary(output_true_seqs)

#change data types
output_true_seqs$ASV <- as.factor(output_true_seqs$ASV)
output_true_seqs$Fish <- as.factor(output_true_seqs$Fish)
output_true_seqs$Output <- as.factor(output_true_seqs$Output)

#reorder factors (should have named outputs 01-10 instead of 1-10)
output_true_seqs$Output <- factor(output_true_seqs$Output, levels=c('output1','output2',
                                                                    'output3','output4',
                                                                    'output5','output6',
                                                                    'output7','output8',
                                                                    'output9','output10',
                                                                    'output11','output12',
                                                                    'output13'))


#check
summary(output_true_seqs)

#find how many unique ASVs
length(unique(o13_true$ASV))
#6483

length(unique(output_true_seqs$ASV))
#7828

#why are there 1400 more ASVs overall than in output 13? 

ggplot(output_true_seqs, aes(x=Sequence_length, fill = Output)) + 
  geom_histogram(binwidth = 3) + 
  theme_linedraw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  facet_grid(.~Output) +
  xlab("Length of ASV sequence") + ylab("Frequency") + 
  ggtitle("Distribution of (consistently detected) ASV lengths retained from different UNIX pipeline parameters") +
  guides(fill="none")


#going to need to make a line graph tracking the reads through each step of the pipeline

#get some summary stats (can use this data if want to add labels to the graph)
avgASVnumtrue <- aggregate(Num.ASVs ~ Output, data = output_true, FUN = mean)

avgASVnumFishtrue <- aggregate(Num.ASVs ~ Output + All.reps + Fish, data = output_true, FUN = mean)

avgdepthtrue <- aggregate(Avg.depth ~ Output + All.reps, data = output_true, FUN = mean)

avgdepthFishtrue <- aggregate(Avg.depth ~ Output + All.reps +Fish, data = output_true, FUN = mean)



#align the seqs----
#make data frame of only ASV and seq then convert to .fasta
UNIXfasta <- output_true_seqs[c(1:2)]
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





