#ASV analysis of SINE UNIX optimization
#1 <-readr::read_tsv("UNIX-optimization/output1/ASV/nochim1.tsv", show_col_types = FALSE)
#2 <-readr::read_tsv("UNIX-optimization/output2/ASV/nochim2.tsv", show_col_types = FALSE)
#3 <-readr::read_tsv("UNIX-optimization/output3/ASV/nochim.tsv", show_col_types = FALSE)
#4 <-readr::read_tsv("UNIX-optimization/output4/ASV/nochim.tsv", show_col_types = FALSE)
#5 <-readr::read_tsv("UNIX-optimization/output5/ASV/nochim.tsv", show_col_types = FALSE)
#6 <-readr::read_tsv("UNIX-optimization/output6/ASV/nochim.tsv", show_col_types = FALSE)
#7 <-readr::read_tsv("UNIX-optimization/output7/ASV/nochim.tsv", show_col_types = FALSE)
#8 <-readr::read_tsv("UNIX-optimization/output8/ASV/nochim.tsv", show_col_types = FALSE)
#9 <-readr::read_tsv("UNIX-optimization/output9/ASV/nochim.tsv", show_col_types = FALSE)
#10 <-readr::read_tsv("UNIX-optimization/output10/ASV/nochim.tsv", show_col_types = FALSE)


#set wd to "/Users/samanthabeal/Documents/MSc/Bioinformatics"
getwd()
#"/Users/samanthabeal"

#setwd("Documents/MSc/Bioinformatics")

#output4----
#load data----
library(readr)
output4 <-readr::read_tsv("UNIX-optimization/output4/ASV/nochim.tsv", show_col_types = FALSE)
#rename columns to only be fish and index
names(output4) = gsub(pattern = "_.*", replacement = "", x = names(output4))

#pull out individual fish----
#190222-A1----
A1190222o4 <- output4[c(1:4)]
names(A1190222o4)[1] <- "ASV"
names(A1190222o4)[2] <- "PCR.1"
names(A1190222o4)[3] <- "PCR.2"
names(A1190222o4)[4] <- "PCR.3"

#determine if the seq is detected in all 3 reps
A1190222o4$Consistently_detected <- ifelse(A1190222o4$PCR.1==0 | A1190222o4$PCR.2== 0 | A1190222o4$PCR.3==0, "false", "true")

#make new df with only the seqs detected in all 3
A1190222o4_true <- subset(A1190222o4, Consistently_detected == "true",
                          select=c("ASV","PCR.1","PCR.2","PCR.3"))
A1190222o4_true$Fish <- "A1190222"
A1190222o4_true$Output <- "output4"

#compare size of two dfs (# ASVs and total reads)
#all ASVs
A1190222o4sum <- as.data.frame(colSums(A1190222o4[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
A1190222o4sum <- as.data.frame(t(A1190222o4sum))

#add column of avg depth across all PCR reps
A1190222o4sum$Avg.depth <-apply(A1190222o4sum,1,mean)
#add columns of total depth
A1190222o4sum$Total.depth <-apply(A1190222o4sum[,1:3],1,sum)

#add column of number of ASVs in the sample (each fish doesnt have reads in all ASVs)
A1190222o4sum$Num.ASVs <- sum(A1190222o4$PCR.1>0 | A1190222o4$PCR.2>0 | A1190222o4$PCR.3>0)

#add column of output folder (so can plot with others for comparison of parameters)
A1190222o4sum$Output <- "output4"

#add column of if in all 3 reps
A1190222o4sum$All.reps <- "No"

#only true ASVs
A1190222o4_truesum <- as.data.frame(colSums(A1190222o4_true[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
A1190222o4_truesum <- as.data.frame(t(A1190222o4_truesum))

#add column of avg depth across all PCR reps
A1190222o4_truesum $Avg.depth <-apply(A1190222o4_truesum,1,mean)

#add columns of total depth
A1190222o4_truesum $Total.depth <-apply(A1190222o4_truesum[,1:3],1,sum)

#add column of number of ASVs in the sample
A1190222o4_truesum$Num.ASVs <- sum(A1190222o4_true$PCR.1>0 | A1190222o4_true$PCR.2>0 | A1190222o4_true$PCR.3>0)

#add column of pipeline (so can plot with UNIX for comparison)
A1190222o4_truesum$Output <- "output4"

#add column of if in all 3 reps
A1190222o4_truesum$All.reps <- "Yes"

#190922-J2 add data frames together
A1190222_output4 <- rbind(A1190222o4_truesum, A1190222o4sum)

#add column of fish ID
A1190222_output4$Fish <- "190922-A1"

#190222-A2----
A2190222o4 <- output4[c(1,5:7)]
names(A2190222o4)[1] <- "ASV"
names(A2190222o4)[2] <- "PCR.1"
names(A2190222o4)[3] <- "PCR.2"
names(A2190222o4)[4] <- "PCR.3"

#determine if the seq is detected in all 3 reps
A2190222o4$Consistently_detected <- ifelse(A2190222o4$PCR.1==0 | A2190222o4$PCR.2== 0 | A2190222o4$PCR.3==0, "false", "true")

#make new df with only the seqs detected in all 3
A2190222o4_true <- subset(A2190222o4, Consistently_detected == "true",
                          select=c("ASV","PCR.1","PCR.2","PCR.3"))
A2190222o4_true$Fish <- "A2190222"
A2190222o4_true$Output <- "output4"

#compare size of two dfs (# ASVs and total reads)
#all ASVs
A2190222o4sum <- as.data.frame(colSums(A2190222o4[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
A2190222o4sum <- as.data.frame(t(A2190222o4sum))

#add column of avg depth across all PCR reps
A2190222o4sum$Avg.depth <-apply(A2190222o4sum,1,mean)
#add columns of total depth
A2190222o4sum$Total.depth <-apply(A2190222o4sum[,1:3],1,sum)

#add column of number of ASVs in the sample (each fish doesnt have reads in all ASVs)
A2190222o4sum$Num.ASVs <- sum(A2190222o4$PCR.1>0 | A2190222o4$PCR.2>0 | A2190222o4$PCR.3>0)

#add column of output folder (so can plot with others for comparison of parameters)
A2190222o4sum$Output <- "output4"

#add column of if in all 3 reps
A2190222o4sum$All.reps <- "No"

#only true ASVs
A2190222o4_truesum <- as.data.frame(colSums(A2190222o4_true[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
A2190222o4_truesum <- as.data.frame(t(A2190222o4_truesum))

#add column of avg depth across all PCR reps
A2190222o4_truesum $Avg.depth <-apply(A2190222o4_truesum,1,mean)

#add columns of total depth
A2190222o4_truesum $Total.depth <-apply(A2190222o4_truesum[,1:3],1,sum)

#add column of number of ASVs in the sample
A2190222o4_truesum$Num.ASVs <- sum(A2190222o4_true$PCR.1>0 | A2190222o4_true$PCR.2>0 | A2190222o4_true$PCR.3>0)

#add column of pipeline (so can plot with UNIX for comparison)
A2190222o4_truesum$Output <- "output4"

#add column of if in all 3 reps
A2190222o4_truesum$All.reps <- "Yes"

#190922-J2 add data frames together
A2190222_output4 <- rbind(A2190222o4_truesum, A2190222o4sum)

#add column of fish ID
A2190222_output4$Fish <- "190922-A2"

#190222-H1----
H1190222o4 <- output4[c(1,8:10)]
names(H1190222o4)[1] <- "ASV"
names(H1190222o4)[2] <- "PCR.1"
names(H1190222o4)[3] <- "PCR.2"
names(H1190222o4)[4] <- "PCR.3"

#determine if the seq is detected in all 3 reps
H1190222o4$Consistently_detected <- ifelse(H1190222o4$PCR.1==0 | H1190222o4$PCR.2== 0 | H1190222o4$PCR.3==0, "false", "true")

#make new df with only the seqs detected in all 3
H1190222o4_true <- subset(H1190222o4, Consistently_detected == "true",
                          select=c("ASV","PCR.1","PCR.2","PCR.3"))
H1190222o4_true$Fish <- "H1190222"
H1190222o4_true$Output <- "output4"

#compare size of two dfs (# ASVs and total reads)
#all ASVs
H1190222o4sum <- as.data.frame(colSums(H1190222o4[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
H1190222o4sum <- as.data.frame(t(H1190222o4sum))

#add column of avg depth across all PCR reps
H1190222o4sum$Avg.depth <-apply(H1190222o4sum,1,mean)
#add columns of total depth
H1190222o4sum$Total.depth <-apply(H1190222o4sum[,1:3],1,sum)

#add column of number of ASVs in the sample (each fish doesnt have reads in all ASVs)
H1190222o4sum$Num.ASVs <- sum(H1190222o4$PCR.1>0 | H1190222o4$PCR.2>0 | H1190222o4$PCR.3>0)

#add column of output folder (so can plot with others for comparison of parameters)
H1190222o4sum$Output <- "output4"

#add column of if in all 3 reps
H1190222o4sum$All.reps <- "No"

#only true ASVs
H1190222o4_truesum <- as.data.frame(colSums(H1190222o4_true[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
H1190222o4_truesum <- as.data.frame(t(H1190222o4_truesum))

#add column of avg depth across all PCR reps
H1190222o4_truesum $Avg.depth <-apply(H1190222o4_truesum,1,mean)

#add columns of total depth
H1190222o4_truesum $Total.depth <-apply(H1190222o4_truesum[,1:3],1,sum)

#add column of number of ASVs in the sample
H1190222o4_truesum$Num.ASVs <- sum(H1190222o4_true$PCR.1>0 | H1190222o4_true$PCR.2>0 | H1190222o4_true$PCR.3>0)

#add column of pipeline (so can plot with UNIX for comparison)
H1190222o4_truesum$Output <- "output4"

#add column of if in all 3 reps
H1190222o4_truesum$All.reps <- "Yes"

#190922-J2 add data frames together
H1190222_output4 <- rbind(H1190222o4_truesum, H1190222o4sum)

#add column of fish ID
H1190222_output4$Fish <- "190922-H1"

#190222-M1----
M1190222o4 <- output4[c(1,11:13)]
names(M1190222o4)[1] <- "ASV"
names(M1190222o4)[2] <- "PCR.1"
names(M1190222o4)[3] <- "PCR.2"
names(M1190222o4)[4] <- "PCR.3"

#determine if the seq is detected in all 3 reps
M1190222o4$Consistently_detected <- ifelse(M1190222o4$PCR.1==0 | M1190222o4$PCR.2== 0 | M1190222o4$PCR.3==0, "false", "true")

#make new df with only the seqs detected in all 3
M1190222o4_true <- subset(M1190222o4, Consistently_detected == "true",
                          select=c("ASV","PCR.1","PCR.2","PCR.3"))
M1190222o4_true$Fish <- "M1190222"
M1190222o4_true$Output <- "output4"


#compare size of two dfs (# ASVs and total reads)
#all ASVs
M1190222o4sum <- as.data.frame(colSums(M1190222o4[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
M1190222o4sum <- as.data.frame(t(M1190222o4sum))

#add column of avg depth across all PCR reps
M1190222o4sum$Avg.depth <-apply(M1190222o4sum,1,mean)
#add columns of total depth
M1190222o4sum$Total.depth <-apply(M1190222o4sum[,1:3],1,sum)

#add column of number of ASVs in the sample (each fish doesnt have reads in all ASVs)
M1190222o4sum$Num.ASVs <- sum(M1190222o4$PCR.1>0 | M1190222o4$PCR.2>0 | M1190222o4$PCR.3>0)

#add column of output folder (so can plot with others for comparison of parameters)
M1190222o4sum$Output <- "output4"

#add column of if in all 3 reps
M1190222o4sum$All.reps <- "No"

#only true ASVs
M1190222o4_truesum <- as.data.frame(colSums(M1190222o4_true[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
M1190222o4_truesum <- as.data.frame(t(M1190222o4_truesum))

#add column of avg depth across all PCR reps
M1190222o4_truesum $Avg.depth <-apply(M1190222o4_truesum,1,mean)

#add columns of total depth
M1190222o4_truesum $Total.depth <-apply(M1190222o4_truesum[,1:3],1,sum)

#add column of number of ASVs in the sample
M1190222o4_truesum$Num.ASVs <- sum(M1190222o4_true$PCR.1>0 | M1190222o4_true$PCR.2>0 | M1190222o4_true$PCR.3>0)

#add column of pipeline (so can plot with UNIX for comparison)
M1190222o4_truesum$Output <- "output4"

#add column of if in all 3 reps
M1190222o4_truesum$All.reps <- "Yes"

#190922-J2 add data frames together
M1190222_output4 <- rbind(M1190222o4_truesum, M1190222o4sum)

#add column of fish ID
M1190222_output4$Fish <- "190922-M1"

#190918-A1----
A1190918o4 <- output4[c(1,14:16)]
names(A1190918o4)[1] <- "ASV"
names(A1190918o4)[2] <- "PCR.1"
names(A1190918o4)[3] <- "PCR.2"
names(A1190918o4)[4] <- "PCR.3"

#determine if the seq is detected in all 3 reps
A1190918o4$Consistently_detected <- ifelse(A1190918o4$PCR.1==0 | A1190918o4$PCR.2== 0 | A1190918o4$PCR.3==0, "false", "true")

#make new df with only the seqs detected in all 3
A1190918o4_true <- subset(A1190918o4, Consistently_detected == "true",
                          select=c("ASV","PCR.1","PCR.2","PCR.3"))
A1190918o4_true$Fish <- "A1190918"
A1190918o4_true$Output <- "output4"
#compare size of two dfs (# ASVs and total reads)
#all ASVs
A1190918o4sum <- as.data.frame(colSums(A1190918o4[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
A1190918o4sum <- as.data.frame(t(A1190918o4sum))

#add column of avg depth across all PCR reps
A1190918o4sum$Avg.depth <-apply(A1190918o4sum,1,mean)
#add columns of total depth
A1190918o4sum$Total.depth <-apply(A1190918o4sum[,1:3],1,sum)

#add column of number of ASVs in the sample (each fish doesnt have reads in all ASVs)
A1190918o4sum$Num.ASVs <- sum(A1190918o4$PCR.1>0 | A1190918o4$PCR.2>0 | A1190918o4$PCR.3>0)

#add column of output folder (so can plot with others for comparison of parameters)
A1190918o4sum$Output <- "output4"

#add column of if in all 3 reps
A1190918o4sum$All.reps <- "No"

#only true ASVs
A1190918o4_truesum <- as.data.frame(colSums(A1190918o4_true[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
A1190918o4_truesum <- as.data.frame(t(A1190918o4_truesum))

#add column of avg depth across all PCR reps
A1190918o4_truesum $Avg.depth <-apply(A1190918o4_truesum,1,mean)

#add columns of total depth
A1190918o4_truesum $Total.depth <-apply(A1190918o4_truesum[,1:3],1,sum)

#add column of number of ASVs in the sample
A1190918o4_truesum$Num.ASVs <- sum(A1190918o4_true$PCR.1>0 | A1190918o4_true$PCR.2>0 | A1190918o4_true$PCR.3>0)

#add column of pipeline (so can plot with UNIX for comparison)
A1190918o4_truesum$Output <- "output4"

#add column of if in all 3 reps
A1190918o4_truesum$All.reps <- "Yes"

#190922-J2 add data frames together
A1190918_output4 <- rbind(A1190918o4_truesum, A1190918o4sum)

#add column of fish ID
A1190918_output4$Fish <- "190918-A1"

#190918-A3----
A3190918o4 <- output4[c(1,17:19)]
names(A3190918o4)[1] <- "ASV"
names(A3190918o4)[2] <- "PCR.1"
names(A3190918o4)[3] <- "PCR.2"
names(A3190918o4)[4] <- "PCR.3"

#determine if the seq is detected in all 3 reps
A3190918o4$Consistently_detected <- ifelse(A3190918o4$PCR.1==0 | A3190918o4$PCR.2== 0 | A3190918o4$PCR.3==0, "false", "true")

#make new df with only the seqs detected in all 3
A3190918o4_true <- subset(A3190918o4, Consistently_detected == "true",
                          select=c("ASV","PCR.1","PCR.2","PCR.3"))
A3190918o4_true$Fish <- "A3190918"
A3190918o4_true$Output <- "output4"

#compare size of two dfs (# ASVs and total reads)
#all ASVs
A3190918o4sum <- as.data.frame(colSums(A3190918o4[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
A3190918o4sum <- as.data.frame(t(A3190918o4sum))

#add column of avg depth across all PCR reps
A3190918o4sum$Avg.depth <-apply(A3190918o4sum,1,mean)
#add columns of total depth
A3190918o4sum$Total.depth <-apply(A3190918o4sum[,1:3],1,sum)

#add column of number of ASVs in the sample (each fish doesnt have reads in all ASVs)
A3190918o4sum$Num.ASVs <- sum(A3190918o4$PCR.1>0 | A3190918o4$PCR.2>0 | A3190918o4$PCR.3>0)

#add column of output folder (so can plot with others for comparison of parameters)
A3190918o4sum$Output <- "output4"

#add column of if in all 3 reps
A3190918o4sum$All.reps <- "No"

#only true ASVs
A3190918o4_truesum <- as.data.frame(colSums(A3190918o4_true[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
A3190918o4_truesum <- as.data.frame(t(A3190918o4_truesum))

#add column of avg depth across all PCR reps
A3190918o4_truesum $Avg.depth <-apply(A3190918o4_truesum,1,mean)

#add columns of total depth
A3190918o4_truesum $Total.depth <-apply(A3190918o4_truesum[,1:3],1,sum)

#add column of number of ASVs in the sample
A3190918o4_truesum$Num.ASVs <- sum(A3190918o4_true$PCR.1>0 | A3190918o4_true$PCR.2>0 | A3190918o4_true$PCR.3>0)

#add column of pipeline (so can plot with UNIX for comparison)
A3190918o4_truesum$Output <- "output4"

#add column of if in all 3 reps
A3190918o4_truesum$All.reps <- "Yes"

#190922-J2 add data frames together
A3190918_output4 <- rbind(A3190918o4_truesum, A3190918o4sum)

#add column of fish ID
A3190918_output4$Fish <- "190918-A3"

#190918-A4----
A4190918o4 <- output4[c(1,20:22)]
names(A4190918o4)[1] <- "ASV"
names(A4190918o4)[2] <- "PCR.1"
names(A4190918o4)[3] <- "PCR.2"
names(A4190918o4)[4] <- "PCR.3"

#determine if the seq is detected in all 3 reps
A4190918o4$Consistently_detected <- ifelse(A4190918o4$PCR.1==0 | A4190918o4$PCR.2== 0 | A4190918o4$PCR.3==0, "false", "true")

#make new df with only the seqs detected in all 3
A4190918o4_true <- subset(A4190918o4, Consistently_detected == "true",
                          select=c("ASV","PCR.1","PCR.2","PCR.3"))
A4190918o4_true$Fish <- "A4190918"
A4190918o4_true$Output <- "output4"


#compare size of two dfs (# ASVs and total reads)
#all ASVs
A4190918o4sum <- as.data.frame(colSums(A4190918o4[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
A4190918o4sum <- as.data.frame(t(A4190918o4sum))

#add column of avg depth across all PCR reps
A4190918o4sum$Avg.depth <-apply(A4190918o4sum,1,mean)
#add columns of total depth
A4190918o4sum$Total.depth <-apply(A4190918o4sum[,1:3],1,sum)

#add column of number of ASVs in the sample (each fish doesnt have reads in all ASVs)
A4190918o4sum$Num.ASVs <- sum(A4190918o4$PCR.1>0 | A4190918o4$PCR.2>0 | A4190918o4$PCR.3>0)

#add column of output folder (so can plot with others for comparison of parameters)
A4190918o4sum$Output <- "output4"

#add column of if in all 3 reps
A4190918o4sum$All.reps <- "No"

#only true ASVs
A4190918o4_truesum <- as.data.frame(colSums(A4190918o4_true[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
A4190918o4_truesum <- as.data.frame(t(A4190918o4_truesum))

#add column of avg depth across all PCR reps
A4190918o4_truesum $Avg.depth <-apply(A4190918o4_truesum,1,mean)

#add columns of total depth
A4190918o4_truesum $Total.depth <-apply(A4190918o4_truesum[,1:3],1,sum)

#add column of number of ASVs in the sample
A4190918o4_truesum$Num.ASVs <- sum(A4190918o4_true$PCR.1>0 | A4190918o4_true$PCR.2>0 | A4190918o4_true$PCR.3>0)

#add column of pipeline (so can plot with UNIX for comparison)
A4190918o4_truesum$Output <- "output4"

#add column of if in all 3 reps
A4190918o4_truesum$All.reps <- "Yes"

#190922-J2 add data frames together
A4190918_output4 <- rbind(A4190918o4_truesum, A4190918o4sum)

#add column of fish ID
A4190918_output4$Fish <- "190918-A4"

#190918-C4----
C4190918o4 <- output4[c(1,23:25)]
names(C4190918o4)[1] <- "ASV"
names(C4190918o4)[2] <- "PCR.1"
names(C4190918o4)[3] <- "PCR.2"
names(C4190918o4)[4] <- "PCR.3"

#determine if the seq is detected in all 3 reps
C4190918o4$Consistently_detected <- ifelse(C4190918o4$PCR.1==0 | C4190918o4$PCR.2== 0 | C4190918o4$PCR.3==0, "false", "true")

#make new df with only the seqs detected in all 3
C4190918o4_true <- subset(C4190918o4, Consistently_detected == "true",
                          select=c("ASV","PCR.1","PCR.2","PCR.3"))
C4190918o4_true$Fish <- "C4190918"
C4190918o4_true$Output <- "output4"


#compare size of two dfs (# ASVs and total reads)
#all ASVs
C4190918o4sum <- as.data.frame(colSums(C4190918o4[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
C4190918o4sum <- as.data.frame(t(C4190918o4sum))

#add column of avg depth across all PCR reps
C4190918o4sum$Avg.depth <-apply(C4190918o4sum,1,mean)
#add columns of total depth
C4190918o4sum$Total.depth <-apply(C4190918o4sum[,1:3],1,sum)

#add column of number of ASVs in the sample (each fish doesnt have reads in all ASVs)
C4190918o4sum$Num.ASVs <- sum(C4190918o4$PCR.1>0 | C4190918o4$PCR.2>0 | C4190918o4$PCR.3>0)

#add column of output folder (so can plot with others for comparison of parameters)
C4190918o4sum$Output <- "output4"

#add column of if in all 3 reps
C4190918o4sum$All.reps <- "No"

#only true ASVs
C4190918o4_truesum <- as.data.frame(colSums(C4190918o4_true[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
C4190918o4_truesum <- as.data.frame(t(C4190918o4_truesum))

#add column of avg depth across all PCR reps
C4190918o4_truesum $Avg.depth <-apply(C4190918o4_truesum,1,mean)

#add columns of total depth
C4190918o4_truesum $Total.depth <-apply(C4190918o4_truesum[,1:3],1,sum)

#add column of number of ASVs in the sample
C4190918o4_truesum$Num.ASVs <- sum(C4190918o4_true$PCR.1>0 | C4190918o4_true$PCR.2>0 | C4190918o4_true$PCR.3>0)

#add column of pipeline (so can plot with UNIX for comparison)
C4190918o4_truesum$Output <- "output4"

#add column of if in all 3 reps
C4190918o4_truesum$All.reps <- "Yes"

#190922-J2 add data frames together
C4190918_output4 <- rbind(C4190918o4_truesum, C4190918o4sum)

#add column of fish ID
C4190918_output4$Fish <- "190918-C4"

#190918-D4----
D4190918o4 <- output4[c(1,26:28)]
names(D4190918o4)[1] <- "ASV"
names(D4190918o4)[2] <- "PCR.1"
names(D4190918o4)[3] <- "PCR.2"
names(D4190918o4)[4] <- "PCR.3"

#determine if the seq is detected in all 3 reps
D4190918o4$Consistently_detected <- ifelse(D4190918o4$PCR.1==0 | D4190918o4$PCR.2== 0 | D4190918o4$PCR.3==0, "false", "true")

#make new df with only the seqs detected in all 3
D4190918o4_true <- subset(D4190918o4, Consistently_detected == "true",
                          select=c("ASV","PCR.1","PCR.2","PCR.3"))
D4190918o4_true$Fish <- "D4190918"
D4190918o4_true$Output <- "output4"


#compare size of two dfs (# ASVs and total reads)
#all ASVs
D4190918o4sum <- as.data.frame(colSums(D4190918o4[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
D4190918o4sum <- as.data.frame(t(D4190918o4sum))

#add column of avg depth across all PCR reps
D4190918o4sum$Avg.depth <-apply(D4190918o4sum,1,mean)
#add columns of total depth
D4190918o4sum$Total.depth <-apply(D4190918o4sum[,1:3],1,sum)

#add column of number of ASVs in the sample (each fish doesnt have reads in all ASVs)
D4190918o4sum$Num.ASVs <- sum(D4190918o4$PCR.1>0 | D4190918o4$PCR.2>0 | D4190918o4$PCR.3>0)

#add column of output folder (so can plot with others for comparison of parameters)
D4190918o4sum$Output <- "output4"

#add column of if in all 3 reps
D4190918o4sum$All.reps <- "No"

#only true ASVs
D4190918o4_truesum <- as.data.frame(colSums(D4190918o4_true[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
D4190918o4_truesum <- as.data.frame(t(D4190918o4_truesum))

#add column of avg depth across all PCR reps
D4190918o4_truesum$Avg.depth <-apply(D4190918o4_truesum,1,mean)

#add columns of total depth
D4190918o4_truesum$Total.depth <-apply(D4190918o4_truesum[,1:3],1,sum)

#add column of number of ASVs in the sample
D4190918o4_truesum$Num.ASVs <- sum(D4190918o4_true$PCR.1>0 | D4190918o4_true$PCR.2>0 | D4190918o4_true$PCR.3>0)

#add column of output (so can plot with UNIX for comparison)
D4190918o4_truesum$Output <- "output4"

#add column of if in all 3 reps
D4190918o4_truesum$All.reps <- "Yes"

#190922-J2 add data frames together
D4190918_output4 <- rbind(D4190918o4_truesum, D4190918o4sum)

#add column of fish ID
D4190918_output4$Fish <- "190918-D4"

#190918-E7----
E7190918o4 <- output4[c(1,29:31)]
names(E7190918o4)[1] <- "ASV"
names(E7190918o4)[2] <- "PCR.1"
names(E7190918o4)[3] <- "PCR.2"
names(E7190918o4)[4] <- "PCR.3"

#determine if the seq is detected in all 3 reps
E7190918o4$Consistently_detected <- ifelse(E7190918o4$PCR.1==0 | E7190918o4$PCR.2== 0 | E7190918o4$PCR.3==0, "false", "true")

#make new df with only the seqs detected in all 3
E7190918o4_true <- subset(E7190918o4, Consistently_detected == "true",
                          select=c("ASV","PCR.1","PCR.2","PCR.3"))
E7190918o4_true$Fish <- "E7190918"
E7190918o4_true$Output <- "output4"

#compare size of two dfs (# ASVs and total reads)
#all ASVs
E7190918o4sum <- as.data.frame(colSums(E7190918o4[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
E7190918o4sum <- as.data.frame(t(E7190918o4sum))

#add column of avg depth across all PCR reps
E7190918o4sum$Avg.depth <-apply(E7190918o4sum,1,mean)
#add columns of total depth
E7190918o4sum$Total.depth <-apply(E7190918o4sum[,1:3],1,sum)

#add column of number of ASVs in the sample (each fish doesnt have reads in all ASVs)
E7190918o4sum$Num.ASVs <- sum(E7190918o4$PCR.1>0 | E7190918o4$PCR.2>0 | E7190918o4$PCR.3>0)

#add column of output folder (so can plot with others for comparison of parameters)
E7190918o4sum$Output <- "output4"

#add column of if in all 3 reps
E7190918o4sum$All.reps <- "No"

#only true ASVs
E7190918o4_truesum <- as.data.frame(colSums(E7190918o4_true[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
E7190918o4_truesum <- as.data.frame(t(E7190918o4_truesum))

#add column of avg depth across all PCR reps
E7190918o4_truesum $Avg.depth <-apply(E7190918o4_truesum,1,mean)

#add columns of total depth
E7190918o4_truesum $Total.depth <-apply(E7190918o4_truesum[,1:3],1,sum)

#add column of number of ASVs in the sample
E7190918o4_truesum$Num.ASVs <- sum(E7190918o4_true$PCR.1>0 | E7190918o4_true$PCR.2>0 | E7190918o4_true$PCR.3>0)

#add column of pipeline (so can plot with UNIX for comparison)
E7190918o4_truesum$Output <- "output4"

#add column of if in all 3 reps
E7190918o4_truesum$All.reps <- "Yes"

#190922-J2 add data frames together
E7190918_output4 <- rbind(E7190918o4_truesum, E7190918o4sum)

#add column of fish ID
E7190918_output4$Fish <- "190918-E7"

#190918-F2----
F2190918o4 <- output4[c(1,32:34)]
names(F2190918o4)[1] <- "ASV"
names(F2190918o4)[2] <- "PCR.1"
names(F2190918o4)[3] <- "PCR.2"
names(F2190918o4)[4] <- "PCR.3"

#determine if the seq is detected in all 3 reps
F2190918o4$Consistently_detected <- ifelse(F2190918o4$PCR.1==0 | F2190918o4$PCR.2== 0 | F2190918o4$PCR.3==0, "false", "true")

#make new df with only the seqs detected in all 3
F2190918o4_true <- subset(F2190918o4, Consistently_detected == "true",
                          select=c("ASV","PCR.1","PCR.2","PCR.3"))
F2190918o4_true$Fish <- "F2190918"
F2190918o4_true$Output <- "output4"

#compare size of two dfs (# ASVs and total reads)
#all ASVs
F2190918o4sum <- as.data.frame(colSums(F2190918o4[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
F2190918o4sum <- as.data.frame(t(F2190918o4sum))

#add column of avg depth across all PCR reps
F2190918o4sum$Avg.depth <-apply(F2190918o4sum,1,mean)
#add columns of total depth
F2190918o4sum$Total.depth <-apply(F2190918o4sum[,1:3],1,sum)

#add column of number of ASVs in the sample (each fish doesnt have reads in all ASVs)
F2190918o4sum$Num.ASVs <- sum(F2190918o4$PCR.1>0 | F2190918o4$PCR.2>0 | F2190918o4$PCR.3>0)

#add column of output folder (so can plot with others for comparison of parameters)
F2190918o4sum$Output <- "output4"

#add column of if in all 3 reps
F2190918o4sum$All.reps <- "No"

#only true ASVs
F2190918o4_truesum <- as.data.frame(colSums(F2190918o4_true[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
F2190918o4_truesum <- as.data.frame(t(F2190918o4_truesum))

#add column of avg depth across all PCR reps
F2190918o4_truesum $Avg.depth <-apply(F2190918o4_truesum,1,mean)

#add columns of total depth
F2190918o4_truesum $Total.depth <-apply(F2190918o4_truesum[,1:3],1,sum)

#add column of number of ASVs in the sample
F2190918o4_truesum$Num.ASVs <- sum(F2190918o4_true$PCR.1>0 | F2190918o4_true$PCR.2>0 | F2190918o4_true$PCR.3>0)

#add column of pipeline (so can plot with UNIX for comparison)
F2190918o4_truesum$Output <- "output4"

#add column of if in all 3 reps
F2190918o4_truesum$All.reps <- "Yes"

#190922-J2 add data frames together
F2190918_output4 <- rbind(F2190918o4_truesum, F2190918o4sum)

#add column of fish ID
F2190918_output4$Fish <- "190918-F2"

#190922-D2----
D2190922o4 <- output4[c(1,35:37)]
names(D2190922o4)[1] <- "ASV"
names(D2190922o4)[2] <- "PCR.1"
names(D2190922o4)[3] <- "PCR.2"
names(D2190922o4)[4] <- "PCR.3"

#determine if the seq is detected in all 3 reps
D2190922o4$Consistently_detected <- ifelse(D2190922o4$PCR.1==0 | D2190922o4$PCR.2== 0 | D2190922o4$PCR.3==0, "false", "true")

#make new df with only the seqs detected in all 3
D2190922o4_true <- subset(D2190922o4, Consistently_detected == "true",
                          select=c("ASV","PCR.1","PCR.2","PCR.3"))
D2190922o4_true$Fish <- "D2190922"
D2190922o4_true$Output <- "output4"

#compare size of two dfs (# ASVs and total reads)
#all ASVs
D2190922o4sum <- as.data.frame(colSums(D2190922o4[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
D2190922o4sum <- as.data.frame(t(D2190922o4sum))

#add column of avg depth across all PCR reps
D2190922o4sum$Avg.depth <-apply(D2190922o4sum,1,mean)
#add columns of total depth
D2190922o4sum$Total.depth <-apply(D2190922o4sum[,1:3],1,sum)

#add column of number of ASVs in the sample (each fish doesnt have reads in all ASVs)
D2190922o4sum$Num.ASVs <- sum(D2190922o4$PCR.1>0 | D2190922o4$PCR.2>0 | D2190922o4$PCR.3>0)

#add column of output folder (so can plot with others for comparison of parameters)
D2190922o4sum$Output <- "output4"

#add column of if in all 3 reps
D2190922o4sum$All.reps <- "No"

#only true ASVs
D2190922o4_truesum <- as.data.frame(colSums(D2190922o4_true[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
D2190922o4_truesum <- as.data.frame(t(D2190922o4_truesum))

#add column of avg depth across all PCR reps
D2190922o4_truesum $Avg.depth <-apply(D2190922o4_truesum,1,mean)

#add columns of total depth
D2190922o4_truesum $Total.depth <-apply(D2190922o4_truesum[,1:3],1,sum)

#add column of number of ASVs in the sample
D2190922o4_truesum$Num.ASVs <- sum(D2190922o4_true$PCR.1>0 | D2190922o4_true$PCR.2>0 | D2190922o4_true$PCR.3>0)

#add column of pipeline (so can plot with UNIX for comparison)
D2190922o4_truesum$Output <- "output4"

#add column of if in all 3 reps
D2190922o4_truesum$All.reps <- "Yes"

#190922-J2 add data frames together
D2190922_output4 <- rbind(D2190922o4_truesum, D2190922o4sum)

#add column of fish ID
D2190922_output4$Fish <- "190922-D2"

#190922-F2----
F2190922o4 <- output4[c(1,38:40)]
names(F2190922o4)[1] <- "ASV"
names(F2190922o4)[2] <- "PCR.1"
names(F2190922o4)[3] <- "PCR.2"
names(F2190922o4)[4] <- "PCR.3"

#determine if the seq is detected in all 3 reps
F2190922o4$Consistently_detected <- ifelse(F2190922o4$PCR.1==0 | F2190922o4$PCR.2== 0 | F2190922o4$PCR.3==0, "false", "true")

#make new df with only the seqs detected in all 3
F2190922o4_true <- subset(F2190922o4, Consistently_detected == "true",
                          select=c("ASV","PCR.1","PCR.2","PCR.3"))
F2190922o4_true$Fish <- "F2190922"
F2190922o4_true$Output <- "output4"

#compare size of two dfs (# ASVs and total reads)
#all ASVs
F2190922o4sum <- as.data.frame(colSums(F2190922o4[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
F2190922o4sum <- as.data.frame(t(F2190922o4sum))

#add column of avg depth across all PCR reps
F2190922o4sum$Avg.depth <-apply(F2190922o4sum,1,mean)
#add columns of total depth
F2190922o4sum$Total.depth <-apply(F2190922o4sum[,1:3],1,sum)

#add column of number of ASVs in the sample (each fish doesnt have reads in all ASVs)
F2190922o4sum$Num.ASVs <- sum(F2190922o4$PCR.1>0 | F2190922o4$PCR.2>0 | F2190922o4$PCR.3>0)

#add column of output folder (so can plot with others for comparison of parameters)
F2190922o4sum$Output <- "output4"

#add column of if in all 3 reps
F2190922o4sum$All.reps <- "No"

#only true ASVs
F2190922o4_truesum <- as.data.frame(colSums(F2190922o4_true[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
F2190922o4_truesum <- as.data.frame(t(F2190922o4_truesum))

#add column of avg depth across all PCR reps
F2190922o4_truesum $Avg.depth <-apply(F2190922o4_truesum,1,mean)

#add columns of total depth
F2190922o4_truesum $Total.depth <-apply(F2190922o4_truesum[,1:3],1,sum)

#add column of number of ASVs in the sample
F2190922o4_truesum$Num.ASVs <- sum(F2190922o4_true$PCR.1>0 | F2190922o4_true$PCR.2>0 | F2190922o4_true$PCR.3>0)

#add column of pipeline (so can plot with UNIX for comparison)
F2190922o4_truesum$Output <- "output4"

#add column of if in all 3 reps
F2190922o4_truesum$All.reps <- "Yes"

#190922-J2 add data frames together
F2190922_output4 <- rbind(F2190922o4_truesum, F2190922o4sum)

#add column of fish ID
F2190922_output4$Fish <- "190922-F2"

#190922-G2----
G2190922o4 <- output4[c(1,41:43)]
names(G2190922o4)[1] <- "ASV"
names(G2190922o4)[2] <- "PCR.1"
names(G2190922o4)[3] <- "PCR.2"
names(G2190922o4)[4] <- "PCR.3"

#determine if the seq is detected in all 3 reps
G2190922o4$Consistently_detected <- ifelse(G2190922o4$PCR.1==0 | G2190922o4$PCR.2== 0 | G2190922o4$PCR.3==0, "false", "true")

#make new df with only the seqs detected in all 3
G2190922o4_true <- subset(G2190922o4, Consistently_detected == "true",
                          select=c("ASV","PCR.1","PCR.2","PCR.3"))
G2190922o4_true$Fish <- "G2190922"
G2190922o4_true$Output <- "output4"

#compare size of two dfs (# ASVs and total reads)
#all ASVs
G2190922o4sum <- as.data.frame(colSums(G2190922o4[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
G2190922o4sum <- as.data.frame(t(G2190922o4sum))

#add column of avg depth across all PCR reps
G2190922o4sum$Avg.depth <-apply(G2190922o4sum,1,mean)
#add columns of total depth
G2190922o4sum$Total.depth <-apply(G2190922o4sum[,1:3],1,sum)

#add column of number of ASVs in the sample (each fish doesnt have reads in all ASVs)
G2190922o4sum$Num.ASVs <- sum(G2190922o4$PCR.1>0 | G2190922o4$PCR.2>0 | G2190922o4$PCR.3>0)

#add column of output folder (so can plot with others for comparison of parameters)
G2190922o4sum$Output <- "output4"

#add column of if in all 3 reps
G2190922o4sum$All.reps <- "No"

#only true ASVs
G2190922o4_truesum <- as.data.frame(colSums(G2190922o4_true[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
G2190922o4_truesum <- as.data.frame(t(G2190922o4_truesum))

#add column of avg depth across all PCR reps
G2190922o4_truesum $Avg.depth <-apply(G2190922o4_truesum,1,mean)

#add columns of total depth
G2190922o4_truesum $Total.depth <-apply(G2190922o4_truesum[,1:3],1,sum)

#add column of number of ASVs in the sample
G2190922o4_truesum$Num.ASVs <- sum(G2190922o4_true$PCR.1>0 | G2190922o4_true$PCR.2>0 | G2190922o4_true$PCR.3>0)

#add column of pipeline (so can plot with UNIX for comparison)
G2190922o4_truesum$Output <- "output4"

#add column of if in all 3 reps
G2190922o4_truesum$All.reps <- "Yes"

#190922-J2 add data frames together
G2190922_output4 <- rbind(G2190922o4_truesum, G2190922o4sum)

#add column of fish ID
G2190922_output4$Fish <- "190922-G2"

#190922-H2----
H2190922o4 <- output4[c(1,44:46)]
names(H2190922o4)[1] <- "ASV"
names(H2190922o4)[2] <- "PCR.1"
names(H2190922o4)[3] <- "PCR.2"
names(H2190922o4)[4] <- "PCR.3"

#determine if the seq is detected in all 3 reps
H2190922o4$Consistently_detected <- ifelse(H2190922o4$PCR.1==0 | H2190922o4$PCR.2== 0 | H2190922o4$PCR.3==0, "false", "true")

#make new df with only the seqs detected in all 3
H2190922o4_true <- subset(H2190922o4, Consistently_detected == "true",
                          select=c("ASV","PCR.1","PCR.2","PCR.3"))
H2190922o4_true$Fish <- "H2190922"
H2190922o4_true$Output <- "output4"

#compare size of two dfs (# ASVs and total reads)
#all ASVs
H2190922o4sum <- as.data.frame(colSums(H2190922o4[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
H2190922o4sum <- as.data.frame(t(H2190922o4sum))

#add column of avg depth across all PCR reps
H2190922o4sum$Avg.depth <-apply(H2190922o4sum,1,mean)
#add columns of total depth
H2190922o4sum$Total.depth <-apply(H2190922o4sum[,1:3],1,sum)

#add column of number of ASVs in the sample (each fish doesnt have reads in all ASVs)
H2190922o4sum$Num.ASVs <- sum(H2190922o4$PCR.1>0 | H2190922o4$PCR.2>0 | H2190922o4$PCR.3>0)

#add column of output folder (so can plot with others for comparison of parameters)
H2190922o4sum$Output <- "output4"

#add column of if in all 3 reps
H2190922o4sum$All.reps <- "No"

#only true ASVs
H2190922o4_truesum <- as.data.frame(colSums(H2190922o4_true[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
H2190922o4_truesum <- as.data.frame(t(H2190922o4_truesum))

#add column of avg depth across all PCR reps
H2190922o4_truesum $Avg.depth <-apply(H2190922o4_truesum,1,mean)

#add columns of total depth
H2190922o4_truesum $Total.depth <-apply(H2190922o4_truesum[,1:3],1,sum)

#add column of number of ASVs in the sample
H2190922o4_truesum$Num.ASVs <- sum(H2190922o4_true$PCR.1>0 | H2190922o4_true$PCR.2>0 | H2190922o4_true$PCR.3>0)

#add column of pipeline (so can plot with UNIX for comparison)
H2190922o4_truesum$Output <- "output4"

#add column of if in all 3 reps
H2190922o4_truesum$All.reps <- "Yes"

#190922-J2 add data frames together
H2190922_output4 <- rbind(H2190922o4_truesum, H2190922o4sum)

#add column of fish ID
H2190922_output4$Fish <- "190922-H2"

#190922-J2----
J2190922o4 <- output4[c(1,47:49)]
names(J2190922o4)[1] <- "ASV"
names(J2190922o4)[2] <- "PCR.1"
names(J2190922o4)[3] <- "PCR.2"
names(J2190922o4)[4] <- "PCR.3"

#determine if the seq is detected in all 3 reps
J2190922o4$Consistently_detected <- ifelse(J2190922o4$PCR.1==0 | J2190922o4$PCR.2== 0 | J2190922o4$PCR.3==0, "false", "true")

#make new df with only the seqs detected in all 3
J2190922o4_true <- subset(J2190922o4, Consistently_detected == "true",
                          select=c("ASV","PCR.1","PCR.2","PCR.3"))
J2190922o4_true$Fish <- "J2190922"
J2190922o4_true$Output <- "output4"


#compare size of two dfs (# ASVs and total reads)
#all ASVs
J2190922o4sum <- as.data.frame(colSums(J2190922o4[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
J2190922o4sum <- as.data.frame(t(J2190922o4sum))

#add column of avg depth across all PCR reps
J2190922o4sum$Avg.depth <-apply(J2190922o4sum,1,mean)
#add columns of total depth
J2190922o4sum$Total.depth <-apply(J2190922o4sum[,1:3],1,sum)

#add column of number of ASVs in the sample (each fish doesnt have reads in all ASVs)
J2190922o4sum$Num.ASVs <- sum(J2190922o4$PCR.1>0 | J2190922o4$PCR.2>0 | J2190922o4$PCR.3>0)

#add column of output folder (so can plot with others for comparison of parameters)
J2190922o4sum$Output <- "output4"

#add column of if in all 3 reps
J2190922o4sum$All.reps <- "No"

#only true ASVs
J2190922o4_truesum <- as.data.frame(colSums(J2190922o4_true[,c("PCR.1", "PCR.2", "PCR.3")]))

#transpose
J2190922o4_truesum <- as.data.frame(t(J2190922o4_truesum))

#add column of avg depth across all PCR reps
J2190922o4_truesum $Avg.depth <-apply(J2190922o4_truesum,1,mean)

#add columns of total depth
J2190922o4_truesum $Total.depth <-apply(J2190922o4_truesum[,1:3],1,sum)

#add column of number of ASVs in the sample
J2190922o4_truesum$Num.ASVs <- sum(J2190922o4_true$PCR.1>0 | J2190922o4_true$PCR.2>0 | J2190922o4_true$PCR.3>0)

#add column of pipeline (so can plot with UNIX for comparison)
J2190922o4_truesum$Output <- "output4"

#add column of if in all 3 reps
J2190922o4_truesum$All.reps <- "Yes"

#190922-J2 add data frames together
J2190922_output4 <- rbind(J2190922o4_truesum, J2190922o4sum)

#add column of fish ID
J2190922_output4$Fish <- "190922-J2"






#put it all all together ----
all_output4 <- rbind(A1190222_output4,A2190222_output4,H1190222_output4,M1190222_output4,
                     A1190918_output4,A3190918_output4,A4190918_output4,C4190918_output4,
                     D4190918_output4,E7190918_output4,F2190918_output4,D2190922_output4,
                     F2190922_output4,G2190922_output4,H2190922_output4,J2190922_output4)

#reorganizae
rownames(all_output4)<-c(1:32)
library(tidyverse)
all_output4 <- all_output4 %>% relocate(Fish, .before = PCR.1)

#save for future analysis
write_csv(all_output4, "UNIX-optimization/organized_ASVs/output4_ASVanalysis.csv")


#change data types
all_output4$Fish <- as.factor(all_output4$Fish)
all_output4$Output <- as.factor(all_output4$Output)
all_output4$All.reps <- as.factor(all_output4$All.reps)


#plot with bars for individual fish
#want to add a label to each bar showing how many ASVs this average was taken from
avgASV <- aggregate(Num.ASVs ~ Fish + All.reps, data = all_output4, FUN = mean)

avgdepth <- aggregate(Avg.depth ~ Fish + All.reps, data = all_output4, FUN = mean)


#plot with bars for individual fish
ggplot(avgdepth, aes(fill=All.reps, y=Avg.depth, x=Fish)) + 
  geom_bar(position="dodge", stat="identity") + 
  theme_linedraw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("Pipeline") + ylab("Average read depth") + 
  labs(fill= "ASV detected in all PCR replicates?") +
  geom_text(aes(label = round(avgASV$Num.ASVs, digits = 0)), position=position_dodge(width=0.9), vjust=-0.25, size=3)

#average across ALL fish
avgallASV <- aggregate(Num.ASVs ~ All.reps, data = all_output4, FUN = mean)

avgalldepth <- aggregate(Avg.depth ~ All.reps, data = all_output4, FUN = mean)

ggplot(avgalldepth, aes(fill=All.reps, y=Avg.depth, x=All.reps)) + 
  geom_bar(position="dodge", stat="identity") + 
  theme_linedraw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("ASV detected in all PCR replicates?") + ylab("Average read depth") + 
  labs(fill= "ASV detected in all PCR replicates?") +
  ggtitle("ASV analysis of output4") +
  guides(fill="none") +
  geom_text(aes(label = round(avgallASV$Num.ASVs, digits = 0)), position=position_dodge(width=0.9), vjust=-0.25, size=3)



#only looking at consistently detected ASVs----
alltrueo4 <- rbind(A1190222o4_true,A2190222o4_true,H1190222o4_true,M1190222o4_true,
                   A1190918o4_true,A3190918o4_true,A4190918o4_true,C4190918o4_true,
                   D4190918o4_true,E7190918o4_true,F2190918o4_true,D2190922o4_true,
                   F2190922o4_true,G2190922o4_true,H2190922o4_true,J2190922o4_true)

#save for future analysis
write_csv(alltrueo4, "UNIX-optimization/organized_ASVs/output4_TrueASVs.csv")

uniqueASVo4 <- unique(alltrueo4$ASV)
length(unique(alltrueo4$ASV))



#1:1574 ASVs seen in all FISH (but not all of these ASVs are present in EACH fish)
#2:1565
#3:1606
#4:1612
#5:1613
#6:1619
#7:1210
#8:1651
#9:1236
#10:2795

summary(alltrueo4)
alltrueo4$ASV <- as.factor(alltrueo4$ASV)
alltrueo4$Fish <- as.factor(alltrueo4$Fish)
alltrueo4$Output <- as.factor(alltrueo4$Output)

#add column of avg depth across all PCR reps
alltrueo4$Avg.depth <-apply(alltrueo4[,2:4],1,mean)
o4trueASVdepth <- aggregate(Avg.depth ~ Fish, data = alltrueo4, FUN = sum)
sum(o4trueASVdepth$Avg.depth)

#avg depth of unique ASVs per fish (avg. across PCR reps) summed across all individuals
#1:62805.67
#2:61970
#3:88145.67
#4:88413
#5:85689
#6:85983.67
#7:72895.33
#8:86018.67
#9:72923.67
#10:97960

#is this the most informative metric?