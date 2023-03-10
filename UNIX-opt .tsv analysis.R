# Install the required package
install.packages("readr")
# Load the installed Package
library(readr)

#be in "/Users/samanthabeal/Documents/MSc/Bioinformatics"
#output1----
setwd("../../..")
getwd()

#derep
derep <-readr::read_tsv("UNIX-optimization/output1/ASV/derep1.tsv", show_col_types = FALSE)
derep2 <- derep[,-1]
rownames(derep2) <- derep[,1]
derep3 <- as.data.frame(colSums(derep2))

#denoise
denoised <-readr::read_tsv("UNIX-optimization/output1/ASV/denoised1.tsv", show_col_types = FALSE)
denoised2 <- denoised[,-1]
rownames(denoised2) <- denoised[,1]
denoised3 <- as.data.frame(colSums(denoised2))

#nochim
nochim <-readr::read_tsv("UNIX-optimization/output1/ASV/nochim1.tsv", show_col_types = FALSE)
nochim2 <- nochim[,-1]
rownames(nochim2) <- nochim[,1]
nochim3 <- as.data.frame(colSums(nochim2))





#output2----
setwd("../../..")
getwd()
setwd("UNIX-optimization/output2/ASV")

#derep
#derep <-readr::read_tsv("derep1.tsv", show_col_types = FALSE)
#derep2 <- derep[,-1]
#rownames(derep2) <- derep[,1]
#derep3 <- as.data.frame(rowSums(derep2))
#colSums(derep3)

#denoise
denoised <-readr::read_tsv("denoised2.tsv", show_col_types = FALSE)
denoised2 <- denoised[,-1]
rownames(denoised2) <- denoised[,1]
denoised3 <- as.data.frame(rowSums(denoised2))
colSums(denoised3)

#nochim
nochim <-readr::read_tsv("nochim1.tsv", show_col_types = FALSE)
nochim2 <- nochim[,-1]
rownames(nochim2) <- nochim[,1]
nochim3 <- as.data.frame(rowSums(nochim2))
colSums(nochim3)


#cut2----
setwd("../../..")
getwd()


#derep
derep <-readr::read_tsv("UNIX-optimization/cut2/ASV/derep.tsv", show_col_types = FALSE)
derep2 <- derep[,-1]
rownames(derep2) <- derep[,1]
derep3 <- as.data.frame(colSums(derep2))

#denoise
denoised <-readr::read_tsv("UNIX-optimization/cut2/ASV/denoised.tsv", show_col_types = FALSE)
denoised2 <- denoised[,-1]
rownames(denoised2) <- denoised[,1]
denoised3 <- as.data.frame(colSums(denoised2))

#nochim
nochim <-readr::read_tsv("UNIX-optimization/cut2/ASV/nochim.tsv", show_col_types = FALSE)
nochim2 <- nochim[,-1]
rownames(nochim2) <- nochim[,1]
nochim3 <- as.data.frame(colSums(nochim2))


#output5----
setwd("../../..")
getwd()

#derep
derep <-readr::read_tsv("UNIX-optimization/output5/ASV/derep.tsv", show_col_types = FALSE)
derep2 <- derep[,-1]
rownames(derep2) <- derep[,1]
derep3 <- as.data.frame(colSums(derep2))

#denoise
denoised <-readr::read_tsv("UNIX-optimization/output5/ASV/denoised.tsv", show_col_types = FALSE)
denoised2 <- denoised[,-1]
rownames(denoised2) <- denoised[,1]
denoised3 <- as.data.frame(colSums(denoised2))

#nochim
nochim <-readr::read_tsv("UNIX-optimization/output5/ASV/nochim.tsv", show_col_types = FALSE)
nochim2 <- nochim[,-1]
rownames(nochim2) <- nochim[,1]
nochim3 <- as.data.frame(colSums(nochim2))


#output6----
setwd("../../..")
getwd()

#derep
derep <-readr::read_tsv("UNIX-optimization/output6/ASV/derep.tsv", show_col_types = FALSE)
derep2 <- derep[,-1]
rownames(derep2) <- derep[,1]
derep3 <- as.data.frame(colSums(derep2))

#denoise
denoised <-readr::read_tsv("UNIX-optimization/output6/ASV/denoised.tsv", show_col_types = FALSE)
denoised2 <- denoised[,-1]
rownames(denoised2) <- denoised[,1]
denoised3 <- as.data.frame(colSums(denoised2))

#nochim
nochim <-readr::read_tsv("UNIX-optimization/output6/ASV/nochim.tsv", show_col_types = FALSE)
nochim2 <- nochim[,-1]
rownames(nochim2) <- nochim[,1]
nochim3 <- as.data.frame(colSums(nochim2))


#output7----
setwd("../../..")
getwd()

#derep
derep <-readr::read_tsv("UNIX-optimization/output7/ASV/derep.tsv", show_col_types = FALSE)
derep2 <- derep[,-1]
rownames(derep2) <- derep[,1]
derep3 <- as.data.frame(colSums(derep2))

#denoise
denoised <-readr::read_tsv("UNIX-optimization/output7/ASV/denoised.tsv", show_col_types = FALSE)
denoised2 <- denoised[,-1]
rownames(denoised2) <- denoised[,1]
denoised3 <- as.data.frame(colSums(denoised2))

#nochim
nochim <-readr::read_tsv("UNIX-optimization/output7/ASV/nochim.tsv", show_col_types = FALSE)
nochim2 <- nochim[,-1]
rownames(nochim2) <- nochim[,1]
nochim3 <- as.data.frame(colSums(nochim2))


#output8----
setwd("../../..")
getwd()

#derep
derep <-readr::read_tsv("UNIX-optimization/output8/ASV/derep.tsv", show_col_types = FALSE)
derep2 <- derep[,-1]
rownames(derep2) <- derep[,1]
derep3 <- as.data.frame(colSums(derep2))

#denoise
denoised <-readr::read_tsv("UNIX-optimization/output8/ASV/denoised.tsv", show_col_types = FALSE)
denoised2 <- denoised[,-1]
rownames(denoised2) <- denoised[,1]
denoised3 <- as.data.frame(colSums(denoised2))


#nochim
nochim <-readr::read_tsv("UNIX-optimization/output8/ASV/nochim.tsv", show_col_types = FALSE)
nochim2 <- nochim[,-1]
rownames(nochim2) <- nochim[,1]
nochim3 <- as.data.frame(colSums(nochim2))
