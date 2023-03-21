#Analysis of .tsv files to get reads/fish
#completed for outputs: 1,2,5,6,7,8,9,10

#last pushed to github: after out10, 23.03.14


# Install the required package
install.packages("readr")
# Load the installed Package
library(readr)

#be in "/Users/samanthabeal/Documents/MSc/Bioinformatics"

#output1----
setwd("../../..")
getwd()

#derep
derep <-readr::read_tsv("UNIX-optimization/output1/ASV/derep.tsv", show_col_types = FALSE)
derep2 <- derep[,-1]
rownames(derep2) <- derep[,1]
derep3 <- as.data.frame(colSums(derep2))

#denoise
denoised <-readr::read_tsv("UNIX-optimization/output1/ASV/denoised.tsv", show_col_types = FALSE)
denoised2 <- denoised[,-1]
rownames(denoised2) <- denoised[,1]
denoised3 <- as.data.frame(colSums(denoised2))

#nochim
nochim <-readr::read_tsv("UNIX-optimization/output1/ASV/nochim.tsv", show_col_types = FALSE)
nochim2 <- nochim[,-1]
rownames(nochim2) <- nochim[,1]
nochim3 <- as.data.frame(colSums(nochim2))





#output2----
setwd("../../..")
getwd()
setwd("UNIX-optimization/output2/ASV")

#derep
#derep <-readr::read_tsv("derep.tsv", show_col_types = FALSE)
#derep2 <- derep[,-1]
#rownames(derep2) <- derep[,1]
#derep3 <- as.data.frame(rowSums(derep2))
#colSums(derep3)

#denoise
denoised <-readr::read_tsv("denoised.tsv", show_col_types = FALSE)
denoised2 <- denoised[,-1]
rownames(denoised2) <- denoised[,1]
denoised3 <- as.data.frame(rowSums(denoised2))
colSums(denoised3)

#nochim
nochim <-readr::read_tsv("nochim.tsv", show_col_types = FALSE)
nochim2 <- nochim[,-1]
rownames(nochim2) <- nochim[,1]
nochim3 <- as.data.frame(rowSums(nochim2))
colSums(nochim3)


#output5 (formerly cut2)----
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


#output9----
setwd("../../..")
getwd()

#derep
derep <-readr::read_tsv("UNIX-optimization/output9/ASV/derep.tsv", show_col_types = FALSE)
derep2 <- derep[,-1]
rownames(derep2) <- derep[,1]
derep3 <- as.data.frame(colSums(derep2))

#denoise
denoised <-readr::read_tsv("UNIX-optimization/output9/ASV/denoised.tsv", show_col_types = FALSE)
denoised2 <- denoised[,-1]
rownames(denoised2) <- denoised[,1]
denoised3 <- as.data.frame(colSums(denoised2))


#nochim
nochim <-readr::read_tsv("UNIX-optimization/output9/ASV/nochim.tsv", show_col_types = FALSE)
nochim2 <- nochim[,-1]
rownames(nochim2) <- nochim[,1]
nochim3 <- as.data.frame(colSums(nochim2))

#output10----

setwd("../../..")
getwd()

#derep
derep <-readr::read_tsv("UNIX-optimization/output10/ASV/derep.tsv", show_col_types = FALSE)
derep2 <- derep[,-1]
rownames(derep2) <- derep[,1]
derep3 <- as.data.frame(colSums(derep2))

#denoise
denoised <-readr::read_tsv("UNIX-optimization/output10/ASV/denoised.tsv", show_col_types = FALSE)
denoised2 <- denoised[,-1]
rownames(denoised2) <- denoised[,1]
denoised3 <- as.data.frame(colSums(denoised2))

#nochim
nochim <-readr::read_tsv("UNIX-optimization/output10/ASV/nochim.tsv", show_col_types = FALSE)
nochim2 <- nochim[,-1]
rownames(nochim2) <- nochim[,1]
nochim3 <- as.data.frame(colSums(nochim2))

#output11----

setwd("../../..")
getwd()

#derep
derep <-readr::read_tsv("UNIX-optimization/output11/ASV/derep.tsv", show_col_types = FALSE)
derep <- as.data.frame(derep)
derep2 <- derep[,-1]
rownames(derep2) <- derep[,1]
derep3 <- as.data.frame(colSums(derep2))

#denoise
denoised <-readr::read_tsv("UNIX-optimization/output11/ASV/denoised.tsv", show_col_types = FALSE)
denoised <- as.data.frame(denoised)
denoised2 <- denoised[,-1]
rownames(denoised2) <- denoised[,1]
denoised3 <- as.data.frame(colSums(denoised2))

#nochim
nochim <-readr::read_tsv("UNIX-optimization/output11/ASV/nochim.tsv", show_col_types = FALSE)
nochim <- as.data.frame(nochim)
nochim2 <- nochim[,-1]
rownames(nochim2) <- nochim[,1]
nochim3 <- as.data.frame(colSums(nochim2))

#output12----
setwd("../../..")
getwd()

#derep
derep <-readr::read_tsv("UNIX-optimization/output12/ASV/derep.tsv", show_col_types = FALSE)
derep <- as.data.frame(derep)
derep2 <- derep[,-1]
rownames(derep2) <- derep[,1]
derep3 <- as.data.frame(colSums(derep2))

#denoise
denoised <-readr::read_tsv("UNIX-optimization/output12/ASV/denoised.tsv", show_col_types = FALSE)
denoised <- as.data.frame(denoised)
denoised2 <- denoised[,-1]
rownames(denoised2) <- denoised[,1]
denoised3 <- as.data.frame(colSums(denoised2))

#nochim
nochim <-readr::read_tsv("UNIX-optimization/output12/ASV/nochim.tsv", show_col_types = FALSE)
nochim <- as.data.frame(nochim)
nochim2 <- nochim[,-1]
rownames(nochim2) <- nochim[,1]
nochim3 <- as.data.frame(colSums(nochim2))

#output13----
setwd("..")
getwd()

#derep
derep <-readr::read_tsv("UNIX-optimization/output13/ASV/derep.tsv", show_col_types = FALSE)
derep <- as.data.frame(derep)
derep2 <- derep[,-1]
rownames(derep2) <- derep[,1]
derep3 <- as.data.frame(colSums(derep2))

#denoise
denoised <-readr::read_tsv("UNIX-optimization/output13/ASV/denoised.tsv", show_col_types = FALSE)
denoised <- as.data.frame(denoised)
denoised2 <- denoised[,-1]
rownames(denoised2) <- denoised[,1]
denoised3 <- as.data.frame(colSums(denoised2))

#nochim
nochim <-readr::read_tsv("UNIX-optimization/output13/ASV/nochim.tsv", show_col_types = FALSE)
nochim <- as.data.frame(nochim)
nochim2 <- nochim[,-1]
rownames(nochim2) <- nochim[,1]
nochim3 <- as.data.frame(colSums(nochim2))
