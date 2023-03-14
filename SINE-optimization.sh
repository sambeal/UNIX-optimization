#!/bin/bash

cd /Users/samanthabeal/Documents/MSc/Bioinformatics/UNIX-optimization

#this directory contains:
#input folder of gDNA from 16 AWF (B2 fish was removed due to failure of 1 PCR rep)
#individual output folders for each set of parameter trials

#output1 = all default
#output2 = cutadapt -e0.15
#output3 = cutadapt -e0.5
#output4 = " and fastq_filter maxEE 2
#cut1 = cutadapt -e0.3
#output5 = cutadapt -e0.4
#putput6 = " and fastq_filter maxEE 2
#output7 = " and unoise alpha 1 (resulted in less reads that output 5)
#output8 = cutadapt -e0.4, fastq_filter maxEE 2, and minsize 4 (8K more reads than output 5)
#output9 = " and alpha 1 (60K reads LESS than output5, 7)
#output10 = cutadapt -e0.4, fastq_filter maxEE 2, minsize 4, alpha 3

################ Input data ################

# list files in folder
# extract part of name before second underscore and find unique
ls input | cut -d_ -f1,2 | sort | uniq

# store list of names in bash variable and check contents
samples=$(ls input | cut -d_ -f1,2 | sort | uniq)
echo $samples


# count number of sequences across all files in folder
cd input
gzcat *.fastq.gz | grep -c "^@M00" 
#2,197,228

# count number of sequences per file and output as list
# must be in the folder where you want to count sequences
#initial input files
#gunzip = g unzip (opposit of this is gzip)

for s in $samples;
do
echo "${s}_L001_R1_001.fastq.gz" \
count=$(gunzip -c ${s}_L001_R1_001.fastq.gz | grep -c "^@M00" ) \
>> SeqNumber_input.txt;
done

# ~20K-70K reads/sample



################ Primer removal ################

cd ../.. #(back to UNIX-optimization)

#do not anchor (^) primers - anchoring tells cutadapt that the primer is at 
#the beginning of ther read
#set minimum overlap parameter (-O) to 18 
#(i.e. must find, in order at least 18 bases of primer seq)
#add 'X' to indicate nothing (i.e. sequence should start at first bp of primer)

#remove primers
# SINE (Smal cor-II) - only R1 sequenced BUT BOTH PRIMERS ARE ON IT DUH
#need to remove the reverse compliment "AAAAGCGTCTGCTAAATGGCA" - two step process
#forward = 20bp
#revcomp = 21 bp

#remove reverse first, put into new directory, move into that directory, remove forward
mkdir output9
mkdir revcompremoved9

for s in $samples;
do
cutadapt -a "AAAAGCGTCTGCTAAATGGCA;e=0.4" \
-o revcompremoved9/${s}_L001_R1_001.fastq.gz --discard-untrimmed \
input/${s}_L001_R1_001.fastq.gz;
done

#now remove forward primer
#expect SINE to be at the beginning of the read but won't anchor it just in case

#look into what cutadapt is doing by default here
#some errors are being allowed, since these sequences are longer 
#and so the default 10% allows 2 errors in these primers

#changing parameters in cutadapt:
#You can do so by adding a semicolon and parameter=value to the end of the adapter sequence
#https://cutadapt.readthedocs.io/en/stable/guide.html#search-parameters

for s in $samples;
do
cutadapt -g "TAGCTCAGCTGGTAGAGCAC;e=0.4" \
-o output9/${s}_L001_R1_001.fastq.gz --discard-untrimmed \
revcompremoved9/${s}_L001_R1_001.fastq.gz;
done

# count number of sequences across all files in folder
cd output9
gzcat *.fastq.gz | grep -c "^@M00" 

for s in $samples;
do
            echo "${s}_L001_R1_001.fastq.gz" \
            count=$(gunzip -c ${s}_L001_R1_001.fastq.gz | grep -c "^@M00" ) \
       >> SeqNumber_primerremoval.txt;
done


################ Sequence quality ################
# make sure in folder 'output'
mkdir qc

# Check sequence quality
# Run fastqc on each file in input
ls *.fastq.gz | parallel 'fastqc {}'

# Move qc outputs
mv /Users/samanthabeal/Documents/MSc/Bioinformatics/UNIX-optimization/output9/*.html /Users/samanthabeal/Documents/MSc/Bioinformatics/UNIX-optimization/output9/qc
mv /Users/samanthabeal/Documents/MSc/Bioinformatics/UNIX-optimization/output9/*.zip /Users/samanthabeal/Documents/MSc/Bioinformatics/UNIX-optimization/output9/qc

cd qc
multiqc .



################ Concatenate Samples ################

cd .. #(back to output)
mkdir ASV #(in output)

gunzip *.fastq.gz
#turns fastq.gz to .fastq
#g unzip

# add sample name to sequences
# combine all sequences in one file
for f in *;
do                         
sed -e "s/\(^@M00.*\) .*$/\1;sample=${f%.*};/" $f \
>> ASV/concatenated.fastq;
done


#count seqs
cd ASV
grep -c "M00" concatenated.fastq


################ Quality Filtering ################

#parameters of concern:
#--fastq_maxee 1 = default

# make sure in output folder

vsearch --fastx_filter concatenated.fastq --fastq_maxee 2 --fastaout concatenated.fasta

#count seqs - this way is not informative as to the depth/fish, only total
grep -c ">M00" concatenated.fasta

#making a .tsv of the seqs takes WAY TOO LONG
#will just look at the derep .tsv as no other sequence manipulation happens between 
#making the concatenated.fasta (after filtering) and dereplicating the data (reorganizing it)

#grep will work for this since each unique seq is it's own line
#won't work downstream of here once the seqs are depreplicated 
#(all unique seqs = smooshed together with a size value)


################ Dereplication ################

# header of each unique sequence records number of copies in dataset

vsearch --derep_fulllength concatenated.fasta --output derep.fasta --sizeout --relabel uniq

#count seqs
vsearch --search_exact concatenated.fasta -db derep.fasta -otutabout derep.tsv


################ Denoising ################

#paramters of concern:
#--minsize 8 = default
#--unoise_alpha 2 = default
#lower unoise = less strict?
vsearch --cluster_unoise derep.fasta --minsize 4 --unoise_alpha 3 --centroids denoised.fasta

#count seqs
vsearch --search_exact concatenated.fasta -db denoised.fasta -otutabout denoised.tsv

################ Chimera Filtering ################

vsearch --uchime3_denovo denoised.fasta --nonchimeras nochim.fasta

#count seqs = following 2 steps of code


################Mapping Reads to ASVs################

vsearch --search_exact concatenated.fasta -db nochim.fasta -otutabout nochim.tsv


################Identifying ASVs with Blast################

# -outmt 6 is standard blast tabulated output
# set strict evalue and perc_identity since expect many closely related species
# -evalue 
# -perc_identity

cd ../..

#blast final file - after chimeras are removed - unsure if i need to do this since
#expect all sequences to be whitefish 

blastn -query output9/ASV/nochim.fasta -subject SINE.fasta -outfmt 6 -out output9/ASV/nochim.txt \
-num_threads 1 -evalue 0.001 -perc_identity 97

#lose most of my reads at this step


