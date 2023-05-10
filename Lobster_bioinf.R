# Make sure you have access to a unix system
# install conda (or miniconda etc)
# create an environment in conda, call it a name and install:

# cutadapt
conda install -c bioconda cutadapt
# sabre
conda install -c bioconda sabre
# dos2unix
conda install -c conda-forge dos2unix
# R
conda install -c r r

# Open the terminal and navigate to your sequencing raw_data folder.
cd /Users/siljeleiknes/Masteroppgave/seqdata/X204SC22122693-Z01-F001/01.RawData

## I made 5 folders - one for each plate in the same directory as all the sequence folders. I made a text file list of each sample in each plate, 
## from your excel sheet, then made sure is is in unix format (dos2unix below  use conda to install) then moved each sample that was in each plate into the appropriate folder
# make directories
mkdir leu_p1
mkdir leu_p2
mkdir leu_p3
mkdir leu_p4
mkdir leu_p5
# Make sure text lists are right format
conda install -c conda-forge dos2unix
dos2unix LEu_P1.txt
dos2unix LEu_P2.txt
dos2unix LEu_p3.txt
dos2unix LEu_P4.txt
dos2unix LEu_P5.txt
## moved sample sequence folders into appropriate plate
for d in $(cat LEu_P1.txt); do mv $d leu_p1/ ; done
for d in $(cat LEu_P2.txt); do mv $d leu_p2/ ; done
for d in $(cat LEu_p3.txt); do mv $d leu_p3/ ; done
for d in $(cat LEu_P4.txt); do mv $d leu_p4/ ; done
for d in $(cat LEu_P5.txt); do mv $d leu_p5/ ; done
## went into each plate folder and extracted all .fq.gz files into the main directory (i.e. plate1, plate2...)

## go into plate1 directory
cd leu_p1
## extract all sequence files from sub-directories into main directory
mv **/*.fq.gz /Users/siljeleiknes/Masteroppgave/seqdata/X204SC22122693-Z01-F001/01.RawData/leu_p1
## exit that directory 
cd ..
## etc
cd leu_p2
mv **/*.fq.gz /Users/siljeleiknes/Masteroppgave/seqdata/X204SC22122693-Z01-F001/01.RawData/leu_p2
cd ..
cd leu_p3
mv **/*.fq.gz /Users/siljeleiknes/Masteroppgave/seqdata/X204SC22122693-Z01-F001/01.RawData/leu_p3
cd ..
cd leu_p4
mv **/*.fq.gz /Users/siljeleiknes/Masteroppgave/seqdata/X204SC22122693-Z01-F001/01.RawData/leu_p4
cd ..
cd leu_p5
mv **/*.fq.gz /Users/siljeleiknes/Masteroppgave/seqdata/X204SC22122693-Z01-F001/01.RawData/leu_p5

#### RUN THROUGH OF PLATE 4 ONLY
## now splitting replicates in each plate - a shell script using sabre
## go into plate4 directory
cd leu_p4

## based on the forward barcode read for plate4 that is designated in another text file in the same directory (plate4_barcode_data.txt)
## make shell script to run sabre 
echo 'for i in *1.fq.gz; do bn=${i/1.fq.gz};
sabre pe -f ${bn}1.fq.gz -r ${bn}2.fq.gz -b /Users/siljeleiknes/Masteroppgave/seqdata/X204SC22122693-Z01-F001/01.RawData/leu_p4/repbarcodes.txt -u ${bn}unassigned1.fq -w ${bn}unassigned.fq;
mv rep1f ${bn}_rep1f.fq;
mv rep1r ${bn}_rep1r.fq;
mv rep2f ${bn}_rep2f.fq;
mv rep2r ${bn}_rep2r.fq;
mv rep3f ${bn}_rep3f.fq;
mv rep3r ${bn}_rep3r.fq;
mv rep4f ${bn}_rep4f.fq;
mv rep4r ${bn}_rep4r.fq;  done' > sabreplate4.sh

# run shell script
bash sabreplate4.sh

## Using cut adapt to remove primers - Doing replicate by replicate and make reverse primer an exact match to the corresponding reverse barcode that was not removed by sabre

#### Trnl base primers (i.e. tag IDs etc) are:
# Euka Fwd: TGGTGCATGGCCGTTCTTAGT
# Euka Rv: CATCTAAGGGCATCACAGACC

#Plate 4 rep 1 reverse barcode is CCACGTC
# Trnlg (fwd): TGGTGCATGGCCGTTCTTAGT
# Trnlh (rcRv): GGTCTGTGATGCCCTTAGATGGACGTGG
# Trnlg (rcfwd): ACTAAGAACGGCCATGCACCA
# Trnlh (Rv): CCACGTCCATCTAAGGGCATCACAGACC

echo 'for i in *rep1f.fq; do bn=${i/rep1f.fq}; cutadapt -a TGGTGCATGGCCGTTCTTAGT...GGTCTGTGATGCCCTTAGATGGACGTGG -A CCACGTCCATCTAAGGGCATCACAGACC...ACTAAGAACGGCCATGCACCA --untrimmed-output ${bn}.rep1out1.fq.gz --untrimmed-paired-output ${bn}.rep1out2.fq.gz -o ${bn}.rep1.trim1.fq.gz -p ${bn}.rep1.trim2.fq.gz ${bn}rep1f.fq ${bn}rep1r.fq; done' > plate4rep1.sh
bash plate4rep1.sh

#Plate 4 rep 2 reverse barcode is TTCTCAGC
# Trnlg (fwd): TGGTGCATGGCCGTTCTTAGT
# Trnlh (rcRv): GGTCTGTGATGCCCTTAGATGGCTGAGAA
# Trnlg (rcfwd): ACTAAGAACGGCCATGCACCA
# Trnlh (Rv): TTCTCAGCCATCTAAGGGCATCACAGACC

echo 'for i in *rep2f.fq; do bn=${i/rep2f.fq}; cutadapt -a TGGTGCATGGCCGTTCTTAGT...GGTCTGTGATGCCCTTAGATGGCTGAGAA -A TTCTCAGCCATCTAAGGGCATCACAGACC...ACTAAGAACGGCCATGCACCA --untrimmed-output ${bn}.rep2out1.fq.gz --untrimmed-paired-output ${bn}.rep2out2.fq.gz -o ${bn}.rep2.trim1.fq.gz -p ${bn}.rep2.trim2.fq.gz ${bn}rep2f.fq ${bn}rep2r.fq; done' > plate4rep2.sh
bash plate4rep2.sh

#Plate 4 rep 3 reverse barcode is CTAGG
# Trnlg (fwd): TGGTGCATGGCCGTTCTTAGT
# Trnlh (rcRv): GGTCTGTGATGCCCTTAGATGCCTAG
# Trnlg (rcfwd): ACTAAGAACGGCCATGCACCA
# Trnlh (Rv): CTAGGCATCTAAGGGCATCACAGACC

echo 'for i in *rep3f.fq; do bn=${i/rep3f.fq}; cutadapt -a TGGTGCATGGCCGTTCTTAGT...GGTCTGTGATGCCCTTAGATGCCTAG -A CTAGGCATCTAAGGGCATCACAGACC...ACTAAGAACGGCCATGCACCA --untrimmed-output ${bn}.rep3out1.fq.gz --untrimmed-paired-output ${bn}.rep3out2.fq.gz -o ${bn}.rep3.trim1.fq.gz -p ${bn}.rep3.trim2.fq.gz ${bn}rep3f.fq ${bn}rep3r.fq; done' > plate4rep3.sh
bash plate4rep3.sh

#Plate 4 rep 4 reverse barcode is TGCTTA
# Trnlg (fwd): TGGTGCATGGCCGTTCTTAGT
# Trnlh (rcRv): GGTCTGTGATGCCCTTAGATGTAAGCA
# Trnlg (rcfwd): ACTAAGAACGGCCATGCACCA
# Trnlh (Rv): TGCTTACATCTAAGGGCATCACAGACC

echo 'for i in *rep4f.fq; do bn=${i/rep4f.fq}; cutadapt -a TGGTGCATGGCCGTTCTTAGT...GGTCTGTGATGCCCTTAGATGTAAGCA -A TGCTTACATCTAAGGGCATCACAGACC...ACTAAGAACGGCCATGCACCA --untrimmed-output ${bn}.rep4out1.fq.gz --untrimmed-paired-output ${bn}.rep4out2.fq.gz -o ${bn}.rep4.trim1.fq.gz -p ${bn}.rep4.trim2.fq.gz ${bn}rep4f.fq ${bn}rep4r.fq; done' > plate4rep4.sh
bash plate4rep4.sh

## moving trimmed files to their own directory for dada2 analysis
mkdir leu_p4trimmed
mv *trim1.fq.gz /Users/siljeleiknes/Masteroppgave/seqdata/X204SC22122693-Z01-F001/01.RawData/leu_p4/leu_p4trimmed
mv *trim2.fq.gz /Users/siljeleiknes/Masteroppgave/seqdata/X204SC22122693-Z01-F001/01.RawData/leu_p4/leu_p4trimmed

#### Moving into R for DADA2 processing... you can either run R from the terminal or just go to Rstudio and point 
#### to the folder where your trimmed sequence files are
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2")
library(dada2)
packageVersion("dada2")

##### PLate 1 in DADA2 ######################################################################################################################
## Making filepath based on where the trimmed files are. Now is a good point to follow the DADA2 tutorial while trying these steps
setwd("/Users/siljeleiknes/Masteroppgave/seqdata/X204SC22122693-Z01-F001/01.RawData/leu_p1/leu_p1trimmed")
path <- "/Users/siljeleiknes/Masteroppgave/seqdata/X204SC22122693-Z01-F001/01.RawData/leu_p1/leu_p1trimmed"
fnFs <- sort(list.files(path, pattern=".trim1.fq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern=".trim2.fq.gz", full.names = TRUE))
###tidying up sample names and replicates
sample.names <-sapply(strsplit(basename(fnFs), "_"), function(x){paste(x[[1]], x[[2]], sep="_")})
reps <- rep(c("r1", "r2", "r3", "r4"))

sample.names <-  paste0(sample.names, reps)

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

## checking some quality plots 
plotQualityProfile(fnFs[10:20]) # 
plotQualityProfile(fnRs[10:20])  # 

## Filtering with standard parameters (did not truncate reads as they are already very short - no quality issues on the tails)
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, # truncLen=c(45,45),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE

head(out, n=50)
## Learning error rates on these data
errF <- learnErrors(filtFs, multithread=TRUE, MAX_CONSIST=20)
errR <- learnErrors(filtRs, multithread=TRUE, MAX_CONSIST=20)

## Having a look at error plots
oldFerrs <- plotErrors(errF, nominalQ=TRUE)
oldFerrs
#dev.off() 
#pdf("rv_error.pdf")
oldRerrs <- plotErrors(errR, nominalQ=TRUE)
oldRerrs

### applying dada2 core inference algorithm - using default of all libraries processed seperately - no pooling
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

###  also applying core inference algorith - but using pooling (all samples pooled together for sample inference
dadaFPPs <- dada(filtFs, err=errF, multithread=TRUE, pool=TRUE)
dadaRPPs <- dada(filtRs, err=errR, multithread=TRUE, pool=TRUE)

###  also applying core inference algorith - but using pseudo pooling (samples processed independent after sharing info
### between samples)
dadaFpsPPs <- dada(filtFs, err=errF, multithread=TRUE, pool="pseudo")
dadaRpsPPs <- dada(filtRs, err=errR, multithread=TRUE, pool="pseudo")

## Merging the paired-ends - no pooling
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

## Merging the paired-ends - true pooling
mergersPP <- mergePairs(dadaFPPs, filtFs, dadaRPPs, filtRs, verbose=TRUE)

### Merging the paired-ends - pseudo-pooling
mergers_psPP <- mergePairs(dadaFpsPPs, filtFs, dadaRpsPPs, filtRs, verbose=TRUE)

## Making sequence table
leu_p1_seqtab <- makeSequenceTable(mergers)
leu_p1_seqtabPP <- makeSequenceTable(mergersPP)
leu_p1_seqtabpsPP <- makeSequenceTable(mergers_psPP)

## keeping track of pipeline loss up to merging tables
#### Tracking read loss through the pipeline - non pooled
getN <- function(x) sum(getUniques(x))
leu_p1_track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN))
colnames(leu_p1_track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
head(leu_p1_track)

#### Tracking read loss through the pipeline - total pooled
leu_p1_trackPP <- cbind(out, sapply(dadaFPPs, getN), sapply(dadaRPPs, getN), sapply(mergersPP, getN))
colnames(leu_p1_trackPP) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
head(leu_p1_trackPP)

#### Tracking read loss through the pipeline - psuedo pooled
leu_p1_trackpsPP <- cbind(out, sapply(dadaFpsPPs, getN), sapply(dadaRpsPPs, getN), sapply(mergers_psPP, getN))
colnames(leu_p1_trackpsPP) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
head(leu_p1_trackpsPP)

##### PLate 2 in DADA2 ######################################################################################################################
## Making filepath based on where the trimmed files are. Now is a good point to follow the DADA2 tutorial while trying these steps
setwd("/Users/siljeleiknes/Masteroppgave/seqdata/X204SC22122693-Z01-F001/01.RawData/leu_p2/leu_p2trimmed")
path <- "/Users/siljeleiknes/Masteroppgave/seqdata/X204SC22122693-Z01-F001/01.RawData/leu_p2/leu_p2trimmed"
fnFs <- sort(list.files(path, pattern=".trim1.fq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern=".trim2.fq.gz", full.names = TRUE))
###tidying up sample names and replicates
sample.names <-sapply(strsplit(basename(fnFs), "_"), function(x){paste(x[[1]], x[[2]], sep="_")})
reps <- rep(c("r1", "r2", "r3", "r4"))

sample.names <-  paste0(sample.names, reps)

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

## checking some quality plots 
plotQualityProfile(fnFs[10:20]) # 
plotQualityProfile(fnRs[10:20])  # 

## Filtering with standard parameters (did not truncate reads as they are already very short - no quality issues on the tails)
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, # truncLen=c(45,45),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE

head(out, n=50)
## Learning error rates on these data
errF <- learnErrors(filtFs, multithread=TRUE, MAX_CONSIST=20)
errR <- learnErrors(filtRs, multithread=TRUE, MAX_CONSIST=20)

## error - some files did not pass trim and need to be removed from lists - which files??
df.fe <-  data.frame(theref = file.exists(filtFs), therer = file.exists(filtRs), filef = filtRs, filer = filtRs)
subset(df.fe, theref == "FALSE") 
subset(df.fe, therer == "FALSE") 

##redefining 'out' matrix so these files are not included
out <- out[file.exists(filtFs),]

### dropping samples that are empty - these numbers will change plate by plate - BE CAREFUL THESE ARE FOR PLATE 1 ONLY
### to tell which numbers to add, loo at the output from the above subset commands.
filtFs <- filtFs[-c(346,359)]
filtRs <- filtRs[-c(346,359)]

errF <- learnErrors(filtFs, multithread=TRUE, MAX_CONSIST=20)
errR <- learnErrors(filtRs, multithread=TRUE, MAX_CONSIST=20)

## Having a look at error plots
oldFerrs <- plotErrors(errF, nominalQ=TRUE)
oldFerrs
#dev.off() 
#pdf("rv_error.pdf")
oldRerrs <- plotErrors(errR, nominalQ=TRUE)
oldRerrs

### applying dada2 core inference algorithm - using default of all libraries processed seperately - no pooling
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

###  also applying core inference algorith - but using pooling (all samples pooled together for sample inference
dadaFPPs <- dada(filtFs, err=errF, multithread=TRUE, pool=TRUE)
dadaRPPs <- dada(filtRs, err=errR, multithread=TRUE, pool=TRUE)

###  also applying core inference algorith - but using pseudo pooling (samples processed independent after sharing info
### between samples)
dadaFpsPPs <- dada(filtFs, err=errF, multithread=TRUE, pool="pseudo")
dadaRpsPPs <- dada(filtRs, err=errR, multithread=TRUE, pool="pseudo")

## Merging the paired-ends - no pooling
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

## Merging the paired-ends - true pooling
mergersPP <- mergePairs(dadaFPPs, filtFs, dadaRPPs, filtRs, verbose=TRUE)

### Merging the paired-ends - pseudo-pooling
mergers_psPP <- mergePairs(dadaFpsPPs, filtFs, dadaRpsPPs, filtRs, verbose=TRUE)

## Making sequence table
leu_p2_seqtab <- makeSequenceTable(mergers)
leu_p2_seqtabPP <- makeSequenceTable(mergersPP)
leu_p2_seqtabpsPP <- makeSequenceTable(mergers_psPP)

## keeping track of pipeline loss up to merging tables
#### Tracking read loss through the pipeline - non pooled
getN <- function(x) sum(getUniques(x))
leu_p2_track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN))
colnames(leu_p2_track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
head(leu_p2_track)

#### Tracking read loss through the pipeline - total pooled
leu_p2_trackPP <- cbind(out, sapply(dadaFPPs, getN), sapply(dadaRPPs, getN), sapply(mergersPP, getN))
colnames(leu_p2_trackPP) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
head(leu_p2_trackPP)

#### Tracking read loss through the pipeline - psuedo pooled
leu_p2_trackpsPP <- cbind(out, sapply(dadaFpsPPs, getN), sapply(dadaRpsPPs, getN), sapply(mergers_psPP, getN))
colnames(leu_p2_trackpsPP) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
head(leu_p2_trackpsPP)

##### PLate 3 in DADA2 ######################################################################################################################
## Making filepath based on where the trimmed files are. Now is a good point to follow the DADA2 tutorial while trying these steps
setwd("/Users/siljeleiknes/Masteroppgave/seqdata/X204SC22122693-Z01-F001/01.RawData/leu_p3/leu_p3trimmed")
path <- "/Users/siljeleiknes/Masteroppgave/seqdata/X204SC22122693-Z01-F001/01.RawData/leu_p3/leu_p3trimmed"
fnFs <- sort(list.files(path, pattern=".trim1.fq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern=".trim2.fq.gz", full.names = TRUE))
###tidying up sample names and replicates
sample.names <-sapply(strsplit(basename(fnFs), "_"), function(x){paste(x[[1]], x[[2]], sep="_")})
reps <- rep(c("r1", "r2", "r3", "r4"))

sample.names <-  paste0(sample.names, reps)

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

## checking some quality plots 
plotQualityProfile(fnFs[10:20]) # 
plotQualityProfile(fnRs[10:20])  # 

## Filtering with standard parameters (did not truncate reads as they are already very short - no quality issues on the tails)
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, # truncLen=c(45,45),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE

head(out, n=50)
## Learning error rates on these data
errF <- learnErrors(filtFs, multithread=TRUE, MAX_CONSIST=20)
errR <- learnErrors(filtRs, multithread=TRUE, MAX_CONSIST=20)

## error - some files did not pass trim and need to be removed from lists - which files??
df.fe <-  data.frame(theref = file.exists(filtFs), therer = file.exists(filtRs), filef = filtRs, filer = filtRs)
subset(df.fe, theref == "FALSE") 
subset(df.fe, therer == "FALSE") 

##redefining 'out' matrix so these files are not included
out <- out[file.exists(filtFs),]

### dropping samples that are empty - these numbers will change plate by plate - BE CAREFUL THESE ARE FOR PLATE 1 ONLY
### to tell which numbers to add, loo at the output from the above subset commands.
filtFs <- filtFs[-c(356)]
filtRs <- filtRs[-c(356)]

errF <- learnErrors(filtFs, multithread=TRUE, MAX_CONSIST=20)
errR <- learnErrors(filtRs, multithread=TRUE, MAX_CONSIST=20)

## Having a look at error plots
oldFerrs <- plotErrors(errF, nominalQ=TRUE)
oldFerrs
#dev.off() 
#pdf("rv_error.pdf")
oldRerrs <- plotErrors(errR, nominalQ=TRUE)
oldRerrs

### applying dada2 core inference algorithm - using default of all libraries processed seperately - no pooling
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

###  also applying core inference algorith - but using pooling (all samples pooled together for sample inference
dadaFPPs <- dada(filtFs, err=errF, multithread=TRUE, pool=TRUE)
dadaRPPs <- dada(filtRs, err=errR, multithread=TRUE, pool=TRUE)

###  also applying core inference algorith - but using pseudo pooling (samples processed independent after sharing info
### between samples)
dadaFpsPPs <- dada(filtFs, err=errF, multithread=TRUE, pool="pseudo")
dadaRpsPPs <- dada(filtRs, err=errR, multithread=TRUE, pool="pseudo")

## Merging the paired-ends - no pooling
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

## Merging the paired-ends - true pooling
mergersPP <- mergePairs(dadaFPPs, filtFs, dadaRPPs, filtRs, verbose=TRUE)

### Merging the paired-ends - pseudo-pooling
mergers_psPP <- mergePairs(dadaFpsPPs, filtFs, dadaRpsPPs, filtRs, verbose=TRUE)

## Making sequence table
leu_p3_seqtab <- makeSequenceTable(mergers)
leu_p3_seqtabPP <- makeSequenceTable(mergersPP)
leu_p3_seqtabpsPP <- makeSequenceTable(mergers_psPP)

## keeping track of pipeline loss up to merging tables
#### Tracking read loss through the pipeline - non pooled
getN <- function(x) sum(getUniques(x))
leu_p3_track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN))
colnames(leu_p3_track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
head(leu_p3_track)

#### Tracking read loss through the pipeline - total pooled
leu_p3_trackPP <- cbind(out, sapply(dadaFPPs, getN), sapply(dadaRPPs, getN), sapply(mergersPP, getN))
colnames(leu_p3_trackPP) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
head(leu_p3_trackPP)

#### Tracking read loss through the pipeline - psuedo pooled
leu_p3_trackpsPP <- cbind(out, sapply(dadaFpsPPs, getN), sapply(dadaRpsPPs, getN), sapply(mergers_psPP, getN))
colnames(leu_p3_trackpsPP) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
head(leu_p3_trackpsPP)

##### PLate 4 in DADA2 ######################################################################################################################
## Making filepath based on where the trimmed files are. Now is a good point to follow the DADA2 tutorial while trying these steps
setwd("/Users/siljeleiknes/Masteroppgave/seqdata/X204SC22122693-Z01-F001/01.RawData/leu_p4/leu_p4trimmed")
path <- "/Users/siljeleiknes/Masteroppgave/seqdata/X204SC22122693-Z01-F001/01.RawData/leu_p4/leu_p4trimmed"
fnFs <- sort(list.files(path, pattern=".trim1.fq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern=".trim2.fq.gz", full.names = TRUE))
###tidying up sample names and replicates
sample.names <-sapply(strsplit(basename(fnFs), "_"), function(x){paste(x[[1]], x[[2]], sep="_")})
reps <- rep(c("r1", "r2", "r3", "r4"))

sample.names <-  paste0(sample.names, reps)

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

## checking some quality plots 
plotQualityProfile(fnFs[10:20]) # 
plotQualityProfile(fnRs[10:20])  # 

## Filtering with standard parameters (did not truncate reads as they are already very short - no quality issues on the tails)
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, # truncLen=c(45,45),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE

head(out, n=50)
## Learning error rates on these data
errF <- learnErrors(filtFs, multithread=TRUE, MAX_CONSIST=20)
errR <- learnErrors(filtRs, multithread=TRUE, MAX_CONSIST=20)

## error - some files did not pass trim and need to be removed from lists - which files??
df.fe <-  data.frame(theref = file.exists(filtFs), therer = file.exists(filtRs), filef = filtRs, filer = filtRs)
subset(df.fe, theref == "FALSE") 
subset(df.fe, therer == "FALSE") 

##redefining 'out' matrix so these files are not included
out <- out[file.exists(filtFs),]

### dropping samples that are empty - these numbers will change plate by plate - BE CAREFUL THESE ARE FOR PLATE 1 ONLY
### to tell which numbers to add, loo at the output from the above subset commands.

## Learning error rates on these data
errF <- learnErrors(filtFs, multithread=TRUE, MAX_CONSIST=20)
errR <- learnErrors(filtRs, multithread=TRUE, MAX_CONSIST=20)

## Having a look at error plots
oldFerrs <- plotErrors(errF, nominalQ=TRUE)
oldFerrs
#dev.off() 
#pdf("rv_error.pdf")
oldRerrs <- plotErrors(errR, nominalQ=TRUE)
oldRerrs

### applying dada2 core inference algorithm - using default of all libraries processed seperately - no pooling
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

###  also applying core inference algorith - but using pooling (all samples pooled together for sample inference

dadaFPPs <- dada(filtFs, err=errF, multithread=TRUE, pool=TRUE)
dadaRPPs <- dada(filtRs, err=errR, multithread=TRUE, pool=TRUE)

###  also applying core inference algorith - but using pseudo pooling (samples processed independent after sharing info
### between samples)

dadaFpsPPs <- dada(filtFs, err=errF, multithread=TRUE, pool="pseudo")
dadaRpsPPs <- dada(filtRs, err=errR, multithread=TRUE, pool="pseudo")

## Merging the paired-ends - no pooling
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

## Merging the paired-ends - true pooling
mergersPP <- mergePairs(dadaFPPs, filtFs, dadaRPPs, filtRs, verbose=TRUE)

### Merging the paired-ends - pseudo-pooling
mergers_psPP <- mergePairs(dadaFpsPPs, filtFs, dadaRpsPPs, filtRs, verbose=TRUE)

## Making sequence table
leu_p4_seqtab <- makeSequenceTable(mergers)
leu_p4_seqtabPP <- makeSequenceTable(mergersPP)
leu_p4_seqtabpsPP <- makeSequenceTable(mergers_psPP)

## keeping track of pipeline loss up to merging tables
#### Tracking read loss through the pipeline - non pooled
getN <- function(x) sum(getUniques(x))
leu_p4_track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN))
colnames(leu_p4_track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
head(leu_p4_track)
getwd()
write.csv(leu_p4_track, "plate4track.csv")
#### Tracking read loss through the pipeline - total pooled
leu_p4_trackPP <- cbind(out, sapply(dadaFPPs, getN), sapply(dadaRPPs, getN), sapply(mergersPP, getN))
colnames(leu_p4_trackPP) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
head(leu_p4_trackPP)

#### Tracking read loss through the pipeline - psuedo pooled
leu_p4_trackpsPP <- cbind(out, sapply(dadaFpsPPs, getN), sapply(dadaRpsPPs, getN), sapply(mergers_psPP, getN))
colnames(leu_p4_trackpsPP) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
head(leu_p4_trackpsPP)

##### PLate 5 in DADA2 ######################################################################################################################
## Making filepath based on where the trimmed files are. Now is a good point to follow the DADA2 tutorial while trying these steps
setwd("/Users/siljeleiknes/Masteroppgave/seqdata/X204SC22122693-Z01-F001/01.RawData/leu_p5/leu_p5trimmed")
path <- "/Users/siljeleiknes/Masteroppgave/seqdata/X204SC22122693-Z01-F001/01.RawData/leu_p5/leu_p5trimmed"
fnFs <- sort(list.files(path, pattern=".trim1.fq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern=".trim2.fq.gz", full.names = TRUE))
###tidying up sample names and replicates
sample.names <-sapply(strsplit(basename(fnFs), "_"), function(x){paste(x[[1]], x[[2]], sep="_")})
reps <- rep(c("r1", "r2", "r3", "r4"))

sample.names <-  paste0(sample.names, reps)

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

## checking some quality plots 
plotQualityProfile(fnFs[10:20]) # 
plotQualityProfile(fnRs[10:20])  # 

## Filtering with standard parameters (did not truncate reads as they are already very short - no quality issues on the tails)
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, # truncLen=c(45,45),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE

head(out, n=50)
## Learning error rates on these data
errF <- learnErrors(filtFs, multithread=TRUE, MAX_CONSIST=20)
errR <- learnErrors(filtRs, multithread=TRUE, MAX_CONSIST=20)

## error - some files did not pass trim and need to be removed from lists - which files??
df.fe <-  data.frame(theref = file.exists(filtFs), therer = file.exists(filtRs), filef = filtRs, filer = filtRs)
subset(df.fe, theref == "FALSE") 
subset(df.fe, therer == "FALSE") 

##redefining 'out' matrix so these files are not included
out <- out[file.exists(filtFs),]

### dropping samples that are empty - these numbers will change plate by plate - BE CAREFUL THESE ARE FOR PLATE 1 ONLY
### to tell which numbers to add, loo at the output from the above subset commands.
filtFs <- filtFs[-c(299)]
filtRs <- filtRs[-c(299)]

## Learning error rates on these data

errF <- learnErrors(filtFs, multithread=TRUE, MAX_CONSIST=20)
errR <- learnErrors(filtRs, multithread=TRUE, MAX_CONSIST=20)

## Having a look at error plots
oldFerrs <- plotErrors(errF, nominalQ=TRUE)
oldFerrs
#dev.off() 
#pdf("rv_error.pdf")
oldRerrs <- plotErrors(errR, nominalQ=TRUE)
oldRerrs

### applying dada2 core inference algorithm - using default of all libraries processed seperately - no pooling
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

###  also applying core inference algorith - but using pooling (all samples pooled together for sample inference
dadaFPPs <- dada(filtFs, err=errF, multithread=TRUE, pool=TRUE)
dadaRPPs <- dada(filtRs, err=errR, multithread=TRUE, pool=TRUE)

###  also applying core inference algorith - but using pseudo pooling (samples processed independent after sharing info
### between samples)
dadaFpsPPs <- dada(filtFs, err=errF, multithread=TRUE, pool="pseudo")
dadaRpsPPs <- dada(filtRs, err=errR, multithread=TRUE, pool="pseudo")

## Merging the paired-ends - no pooling
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

## Merging the paired-ends - true pooling
mergersPP <- mergePairs(dadaFPPs, filtFs, dadaRPPs, filtRs, verbose=TRUE)

### Merging the paired-ends - pseudo-pooling
mergers_psPP <- mergePairs(dadaFpsPPs, filtFs, dadaRpsPPs, filtRs, verbose=TRUE)

## Making sequence table
leu_p5_seqtab <- makeSequenceTable(mergers)
leu_p5_seqtabPP <- makeSequenceTable(mergersPP)
leu_p5_seqtabpsPP <- makeSequenceTable(mergers_psPP)

## keeping track of pipeline loss up to merging tables
#### Tracking read loss through the pipeline - non pooled
getN <- function(x) sum(getUniques(x))
leu_p5_track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN))
colnames(leu_p5_track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
head(leu_p5_track)

#### Tracking read loss through the pipeline - total pooled
leu_p5_trackPP <- cbind(out, sapply(dadaFPPs, getN), sapply(dadaRPPs, getN), sapply(mergersPP, getN))
colnames(leu_p5_trackPP) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
head(leu_p5_trackPP)

#### Tracking read loss through the pipeline - psuedo pooled
leu_p5_trackpsPP <- cbind(out, sapply(dadaFpsPPs, getN), sapply(dadaRpsPPs, getN), sapply(mergers_psPP, getN))
colnames(leu_p5_trackpsPP) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
head(leu_p5_trackpsPP)

library(stringr)
library(abind)
library(tidyverse)
library(devtools)
install_github("tobiasgf/lulu")
library(lulu)
library(textshape)
library(Hmisc)
library("myTAI")
library(dplyr)
library(data.table)
library(phyloseq)
library(DECIPHER)
library(phangorn)
library(metacoder)
library(vegan)
library(ggplot2)
library(gridExtra)
library("taxonomizr")
library("taxize")
library("taxizedb")
library(phylosmith)

## making unique names for controls - example
rownames(leu_p1_seqtab) <- str_replace_all(rownames(leu_p1_seqtab), "lake", "lakep1")

#### making lists of all the matrices for each run (no pooling, pooling, Pseudo-pooling)
npoollist <- list(leu_p1_seqtab, leu_p2_seqtab, leu_p3_seqtab, leu_p4_seqtab, leu_p5_seqtab)
poollist <- list(leu_p1_seqtabPP, leu_p2_seqtabPP, leu_p3_seqtabPP, leu_p4_seqtabPP, leu_p5_seqtabPP)
pspoollist <- list(leu_p1_seqtabpsPP, leu_p2_seqtabpsPP, leu_p3_seqtabpsPP, leu_p4_seqtabpsPP, leu_p5_seqtabpsPP)

## Tidying up names
name.change <- function(x){
  rownames(x) <-  gsub("_F_filt.fastq.gz", "", rownames(x), ignore.case =FALSE, perl = TRUE)
  rownames(x) <-  gsub("_EKDL220017784-1A", "", rownames(x), ignore.case =FALSE, perl = TRUE)
  return(x)
}

index.info <- function(x){
  y <- data.frame(matrix(NA,    # Create empty data frame
                         nrow = nrow(x),
                         ncol = 0))
  y$rep <- str_sub(rownames(x), start= -1)
  y$sample <- str_sub(rownames(x),1,nchar(rownames(x))-2)
  y$sample <- ifelse(startsWith(y$sample, "L2")==TRUE, substring(y$sample, 2), y$sample)
  y$full <- rownames(x)
  y$totseq <- rowSums(x)
  return(y)
}
##tidying up names and getting indexed summary stats - total sequences per replicate pcr per sample
npoolred <- lapply(npoollist, name.change)
names(npoolred) <- c("leu1_seqtab", "leu2_seqtab", "leu3_seqtab", "leu4_seqtab", "leu5_seqtab")
rownames(npoolred$leu1_seqtab)
npoolredind <- lapply(npoolred, index.info)
names(npoolredind) <- c("leu1_seqtab", "leu2_seqtab", "leu3_seqtab", "leu4_seqtab", "leu5_seqtab")

poolred <- lapply(poollist, name.change)
names(poolred) <- c("leu1_seqtab", "leu2_seqtab", "leu3_seqtab", "leu4_seqtab", "leu5_seqtab")
poolredind <- lapply(poolred, index.info)
names(poolredind) <- c("leu1_seqtabPP", "leu2_seqtabPP", "leu3_seqtabPP", "leu4_seqtabPP", "leu5_seqtabPP")

pspoolred <- lapply(poollist, name.change)
names(pspoolred) <- c("leu1_seqtab", "leu2_seqtab", "leu3_seqtab", "leu4_seqtab", "leu5_seqtab")
rownames(pspoolred$leu2_seqtab)
pspoolredind <- lapply(pspoolred, index.info)
names(pspoolredind) <- c("leu1_seqtabpsPP", "leu2_seqtabpsPP", "leu3_seqtabpsPP", "leu4_seqtabpsPP", "leu5_seqtabpsPP")

#### Now that all plates have been processed seperately, joining runs together for chimera removal
merge.npoolred <- mergeSequenceTables(npoolred$leu1_seqtab, npoolred$leu2_seqtab, npoolred$leu3_seqtab, npoolred$leu4_seqtab, npoolred$leu5_seqtab)
merge.poolred <- mergeSequenceTables(poolred$leu1_seqtab, poolred$leu2_seqtab, poolred$leu3_seqtab, poolred$leu4_seqtab, poolred$leu5_seqtab)
merge.pspoolred <- mergeSequenceTables(pspoolred$leu1_seqtab, pspoolred$leu2_seqtab, pspoolred$leu3_seqtab, pspoolred$leu4_seqtab, pspoolred$leu5_seqtab)

merge.npoolred.nochim <- removeBimeraDenovo(merge.npoolred, method="consensus", multithread=TRUE, verbose=TRUE)
dim(merge.npoolred.nochim)
sum(merge.npoolred.nochim)/sum(merge.npoolred)

merge.poolred.nochim <- removeBimeraDenovo(merge.poolred, method="consensus", multithread=TRUE, verbose=TRUE)
dim(merge.poolred.nochim)
sum(merge.poolred.nochim)/sum(merge.poolred)

merge.pspoolred.nochim <- removeBimeraDenovo(merge.pspoolred, method="consensus", multithread=TRUE, verbose=TRUE)
dim(merge.pspoolred.nochim)
sum(merge.pspoolred.nochim)/sum(merge.pspoolred)

str(merge.pspoolred.nochim)

## making single chimera table back into list of individual plates for NTC extraction
npoolred <- as.data.frame(merge.npoolred.nochim)
sample.info.nopool <- bind_rows(npoolredind, .id = "column_label")
rownames(sample.info.nopool) <- sample.info.nopool$full
npoolred <-  merge(npoolred, sample.info.nopool ,by =  'row.names', all.x=TRUE)
rownames(npoolred) <- npoolred$Row.names
npoolred <- split(npoolred, npoolred$column_label)
colnames(npoolred$leu1_seqtab)
dropnames <- colnames(npoolred$leu1_seqtab[, c(which(nchar(colnames(npoolred$leu1_seqtab))< 20))])
npoolred <- lapply(npoolred, function(x) x[!(names(x) %in% dropnames)])

poolred <- as.data.frame(merge.poolred.nochim)
sample.info.nopool <- bind_rows(poolredind, .id = "column_label")
rownames(sample.info.nopool) <- sample.info.nopool$full
poolred <-  merge(poolred, sample.info.nopool ,by =  'row.names', all.x=TRUE)
rownames(poolred) <- poolred$Row.names
poolred <- split(poolred, poolred$column_label)
colnames(poolred$leu1_seqtabPP)
dropnames <- colnames(poolred$leu1_seqtabPP[, c(which(nchar(colnames(poolred$leu1_seqtabPP))< 20))])
poolred <- lapply(poolred, function(x) x[!(names(x) %in% dropnames)])

pspoolred <- as.data.frame(merge.pspoolred.nochim)
sample.info.nopool <- bind_rows(pspoolredind, .id = "column_label")
rownames(sample.info.nopool) <- sample.info.nopool$full
pspoolred <-  merge(pspoolred, sample.info.nopool ,by =  'row.names', all.x=TRUE)
rownames(pspoolred) <- pspoolred$Row.names
pspoolred <- split(pspoolred, pspoolred$column_label)
colnames(pspoolred$leu1_seqtabpsPP)
dropnames <- colnames(pspoolred$leu1_seqtabpsPP[, c(which(nchar(colnames(pspoolred$leu1_seqtabpsPP))< 20))])
pspoolred <- lapply(pspoolred, function(x) x[!(names(x) %in% dropnames)])

## Moving to LULU curation
library(devtools)
install_github("tobiasgf/lulu")  
library(lulu)

str(merge.npoolred.nochim)
uniquesToFasta(merge.npoolred.nochim, fout="/Users/siljeleiknes/Masteroppgave/seqdata/X204SC22122693-Z01-F001/Nopool.nochim.fasta", ids=paste0("OTU",  seq(length(getSequences(merge.npoolred.nochim)))))
uniquesToFasta(merge.poolred.nochim, fout="/Users/siljeleiknes/Masteroppgave/seqdata/X204SC22122693-Z01-F001/Pooled.nochim.fasta", ids=paste0("OTU",seq(length(getSequences(merge.poolred.nochim)))))
uniquesToFasta(merge.pspoolred.nochim, fout="/Users/siljeleiknes/Masteroppgave/seqdata/X204SC22122693-Z01-F001/psPooled.nochim.fasta", ids=paste0("OTU", seq(length(getSequences(merge.pspoolred.nochim)))))

## Make LULU OTU tables (OTUs: rows, samples: columns)
npool.lulu <- merge.npoolred.nochim
colnames(npool.lulu) <- paste0("OTU", seq(length(getSequences(merge.npoolred.nochim))))
npool.lulu <- t(npool.lulu)

pool.lulu <- merge.poolred.nochim
colnames(pool.lulu) <- paste0("OTU", seq(length(getSequences(merge.poolred.nochim))))
pool.lulu <- t(pool.lulu)

pspool.lulu <- merge.pspoolred.nochim
colnames(pspool.lulu) <- paste0("OTU", seq(length(getSequences(merge.pspoolred.nochim))))
pspool.lulu <- t(pspool.lulu)

########### THIS NEXT PART IN THE BASH TERMINAL WITH BLAST INSTALLED IN PATH (or conda environment)

#First produce a blast databases with the OTUs
makeblastdb -in Nopool.nochim.fasta -parse_seqids -dbtype nucl
makeblastdb -in Pooled.nochim.fasta -parse_seqids -dbtype nucl
makeblastdb -in psPooled.nochim.fasta -parse_seqids -dbtype nucl
# Then blast the OTUs against the database
blastn -db Nopool.nochim.fasta -outfmt '6 qseqid sseqid pident' -out NoPool_match_list.txt -qcov_hsp_perc 80 -perc_identity 84 -query Nopool.nochim.fasta
blastn -db Pooled.nochim.fasta -outfmt '6 qseqid sseqid pident' -out Pool_match_list.txt -qcov_hsp_perc 80 -perc_identity 84 -query Pooled.nochim.fasta
blastn -db psPooled.nochim.fasta -outfmt '6 qseqid sseqid pident' -out psPool_match_list.txt -qcov_hsp_perc 80 -perc_identity 84 -query psPooled.nochim.fasta

##### Back into R
## RUNNING Lulu with default parameters
setwd("/Users/siljeleiknes/Masteroppgave/seqdata/X204SC22122693-Z01-F001")
NoPool_match_list.txt <- read.table("NoPool_match_list.txt")
str(NoPool_match_list.txt)
str(npool.lulu)
Pool_match_list.txt <- read.table("Pool_match_list.txt")
psPool_match_list.txt <- read.table("psPool_match_list.txt")
nopool.nochim.curated_result <- lulu(as.data.frame(npool.lulu), NoPool_match_list.txt)
pool.nochim.curated_result <- lulu(as.data.frame(pool.lulu), Pool_match_list.txt)
pspool.nochim.curated_result <- lulu(as.data.frame(pspool.lulu), psPool_match_list.txt)

## Check out how many OTUs were collapsed:
print(paste0("Not Pooled: ", "OTUs after Lulu: ", nopool.nochim.curated_result$curated_count, " --- ", "OTUs before Lulu: ", nrow(nopool.nochim.curated_result$original_table)))
print(paste0("Pooled: ", "OTUs after Lulu: ", pool.nochim.curated_result$curated_count, " --- ", "OTUs before Lulu: ", nrow(pool.nochim.curated_result$original_table)))
print(paste0("Pooled: ", "OTUs after Lulu: ", pspool.nochim.curated_result$curated_count, " --- ", "OTUs before Lulu: ", nrow(pspool.nochim.curated_result$original_table)))

## Making sequence tables compatible with summary routines below
## Making vector of row numbers of OTUs kept after curation - correspond to column numbers from the precurated data
rownames(nopool.nochim.curated_result$curated_table)

nopool.kept.otus <-  as.numeric(gsub("OTU", "", rownames(nopool.nochim.curated_result$curated_table), ignore.case =FALSE, perl = TRUE))
pool.kept.otus <-  as.numeric(gsub("OTU", "", rownames(pool.nochim.curated_result$curated_table), ignore.case =FALSE, perl = TRUE))
pspool.kept.otus <-  as.numeric(gsub("OTU", "", rownames(pspool.nochim.curated_result$curated_table), ignore.case =FALSE, perl = TRUE))

#t() -> transpose 
nopool.lulu <- t(nopool.nochim.curated_result$curated_table)
colnames(nopool.lulu) <- colnames(merge.npoolred.nochim[, nopool.kept.otus])
pool.lulu <- t(pool.nochim.curated_result$curated_table)
colnames(pool.lulu) <- colnames(merge.poolred.nochim[, pool.kept.otus])
pspool.lulu <- t(pspool.nochim.curated_result$curated_table)
colnames(pspool.lulu) <- colnames(merge.pspoolred.nochim[, pspool.kept.otus])

########### All steps from here will require the use of the above merge.xxxx.nochim tables for pre lulu treatment and xxx.lulu for post lulu treatment
## For raw data - removing non-biological samples (blanks, ntc, positive controls and tagcatches) and making lists of dataframes by samples for later

## Writing function:
raw.to.sample <-  function (x){
  cv <-  as.data.frame(x[!grepl("blank|NTC|lake|tagcatch|L20_30|L20_60|L20_90", rownames(x)),])
  sample.info.nopool <- bind_rows(npoolredind, .id = "column_label")
  rownames(sample.info.nopool) <- sample.info.nopool$full
  cv <-  merge(cv, sample.info.nopool, by = 'row.names', all.x=TRUE)
  rownames(cv) <- cv$Row.names
  cv <- split(cv, cv$sample)
  dropnames <- colnames(cv[[1]][, c(which(nchar(colnames(cv[[1]]))< 20))])
  cv <- lapply(cv, function(y) y[!(names(y) %in% dropnames)])
  cv <- lapply(cv,as.matrix)
  return(cv)
}

## converting raw data into sample-wise matrices with only biological samples
npoolred.sample <- raw.to.sample(merge.npoolred.nochim)
poolred.sample <- raw.to.sample(merge.poolred.nochim)
pspoolred.sample <- raw.to.sample(merge.pspoolred.nochim)

nopool.lulu.sample <- raw.to.sample(nopool.lulu)
pool.lulu.sample <- raw.to.sample(pool.lulu)
pspool.lulu.sample <- raw.to.sample(pspool.lulu)

npoolred.sample <-  merge.npoolred.nochim[!grepl("blank|NTC|lake|tagcatch|L20_30|L20_60|L20_90", rownames(merge.npoolred.nochim)),]

#### Now all raw data by sample stored in npoolred.sample, poolred.sample & pspoolred.sample
### all lulu curated data stored in nopool.lulu.sample , pool.lulu.sample, pspool.lulu.sample
### Get extraction data
getwd()
setwd("/Users/siljeleiknes/Masteroppgave/seqdata/X204SC22122693-Z01-F001")
exdata <- read.csv("extractions.csv", sep=";")
str(exdata)
length(unique(exdata$blank))
unique(exdata$blank)

## Subtracting the amount of sequences across all replicate NTC libraries within each OTU for each plate
## making function:
ntc.change <- function(x){
  mind <- apply(x[grep("NTC", rownames(x)), ], 2, function(y) max(y, na.rm = TRUE))
  x1 <- sweep(x, 2, mind)
  x1 <- pmax(x1,0)
  return(x1)
}

##### removing NTC counts - converting back to matrix first to ease later pain
npoolred <- lapply(npoolred, as.matrix)  
poolred <- lapply(poolred, as.matrix)  
pspoolred <- lapply(pspoolred, as.matrix) 

minusNTCnpool <- lapply(npoolred, ntc.change)
minusNTCpool <- lapply(poolred, ntc.change)
minusNTCpspool <- lapply(pspoolred, ntc.change)

# NOTE: To ntc control for lulu samples, need to reorganise samples into 4 lists for plate group samples were originally PCR'ed in
## function to make 4 plate list from current matrix:
raw.to.plate <-  function (x){
  ef <- merge(x, sample.info.nopool ,by =  'row.names', all.x=TRUE)
  rownames(ef) <- ef$Row.names
  ef <- split(ef, ef$column_label)
  dropnames <- colnames(ef[[1]][, c(which(nchar(colnames(ef[[1]]))< 20))])
  ef <- lapply(ef, function(y) y[!(names(y) %in% dropnames)])
  ef <- lapply(ef,as.matrix)
  return(ef)
}

plate.nopool.lulu <- raw.to.plate(nopool.lulu)
plate.pool.lulu <- raw.to.plate(pool.lulu)
plate.pspool.lulu <- raw.to.plate(pspool.lulu)

lulu.minusNTCnpool <- lapply(plate.nopool.lulu, ntc.change)
lulu.minusNTCpool <- lapply(plate.pool.lulu, ntc.change)
lulu.minusNTCpspool <- lapply(plate.pspool.lulu, ntc.change)

## For NTC controlled data - removing non-biological samples (blanks, ntc, positive controls and tagcatches) 
## and making lists of matrices by samples and by blank for later
## FUNCTION: ntc.to.sample removing non-biological samples and making list of separate matrix per sample for ntc control summaries
## FUNCTION: ntc.to.blankcontrol: making matrices out of NTC controlled data with each matrix 
## corresponding to an extraction group to control for blank extract

ntc.to.sample <-  function (x){
  z <- lapply(x, function(y) y[!grepl("blank|NTC|lake|tagcatch|L20_30|L20_60|L20_90", rownames(y)),] )
  z <- mergeSequenceTables(z[[1]], z[[2]], z[[3]], z[[4]], z[[5]])
  z <- as.data.frame(z)
  sample.info.nopool <- bind_rows(npoolredind, .id = "column_label")
  rownames(sample.info.nopool) <- sample.info.nopool$full
  z <-  merge(z, sample.info.nopool ,by =  'row.names', all.x=TRUE)
  rownames(z) <- z$Row.names
  z <- split(z, z$sample)
  dropnames <- colnames(z[[1]][, c(which(nchar(colnames(z[[1]]))< 20))])
  z <- lapply(z, function(x) x[!(names(x) %in% dropnames)])
  z <- lapply(z,as.matrix)
  return(z)
}

minusNTCnpool.sample <- ntc.to.sample(minusNTCnpool)
minusNTCpool.sample <- ntc.to.sample(minusNTCpool)
minusNTCpspool.sample <- ntc.to.sample(minusNTCpspool)

lulu.minusNTCnpool.sample <- ntc.to.sample(lulu.minusNTCnpool)
lulu.minusNTCpool.sample <- ntc.to.sample(lulu.minusNTCpool)
lulu.minusNTCpspool.sample <- ntc.to.sample(lulu.minusNTCpspool)

library(textshape)
## function to convert data to list of matrices split by blank extract batch
ntc.to.blankcontrol <-  function (x){
  l <- lapply(x, as.matrix)
  z <- mergeSequenceTables(l[[1]], l[[2]], l[[3]], l[[4]], l[[5]])
  sample.info.nopool <- bind_rows(npoolredind, .id = "column_label")
  row.names(sample.info.nopool) = sample.info.nopool$full
  z = merge(z, sample.info.nopool, by =  'row.names', all.x=TRUE)
  z = merge(z, exdata, by.x="sample", by.y="sample1", all.x=TRUE)
  z = split(z, z$blank)
  z = lapply(z, function(y) textshape::column_to_rownames(y, loc = 2))
  dropnames <- colnames(z[[1]][, c(which(nchar(colnames(z[[1]]))< 20))])
  z = lapply(z, function(x) x[!(names(x) %in% dropnames)])
  z = lapply(z, as.matrix)
  return(z)
}

minusNTCBLANKnpool <- ntc.to.blankcontrol(minusNTCnpool)
minusNTCBLANKpool <- ntc.to.blankcontrol(minusNTCpool)
minusNTCBLANKpspool <- ntc.to.blankcontrol(minusNTCpspool)
lulu.minusNTCBLANKnpool <- ntc.to.blankcontrol(lulu.minusNTCnpool)
lulu.minusNTCBLANKpool <- ntc.to.blankcontrol(lulu.minusNTCpool)
lulu.minusNTCBLANKpspool <- ntc.to.blankcontrol(lulu.minusNTCpspool)

## Blank subtraction function same as NTC subtraction function
blank.change <- function(x){
  mind <- apply(x[grep("blank|L20_30|L20_60|L20_90", rownames(x)), ], 2, function(y) max(y, na.rm = TRUE))
  x1 <- sweep(x, 2, mind)
  x1 <- pmax(x1,0)
  return(x1)
}

minusNTCBLANKnpool <- lapply(minusNTCBLANKnpool, blank.change)
minusNTCBLANKpool <- lapply(minusNTCBLANKpool, blank.change)
minusNTCBLANKpspool <- lapply(minusNTCBLANKpspool, blank.change)
lulu.minusNTCBLANKnpool <- lapply(lulu.minusNTCBLANKnpool, blank.change)
lulu.minusNTCBLANKpool <- lapply(lulu.minusNTCBLANKpool, blank.change)
lulu.minusNTCBLANKpspool <- lapply(lulu.minusNTCBLANKpspool, blank.change)

## function to remove control samples (blanks, ntc, positive controls and tagcatches) and making lists of dataframes by samples for later
ntcblank.to.samplelist <-  function (x){
  minusNTCBLANK.sample <- mergeSequenceTables(x$L20_30, x$L20_60, x$L20_90, x$L20_blank4,
                                              x$L22_blank1, x$L22_blank2, x$L22_blank3, x$L22_blank4,
                                              x$L22_blank5, x$L22_blank6, x$L22_blank7, x$L22_blank8, 
                                              x$L22_blank9, x$L22_blank10, x$L22_blank11, x$L22_blank12,
                                              x$L22_blank13, x$L22_blank14)
  minusNTCBLANK.sample <- minusNTCBLANK.sample[!grepl("blank|NTC|lake|tagcatch|L20_30|L20_60|L20_90", rownames(minusNTCBLANK.sample)),]
  sample.info.nopool <- bind_rows(npoolredind, .id = "column_label")
  rownames(sample.info.nopool) <- sample.info.nopool$full
  minusNTCBLANK.sample <-  merge(minusNTCBLANK.sample, sample.info.nopool ,by =  'row.names', all.x=TRUE)
  rownames(minusNTCBLANK.sample) <- minusNTCBLANK.sample$Row.names
  minusNTCBLANK.sample <- split(minusNTCBLANK.sample, minusNTCBLANK.sample$sample)
  dropnames <- colnames(minusNTCBLANK.sample[[1]][, c(which(nchar(colnames(minusNTCBLANK.sample[[1]]))< 20))])
  minusNTCBLANK.sample <- lapply(minusNTCBLANK.sample, function(x) x[!(names(x) %in% dropnames)])
  minusNTCBLANK.sample <- lapply(minusNTCBLANK.sample, as.matrix)
  return(minusNTCBLANK.sample)
}

minusNTCBLANKnpool.sample <- ntcblank.to.samplelist(minusNTCBLANKnpool)
minusNTCBLANKpool.sample <- ntcblank.to.samplelist(minusNTCBLANKpool)
minusNTCBLANKpspool.sample <- ntcblank.to.samplelist(minusNTCBLANKpspool)
lulu.minusNTCBLANKnpool.sample <- ntcblank.to.samplelist(lulu.minusNTCBLANKnpool)
lulu.minusNTCBLANKpool.sample <- ntcblank.to.samplelist(lulu.minusNTCBLANKpool)
lulu.minusNTCBLANKpspool.sample <- ntcblank.to.samplelist(lulu.minusNTCBLANKpspool)

#### Replicate control 
## functions to only keep columns where ASV occurs in => n (i.e. 2,3,4 replicates)

## First Removing samples where any replicate PCR has <1000 sequences
minusNTCBLANKnpool.sample.1k <- minusNTCBLANKnpool.sample[sapply(minusNTCBLANKnpool.sample, function(x) !any(rowSums(x) < 1000))]
minusNTCBLANKpool.sample.1k <- minusNTCBLANKpool.sample[sapply(minusNTCBLANKpool.sample, function(x) !any(rowSums(x) < 1000))]
minusNTCBLANKpspool.sample.1k <- minusNTCBLANKpspool.sample[sapply(minusNTCBLANKpspool.sample, function(x) !any(rowSums(x) < 1000))]
lulu.minusNTCBLANKnpool.sample.1k <- lulu.minusNTCBLANKnpool.sample[sapply(lulu.minusNTCBLANKnpool.sample, function(x) !any(rowSums(x) < 1000))]
lulu.minusNTCBLANKpool.sample.1k <- lulu.minusNTCBLANKpool.sample[sapply(lulu.minusNTCBLANKpool.sample, function(x) !any(rowSums(x) < 1000))]
lulu.minusNTCBLANKpspool.sample.1k <- lulu.minusNTCBLANKpspool.sample[sapply(lulu.minusNTCBLANKpspool.sample, function(x) !any(rowSums(x) < 1000))]

setwd("/Users/siljeleiknes")
getwd()
write.csv(lulu.minusNTCBLANKnpool.sample.1k, "lulu_1k.csv")

rep.groups2 <- function(x){
  r2 <- apply(x, 2, function(c) replace(c, sum(c!=0)<2, 0))
  return(r2)
}
rep.groups3 <- function(x){
  r3 <- apply(x, 2, function(c) replace(c, sum(c!=0)<3, 0))
  return(r3)
}
rep.groups4 <- function(x){
  r4 <- apply(x, 2, function(c) replace(c, sum(c!=0)<4, 0))
  return(r4)
}

rg2.lulu.minusNTCBLANKnpool.sample.1k<- lapply(lulu.minusNTCBLANKnpool.sample.1k, rep.groups2)
rg3.lulu.minusNTCBLANKnpool.sample.1k<- lapply(lulu.minusNTCBLANKnpool.sample.1k, rep.groups3)
rg4.lulu.minusNTCBLANKnpool.sample.1k<- lapply(lulu.minusNTCBLANKnpool.sample.1k, rep.groups4)

rg2.lulu.minusNTCBLANKpool.sample.1k<- lapply(lulu.minusNTCBLANKpool.sample.1k, rep.groups2)
rg3.lulu.minusNTCBLANKpool.sample.1k<- lapply(lulu.minusNTCBLANKpool.sample.1k, rep.groups3)
rg4.lulu.minusNTCBLANKpool.sample.1k<- lapply(lulu.minusNTCBLANKpool.sample.1k, rep.groups4)

rg2.lulu.minusNTCBLANKpspool.sample.1k<- lapply(lulu.minusNTCBLANKpspool.sample.1k, rep.groups2)
rg3.lulu.minusNTCBLANKpspool.sample.1k<- lapply(lulu.minusNTCBLANKpspool.sample.1k, rep.groups3)
rg4.lulu.minusNTCBLANKpspool.sample.1k<- lapply(lulu.minusNTCBLANKpspool.sample.1k, rep.groups4)

rg2.nopool <- lapply(minusNTCBLANKnpool.sample, rep.groups2)
rg2.pooled <- lapply(minusNTCBLANKpool.sample, rep.groups2)
rg2.pspooled <- lapply(minusNTCBLANKpspool.sample, rep.groups2)
## doing LULU replicate control on 1k plus samples
rg2.lulu.nopool.1k <- lapply(lulu.minusNTCBLANKnpool.sample.1k, rep.groups2)
#rg2.lulu.nopool.1k <- lapply(rg2.lulu.nopool.1k, function (x) x[is.nan(x)] <-0) 
rg2.lulu.pooled.1k <- lapply(lulu.minusNTCBLANKpool.sample.1k, rep.groups2)
#rg2.lulu.pooled.1k <- lapply(rg2.lulu.pooled.1k, function (x) x[is.nan(x)] <-0) 
rg2.lulu.pspooled.1k <- lapply(lulu.minusNTCBLANKpspool.sample.1k, rep.groups2)
#rg2.lulu.pspooled.1k <- lapply(rg2.lulu.pspooled.1k, function (x) x[is.nan(x)] <-0)

rg3.nopool <- lapply(minusNTCBLANKnpool.sample, rep.groups3)
rg3.pooled <- lapply(minusNTCBLANKpool.sample, rep.groups3)
rg3.pspooled <- lapply(minusNTCBLANKpspool.sample, rep.groups3)
## doing LULU replicate control on 1k plus samples
rg3.lulu.nopool.1k <- lapply(lulu.minusNTCBLANKnpool.sample.1k, rep.groups3)
#rg3.lulu.nopool.1k <- lapply(rg3.lulu.nopool.1k, function (x) x[is.nan(x)] <-0) 
rg3.lulu.pooled.1k <- lapply(lulu.minusNTCBLANKpool.sample.1k, rep.groups3)
#rg3.lulu.pooled.1k <- lapply(rg3.lulu.pooled.1k, function (x) x[is.nan(x)] <-0) 
rg3.lulu.pspooled.1k <- lapply(lulu.minusNTCBLANKpspool.sample.1k, rep.groups3)
#rg3.lulu.pspooled.1k <- lapply(rg3.lulu.pspooled.1k, function (x) x[is.nan(x)] <-0)

rg4.nopool <- lapply(minusNTCBLANKnpool.sample, rep.groups4)
rg4.pooled <- lapply(minusNTCBLANKpool.sample, rep.groups4)
rg4.pspooled <- lapply(minusNTCBLANKpspool.sample, rep.groups4)
## doing LULU replicate control on 1k plus samples
rg4.lulu.nopool.1k <- lapply(lulu.minusNTCBLANKnpool.sample.1k, rep.groups4)
#rg4.lulu.nopool.1k <- lapply(rg4.lulu.nopool.1k, function (x) x[is.nan(x)] <-0) 
rg4.lulu.pooled.1k <- lapply(lulu.minusNTCBLANKpool.sample.1k, rep.groups4)
#rg4.lulu.pooled.1k <- lapply(rg4.lulu.pooled.1k, function (x) x[is.nan(x)] <-0) 
rg4.lulu.pspooled.1k <- lapply(lulu.minusNTCBLANKpspool.sample.1k, rep.groups4)
#rg4.lulu.pspooled.1k <- lapply(rg4.lulu.pspooled.1k, function (x) x[is.nan(x)] <-0)

### Finally, assigning taxonomy to OTUs - Only doing this with LULU OTUs as that is the dataset we will ultimately use for further analyses
getwd()
setwd("/Users/siljeleiknes/Masteroppgave/seqdata/X204SC22122693-Z01-F001")

minusNTCBLANKnpool.sample <- lapply(minusNTCBLANKnpool.sample, as.matrix)

uniquesToFasta(minusNTCBLANKnpool.sample$`20_1`, fout='/Users/siljeleiknes/Masteroppgave/seqdata/X204SC22122693-Z01-F001/raw.nopool.fasta', ids=paste0("OTU", seq(length(getSequences(minusNTCBLANKnpool.sample$`20_1`)))))
uniquesToFasta(minusNTCBLANKpool.sample$`20_1`, fout='/Users/siljeleiknes/Masteroppgave/seqdata/X204SC22122693-Z01-F001/raw.pool.fasta', ids=paste0("OTU", seq(length(getSequences(minusNTCBLANKpool.sample$`20_1`)))))
uniquesToFasta(minusNTCBLANKpspool.sample$`20_1`, fout='/Users/siljeleiknes/Masteroppgave/seqdata/X204SC22122693-Z01-F001/raw.pspool.fasta', ids=paste0("OTU", seq(length(getSequences(minusNTCBLANKpspool.sample$`20_1`)))))

uniquesToFasta(lulu.minusNTCBLANKnpool.sample$`20_1`, fout='/Users/siljeleiknes/Masteroppgave/seqdata/X204SC22122693-Z01-F001/lulu.nopool.fasta', ids=paste0("OTU", seq(length(getSequences(lulu.minusNTCBLANKnpool.sample$`20_1`)))))
uniquesToFasta(lulu.minusNTCBLANKpool.sample$`20_1`, fout='/Users/siljeleiknes/Masteroppgave/seqdata/X204SC22122693-Z01-F001/lulu.pool.fasta', ids=paste0("OTU", seq(length(getSequences(lulu.minusNTCBLANKpool.sample$`20_1`)))))
uniquesToFasta(lulu.minusNTCBLANKpspool.sample$`20_1`, fout='/Users/siljeleiknes/Masteroppgave/seqdata/X204SC22122693-Z01-F001/lulu.pspool.fasta', ids=paste0("OTU", seq(length(getSequences(lulu.minusNTCBLANKpspool.sample$`20_1`)))))

### Externally blasting sequence files that are not being assigned to species with custom database
### performing Megan analysis LCA analysis (naive algorithm with default parameters)
### Exporting Megan results in csv format with OTU name to assign taxonomy back to
### Set column names OTU for OTUs and ID for the names, make OTU column into rownames
getwd()
setwd("/Users/siljeleiknes/Masteroppgave/SiljeFASTAfiles")
rawnopoolmegan <- read.csv("raw.nopool.fasta_exNAME.txt", header=FALSE)
colnames(rawnopoolmegan) <- c("OTU","ID")
rawnopoolmegan <- rawnopoolmegan %>% remove_rownames %>% column_to_rownames

rawpoolmegan <- read.csv("raw.pool.fasta_exNAME.txt", header=FALSE)
colnames(rawpoolmegan) <- c("OTU","ID")
rawpoolmegan <- rawpoolmegan %>% remove_rownames %>% column_to_rownames

rawpspoolmegan <- read.csv("raw.pspool.fasta_exNAME.txt", header=FALSE)
colnames(rawpspoolmegan) <- c("OTU","ID")
rawpspoolmegan <- rawpspoolmegan %>% remove_rownames %>% column_to_rownames

lulunopoolmegan <- read.csv("lulu.nopool.fasta_exNAME.txt", header=FALSE)
colnames(lulunopoolmegan) <- c("OTU","ID")
lulunopoolmegan <- lulunopoolmegan %>% remove_rownames %>% column_to_rownames

lulupoolmegan <- read.csv("lulu.pool.fasta_exNAME.txt", header=FALSE)
colnames(lulupoolmegan) <- c("OTU","ID")
lulupoolmegan <- lulupoolmegan %>% remove_rownames %>% column_to_rownames

lulupspoolmegan <- read.csv("lulu.pspool.fasta_exNAME.txt", header=FALSE)
colnames(lulupspoolmegan) <- c("OTU","ID")
lulupspoolmegan <- lulupspoolmegan %>% remove_rownames %>% column_to_rownames

### Using Lulu curated, DADA2 not pooled data for further analysis::
lunptax <-  read.csv("lulu.nopool.fasta.noenv-exTaxIDs.txt", header =FALSE)

#function to prepare taxID tables and return taxonomy
colrowtax <- function(x){
  lahs <- x
  colnames(lahs) <- c("OTU","taxID")
  lahs1 <- lahs
  lahs<-getTaxonomy(lahs[,2],'accessionTaxa.sql', desiredTaxa = c("superkingdom", "phylum", "subphylum", "class", "order", "family", "genus", "species"))
  rownames(lahs) <- lahs1$OTU
  return(lahs)
}

lunptax <-  colrowtax(lunptax)

## getting sequences back in
library(seqinr)
lulunpseq <- read.fasta("lulu.nopool.fasta")
row.names(lulunpseq) <- lulunpseq$seq.name
str(lulunpseq)

## Getting raw lulu matrix to do controls with 
str(nopool.lulu)

remotes::install_github("jhbadger/cipphyloseq")
remotes::install_github("helixcn/phylotools", build_vignettes = TRUE)

library(phylotools)

## Making sure that sequences in fasta file (that were used for taxonomy are allocated to the correct OTU name)
s = readDNAStringSet("lulu.nopool.fasta")
dss2df <- function(dss) data.frame(width=width(dss), seq=as.character(dss), names=names(dss))
sdf <- dss2df(s)
sdf$seqw2 <- sdf$seq

k = as.data.frame(colnames(nopool.lulu))
k$seqw2 <- k$`colnames(nopool.lulu)`
k$id  <- 1:nrow(k)

k1 <- sdf
str(k)
k3 <- merge(k,k1, by="seqw2")
k3 <- k3[order(k3$id), ]

## changing colnames of lulu object so they line up with fasta file and taxonomy
colnames(nopool.lulu) <- k3$names

##sample data
rawrepsamp <- read.csv("lulu_np_repsampledata.csv")
str(rawrepsamp)
rownames(rawrepsamp) = rawrepsamp$rep

taxa_names(controlcheck)
str(controlcheck)

controlcheck <- phyloseq(otu_table(nopool.lulu, taxa_are_rows = FALSE), sample_data(rawrepsamp), tax_table(lunptax), s)

## using callaghan decontam package to control for contamination - inspecting individual PCR library sizes
#BiocManager::install("decontam")
library(decontam)
df <- as.data.frame(sample_data(controlcheck)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(controlcheck)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=group)) + geom_point() + ylim(0,20000)

## using the prevalence of contaminant function 
sample_data(controlcheck)$is.neg <- sample_data(controlcheck)$group == "NTC" | sample_data(controlcheck)$group == "BLANKextract"| sample_data(controlcheck)$group == "External_control"

contamdf.prev <- isContaminant(controlcheck, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)

ps.pa <- transform_sample_counts(controlcheck, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(controlcheck)$group == "NTC" | sample_data(controlcheck)$group == "BLANKextract"| sample_data(controlcheck)$group == "External_control", ps.pa)
ps.pa.pos <- prune_samples(sample_data(controlcheck)$group == "Faeces_sample", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + 
  geom_point() +
  xlab("Prevalence (Negative Controls)") + 
  ylab("Prevalence (True Samples)") ## splits out samples nicely

## now removing contaminants from samples
remo <- subset(contamdf.prev, contaminant == TRUE)
remo[order(remo$prev), ]

badTaxa = row.names(remo)
goodTaxa <- setdiff(taxa_names(controlcheck), badTaxa)
controlled <- prune_taxa(goodTaxa, controlcheck)

#### Now decontaminated data is in the object that "controlled" - checking for cross-contamination of positive control
controlled

## exported to excel
write.csv(otu_table(controlled), "check.csv")
#Removed any sample with an average of postive control sequences across PCR replicates >0.01% of total sequences & in more than 2 PCR replicates

rmsamps <- c("L20_110",
             "L20_16",
             "L20_21",
             "L20_23",
             "L20_26",
             "L20_33",
             "L20_39",
             "L20_4",
             "L20_41",
             "L20_48",
             "L20_49",
             "L20_52",
             "L20_53",
             "L20_58",
             "L20_6",
             "L20_63",
             "L20_65",
             "L20_66",
             "L20_7",
             "L20_71",
             "L20_72",
             "L20_73",
             "L20_77",
             "L20_78",
             "L20_84",
             "L20_85",
             "L20_92",
             "L22_102",
             "L22_114",
             "L22_118",
             "L22_121",
             "L22_122",
             "L22_131",
             "L22_137",
             "L22_16",
             "L22_64",
             "L22_68",
             "L22_92")

## Removing samples
goodTaxa <- setdiff(sample_data(controlled)$samp, rmsamps)

controlled2  <- controlled %>%
  subset_samples(samp %in% goodTaxa)

library(microViz)
controlled2  %>% 
  tax_transform("clr") %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  ord_calc() %>% 
  ord_plot(color = "group") +
  scale_colour_brewer(palette = "Dark2")

### Note that positive control OTU is OTU3 - remove for later.
### Now re-extracting OTU table again to go through NTC & extract blank controls
lulunpcont <- otu_table(controlled2)

## Subtracting the amount of sequences across all replicate NTC libraries within each OTU for each plate
## making function:
ntc.change <- function(x){
  mind <- apply(x[grep("NTC", rownames(x)), ], 2, function(y) max(y, na.rm = TRUE))
  x1 <- sweep(x, 2, mind)
  x1 <- pmax(x1,0)
  return(x1)
}

raw.to.plate <-  function (x){
  ef <- merge(x, sample.info.nopool ,by =  'row.names', all.x=TRUE)
  rownames(ef) <- ef$Row.names
  ef <- split(ef, ef$column_label)
  ef <- lapply(ef, function(y) select(y, contains("OTU")))
  ef <- lapply(ef,as.matrix)
  return(ef)
}
library(dada2)
library(dplyr)
## converting data to platewide data for NTC control
plate.nopool.lulu <- raw.to.plate(lulunpcont)
## subtracting the largest NTC replicate count in each OTU from all samples
lulu.minusNTCnpool <- lapply(plate.nopool.lulu, ntc.change)
## now making lists of samples according to the blank they were extracted with
lulu.minusNTCBLANKnpool <- ntc.to.blankcontrol(lulu.minusNTCnpool)
## Blank control as per NTC control (subtraction of largest amount)
lulu.minusNTCBLANKnpool <- lapply(lulu.minusNTCBLANKnpool, blank.change)

## Now removing all control & non-faeces samples and making lists of samples for replicate control
lulu.minusNTCBLANKnpool.sample <- ntcblank.to.samplelist(lulu.minusNTCBLANKnpool)

### Now doing replicate control
rg2.lulu.nopool <- lapply(lulu.minusNTCBLANKnpool.sample, rep.groups2)
rg3.lulu.nopool <- lapply(lulu.minusNTCBLANKnpool.sample, rep.groups3)
rg4.lulu.nopool <- lapply(lulu.minusNTCBLANKnpool.sample, rep.groups4)

to.one.matrix <- function(x){
  lah <- do.call(rbind.data.frame, x)
  rownames(lah) <- names(x)
  colnames(lah) <- names(x[[1]])
  lah <- as.matrix(lah)
  return(lah)
}

make.phylo <- function (x, y, z, k){ # x = OTU data to this point, y = taxonomic data, z = sample data, k = sequence data
  ## Summing all the replicates to get a single row for each sample
  single <- lapply(x, function(w) colSums(w) )
  ### getting a single matrix from list of samples (with single row)
  test <- to.one.matrix(single)
  ## Giving the columns the OTU name instead of the sequence 
  # colnames(test) <- c(k$seq.name)
  ### Making phyloseq initial object
  wanted <- phyloseq(otu_table(test, taxa_are_rows = FALSE), tax_table(y), sample_data(z), k)
  return(wanted)
}

nofil <- make.phylo(lulu.minusNTCBLANKnpool.sample, tax_table(controlled2), ex1, refseq(controlled2))
rep2fil <- make.phylo(rg2.lulu.nopool, tax_table(controlled2), ex1, refseq(controlled2))
rep3fil <- make.phylo(rg3.lulu.nopool, tax_table(controlled2), ex1, refseq(controlled2))
rep4fil <- make.phylo(rg4.lulu.nopool, tax_table(controlled2), ex1, refseq(controlled2))

lulu.np.list <- list(nofil,rep2fil,rep3fil,rep4fil)

## Tidying up the phyloseqs object with phyloseq_validate (removing taxa with zero reads)
lulu.np.list <- lapply(lulu.np.list, function (x) phyloseq_validate(x, remove_undetected = TRUE))

## take a look at the data
lulu.np.list
lulu.np.list[[1]]
tax_table(lulu.np.list[[1]])
otu_table(lulu.np.list[[1]])

## fixing up taxa table so that it fills out all the way to species level, even if there is no ID out that far
lulu.np.list[[1]] <- lulu.np.list[[1]] %>%
  tax_fix(
    min_length = 4,
    unknowns = c(""),
    sep = " ", anon_unique = TRUE,
    suffix_rank = "classified")
lulu.np.list[[2]] <- lulu.np.list[[2]] %>%
  tax_fix(
    min_length = 4,
    unknowns = c(""),
    sep = " ", anon_unique = TRUE,
    suffix_rank = "classified")
lulu.np.list[[3]] <- lulu.np.list[[3]] %>%
  tax_fix(
    min_length = 4,
    unknowns = c(""),
    sep = " ", anon_unique = TRUE,
    suffix_rank = "classified")
lulu.np.list[[4]] <- lulu.np.list[[4]] %>%
  tax_fix(
    min_length = 4,
    unknowns = c(""),
    sep = " ", anon_unique = TRUE,
    suffix_rank = "classified")

## removing OTU1 & OTU18 - Why? They are lobster DNA (OTU1) and lobster DNA with INDEL (18) - OTU3 also as it is positive control DNA
badTaxa = c("OTU1", "OTU3", "OTU18")
goodTaxa <- setdiff(taxa_names(lulu.np.list[[1]] ), badTaxa)

lulu.np.list <- lapply(lulu.np.list, function (x) prune_taxa(goodTaxa, x))

## removing samples with >200 identifiable to phylum prey reads
lulu.np.list100 <- lapply(lulu.np.list, function (x) prune_samples(sample_sums(x)>=100, x))
lulu.np.list200 <- lapply(lulu.np.list, function (x) prune_samples(sample_sums(x)>=200, x))
lulu.np.list500 <- lapply(lulu.np.list, function (x) prune_samples(sample_sums(x)>=500, x))
lulu.np.list1k <- lapply(lulu.np.list, function (x) prune_samples(sample_sums(x)>=1000, x))

## installing a package with handy phyloseq functions
install.packages("remotes")
remotes::install_github("cpauvert/psadd")
library(psadd)
## Tidying up phyloseq object
lulu.np.list500 <- lapply(lulu.np.list500, function (x) subset_samples_no_zero(x))
lulu.np.list1k <- lapply(lulu.np.list1k, function (x) subset_samples_no_zero(x))

## Now just getting metazoan diet taxa (not fungi, protists parasites etc) . NOTE - if, for example interested in eukaryotic microbiome instead of diet, filter different phylums
lulu.np.diet100<- lapply(lulu.np.list100, function (x) subset_taxa(x, phylum == "Annelida"| phylum =="Arthropoda"|phylum =="Bryozoa"|phylum =="Chlorophyta"|
                                                                     phylum =="Chordata"|phylum =="Cnidaria"|phylum =="Echinodermata"|phylum =="Haptophyta"|phylum =="Kinorhyncha"|
                                                                     phylum =="Mollusca"|phylum =="Nematoda"|phylum =="Nemertea"|phylum =="Porifera"|phylum =="Rhodophyta"|
                                                                     phylum =="Streptophyta"))
lulu.np.diet200 <- lapply(lulu.np.list200, function (x) subset_taxa(x, phylum == "Annelida"| phylum =="Arthropoda"|phylum =="Bryozoa"|phylum =="Chlorophyta"|
                                                                      phylum =="Chordata"|phylum =="Cnidaria"|phylum =="Echinodermata"|phylum =="Haptophyta"|phylum =="Kinorhyncha"|
                                                                      phylum =="Mollusca"|phylum =="Nematoda"|phylum =="Nemertea"|phylum =="Porifera"|phylum =="Rhodophyta"|
                                                                      phylum =="Streptophyta"))
lulu.np.diet500 <- lapply(lulu.np.list500, function (x) subset_taxa(x, phylum == "Annelida"| phylum =="Arthropoda"|phylum =="Bryozoa"|phylum =="Chlorophyta"|
                                                                      phylum =="Chordata"|phylum =="Cnidaria"|phylum =="Echinodermata"|phylum =="Haptophyta"|phylum =="Kinorhyncha"|
                                                                      phylum =="Mollusca"|phylum =="Nematoda"|phylum =="Nemertea"|phylum =="Porifera"|phylum =="Rhodophyta"|
                                                                      phylum =="Streptophyta"))
lulu.np.diet1k <- lapply(lulu.np.list1k, function (x) subset_taxa(x, phylum == "Annelida"| phylum =="Arthropoda"|phylum =="Bryozoa"|phylum =="Chlorophyta"|
                                                                    phylum =="Chordata"|phylum =="Cnidaria"|phylum =="Echinodermata"|phylum =="Haptophyta"|phylum =="Kinorhyncha"|
                                                                    phylum =="Mollusca"|phylum =="Nematoda"|phylum =="Nemertea"|phylum =="Porifera"|phylum =="Rhodophyta"|
                                                                    phylum =="Streptophyta"))
## Now removing samples with >500 & 1000 DIET reads
lulu.np.diet100 <- lapply(lulu.np.diet100, function (x) prune_samples(sample_sums(x)>=100, x))
lulu.np.diet200 <- lapply(lulu.np.diet200, function (x) prune_samples(sample_sums(x)>=200, x))
lulu.np.diet500 <- lapply(lulu.np.diet500, function (x) prune_samples(sample_sums(x)>=500, x))
lulu.np.diet1k <- lapply(lulu.np.diet1k, function (x) prune_samples(sample_sums(x)>=1000, x))

# making sure there are no samples with zero reads
lulu.np.diet100 <- lapply(lulu.np.diet100, function (x) subset_samples_no_zero(x))
lulu.np.diet200 <- lapply(lulu.np.diet200, function (x) subset_samples_no_zero(x))
lulu.np.diet500 <- lapply(lulu.np.diet500, function (x) subset_samples_no_zero(x))
lulu.np.diet1k <- lapply(lulu.np.diet1k, function (x) subset_samples_no_zero(x))
## after previous filtering, again, removing taxa with zero reads
lulu.np.diet100 <- lapply(lulu.np.diet100, function (x) phyloseq_validate(x, remove_undetected = TRUE))
lulu.np.diet200 <- lapply(lulu.np.diet200, function (x) phyloseq_validate(x, remove_undetected = TRUE))
lulu.np.diet500 <- lapply(lulu.np.diet500, function (x) phyloseq_validate(x, remove_undetected = TRUE))
lulu.np.diet1k <- lapply(lulu.np.diet1k, function (x) phyloseq_validate(x, remove_undetected = TRUE))

lulu.np.diet100 ## all 4 replicate controls with samples with <100 prey sequences discarded
lulu.np.diet200 ## all 4 replicate controls with samples with <200 prey sequences discarded
lulu.np.diet500 ## all 4 replicate controls with samples with <500 prey sequences discarded
lulu.np.diet1k ## all 4 replicate controls with samples with <1000 prey sequences discarded

## plots to look at library sizes, number of OTUs per sample etc 
df <- as.data.frame(sample_data(lulu.np.diet500[[3]])) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(lulu.np.diet500[[3]])
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
boxplot(df$LibrarySize~df$yrtrt)
ggplot(data=df, aes(x=Index, y=LibrarySize, color=trt)) + geom_point() + ylim(0, 20000)

## counts of OTUs per sample
remotes::install_github("vmikk/metagMisc")
library(metagMisc)
OTU_counts <- phyloseq_ntaxa_by_tax(lulu.np.diet500[[3]], TaxRank = "superkingdom")
OTU_counts$lsize <- sample_sums(lulu.np.diet500[[3]])

## relationship between sequencing depth and diet diversity
ggplot(data=OTU_counts, aes(x=N.OTU, y=lsize, color=yrtrt)) + geom_point() # + ylim(0, 20000)

## Creating a phylogenetic tree for further analysis
library(DECIPHER)

## Aligning  - doing this without lapply.. too complicated running through 100, 200, 500 & 1k seperately
lulu.np.all <- AlignSeqs(DNAStringSet(refseq(lulu.np.diet100[[1]])), anchor=NA)
lulu.np.r2 <- AlignSeqs(DNAStringSet(refseq(lulu.np.diet100[[2]])), anchor=NA)
lulu.np.r3 <- AlignSeqs(DNAStringSet(refseq(lulu.np.diet100[[3]])), anchor=NA)
lulu.np.r4 <- AlignSeqs(DNAStringSet(refseq(lulu.np.diet100[[4]])), anchor=NA)

## making maximum likelihood trees from each alignment
library(phangorn)
phang.align <- phyDat(as(lulu.np.all, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)
fitGTR <- update(fit, k=4, inv=0.2)
lulu.np.all.fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                                rearrangement = "stochastic", control = pml.control(trace = 0))

phang.align <- phyDat(as(lulu.np.r2, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)
fitGTR <- update(fit, k=4, inv=0.2)
lulu.np.r2.fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                               rearrangement = "stochastic", control = pml.control(trace = 0))

phang.align <- phyDat(as(lulu.np.r3, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)
fitGTR <- update(fit, k=4, inv=0.2)
lulu.np.r3.fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                               rearrangement = "stochastic", control = pml.control(trace = 0))

phang.align <- phyDat(as(lulu.np.r4, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)
fitGTR <- update(fit, k=4, inv=0.2)
lulu.np.r4.fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                               rearrangement = "stochastic", control = pml.control(trace = 0))

## Adding trees

lulu.np.diet100[[1]]<- merge_phyloseq(lulu.np.diet100[[1]], lulu.np.all.fitGTR$tree)
lulu.np.diet100[[2]]<- merge_phyloseq(lulu.np.diet100[[2]], lulu.np.r2.fitGTR$tree)
lulu.np.diet100[[3]]<- merge_phyloseq(lulu.np.diet100[[3]], lulu.np.r3.fitGTR$tree)
lulu.np.diet100[[4]]<- merge_phyloseq(lulu.np.diet100[[4]], lulu.np.r4.fitGTR$tree)

####################
lulu.np.all <- AlignSeqs(DNAStringSet(refseq(lulu.np.diet200[[1]])), anchor=NA)
lulu.np.r2 <- AlignSeqs(DNAStringSet(refseq(lulu.np.diet200[[2]])), anchor=NA)
lulu.np.r3 <- AlignSeqs(DNAStringSet(refseq(lulu.np.diet200[[3]])), anchor=NA)
lulu.np.r4 <- AlignSeqs(DNAStringSet(refseq(lulu.np.diet200[[4]])), anchor=NA)

## making maximum likelihood trees from each alignment 
phang.align <- phyDat(as(lulu.np.all, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)
fitGTR <- update(fit, k=4, inv=0.2)
lulu.np.all.fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                                rearrangement = "stochastic", control = pml.control(trace = 0))

phang.align <- phyDat(as(lulu.np.r2, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)
fitGTR <- update(fit, k=4, inv=0.2)
lulu.np.r2.fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                               rearrangement = "stochastic", control = pml.control(trace = 0))

phang.align <- phyDat(as(lulu.np.r3, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)
fitGTR <- update(fit, k=4, inv=0.2)
lulu.np.r3.fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                               rearrangement = "stochastic", control = pml.control(trace = 0))

phang.align <- phyDat(as(lulu.np.r4, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)
fitGTR <- update(fit, k=4, inv=0.2)
lulu.np.r4.fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                               rearrangement = "stochastic", control = pml.control(trace = 0))


## Adding trees
lulu.np.diet200[[1]]<- merge_phyloseq(lulu.np.diet200[[1]], lulu.np.all.fitGTR$tree)
lulu.np.diet200[[2]]<- merge_phyloseq(lulu.np.diet200[[2]], lulu.np.r2.fitGTR$tree)
lulu.np.diet200[[3]]<- merge_phyloseq(lulu.np.diet200[[3]], lulu.np.r3.fitGTR$tree)
lulu.np.diet200[[4]]<- merge_phyloseq(lulu.np.diet200[[4]], lulu.np.r4.fitGTR$tree)

##################
## Aligning  - doing this without lapply.. too complicated running through 500 
lulu.np.all <- AlignSeqs(DNAStringSet(refseq(lulu.np.diet500[[1]])), anchor=NA)
lulu.np.r2 <- AlignSeqs(DNAStringSet(refseq(lulu.np.diet500[[2]])), anchor=NA)
lulu.np.r3 <- AlignSeqs(DNAStringSet(refseq(lulu.np.diet500[[3]])), anchor=NA)
lulu.np.r4 <- AlignSeqs(DNAStringSet(refseq(lulu.np.diet500[[4]])), anchor=NA)

## making maximum likelihood trees from each alignment 
phang.align <- phyDat(as(lulu.np.all, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)
fitGTR <- update(fit, k=4, inv=0.2)
lulu.np.all.fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                                rearrangement = "stochastic", control = pml.control(trace = 0))

phang.align <- phyDat(as(lulu.np.r2, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)
fitGTR <- update(fit, k=4, inv=0.2)
lulu.np.r2.fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                               rearrangement = "stochastic", control = pml.control(trace = 0))

phang.align <- phyDat(as(lulu.np.r3, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)
fitGTR <- update(fit, k=4, inv=0.2)
lulu.np.r3.fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                               rearrangement = "stochastic", control = pml.control(trace = 0))

phang.align <- phyDat(as(lulu.np.r4, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)
fitGTR <- update(fit, k=4, inv=0.2)
lulu.np.r4.fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                               rearrangement = "stochastic", control = pml.control(trace = 0))


## Adding trees
lulu.np.diet500[[1]]<- merge_phyloseq(lulu.np.diet500[[1]], lulu.np.all.fitGTR$tree)
lulu.np.diet500[[2]]<- merge_phyloseq(lulu.np.diet500[[2]], lulu.np.r2.fitGTR$tree)
lulu.np.diet500[[3]]<- merge_phyloseq(lulu.np.diet500[[3]], lulu.np.r3.fitGTR$tree)
lulu.np.diet500[[4]]<- merge_phyloseq(lulu.np.diet500[[4]], lulu.np.r4.fitGTR$tree)


##################
## Aligning  - doing this without lapply.. too complicated running through 1k
lulu.np.all <- AlignSeqs(DNAStringSet(refseq(lulu.np.diet1k[[1]])), anchor=NA)
lulu.np.r2 <- AlignSeqs(DNAStringSet(refseq(lulu.np.diet1k[[2]])), anchor=NA)
lulu.np.r3 <- AlignSeqs(DNAStringSet(refseq(lulu.np.diet1k[[3]])), anchor=NA)
lulu.np.r4 <- AlignSeqs(DNAStringSet(refseq(lulu.np.diet1k[[4]])), anchor=NA)

## making maximum likelihood trees from each alignment 
phang.align <- phyDat(as(lulu.np.all, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)
fitGTR <- update(fit, k=4, inv=0.2)
lulu.np.all.fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                                rearrangement = "stochastic", control = pml.control(trace = 0))

phang.align <- phyDat(as(lulu.np.r2, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)
fitGTR <- update(fit, k=4, inv=0.2)
lulu.np.r2.fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                               rearrangement = "stochastic", control = pml.control(trace = 0))

phang.align <- phyDat(as(lulu.np.r3, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)
fitGTR <- update(fit, k=4, inv=0.2)
lulu.np.r3.fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                               rearrangement = "stochastic", control = pml.control(trace = 0))

phang.align <- phyDat(as(lulu.np.r4, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)
fitGTR <- update(fit, k=4, inv=0.2)
lulu.np.r4.fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                               rearrangement = "stochastic", control = pml.control(trace = 0))

detach("package:phangorn", unload=TRUE)

## Adding trees
lulu.np.diet1k[[1]]<- merge_phyloseq(lulu.np.diet1k[[1]], lulu.np.all.fitGTR$tree)
lulu.np.diet1k[[2]]<- merge_phyloseq(lulu.np.diet1k[[2]], lulu.np.r2.fitGTR$tree)
lulu.np.diet1k[[3]]<- merge_phyloseq(lulu.np.diet1k[[3]], lulu.np.r3.fitGTR$tree)
lulu.np.diet1k[[4]]<- merge_phyloseq(lulu.np.diet1k[[4]], lulu.np.r4.fitGTR$tree)

# Move on to statistical analysis and data visualization from here
