setwd("/home/jnguy254/projects/def-ccastel/SharedResources/CLSA_private/data/epigenetics")

#Load packages
suppressMessages(library(minfi))
suppressMessages(library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19))

############################
### Load data into minfi ###
############################

#Read in the sample sheet for the experiment
dataDirectory <-  "/home/jnguy254/projects/def-ccastel/SharedResources/CLSA_private/data/epigenetics"
pdata <- read.metharray.sheet(dataDirectory, pattern="CLSA_metadata_FINAL_Apr2021-deidentified.csv")
pdata$Basename <- gsub("_Red.idat", "", pdata$IDAT_FILE_NAME_Red_Cy5)


#Read in the raw data from the IDAT files
rgSet <- read.metharray.exp(base = dataDirectory, targets=pdata)
sampleNames(rgSet) <- pdata$ADM_EPIGEN2_COM

#Calculate the detection p-values
detP <- suppressMessages(detectionP(rgSet))

############################
#Track how many samples were removed due to poor quality
samples_before <- dim(rgSet)[2]

#Remove poor quality samples
keep <- colMeans(detP) < 0.05
rgSet <- rgSet[,keep]

#Print out number of samples removed due to poor quality
samples_removed <- samples_before - dim(detP)[2]
message("----- ", samples_removed, " sample(s) removed due to poor quality") #----- 0 sample(s) removed due to poor quality


############################
#Normalize the data
mSetSq <- suppressMessages(preprocessNoob(rgSet))
mSetSq <- ratioConvert(mSetSq)


############################
#Ensure probes are in the same order in the mSetSq and detP objects
detP <- detP[match(featureNames(mSetSq),rownames(detP)),]


#Remove any probes that have failed in one or more samples
probes_before <- dim(mSetSq)[1]
keep <- rowSums(detP < 0.01) == ncol(mSetSq)
mSetSqFlt <- mSetSq[keep,]

#Print out number of probes removed for failing in one or more samples
probes_removed <- probes_before - dim(mSetSqFlt)[1]
message("----- ", probes_removed, " probe(s) removed for failing in one or more samples") #----- 35533 probe(s) removed for failing in one or more samples
probes_before <- dim(mSetSqFlt)[1]

############################
#Remove probes with SNPs at CpG site
mSetSqFlt <- dropLociWithSnps(mSetSqFlt)

#Print out number of probes removed for having SNPs at CpG site
probes_removed <- probes_before - dim(mSetSqFlt)[1]
message("----- ", probes_removed, " probe(s) removed for having SNPs at CpG site") #----- 25504 probe(s) removed for having SNPs at CpG site
probes_before <- dim(mSetSqFlt)[1]

############################
#Exclude cross reactive probes
xReactiveProbes <- read.csv(file=paste(dataDirectory,
                                       "PidsleyCrossReactiveEPIC.csv",
                                       sep="/"), stringsAsFactors=FALSE)
keep <- !(featureNames(mSetSqFlt) %in% xReactiveProbes$TargetID)
mSetSqFlt <- mSetSqFlt[keep,]

#Print out number of probes removed for being cross reactive
probes_removed <- probes_before - dim(mSetSqFlt)[1]
message("----- ", probes_removed, " probe(s) removed for being cross reactive") #----- 24846 probe(s) removed for being cross reactive
probes_before <- dim(mSetSqFlt)[1]


############################
#Remove Sex Probes
ann <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
sexProbes <- ann[which(ann$chr %in% c("chrX", "chrY")),]
keep <- !(featureNames(mSetSqFlt) %in% sexProbes$Name)
mSetSqFlt <- mSetSqFlt[keep,]

#Print out number of probes removed for being cross reactive
probes_removed <- probes_before - dim(mSetSqFlt)[1]
message("----- ", probes_removed, " probe(s) removed for being on sex chromosomes") #----- 16805 probe(s) removed for being on sex chromosomes



#Print out the number of probes remaining
message("----- ", dim(mSetSqFlt)[1], " probe(s) remaining for analysis") #----- 763171 probe(s) remaining for analysis


############################
#Calculate methylation beta values
bVals <- getBeta(mSetSqFlt)

write.csv(bVals, file = "../../results/epigenetics/bVals_Noob.csv", quotes = FALSE, rownames = TRUE)


