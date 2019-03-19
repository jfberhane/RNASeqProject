################################################################################
### R script to compare several conditions with the SARTools and DESeq2 packages
### Hugo Varet
### March 20th, 2018
### designed to be executed with SARTools 1.6.6
################################################################################

################################################################################
###                parameters: to be modified by the user                    ###
################################################################################
rm(list=ls())                                        # remove all the objects from the R session

workDir <- setwd("~/Documents/RNASeqProject/SARTools.DESeq2.transcripts.batch")      # working directory for the R session

projectName <- "SARTools.DESeq2.genes"                         # name of the project
author <- "Julie Berhane"                                # author of the statistical analysis/report

targetFile <- "./transcripts.target.txt"                           # path to the design/target file
rawDir <- "./"                                      # path to the directory containing raw counts files
featuresToRemove <- c("NULL")        # names of the features to be removed
                      # (specific HTSeq-count information and rRNA for example)
                      # NULL if no feature to remove

varInt <- "Treatment"                                    # factor of interest
condRef <- "Untreated"                                      # reference biological condition
batch <- "batch"                                        # blocking factor: NULL (default) or "batch" for example

idColumn = 1                                         # column with feature Ids (usually 1)
countColumn = 5                                      # column with counts  (2 for htseq-count, 7 for featurecounts, 5 for RSEM/Salmon, 4 for kallisto)
rowSkip = 0                                          # rows to skip (not including header) 

fitType <- "parametric"                              # mean-variance relationship: "parametric" (default), "local" or "mean"
cooksCutoff <- TRUE                                  # TRUE/FALSE to perform the outliers detection (default is TRUE)
independentFiltering <- TRUE                         # TRUE/FALSE to perform independent filtering (default is TRUE)
alpha <- 0.05                                        # threshold of statistical significance
pAdjustMethod <- "BH"                                # p-value adjustment method: "BH" (default) or "BY"

typeTrans <- "VST"                                   # transformation for PCA/clustering: "VST" or "rlog"
locfunc <- "median"                                  # "median" (default) or "shorth" to estimate the size factors


colors <- c("dodgerblue","firebrick1",               # vector of colors of each biological condition on the plots
            "MediumVioletRed","SpringGreen")

forceCairoGraph <- FALSE

################################################################################
###                             running script                               ###
################################################################################
setwd(workDir)

if (!require("BiocManager")) install.packages("BiocManager"); library(BiocManager)
if (!require("DESeq2")) BiocManager::install("DESeq2"); library(DESeq2)
if (!require("edgeR")) BiocManager::install("edgeR"); library(edgeR)
if (!require("genefilter")) BiocManager::install("genefilter"); library(genefilter)

# PC Users only, install Rtools https://cran.r-project.org/bin/windows/Rtools/

if (!require("devtools")) install.packages("devtools"); library(devtools)
if (!require("SARTools")) install_github("KField-Bucknell/SARTools", build_vignettes=TRUE, force=TRUE); library(SARTools)
if (forceCairoGraph) options(bitmapType="cairo")

# checking parameters
checkParameters.DESeq2(projectName=projectName,author=author,targetFile=targetFile,
                       rawDir=rawDir,featuresToRemove=featuresToRemove,varInt=varInt,
                       condRef=condRef,batch=batch,fitType=fitType,cooksCutoff=cooksCutoff,
                       independentFiltering=independentFiltering,alpha=alpha,pAdjustMethod=pAdjustMethod,
                       typeTrans=typeTrans,locfunc=locfunc,colors=colors)

# loading target file
target <- loadTargetFile(targetFile=targetFile, varInt=varInt, condRef=condRef, batch=batch)

# loading counts
counts <- loadCountData(target=target, rawDir=rawDir, featuresToRemove=featuresToRemove, 
                        skip=rowSkip, idColumn=idColumn, countColumn=countColumn)

# description plots
majSequences <- descriptionPlots(counts=counts, group=target[,varInt], col=colors)

# analysis with DESeq2
out.DESeq2 <- run.DESeq2(counts=counts, target=target, varInt=varInt, batch=batch,
                         locfunc=locfunc, fitType=fitType, pAdjustMethod=pAdjustMethod,
                         cooksCutoff=cooksCutoff, independentFiltering=independentFiltering, alpha=alpha)

# PCA + clustering
exploreCounts(object=out.DESeq2$dds, group=target[,varInt], typeTrans=typeTrans, col=colors)

# summary of the analysis (boxplots, dispersions, diag size factors, export table, nDiffTotal, histograms, MA plot)
summaryResults <- summarizeResults.DESeq2(out.DESeq2, group=target[,varInt], col=colors,
                                          independentFiltering=independentFiltering,
                                          cooksCutoff=cooksCutoff, alpha=alpha)

# save image of the R session
save.image(file=paste0(projectName, ".RData"))

# generating HTML report
writeReport.DESeq2(target=target, counts=counts, out.DESeq2=out.DESeq2, summaryResults=summaryResults,
                   majSequences=majSequences, workDir=workDir, projectName=projectName, author=author,
                   targetFile=targetFile, rawDir=rawDir, featuresToRemove=featuresToRemove, varInt=varInt,
                   condRef=condRef, batch=batch, fitType=fitType, cooksCutoff=cooksCutoff,
                   independentFiltering=independentFiltering, alpha=alpha, pAdjustMethod=pAdjustMethod,
                   typeTrans=typeTrans, locfunc=locfunc, colors=colors)

