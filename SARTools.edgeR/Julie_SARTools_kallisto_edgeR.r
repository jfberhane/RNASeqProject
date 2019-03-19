################################################################################
### R script to compare several conditions with the SARTools and edgeR packages
### Hugo Varet
### March 20th, 2018
### designed to be executed with SARTools 1.6.6
################################################################################

################################################################################
###                parameters: to be modified by the user                    ###
################################################################################
rm(list=ls())                                        # remove all the objects from the R session

workDir <- setwd("~/Documents/RNASeqProject/SARTools.edgeR")      # working directory for the R session

projectName <- "SARTools.kallisto.edgeR"                         # name of the project
author <- "Julie Berhane"                                # author of the statistical analysis/report

targetFile <- "~/Documents/RNASeqProject/SARTools.edgeR/edgeR.target.txt"                           # path to the design/target file
rawDir <- "~/Documents/RNASeqProject/SARTools.edgeR"                                      # path to the directory containing raw counts files
featuresToRemove <- c("NULL")# NULL if no feature to remove

varInt <- "Treatment"                                    # factor of interest
condRef <- "Untreated"                                      # reference biological condition
batch <- NULL                                        # blocking factor: NULL (default) or "batch" for example

idColumn = 1                                         # column with feature Ids (usually 1)
countColumn = 4                                      # column with counts  (2 for htseq-count, 7 for featurecounts, 5 for RSEM/Salmon, 4 for kallisto)
rowSkip = 0                                          # rows to skip (not including header) 

alpha <- 0.05                                        # threshold of statistical significance
pAdjustMethod <- "BH"                                # p-value adjustment method: "BH" (default) or "BY"

cpmCutoff <- 1                                       # counts-per-million cut-off to filter low counts
gene.selection <- "pairwise"                         # selection of the features in MDSPlot
normalizationMethod <- "TMM"                         # normalization method: "TMM" (default), "RLE" (DESeq) or "upperquartile"

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
checkParameters.edgeR(projectName=projectName,author=author,targetFile=targetFile,
                      rawDir=rawDir,featuresToRemove=featuresToRemove,varInt=varInt,
                      condRef=condRef,batch=batch,alpha=alpha,pAdjustMethod=pAdjustMethod,
                      cpmCutoff=cpmCutoff,gene.selection=gene.selection,
                      normalizationMethod=normalizationMethod,colors=colors)

# loading target file
target <- loadTargetFile(targetFile=targetFile, varInt=varInt, condRef=condRef, batch=batch)

# loading counts
counts <- loadCountData(target=target, rawDir=rawDir, featuresToRemove=featuresToRemove, 
                        skip=rowSkip, idColumn=idColumn, countColumn=countColumn)

# description plots
majSequences <- descriptionPlots(counts=counts, group=target[,varInt], col=colors)

# edgeR analysis
out.edgeR <- run.edgeR(counts=counts, target=target, varInt=varInt, condRef=condRef,
                       batch=batch, cpmCutoff=cpmCutoff, normalizationMethod=normalizationMethod,
                       pAdjustMethod=pAdjustMethod)

# MDS + clustering
exploreCounts(object=out.edgeR$dge, group=target[,varInt], gene.selection=gene.selection, col=colors)

# summary of the analysis (boxplots, dispersions, export table, nDiffTotal, histograms, MA plot)
summaryResults <- summarizeResults.edgeR(out.edgeR, group=target[,varInt], counts=counts, alpha=alpha, col=colors)

# save image of the R session
save.image(file=paste0(projectName, ".RData"))

# generating HTML report
writeReport.edgeR(target=target, counts=counts, out.edgeR=out.edgeR, summaryResults=summaryResults,
                  majSequences=majSequences, workDir=workDir, projectName=projectName, author=author,
                  targetFile=targetFile, rawDir=rawDir, featuresToRemove=featuresToRemove, varInt=varInt,
                  condRef=condRef, batch=batch, alpha=alpha, pAdjustMethod=pAdjustMethod, cpmCutoff=cpmCutoff,
                  colors=colors, gene.selection=gene.selection, normalizationMethod=normalizationMethod)
