# My R Script
# Author: Pavel Artemov
# Date: 2024-02-05
# Description: This script divides pchic rds into subfiles and also exports the distance file for imputation
# Usage: Rscript chr_split.R pchic.rds design/ distout.rds split_pchic_filename_prefix_
# arq1 is pchic rds file input
# arg2 is design directory needs to end in /
# arg3 is distance file out .rds
# arg4 is filename prefix of file output 

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=4) {
  stop("Wrong number of arguements", call.=FALSE) } else if (length(args)==4) {
  
  print('Arguements supplied succesfully to pchic split script succefully')}

#pchic_rds= '~/../projects/lms-spivakov-analysis/live/artemov/ILC3/K562/K562_5kb_solbaits_merged_reweighted_Step2.Rds'
#design= '~/../projects/lms-spivakov-analysis/live/Design/Human_hg38_bin5K_sol_baits/'
#dist_out= '~/../projects/lms-spivakov-analysis/live/artemov/ILC3/K562/PCHIC_5kb/K562_dist_5kb.rds'
#split_pchic='~/../projects/lms-spivakov-analysis/live/artemov/ILC3/K562/PCHIC_5kb/K562_5kb_'

library(dplyr)
library(stringr)
library(data.table)
#library(devtools)
library(Chicago)

args = commandArgs(trailingOnly=TRUE)
pchic <- file.path(args[1])
testDesignDir <- file.path(args[2])
bmap = list.files(testDesignDir,pattern = '*.baitmap$')[1]
rmap = list.files(testDesignDir,pattern = '*.rmap$')[1]
dist_out <- args[3]
split_pchic <- args[4]
print('reading PCHIC')
pchic <- readRDS(pchic)


pchic <- modifySettings(pchic, designDir=testDesignDir)
print("PCHiC settings modified")

if(is.null(pchic@params$distFunParams)){
  print('Estimating Distance Parameters')
  pchic <- estimateDistFun(pchic,plot=FALSE)
}else{
  print('Distance Parameters already estimated')
}

exclude_cols <- c("numPairs", "nTrans", "otherEndChr", "baitChr", "N.3", "transBaitLen")

if(("otherEndChr" %in% names(pchic@x))){
  print('Removing trans interactions')
  pchic@x <- pchic@x[otherEndChr==baitChr]
}else{
  print('No trans interactions')
}

if(sum(exclude_cols %in% names(pchic@x))>0){
  print('Removing unneded columns')
  desired_cols <- setdiff(colnames(pchic@x), exclude_cols)
  pchic@x <- pchic@x[, ..desired_cols, with = FALSE]
}else{
  print('No unneeded columns in CHiCAGO RDS')
}


print('The distance parameters are:')
print(pchic@params$distFunParams)
print(paste0('Saving them at ',dist_out))
saveRDS(pchic@params$distFunParams,dist_out)
print('Distance parameters saved.')

print(paste('The baitmap is',bmap))
bmap <- fread(paste0(testDesignDir,bmap))
merged <- merge(pchic@x,bmap, by.y = 'V4', by.x = 'baitID', all.x = TRUE)
setnames(merged,'V1','baitChr')
print(paste('The rmap is',rmap))
rmap <- fread(paste0(testDesignDir,rmap))
merged2 <- merge(merged,rmap, by.y = 'V4', by.x = 'otherEndID', all.x = TRUE)
setnames(merged2,'V1','otherEndChr')
merged2 <- merged2[baitChr == otherEndChr]
chromosomes <- unique(merged$baitChr)
for (i in chromosomes){
  small <- merged2[baitChr == i]
  print(paste('Saving pchic chromosome',i))
  saveRDS(small, paste0(split_pchic,i,'.rds'))
  print('Moving onto next chromosome.')
}
