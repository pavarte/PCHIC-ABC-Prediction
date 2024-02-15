# My R Script
# Author: Pavel Artemov
# Date: 2024-02-05
# Description: This script divides pchic rds into subfiles and also exports the distance file for imputation

# Usage: Rscript imputation_script.R ${pchic.rds} ${design} ${enhancerdir/} ${split_pchic} ${distout.rds} ${imputed_pchic_prefix}
# arg1 is pchic.rds, RDS file from step 2 of running 
# arg2 is testDesignDir, needs to end in / 
# arq3 is input_cand_dir, needs to end in /
# arg4 is split_pchic, filename prefix only
# arg5 is RDS file, distance parameters
# arg6 is output directory, needs to end in /  

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=6) {
  stop("Wrong number of arguements", call.=FALSE) } else if (length(args)==4) {

  print('Arguements supplied succesfully to imputation script')}


library(dplyr)
library(stringr)
library(data.table)
library(Chicago)

##### Input some arguements for imputation

maxWindow = 5000000
minWindow = 0
chromosomes = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X")
pchic.rds = args[1]
testDesignDir = args[2]
input_cand_dir = args[3]
split_pchic = args[4]
distout.rds = args[5]
pchic_out_prefix = file.path(args[6])

##### Splitting Enhancer and Gene Lists by chromosome
chr_split_command <- paste0("Rscript chr_split.r ",input_cand_dir)
# Execute the command
system(chr_split_command)
message("Succefully finished chromosome split")
##### Splitting PCHiC by chromosome
chic_split_command <- paste("Rscript chic_split.r",pchic.rds, testDesignDir, distout.rds, split_pchic)
# Execute the command
system(chic_split_command)
message("Succefully finished PCHiC split")
##### Examples of arguments
#pchic.rds <- readRDS("~/../projects/lms-spivakov-analysis/live/artemov/ILC3/ILC3/input_data/contact/ILC/hILC3_20K_merged_Step2.Rds")
#testDesignDir <- file.path("~/../projects/lms-spivakov-analysis/live/artemov/Design_for_PCHIC/Human_hg38_bin5K_sol_baits/")
#input_cand_dir = file.path('~/../projects/lms-spivakov-analysis/live/artemov/ILC3/ILC3/ABC_output_neighb/')
#split_pchic = file.path('~/../projects/lms-spivakov-analysis/live/artemov/CD4_for_swap/5kbp/cd4_5kbp_')
#distout.rds = readRDS('~/../projects/lms-spivakov-analysis/live/artemov/CD4_for_swap/5kbp/cd4_5kbp_dist.rds')
#pchic_out_prefix = file.path('~/../projects/lms-spivakov-analysis/live/artemov/abc_swap/imputed_data/Cc_Ia_5/')
##### Creating candidate pairs from Candidate Enhancers and Candidate Genes

#reading in the distance parameters
distout.rds = readRDS(args[5])

# preparing folders for the input
for (i in chromosomes){
  baitmap=list.files(path = testDesignDir,pattern = "\\.baitmap$")[1]
  rmap=list.files(path = testDesignDir, pattern = "\\.rmap$")[1]
  candidate_enhancers <- fread(paste0(input_cand_dir,'EnhancerList_chr',i,'.txt'))
  candidate_genes <- fread(paste0(input_cand_dir,'GeneList_chr',i,'.txt'))
  setDT(candidate_enhancers)
  setDT(candidate_genes)
  candidate_enhancers[,chr:=str_split(chr,'chr', simplify = TRUE)[,2]]
  candidate_genes[,chr:=str_split(chr,'chr', simplify = TRUE)[,2]]
  candidate_enhancers[,chr:=as.character(chr)]
  candidate_genes[,chr:=as.character(chr)] 
  print(head(candidate_genes))
  print(head(candidate_enhancers))
  
  
  #### Reading baitmap and rmap
  baitmap=list.files(path = testDesignDir, pattern = "\\.baitmap$")[1]
  rmap=list.files(path = testDesignDir, pattern = "\\.rmap$")[1]
  baitmap <- fread(paste0(testDesignDir,baitmap),sep ='\t',header=FALSE)
  rmap <- fread(paste0(testDesignDir,rmap),sep ='\t',header=FALSE)
  id_space = merge(baitmap,rmap,by = c('V1','V2','V3','V4'), all.y = TRUE)
  id_space[is.na(V5),V5 := '.']
  id_space[,V1:= as.character(V1)]
  setnames(rmap, c('otherEndChr','oestart','oeend','otherEndID'))
  setnames(baitmap, c('baitChr','bstart','bend','baitID','bait_name'))
  print(head(id_space))
  
  # Merging enhancers and genes 
  if("DHS.RPKM.quantile" %in% names(candidate_enhancers)) {
      ce_test <- candidate_enhancers[, .(chr, start, end, DHS.RPKM.quantile, H3K27ac.RPKM.quantile, cellType, class, genicSymbol, name, activity_base)]
  } else if("ATAC.RPKM.quantile" %in% names(candidate_enhancers)) {
      ce_test <- candidate_enhancers[, .(chr, start, end, ATAC.RPKM.quantile, H3K27ac.RPKM.quantile, cellType, class, genicSymbol, name, activity_base)]
  } else {
      stop("Neither DHS.RPKM.quantile nor ATAC.RPKM.quantile found in candidate_enhancers")
  } 
  cg_test <- candidate_genes[,.(chr,start,end,name,score,strand,symbol,tss,Expression,is_ue,cellType,PromoterActivityQuantile)]
  setkey(ce_test,chr)
  setkey(cg_test,chr)
  
  cand_pairs = merge(cg_test,ce_test,by='chr',allow.cartesian = TRUE)
  cand_pairs[,distSign := abs(0.5*(start.y+end.y)-tss)]
  rm(ce_test)
  rm(cg_test)
  
  # Filtering based on maxWindow and minWindow
  cand_pairs = cand_pairs[distSign <= maxWindow & distSign > minWindow]
  print(paste("Dimension of enh-gene map before overlap",dim(cand_pairs)))
  
  # Overlapping ID with gene coordinate
  print('Overlapping ID with gene coordinates')
  cand_pairs[,chr:= as.character(chr)]
  cand_pairs[,tss1:=tss+1]
  setkey(cand_pairs,chr,tss,tss1)
  setnames(id_space,  c('baitChr','bstart','bend','baitID','bname'))
  setkey(id_space,baitChr,bstart,bend,baitID)
  cand_pairs <- foverlaps(cand_pairs,id_space,by.x = c('chr','tss','tss1'),by.y=c('baitChr','bstart','bend'), type = 'any')
  print(paste("Nrows of overlap of genes with id space",dim(cand_pairs)[1]))
  print(head(cand_pairs))
  
  # Overlapping ID with enh coordinates
  print('Overlapping ID with enh coordinates')
  setkey(cand_pairs,chr,start.y,end.y)
  id_space = id_space[,c('baitChr','bstart','bend','baitID','bname')]
  setnames(id_space, c('otherEndChr','oestart','oeend','otherEndID','otherEnd_name'))
  setkey(id_space,otherEndChr,oestart,oeend,otherEndID)
  cand_pairs[, midstart.y:=floor((start.y+end.y)/2)]
  cand_pairs[, midend.y:=midstart.y+1]
  ### overlaping by enhancer midpoint like in the ABC
  cand_pairs <- foverlaps(cand_pairs,id_space,by.x = c('chr','midstart.y','midend.y'),by.y=c('otherEndChr','oestart','oeend'), type = 'any')
  
  print(paste("Nrows of overlap of genes with id space",dim(cand_pairs)[1]))
  cand_pairs[,baitChr:=chr]
  cand_pairs[,otherEndChr:=chr]
  print('Overlapping done, created copies of chromosome columns.')
  cand_pairs <- cand_pairs[,c('otherEndID','otherEndChr','baitID','baitChr','start.x','end.x','start.y','end.y','strand','tss','bstart','bend','oestart','oeend','otherEnd_name','distSign')]
  cand_pairs <- unique(cand_pairs) # duplicated rows will have arisen due to diff isoforms with the same TSS
  chromosomes = unique(cand_pairs[,baitChr]) # do i need this?
  
  #very important
  baitmap = rmap
  setnames(baitmap, c('baitChr','bstart','bend','baitID'))
  print(paste('Starting import of chicago object for chr', i))
  print(colnames(baitmap))
  print(colnames(rmap))

  ##### Imputting data based on chicago object
  print('Imputting data based on chicago object')
  pchic_data <- readRDS(paste0(split_pchic,i,'.rds'))
  setDT(pchic_data)
  print(paste('Starting chromosome',i))
  print('Check the chr!')
  setkey(pchic_data,otherEndID,otherEndChr,baitID,baitChr)
  print(paste('Line 110',i))
  cand_pairs_chr = cand_pairs[baitChr == i,]
  print("Candidate pairs head is ")
  print(head(cand_pairs_chr))
  print("no Chr error")
  rm(cand_pairs)
  pchic_data[,N_imp.y:= pmax(N,Bmean)/(s_i*s_j)] #check the formula before running (note plotBaits plots unscaled N's and Bmeans, N_imp in contrast is on the normalised scale, comparable with the distance function)
  print('No N_imp error')
  print(head(cand_pairs_chr))
  cand_pairs_chr[,N_imp.x := Chicago:::.distFun(abs(distSign),distout.rds)]
  print('No Chicago error')
  cand_pairs_chr <- cand_pairs_chr[abs(distSign) > 0,]
  print('No other imp error')
  
  print(head(cand_pairs_chr))
  print(head(pchic_data))
  setkey(pchic_data,baitID,otherEndID)
  setkey(cand_pairs_chr,baitID, otherEndID)
  
  print('Imputing')
  imputed_data <- merge(cand_pairs_chr, pchic_data, by = c('otherEndID','baitID'), all.x = TRUE)
  imputed_data[, contact:=pmax(N_imp.x, N_imp.y, na.rm = T)]
  toplim<-quantile(imputed_data$contact,0.999)
  imputed_data[contact>toplim, contact:= toplim]
  #imputed_data[contact>100, contact:=100]
  print(head(imputed_data))
  
  
  ##### Exporting the results
  print('Exporting the results')
  imputed_data[,tss1 := tss-1]
  imputed_data[,tss2 := tss+1]
  out = imputed_data[,c('baitChr.x','tss1','tss2','otherEndChr.x','start.y','end.y','otherEnd_name','contact','score')] # todo: see if we want to retain any other columns here
  out[,baitChr.x := paste0('chr',baitChr.x)]
  out[,otherEndChr.x := paste0('chr',otherEndChr.x)]
  setnames(out, c("bait_chr", "bait_start", "bait_end","otherEnd_chr", "otherEnd_start", "otherEnd_end", "otherEnd_name", 'contact','CG_score'))
  print("here is the imputed data head")
  print(head(out))
  system(paste0('mkdir ',pchic_out_prefix,'chr',i))
  fwrite(out, paste0(pchic_out_prefix,'chr',i,'/CG_score_chr',i,'.bedpe'), sep="\t", quote=FALSE, row.names=FALSE,col.names = FALSE)
  fwrite(out[,1:8], paste0(pchic_out_prefix,'chr',i,'/chr',i,'.bedpe'), sep="\t", quote=FALSE, row.names=FALSE,col.names = FALSE)
  system(paste0('gzip ',pchic_out_prefix,'chr',i,'/chr',i,'.bedpe'))
  print(paste('Finished chromosome',i))
  rm(out)
  rm(imputed_data)
  system(paste0('rm -r ',input_cand_dir,'EnhancerList_chr',i,'.txt'))
  system(paste0('rm -r ',input_cand_dir,'GeneList_chr',i,'.txt'))
  system(paste0('rm -r ',split_pchic,i,'.rds'))
}




