library(dplyr)
library(stringr)
library(data.table)
##### Creating candidate pairs from Candidate Enhancers and Candidate Genes
args = commandArgs(trailingOnly=TRUE)

chromosomes = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")
maxWindow = args[1]
minWindow = args[2]

##### Reading params and design dir set up 
#params <- readRDS("~/ilc3_0922/contact/raw/dist_params/ILC3_5K_dist.rds")
params = args[3]
#testDesignDir <- file.path("~/../projects/lms-spivakov-analysis/live/Design/Human_hg38_bin5K_sol_baits/")
testDesignDir = args[4]
#input_cand_dir = file.path('~/ilc3_0922/ABC_output_neigh/')
input_cand_dir = args[5]
#input_pchic = file.path('~/ilc3_0922/contact/raw/5K/ILC3_5K_raw_')
input_pchic = args[6]
#output_dir = file.path('~/ilc3_0922/contact/imputed/5K_corr/')
output_dir = args[6]
baitmap=args[7] #human_DpnII_5K_sol_baits.baitmap
rmap=args[8] #"human_DpnII_5K_sol_baits.rmap
##### Starting imputation for loop 
for (i in chromosomes){
  candidate_enhancers <- fread(paste0(input_cand_dir,'EnhancerList_chr',i,'.txt'))
  candidate_genes <- fread(paste0(input_cand_dir,'GeneList_chr',i,'.txt'))
  setDT(candidate_enhancers)
  setDT(candidate_genes)
  candidate_enhancers[,chr:=str_split(chr,'chr', simplify = TRUE)[,2]]
  candidate_genes[,chr:=str_split(chr,'chr', simplify = TRUE)[,2]]
  candidate_enhancers[,chr:=as.character(chr)]
  candidate_genes[,chr:=as.character(chr)]
  
  
  
  #### Reading baitmap and rmap
  baitmap <- fread(paste0(testDesignDir,baitmap),sep ='\t',header=FALSE)
  rmap <- fread(paste0(testDesignDir,rmap),sep ='\t',header=FALSE)
  id_space = merge(baitmap,rmap,by = c('V1','V2','V3','V4'), all.y = TRUE)
  id_space[is.na(V5),V5 := '.']
  id_space[,V1:= as.character(V1)]
  setnames(rmap, c('otherEndChr','oestart','oeend','otherEndID'))
  setnames(baitmap, c('baitChr','bstart','bend','baitID','bait_name'))
  
  # Merging enhancers and genes 
  ce_test <- candidate_enhancers[,.(chr,start,end,ATAC.RPKM.quantile,H3K27ac.RPKM.quantile,cellType,class,genicSymbol,name,activity_base)]
  cg_test <- candidate_genes[,.(chr,start,end,name,score,strand,symbol,tss,Expression,Expression.quantile,is_ue,cellType,PromoterActivityQuantile)]
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
  cand_pairs[,chr:= as.character(chr)]
  cand_pairs[,tss1:=tss+1]
  setkey(cand_pairs,chr,tss,tss1)
  setnames(id_space,  c('baitChr','bstart','bend','baitID','bname'))
  setkey(id_space,baitChr,bstart,bend,baitID)
  cand_pairs <- foverlaps(cand_pairs,id_space,by.x = c('chr','tss','tss1'),by.y=c('baitChr','bstart','bend'), type = 'any')
  print(paste("Nrows of overlap of genes with id space",dim(cand_pairs)[1]))
  
  
  # Overlapping ID with enh coordinates
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
  pchic <- readRDS(paste0(input_pchic,i,'.rds'))
  setDT(pchic)
  print(paste('Starting chromosome',i))
  pchic_data <- pchic#@x remove hash if it's the whole file  
  rm(pchic)
  setkey(pchic_data,otherEndID,otherEndChr,baitID,baitChr)
  cand_pairs_chr = cand_pairs[baitChr == i]
  rm(cand_pairs)
  pchic_data[,N_imp.y:= pmax(N,Bmean)/(s_i*s_j)] #check the formula before running (note plotBaits plots unscaled N's and Bmeans, N_imp in contrast is on the normalised scale, comparable with the distance function)
  cand_pairs_chr[,N_imp.x := Chicago:::.distFun(abs(distSign),params)]
  cand_pairs_chr <- cand_pairs_chr[abs(distSign) > 0,]
  setkey(pchic_data,baitID,otherEndID)
  setkey(cand_pairs_chr,baitID, otherEndID)
  
  imputed_data <- merge(cand_pairs_chr,pchic_data,by = c('otherEndID','baitID'),all.x = TRUE)
  imputed_data[, contact:=pmax(N_imp.x, N_imp.y, na.rm = T)]
  toplim<-quantile(imputed_data$contact,0.999)
  imputed_data[contact>toplim, contact:= toplim]
  #imputed_data[contact>100, contact:=100]
  print(head(imputed_data))
  
  
  ##### Exporting the results
  imputed_data[,tss1 := tss-1]
  imputed_data[,tss2 := tss+1]
  out = imputed_data[,c('baitChr.x','tss1','tss2','otherEndChr.x','start.y','end.y','otherEnd_name','contact','score')] # todo: see if we want to retain any other columns here
  out[,baitChr.x := paste0('chr',baitChr.x)]
  out[,otherEndChr.x := paste0('chr',otherEndChr.x)]
  setnames(out, c("bait_chr", "bait_start", "bait_end","otherEnd_chr", "otherEnd_start", "otherEnd_end", "otherEnd_name", 'contact','CG_score'))
  print("here is the imputed data head")
  print(head(out))
  fwrite(out, paste0(output_dir,'chr',i,'/CG_score_chr',i,'.bedpe'), sep="\t", quote=FALSE, row.names=FALSE,col.names = FALSE)
  fwrite(out[,1:8], paste0(output_dir,'chr',i,'/chr',i,'.bedpe'), sep="\t", quote=FALSE, row.names=FALSE,col.names = FALSE)
  system(paste0('gzip ',output_dir,'chr',i,'/chr',i,'.bedpe'))
  print(paste('Finished chromosome',i))
  rm(out)
  rm(imputed_data)
}





