# My R Script
# Author: Pavel Artemov
# Date: 2024-10
# Description: This script divides pchic rds into subfiles and also exports the distance file for imputation

# Usage: Rscript imputation_script.R ${pchic.rds} ${design} ${enhancerdir/} ${split_pchic} ${distout.rds} ${imputed_pchic_prefix}
# arg1 is pchic.rds, RDS file from step 2 of running 
# arg2 is testDesignDir, needs to end in / 
# arq3 is input_cand_dir, needs to end in /
# arg4 is split_pchic, filename prefix only
# arg5 is RDS file, distance parameters
# arg6 is output directory, needs to end in /  

# Below are some examples of arguments for input 
# args = vector("character")
# args[1] = "/rds/general/project/lms-spivakov-analysis/live/artemov/ABC_MASTER/ABC_data/input/K562_roadmap_Ery_MS/Ery_Step2_Merged.rds"
# args[2] = "/rds/general/project/lms-spivakov-analysis/live/Design/Human_hg19/"
# args[3] = "/rds/general/project/lms-spivakov-analysis/live/artemov/ABC_MASTER/ABC_data/input/K562_roadmap_Ery_MS/"
# args[4] = "/rds/general/project/lms-spivakov-analysis/live/artemov/ABC_MASTER/ABC_data/input/K562_roadmap_Ery_MS/K562_roadmap_Ery_MS_"
# args[5] = "~/analysis/artemov/ABC_MASTER/ABC_data/input/K562_roadmap_Ery_MS/K562_roadmap_Ery_MS_dist.rds"
# args[6] = "/rds/general/project/lms-spivakov-analysis/live/artemov/ABC_MASTER/ABC_data/input/K562_roadmap_Ery_MS/imputed_contact/"
# args[7] = 1,2 – number of fragment steps to be used for imputation

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=7) {
  stop("imputation_script.R: Wrong number of arguments", call.=FALSE) 
} else {
  print('Arguments supplied succesfully to imputation script')
}


library(dplyr)
library(stringr)
library(data.table)
library(Chicago)
library(ggplot2)
library(ggpubr)

##### Input some arguements for imputation
maxWindow = 5000000
minWindow = 0
pchic.rds = args[1]
testDesignDir = args[2]
input_cand_dir = args[3]
split_pchic = args[4]
distout.rds = args[5]
pchic_out_prefix = file.path(args[6])

##### Splitting Enhancer and Gene Lists by chromosome
chr_split_command <- paste0("Rscript chr_split.r ",input_cand_dir)
# Execute the command
if(!system(chr_split_command)){
  message("Successfully finished chromosome split")
}else{
  message("Error in the chromosomes split")
}
##### Splitting PCHiC by chromosome
chic_split_command <- paste("Rscript chic_split.r",pchic.rds, testDesignDir, distout.rds, split_pchic)
# Execute the command
if(!system(chic_split_command)){
  message("Succesfully finished PCHiC split")
}else{
  message("Error in the PCHiC split")
}

enhancers=fread(paste0(input_cand_dir,'/EnhancerList.txt'))
enhancers=enhancers[,chr:=str_split(chr,'chr', simplify = TRUE)[,2]]
chromosomes=unique(enhancers[,chr])
print("Using chromosomes in enhancers:")
print(chromosomes)
rm(enhancers)
genes=fread(paste0(input_cand_dir,'/GeneList.txt'))
genes=genes[,chr:=str_split(chr,'chr', simplify = TRUE)[,2]]
chromosomes=unique(genes[,chr])
print("Using chromosomes in genes:")
print(chromosomes)
rm(genes)
chromosomes=chromosomes[chromosomes!='Y']
#reading in the distance parameters
print("Reading distance file")
print(args[5])
distout.rds = readRDS(args[5])
# preparing folders for the input
#chromosomes=1
imputed_data_all <- list()

for (i in chromosomes){
message(paste0("Imputing for chromosome",i))
baitmap=list.files(path = testDesignDir, pattern = "\\.baitmap$")[1]
rmap=list.files(path = testDesignDir, pattern = "\\.rmap$")[1]
baitmap <- fread(paste0(testDesignDir,baitmap),sep ='\t',header=FALSE)
rmap <- fread(paste0(testDesignDir,rmap),sep ='\t',header=FALSE)
id_space = merge(baitmap,rmap,by = c('V1','V2','V3','V4'), all.y = TRUE)
id_space[is.na(V5),V5 := '.']
id_space[,V1:= as.character(V1)]
setnames(rmap, c('otherEndChr','oestart','oeend','otherEndID'))
setnames(baitmap, c('baitChr','bstart','bend','baitID','bait_name'))
print("Resulting ID Space head")
print(head(id_space))
print(unique(id_space[,V1]))
print("Unique chromosomes in rmap")
print(unique(rmap[,otherEndChr]))
print("Unique chromosomes in baitmap")
print(unique(baitmap[,baitChr]))

bmap=baitmap
candidate_enhancers <- fread(paste0(input_cand_dir,'EnhancerList_chr',i,'.txt'))
candidate_genes <- fread(paste0(input_cand_dir,'GeneList_chr',i,'.txt'))
setDT(candidate_enhancers)
setDT(candidate_genes)
candidate_enhancers[,chr:=str_split(chr,'chr', simplify = TRUE)[,2]]
candidate_genes[,chr:=str_split(chr,'chr', simplify = TRUE)[,2]]
print("Candidate enhancers chromosomes")
print(unique(candidate_enhancers[,chr]))
print("Candidate genes chromosomes")
print(unique(candidate_genes[,chr]))
candidate_enhancers[,chr:=as.character(chr)]
candidate_genes[,chr:=as.character(chr)] 


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
cand_pairs[,tss1:=tss]
cand_pairs[,tss:=tss]

setkey(cand_pairs,chr,tss,tss1)
setnames(id_space,  c('baitChr','bstart','bend','baitID','bname'))
setkey(id_space,baitChr,bstart,bend,baitID)
cand_pairs <- foverlaps(cand_pairs,id_space,by.x = c('chr','tss','tss1'),by.y=c('baitChr','bstart','bend'), type = 'any')
print(paste("Nrows of overlap of genes with id space",dim(cand_pairs)[1]))

# Overlapping ID with enh coordinates
print('Overlapping ID with enh coordinates')
setkey(cand_pairs,chr,start.y,end.y)
id_space = id_space[,c('baitChr','bstart','bend','baitID','bname')]
setnames(id_space, c('otherEndChr','oestart','oeend','otherEndID','otherEnd_name'))
setkey(id_space,otherEndChr,oestart,oeend,otherEndID)

cand_pairs <- foverlaps(cand_pairs,id_space,by.x = c('chr','start.y','end.y'),by.y=c('otherEndChr','oestart','oeend'), type = 'any') ## MS cand_pairs examples
cand_pairs <- cand_pairs[,c('otherEndID','baitID','start.x','end.x','start.y','end.y','strand','tss','bstart','bend','oestart','oeend','otherEnd_name','distSign','chr', 'name.y') ]
cand_pairs[,baitChr:=chr]
cand_pairs[,otherEndChr:=chr]

baitmap = rmap
setnames(baitmap, c('baitChr','bstart','bend','baitID'))
print(paste('Starting import of chicago object for chr', i))


##### Imputting data based on chicago object
message('Imputing data based on chicago object')
print('Imputing data based on chicago object')
pchic_data <- readRDS(paste0(split_pchic,i,'.rds'))
setDT(pchic_data)
print(dim(pchic_data))
print("PCHiC data chromosomes")
print(unique(pchic_data[,baitChr]))
print(unique(pchic_data[,otherEndChr]))
message(paste('Starting chromosome',i))
print('Check the chr!')
setkey(pchic_data,otherEndID,otherEndChr,baitID,baitChr)
cand_pairs_chr = cand_pairs[baitChr == i] # note this is the cand_pairs_ms approach
print("no Chr error")
#rm(cand_pairs) to save memory
# Note that Bmean/(s_i*s_j) should be equivalent to how N_imp.x is defined below
pchic_data[,N_imp.y:= pmax(N,Bmean)/(s_i*s_j)] #check the formula before running (note plotBaits plots unscaled N's and Bmeans, N_imp in contrast is on the normalised scale, comparable with the distance function)

directory_path <- input_cand_dir
pchic_data <- pchic_data %>%
  mutate(group = case_when(
    N > Bmean ~ "Observed, N > Bmean",
    Bmean >= N ~ "Observed, Bmean >= N"
  ))
pchic_data_mut = pchic_data
pchic_data_mut = pchic_data_mut[abs(pchic_data_mut$distSign) < 200000]
pchic_data_mut$distSign_bin <- cut(abs(pchic_data_mut$distSign), breaks = 10)

gg0=ggplot(pchic_data_mut, aes(x = distSign_bin, y = asinh(N_imp.y), fill = group)) +
  geom_boxplot() +
  labs(title = paste("Imputed contact distribution against distance from TSS chr",i),
       x = "Distance bin to TSS",
       y = "Contact Frequency (asinh scale)",
       fill = "Type of PCHiC contact") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(gg0, filename = paste0(directory_path, "/observed.only.distr_N_v_Bmean_for_chr", i, ".pdf"), width = 16, height = 8)


# INCORPORATE INTO THE CODE ##
message("Percent contacts imputed from Bmean: ", round(nrow(pchic_data[N_imp.y==Bmean/(s_i*s_j)])/nrow(pchic_data)*100,2))

print('No N_imp error')
cand_pairs_chr[,N_imp.x := Chicago:::.distFun(abs(distSign),distout.rds)]
print('No Chicago error')
cand_pairs_chr <- cand_pairs_chr[abs(distSign) > 0,]
print('No other imp error')

setkey(pchic_data,baitID,otherEndID)
setkey(cand_pairs_chr,baitID, otherEndID)
print('Imputing')
imputed_data1 <- merge(cand_pairs_chr, pchic_data, by = c('otherEndID','baitID'), all.x = TRUE,allow.cartesian = TRUE)
imputed_data1[!is.na(N_imp.y)]
##swap pchic name and then do again 
setnames(pchic_data, old = c("otherEndID", "V2.y", "V3.y"),
         new = c("temp_baitID", "temp_V2.x", "temp_V3.x"))
setnames(pchic_data, old = c("baitID", "V2.x", "V3.x"),
         new = c("otherEndID", "V2.y", "V3.y"))
setnames(pchic_data, old = c("temp_baitID", "temp_V2.x", "temp_V3.x"),
         new = c("baitID", "V2.x", "V3.x"))
## pchic data at this point has flipped IDs of bait and other end so it needs to be flipped back if you decide to use it downstream
imputed_data <- merge(imputed_data1, pchic_data, by = c('otherEndID','baitID'), all.x = TRUE,allow.cartesian = TRUE)
imputed_data[!is.na(N_imp.y.y)]

imputed_data[is.na(N_imp.x),N_imp.x:=0]
imputed_data[is.na(N_imp.y.x),N_imp.y.x:=0]
imputed_data[is.na(N_imp.y.y),N_imp.y.y:=0]
imputed_data[, contact:=pmax(N_imp.x, N_imp.y.x,N_imp.y.y, na.rm = T)]
imputed_data
# Changing the algorithm! Not using the midpoint for mapping to the restriction fragment but the fragment with the highest contact.
# super important pay attention here 
imputed_data_collapsed = imputed_data[, {sel = which(contact==max(contact))[1]; 
.(baitID=baitID[sel], baitChr.x=baitChr.x[sel], baitChr.y=baitChr.y[sel],
  otherEndID=otherEndID[sel], otherEndChr.x=otherEndChr.x[sel],otherEndChr.y=otherEndChr.y[sel], otherEnd_name=otherEnd_name[sel],
  V2.x.y=V2.x.y[sel], V3.x.y=V3.y.y[sel],V2.x.x=V2.x.x[sel], V3.x.x=V3.x.x[sel],
  V2.y.y=V2.y.y[sel], V3.y.y=V3.y.y[sel],V2.y.x=V2.y.x[sel], V3.y.x=V3.y.x[sel],
  s_j.x=s_j.x[sel], s_j.y=s_j.y[sel], 
  start.y=start.y[sel], end.y=end.y[sel], group.x=group.x[sel],group.y=group.y[sel],score.x=score.x[sel], score.y=score.y[sel],distbin.x=distbin.x[sel], distSign.x=distSign.x[sel],
  contact=max(contact))}, by=c("name.y", "tss")]

imputed_data0 <- imputed_data ## saving just in case 
imputed_data <- imputed_data_collapsed
real_baitmap = fread(paste0(testDesignDir,list.files(path = testDesignDir, pattern = "\\.baitmap$")[1]),sep ='\t',header=FALSE)

setnames(pchic_data, old = c("otherEndID", "V2.y", "V3.y"),
         new = c("temp_baitID", "temp_V2.x", "temp_V3.x"))
setnames(pchic_data, old = c("baitID", "V2.x", "V3.x"),
         new = c("otherEndID", "V2.y", "V3.y"))
setnames(pchic_data, old = c("temp_baitID", "temp_V2.x", "temp_V3.x"),
         new = c("baitID", "V2.x", "V3.x"))
message("Imputed based on distance: ", round(nrow(imputed_data[is.na(s_j.x)&is.na(s_j.y)])/nrow(imputed_data)*100,2))
message("Of the imputed based on distance, have promoter-containing fragments not in baitmap (baitmap/capture design issue or candidate pairs issue): ", round(nrow(imputed_data[is.na(s_j.x)&is.na(s_j.y) & !baitID%in%real_baitmap$V4])/nrow(imputed_data[is.na(s_j.x)&is.na(s_j.y)])*100,2))
message("Of the imputed based on distance with promoter-containing fragments in baitmap, have baits that are not in the chicago object (coverage issue): ", 
        round(nrow(imputed_data[is.na(s_j.x)&is.na(s_j.y) & baitID%in%real_baitmap$V4 & !baitID%in%pchic_data$baitID])/nrow(imputed_data[is.na(s_j.x)&is.na(s_j.y) & baitID%in%real_baitmap$V4])*100,2))
message("Of the imputed based on distance with promoter-containing fragments in the chicago object, have enhancer-containing fragments that are not in the chicago object (coverage issue): ", 
        round(nrow(imputed_data[is.na(s_j.x)&is.na(s_j.y) & baitID%in%pchic_data$baitID & !otherEndID %in%pchic_data$otherEndID])/nrow(imputed_data[is.na(s_j.x)&is.na(s_j.y) & baitID%in%pchic_data$baitID])*100,2))
### 
pchic_data[, id:=paste(baitID,otherEndID, sep="_")]
pchic_data[, reverseid:=paste(otherEndID,baitID, sep="_")]
imputed_data[, id:=paste(baitID,otherEndID, sep="_")] ##potentially need to flip it 

message("Of the imputed based on distance with promoter- and enhancer-containing fragments in the chicago object, have the respective bait-otherEnd contacts not in the chicago object (coverage issue): ", 
        round(nrow(imputed_data[is.na(s_j.x)&is.na(s_j.y) & baitID%in%pchic_data$baitID & otherEndID %in%pchic_data$otherEndID & !id%in%pchic_data$id ])/nrow(imputed_data[is.na(s_j.x)&is.na(s_j.y) & baitID%in%pchic_data$baitID & otherEndID%in% pchic_data$otherEndID])*100,2))

message("Total: ", nrow(imputed_data))
message("Total_imputed_based_on_distance(IBOD): ", nrow(imputed_data[is.na(s_j.x)&is.na(s_j.y)]))
message("IBOD_promoter_not_in_baitmap:",nrow(imputed_data[is.na(s_j.x)&is.na(s_j.y) & !baitID%in%real_baitmap$V4]))
message("IBOD_promoter_in_baitmap_but_not_in_chicago_object:",nrow(imputed_data[is.na(s_j.x)&is.na(s_j.y) & baitID%in%real_baitmap$V4 & !baitID%in%pchic_data$baitID]))
message("IBOD_promoter_in_chicago_object_enhancer_not_in_chicago_object:", nrow(imputed_data[is.na(s_j.x)&is.na(s_j.y) & baitID%in%pchic_data$baitID & !otherEndID %in% pchic_data$otherEndID]))
message("IBOD_contact_not_in_chicago_object: ",nrow(imputed_data[is.na(s_j.x)&is.na(s_j.y) & baitID%in%pchic_data$baitID & otherEndID %in%pchic_data$otherEndID & !id%in%c(pchic_data$id, pchic_data$reverseid) ]))
# according to MS, even though for these cases immediately above we have access to Bmean, s_i and s_j, there's no point looking them up to compute N_imp.y and doing it from the distance function directly as 
# it is done for N_imp.x should produce the same result
## plot for distance profile
## ask the opposite question - which ones are in PCHiC data not in the other one 
column_names <- names(imputed_data)
repeated_names <- substr(column_names, 1, nchar(column_names) - 2)
prefixes <- unique(repeated_names[duplicated(repeated_names)])

print(column_names)
print(prefixes)
##
#prefixes <- unique(gsub("\\.(x|y)$", "", columns_with_x_y))

for (col in prefixes) {
  col_x <- paste0(col, ".x")
  col_y <- paste0(col, ".y")

  # Check if both .x and .y versions exist in the data.table
  if (col_x %in% names(imputed_data) & col_y %in% names(imputed_data)) {
    print(col)

    # Get column types for both .x and .y columns
    type_x <- typeof(imputed_data[[col_x]])
    type_y <- typeof(imputed_data[[col_y]])

    # Coalesce based on the type of columns
    if (type_x == "double" & type_y == "double") {
      # Both are numeric
      imputed_data[, (col) := fcoalesce(as.double(get(col_x)), as.double(get(col_y)))]
    } else if (type_x == "character" & type_y == "character") {
      # Both are strings
      imputed_data[, (col) := fcoalesce(as.character(get(col_x)), as.character(get(col_y)))]
    } else {
      # Handle mixed types
      imputed_data[, (col) := fcoalesce(as.character(get(col_x)), as.character(get(col_y)))]
    }

    # Optionally, remove the original .x and .y columns
      imputed_data[, c(col_x, col_y) := NULL]
  }
}
#print("Saving imputed data for tests")
#fwrite(imputed_data,"/home/pa2915/analysis/artemov/ABC_MASTER/scripts/imputed_data_for_test.csv",sep='\t')
##### Some modifications for downstream export
imputed_data[,tss1 := tss]
imputed_data[,tss2 := tss]
imputed_data[,distSign:=0.5*(start.y+end.y)-tss]
imputed_data[,enh_midpoint:=0.5*(start.y+end.y)]
imputed_data[,start.y:=floor(enh_midpoint)]
imputed_data[,end.y:=floor(enh_midpoint)+1]
#imputed_data[,start.y:=floor(enh_midpoint)-1]
#imputed_data[,start.y:=floor(enh_midpoint)+1]
##### Capping N based on nearest fragment and other end
rmap=list.files(path = testDesignDir, pattern = "\\.rmap$")[1]
rmap <- fread(paste0(testDesignDir,rmap),sep ='\t',header=FALSE)
rmap=rmap[V1==imputed_data$baitChr[1]]
nuber_of_steps = as.numeric(args[7])
rmap[,V2:=as.numeric(V2)]
rmap[,V3:=as.numeric(V3)]
imputed_data[,tss:=as.numeric(tss)]
imputed_data[,distSign:=as.numeric(distSign)]
setnames(rmap,c("baitChr","next_frag_boundary","next_frag_outer_boundary","enh_ID_plus"))
rmap[,frag_length:=abs(next_frag_outer_boundary-next_frag_boundary)]
if(median(rmap[,frag_length])>1500){
#median_frag_length=5210
median_frag_length=1500

}else{
median_frag_length=1500
}


## Potentially change one selection just to be based on one of the boundaries 
message(paste0("Number of contacts, N capped because distance of enhancer is closer than next and prev fragment from TSS: ",
       nrow(imputed_data[abs(distSign) <= median_frag_length,])))
imputed_data[!is.na(group), group := "Contact frequency observed"]
imputed_data[is.na(group), group := "Imputed (using CHiCAGO distance function)"]
imputed_data[abs(distSign) <= median_frag_length & contact>Chicago:::.distFun(median_frag_length,distout.rds), group:='Capped']
imputed_data[abs(distSign) <= median_frag_length & contact>Chicago:::.distFun(median_frag_length,distout.rds), contact:=Chicago:::.distFun(median_frag_length,distout.rds)]


##### Plots for diagnostics
print('Creating plots for the distance and imputed contact distribution')

imputed_data_mut = imputed_data[abs(imputed_data$distSign) <200000]
imputed_data_mut$distSign_bin <- cut(abs(imputed_data_mut$distSign), breaks = 10)
print(unique(imputed_data_mut[,group]))
imputed_data_mut[is.na(group), group := "Imputed (using CHiCAGO distance function)"]
gg01=ggplot(imputed_data_mut, aes(x = distSign_bin, y = asinh(contact), fill = group)) +
  geom_boxplot() +
  labs(title = paste("Imputed contact distribution against distance from TSS chr",i),
       x = "Distance bin to TSS",
       y = "Contact Frequency (asinh scale)",
       fill = "Type of PCHiC contact") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(gg01, filename = paste0(directory_path, "/distr_imputed.data_N_v_Bmean_for_chr", i, ".pdf"), width = 16, height = 8)

imputed_data_mut = imputed_data[abs(imputed_data$distSign) <25000]
imputed_data_mut$distSign_bin <- cut(abs(imputed_data_mut$distSign), breaks = 10)
print(unique(imputed_data_mut[,group]))
imputed_data_mut[is.na(group), group := "Imputed (using CHiCAGO distance function)"]

gg02=ggplot(imputed_data_mut, aes(x = distSign_bin, y = asinh(contact), fill = group)) +
  geom_boxplot() +
  labs(title = paste("Imputed contact distribution against distance from TSS chr",i),
       x = "Distance bin to TSS",
       y = "Contact Frequency (asinh scale)",
       fill = "Type of PCHiC contact") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(gg02, filename = paste0(directory_path, "/distr_imputed.data_v.close.range_N_v_Bmean_for_chr", i, ".pdf"), width = 16, height = 8)

imputed_data_mut = imputed_data[abs(imputed_data$distSign) <5000]
imputed_data_mut$distSign_bin <- cut(abs(imputed_data_mut$distSign), breaks = 11)
print(unique(imputed_data_mut[,group]))
imputed_data_mut[is.na(group), group := "Imputed (using CHiCAGO distance function)"]

gg02=ggplot(imputed_data_mut, aes(x = distSign_bin, y = asinh(contact), fill = group)) +
  geom_boxplot() +
  labs(title = paste("Imputed contact distribution against distance from TSS chr",i),
       x = "Distance bin to TSS",
       y = "Contact Frequency (asinh scale)",
       fill = "Type of PCHiC contact") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(gg02, filename = paste0(directory_path, "/distr_imputed.data_v.v.close.range_N_v_Bmean_for_chr", i, ".pdf"), width = 16, height = 8)


##
data_for_plot <- rbind(
  data.frame(abs_distSign = abs(imputed_data[!is.na(imputed_data$s_j), "distSign"]),
             group = "Contact inferred from observed contact frequencies (total)"),
  data.frame(abs_distSign = abs(imputed_data[is.na(imputed_data$s_j), "distSign"]),
             group = "Imputed using CHiCAGO Distance Function (not observed in PCHiC)")
)
##
gg1 <- ggplot(data_for_plot, aes(x = distSign, fill = group)) +
  geom_histogram(position = "identity", alpha = 0.25, bins = 30) +
  labs(title = paste("Distance distributions of imputed contacts",i),
       x = "abs(distance in bp)",
       y = "Count") +
  theme_minimal() +
  scale_fill_manual(values = c("red", "blue","grey","green")) +
  theme(legend.title = element_blank())

gg2 <- ggplot(data_for_plot, aes(x = distSign, fill = group)) +
  geom_histogram(aes(y = after_stat(density)), position = "identity", alpha = 0.25, bins = 30) +
  labs(title = paste("Distance distributions of imputed contacts chr",i),
       x = "abs(distance in bp)",
       y = "Density") +
  theme_minimal() +
  scale_fill_manual(values = c("red", "blue","grey","green")) +
  theme(legend.title = element_blank())
gg3 <- ggarrange(gg1,gg2)
##
imputed_data_mut <- imputed_data %>%
  mutate(group = case_when(
    !is.na(s_j) ~ "Observed contact frequencies",
    is.na(s_j) ~ "Imputed using CHiCAGO Distance Function (not observed in PCHiC)"
  ))

# Bin the distSign values (you can adjust the number of bins or the bin width)
imputed_data_mut$distSign_bin <- cut(abs(imputed_data_mut$distSign), breaks = 10)

gg4 <- ggplot(imputed_data_mut, aes(x = distSign_bin, y = asinh(contact), fill = group)) +
  geom_boxplot() +
  labs(title = paste("Contact frequency distribution against distance from TSS for chr",i),
       x = "abs(distance in bp)",
       y = "Contact Frequency (asinh scale)",
       fill = "Group") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotating x-axis labels for better readability

imputed_data_mut <- imputed_data_mut[abs(imputed_data_mut$distSign)<1000000,]
imputed_data_mut$distSign_bin <- cut(abs(imputed_data_mut$distSign), breaks = 10)

gg5 <- ggplot(imputed_data_mut, aes(x = distSign_bin, y = asinh(contact), fill = group)) +
  geom_boxplot() +
  labs(title = paste("Contact frequency distribution against distance from TSS for chr",i),
       x = "abs(distance in bp)",
       y = "Contact Frequency (asinh scale)",
       fill = "Group") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotating x-axis labels for better readability
imputed_data_mut <- imputed_data_mut[abs(imputed_data_mut$distSign)<100000,]
imputed_data_mut$distSign_bin <- cut(abs(imputed_data_mut$distSign), breaks = 10)
gg6 <- ggplot(imputed_data_mut, aes(x = distSign_bin, y = asinh(contact), fill = group)) +
  geom_boxplot() +
  labs(title = paste("Contact frequency distribution against distance from TSS for chr",i),
       x = "abs(distance in bp)",
       y = "Contact Frequency (asinh scale)",
       fill = "Group") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
# Define the directory path
directory_path <- input_cand_dir

# Now save the plots to the created directory
ggsave(gg3, filename = paste0(directory_path, "/distance_distr_imputed_contacts_for_chr", i, ".pdf"), width = 16, height = 8)
ggsave(gg4, filename = paste0(directory_path, "/count_distr_imputed_contacts_for_chr", i, ".pdf"), width = 16, height = 8)
ggsave(gg5, filename = paste0(directory_path, "/count_distr_close_range_imputed_contacts_for_chr", i, ".pdf"), width = 16, height = 8)
ggsave(gg6, filename = paste0(directory_path, "/count_distr_v_close_range_imputed_contacts_for_chr", i, ".pdf"), width = 16, height = 8)

print("Chromosome names of bait chr of exported data")
print(unique(imputed_data[,baitChr]))
print("Chromosome names of other end chr of exported data")
print(unique(imputed_data[,otherEndChr]))
### Actual export and saving
out = imputed_data[,c('baitChr','tss1','tss2','otherEndChr','start.y','end.y','otherEnd_name','contact','score','otherEndID','baitID')] # todo: see if we want to retain any other columns here
setnames(out, c("bait_chr", "bait_start", "bait_end","otherEnd_chr", "otherEnd_start", "otherEnd_end", "otherEnd_name", 'contact','CG_score','otherEndID', 'baitID'))
out[,bait_chr:=paste0('chr',bait_chr)]
out[,otherEnd_chr:=paste0('chr',otherEnd_chr)]
setkey(out,bait_chr,bait_start,bait_end,otherEnd_chr,otherEnd_start,otherEnd_end,otherEnd_name,otherEndID,baitID)
print("Saving the exported file")
system(paste0('rm -r ',pchic_out_prefix,'/chr',i))
system(paste0('mkdir ',pchic_out_prefix,'/chr',i))
fwrite(out, paste0(pchic_out_prefix,'/chr',i,'/CG_score_chr',i,'.bedpe'), sep="\t", quote=FALSE, row.names=FALSE,col.names = FALSE)
fwrite(out[,1:8], paste0(pchic_out_prefix,'/chr',i,'/chr',i,'.bedpe'), sep="\t", quote=FALSE, row.names=FALSE,col.names = FALSE)
system(paste0('gzip ',pchic_out_prefix,'/chr',i,'/chr',i,'.bedpe'))
print(paste('Finished chromosome',i))
rm(out)
# === NEW ===
imputed_data_all[[i]] <- imputed_data
print("test of additional plot")
print("testing imputed_data")
print(head(imputed_data))
print("test of imputed_data_all")
print(str(imputed_data_all))
print(head(imputed_data_all))
rm(imputed_data)
system(paste0('rm -r ',input_cand_dir,'/EnhancerList_chr',i,'.txt'))
system(paste0('rm -r ',input_cand_dir,'/GeneList_chr',i,'.txt'))
system(paste0('rm -r ',split_pchic,i,'.rds'))
}
# === GENOME-WIDE PLOTTING ===

library(ggplot2)
library(dplyr)
library(ggpubr)

print("Combining imputed data across all chromosomes")

# Combine all chromosome imputed data
imputed_data_genome <- rbindlist(imputed_data_all, use.names = TRUE, fill = TRUE)

# Assign group labels
imputed_data_genome[!is.na(s_j), group := "Observed contact frequencies"]
imputed_data_genome[is.na(s_j), group := "Imputed using CHiCAGO Distance Function (not observed in PCHiC)"]

# Convert types safely
imputed_data_genome$distSign <- as.numeric(imputed_data_genome$distSign)
imputed_data_genome$contact <- as.numeric(imputed_data_genome$contact)

# Filter and bin
imputed_data_mut_genome <- imputed_data_genome[abs(distSign) < 1000000]
imputed_data_mut_genome$distSign_bin <- cut(abs(imputed_data_mut_genome$distSign), breaks = 10)
print(head(imputed_data_mut_genome))
# Boxplot
gg_genome <- ggplot(imputed_data_mut_genome, aes(x = distSign_bin, y = asinh(contact), fill = group)) +
  geom_boxplot() +
  labs(title = "Genome-wide Contact Frequency vs Distance",
       x = "Distance bin to TSS",
       y = "Contact Frequency (asinh scale)",
       fill = "Group") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the genome-wide plot
ggsave(gg_genome, filename = paste0(input_cand_dir, "/genomewide_contact_distribution.pdf"), width = 16, height = 8)
print("Saved genome-wide contact distribution plot")

# Additional Genome-wide Diagnostic Plots

# ~~~ 1. CLOSE RANGE CONTACT FREQUENCY BOXPLOTS ~~~

# <100kb
imputed_data_mut_genome_100k <- imputed_data_genome[abs(distSign) < 100000]
imputed_data_mut_genome_100k$distSign_bin <- cut(abs(imputed_data_mut_genome_100k$distSign), breaks = 10)

gg_g100k <- ggplot(imputed_data_mut_genome_100k, aes(x = distSign_bin, y = asinh(contact), fill = group)) +
  geom_boxplot() +
  labs(title = "Genome-wide Contact Frequency (±100kb)",
       x = "Distance bin to TSS",
       y = "Contact Frequency (asinh scale)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(gg_g100k, filename = paste0(input_cand_dir, "/genomewide_contact_distribution_100kb.pdf"), width = 16, height = 8)

# <25kb
imputed_data_mut_genome_25k <- imputed_data_genome[abs(distSign) < 25000]
imputed_data_mut_genome_25k$distSign_bin <- cut(abs(imputed_data_mut_genome_25k$distSign), breaks = 10)

gg_g25k <- ggplot(imputed_data_mut_genome_25k, aes(x = distSign_bin, y = asinh(contact), fill = group)) +
  geom_boxplot() +
  labs(title = "Genome-wide Contact Frequency (±25kb)",
       x = "Distance bin to TSS",
       y = "Contact Frequency (asinh scale)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(gg_g25k, filename = paste0(input_cand_dir, "/genomewide_contact_distribution_25kb.pdf"), width = 16, height = 8)

# <5kb
imputed_data_mut_genome_5k <- imputed_data_genome[abs(distSign) < 5000]
imputed_data_mut_genome_5k$distSign_bin <- cut(abs(imputed_data_mut_genome_5k$distSign), breaks = 11)

gg_g5k <- ggplot(imputed_data_mut_genome_5k, aes(x = distSign_bin, y = asinh(contact), fill = group)) +
  geom_boxplot() +
  labs(title = "Genome-wide Contact Frequency (±5kb)",
       x = "Distance bin to TSS",
       y = "Contact Frequency (asinh scale)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(gg_g5k, filename = paste0(input_cand_dir, "/genomewide_contact_distribution_5kb.pdf"), width = 16, height = 8)

# ~~~ 2. HISTOGRAMS ~~~

data_for_plot <- rbind(
  data.frame(distSign = abs(imputed_data_genome[!is.na(imputed_data_genome$s_j), "distSign"]),
             group = "Contact inferred from observed contact frequencies (total)"),
  data.frame(distSign = abs(imputed_data_genome[is.na(imputed_data_genome$s_j), "distSign"]),
             group = "Imputed using CHiCAGO Distance Function (not observed in PCHiC)")
)

# Raw count
gg_g_hist1 <- ggplot(data_for_plot, aes(x = distSign, fill = group)) +
  geom_histogram(position = "identity", alpha = 0.25, bins = 30) +
  labs(title = "Genome-wide Imputed Contact Distance Distribution",
       x = "abs(distance in bp)", y = "Count") +
  theme_minimal() +
  theme(legend.title = element_blank())

ggsave(gg_g_hist1, filename = paste0(input_cand_dir, "/genomewide_distance_histogram_count.pdf"), width = 16, height = 8)

# Density
gg_g_hist2 <- ggplot(data_for_plot, aes(x = distSign, fill = group)) +
  geom_histogram(aes(y = after_stat(density)), position = "identity", alpha = 0.25, bins = 30) +
  labs(title = "Genome-wide Imputed Contact Distance Distribution (Density)",
       x = "abs(distance in bp)", y = "Density") +
  theme_minimal() +
  theme(legend.title = element_blank())

ggsave(gg_g_hist2, filename = paste0(input_cand_dir, "/genomewide_distance_histogram_density.pdf"), width = 16, height = 8)


