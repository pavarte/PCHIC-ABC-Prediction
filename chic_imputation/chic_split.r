library(dplyr)
library(stringr)
library(data.table)
library(Chicago)
args = commandArgs(trailingOnly=TRUE)
pchic <- args[1]
testDesignDir <- file.path(args[2])
bmap <- args[3]
rmap <- args[4]
dist_out <- args[5]
pchic_out <- args[6]
pchic <- readRDS(pchic)
pchic <- modifySettings(pchic, designDir=testDesignDir)
print('The distance parameters are:')
print(pchic@params$distFunParams)
print(paste0('Saving them at ',dist_out))
#saveRDS(pchic@params$distFunParams,dist_out)
print('Distance parameters saved.')

print(paste('The baitmap is',bmap))
bmap <- fread(bmap)
merged <- merge(pchic@x,bmap, by.y = 'V4', by.x = 'baitID', all.x = TRUE)
setnames(merged,'V1','baitChr')
print(paste('The rmap is',rmap))
rmap <- fread(rmap)
merged2 <- merge(merged,rmap, by.y = 'V4', by.x = 'otherEndID', all.x = TRUE)
setnames(merged2,'V1','otherEndChr')
merged2 <- merged2[baitChr == otherEndChr]
chromosomes <- unique(merged$baitChr)
for (i in chromosomes){
  small <- merged2[baitChr == i]
  print(paste('Saving pchic chromosome',i))
  saveRDS(small, paste0(pchic_out,i,'.rds'))
  print('Moving onto next chromosome.')
}



