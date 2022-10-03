library(data.table)
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {

print('Arguements supplied succesfully')

path = args[1]
file = fread(path)  
print(paste('There are', dim(file)[1],'regions'))
chromosomes <- unique(file[,chr])
write_path = strsplit(path,'/')[[1]]


for(i in chromosomes){
  print(paste('Exporting', i,'chromosomes.'))
  write_path_i = write_path
  write_path_i[length(write_path_i)] = paste0(strsplit(write_path_i[length(write_path_i)],".txt")[[1]],'_',i,'.txt')
  write_path_i = paste0(write_path_i, collapse='/')
  print(paste('New path is ', write_path_i))
  fwrite(file[chr == i,],write_path_i,sep='\t')
}
}
