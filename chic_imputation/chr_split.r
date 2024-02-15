library(data.table)
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=1) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  
  #
  #args[1] = "/rds/general/user/pa2915/projects/lms-spivakov-analysis/live/artemov/ABC_MASTER/ABC_data/K562_roadmap/"
  print('Arguements supplied succesfully to Enhancer and Gene split script succefully')}
  #path = file.path(paste0(","EnhancerList.txt"))
  
  ## Enhancer List export
  path = file.path(paste0(args[1],"EnhancerList.txt"))
  e.file = fread(path)
  print(paste('There are', dim(e.file)[1],'regions in enhancer list'))
  chromosomes <- unique(e.file[,chr])
  write_path = path
  
  for(i in chromosomes){
    print(paste('Exporting genes', i,'chromosomes.'))
    write_path_i = write_path
    write_path_i[length(write_path_i)] = paste0(strsplit(write_path_i[length(write_path_i)],".txt")[[1]],'_',i,'.txt')
    print(paste('New enhancer path is ', write_path_i))
    fwrite(e.file[chr == i,],write_path_i,sep='\t')
  }
  
  
  ## GeneList Export 
  path = file.path(paste0(args[1],"GeneList.txt"))
  g.file = fread(path)
  print(paste('There are', dim(g.file)[1],'regions in gene list'))
  chromosomes <- unique(g.file[,chr])
  write_path = path
  
  for(i in chromosomes){
    print(paste('Exporting genes', i,'chromosomes.'))
    write_path_i = write_path
    write_path_i[length(write_path_i)] = paste0(strsplit(write_path_i[length(write_path_i)],".txt")[[1]],'_',i,'.txt')
    print(paste('New enhancer path is ', write_path_i))
    fwrite(g.file[chr == i,],write_path_i,sep='\t')
  }
  
