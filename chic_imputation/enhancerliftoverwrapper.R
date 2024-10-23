library(data.table)
library(GenomicRanges)
library(rtracklayer)

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=2) {
  stop("Wrong number of arguements: Rscript enhancerliftoverwrapper.R <abcenhancerlist>.txt <abcgenelist>.txt", call.=FALSE) } else if (length(args)==2) {
    
    print('Arguements supplied succesfully to the liftover script')}

enhancer_conversion <- function(file_path) {
  # Read the file
  enhancer <- fread(file_path)
  
  # Rest of the code remains the same
  ## conversion to GRanges
  enhancer[,id:=rownames(enhancer)]
  enhancer.gr <- GRanges(seqnames=enhancer$chr, 
                         ranges=IRanges(start=enhancer$start,
                                        end=enhancer$end),id = enhancer$id)
  ## Perform the liftOver
  path = '/home/pa2915/analysis/artemov/ABC_MASTER/scripts/hg19ToHg38.over.chain'
  ch = import.chain(path)
  result <- liftOver(enhancer.gr, ch)
  combinedGR <- do.call(c, result)
  
  ## Duplicated ranges in hg38 
  dups <- combinedGR[duplicated(combinedGR$id),]
  unique <- combinedGR[!(combinedGR$id %in% dups$id),]
  dups$width <- width(dups)
  dups_sorted <- dups[order(dups$id, -width(dups))]
  dups_unique <- dups_sorted[!duplicated(dups_sorted$id)]
  dups_unique <- dups_unique[,"id"]
  message("The length of original is ",length(combinedGR), ", the length of recombined is ", length(dups_unique)+length(unique))
  combinedGR.dedups <- c(unique, dups_unique)
  combinedGR.dedups <- as.data.table(combinedGR.dedups[order(as.numeric(combinedGR.dedups$id)),])
  
  enhancer.hg38 <- combinedGR.dedups[,c("seqnames","start","end","id")]
  colnames(enhancer.hg38) <- c("chr.hg38","start.hg38","end.hg38","id")
  new.enhancer <- merge(enhancer, enhancer.hg38, by.x="id", by.y="id")
  new.enhancer <- new.enhancer[,-c("id","chr","start","end")]
  setnames(new.enhancer,c("chr.hg38","start.hg38","end.hg38"),c("chr","start","end"))
  
  # Save the output
  fwrite(new.enhancer, gsub("\\.txt$", "hg38.txt", file_path),sep='\t')
}

enhancer_conversion(args[1])
enhancer_conversion(args[2])
