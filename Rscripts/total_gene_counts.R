library(GenomicAlignments)
library(BiocParallel)
register(MulticoreParam(16))

#source('/home/ev250/Cincinatti/Functions/various.R')

# ebg = "/home/kgoldmann/Documents/PEAC_eqtl/Outputs/b37_ebg.rds"
# in.files = output.df$al.file[1]
# out.files = paste("/home/kgoldmann/Documents/PEAC_eqtl/Outputs/RNA_counts/", output.df$SampleID..QMUL.or.Genentech.,".txt", sep="")
# mode="Union"
# ignore.strand=TRUE

# get args from snakemake
ebg <-  snakemake@input[[1]]
in.files = list.files("/home/kgoldmann/Documents/PEAC_eqtl/Outputs/STAR/2", recursive = TRUE, full.names = TRUE) #snakemake@input[[2]]
in.files = in.files[grepl("Aligned.sortedByCoord.out.bam", in.files)]
#out.files = snakemake@output
mode = snakemake@params[['mode']]
ignore.strand = as.logical(snakemake@params[['ignore_strand']])

ebg <- readRDS(ebg)

for(input in in.files){
  print(input)
  output = paste("/home/kgoldmann/Documents/PEAC_eqtl/Outputs/RNA_counts/", sapply(strsplit(input, "/"), '[[', 9), ".txt", sep="")
  path <-  dirname(input)
  bam.name <- basename(input)
  
  ## transform args to feed function
  sample <- basename(path)
  reads <- basename(dirname(path))
  singleEnd <- ifelse(reads=="Paired", FALSE, TRUE)
  
  ## Prepare matrix of counts per gene:
  
  ##cat(class(ebg), path, bam.name, mode, singleEnd, ignore.strand, out)
  se <- summarizeOverlaps(features=ebg, reads=paste(path, bam.name, sep="/"), mode=mode, ignore.strand=ignore.strand)
  
  counts = assay(se)
  colnames(counts) = sapply(strsplit(input, '/'), '[[', 9)
  #counts <- counts_sample_sub(ebg,path,bam.name, mode,singleEnd, ignore.strand)
  
  #names(counts)[2:ncol(counts)] <- sample
  
  write.table(counts, file=output, row.names=T)
}
