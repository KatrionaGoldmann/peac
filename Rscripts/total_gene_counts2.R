require(GenomicAlignments, quietly=TRUE)
require(BiocParallel, quietly=TRUE)

#register(MulticoreParam(16))

#source('/Functions/various.R')

# ebg = "/home/kgoldmann/Documents/PEAC_eqtl/Outputs/b37_ebg.rds"
# mode="Union"
# ignore.strand=TRUE

# get args from snakemake
ebg <-  snakemake@input[[1]]
in.files = snakemake@input[[2]] # in.files = "/media/d1/STAR/Paired/SAM24297981/Aligned.sortedByCoord.out.bam"
out.files = snakemake@output # out.files = "/home/kgoldmann/Documents/PEAC_eqtl/Outputs/RNA_counts/SAM24297981.txt"
mode = snakemake@params[['mode']] 
ignore.strand = as.logical(snakemake@params[['ignore_strand']]) 

ebg <- readRDS(ebg)

singleEnd = FALSE

fun.runner = function(input){
  print(input)
  #output = paste("/home/kgoldmann/Documents/PEAC_eqtl/Outputs/RNA_counts/", sapply(strsplit(input, "/"), '[[', 6), ".txt", sep="")
  path <-  dirname(input)
  bam.name <- basename(input)
  
  ## transform args to feed function
  sample <- basename(path)
  reads <- basename(dirname(path))
  singleEnd <- ifelse(reads=="Paired", FALSE, TRUE)
  
  ## Prepare matrix of counts per gene:
  print(paste(sample, "1"))
  ##cat(class(ebg), path, bam.name, mode, singleEnd, ignore.strand, out)
  se <- summarizeOverlaps(features=ebg, reads=paste(path, bam.name, sep="/"), mode=mode, ignore.strand=ignore.strand)
  print(paste(sample, "2"))
  
  counts = assay(se)
  colnames(counts) = sapply(strsplit(input, '/'), '[[', 6)
  #counts <- counts_sample_sub(ebg,path,bam.name, mode,singleEnd, ignore.strand)
  print(paste(sample, "3"))
  #names(counts)[2:ncol(counts)] <- sample
  write.table(counts, file=as.character(out.files), quote=FALSE)
  
  counts
}

counts = fun.runner(in.files)
print(paste("run for file", out.files))

#zz <- file(description=out.files,"w")

#close(zz)


print("complete")


# set.seed( 123, kind = "L'Ecuyer-CMRG" )
# mclapply(in.files, fun.runner(x), mc.cores=16)
