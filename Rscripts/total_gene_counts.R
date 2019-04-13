library(GenomicAlignments, quietly = TRUE, warn.conflicts = FALSE)
library(BiocParallel, quietly = TRUE, warn.conflicts = FALSE)
register(MulticoreParam(16))

#source('/home/ev250/Cincinatti/Functions/various.R')

# ebg = "/home/kgoldmann/Documents/PEAC_eqtl/Outputs/b37_ebg.rds"
# in.files = output.df$al.file[1]
# out.files = paste("/home/kgoldmann/Documents/PEAC_eqtl/Outputs/RNA_counts/", output.df$SampleID..QMUL.or.Genentech.,".txt", sep="")
# mode="Union"
# ignore.strand=TRUE

# get args from snakemake
ebg <-  snakemake@input[[1]] # /home/kgoldmann/Documents/PEAC_eqtl/Outputs/b37_ebg.rds
in.files = snakemake@input[[2]] #list.files("/mnt/volume1/STAR", recursive = TRUE, full.names = TRUE)[1]
#in.files = in.files[grepl("Aligned.sortedByCoord.out.bam", in.files)]
#out.files = snakemake@output
mode = snakemake@params[['mode']]
ignore.strand = as.logical(snakemake@params[['ignore_strand']])
output = snakemake@output[[1]] #paste("/home/kgoldmann/Documents/PEAC_eqtl/Outputs/RNA_counts/", sapply(strsplit(input, "/"), '[[', 7), ".txt", sep="")

ebg <- readRDS(ebg)


print(paste(Sys.time(), ":" , output))

path <-  dirname(in.files)
bam.name <- basename(in.files)

## transform args to feed function
sample <- basename(path)
reads <- basename(dirname(path))
singleEnd <- ifelse(reads=="Paired", FALSE, TRUE)

## Prepare matrix of counts per gene:
print(paste(Sys.time(), ": calculating se"))
se <- summarizeOverlaps(features=ebg, reads=paste(path, bam.name, sep="/"), mode=mode, ignore.strand=ignore.strand)
print(paste(Sys.time(), ": se calculated"))

counts = assay(se)
colnames(counts) = sapply(strsplit(in.files, '/'), '[[', 7)


write.table(counts, file=output, row.names=T)


