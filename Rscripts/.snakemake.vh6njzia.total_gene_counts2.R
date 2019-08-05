
######## Snakemake header ########
library(methods)
Snakemake <- setClass(
    "Snakemake",
    slots = c(
        input = "list",
        output = "list",
        params = "list",
        wildcards = "list",
        threads = "numeric",
        log = "list",
        resources = "list",
        config = "list",
        rule = "character"
    )
)
snakemake <- Snakemake(
    input = list('/home/kgoldmann/Documents/PEAC_eqtl/Outputs/b37_ebg.rds', '/media/d1/STAR/Paired/SAM24298030/Aligned.sortedByCoord.out.bam'),
    output = list('/home/kgoldmann/Documents/PEAC_eqtl/Outputs/RNA_counts/SAM24298030.txt'),
    params = list('Union', 'TRUE', "mode" = 'Union', "ignore_strand" = 'TRUE'),
    wildcards = list('SAM24298030', "sample" = 'SAM24298030'),
    threads = 1,
    log = list(),
    resources = list(),
    config = list("scripts_dir" = '/home/kgoldmann/Documents/peac/', "output_dir" = '/home/kgoldmann/Documents/PEAC_eqtl/Outputs', "STAR" = '/home/kgoldmann/Applications/STAR/source/STAR', "ref_fasta" = '/home/kgoldmann/Documents/PEAC_eqtl/Data/Reference/GRCh37/Homo_sapiens.GRCh37.dna.primary_assembly.fa', "ref_gtf" = '/home/kgoldmann/Documents/PEAC_eqtl/Data/Reference/GRCh37/Homo_sapiens.GRCh37.87.gtf', "indices" = '/home/kgoldmann/Documents/PEAC_eqtl/Outputs/indices', "ebg" = '/home/kgoldmann/Documents/PEAC_eqtl/Outputs/b37_ebg.rds', "filter" = 100, "sample_file" = '/home/kgoldmann/Documents/PEAC_eqtl/Data/RA.csv', "samples" = '/home/kgoldmann/NAS/RNASEQ/PEAC_RawData/SAM9103822_R1.fastq.gz', "sample_meta" = '/home/kgoldmann/Documents/PEAC_eqtl/Outputs/PEAC_eth.txt', "geno_vcf" = '/home/kgoldmann/Documents/PEAC_eqtl/Data/vcf_list.txt', "geno_vcf2" = '/home/kgoldmann/Documents/PEAC_eqtl/Data/vcf_list2.txt', "ref_bcf" = '/home/kgoldmann/Documents/PEAC_eqtl/Data/RP/refPanel_list.txt', "ref_legend" = '/home/kgoldmann/Documents/PEAC_eqtl/Data/RP/refPanel_leg_list.txt', "N factors" = 10),
    rule = 'total_gene_counts'
)
######## Original script #########
library(GenomicAlignments)
library(BiocParallel)
#register(MulticoreParam(16))

#source('/Functions/various.R')

# ebg = "/home/kgoldmann/Documents/PEAC_eqtl/Outputs/b37_ebg.rds"
# in.files = output.df$al.file[1]
# out.files = paste("/home/kgoldmann/Documents/PEAC_eqtl/Outputs/RNA_counts/", output.df$SampleID..QMUL.or.Genentech.,".txt", sep="")
# mode="Union"
# ignore.strand=TRUE

# get args from snakemake
ebg <-  snakemake@input[[1]]
#ebg = "/home/kgoldmann/Documents/PEAC_eqtl/Outputs/b37_ebg.rds"
in.files = snakemake@input[[2]] # "/media/d1/STAR/Paired/SAM24297981/Aligned.sortedByCoord.out.bam"
out.files = snakemake@output # "/home/kgoldmann/Documents/PEAC_eqtl/Outputs/RNA_counts/SAM24297981.txt"
mode = snakemake@params[['mode']] 
#mode = "Union"
ignore.strand = as.logical(snakemake@params[['ignore_strand']]) 
# ignore.strand = TRUE



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
  counts
}

counts = fun.runner(in.files)
print(paste("run for file", out.files))

zz <- file(description=out.files,"w")
write.table(counts, file=zz, quote=FALSE)
close(zz)


print("complete")


# set.seed( 123, kind = "L'Ecuyer-CMRG" )
# mclapply(in.files, fun.runner(x), mc.cores=16)
