
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
        rule = "character",
        bench_iteration = "numeric",
        scriptdir = "character",
        source = "function"
    )
)
snakemake <- Snakemake(
    input = list('/home/kgoldmann/Documents/PEAC_eqtl/Outputs/b37_ebg.rds', '/home/kgoldmann/Documents/PEAC_eqtl/Data/FASTQ/LIB12092_SAM9103802_R1.fastq.gz', '/home/kgoldmann/Documents/PEAC_eqtl/Data/FASTQ/LIB12093_SAM9103803_R1.fastq.gz', '/home/kgoldmann/Documents/PEAC_eqtl/Data/FASTQ/LIB12094_SAM9103804_R1.fastq.gz', '/home/kgoldmann/Documents/PEAC_eqtl/Data/FASTQ/LIB12095_SAM9103805_R1.fastq.gz', '/home/kgoldmann/Documents/PEAC_eqtl/Data/FASTQ/LIB12096_SAM9103806_R1.fastq.gz', '/home/kgoldmann/Documents/PEAC_eqtl/Data/FASTQ/LIB12097_SAM9103807_R1.fastq.gz', '/home/kgoldmann/Documents/PEAC_eqtl/Data/FASTQ/LIB12098_SAM9103808_R1.fastq.gz', '/home/kgoldmann/Documents/PEAC_eqtl/Data/FASTQ/LIB12099_SAM9103809_R1.fastq.gz', '/home/kgoldmann/Documents/PEAC_eqtl/Data/FASTQ/LIB12100_SAM9103810_R1.fastq.gz', '/home/kgoldmann/Documents/PEAC_eqtl/Data/FASTQ/LIB12101_SAM9103811_R1.fastq.gz', '/home/kgoldmann/Documents/PEAC_eqtl/Data/FASTQ/LIB12102_SAM9103812_R1.fastq.gz', '/home/kgoldmann/Documents/PEAC_eqtl/Data/FASTQ/LIB12103_SAM9103813_R1.fastq.gz', '/home/kgoldmann/Documents/PEAC_eqtl/Data/FASTQ/LIB12104_SAM9103814_R1.fastq.gz', '/home/kgoldmann/Documents/PEAC_eqtl/Data/FASTQ/LIB12105_SAM9103815_R1.fastq.gz', '/home/kgoldmann/Documents/PEAC_eqtl/Data/FASTQ/LIB12106_SAM9103816_R1.fastq.gz', '/home/kgoldmann/Documents/PEAC_eqtl/Data/FASTQ/LIB12107_SAM9103817_R1.fastq.gz', '/home/kgoldmann/Documents/PEAC_eqtl/Data/FASTQ/LIB12108_SAM9103818_R1.fastq.gz', '/home/kgoldmann/Documents/PEAC_eqtl/Data/FASTQ/LIB12109_SAM9103819_R1.fastq.gz', '/home/kgoldmann/Documents/PEAC_eqtl/Data/FASTQ/LIB12110_SAM9103820_R1.fastq.gz', '/home/kgoldmann/Documents/PEAC_eqtl/Data/FASTQ/LIB12111_SAM9103821_R1.fastq.gz', '/home/kgoldmann/Documents/PEAC_eqtl/Data/FASTQ/LIB12112_SAM9103822_R1.fastq.gz', '/home/kgoldmann/Documents/PEAC_eqtl/Data/FASTQ/LIB12113_SAM9103823_R1.fastq.gz', '/home/kgoldmann/Documents/PEAC_eqtl/Data/FASTQ/LIB12114_SAM9103824_R1.fastq.gz', '/home/kgoldmann/Documents/PEAC_eqtl/Data/FASTQ/LIB12115_SAM9103825_R1.fastq.gz', '/home/kgoldmann/Documents/PEAC_eqtl/Data/FASTQ/LIB12116_SAM9103826_R1.fastq.gz', '/home/kgoldmann/Documents/PEAC_eqtl/Data/FASTQ/LIB12117_SAM9103827_R1.fastq.gz', '/home/kgoldmann/Documents/PEAC_eqtl/Data/FASTQ/LIB12118_SAM9103828_R1.fastq.gz', '/home/kgoldmann/Documents/PEAC_eqtl/Data/FASTQ/LIB12119_SAM9103829_R1.fastq.gz', '/home/kgoldmann/Documents/PEAC_eqtl/Data/FASTQ/LIB12120_SAM9103830_R1.fastq.gz', '/home/kgoldmann/Documents/PEAC_eqtl/Data/FASTQ/LIB12121_SAM9103832_R1.fastq.gz', '/home/kgoldmann/Documents/PEAC_eqtl/Data/FASTQ/LIB12122_SAM9103833_R1.fastq.gz', '/home/kgoldmann/Documents/PEAC_eqtl/Data/FASTQ/LIB12123_SAM9103834_R1.fastq.gz', '/home/kgoldmann/Documents/PEAC_eqtl/Data/FASTQ/LIB12124_SAM9103835_R1.fastq.gz', '/home/kgoldmann/Documents/PEAC_eqtl/Data/FASTQ/LIB12125_SAM9103836_R1.fastq.gz', '/home/kgoldmann/Documents/PEAC_eqtl/Data/FASTQ/LIB12126_SAM9103837_R1.fastq.gz', '/home/kgoldmann/Documents/PEAC_eqtl/Data/FASTQ/LIB12127_SAM9103838_R1.fastq.gz', '/home/kgoldmann/Documents/PEAC_eqtl/Data/FASTQ/LIB12128_SAM9103839_R1.fastq.gz'),
    output = list('/home/kgoldmann/Documents/PEAC_eqtl/Outputs/RNA_counts/SAM9103802.txt', '/home/kgoldmann/Documents/PEAC_eqtl/Outputs/RNA_counts/SAM9103803.txt', '/home/kgoldmann/Documents/PEAC_eqtl/Outputs/RNA_counts/SAM9103804.txt', '/home/kgoldmann/Documents/PEAC_eqtl/Outputs/RNA_counts/SAM9103805.txt', '/home/kgoldmann/Documents/PEAC_eqtl/Outputs/RNA_counts/SAM9103806.txt', '/home/kgoldmann/Documents/PEAC_eqtl/Outputs/RNA_counts/SAM9103807.txt', '/home/kgoldmann/Documents/PEAC_eqtl/Outputs/RNA_counts/SAM9103808.txt', '/home/kgoldmann/Documents/PEAC_eqtl/Outputs/RNA_counts/SAM9103809.txt', '/home/kgoldmann/Documents/PEAC_eqtl/Outputs/RNA_counts/SAM9103810.txt', '/home/kgoldmann/Documents/PEAC_eqtl/Outputs/RNA_counts/SAM9103811.txt', '/home/kgoldmann/Documents/PEAC_eqtl/Outputs/RNA_counts/SAM9103812.txt', '/home/kgoldmann/Documents/PEAC_eqtl/Outputs/RNA_counts/SAM9103813.txt', '/home/kgoldmann/Documents/PEAC_eqtl/Outputs/RNA_counts/SAM9103814.txt', '/home/kgoldmann/Documents/PEAC_eqtl/Outputs/RNA_counts/SAM9103815.txt', '/home/kgoldmann/Documents/PEAC_eqtl/Outputs/RNA_counts/SAM9103816.txt', '/home/kgoldmann/Documents/PEAC_eqtl/Outputs/RNA_counts/SAM9103817.txt', '/home/kgoldmann/Documents/PEAC_eqtl/Outputs/RNA_counts/SAM9103818.txt', '/home/kgoldmann/Documents/PEAC_eqtl/Outputs/RNA_counts/SAM9103819.txt', '/home/kgoldmann/Documents/PEAC_eqtl/Outputs/RNA_counts/SAM9103820.txt', '/home/kgoldmann/Documents/PEAC_eqtl/Outputs/RNA_counts/SAM9103821.txt', '/home/kgoldmann/Documents/PEAC_eqtl/Outputs/RNA_counts/SAM9103822.txt', '/home/kgoldmann/Documents/PEAC_eqtl/Outputs/RNA_counts/SAM9103823.txt', '/home/kgoldmann/Documents/PEAC_eqtl/Outputs/RNA_counts/SAM9103824.txt', '/home/kgoldmann/Documents/PEAC_eqtl/Outputs/RNA_counts/SAM9103825.txt', '/home/kgoldmann/Documents/PEAC_eqtl/Outputs/RNA_counts/SAM9103826.txt', '/home/kgoldmann/Documents/PEAC_eqtl/Outputs/RNA_counts/SAM9103827.txt', '/home/kgoldmann/Documents/PEAC_eqtl/Outputs/RNA_counts/SAM9103828.txt', '/home/kgoldmann/Documents/PEAC_eqtl/Outputs/RNA_counts/SAM9103829.txt', '/home/kgoldmann/Documents/PEAC_eqtl/Outputs/RNA_counts/SAM9103830.txt', '/home/kgoldmann/Documents/PEAC_eqtl/Outputs/RNA_counts/SAM9103832.txt', '/home/kgoldmann/Documents/PEAC_eqtl/Outputs/RNA_counts/SAM9103833.txt', '/home/kgoldmann/Documents/PEAC_eqtl/Outputs/RNA_counts/SAM9103834.txt', '/home/kgoldmann/Documents/PEAC_eqtl/Outputs/RNA_counts/SAM9103835.txt', '/home/kgoldmann/Documents/PEAC_eqtl/Outputs/RNA_counts/SAM9103836.txt', '/home/kgoldmann/Documents/PEAC_eqtl/Outputs/RNA_counts/SAM9103837.txt', '/home/kgoldmann/Documents/PEAC_eqtl/Outputs/RNA_counts/SAM9103838.txt', '/home/kgoldmann/Documents/PEAC_eqtl/Outputs/RNA_counts/SAM9103839.txt'),
    params = list('Union', 'TRUE', "mode" = 'Union', "ignore_strand" = 'TRUE'),
    wildcards = list(),
    threads = 1,
    log = list(),
    resources = list(),
    config = list("scripts_dir" = '/home/kgoldmann/Documents/Git/peac/', "output_dir" = '/home/kgoldmann/Documents/PEAC_eqtl/Outputs', "STAR" = '/home/kgoldmann/Applications/STAR-2.7.0d/source/STAR', "ref_fasta" = '/home/kgoldmann/Documents/PEAC_eqtl/Data/Homo_sapiens.GRCh37.dna.primary_assembly.fa', "ref_gtf" = '/home/kgoldmann/Documents/PEAC_eqtl/Data/Homo_sapiens.GRCh37.87.gtf', "indices" = '/home/kgoldmann/Documents/PEAC_eqtl/Outputs/indices', "ebg" = '/home/kgoldmann/Documents/PEAC_eqtl/Outputs/b37_ebg.rds', "filter" = 100, "sample_file" = '/home/kgoldmann/Documents/PEAC_eqtl/Data/RA.csv', "samples" = '/home/kgoldmann/NAS/RNASEQ/PEAC_RawData/SAM9103822_R1.fastq.gz', "sample_meta" = '/mrc-bsu/scratch/ev250/emedlab/RNA/objects/PEAC_eth.txt', "geno_vcf" = '/home/kgoldmann/Documents/PEAC_eqtl/Data/vcf_list.txt', "ref_bcf" = '/home/kgoldmann/Documents/PEAC_eqtl/Data/RP/refPanel_list.txt', "ref_legend" = '/home/kgoldmann/Documents/PEAC_eqtl/Data/RP/refPanel_leg_list.txt', "N factors" = 10),
    rule = 'total_gene_counts',
    bench_iteration = as.numeric(NA),
    scriptdir = '/home/kgoldmann/Documents/Git/peac/Rscripts',
    source = function(...){
        wd <- getwd()
        setwd(snakemake@scriptdir)
        source(...)
        setwd(wd)
    }
)

######## Original script #########
#source('/home/ev250/Cincinatti/Functions/various.R')


## get args from snakemake
ebg <-  snakemake@input[[1]]
in.files = read.csv(snakemake@input[[2]], header=FALSE)
out.files = paste("/home/kgoldmann/Documents/PEAC_eqtl/Outputs/RNA_counts/", in.files$V1 ,".txt", sep="")

for(i in 1:nrow(in.files)){
  file = as.character(in.files$V14)[i]
  path <-  dirname(file)
  bam.name <- basename(file)
  
  mode <- snakemake@params[['mode']]
  ignore.strand <- as.logical(snakemake@params[['ignore_strand']])
  
  out <- out.files[i]
  
  ## transform args to feed function
  
  sample <- basename(path)
  reads <- basename(dirname(path))
  
  ebg <- readRDS(ebg)
  
  singleEnd <- ifelse(reads=="Paired", FALSE, TRUE)
  
  ## Prepare matrix of counts per gene:
  
  ##cat(class(ebg), path, bam.name, mode, singleEnd, ignore.strand, out)
  
  counts <- counts_sample_sub(ebg,path,bam.name, mode,singleEnd, ignore.strand)
  
  names(counts)[2:ncol(counts)] <- sample
  
  write.table(counts, file=out , row.names=F)
}