
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
    input = list('/home/kgoldmann/Documents/PEAC_eqtl/Outputs/RNA_counts/groups/Genentech.txt', '/home/kgoldmann/Documents/PEAC_eqtl/Outputs/gene_coord.txt', '/home/kgoldmann/NAS/GCPEAC/Kevin/PEAC-eQTL/KevinAnalysis/ApocritaHPC/PEAC/imputed_vcf/PEAC_1000G_Phase3_Oct14_chr21.vcf.gz', '/home/kgoldmann/NAS/GCPEAC/Kevin/PEAC-eQTL/KevinAnalysis/ApocritaHPC/PEAC/imputed_vcf/PEAC_1000G_Phase3_Oct14_chr22.vcf.gz', "counts" = c('/home/kgoldmann/Documents/PEAC_eqtl/Outputs/RNA_counts/groups/Genentech.txt'), "genecoord" = '/home/kgoldmann/Documents/PEAC_eqtl/Outputs/gene_coord.txt', "vcf" = c('/home/kgoldmann/NAS/GCPEAC/Kevin/PEAC-eQTL/KevinAnalysis/ApocritaHPC/PEAC/imputed_vcf/PEAC_1000G_Phase3_Oct14_chr21.vcf.gz', '/home/kgoldmann/NAS/GCPEAC/Kevin/PEAC-eQTL/KevinAnalysis/ApocritaHPC/PEAC/imputed_vcf/PEAC_1000G_Phase3_Oct14_chr22.vcf.gz')),
    output = list('/home/kgoldmann/Documents/PEAC_eqtl/Outputs/deseq2/inputs/ENSG00000015475.rds', "in_deseq" = '/home/kgoldmann/Documents/PEAC_eqtl/Outputs/deseq2/inputs/ENSG00000015475.rds'),
    params = list('22', -55, 5, 0.9, 5, '/home/kgoldmann/Documents/PEAC_eqtl/Outputs/deseq2/inputs', "chrom" = '22', "snps" = -55, "nhets" = 5, "tag" = 0.9, "missing" = 5, "out" = '/home/kgoldmann/Documents/PEAC_eqtl/Outputs/deseq2/inputs'),
    wildcards = list('ENSG00000015475', "gene" = 'ENSG00000015475'),
    threads = 1,
    log = list(),
    resources = list(),
    config = list("scripts_dir" = '/home/kgoldmann/Documents/Git/peac/', "output_dir" = '/home/kgoldmann/Documents/PEAC_eqtl/Outputs', "STAR" = '/home/kgoldmann/Applications/STAR-2.7.0d/source/STAR', "ref_fasta" = '/home/kgoldmann/Documents/PEAC_eqtl/Data/Homo_sapiens.GRCh37.dna.primary_assembly.fa', "ref_gtf" = '/home/kgoldmann/Documents/PEAC_eqtl/Data/Homo_sapiens.GRCh37.87.gtf', "indices" = '/home/kgoldmann/Documents/PEAC_eqtl/Outputs/indices', "ebg" = '/home/kgoldmann/Documents/PEAC_eqtl/Outputs/b37_ebg.rds', "filter" = 100, "sample_file" = '/home/kgoldmann/Documents/PEAC_eqtl/Data/RA.csv', "samples" = '/home/kgoldmann/NAS/RNASEQ/PEAC_RawData/SAM9103822_R1.fastq.gz', "sample_meta" = '/home/kgoldmann/Documents/PEAC_eqtl/Outputs/PEAC_eth.txt', "geno_vcf" = '/home/kgoldmann/Documents/PEAC_eqtl/Data/vcf_list.txt', "ref_bcf" = '/home/kgoldmann/Documents/PEAC_eqtl/Data/RP/refPanel_list.txt', "ref_legend" = '/home/kgoldmann/Documents/PEAC_eqtl/Data/RP/refPanel_leg_list.txt', "N factors" = 10),
    rule = 'Deseq2_inputs',
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
source("/home/kgoldmann/Documents/Git/peac/Rfunctions/inputs.eQTL.dseq2.R")

#genes = snakemake@params[['gene']]

for(i in 1:length(snakemake@params[['gene']])){
in.deseq2(gene=snakemake@params[['gene']][i],
          chr=as.numeric(snakemake@params[['chrom']]),
          snps=as.numeric(snakemake@params[['snps']]),
          counts.f=snakemake@input[['counts']],
          gene.coord=snakemake@input[['genecoord']],
          vcf=snakemake@input[['vcf']],
          nhets=as.numeric(snakemake@params[['nhets']]),
          tag.threshold=as.numeric(snakemake@params[['tag']]),
          out=snakemake@params[['out']],
          missing=as.numeric(snakemake@params[['missing']])
          )
}

#gene, chr, snps=5*10^5,counts.f,gene.coord,vcf, nhets=5,tag.threshold=.9, out=".", prefix=NULL, missing=5

#gns = read.table("/home/kgoldmann/Documents/PEAC_eqtl/Outputs/gene_coord.txt", header=T)
# gene = c("ENSG00000008735", "ENSG00000015475", "ENSG00000025708", "ENSG00000025770", "ENSG00000040608")
# 
# chr=22
# 
# snps=5e-5
# counts.f="/home/kgoldmann/Documents/PEAC_eqtl/Outputs/RNA_counts/groups/Genentech.txt"
# gene.coord="/home/kgoldmann/Documents/PEAC_eqtl/Outputs/gene_coord.txt"
# vcf= '/home/kgoldmann/NAS/GCPEAC/Kevin/PEAC-eQTL/KevinAnalysis/ApocritaHPC/PEAC/imputed_vcf/PEAC_1000G_Phase3_Oct14_chr22.vcf.gz'
# nhets=5
# tag.threshold=0.9
# prefix=NULL
# missing=5
# out='/home/kgoldmann/Documents/PEAC_eqtl/Outputs/deseq2/inputs'
# 
# for (i in gene){
# in.deseq2(gene=i,
#           chr=chr,
#           snps=snps,
#           counts.f=counts.f,
#           gene.coord=gene.coord,
#           vcf=vcf,
#           nhets=nhets,
#           tag.threshold=tag.threshold,
#           out=out,
#           missing=missing)
# }
# 
