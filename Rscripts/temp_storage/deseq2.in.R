source("/home/kgoldmann/Documents/Git/peac/Rfunctions/inputs.eQTL.dseq2.R")

print(paste0("gene=", snakemake@wildcards[['gene']]))
print(paste0("chr=", snakemake@params[['chrom']]))
#vcfcheat = paste0("/home/kgoldmann/NAS/GCPEAC/Kevin/PEAC-eQTL/KevinAnalysis/ApocritaHPC/PEAC/imputed_vcf/PEAC_1000G_Phase3_Oct14_chr", snakemake@params[['chrom']],".vcf.gz")
print(paste0("vcf=", snakemake@input[['vcf']]))

library(purrr)
in.deseq3 = possibly(in.deseq2, otherwise=NA)

in.deseq3(gene=snakemake@wildcards[['gene']],
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



# 
# #gene, chr, snps=5*10^5,counts.f,gene.coord,vcf, nhets=5,tag.threshold=.9, out=".", prefix=NULL, missing=5
# 
# #gns = read.table("/home/kgoldmann/Documents/PEAC_eqtl/Outputs/gene_coord.txt", header=T)
# gene = "ENSG00000273091"# c("ENSG00000008735", "ENSG00000015475", "ENSG00000025708", "ENSG00000025770", "ENSG00000040608")
# 
# chr=21
# 
# snps=5e-5
# counts.f="/home/kgoldmann/Documents/PEAC_eqtl/Outputs/RNA_counts/groups/Genentech.txt"
# gene.coord="/home/kgoldmann/Documents/PEAC_eqtl/Outputs/gene_coord.txt"
# vcf= '/home/kgoldmann/NAS/GCPEAC/Kevin/PEAC-eQTL/KevinAnalysis/ApocritaHPC/PEAC/imputed_vcf/PEAC_1000G_Phase3_Oct14_chr21.vcf.gz'
# nhets=5
# tag.threshold=0.9
# prefix=NULL
# missing=5
# out='/home/kgoldmann/Documents/PEAC_eqtl/Outputs/deseq2/inputs'
# 
# 
# in.deseq2(gene=gene,
#           chr=chr,
#           snps=snps,
#           counts.f=counts.f,
#           gene.coord=gene.coord,
#           vcf=vcf,
#           nhets=nhets,
#           tag.threshold=tag.threshold,
#           out=out,
#           missing=missing)

# #
