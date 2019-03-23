library(gdsfmt)
library(SNPRelate)
library(parallel)

## Reads data from vcf files and convert into SNP GDS Format, calls snpgdsVCF2GDS from SNPRelate package

in.list = c("/home/kgoldmann/Documents/PEAC_eqtl/Outputs/DNA/PEAC_chr22_sub.vcf.gz",
            "/home/kgoldmann/Documents/PEAC_eqtl/Outputs/DNA/RP_chr22_sub.vcf.gz") #expand(config['output_dir'] + "/DNA/RP_chr{chrom}_sub.vcf.gz", chrom=vcf(config["geno_vcf"]).keys())
len = length(in.list)
method="biallelic.only"
##prefix=config['output_dir'] + "/DNA/"
out.list = list("peac"="/home/kgoldmann/Documents/PEAC_eqtl/Outputs/DNA/PEAC_PCA.gds", "rp"="/home/kgoldmann/Documents/PEAC_eqtl/Outputs/DNA/RP_PCA.gds")
vcf.in = list(in.list[1:(len/2)], in.list[(len/2+1):len])

len=length(snakemake@input)
in.list = snakemake@input 
vcf.in = list(in.list[1:(len/2)], in.list[(len/2+1):len])
out.list = snakemake@output


# snpgdsVCF2GDS(vcf.fn=unlist(vcf.in[[1]]),
#               method="biallelic.only",
#               out.fn=out.list[[1]])

set.seed(1234, kind = "L'Ecuyer-CMRG" )
mclapply(1:2 ,function(i)  snpgdsVCF2GDS(vcf.fn=unlist(vcf.in[[i]]), method=snakemake@params[['method']],out.fn=snakemake@output[[i]]), mc.cores = snakemake@threads)

