library(gdsfmt)
library(SNPRelate)
library(parallel)

## Reads data from vcf files and convert into SNP GDS Format, calls snpgdsVCF2GDS from SNPRelate package



len=length(snakemake@input)
in.list = snakemake@input 
vcf.in = list(in.list[1:(len/2)], in.list[(len/2+1):len])
out.list = snakemake@output

mclapply(1:2, function(i) snpgdsVCF2GDS(vcf.fn=unlist(vcf.in[[i]]),
              method="biallelic.only",
              out.fn=out.list[[i]]), mc.cores=2)

# in.list = c("/home/kgoldmann/Documents/PEAC_eqtl/Outputs/DNA/PEAC_chr20_sub.vcf.gz",
#             "/home/kgoldmann/Documents/PEAC_eqtl/Outputs/DNA/RP_chr20_sub.vcf.gz") #expand(config['output_dir'] + "/DNA/RP_chr{chrom}_sub.vcf.gz", chrom=vcf(config["geno_vcf"]).keys())
# len = length(in.list)
# method="biallelic.only"
# # ##prefix=config['output_dir'] + "/DNA/"
# out.list = list("peac"="/home/kgoldmann/Documents/PEAC_eqtl/Outputs/DNA/PEAC_PCA.gds", "rp"="/home/kgoldmann/Documents/PEAC_eqtl/Outputs/DNA/RP_PCA_trial.gds")
# vcf.in = list(in.list[1:(len/2)], in.list[(len/2+1):len])
# 
# 
# set.seed(1234, kind = "L'Ecuyer-CMRG" )
# library(parallel)
# 
# library(GWASTools)
# ?convertVcfGds
# library(SNPRelate)
# snpgdsVCF2GDS(unlist(vcf.in[[2]]), out.list[[2]])
# 
# snpgdsVCF2GDS(vcf.fn=unlist(vcf.in[[2]]), method="biallelic.only", out.fn=out.list[[2]])
# 
# mclapply(1:2 ,function(i)  snpgdsVCF2GDS(vcf.fn=unlist(vcf.in[[i]]), method='copy.num.of.ref',out.fn=out.list[[i]]), mc.cores = 2)
# library(vcfR)
# x = read.vcfR(in.list[2], verbose=FALSE)
# get = extract.gt(x, "GT")
# 
# # find the snps with more alts
# which(apply(get[1:1000, ], 1, function(y) "2|0" %in% y))
# table(get["rs188248858", ])
# 
# # what do those info look like
# temp = getFIX(x)
# temp[temp[, "ID"] == "rs188248858", ]
# any(grepl(",", temp[, "ALT"]))
# 
# # Error in snpgdsVCF2GDS(unlist(vcf.in[[2]]), out.list[[2]]) : 
# #   Genotype code is out of range "2".
# # FILE: /home/kgoldmann/Documents/PEAC_eqtl/Outputs/DNA/RP_chr9_trial.vcf.gz.gz
# # LINE: 380, COLUMN: 119, 0|2
# 
# which(is.biallelic(x)!=TRUE)
# 
# get2 = extract.info(x, "GT")
# table()
# 
