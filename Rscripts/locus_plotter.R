# Data visualisations
# http://bioconductor.org/packages/release/bioc/vignettes/ggbio/inst/doc/ggbio.pdf

library(data.table)
library(ggplot2)
library(ggrepel)
library(ggbio)
library(ggpubr)
library(EnsDb.Hsapiens.v75)
data(genesymbol, package = "biovizBase")

# Load the script for locuszoom
source("/home/kgoldmann/Documents/peac/Rfunctions/locuszoom.R")


covariates = "EV1-PEER1"

# Synovium
############
file.in = paste0("/home/kgoldmann/Documents/PEAC_eqtl/Results/Synovium/", covariates, "/Synovium_", covariates, ".csv")
SNP_file = "/media/d1/KG_Outputs/Syn_out_KG/matqtl/inputs/genotype.txt"
gene.pos = "/media/d1/KG_Outputs/Syn_out_KG//matqtl/inputs/gene_location.txt"
exp.in = fread("/media/d1/KG_Outputs/Syn_out_KG/matqtl/inputs/gene_expression_cqn.txt")
window=1e6

# Plot for only the significant genes
all = fread(file.in, nrows=20000)
all = all[all$symbol != "" & all$pvalue < 1e-8]
all$count = table(all$symbol)[match(all$symbol, names(table(all$symbol)))]
all = all[rev(order(all$count)), ]

# Output the information for each gene
genes.info=data.frame("gene"=unique(all$symbol), "Ensembl"=unique(all$gene))
genes.info$"sig.snps"= all$count[match(genes.info$gene, all$symbol)]
genes.info$top.snp = NA
genes.info$top.p = NA
for(i in 1:nrow(genes.info)){
  genes.info$top.snp[i] = all$rs.id[all$symbol == genes.info$gene[i]][1]
  genes.info$top.p[i] = all$pvalue[all$symbol == genes.info$gene[i]][1]
}
genes.info = genes.info[order(genes.info$top.p), ]
write.csv(genes.info, paste0("/home/kgoldmann/Documents/PEAC_eqtl/Results/Synovium/", covariates, "/sig_genes.csv"), quote=F, row.names=F)


common.genes = unique(all$symbol[all$count > 1])[1:50]
top.genes = as.character(data.frame(fread(file.in, select="symbol", nrows=100))$symbol)
top.genes = top.genes[top.genes != ""][1:50]
all.genes = unique(c(common.genes, top.genes))
all.genes = unique(c(all.genes[! is.na(all.genes)], c("ERAP1", "ERAP2")))

load("/home/kgoldmann/Documents/gcpeac/Chris/PEAC_data/PEACmetadataV12.RData")
m1 = metadata$baseline
m1 = m1[match(colnames(exp.in)[2:ncol(exp.in)], m1$GenentechID), ]

snp.cors = data.frame()

invisible(file.remove(list.files(paste0("/home/kgoldmann/Documents/PEAC_eqtl/Results/Synovium/", covariates, "/locuszoom/"), full.names = T)))
for(gene in all.genes){
  print(paste(gene, which(all.genes == gene), "/", length(all.genes)))
  
  # check which genes.df within the plot neighbourhood for plotting
  genes.df = fread(gene.pos)
  gp = genes.df$start[genes.df$symbol == gene] + (genes.df$end[genes.df$symbol == gene] - genes.df$start[genes.df$symbol == gene])/2
  gene.range = c(gp-window/2, gp+window/2)
  chrom = genes.df$chrom[genes.df$symbol == gene]
  g.start = genes.df$start[genes.df$symbol == gene]
  g.end = genes.df$end[genes.df$symbol == gene]
  
  genes.df = genes.df[genes.df$start >= gene.range[1] & genes.df$end <= gene.range[2] & genes.df$chrom == genes.df$chrom[genes.df$symbol == gene] & genes.df$symbol != "", ]
  genes.df=genes.df[order(genes.df$start), ]
  
  plots=local_plotter(sym=gene, chr=chrom, width=window, file=file.in, col = "symbol", SNP_file_name = SNP_file, 
                      ld.calc=TRUE, genes.plot=genes.df$symbol, highlight.gene = T, g.start=g.start, g.end=g.end, 
                      exp = exp.in, meta=m1)
  

  
  # Plot if snps comply with distribution rules
  if(! "reason" %in% names(plots)){
    snp.cors = rbind(snp.cors, plots$snp.corr)
    
    p1 = ggarrange(plots$snps, plots$genes, nrow=2, ncol=1, align="v", legend="bottom", common.legend=T, heights = c(0.6, 0.4), widths=c(0.6, 0.4))#, nrow=2, heights=c(0.1, 0.9))
    p2 = ggarrange(plots$bp, plots$sp, nrow=2, ncol=1, heights=c(0.5, 1))
    
    p = ggarrange(p1, p2, nrow=1, ncol=2, widths=c(0.6, 0.4))
    
    pdf(paste0("/home/kgoldmann/Documents/PEAC_eqtl/Results/Synovium/", covariates, "/locuszoom/Synovium_", gene, ".pdf"), onefile=FALSE, width=10)
    print(annotate_figure(p, top = paste("SNPs correlated with", gene)))
    dev.off() 
    
    
    
  }
}

write.csv(snp.cors[snp.cors$Significant == T, ], paste0("/home/kgoldmann/Documents/PEAC_eqtl/Results/Synovium/", covariates, "/snp_clin_corrs.csv"), quote=F, row.names = F)

# sudo rm /home/kgoldmann/Documents/gcpeac/Katriona/PEAC_analysis/PEAC_eQTL/Results/Synovium/EV1-PEER1/locuszoom/*
# sudo head -1000 /home/kgoldmann/Documents/PEAC_eqtl/Results/Synovium/EV1-PEER1/Synovium_EV1-PEER1.csv > /home/kgoldmann/Documents/PEAC_eqtl/Results/Synovium/EV1-PEER1/Synovium_EV1-PEER1_top1000.csv
# sudo cp -r /home/kgoldmann/Documents/PEAC_eqtl/Results/Synovium/EV1-PEER1/Synovium_EV1-PEER1_top1000.csv /home/kgoldmann/Documents/gcpeac/Katriona/PEAC_analysis/PEAC_eQTL/Results/Synovium/EV1-PEER1/Synovium_EV1-PEER1_top1000.csv





# Blood
############
file.in = "/home/kgoldmann/Documents/PEAC_eqtl/Results/Blood/EV1-PEER4/Blood_EV1-PEER4.csv"
SNP_file = "/home/kgoldmann/Documents/PEAC_eqtl/Outputs_Blood/matqtl/inputs/genotype.txt"
gene.pos = "/home/kgoldmann/Documents/PEAC_eqtl/Outputs_Blood/matqtl/inputs/gene_location.txt"
exp.in = fread("/home/kgoldmann/Documents/PEAC_eqtl/Outputs_Blood/matqtl/inputs/gene_expression_cqn.txt")

window=1e6

all = fread(file.in, nrows=20000)
all = all[all$symbol != "" & all$pvalue < 1e-8]
all$count = table(all$symbol)[match(all$symbol, names(table(all$symbol)))]
all = all[rev(order(all$count)), ]

# Output the information for each gene
genes.info=data.frame("gene"=unique(all$symbol), "Ensembl"=unique(all$gene))
genes.info$"sig.snps"= all$count[match(genes.info$gene, all$symbol)]
genes.info$top.snp = NA
genes.info$top.p = NA
for(i in 1:nrow(genes.info)){
  genes.info$top.snp[i] = all$rs.id[all$symbol == genes.info$gene[i]][1]
  genes.info$top.p[i] = all$pvalue[all$symbol == genes.info$gene[i]][1]
}
genes.info = genes.info[order(genes.info$top.p), ]
write.csv(genes.info, paste0("/home/kgoldmann/Documents/PEAC_eqtl/Results/Blood/", covariates, "/sig_genes.csv"), quote=F, row.names=F)


common.genes = unique(all$symbol[all$count > 1])[1:50]
top.genes = as.character(data.frame(fread(file.in, select="symbol", nrows=100))$symbol)
top.genes = top.genes[top.genes != ""][1:50]
all.genes = unique(c(common.genes, top.genes))
all.genes = unique(c(all.genes[! is.na(all.genes)], c("ERAP1", "ERAP2")))


for(gene in all.genes){
  print(paste(gene, which(all.genes == gene), "/", length(all.genes)))
  # check which genes within the plot neighbourhood for plotting
  genes.df = fread(gene.pos)
  gp = genes.df$start[genes.df$symbol == gene] + (genes.df$end[genes.df$symbol == gene] - genes.df$start[genes.df$symbol == gene])/2
  gene.range = c(gp-window/2, gp+window/2)
  chrom = genes.df$chrom[genes.df$symbol == gene]
  g.start = genes.df$start[genes.df$symbol == gene]
  g.end = genes.df$end[genes.df$symbol == gene]
  
  genes.df = genes.df[genes.df$start >= gene.range[1] & genes.df$end <= gene.range[2] & genes.df$chrom == genes.df$chrom[genes.df$symbol == gene] & genes.df$symbol != "", ]
  genes.df=genes.df[order(genes.df$start), ]
  
  plots=local_plotter(sym=gene, chr=chrom, width=window, file=file.in, col = "symbol", SNP_file_name = SNP_file, 
                      ld.calc=TRUE, genes.plot=genes.df$symbol, highlight.gene = T, g.start=g.start, g.end=g.end, 
                      exp = exp.in)
  p1 = ggarrange(plots$snps, plots$genes, nrow=2, ncol=1, align="v", legend="bottom", common.legend=T, heights = c(0.6, 0.4), widths=c(0.6, 0.4))#, nrow=2, heights=c(0.1, 0.9))
  p2 = ggarrange(plots$bp, plots$table, nrow=3, ncol=1, heights=c(0.5, 0.3, 0.7))
  p = ggarrange(p1, p2, nrow=1, ncol=2, widths=c(0.7, 0.3))
  
  pdf(paste0("/home/kgoldmann/Documents/PEAC_eqtl/Results/Blood/EV1-PEER4/locuszoom/Blood_", gene, ".pdf"), onefile=FALSE)
  print(annotate_figure(p, top = paste("SNPs correlated with", gene)))
  dev.off()
}

# sudo cp -r /home/kgoldmann/Documents/PEAC_eqtl/Results/Blood/EV1-PEER4/locuszoom/* /home/kgoldmann/Documents/gcpeac/Katriona/PEAC_analysis/PEAC_eQTL/Results/Blood/EV1-PEER4/locuszoom/

