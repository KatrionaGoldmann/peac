# Data visualisations
# http://bioconductor.org/packages/release/bioc/vignettes/ggbio/inst/doc/ggbio.pdf

library(data.table)
library(ggplot2)
library(ggrepel)
library(ggbio)
library(gdata)
library(ggpubr)
library(EnsDb.Hsapiens.v75)
# data(genesymbol, package = "biovizBase")

# Load the script for locuszoom
source("/home/kgoldmann/Documents/peac/Rfunctions/locuszoom.R")

############
# Synovium
############
SNP_file = "/media/d1/KG_Outputs/Syn_out_KG/matqtl/inputs/genotype.txt"
gene.pos = "/media/d1/KG_Outputs/Syn_out_KG/matqtl/inputs/gene_location.txt"
exp.in = fread("/media/d1/KG_Outputs/Syn_out_KG/matqtl/inputs/gene_expression_cqn.txt")
window=1e6
meta = read.table("/home/kgoldmann/Documents/PEAC_eqtl/Data/PEAC/PEAC_eth_syn.txt")

cvs = list.files("/home/kgoldmann/Documents/PEAC_eqtl/Results/Synovium", full.names = FALSE, recursive = T)
cvs = cvs[! grepl("summary|locuszoom|sig_genes|snp_clin_corrs", cvs)]
cvs = gsub("\\/.*", "", cvs)

for (covariates in cvs){
  
  file.in = paste0("/home/kgoldmann/Documents/PEAC_eqtl/Results/Synovium/", covariates, "/Synovium_", covariates, ".csv")

  # Plot for only the significant genes
  all = fread(file.in, nrows=20000)
  all = all[all$symbol != "" & all$pvalue < 5e-8]
  all$count = table(all$gene)[match(all$gene, names(table(all$gene)))]
  all = all[rev(order(all$count)), ]
  # Make note of the haplotypic genes
  hap.genes = all[! duplicated(paste(all$symbol, all$gene))]
  hap.genes = sapply(names(which(table(hap.genes$symbol) > 1)), function(x) unique(hap.genes$gene[hap.genes$symbol == x]))

  # Output the information for each gene
  genes.info=data.frame("gene"=all$symbol[match(unique(all$gene), all$gene)], "Ensembl"=unique(all$gene))
  genes.info$"sig.snps"= all$count[match(genes.info$gene, all$symbol)]
  genes.info$top.snp = NA
  genes.info$top.p = NA
  for(i in 1:nrow(genes.info)){
    genes.info$top.snp[i] = all$rs.id[all$symbol == genes.info$gene[i]][1]
    genes.info$top.p[i] = all$pvalue[all$symbol == genes.info$gene[i]][1]
  }
  genes.info = genes.info[order(genes.info$top.p), ]
  write.csv(genes.info, file=paste0("/home/kgoldmann/Documents/PEAC_eqtl/Results/Synovium/", covariates, "/sig_genes.csv"), quote=F, row.names=F)
  
  common.genes = unique(all$gene[all$count > 1 & all$maf > 0.05])[1:50]
  top.genes = as.character(data.frame(fread(file.in, select="gene", nrows=100))$gene)
  top.genes = top.genes[top.genes != ""][1:50]
  all.genes = unique(c(common.genes, top.genes))
  all.genes = unique(c(all.genes[! is.na(all.genes)]))
  #names(hap.genes)[sapply(hap.genes, function(x) length(which(x %in% all.genes)) > 1)]
  #all.genes = all.genes[! all.genes %in% hap.genes]
  
  pass = all[all$gene %in% all.genes, ]
  pass = pass %>%  distinct(gene, symbol, .keep_all = T)
  pass = sapply(unique(pass$gene), function(x) unique(pass$symbol[pass$gene == x]))
  
  #load("/home/kgoldmann/Documents/gcpeac/Chris/PEAC_data/PEACmetadataV12.RData")
  m1 = meta
  m1 = m1[match(colnames(exp.in)[2:ncol(exp.in)], m1$samp.et), ]
  
  snp.cors = data.frame()
  
  # load("/media/gcpeac/Chris/PEAC_data/PEACmetadataV12.RData")
  # load("/media/gcpeac/Chris/PEAC_data/PEACexpressiondataV2.RData")
  # vst = filteredexpression$vst
  # remove(filteredexpression, metadata, txiraw)
  vst = readRDS("/media/d1/KG_Outputs/Syn_out_KG/matqtl/inputs/vst.rds")
  
  
  dir.create(paste0("/home/kgoldmann/Documents/PEAC_eqtl/Results/Synovium/", covariates, "/locuszoom/"), showWarnings = FALSE)
  #invisible(file.remove(list.files(paste0("/home/kgoldmann/Documents/PEAC_eqtl/Results/Synovium/", covariates, "/locuszoom/"), full.names = T)))
  for(k in 1:length(pass)){
    gene = pass[k]
    print(paste0(gene, "   ", k, "/", length(pass)))
    
    # check which genes.df within the plot neighbourhood for plotting
    genes.df = fread(gene.pos)
    gp = genes.df$start[genes.df$gene_id == names(pass)[k] ] + (genes.df$end[genes.df$gene_id ==  names(pass)[k]  ] - genes.df$start[genes.df$gene_id == names(pass)[k] ])/2
    gene.range = c(gp-window/2, gp+window/2)
    chrom = genes.df$chrom[genes.df$gene_id ==  names(pass)[k] ]
    g.start = genes.df$start[genes.df$gene_id ==  names(pass)[k] ]
    g.end = genes.df$end[genes.df$gene_id ==  names(pass)[k] ]
    
    genes.df = genes.df[genes.df$start >= gene.range[1] & genes.df$end <= gene.range[2] & genes.df$chrom == genes.df$chrom[genes.df$gene_id %in%  names(pass)[k] ] & genes.df$symbol != "", ]
    genes.df=genes.df[order(genes.df$start), ]
    
    plots=local_plotter(sym=gene, chr=chrom, width=window, file=file.in, col = "symbol", SNP_file_name = SNP_file, 
                        ld.calc=TRUE, genes.plot=genes.df$gene_id, highlight.gene = T, g.start=g.start, g.end=g.end, 
                        exp = exp.in, meta=m1, vst=vst)
    
    
    # Plot if snps comply with distribution rules
    if(! "reason" %in% names(plots)){
      snp.cors = rbind(snp.cors, plots$snp.corr)
      
      p1 = ggarrange(plots$snps, plots$legend, plots$genes, nrow=2, ncol=2, align="v", heights = c(0.6, 0.4), widths=c(0.6, 0.15))#, nrow=2, heights=c(0.1, 0.9))
      p2 = ggarrange(plots$bp, plots$sp, nrow=2, ncol=1, heights=c(0.5, 1))
      
      p = ggarrange(p1, p2, nrow=1, ncol=2, widths=c(0.6, 0.4))
      
      
      
      ggsave(print(annotate_figure(p, top = paste("SNPs correlated with", gene))), 
             filename =paste0("/home/kgoldmann/Documents/PEAC_eqtl/Results/Synovium/", covariates, "/locuszoom/Synovium_", gene, "_", names(pass)[k], ".pdf"), 
             device = cairo_pdf)
      
      # pdf(paste0("/home/kgoldmann/Documents/PEAC_eqtl/Results/Synovium/", covariates, "/locuszoom/Synovium_", gene, "_", names(pass)[k], ".pdf"), onefile=FALSE, width=10)
      # print(annotate_figure(p, top = paste("SNPs correlated with", gene)))
      # dev.off() 
      # 
      
      
    }
  }
  
  write.csv(snp.cors[snp.cors$Significant == T, ], paste0("/home/kgoldmann/Documents/PEAC_eqtl/Results/Synovium/", covariates, "/snp_clin_corrs.csv"), quote=F, row.names = F)
  
}


# sudo rm -r /home/kgoldmann/Documents/gcpeac/Katriona/PEAC_analysis/PEAC_eQTL/Results/Synovium/EV1234_PEER134/
# sudo cp -r /home/kgoldmann/Documents/PEAC_eqtl/Results/Synovium/EV1234_PEER1234/locuszoom/ /home/kgoldmann/Documents/gcpeac/Katriona/PEAC_analysis/PEAC_eQTL/Results/Synovium/EV1234_PEER134/



############
# Blood
############
SNP_file = "/media/d1/KG_Outputs/Bld_out_KG/matqtl/inputs/genotype.txt"
gene.pos = "/media/d1/KG_Outputs/Bld_out_KG/matqtl/inputs/gene_location.txt"    
exp.in = fread("/media/d1/KG_Outputs/Bld_out_KG/matqtl/inputs/gene_expression_cqn.txt")
window=1e6
meta = read.table("/home/kgoldmann/Documents/PEAC_eqtl/Data/PEAC/PEAC_eth_bld.txt")

cvs = list.files("/home/kgoldmann/Documents/PEAC_eqtl/Results/Synovium", full.names = FALSE, recursive = T)
cvs = cvs[! grepl("summary|locuszoom|sig_genes|snp_clin_corrs", cvs)]
cvs = gsub("\\/.*", "", cvs)

for (covariates in cvs){
  print(paste("###", covariates))
  file.in = paste0("/home/kgoldmann/Documents/PEAC_eqtl/Results/Blood/", covariates, "/Blood_", covariates, ".csv")
  
  # Plot for only the significant genes
  all = fread(file.in, nrows=20000)
  all = all[all$symbol != "" & all$pvalue < 5e-8]
  all$count = table(all$gene)[match(all$gene, names(table(all$gene)))]
  all = all[rev(order(all$count)), ]
  # Make note of the haplotypic genes
  hap.genes = all[! duplicated(paste(all$symbol, all$gene))]
  hap.genes = sapply(names(which(table(hap.genes$symbol) > 1)), function(x) unique(hap.genes$gene[hap.genes$symbol == x]))
  
  # Output the information for each gene
  genes.info=data.frame("gene"=all$symbol[match(unique(all$gene), all$gene)], "Ensembl"=unique(all$gene))
  genes.info$"sig.snps"= all$count[match(genes.info$gene, all$symbol)]
  genes.info$top.snp = NA
  genes.info$top.p = NA
  for(i in 1:nrow(genes.info)){
    genes.info$top.snp[i] = all$rs.id[all$symbol == genes.info$gene[i]][1]
    genes.info$top.p[i] = all$pvalue[all$symbol == genes.info$gene[i]][1]
  }
  genes.info = genes.info[order(genes.info$top.p), ]
  write.csv(genes.info, paste0("/home/kgoldmann/Documents/PEAC_eqtl/Results/Blood/", covariates, "/sig_genes.csv"), quote=F, row.names=F)
  
  common.genes = unique(all$gene[all$count > 1 & all$maf > 0.05])[1:50]
  top.genes = as.character(data.frame(fread(file.in, select="gene", nrows=100))$gene)
  top.genes = top.genes[top.genes != ""][1:50]
  all.genes = unique(c(common.genes, top.genes))
  all.genes = unique(c(all.genes[! is.na(all.genes)]))
  #names(hap.genes)[sapply(hap.genes, function(x) length(which(x %in% all.genes)) > 1)]
  #all.genes = all.genes[! all.genes %in% hap.genes]
  
  pass = all[all$gene %in% all.genes, ]
  pass = pass %>% 
    distinct(gene, symbol, .keep_all = T)
  pass = sapply(unique(pass$gene), function(x) unique(pass$symbol[pass$gene == x]))
  
  m1 = meta
  m1 = m1[match(colnames(exp.in)[2:ncol(exp.in)], m1$samp.et), ]
  
  
  snp.cors = data.frame()
  vst = readRDS("/media/d1/KG_Outputs/Bld_out_KG/matqtl/inputs/vst.rds")
  
  
  dir.create(paste0("/home/kgoldmann/Documents/PEAC_eqtl/Results/Blood/", covariates, "/locuszoom/"), showWarnings = FALSE)
  #invisible(file.remove(list.files(paste0("/home/kgoldmann/Documents/PEAC_eqtl/Results/Blood/", covariates, "/locuszoom/"), full.names = T)))
  for(k in 1:length(pass)){
    gene = pass[k]
    print(paste0(gene, "   ", k, "/", length(pass)))
    
    # check which genes.df within the plot neighbourhood for plotting
    genes.df = fread(gene.pos)
    gp = genes.df$start[genes.df$gene_id == names(pass)[k] ] + (genes.df$end[genes.df$gene_id ==  names(pass)[k]  ] - genes.df$start[genes.df$gene_id == names(pass)[k] ])/2
    gene.range = c(gp-window/2, gp+window/2)
    chrom = genes.df$chrom[genes.df$gene_id ==  names(pass)[k] ]
    g.start = genes.df$start[genes.df$gene_id ==  names(pass)[k] ]
    g.end = genes.df$end[genes.df$gene_id ==  names(pass)[k] ]
    
    genes.df = genes.df[genes.df$start >= gene.range[1] & genes.df$end <= gene.range[2] & genes.df$chrom == genes.df$chrom[genes.df$gene_id %in%  names(pass)[k] ] & genes.df$symbol != "", ]
    genes.df=genes.df[order(genes.df$start), ]
    
    plots=local_plotter(sym=gene, chr=chrom, width=window, file=file.in, col = "symbol", SNP_file_name = SNP_file, 
                        ld.calc=TRUE, genes.plot=genes.df$gene_id, highlight.gene = T, g.start=g.start, g.end=g.end, 
                        exp = exp.in, meta=m1, vst=vst)
    
    saveRDS(plots, paste0("/home/kgoldmann/Documents/PEAC_eqtl/Results/Blood/", covariates, "/locuszoom/Blood_", gene, "_", names(pass)[k], ".rds"))
    
    # Plot if snps comply with distribution rules
    if(! "reason" %in% names(plots)){
      #snp.cors = rbind(snp.cors, plots$snp.corr)
      
      p1 = ggarrange(plots$snps, plots$legend, plots$genes, nrow=2, ncol=2, align="v", heights = c(0.6, 0.4), widths=c(0.6, 0.15))#, nrow=2, heights=c(0.1, 0.9))
      p2 = ggarrange(plots$bp, plots$sp, nrow=2, ncol=1, heights=c(0.5, 1))
      
      p = ggarrange(p1, p2, nrow=1, ncol=2, widths=c(0.75, 0.3))
      
      ggsave(print(annotate_figure(p, top = paste("SNPs correlated with", gene))), width=12,
             filename =paste0("/home/kgoldmann/Documents/PEAC_eqtl/Results/Blood/", covariates, "/locuszoom/Blood_", gene, "_", names(pass)[k], ".pdf"), 
             device = cairo_pdf)
      
      
      # pdf(paste0("/home/kgoldmann/Documents/PEAC_eqtl/Results/Blood/", covariates, "/locuszoom/Blood_", gene, "_", names(pass)[k], ".pdf"), onefile=FALSE, width=10)
      # print(annotate_figure(p, top = paste("SNPs correlated with", gene)))
      # dev.off() 
      
      
      
    }
  }
  
  write.csv(snp.cors[snp.cors$Significant == T, ], paste0("/home/kgoldmann/Documents/PEAC_eqtl/Results/Blood/", covariates, "/snp_clin_corrs.csv"), quote=F, row.names = F)
  keep(local_plotter, sig, SNP_file, gene.pos, exp.in , window=1e6, cvs, sure=T)
}



# sudo cp -r /home/kgoldmann/Documents/PEAC_eqtl/Results/* ~/Desktop/


