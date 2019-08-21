# Data visualisations

library(data.table)
library(ggplot2)
library(ggrepel)
library(ggbio)
library(ggpubr)
library(EnsDb.Hsapiens.v75)


data(genesymbol, package = "biovizBase")

p.ideo <- Ideogram(genome = "hg19", aspect.ratio=1/50, fill="green", color="blue")

local_plotter = function(sym, width=1e6, file=file.in, col = "symbol", SNP_file_name = SNP_file, ld.calc=TRUE, genes.plot=c("CHURC1"), ...){
  print("test1")
  mat.df = data.frame(fread(cmd=paste0("grep -e ", sym, " -e statistic " , " ",  file)))
  mat.df = mat.df[mat.df[, col] == sym, ]
  
  gene.centre = mat.df$g_start[1] + (mat.df$g_end[1] - mat.df$g_start[1])/2
  print("test2")
  # Calc the LD values
  if(ld.calc == TRUE){
    ld.df = data.frame(fread(cmd=paste0("grep -e ", paste0(unique(mat.df$snps), collapse=" -e "), " ",  SNP_file_name)))
    rownames(ld.df) = ld.df[,1]
    ld.df = ld.df[, 2:ncol(ld.df)]
    mat.df$ld = NA
    for(i in 1:nrow(mat.df)){
      c = cor(as.numeric(ld.df[mat.df$snps[i], ]), 
              as.numeric(ld.df[mat.df$snps[which(-log10(mat.df$pvalue) == max(-log10(mat.df$pvalue), na.rm=T))], ]), 
              use="complete.obs")
      mat.df$ld[i] = c
    }
  } else{mat.df$ld = 1}
  print("test3")
  le = ggplot(mat.df, aes(x=pos, y=-log10(pvalue), color=ld)) + 
    geom_point() + xlab("") + 
    theme_classic() +
    geom_text_repel(data=mat.df[which(mat.df$pvalue == min(mat.df$pvalue, na.rm=T)), ], color="black",
                    aes(x=pos, y=-log10(pvalue), label=rs.id)) + 
    scale_color_continuous(type = "viridis") + 
    labs(color = "Linkage Disequilibrium") + 
    theme(plot.title = element_text(hjust = 0.5)) 
  
  print("test4")
  #retrieve information of the gene of interest
  p.txdb = autoplot(EnsDb.Hsapiens.v75, ~ (symbol %in% genes.plot), names.expr="gene_name", ...)
  print("test5")
  g.plot =  attr(p.txdb, 'ggplot') +  theme_classic() + theme(axis.line.y=element_blank())
  print("test6")
  plot.range = c(min(c(layer_scales(le)$x$range$range[1], layer_scales(g.plot)$x$range$range[1])), 
                 max(c(layer_scales(le)$x$range$range[2], layer_scales(g.plot)$x$range$range[2])))
  print("test7")
  return(list("snps" = le + xlim(plot.range), "genes"=g.plot + xlim(plot.range), snp.range=layer_scales(le)$x$range$range))
}


# Synovium
############
file.in = "/home/kgoldmann/Documents/PEAC_eqtl/Results/Synovium/EV1-PEER4/Synovium_EV1-PEER4.csv"
SNP_file = "/media/d1/Syn_out_KG/matqtl/inputs/genotype.txt"
gene.pos = "/media/d1/Syn_out_KG//matqtl/inputs/gene_location.txt"
window=1e6

top.genes = as.character(data.frame(fread(file.in, select="symbol", nrows=100))$symbol)
top.genes = top.genes[top.genes != ""][1:20]

for(gene in top.genes[! top.genes %in% c("HCG4P7", "NUPR1L", "PSPHP1")]){
  print(gene)
  # check which genes.df within the plot neighbourhood for plotting
  genes.df = fread(gene.pos)
  gp = genes.df$start[genes.df$symbol == gene] + (genes.df$end[genes.df$symbol == gene] - genes.df$start[genes.df$symbol == gene])/2
  gene.range = c(gp-window/2, gp+window/2)
  chr.plot = p.ideo + xlim(GRanges(paste0("chr", genes.df$chrom[genes.df$symbol == gene]), IRanges(gene.range[1], gene.range[2]))) + labs(title=paste("Snps correlated with", gene))
  genes.df = genes.df[genes.df$symbol %in% select(EnsDb.Hsapiens.v75, keys=unique(as.character(genes.df$chrom[genes.df$symbol == gene])), columns=c("SYMBOL"), keytype="SEQNAME")$SYMBOL, ]
  genes.df = genes.df[genes.df$start >= gene.range[1] & genes.df$end <= gene.range[2] & genes.df$chrom == genes.df$chrom[genes.df$symbol == gene] & genes.df$symbol != "", ]
  genes.df=genes.df[order(genes.df$start), ]

  plots=local_plotter(sym=gene, width=window, file=file.in, col = "symbol", SNP_file_name = SNP_file, ld.calc=TRUE, genes.plot=genes.df$symbol)
  p = ggarrange(attr(chr.plot, 'ggplot'), ggarrange(plots$snps, plots$genes,  nrow=2, ncol=1, align="v", legend="bottom", common.legend=T, heights = c(0.6, 0.4)), nrow=2, heights=c(0.1, 0.9))
  
  pdf(paste0("/home/kgoldmann/Documents/PEAC_eqtl/Results/Synovium/EV1-PEER4/locuszoom/Synovium_", gene, ".pdf"), onefile=FALSE)
  print(annotate_figure(p, top = paste("SNPs correlated with", gene)))
  dev.off() 

}  




# Blood
############
file.in = "/home/kgoldmann/Documents/PEAC_eqtl/Results/Blood/EV1-PEER4/Blood_EV1-PEER4.csv"
SNP_file = "/home/kgoldmann/Documents/PEAC_eqtl/Outputs_Blood/matqtl/inputs/genotype.txt"
gene.pos = "/home/kgoldmann/Documents/PEAC_eqtl/Outputs_Blood/matqtl/inputs/gene_location.txt"
dir.create("/home/kgoldmann/Documents/PEAC_eqtl/Results/Blood/EV1-PEER4/locuszoom/", showWarnings = F)
window=1e6

top.genes = as.character(data.frame(fread(file.in, select="symbol", nrows=100))$symbol)
top.genes = top.genes[top.genes != ""][1:20]

for(gene in top.genes[! top.genes %in% c("PSPHP1", "C4BPA")]){
  print(gene)
  # check which genes within the plot neighbourhood for plotting
  genes.df = fread(gene.pos)
  gp = genes.df$start[genes.df$symbol == gene] + (genes.df$end[genes.df$symbol == gene] - genes.df$start[genes.df$symbol == gene])/2
  gene.range = c(gp-window/2, gp+window/2)
  chr.plot = p.ideo + xlim(GRanges(paste0("chr", genes.df$chrom[genes.df$symbol == gene]), IRanges(gene.range[1], gene.range[2]))) + labs(title=paste("Snps correlated with", gene))
  genes.df = genes.df[genes.df$symbol %in% select(EnsDb.Hsapiens.v75, keys=unique(as.character(genes.df$chrom[genes.df$symbol == gene])), columns=c("SYMBOL"), keytype="SEQNAME")$SYMBOL, ]
  genes.df = genes.df[genes.df$start >= gene.range[1] & genes.df$end <= gene.range[2] & genes.df$chrom == genes.df$chrom[genes.df$symbol == gene] & genes.df$symbol != "", ]
  genes.df=genes.df[order(genes.df$start), ]
  
  plots = local_plotter(sym=gene, width=window, file=file.in, col = "symbol", SNP_file_name = SNP_file, ld.calc=TRUE, genes.plot=genes.df$symbol)
  p = ggarrange(attr(chr.plot, 'ggplot'), ggarrange(plots$snps, plots$genes,  nrow=2, ncol=1, align="v", legend="bottom", common.legend=T, heights = c(0.6, 0.4)), nrow=2, heights=c(0.1, 0.9))
  
  pdf(paste0("/home/kgoldmann/Documents/PEAC_eqtl/Results/Blood/EV1-PEER4/locuszoom/Blood_", gene, ".pdf"), onefile=FALSE)
  print(annotate_figure(p, top = paste("SNPs correlated with", gene)))
  dev.off()
}

