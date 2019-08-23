# Data visualisations
# http://bioconductor.org/packages/release/bioc/vignettes/ggbio/inst/doc/ggbio.pdf

library(data.table)
library(ggplot2)
library(ggrepel)
library(ggbio)
library(ggpubr)
library(EnsDb.Hsapiens.v75)

data(genesymbol, package = "biovizBase")

p.ideo <- Ideogram(genome = "hg19", aspect.ratio=1/50, fill="green", color="blue")

local_plotter = function(sym, chr=NULL, width=1e6, file=file.in, col = "symbol", SNP_file_name = SNP_file, 
                         ld.calc=TRUE, genes.plot=c("CHURC1"), highlight.gene=F, g.start=NULL, g.end=NULL, exp=NULL, ...){

  mat.df = data.frame(fread(cmd=paste0("grep -e ", sym, " -e statistic " , " ",  file)))
  mat.df = mat.df[mat.df[, col] == sym, ]
  
  gene.centre = mat.df$g_start[1] + (mat.df$g_end[1] - mat.df$g_start[1])/2

  # Calc the LD values
  ld.df = data.frame(fread(cmd=paste0("grep -e ", paste0(unique(mat.df$snps), collapse=" -e "), " ",  SNP_file_name)))
  rownames(ld.df) = ld.df[,1]
  ld.df = ld.df[, 2:ncol(ld.df)]
  top.snp = as.numeric(ld.df[mat.df$snps[which(-log10(mat.df$pvalue) == max(-log10(mat.df$pvalue), na.rm=T))], ])
  snp.id = mat.df$snps[which(-log10(mat.df$pvalue) == max(-log10(mat.df$pvalue), na.rm=T))]
  if(ld.calc == TRUE){
    mat.df$ld = NA
    for(i in 1:nrow(mat.df)){
      c = cor(as.numeric(ld.df[mat.df$snps[i], ]), top.snp, use="complete.obs")
      mat.df$ld[i] = c
    }
  } else{mat.df$ld = 1}

  gene.exp = as.numeric(as.character(exp[which(exp$V1 == unique(mat.df$gene[mat.df$symbol == sym])), 2:ncol(exp)]))
  
  dat = data.frame(x =top.snp, y=gene.exp)
  dat = dat[! is.na(dat$x), ]
  factor.match = setNames(c(paste0(strsplit(snp.id, ":")[[1]][3], "/", strsplit(snp.id, ":")[[1]][3]), 
                   paste0(strsplit(snp.id, ":")[[1]][3], "/", strsplit(snp.id, ":")[[1]][4]), 
                   paste0(strsplit(snp.id, ":")[[1]][4], "/", strsplit(snp.id, ":")[[1]][4])), c(0, 1, 2))
  dat$snps = factor(dat$x, labels=factor.match[match(levels(factor(dat$x)), names(factor.match))])
  
  snp.box = ggplot(dat, aes(x = snps, y=y, color=snps, fill=snps)) + 
    geom_boxplot(alpha=0.3, outlier.shape=NA) + 
    geom_jitter(width=0.25) + 
    theme_classic() + 
    theme(legend.position = "none") + 
    labs(x=unique(mat.df$rs.id[mat.df$snps == snp.id]), y=sym)
  
  
  le = ggplot(mat.df, aes(x=pos, y=-log10(pvalue), color=ld)) + 
    geom_point() + xlab("") + 
    theme_classic() +
    geom_text_repel(data=mat.df[which(mat.df$pvalue == min(mat.df$pvalue, na.rm=T)), ], color="black",
                    aes(x=pos, y=-log10(pvalue), label=rs.id)) + 
    scale_color_continuous(type = "viridis") + 
    labs(color = "Linkage Disequilibrium") + 
    theme(plot.title = element_text(hjust = 0.5)) 
  
  if(highlight.gene == T){le = le + annotate("rect", xmin=g.start, xmax=g.end, ymin=-Inf, ymax=Inf, alpha=0.3, fill="salmon")}
  
  print("test1")
  #retrieve information of the gene of interest (longest protein coding)
  TX = transcripts(EnsDb.Hsapiens.v75, filter = ~ (symbol %in% genes.plot & seq_name == chr))# & tx_biotype == "protein_coding" & gene_biotype == "protein_coding"))
  lengths = data.frame("symbol"=TX$symbol, "length"=(TX$tx_cds_seq_end - TX$tx_cds_seq_start), "start"=TX$tx_cds_seq_start, "end"=TX$tx_cds_seq_end, "id"=TX$tx_id)
  lengths = lengths[! is.na(lengths$length), ]
  l = list()
  for(i in unique(lengths$symbol)){
    l[[i]] =lengths$id[as.numeric(as.character(lengths$length)) == max(lengths$length[lengths$symbol == i], na.rm=T)][1]
  }
  #   & tx_biotype == "protein_coding" & gene_biotype == "protein_coding"
  p.txdb = autoplot(EnsDb.Hsapiens.v75, ~ (symbol %in% genes.plot & seq_name == chr & tx_id %in% unlist(l)), names.expr="gene_name", ...)
  
 
  print("test2")
  g.plot =  attr(p.txdb, 'ggplot') +  theme_classic() + theme(axis.line.y=element_blank())
  if(highlight.gene == T){g.plot =g.plot + annotate("rect", xmin=g.start, xmax=g.end, ymin=-Inf, ymax=Inf, alpha=0.3, fill="salmon")}

  plot.range = c(min(c(layer_scales(le)$x$range$range[1], layer_scales(g.plot)$x$range$range[1])), 
                 max(c(layer_scales(le)$x$range$range[2], layer_scales(g.plot)$x$range$range[2])))

  return(list("snps" = le + xlim(plot.range), "genes"=g.plot + xlim(plot.range), "bp"=snp.box, snp.range=layer_scales(le)$x$range$range))
}


# Synovium
############
file.in = "/home/kgoldmann/Documents/PEAC_eqtl/Results/Synovium/EV1-PEER4/Synovium_EV1-PEER4.csv"
SNP_file = "/media/d1/Syn_out_KG/matqtl/inputs/genotype.txt"
gene.pos = "/media/d1/Syn_out_KG//matqtl/inputs/gene_location.txt"
exp.in = fread("/media/d1/Syn_out_KG/matqtl/inputs/gene_expression_cqn.txt")
window=1e6


all = fread(file.in)
all = all[all$pvalue < 0.001 & all$symbol != ""]
all$count = table(all$symbol)[match(all$symbol, names(table(all$symbol)))]
all = all[rev(order(all$count)), ]
common.genes = unique(all$symbol[all$count > 1])[1:20]
top.genes = as.character(data.frame(fread(file.in, select="symbol", nrows=100))$symbol)
top.genes = top.genes[top.genes != ""][1:20]
all.genes = unique(c(common.genes, top.genes))

invisible(file.remove(list.files("/home/kgoldmann/Documents/PEAC_eqtl/Results/Synovium/EV1-PEER4/locuszoom/", full.names = T)))
for(gene in all.genes){
  print(paste(gene, which(all.genes == gene), "/", length(all.genes)))
  # check which genes.df within the plot neighbourhood for plotting
  genes.df = fread(gene.pos)
  gp = genes.df$start[genes.df$symbol == gene] + (genes.df$end[genes.df$symbol == gene] - genes.df$start[genes.df$symbol == gene])/2
  gene.range = c(gp-window/2, gp+window/2)
  chrom = genes.df$chrom[genes.df$symbol == gene]
  g.start = genes.df$start[genes.df$symbol == gene]
  g.end = genes.df$end[genes.df$symbol == gene]
  #chr.plot = p.ideo + xlim(GRanges(paste0("chr", chrom), IRanges(gene.range[1], gene.range[2]))) + labs(title=paste("Snps correlated with", gene))
  #genes.df = genes.df[genes.df$symbol %in% select(EnsDb.Hsapiens.v75, keys=unique(as.character(genes.df$chrom[genes.df$symbol == gene])), columns=c("SYMBOL"), keytype="SEQNAME")$SYMBOL, ]
  genes.df = genes.df[genes.df$start >= gene.range[1] & genes.df$end <= gene.range[2] & genes.df$chrom == genes.df$chrom[genes.df$symbol == gene] & genes.df$symbol != "", ]
  genes.df=genes.df[order(genes.df$start), ]

  plots=local_plotter(sym=gene, chr=chrom, width=window, file=file.in, col = "symbol", SNP_file_name = SNP_file, 
                      ld.calc=TRUE, genes.plot=genes.df$symbol, highlight.gene = T, g.start=g.start, g.end=g.end, 
                      exp = exp.in)
  #p = ggarrange(attr(chr.plot, 'ggplot'), 
  p1 = ggarrange(plots$snps, plots$genes, nrow=2, ncol=1, align="v", legend="bottom", common.legend=T, heights = c(0.6, 0.4), widths=c(0.6, 0.4))#, nrow=2, heights=c(0.1, 0.9))
  p2 = ggarrange(plots$bp, nrow=2, ncol=1, heights=c(0.3, 0.7))
  p = ggarrange(p1, p2, nrow=1, ncol=2, widths=c(0.7, 0.3))
  
  pdf(paste0("/home/kgoldmann/Documents/PEAC_eqtl/Results/Synovium/EV1-PEER4/locuszoom/Synovium_", gene, ".pdf"), onefile=FALSE, width=10)
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

all = fread(file.in)
all = all[all$pvalue < 0.001 & all$symbol != ""]
all$count = table(all$symbol)[match(all$symbol, names(table(all$symbol)))]
all = all[all$chr != 6, ]
all = all[rev(order(all$count)), ]
common.genes = unique(all$symbol[all$count > 1])[1:20]
top.genes = as.character(data.frame(fread(file.in, select="symbol", nrows=100))$symbol)
top.genes = top.genes[top.genes != ""][1:20]
all.genes = unique(c(common.genes, top.genes))

for(gene in all.genes){
  print(gene)
  # check which genes within the plot neighbourhood for plotting
  genes.df = fread(gene.pos)
  gp = genes.df$start[genes.df$symbol == gene] + (genes.df$end[genes.df$symbol == gene] - genes.df$start[genes.df$symbol == gene])/2
  gene.range = c(gp-window/2, gp+window/2)
  chrom = genes.df$chrom[genes.df$symbol == gene]
  chr.plot = p.ideo + xlim(GRanges(paste0("chr",  chrom), IRanges(gene.range[1], gene.range[2]))) + labs(title=paste("Snps correlated with", gene))
  #genes.df = genes.df[genes.df$symbol %in% select(EnsDb.Hsapiens.v75, keys=unique(as.character(genes.df$chrom[genes.df$symbol == gene])), columns=c("SYMBOL"), keytype="SEQNAME")$SYMBOL, ]
  genes.df = genes.df[genes.df$start >= gene.range[1] & genes.df$end <= gene.range[2] & genes.df$chrom == genes.df$chrom[genes.df$symbol == gene] & genes.df$symbol != "", ]
  genes.df=genes.df[order(genes.df$start), ]
  
  plots=local_plotter(sym=gene, chr=chrom, width=window, file=file.in, col = "symbol", SNP_file_name = SNP_file, 
                      ld.calc=TRUE, genes.plot=genes.df$symbol, highlight.gene = T, g.start=g.start, g.end=g.end)
  p = ggarrange(attr(chr.plot, 'ggplot'), ggarrange(plots$snps, plots$genes,  nrow=2, ncol=1, align="v", legend="bottom", common.legend=T, heights = c(0.6, 0.4)), nrow=2, heights=c(0.1, 0.9))
  
  pdf(paste0("/home/kgoldmann/Documents/PEAC_eqtl/Results/Blood/EV1-PEER4/locuszoom/Blood_", gene, ".pdf"), onefile=FALSE)
  print(annotate_figure(p, top = paste("SNPs correlated with", gene)))
  dev.off()
}

# sudo cp -r /home/kgoldmann/Documents/PEAC_eqtl/Results/Blood/EV1-PEER4/* /home/kgoldmann/Documents/gcpeac/Katriona/PEAC_analysis/PEAC_eQTL/Results/Blood/

