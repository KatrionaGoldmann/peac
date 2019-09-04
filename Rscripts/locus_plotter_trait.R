# Data visualisations
# http://bioconductor.org/packages/release/bioc/vignettes/ggbio/inst/doc/ggbio.pdf

library(data.table)
library(ggplot2)
library(ggrepel)
library(ggbio)
library(ggpubr)
#install_github("drveera/ggman")
library(ggman)

#data(genesymbol, package = "biovizBase")

#p.ideo <- Ideogram(genome = "hg19", aspect.ratio=1/50, fill="green", color="blue")

local_plotter_trait = function(sym, file=file.in, SNP_file_name = SNP_file, ld.calc=FALSE, exp=NULL, ...){

  mat.df = data.frame(fread(cmd=paste0("grep -e '", sym, "' -e statistic " , " ",  file)))
  mat.df = mat.df[mat.df$gene == sym, ]

  # Calc the LD values
  snp.id = mat.df$snps[which(-log10(mat.df$pvalue) == max(-log10(mat.df$pvalue), na.rm=T))[1]]
  ld.df = data.frame(fread(cmd=paste0("grep -e ", paste0(snp.id, collapse=" -e "), " ",  SNP_file_name)))
  rownames(ld.df) = ld.df[,1]
  ld.df = ld.df[, 2:ncol(ld.df)]
  top.snp = as.numeric(ld.df) #[mat.df$snps[which(-log10(mat.df$pvalue) == max(-log10(mat.df$pvalue), na.rm=T))], ])
  
  if(ld.calc == TRUE){
    mat.df$ld = NA
    for(i in 1:nrow(mat.df)){
      c = cor(as.numeric(ld.df[mat.df$snps[i], ]), top.snp, use="complete.obs")
      mat.df$ld[i] = c
    }
  } else{mat.df$ld = 1}

  cols = unique(mat.df$gene[mat.df$gene == sym])
  mod.exp = as.numeric(exp[, get(names(exp)[34])])

  dat = data.frame(x =top.snp, y=mod.exp)
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
  
  
  le = ggman(mat.df, snp="snps", bp="snp_pos", chrom="chr", pvalue="pvalue", pointSize = 2, title="") + theme_classic() 
  le = ggmanLabel(le, labelDfm = mat.df[-log10(mat.df$pvalue) == max(-log10(mat.df$pvalue), na.rm=T), ], snp = "snps", label = "rs.id", type = "text", size = 6)

  # le = ggplot(mat.df, aes(x=snp_pos, y=-log10(pvalue), color=ld)) + 
  #   geom_point() + xlab("") + 
  #   theme_classic() +
  #   geom_text_repel(data=mat.df[mat.df$snps == snp.id, ], color="black", aes(x=snp_pos, y=-log10(pvalue), label=rs.id)) + 
  #   scale_color_continuous(type = "viridis") + 
  #   labs(color = "Linkage Disequilibrium") + 
  #   theme(plot.title = element_text(hjust = 0.5), legend.position = "none", ...) 

  return(list("snps" = le, snp.range=layer_scales(le)$x$range$range, "bp"=snp.box))
}


# Synovium
############
file.in = "/home/kgoldmann/Documents/PEAC_eqtl/Results/Synovium/EV1-PEER4/Synovium_EV1-PEER4_fantom5mods.csv"
SNP_file = "/media/d1/Syn_out_KG/matqtl/inputs/genotype.txt"
exp.in = fread("/media/d1/Syn_out_KG/F5_expression.csv")
#window=1e6


all = fread(file.in, nrows=20000)
all.genes = unique(all$gene)

invisible(file.remove(list.files("/home/kgoldmann/Documents/PEAC_eqtl/Results/Synovium/EV1-PEER4/locuszoom_trait/", full.names = T)))
for(gene in all.genes){
  print(paste(gene, which(all.genes == gene), "/", length(all.genes)))


  plots=local_plotter_trait(sym=gene, file=file.in, SNP_file_name = SNP_file, ld.calc=FALSE,  exp = exp.in)
  #p = ggarrange(attr(chr.plot, 'ggplot'), 
  p1 = ggarrange(plots$snps, plots$bp, nrow=1, ncol=2, legend="none")

  pdf(paste0("/home/kgoldmann/Documents/PEAC_eqtl/Results/Synovium/EV1-PEER4/locuszoom_trait/Synovium_", gene, ".pdf"), onefile=FALSE, width=10)
  print(annotate_figure(p1, top = paste("SNPs correlated with", gene)))
  dev.off() 

}  




# Blood
############
file.in = "/home/kgoldmann/Documents/PEAC_eqtl/Results/Blood/EV1-PEER4/Blood_EV1-PEER4.csv"
SNP_file = "/home/kgoldmann/Documents/PEAC_eqtl/Outputs_Blood/matqtl/inputs/genotype.txt"
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

