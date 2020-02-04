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

local_plotter_trait = function(sym, file=file.in, SNP_file_name = SNP_file, ld.calc=FALSE, exp=NULL, vst, ...){
  
  mat.df = data.frame(fread(cmd=paste0("grep -e '", sym, "' -e statistic " , " ",  file)))
  mat.df = mat.df[mat.df$Fantom.5.Module == sym, ]
  mat.df = mat.df[mat.df$maf >= 0.05, ]
  
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
  
  cols = unique(mat.df$Fantom.5.Module[mat.df$Fantom.5.Module == sym])
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
  
  
  le = ggman(mat.df, snp="snps", bp="pos", chrom="chr", pvalue="pvalue", pointSize = 2, title="") + theme_classic() +lims(y=c(2, max(-log10(mat.df$pvalue), na.rm=T)))
  le = ggmanLabel(le, labelDfm = mat.df[-log10(mat.df$pvalue) == max(-log10(mat.df$pvalue), na.rm=T), ], 
                  snp = "snps", label = "rs.id", type = "text", size = 4, color="black")
  
  return(list("snps" = le, snp.range=layer_scales(le)$x$range$range, "bp"=snp.box))
}


# Synovium
############
file.in = list.files("/home/kgoldmann/Documents/PEAC_eqtl/Results/Synovium", recursive=T, full.names = T) 
file.in = file.in[grepl("fantom5mods.csv", file.in)]

SNP_file = "/media/d1/KG_Outputs/Syn_out_KG/matqtl/inputs/genotype.txt"
exp.in = fread("/media/d1/KG_Outputs/Syn_out_KG/F5_expression.csv")
#window=1e6

for(file in file.in){
  all = fread(file, select ="Fantom 5 Module" )
  all.mods = unique(all$`Fantom 5 Module`)
  
  cv = strsplit(file, split="/")[[1]][8]
  dir.create(paste0("/home/kgoldmann/Documents/PEAC_eqtl/Results/Synovium/", cv, "/locuszoom_trait/"), showWarnings = FALSE)
  invisible(file.remove(list.files(paste0("/home/kgoldmann/Documents/PEAC_eqtl/Results/Synovium/", cv, "/locuszoom_trait/"), full.names = T)))
  for(mod in all.mods){
    print(paste(mod, which(all.mods == mod), "/", length(all.mods)))
    plots=local_plotter_trait(sym=mod, file=file.in, SNP_file_name = SNP_file, ld.calc=FALSE,  exp = exp.in)
    p1 = ggarrange(plots$snps, plots$bp, nrow=1, ncol=2, legend="none")
    
    pdf(paste0("/home/kgoldmann/Documents/PEAC_eqtl/Results/Synovium/", cv, "/locuszoom_trait/Synovium_", mod, ".pdf"), onefile=FALSE, width=10)
    print(annotate_figure(p1, top = paste("SNPs correlated with", mod)))
    dev.off() 
  }  
}





# Blood
############
file.in = list.files("/home/kgoldmann/Documents/PEAC_eqtl/Results/Blood", recursive=T, full.names = T) 
file.in = file.in[grepl("fantom5mods.csv", file.in)]
SNP_file = "/media/d1/KG_Outputs/Bld_out_KG/matqtl/inputs/genotype.txt"
exp.in = fread("/media/d1/KG_Outputs/Bld_out_KG/F5_expression.csv")
window=1e6



for(file in file.in){
  all = fread(file, nrows=50000)
  all.mods = unique(all$`Fantom 5 Module`)
  
  cv = strsplit(file, split="/")[[1]][8]
  dir.create(paste0("/home/kgoldmann/Documents/PEAC_eqtl/Results/Blood/", cv, "/locuszoom_trait/"), showWarnings = FALSE)
  invisible(file.remove(list.files(paste0("/home/kgoldmann/Documents/PEAC_eqtl/Results/Blood/", cv, "/locuszoom_trait/", full.names = T))))
  for(mod in all.mods){
    print(paste(mod, which(all.mods == mod), "/", length(all.mods)))
    plots=local_plotter_trait(sym=mod, file=file, SNP_file_name = SNP_file, ld.calc=FALSE,  exp = exp.in, pcutoff=2)
    p1 = ggarrange(plots$snps, plots$bp, nrow=1, ncol=2, legend="none")
    
    pdf(paste0("/home/kgoldmann/Documents/PEAC_eqtl/Results/Blood/", cv, "/locuszoom_trait/Blood_", mod, ".pdf"), onefile=FALSE, width=10)
    print(annotate_figure(p1, top = paste("SNPs correlated with", mod)))
    dev.off() 
  }  
}

# sudo cp -r /home/kgoldmann/Documents/PEAC_eqtl/Results/Blood/EV1234/* /home/kgoldmann/Documents/gcpeac/Katriona/PEAC_analysis/PEAC_eQTL/Results/Blood/

