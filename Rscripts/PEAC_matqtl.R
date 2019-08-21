library(MatrixEQTL)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggpubr)

#' cis-eqtl analysis using Matrixeqtl
#' http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/runit.html
#' calling Matrix_eQTL main function

gg.manhattan <- function(df, threshold, hlight, col, ylims, title){
  # format df
  df.tmp <- df %>% 
    
    # Compute chromosome size
    group_by(CHR) %>% 
    summarise(chr_len=max(BP)) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(df, ., by=c("CHR"="CHR")) %>%
    
    # Add a cumulative position of each SNP
    arrange(CHR, BP) %>%
    mutate( BPcum=BP+tot) %>%
    
    # Add highlight and annotation information
    mutate( is_highlight=ifelse(SNP %in% hlight, "yes", "no")) %>%
    mutate( is_annotate=ifelse(P < threshold, "yes", "no"))
  
  # get chromosome center positions for x-axis
  axisdf <- df.tmp %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  
  ggplot(df.tmp, aes(x=BPcum, y=-log10(P))) +
    # Show all points
    geom_point(aes(color=as.factor(CHR)), alpha=0.8, size=2) +
    scale_color_manual(values = rep(col, 22 )) +
    
    # custom X axis:
    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0), limits = ylims) + # expand=c(0,0)removes space between plot area and x axis 
    
    # add plot and axis titles
    ggtitle(paste0(title)) +
    labs(x = "Chromosome") +
    
    # add genome-wide sig and sugg lines
    geom_hline(yintercept = -log10(sig)) +
    geom_hline(yintercept = -log10(sugg), linetype="dashed") +
    
    # Add highlighted points
    #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
    
    # Add label using ggrepel to avoid overlapping
    #geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(SNP), alpha=0.7), size=5, force=1.3) +
    
    # Custom the theme:
    theme_bw(base_size = 22) +
    theme( 
      plot.title = element_text(hjust = 0.5, size=14),
      legend.position="none",
      panel.border = element_blank(),
      axis.title.x = element_text(size=12),
      axis.title.y = element_text(size=12),
      axis.text.x = element_text(size=10),
      axis.text.y = element_text(size=10),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
}

dir.create("/home/kgoldmann/Documents/PEAC_eqtl/Results", showWarnings = F)

output.creater = function(cv.df, results.dir, tissue, cv.names, rs.plotter){
  useModel = modelLINEAR
  
  # snp data
  SNP_file_name = paste0(results.dir, "/matqtl/inputs/genotype.txt")
  snps_location_file_name = paste0(results.dir, "/matqtl/inputs/snp_location.txt")
  #snps = fread(cmd=paste("head", SNP_file_name))
  
  # Gene expression file name
  expression_file_name = paste0(results.dir, "/matqtl/inputs/gene_expression_cqn.txt")
  gene_location_file_name = paste0(results.dir, "/matqtl/inputs/gene_location.txt")
  #gene = fread(cmd=paste("head", expression_file_name))
  
  # Covariate matrix
  if (nrow(cv.df) != 0){
    cv.df = cv.df[cv.names, ]
    cv.df = as.matrix(cv.df)
    cv.df = t(apply(cv.df, 1, as.numeric))
  }
  
  
  # df = data.frame("expression"=colnames(gene)[colnames(gene) != "V1"],
  #                 "snps"=colnames(snps)[colnames(snps) != "V1"],
  #                 "covars"=colnames(cv.df))
  # 
  # load("~/Documents/gcpeac/PEAC/PEAC main 081216.rdata")
  # df$conv = gsub("b", "", metadata$SampleID..QMUL.ID.only.[match(df$expression, metadata$SampleID..QMUL.or.Genentech.)])
  # # Sanity check!!
  # if(all(identical(as.character(df$conv), as.character(df$covars)),
  #     identical(as.character(df$snps), as.character(df$covars))) == FALSE) print("ERROR: check alignment")
  
  # Output file name
  covars.info = cv.names
  dir.create(paste0("/home/kgoldmann/Documents/PEAC_eqtl/Results/", tissue), showWarnings = F)
  dir.create(paste0("/home/kgoldmann/Documents/PEAC_eqtl/Results/", tissue, "/", covars.info[1], 
                    "-", covars.info[length(covars.info)]), showWarnings = F)
  output_file_name_cis = paste0("/home/kgoldmann/Documents/PEAC_eqtl/Results/", tissue, "/", 
                                covars.info[1], "-", covars.info[length(covars.info)], "/", 
                                tissue, "_", covars.info[1], "-", covars.info[length(covars.info)], 
                                ".csv")
  
  
  pvOutputThreshold_cis = 1 # Only associations significant at this level will be saved
  #pvOutputThreshold_tra = 0; ## only cis-eqtl
  errorCovariance = numeric(); # Error covariance matrix # Set to numeric() for identity.
  cisDist = 5e5; # Distance for local gene-SNP pairs
  
  ## Load genotype data
  print("## SNP organiser")
  snps = SlicedData$new();
  snps$fileDelimiter = " ";      # the ' ' character
  snps$fileOmitCharacters = "NA"; # denote missing values;
  snps$fileSkipRows = 1;          # one row of column labels
  snps$fileSkipColumns = 1;       # one column of row labels
  snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  snps$LoadFile(SNP_file_name)
  
  ## Load gene expression data
  print("## Gene organiser")
  gene = SlicedData$new();
  gene$fileDelimiter = " ";      # the  character
  gene$fileOmitCharacters = "NA"; # denote missing values;
  gene$fileSkipRows = 1;          # one row of column labels
  gene$fileSkipColumns = 1;       # one column of row labels
  gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  gene$LoadFile(expression_file_name);
  
  ## Load covariates
  #write.table(cv.df, "/home/kgoldmann/Documents/tmp_cv.txt")
  print("## Covariate organiser")
  if(nrow(cv.df) != 0) {
    cvrt = SlicedData$new();
    cvrt$CreateFromMatrix(cv.df)
  } else{
    cvrt = SlicedData$new()
    cvrt$fileDelimiter = " ";      # the TAB character
    cvrt$fileOmitCharacters = "NA"; # denote missing values;
    cvrt$fileSkipRows = 0;          # one row of column labels
    cvrt$fileSkipColumns = 1;       # one column of row labels
  }
  
  ## Run the analysis
  snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
  snpspos2 = snpspos[, 1:3]
  
  genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);
  genepos2 = genepos[, 1:4]
  
  print("## matrix EQTL")
  me = Matrix_eQTL_main(
    snps = snps, gene = gene, cvrt = cvrt, output_file_name = NULL, ## for trans-eqtl
    pvOutputThreshold = 0, ## only cis-eqtl
    useModel = useModel, errorCovariance = errorCovariance, verbose = TRUE,
    output_file_name.cis = NULL, pvOutputThreshold.cis = pvOutputThreshold_cis,
    snpspos = snpspos2, genepos = genepos2, cisDist = cisDist,
    pvalue.hist = "qqplot", min.pv.by.genesnp = FALSE, noFDRsaveMemory = FALSE);
  
  out <- me$cis$eqtls
  out$rs.id = snpspos$rs.id[match(out$snps, snpspos$snp)]
  out$symbol = genepos$symbol[match(out$gene, genepos$gene_id)]
  out$g_start = genepos$start[match(out$gene, genepos$gene_id)]
  out$g_end = genepos$end[match(out$gene, genepos$gene_id)]
  out$chr = gsub("\\:.*", "", out$snps)
  out$pos = as.numeric(as.character(unlist(lapply(out$snps, function(x) unlist(lapply(strsplit(as.character(x), split=":"), "[[", 2))))))
  
  print("## Output")
  fwrite(list(paste("#", tissue), paste("\n#", paste(covars.info, collapse=",")), 
              paste("\n#", "p <", pvOutputThreshold_cis, "#")), 
              row.names=F, file=output_file_name_cis, col.names=F, append=F, quote=F)
  fwrite(out, file=output_file_name_cis, append=T, col.names = T)
  
  
  # # Do a local explorer
  # if(length(rs.plotter) > 0){
  #   central.snp = rs.plotter
  # } else{central.snp = out$rs.id[out$pvalue == min(out$pvalue, na.rm=T)]}
  # window = 1e6
  # 
  # all.snps = out[out$pos <= unique(out$pos[out$rs.id == central.snp]) + window & 
  #                  out$pos >= unique(out$pos[out$rs.id == central.snp]) - window & 
  #                  out$chr == unique(out$chr[out$rs.id == central.snp]), ]
  # all.snps$pvalue = as.numeric(as.character(all.snps$pvalue))
  # 
  # # Calc the LD values
  # temp = fread(cmd=paste0("grep -e ", paste0(all.snps$snps, collapse=" -e "), " ",  SNP_file_name))
  # temp = data.frame(temp)
  # rownames(temp) = temp[,1]
  # temp = temp[, 2:ncol(temp)]
  # all.snps$ld = NA
  # for(i in 1:nrow(all.snps)){
  #   c = cor(as.numeric(temp[as.character(all.snps$snps[i]), ]), 
  #           as.numeric(temp[unique(as.character(out$snps[out$rs.id == central.snp])), ]), 
  #           use="complete.obs")
  #   all.snps$ld[i] = c
  # }
  # 
  # le = ggplot(all.snps, aes(x=pos, y=-log10(pvalue), color=ld)) + 
  #   geom_point() + xlab("") + 
  #   theme_classic() +
  #   geom_text_repel(data=all.snps[-log10(all.snps$pvalue) > 8, ], color="black",
  #                   aes(x=pos, y=-log10(pvalue), label=paste0(rs.id, " vs ", symbol))) + 
  #   scale_color_gradientn(colours = rainbow(5))
  # 
  # p0 = all.snps[! duplicated(all.snps$gene), ]
  # p0 = p0[p0$symbol != "", ]
  # p0 = p0[order(as.numeric(as.character(p0$g_start))), ]
  # 
  # key = setNames(unique(p0$symbol), rep(1:(length(unique(p0$symbol))/2), 2))
  # 
  # p1 = p0[, colnames(p0)[colnames(p0) != "g_end"]]
  # p2 = p0[, colnames(p0)[colnames(p0) != "g_start"]]
  # colnames(p1)[colnames(p1) == "g_start"] = colnames(p2)[colnames(p2) == "g_end"] = "position"
  # plot.df = rbind(p1, p2)
  # 
  # plot.df$symbol = as.character(plot.df$symbol)
  # 
  # labs = data.frame("gene"=unique(plot.df$symbol))
  # labs$xpos = NA
  # labs = labs[! grepl("ENSG", as.character(labs$gene)), ]
  # for(i in 1:nrow(labs)){
  #   labs$xpos[i] = min(plot.df$position[plot.df$symbol == labs$gene[i]], na.rm=T) + (max(plot.df$position[plot.df$symbol == labs$gene[i]], na.rm=T) - min(plot.df$position[plot.df$symbol == labs$gene[i]], na.rm=T))/2
  # }
  # 
  # labs = labs[order(as.numeric(as.character(labs$xpos))), ]
  # plot.df$y = as.numeric(as.character(names(key)[match(plot.df$symbol, key)]))
  # 
  # plot.df = plot.df[! is.na(plot.df$position) & ! is.na(plot.df$y), ]
  # 
  # genes = ggplot(plot.df, aes(x=position, y=y, group=symbol)) + 
  #   geom_line(color="red", size=5) + 
  #   annotate("text", x=as.numeric(as.character(labs$xpos)),
  #            y=0.2+rep(1:(length(labs$xpos)/2), 3)[1:nrow(labs)], 
  #            label=as.character(labs$gene)) + 
  #   theme(panel.background=element_blank(), axis.line.x = element_line(color="black", size = 0.5), 
  #         axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank())
  # 
  # dims = c(min(plot.df$position, na.rm=T), max(plot.df$position, na.rm=T))
  # 
  # f = list.files(gsub(paste0(tissue, "_", covars.info[1], "-", covars.info[length(covars.info)], ".csv"), "", 
  #                     output_file_name_cis))
  # f = f[! grepl("manhattan", f) & grepl("pdf", f)]
  # 
  # pdf(gsub(".csv", paste0("_", length(f) + 1, ".pdf"), output_file_name_cis))
  # print(ggarrange(le + xlim(dims), genes + xlim(dims), nrow=2, ncol=1, heights=c(1, 0.5), 
  #                 align="v", legend="bottom", common.legend=T))
  # invisible(dev.off())
  # 
  # # Create a manhattan plots
  # man.df = out[out$pvalue < 0.05, ]
  # man.df$POS = unlist(lapply(man.df$snps, function(x) unlist(lapply(strsplit(as.character(x), split=":"), "[[", 2))))
  # 
  # man.df$chr = as.numeric(as.character(man.df$chr))
  # man.df$pos = as.numeric(as.character(man.df$pos))
  # man.df$pvalue = as.numeric(as.character(man.df$pvalue))
  # 
  # colnames(man.df) = c("snps", "gene", "statistic", "P", "FDR", "beta",
  #                      "SNP", "symbol", "g_start", "g_end", "CHR", "BP", "POS" )
  # 
  # man.df$CHR = as.numeric(as.character(man.df$CHR))
  # 
  # 
  # mp = gg.manhattan(man.df, hlight=c(""), threshold=1e-6, col=mypalette, ylims=c(0,10), title=paste(c(tissue, ":", paste(covars.info, collapse=", ")), collapse=" "))
  # 
  # pdf(gsub(".csv", "_locuszoom.pdf", output_file_name_cis))
  # print(mp)
  # invisible(dev.off())
  
}

###############
# Synovium
###############
# Covariates file name
covariates_mat =  read.table(paste0("/media/d1/Syn_out_KG/matqtl/inputs/PCA4.PEER4.txt"))
rownames(covariates_mat) = covariates_mat[, "rn"]
covariates_mat = as.matrix(covariates_mat[, colnames(covariates_mat) != "rn"]) 

meta = read.table("/home/kgoldmann/Documents/PEAC_eqtl/Data/PEAC/PEAC_eth.txt")
covariates_mat = rbind(covariates_mat, t(meta[match(colnames(covariates_mat), meta$vcf_id), c("Ethnicity", "Gender")]))
covariates_mat[c("Ethnicity", "Gender"),] = apply(covariates_mat[c("Ethnicity", "Gender"),], 1, function(x) factor(x, labels=1:length(levels(factor(x)))))
covariates_mat2 = t(apply(covariates_mat, 1, as.numeric))
dimnames(covariates_mat2) = dimnames(covariates_mat)
covariates_mat = rbind(covariates_mat2, "none"=1)

mypalette <- c("#E2709A", "#CB4577", "#BD215B", "#970F42", "#75002B")
sig = 5e-8 # significant threshold line
sugg = 1e-6

opts = list(c(paste0("EV", 1:4), paste0("PEER", 1:4)), c("none"), paste0("EV", 1:4), paste0("PEER", 1:4))

for(i in 1:2){
  cv = as.character(opts[[i]])
  print(c("###", cv))
  cm = covariates_mat
  if (length(cv) == 1) cm = matrix(ncol = ncol(covariates_mat), nrow=0)
  output.creater(cv.df=cm, results.dir="/media/d1/Syn_out_KG/", tissue="Synovium", 
                 cv.names=cv, rs.plotter=c())
  
}




###############
# Blood
###############
# Covariates file name
covariates_mat =  read.table(paste0("/home/kgoldmann/Documents/PEAC_eqtl/Outputs_Blood/matqtl/inputs/PCA4.PEER4.txt"))
rownames(covariates_mat) = covariates_mat[, "rn"]
covariates_mat = as.matrix(covariates_mat[, colnames(covariates_mat) != "rn"]) 

meta = read.table("/home/kgoldmann/Documents/PEAC_eqtl/Data/PEAC/PEAC_eth.txt")
covariates_mat = rbind(covariates_mat, t(meta[match(colnames(covariates_mat), meta$vcf_id), c("Ethnicity", "Gender")]))
covariates_mat[c("Ethnicity", "Gender"),] = apply(covariates_mat[c("Ethnicity", "Gender"),], 1, function(x) factor(x, labels=1:length(levels(factor(x)))))
covariates_mat2 = t(apply(covariates_mat, 1, as.numeric))
dimnames(covariates_mat2) = dimnames(covariates_mat)

opts = list(c(paste0("EV", 1:4), paste0("PEER", 1:4)), c("none"), paste0("EV", 1:4), paste0("PEER", 1:4))
for(i in 1:length(opts)){
  cv = as.character(opts[[i]])
  print(c("###", cv))
  cm = covariates_mat
  if (length(cv) == 1) cm = matrix(ncol = ncol(covariates_mat), nrow=0)  
  output.creater(cv.df=cm, results.dir="/home/kgoldmann/Documents/PEAC_eqtl/Outputs_Blood/", 
                 tissue="Blood", cv.names=cv, rs.plotter=c())
  
}





###############
# Sanity test
###############
#lets investigate one of the top hits in blood with covariates as the PCAs

snps = read.table(SNP_file_name)
snp = snps["6:32194854:G:A", ]

genes = read.table(expression_file_name)
gene = genes["ENSG00000198502", ]

df = data.frame("snp"=as.numeric(as.character(snp)), "gene"=as.numeric(as.character(gene)))


ggplot(df, aes(x=snp, y=gene)) + geom_point()
cor.test(x=df$snp, y=df$gene, use="complete.obs")      
