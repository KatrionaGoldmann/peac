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
