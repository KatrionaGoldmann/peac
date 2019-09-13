library(MatrixEQTL)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggpubr)
library(gdata)
library(data.table)

dir.create("/home/kgoldmann/Documents/PEAC_eqtl/Results", showWarnings = F)

module.function <- function(modulegenes, new_exp){
  # check how many new_exp genes are in modulegenes
  checkgenes <- lapply(modulegenes, function(x) {
    colnames(new_exp)[colnames(new_exp) %in% x]
  })
  
  modulelengths <- sapply(checkgenes, length)
  # Remove those modules with 0 genes in new_exp
  modulegenes <- modulegenes[!(modulelengths==0)]
  namelist <- colnames(new_exp[,1:ncol(new_exp)-1])
  
  module.output3 <- sapply(modulegenes, function(x) {
    data <- as.matrix(new_exp[,colnames(new_exp) %in% x])  # Catagorises the data into gene groups
    colsd <- apply(data, 2, sd)
    if (sum(colsd==0)>0) print(colnames(data)[colsd==0])
    scaledata <- scale(data[, colsd!=0])
    
    #scaledata <- scaledata[,colSums(scaledata != 0) != 0]
    svd1 <- svd(t(scaledata), nu=0, nv=3)
    #svd1 <- svd(t(data), nu=0, nv=3)       ### I CHANGED THIS LINE FROM THE ABOVE LINE TO EXPRESS THE SCALED BLOOD BETTER.
    ret <- svd1$v[,1:3]
    sign1 <- sign(sum(cor(ret[,1],  data[,colsd!=0] )))
    if (sign1 != 0)  ret <- sign1* ret
    dimnames(ret) <- list(rownames(data), paste0("ALL", 1:3))
    ret[,1]
  })
}
output.creater = function(exp.df, cv.df, results.dir, tissue, cv.names, rs.plotter, pv=0.001){
  useModel = modelLINEAR
  
  # snp data
  SNP_file_name = paste0(results.dir, "/matqtl/inputs/genotype.txt")
  snps_location_file_name = paste0(results.dir, "/matqtl/inputs/snp_location.txt")

  
  # Covariate matrix
  if (nrow(cv.df) != 0){
    cv.df = cv.df[cv.names, ]
    cv.df = as.matrix(cv.df)
    cv.df = t(apply(cv.df, 1, as.numeric))
  }
  
  # Output file name
  covars.info = cv.names
  dir.create(paste0("/home/kgoldmann/Documents/PEAC_eqtl/Results/", tissue), showWarnings = F)
  dir.create(paste0("/home/kgoldmann/Documents/PEAC_eqtl/Results/", tissue, "/", covars.info[1], 
                    "-", covars.info[length(covars.info)]), showWarnings = F)
  output_file_name_cis = paste0("/home/kgoldmann/Documents/PEAC_eqtl/Results/", tissue, "/", 
                                covars.info[1], "-", covars.info[length(covars.info)], "/", 
                                tissue, "_", covars.info[1], "-", covars.info[length(covars.info)], "_p", pv,  
                                "_fantom5mods.csv")
  
  
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
  gene$CreateFromMatrix(exp.df);
  
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
  
  genepos = data.frame("gene_id" = rownames(exp.df), "chrom"=1, "start"=1001:(1000+nrow(exp.df)), "end" = 2001:(2000+nrow(exp.df)))
  genepos2 = genepos[, 1:4]
  
  print("## matrix EQTL")
  me = Matrix_eQTL_main(
    snps = snps, gene = gene, cvrt = cvrt, output_file_name = NULL, ## for trans-eqtl
    pvOutputThreshold = pv, ## inc trans cis-eqtl
    useModel = useModel, errorCovariance = errorCovariance, verbose = TRUE,
    output_file_name.cis = NULL, pvOutputThreshold.cis = pvOutputThreshold_cis,
    snpspos = snpspos2, genepos = genepos2, cisDist = 0,
    pvalue.hist = "qqplot", min.pv.by.genesnp = FALSE, noFDRsaveMemory = FALSE);
  
  out <- me$trans$eqtls
  out$rs.id = snpspos$rs.id[match(out$snps, snpspos$snp)]
  out$chr = snpspos$chr[match(out$snps, snpspos$snp)]
  out$snp_pos = as.numeric(as.character(unlist(lapply(out$snps, function(x) unlist(lapply(strsplit(as.character(x), split=":"), "[[", 2))))))
  
  if(pv == 1) saveRDS(out, paste0("/home/kgoldmann/Documents/PEAC_eqtl/Results/", tissue, ".rds"))
  
  print("## Output")
  fwrite(list(paste("#", tissue), paste("\n#", paste(covars.info, collapse=",")), 
              paste("\n#", "p <", pv, "#")), 
              row.names=F, file=output_file_name_cis, col.names=F, append=F, quote=F)
  fwrite(out, file=output_file_name_cis, append=T, col.names = T)
  
}

###############
# Synovium
###############
# Covariates file name
covariates_mat =  read.table(paste0("/media/d1/KG_Outputs/Syn_out_KG/matqtl/inputs/PCA4.PEER4.txt"))
rownames(covariates_mat) = covariates_mat[, "rn"]
covariates_mat = as.matrix(covariates_mat[, colnames(covariates_mat) != "rn"]) 

# Get the expression data in matrix format
exp.in = data.frame(fread("/media/d1/KG_Outputs/Syn_out_KG/matqtl/inputs/gene_expression_cqn.txt"))
rownames(exp.in) = exp.in$V1
exp.in = exp.in[, colnames(exp.in) != "V1"]
temp = data.frame(fread("/home/kgoldmann/Documents/PEAC_eqtl/Results/Synovium/EV1-PEER4/Synovium_EV1-PEER4.csv"))
exp.in$gn = temp$symbol[match(rownames(exp.in), temp$gene)]
exp.in = exp.in[exp.in$gn != "", ]
rownames(exp.in) = make.names(exp.in$gn, unique=TRUE)
exp.in = exp.in[, colnames(exp.in) != "gn"]

# Calculate the F5 scores
load("/home/kgoldmann/Documents/gcpeac/Sharmila/Fantom_5/fantom5_ambry ML 05.05.16.RData")
keep(exp.in, module.function, covariates_mat, output.creater, modulegenes, sure=T)
f5 = module.function(modulegenes, new_exp=t(exp.in))

exp.vals=f5[,c("plasma cells", "Natural Killer Cells", "CD14-CD16+ Monocytes", "CD14+ monocyte derived endothelial progenitor cells",
               "CD14+ Monocytes","CD14+CD16- Monocytes","CD14+CD16+ Monocytes","CD19+ B Cells","CD34+ Progenitors",
               "CD4+ T Cells","CD4+CD25-CD45RA- memory conventional T cells","CD4+CD25-CD45RA+ naive conventional T cells",
               "CD4+CD25+CD45RA- memory regulatory T cells","CD4+CD25+CD45RA+ naive regulatory T cells",
               "CD8+ T Cells")]

fwrite(f5, "/media/d1/KG_Outputs/Syn_out_KG/F5_expression.csv", row.names = T)

# Lets see how the PEERs and f5 mods correlate
covar_drivers<-function(pca, clin, index='Sample', block=NULL, unblock=NULL, kernel=NULL, 
                        kpar=NULL, top=NULL, n.pc=5L, label=FALSE, alpha=NULL,
                        p.adj=NULL, title='VariationByFeature', legend='right', hover=FALSE){
  
  # Preliminaries
  loc <- c('bottom', 'left', 'top', 'right','bottomright', 'bottomleft', 'topleft', 'topright')
  if (!legend %in% loc) {stop('legend must be one of ', stringify(loc, 'or'), '.')}
  
  tibble(Feature = colnames(clin)) %>% print(n = nrow(.))
  
  sig <- function(j, pc) {                       # p-val fn
    if (block %>% is.null || j %in% unblock || j == block) {
      mod <- lm(pca[, pc] ~ clin[[j]])
    } else {
      mod <- lm(pca[, pc] ~ clin[[j]] + clin[[block]])
    }
    if_else(clin[[j]] %>% is.numeric,
            summary(mod)$coef[2L, 4L], anova(mod)[1L, 5L])
  }
  
  df <- expand.grid(Feature = colnames(clin), PC = colnames(pca)) %>%
    rowwise(.) %>%
    #mutate(Association = sig(Feature, PC)) %>%   # Populate
    ungroup(.)
  
  df$Association=NA
  for(i in 1:nrow(df)){
    df$Association[i] = cor.test(clin[, df$Feature[i]], pca[, df$PC[i]])$p.value
  }
  
  if (!p.adj %>% is.null) {
    df <- df %>% mutate(Association = p.adjust(Association, method = p.adj))
  }
  alpha = 0.05
  df <- df %>% 
    mutate(Significant = if_else(Association <= alpha, TRUE, FALSE),
           Association = -log(Association))
  
  # Build plot
  if (!p.adj %>% is.null && p.adj %in% c('fdr', 'BH', 'BY')) {
    leg_lab <- expression(~-log(italic(q)))
  } else {
    leg_lab <- expression(~-log(italic(p)))
  }
  p <- ggplot(df, aes(PC, Feature, fill = Association, text = Association, color = Significant)) +
    geom_tile(size = 1L, width = 0.9, height = 0.9) +
    #coord_equal() +
    scale_fill_gradientn(colors = c('white', 'pink', 'orange', 'red', 'darkred'),
                         name = leg_lab) +
    scale_color_manual(values = c('grey90', 'black')) +
    #scale_x_discrete(labels = paste0(unique(df$PC), pve)) +
    guides(color = FALSE) +
    labs(title = title, x = 'Principal Component') +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  if (label) {
    p <- p + geom_text(aes(label = round(Association, 2L)))
  }
  
  # Output
  list("plot"=p, "df" = df)
}
syn = covar_drivers(clin = exp.vals, pca = t(covariates_mat))
pdf("/home/kgoldmann/Documents/PEAC_eqtl/Results/mod_drivers_syn.pdf")
syn$plot
dev.off()
# So lets remove PEER 2->4 from the trait analysis


cm = covariates_mat[1:5, ]
output.creater(exp.df = t(exp.vals), cv.df=cm, results.dir="/media/d1/KG_Outputs/Syn_out_KG/", tissue="Synovium", 
                 cv.names=c(paste0("EV", 1:4), paste0("PEER", 1)), rs.plotter=c())
  

###############
# Blood
###############
# Covariates file name
covariates_mat =  read.table(paste0("/home/kgoldmann/Documents/PEAC_eqtl/Outputs_Blood/matqtl/inputs/PCA4.PEER4.txt"))
rownames(covariates_mat) = covariates_mat[, "rn"]
covariates_mat = as.matrix(covariates_mat[, colnames(covariates_mat) != "rn"]) 

# Get the expression data in matrix format
exp.in = data.frame(fread("/home/kgoldmann/Documents/PEAC_eqtl/Outputs_Blood/matqtl/inputs/gene_expression_cqn.txt"))
rownames(exp.in) = exp.in$V1
exp.in = exp.in[, colnames(exp.in) != "V1"]
temp = data.frame(fread("/home/kgoldmann/Documents/PEAC_eqtl/Results/Blood/EV1-PEER4/Blood_EV1-PEER4.csv"))
exp.in$gn = temp$symbol[match(rownames(exp.in), temp$gene)]
exp.in = exp.in[exp.in$gn != "", ]
rownames(exp.in) = make.names(exp.in$gn, unique=TRUE)
exp.in = exp.in[, colnames(exp.in) != "gn"]

# Calculate the F5 scores
load("/home/kgoldmann/Documents/gcpeac/Sharmila/Fantom_5/fantom5_ambry ML 05.05.16.RData")
keep(exp.in, module.function, covariates_mat, output.creater, modulegenes, sure=T)
f5 = module.function(modulegenes, new_exp=t(exp.in))

exp.vals=f5[,c("plasma cells", "Natural Killer Cells", "CD14-CD16+ Monocytes", "CD14+ monocyte derived endothelial progenitor cells",
               "CD14+ Monocytes","CD14+CD16- Monocytes","CD14+CD16+ Monocytes","CD19+ B Cells","CD34+ Progenitors",
               "CD4+ T Cells","CD4+CD25-CD45RA- memory conventional T cells","CD4+CD25-CD45RA+ naive conventional T cells",
               "CD4+CD25+CD45RA- memory regulatory T cells","CD4+CD25+CD45RA+ naive regulatory T cells",
               "CD8+ T Cells")]
fwrite(f5, "/home/kgoldmann/Documents/PEAC_eqtl/Outputs_Blood/F5_expression.csv", row.names = T)

cm = covariates_mat
output.creater(exp.df = t(exp.vals), cv.df=cm, results.dir="/home/kgoldmann/Documents/PEAC_eqtl/Outputs_Blood/", tissue="Blood", 
               cv.names=c(paste0("EV", 1:4), paste0("PEER", 1:4)), rs.plotter=c(), pv=0.001)
