output.creater = function(cv.df, results.dir, tissue, cv.names, rs.plotter){
  
  suppressMessages(require(MatrixEQTL))
  suppressMessages(require(ggrepel))
  suppressMessages(require(ggplot2))
  suppressMessages(require(dplyr))
  suppressMessages(require(RColorBrewer))
  suppressMessages(require(ggpubr))
  suppressMessages(require(data.table))
  suppressMessages(require(gdata))
  
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
  
  # Output file name
  covars.info = cv.names
  dir.create(paste0("/home/kgoldmann/Documents/PEAC_eqtl/Results/", tissue), showWarnings = F)
  abr = "none"
  if(any(grepl("EV", covars.info)) & ! any(grepl("PEER", covars.info))  ) abr = paste0(c("EV", gsub("EV", "", covars.info[grepl("EV", covars.info)])), collapse="")
  if(any(grepl("PEER", covars.info)) & ! any(grepl("EV", covars.info))  ) abr = paste0(c("PEER", gsub("PEER", "", covars.info[grepl("PEER", covars.info)])), collapse="")
  if(any(grepl("PEER", covars.info)) & any(grepl("EV", covars.info))  )  abr = paste0(c("EV", gsub("EV", "", covars.info[grepl("EV", covars.info)]), "_", "PEER", gsub("PEER", "", covars.info[grepl("PEER", covars.info)])), collapse="")
  
  output_file_name_cis = paste0("/home/kgoldmann/Documents/PEAC_eqtl/Results/", tissue, "/", 
                                abr, "/", tissue, "_", abr, ".csv")
  
  
  pvOutputThreshold_cis = 1 # Only associations significant at this level will be saved
  #pvOutputThreshold_tra = 0; ## only cis-eqtl
  errorCovariance = numeric(); # Error covariance matrix # Set to numeric() for identity.
  cisDist = 5e5; # Distance for local gene-SNP pairs
  
  ## Load genotype data
  print(paste("##", strsplit(as.character(Sys.time()), "\\s+")[[1]][2], ": SNP organiser (1 or 2 mins)"))
  snps = SlicedData$new();
  snps$fileDelimiter = " ";      # the ' ' character
  snps$fileOmitCharacters = "NA"; # denote missing values;
  snps$fileSkipRows = 1;          # one row of column labels
  snps$fileSkipColumns = 1;       # one column of row labels
  snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  suppressMessages(snps$LoadFile(SNP_file_name))
  
  # Calculate the MAF for each snp
  print(paste("##", strsplit(as.character(Sys.time()), "\\s+")[[1]][2], ": Calculate MAF (seconds)"))
  maf.list = vector('list', length(snps))
  for(sl in 1:length(snps)) {
    slice = snps[[sl]];
    maf.list[[sl]] = rowMeans(slice,na.rm=TRUE)/2;
    maf.list[[sl]] = pmin(maf.list[[sl]],1-maf.list[[sl]]);
  }
  maf = setNames(unlist(maf.list), unlist(snps$rowNameSlices))
  
  ## Load gene expression data
  print(paste("##", strsplit(as.character(Sys.time()), "\\s+")[[1]][2], ": Gene organiser (seconds)"))
  gene = SlicedData$new();
  gene$fileDelimiter = " ";      # the  character
  gene$fileOmitCharacters = "NA"; # denote missing values;
  gene$fileSkipRows = 1;          # one row of column labels
  gene$fileSkipColumns = 1;       # one column of row labels
  gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  suppressMessages(gene$LoadFile(expression_file_name));
  
  ## Load covariates
  #write.table(cv.df, "/home/kgoldmann/Documents/tmp_cv.txt")
  print(paste("##", strsplit(as.character(Sys.time()), "\\s+")[[1]][2], ": Covariate organiser (1 or 2 mins)"))
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
  
  print(paste("##", strsplit(as.character(Sys.time()), "\\s+")[[1]][2], ": matrix EQTL (exp around 3 mins)"))
  
  me = suppressMessages(Matrix_eQTL_main(
    snps = snps, gene = gene, cvrt = cvrt, output_file_name = NULL, ## for trans-eqtl
    pvOutputThreshold = 0, ## only cis-eqtl
    useModel = useModel, errorCovariance = errorCovariance, verbose = TRUE,
    output_file_name.cis = NULL, pvOutputThreshold.cis = pvOutputThreshold_cis,
    snpspos = snpspos2, genepos = genepos2, cisDist = cisDist,
    pvalue.hist = "qqplot", min.pv.by.genesnp = FALSE, noFDRsaveMemory = FALSE));
  
  
  print(paste("##",  strsplit(as.character(Sys.time()), "\\s+")[[1]][2], ": Match gene/snp info (exp around 5 mins)"))
  out <- me$cis$eqtls
  out$rs.id = snpspos$rs.id[match(out$snps, snpspos$snp)]
  out$maf = maf[match(out$snps, names(maf))]  # couple of seconds
  out$chr = gsub("\\:.*", "", out$snps)
  out = cbind(out, genepos[match(out$gene, genepos$gene_id), c("symbol", "start", "end", "gene_description")])
  out$pos = gsub('.*:([0-9]+):.*', '\\1', out$snps)
  
  print(paste("##", strsplit(as.character(Sys.time()), "\\s+")[[1]][2], ": Output data - started (exp < 1 min)"))
  dir.create(paste0("/home/kgoldmann/Documents/PEAC_eqtl/Results/", tissue, "/", abr), showWarnings = F)
 
  fwrite(list(paste("#", tissue), paste("\n#", paste(covars.info, collapse=",")), 
              paste("\n#", "p <", pvOutputThreshold_cis, "#")), 
         row.names=F, file=output_file_name_cis, col.names=F, append=F, quote=F)
  fwrite(out, file=output_file_name_cis, append=T, col.names = T)
  
  df_summary = data.frame("Cutoff" = c("those with p < 5e-8", "those with p < 5e-8 and maf > 0.05", "unique snps with p < 5e-8 and maf > 0.05", "unique genes with p < 5e-8 and maf > 0.05", "", ""), 
                          "Number passing" = c(length(out$snps[out$pvalue < 5e-8]), length(out$snps[out$pvalue < 5e-8 & out$maf > 0.05]), length(unique(out$snps[out$pvalue < 5e-8 & out$maf > 0.05])), length(unique(out$gene[out$pvalue < 5e-8 & out$maf > 0.05])), "", "" ))
  
  temp = table(out$symbol[out$pvalue < 5e-8])
  temp = data.frame(sort(temp[temp > 1], decreasing=T))


  fwrite(df_summary, file = gsub(".csv", "_summary.csv", output_file_name_cis), row.names = FALSE)
  fwrite(temp, file = gsub(".csv", "_summary.csv", output_file_name_cis), row.names = FALSE, append=T)
  #fwrite(out[out$FDR < 0.05, ], file=gsub(".csv", "_fdr_0.05.csv", output_file_name_cis), append=F, col.names = T)
  #fwrite(out[out$pvalue < 1e-5, ], file=gsub(".csv", "_supp_p_1e-5.csv", output_file_name_cis), append=F, col.names = T)
}