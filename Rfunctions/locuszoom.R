# Script to output locuszoom plots for a given gene

# sym=gene
# chr=chrom
# width=window
# file=file.in
# col = "symbol"
# SNP_file_name = SNP_file
# ld.calc=TRUE
# genes.plot=genes.df$symbol
# highlight.gene = T
# g.start=g.start
# g.end=g.end
# exp = exp.in
# meta=m1

sig <- function(j, pc, snp.exp, clin) {   
  # p-val fn
  block = NULL
  if (block %>% is.null || j %in% unblock || j == block) {
    mod <- lm(as.numeric(snp.exp[pc, ]) ~ clin[, j])
  } else {
    mod <- lm(as.numeric(snp.exp[pc, ]) ~ clin[, j] + meta[block, ])
  }
  if_else(clin[, j] %>% is.numeric, summary(mod)$coef[2L, 4L], anova(mod)[1L, 5L])
}


library(RcmdrMisc)
local_plotter = function(sym, 
                         chr=NULL, 
                         width=1e6, 
                         file=file.in, 
                         col = "symbol", 
                         SNP_file_name = SNP_file, 
                         ld.calc=FALSE, 
                         genes.plot=c("CHURC1"), 
                         highlight.gene=F, 
                         g.start=NULL, 
                         g.end=NULL, 
                         exp=NULL, 
                         meta=NULL,
                         ...){
  
  # Load in the snp toptable (pvalues with the gene sym)
  mat.df = data.frame(fread(cmd=paste0("grep -e ", sym, " -e statistic " , " ",  file)))
  
  # Define the centre of the gene for window purposes
  gene.centre = mat.df$g_start[1] + (mat.df$g_end[1] - mat.df$g_start[1])/2
  
  # Load in the snp "expression data"
  if(length(unique(mat.df$snps)) > 7000){ # grepl falls over if too many snps so just split in half
    snp.df1 = data.frame(fread(cmd=paste0("grep -e ", paste0(unique(mat.df$snps)[1:3500], collapse=" -e "), " ",  SNP_file_name)))
    snp.df2 = data.frame(fread(cmd=paste0("grep -e ", paste0(unique(mat.df$snps)[3501:length(unique(mat.df$snps))], collapse=" -e "), " ",  SNP_file_name)))
    snp.df = rbind(snp.df1, snp.df2)
  } else{ snp.df = data.frame(fread(cmd=paste0("grep -e ", paste0(unique(mat.df$snps)[1:7000], collapse=" -e "), " ",  SNP_file_name))) } 
  rownames(snp.df) = snp.df[,1]
  snp.df = data.frame(snp.df[, 2:ncol(snp.df)])
  snp.df = snp.df[match(mat.df$snps, rownames(snp.df)), ] # Reorder by pvalue
  
  # Calculate the distribution of for each snp
  l = (apply(snp.df, 1, function(x) table(factor(x, levels=c(0:2)), exclude="none")))
  snp.count = data.frame()
  for(x in 1:length(l)) {
    i = l[[x]]
    i = data.frame(i)[! is.na(data.frame(i)$Var1), ]
    snp.count = rbind(snp.count, t(data.frame(i$Freq)))  
  }
  rownames(snp.count) = names(l)
  
  # Extract the most significant snp (which do not just have two levels OR do not contain only 1 individual in the high risk group)
  # remove snps only two levels
  if(any(snp.count[1, ]== 0)){ reason = "Only two levels"}
  snp.count = snp.count[apply(snp.count, 1, function(x) all(x != 0)), ]
  # remove snps with only 1 person in the high risk group
  if(snp.count[1, 3] == 1){ reason = "Only one person in high risk group"}
  snp.count = snp.count[apply(snp.count, 1, function(x) x[3] != 1), ]
  # The top significant snp, NA if none sig.    
  snp.id = rownames(snp.count)[1] 
  if(mat.df$pvalue[mat.df$snps == snp.id] > 1e-8) {reason = "None significant"}
  
  
  
  # So should we plot? Only if after snp removal the top one is still significanct
  if(mat.df$pvalue[mat.df$snps == snp.id] <= 1e-8){
    
    # Extract expression of the top snp
    top.snp = as.numeric(snp.df[snp.id, ]) 
    gene.exp = as.numeric(as.character(exp[which(exp$V1 == unique(mat.df$gene[mat.df$symbol == sym])), 2:ncol(exp)]))
    
    dat = data.frame(x =top.snp, y=gene.exp)
    dat = dat[! is.na(dat$x), ]
    factor.match = setNames(c(paste0(strsplit(snp.id, ":")[[1]][3], "/", strsplit(snp.id, ":")[[1]][3]), 
                              paste0(strsplit(snp.id, ":")[[1]][3], "/", strsplit(snp.id, ":")[[1]][4]), 
                              paste0(strsplit(snp.id, ":")[[1]][4], "/", strsplit(snp.id, ":")[[1]][4])), c(0, 1, 2))
    dat$snps = factor(dat$x, labels=factor.match[match(levels(factor(dat$x)), names(factor.match))])
    
    table.count <- ggtexttable(data.frame(table(dat$snps)), theme = ttheme("mOrange"), rows=NULL, cols = c("Var", "Count"))
    
    # Calc the LD values
    if(ld.calc == TRUE){
      mat.df$ld = unlist(lapply(1:nrow(mat.df), function(x) {
        cor(as.numeric(  snp.df[rownames(snp.df) == mat.df$snps[x], ]), top.snp, use="complete.obs")^2
      }))
      mat.df$ld[is.na(mat.df$ld)] = 0
    } else{mat.df$ld = 1}
    
    snp.box = ggplot(dat, aes(x = snps, y=y, color=snps, fill=snps)) + 
      geom_boxplot(alpha=0.3, outlier.shape=NA) + 
      geom_jitter(width=0.25) + 
      theme_classic() + 
      theme(legend.position = "none") + 
      labs(x=unique(mat.df$rs.id[mat.df$snps == snp.id]), y=sym)
    
    # Lets look at the top snps where ld != 1 and see how they correlate with clinical parameters
    top10 = mat.df$snps[mat.df$ld != 1 | mat.df$snps == snp.id]
    top10 = snp.df[rownames(snp.df) %in% top10[1:10], ]
    if(! all(mat.df$rs.id[match(rownames(top10), mat.df$snps)] == mat.df$rs.id[match(rownames(top10), mat.df$snps)][1])){
      rownames(top10) = mat.df$rs.id[match(rownames(top10), mat.df$snps)]
    }
    meta = meta[, c("CCP", "RF", "HAQ", "CRP", "ESR", "VAS", "Tender", "Swollen", "DAS28", 
                    "CD3.max", "CD20.max", "CD68L.max", "CD68SL.max", "CD138.max", "DAS28.real.6M", "Delta.DAS", 
                    "Pathotype", "EULAR3response",  "erosionstatus")]
    colnames(meta) = c("CCP", "RF", "HAQ", "CRP", "CCP", "ESR", "VAS", "Tender", "Swollen", "DAS28", 
                       "CD3", "CD20", "CD68L", "CD68SL", "CD138", "Six month DAS28", "Delta DAS28", 
                       "Pathotype",  "EULAR3 Response",  "Erosion Status")

    num.cols <- which(apply(meta, 2, function(x) all(check.numeric(x) == T) ))
    meta[num.cols] <- sapply(meta[num.cols], as.numeric)
    meta[! num.cols] <- sapply(meta[! num.cols], function(x) droplevels(factor(x)))
    
    rownames(meta) = colnames(top10)
    block=NULL
    alpha=0.05
    R.df <- expand.grid(Clinical = colnames(meta), SNP = rownames(top10)) %>%
      rowwise(.) %>%
      mutate(Association = sig(Clinical, SNP, top10, meta)) %>%   # Populate
      mutate(Padj = p.adjust(Association, method = "holm")) %>%
      ungroup(.) %>%
      mutate(Significant = if_else(Association <= alpha, TRUE, FALSE), Association = -log(Association))
    
    # Create a correlation plot for the snps vs clinical params
    sig.plot = ggplot(R.df, aes(SNP, Clinical, fill = Association, text = Association, color=Significant)) +
      geom_tile(size = 1L, width = 0.9, height = 0.9) + 
      scale_fill_gradientn(colors = c('white', 'pink', 'orange', 'red', 'darkred'), values=c(0, 0.1, 0.5, -log10(0.05),  8), name = "-log10(P)", limits = c(0,8), labels=format(c(0, 0.1, 0.5, -log10(0.05),  8), digits=2)) +
      scale_color_manual(values = c('grey90', 'black')) +
      guides(color = FALSE) +
      labs(title = "", x = '', y="") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 315, hjust=0), legend.position="right", 
            plot.margin = unit(c(0,0,0,0), "cm"),  legend.key.width = unit(0.2, "cm"), legend.title = element_text(size = 8))
    R.df$gene = sym
    
    le = ggplot(mat.df, aes(x=pos, y=-log10(pvalue), color=ld)) + 
      geom_point() + xlab("") + 
      theme_classic() +
      geom_text_repel(data=mat.df[mat.df$snps == snp.id, ], color="black", aes(x=pos, y=-log10(pvalue), label=rs.id)) + 
      scale_color_continuous(type = "viridis") + 
      labs(color = expression(paste(r^{2}))) + 
      theme(plot.title = element_text(hjust = 0.5), legend.position = "none") 
    
    if(highlight.gene == T){le = le + annotate("rect", xmin=g.start, xmax=g.end, ymin=-Inf, ymax=Inf, alpha=0.3, fill="salmon")}
    
    #retrieve information of the gene of interest (longest protein coding)
    TX = transcripts(EnsDb.Hsapiens.v75, filter = ~ (symbol %in% genes.plot & seq_name == chr))# & tx_biotype == "protein_coding" & gene_biotype == "protein_coding"))
    TX = data.frame(TX)
    lengths = data.frame("symbol"=TX$symbol, "width"=TX$width, "length"=(TX$tx_cds_seq_end - TX$tx_cds_seq_start), "start"=TX$tx_cds_seq_start, "end"=TX$tx_cds_seq_end, "id"=TX$tx_id)
    lengths = lengths[! is.na(lengths$width), ]
    l = list()
    for(i in unique(lengths$symbol)){
      l[[i]] = lengths$id[as.numeric(as.character(lengths$width)) == max(lengths$width[lengths$symbol == i], na.rm=T)][1]
    }
    #   & tx_biotype == "protein_coding" & gene_biotype == "protein_coding"
    p.txdb = autoplot(EnsDb.Hsapiens.v75, ~ (symbol %in% genes.plot & seq_name == chr & tx_id %in% unlist(l)), names.expr="gene_name", ...)
    
    
    g.plot =  attr(p.txdb, 'ggplot') +  theme_classic() + theme(axis.line.y=element_blank())
    if(highlight.gene == T){g.plot =g.plot + annotate("rect", xmin=g.start, xmax=g.end, ymin=-Inf, ymax=Inf, alpha=0.3, fill="salmon")}
    
    plot.range = c(min(c(layer_scales(le)$x$range$range[1], layer_scales(g.plot)$x$range$range[1])), 
                   max(c(layer_scales(le)$x$range$range[2], layer_scales(g.plot)$x$range$range[2])))
    
    
    return(list("snps" = le + xlim(plot.range), "genes"=g.plot + xlim(plot.range), "bp"=snp.box, snp.range=layer_scales(le)$x$range$range, "table"=table.count, "sp"=sig.plot, "snp.corr"=R.df))
  } else{
    print(paste("For", sym, "no appropriate significant snps found"))
    return(list("reason"=reason))}
}