# Script to output locuszoom plots for a given gene

# sym=gene
# chr=chrom
# width=window
# file=file.in
# col = "symbol"
# SNP_file_name = SNP_file
# ld.calc=TRUE
# genes.plot=genes.df$gene_id
# highlight.gene = T
# g.start=g.start
# g.end=g.end
# exp = exp.in
# meta=m1

sig <- function(j, pc, snp.exp, clin) {   
  # p-val fn
  print(j)
  print(pc)
  df = data.frame("x"=as.numeric(snp.exp[pc, ]), "y"=clin[, j])
  df = df[! is.na(df$x), ]
  df = df[! is.na(df$y), ]
  if(nrow(df) > 2){
    mod <- lm(df$x ~ df$y)
    out <- if_else(clin[, j] %>% is.numeric, summary(mod)$coef[2L, 4L], anova(mod)[1L, 5L])
  } else{ out = NA}
  out
}

# R.df <- expand.grid(Clinical = colnames(meta)[which(! apply(meta, 2, function(x) length(levels(factor(x)))) %in% c(0, 1, nrow(meta)))][6], SNP = rownames(top10)) %>%
#   rowwise(.) %>%
#   dplyr::mutate(Association = sig(Clinical, SNP, top10, meta)) %>%   # Populate
#   dplyr::mutate(Padj = p.adjust(Association, method = "holm")) %>%
#   ungroup(.) %>%
#   dplyr::mutate(Significant = if_else(Association <= 0.05, TRUE, FALSE), Association = -log(Association))
# R.df$SNP = gsub("X", "", R.df$SNP)


library(RcmdrMisc)
library(dplyr)
library(varhandle)
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
                         vst=NULL,
                         ...){
  
  # Load in the snp toptable (pvalues with the gene sym)
  mat.df = data.frame(fread(cmd=paste0("grep -e ", names(sym), " -e statistic " , " ",  file)))
  mat.df$rs.id[mat.df$rs.id == chr] = mat.df$snps[mat.df$rs.id == chr]
  
  # Define the centre of the gene for window purposes
  gene.centre = mat.df$g_start[1] + (mat.df$g_end[1] - mat.df$g_start[1])/2
  
  # Load in the snp "expression data"
  if(length(unique(mat.df$snps)) > 7000){ # grepl falls over if too many snps so just split in half
    in.get = list
    snp.df1 = data.frame(fread(cmd=paste0("grep -e ", paste0(unique(mat.df$snps)[1:3500], collapse=" -e "), " ",  SNP_file_name)))
    snp.df2 = data.frame(fread(cmd=paste0("grep -e ", paste0(unique(mat.df$snps)[3501:length(unique(mat.df$snps))], collapse=" -e "), " ",  SNP_file_name)))
    snp.df = rbind(snp.df1, snp.df2)
  } else{ snp.df = data.frame(fread(cmd=paste0("grep -e ", paste0(unique(mat.df$snps)[1:7000], collapse=" -e "), " ",  SNP_file_name))) } 
  rownames(snp.df) = snp.df[,1]
  snp.df = data.frame(snp.df[, 2:ncol(snp.df)])
  
  
  mat.df = mat.df[mat.df$snps %in% rownames(snp.df), ]
  snp.df = snp.df[match(mat.df$snps, rownames(snp.df)), ] # Reorder by pvalue
  
  
  if(all(mat.df$maf < 0.05)){reason="none reaching maf cutoff"}
  snp.count = mat.df[mat.df$maf > 0.05, ]
  snp.count = snp.count[order(snp.count$pvalue), ]
  snp.id = snp.count$snps[1]
  if(! is.na(snp.id) & all(mat.df$pvalue[mat.df$snps == snp.id] > 5e-8)){reason="none reaching genome wide significance"}
  
  
  # So should we plot? Only if after snp removal the top one is still significanct
  if(! is.na(snp.id) & any(mat.df$pvalue[mat.df$snps == snp.id] <= 5e-8)){
    
    # Extract expression of the top snp
    top.snp = as.numeric(snp.df[snp.id, ]) 
    gene.exp = as.numeric(as.character(exp[which(exp$V1 == unique(mat.df$gene[mat.df$symbol == sym])), 2:ncol(exp)]))
    
    dat = data.frame(x =top.snp, y=gene.exp, id= colnames(exp)[ 2:ncol(exp)], vst = vst[names(sym), match(colnames(exp)[ 2:ncol(exp)], colnames(vst))])
    
    dat = dat[! is.na(dat$x), ]
    factor.match = setNames(c(paste0(strsplit(snp.id, ":")[[1]][3], "/", strsplit(snp.id, ":")[[1]][3]), 
                              paste0(strsplit(snp.id, ":")[[1]][3], "/", strsplit(snp.id, ":")[[1]][4]), 
                              paste0(strsplit(snp.id, ":")[[1]][4], "/", strsplit(snp.id, ":")[[1]][4])), c(0, 1, 2))
    dat$snps = factor(dat$x, labels=factor.match[match(levels(factor(dat$x)), names(factor.match))])
    
    table.count <- ggtexttable(data.frame(table(dat$snps)), theme = ttheme("mOrange"), rows=NULL, cols = c("Var", "Count"))
    
    # Calc the LD values
    #snp.df = snp.df[rowSums(is.na(snp.df)) != ncol(snp.df), ]
    if(ld.calc == TRUE){
      mat.df$ld = unlist(lapply(1:nrow(mat.df), function(x) {
        cor(as.numeric(  snp.df[rownames(snp.df) == mat.df$snps[x], ]), top.snp, use="complete.obs")^2
      }))
      mat.df$ld[is.na(mat.df$ld)] = 0
    } else{mat.df$ld = 1}
    
    snp.box = ggplot(dat[! is.na(dat$vst), ], aes(x = snps, y=as.numeric(vst), color=snps, fill=snps)) + 
      geom_boxplot(alpha=0.3, outlier.shape=NA) + 
      geom_jitter(width=0.25) + 
      theme_classic() + 
      theme(legend.position = "none", axis.text=element_text(colour="black"), aspect.ratio=1) + 
      labs(x=unique(mat.df$rs.id[mat.df$snps == snp.id]), y=sym, title="VST") 
    
    # Lets look at the top snps where ld != 1 and see how they correlate with clinical parameters
    top10 = mat.df$snps[mat.df$ld != 1 | mat.df$snps == snp.id]
    top10 = snp.df[rownames(snp.df) %in% top10[1:10], ]
    if(! all(mat.df$rs.id[match(rownames(top10), mat.df$snps)] == mat.df$rs.id[match(rownames(top10), mat.df$snps)][1])){
      rownames(top10) = make.names(mat.df$rs.id[match(rownames(top10), mat.df$snps)], unique=T)
    }
    
    num.cols <- which(apply(meta, 2, function(x) all(check.numeric(x[! is.na(x)]) == T) ))
    meta[num.cols] <- sapply(meta[num.cols], as.numeric)
    meta[! num.cols] <- sapply(meta[! num.cols], function(x) droplevels(factor(x)))
    
    rownames(meta) = colnames(top10)
    
    R.df <- expand.grid(Clinical = colnames(meta)[which(! apply(meta, 2, function(x) length(levels(factor(x)))) %in% c(0, 1, nrow(meta)))], SNP = rownames(top10)) %>%
      rowwise(.) #%>%
    
    R.df$Association = NA
    for(i in 1:nrow(R.df)){
      df = data.frame("x"=as.numeric(top10[ as.character(R.df$SNP[i]), ]), "y"=meta[, as.character(R.df$Clinical[i])] )
      df = df[! is.na(df$x), ]
      df = df[! is.na(df$y), ]
      if(nrow(df) > 2){
        mod <- lm(df$x ~ df$y)
        out <- if_else(meta[, as.character(R.df$Clinical[i])] %>% is.numeric, summary(mod)$coef[2L, 4L], anova(mod)[1L, 5L])
      } else{ out = NA}
      R.df$Association[i] = out
    }
    R.df = R.df %>%  # Populate
      dplyr::mutate(Padj = p.adjust(Association, method = "holm")) %>%
      ungroup(.) %>%
      dplyr::mutate(Significant = if_else(Association <= 0.05, TRUE, FALSE), Association = -log(Association))
    R.df$SNP = gsub("X", "", R.df$SNP)
    
    # Create a correlation plot for the snps vs clinical params
    sig.plot = ggplot(R.df, aes(SNP, Clinical, fill = Association, text = Association, color=Significant)) +
      geom_tile(size = 1L, width = 0.9, height = 0.9) + 
      scale_fill_gradientn(colors = c('#440154FF','#33638DFF', '#95D840FF', '#FDE725FF'), #c('white',  'dodgerblue1', 'dodgerblue3', 'dodgerblue4'), #, values=c(0, 0.1, -log10(0.05),  8), 
                           name = "-log10(P)", limits = c(0,8), labels=format(c(0, 0.1, 0.5, -log10(0.05),  8), digits=2)) +
      scale_color_manual(values = c('grey90', 'black')) +
      guides(color = FALSE) +
      labs(title = "", x = '', y="") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 315, hjust=0), legend.position="right", 
            plot.margin = unit(c(0,0,0,0), "cm"),  legend.key.width = unit(0.2, "cm"), legend.title = element_text(size = 8), 
            axis.text=element_text(colour="black"))
    R.df$gene = sym
    
    # Discretize the color scheme
    lz_scale = c("#322A82", "#80C2E9", "#64BB48", "#FBA525", "#EF3122")
    lz_labs = c(paste0("(", (1:4/5)-0.2, ", ", (1:4/5), "]"), paste0("(", (5/5)-0.2, ", ", (5/5), ")"))
    mat.df$ld_scale = NA
    for(br in 1:5) {
      max_break = (1:5/5)[br]
      min_break = max_break - 1/5
      mat.df$ld_scale[mat.df$ld >= min_break & mat.df$ld < max_break] = lz_labs[br]
    }
    mat.df$ld_scale[mat.df$ld == 1] = "(0.8, 1)"
    
    le = ggplot(mat.df, aes(x=pos, y=-log10(pvalue), fill=ld_scale)) + 
      geom_point(colour="black", shape=21, size=2.7) + xlab("") + 
      theme_classic() +
      geom_text_repel(data=mat.df[mat.df$snps == snp.id, ], color="black", aes(x=pos, y=-log10(pvalue), label=rs.id)) + 
      scale_fill_manual(values=lz_scale[lz_labs %in% unique(mat.df$ld_scale)], labels=lz_labs[lz_labs %in% unique(mat.df$ld_scale)]) + 
      labs(color = expression(paste(r^{2}))) + 
      theme(plot.title = element_text(hjust = 0.5), legend.position = "none", axis.text=element_text(colour="black"))
    
    # Create legend to mirror the locuszoom style
    dummy_plot = ggplot(data.frame("x"=lz_scale), aes(x, fill=x)) + geom_bar() + scale_fill_manual(values=rev(lz_scale), labels=rev(lz_labs), name=expression(paste(r^{2})))
    dummy_legend = as_ggplot(get_legend(dummy_plot))
    
    if(highlight.gene == T){le = le + annotate("rect", xmin=g.start, xmax=g.end, ymin=-Inf, ymax=Inf, alpha=0.3, fill="salmon")}
    
    #retrieve information of the gene of interest (longest protein coding)
    TX = transcripts(EnsDb.Hsapiens.v75, filter = ~ (gene_id %in% genes.plot & seq_name == chr))# & tx_biotype == "protein_coding" & gene_biotype == "protein_coding"))
    TX = data.frame(TX)
    lengths = data.frame("symbol"=TX$gene_id, "width"=TX$width, "length"=(TX$tx_cds_seq_end - TX$tx_cds_seq_start), "start"=TX$tx_cds_seq_start, "end"=TX$tx_cds_seq_end, "id"=TX$tx_id)
    lengths = lengths[! is.na(lengths$width), ]
    l = list()
    for(i in unique(lengths$symbol)){
      l[[i]] = lengths$id[as.numeric(as.character(lengths$width)) == max(lengths$width[lengths$symbol == i], na.rm=T)][1]
    }
    p.txdb = suppressMessages(autoplot(EnsDb.Hsapiens.v75, ~ (gene_id %in% genes.plot & seq_name == chr & tx_id %in% unlist(l) ), names.expr="gene_name"))
    
    g.plot =  attr(p.txdb, 'ggplot') +  theme_classic() + theme(axis.line.y=element_blank(), axis.text=element_text(colour="black"))
    if(highlight.gene == T){g.plot =g.plot + annotate("rect", xmin=g.start, xmax=g.end, ymin=-Inf, ymax=Inf, alpha=0.3, fill="salmon")}
    
    plot.range = c(min(c(layer_scales(le)$x$range$range[1], layer_scales(g.plot)$x$range$range[1])), 
                   max(c(layer_scales(le)$x$range$range[2], layer_scales(g.plot)$x$range$range[2])))
    
    
    return(list("snps" = le + xlim(plot.range), "genes"=g.plot + xlim(plot.range), "bp"=snp.box, snp.range=layer_scales(le)$x$range$range, "table"=table.count, "sp"=sig.plot, "snp.corr"=R.df, 
                "legend"=dummy_legend))
  } else{
    print(paste("For", sym, "no appropriate significant snps found"))
    return(list("reason"=reason))}
}
