#######################
# Look closer at the ZNF266 snps to see why some correlate with Batch
#######################
# run locuszoom analysis to get top10 for gene ZNF266
plot.list=list()
for(i in  1:nrow(top10)){
  df = data.frame(x=meta$Batch, y=as.numeric(top10[i, ]))
  meta$Batch = gsub("Genentech", "", meta$Batch)
  df = df[! is.na(df$x), ]
  
  
  plot.list[[i]] = ggplot(df, aes(x, y, fill=x)) + geom_boxplot() + geom_jitter(width=0.25) + theme_classic() + labs(y=rownames(top10)[i], x="")+ theme(legend.position = "none") + stat_compare_means()
  
}

pdf("/home/kgoldmann/Documents/PEAC_eqtl/Results/Specific_Analysis/ZNF266_correlation_check.pdf", width=10)
annotate_figure(ggarrange(plotlist = plot.list), top = "ZNF266 Batch correlation check, Synovium, EV1:4&PEER1")
dev.off()

#######################
# Check why HLA-V such negative expression
#######################

exp.in = fread("/media/d1/KG_Outputs/Syn_out_KG/matqtl/inputs/gene_expression_cqn.txt")

range(as.numeric(exp.in[exp.in$V1 == "ENSG00000181126", 2:ncol(exp.in)]), na.rm=T)

counts = fread("/media/d1/KG_Outputs/Syn_out_KG/RNA_counts/groups/Genentech.txt")

range(as.numeric(counts[counts$V1 == "ENSG00000181126", 2:ncol(counts)]))

pdf("/home/kgoldmann/Documents/PEAC_eqtl/Results/Specific_Analysis/HLA-V_counts.pdf")
hist(as.numeric(counts[counts$V1 == "ENSG00000181126", 2:ncol(counts)]), main="HLA-V Counts", xlab = "Counts", breaks=100)
hist(as.numeric(counts[counts$V1 == "ENSG00000196301", 2:ncol(counts)]), main="HLA-DRB9 Counts", xlab = "Counts", breaks=100)
hist(as.numeric(counts[counts$V1 == "ENSG00000164308", 2:ncol(counts)]), main="ERAP2 Counts", xlab = "Counts", breaks=100)
dev.off()

# Just has very low counts

#######################
# Drivers of the EV and PEER values with clinical params
#######################

# Drivers plot, adapted from bioplotr package to include not just pca vals
# Katriona Goldmann

library(gridExtra)
library(grid)
library(dplyr)
library(purrr)
library(ggplot2)
library(ggpubr)

source("/home/kgoldmann/Documents/peac/Rfunctions/covariate_drivers.R")

# Synovium
##############
covariates_mat =  read.table(paste0("/media/d1/KG_Outputs/Syn_out_KG/matqtl/inputs/PCA5.PEER5.txt"))
mp = read.table("/home/kgoldmann/Documents/PEAC_eqtl/Data/PEAC/PEAC_eth_syn.txt")

rownames(covariates_mat) = covariates_mat$rn
covariates_mat = covariates_mat[, colnames(covariates_mat) != "rn"]

meta = readRDS("/home/kgoldmann/Documents/gcpeac/PEAC/PEAC_Imputed_data.rds")
m = meta[match(colnames(covariates_mat), meta$Sample.final), ]
mp = mp[match(colnames(covariates_mat), mp$vcf_id), ]

all(identical(colnames(covariates_mat), as.character(m$Sample.final)), 
    identical(colnames(covariates_mat), as.character(mp$vcf_id)))
df.m = cbind(mp[, c( "Gender", "Batch")],
             m[, c("Ethnicity", "Pathotype", "Age", "CCP", "Inflammatory.score", "CRP", "ESR", "Tender", 
                   "Swollen", "VAS","DAS28.ESR", "DAS28.CRP", "CD3.max", "CD20.max", "CD68L.max", "CD68SL.max", "CD138.max" )])

colnames(df.m) = c("Gender", "Batch", "Ethnicity", "Pathotype", "Age", "CCP", "Inflammatory Score", "CRP", "ESR", "Tender", 
                   "Swollen", "VAS","DAS28 ESR", "DAS28 CRP", "CD3", "CD20", "CD68L", "CD68SL", "CD138" )

df.m$Ethnicity = tolower(df.m$Ethnicity)
df.m$Ethnicity[df.m$Ethnicity == "caucasian"] = "white"
df.m$Ethnicity[df.m$Ethnicity == "pakistani-"] = "pakistani"

syn = covar_drivers_specific(pca = t(covariates_mat), clin = df.m, alpha = 0.01)


pdf("/home/kgoldmann/Documents/PEAC_eqtl/Results/Specific_Analysis/drivers_syn_all.pdf", height=10)
syn$plot + theme(legend.position = "top", axis.text.x = element_text(angle = 315, hjust=0))
print(grid.arrange(tableGrob(syn$df[syn$df$Significant == TRUE, ], theme=ttheme_default(base_size=10, padding = unit(c(1, 1), "mm")))))
dev.off()

# Lets investigate this further
plot.list = list()
for(i in which(syn$df$Significant == T)){
  ev = as.character(syn$df$PC[i])
  var = as.character(syn$df$Feature[i])
  df = data.frame(e = as.numeric(covariates_mat[ev, ]), m = df.m[, var])
  if(class(df$m) == "factor"){
    p2 = ggplot(df, aes(x=m, y=e, color=m, fill=m)) + geom_boxplot(alpha=0.3, outlier.shape=NA) + geom_jitter(width=0.25) + labs(x=var, y=ev) + theme_classic() + theme(legend.position = "none")
  } else{
    p2 = ggplot(df, aes(x=m, y=e)) + geom_point() + labs(x=var, y=ev) + theme_classic()
  }
  plot.list[[length(plot.list)+1]] = p2
}
length(plot.list)

pdf("/home/kgoldmann/Documents/PEAC_eqtl/Results/Specific_Analysis/drivers_syn_closer.pdf", onefile = T)
ggarrange(plotlist = plot.list, nrow=2, ncol=2)
dev.off()


# Blood
##############
covariates_mat =  read.table(paste0("/home/kgoldmann/Documents/PEAC_eqtl/Outputs_Blood/matqtl/inputs/PCA10.PEER10.txt"))
mp = read.table("/home/kgoldmann/Documents/PEAC_eqtl/Data/PEAC/PEAC_eth_blood.txt")

rownames(covariates_mat) = covariates_mat$rn
covariates_mat = covariates_mat[, colnames(covariates_mat) != "rn"]

meta = readRDS("/home/kgoldmann/Documents/gcpeac/PEAC/PEAC_Imputed_data.rds")
m = meta[match(colnames(covariates_mat), meta$Sample.final), ]
mp = mp[match(colnames(covariates_mat), mp$vcf_id), ]

all(identical(colnames(covariates_mat), as.character(m$Sample.final)), 
    identical(colnames(covariates_mat), as.character(mp$vcf_id)))
df.m = cbind(mp[, c( "Gender")],
             m[, c("Ethnicity", "Pathotype", "Age", "CCP", "Inflammatory.score", "CRP", "ESR", "Tender", 
                   "Swollen", "VAS","DAS28.ESR", "DAS28.CRP", "CD3.max", "CD20.max", "CD68L.max", "CD68SL.max", "CD138.max" )])

colnames(df.m) = c("Gender", "Ethnicity", "Pathotype", "Age", "CCP", "Inflammatory Score", "CRP", "ESR", "Tender", 
                   "Swollen", "VAS","DAS28 ESR", "DAS28 CRP", "CD3", "CD20", "CD68L", "CD68SL", "CD138" )

bld = covar_drivers_specific(pca = t(covariates_mat), clin = df.m, alpha = 0.05)

pdf("/home/kgoldmann/Documents/PEAC_eqtl/Results/Specific_Analysis/drivers_bld_all.pdf", height=10)
bld$plot + theme(legend.position = "top", axis.text.x = element_text(angle = 315, hjust=0))
print(grid.arrange(tableGrob(bld$df[bld$df$Significant == TRUE, ], theme=ttheme_default(base_size=6, padding = unit(c(1, 1), "mm")))))
dev.off()


