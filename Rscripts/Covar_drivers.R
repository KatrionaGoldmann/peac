# Drivers plot, adapted from bioplotr package to include not just pca vals
# Katriona Goldmann

library(gridExtra)
library(grid)
library(dplyr)
library(purrr)
library(ggplot2)

covar_drivers<-function(pca, clin, index='Sample', block=NULL, unblock=NULL, kernel=NULL, kpar=NULL, top=NULL, n.pc=5L, label=FALSE, alpha=NULL,
                        p.adj=NULL, title='VariationByFeature', legend='right', hover=FALSE){
  
  # Preliminaries
  if (!alpha %>% is.null) {if (alpha <= 0L || alpha >= 1L) {stop('alpha must be numeric on (0, 1).') }} else { alpha <- 0L}
  if (!p.adj %>% is.null) {
    p_adj <- c('holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr')
    if (!p.adj %in% p_adj) {stop('p.adj must be one of ', stringify(p_adj, 'or'), '. See ?p.adjust.') }
  }
  loc <- c('bottom', 'left', 'top', 'right','bottomright', 'bottomleft', 'topleft', 'topright')
  if (!legend %in% loc) {stop('legend must be one of ', stringify(loc, 'or'), '.')}
  
  tibble(Feature = colnames(clin), Class = clin %>% map_chr(class)) %>% print(n = nrow(.))
  
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
    mutate(Association = sig(Feature, PC)) %>%   # Populate
    ungroup(.)
  if (!p.adj %>% is.null) {
    df <- df %>% mutate(Association = p.adjust(Association, method = p.adj))
  }
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
    coord_equal() +
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





##############
# Synovium
##############
covariates_mat =  read.table(paste0("/media/d1/Syn_out_KG/matqtl/inputs/PCA10.PEER10.txt"))
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
                   "Swollen", "VAS","DAS28.ESR", "DAS28.CRP")])

syn = covar_drivers(pca = t(covariates_mat), clin = df.m, alpha = 0.05)
ss = tableGrob(syn$df[syn$df$Association > -log10(0.05), ])

pdf("/home/kgoldmann/Documents/PEAC_eqtl/Results/drivers_syn_all.pdf", height=10)
syn$plot
print(grid.arrange(ss))
dev.off()

##############
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
                   "Swollen", "VAS","DAS28.ESR", "DAS28.CRP")])

bld = covar_drivers(pca = t(covariates_mat), clin = df.m, alpha = 0.05)
bs = tableGrob(bld$df[bld$df$Association > -log10(0.05), ])


pdf("/home/kgoldmann/Documents/PEAC_eqtl/Results/drivers_bld_all.pdf", height=10)
bld$plot
print(grid.arrange(bs))
dev.off()

