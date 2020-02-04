library(dplyr)

plot_drivers <- function(pcs,
                         clin,
                         block = NULL,
                         unblock = NULL,
                         kernel = NULL,
                         kpar = NULL,
                         top = NULL,
                         n.pc = 5L,
                         label = FALSE,
                         alpha = 0.05,
                         p.adj = NULL,
                         title = 'Variation By Feature',
                         legend = 'right',
                         hover = FALSE) {
  
  
  sig <- function(j, pc) {                       # p-val fn
    mod <- lm(pcs[, pc] ~ clin[[j]])
    if_else(clin[[j]] %>% is.numeric, summary(mod)$coef[2L, 4L], anova(mod)[1L, 5L])
  }
  
  df <- expand.grid(Feature = colnames(clin), PC = colnames(pcs)) %>%
    rowwise(.) %>%
    dplyr::mutate(Association = sig(Feature, PC)) %>%   # Populate
    ungroup(.)
  if (!p.adj %>% is.null) {
    df <- df %>% mutate(Association = p.adjust(Association, method = p.adj))
  }
  df <- df %>% 
    mutate(Significant = if_else(Association <= alpha, TRUE, FALSE), Association = -log(Association))
  
  # Build plot
  if (!p.adj %>% is.null && p.adj %in% c('fdr', 'BH', 'BY')) {
    leg_lab <- expression(~-log(italic(q)))
  } else {
    leg_lab <- expression(~-log(italic(p)))
  }
  p <- ggplot(df, aes(PC, Feature, fill = Association, text = Association,
                      color = Significant)) +
    geom_tile(size = 1L, width = 0.9, height = 0.9) +
    coord_equal() +
    scale_fill_gradientn(colors = c('white',  'dodgerblue1', 'dodgerblue3', 'dodgerblue4'), name = leg_lab) +
    scale_color_manual(values = c('grey90', 'black')) +
    guides(color = FALSE) +
    labs(title = title, x = '', y='') +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 315, hjust = 0))
  if (label) {
    p <- p + geom_text(aes(label = round(Association, 2L)))
  }
  return(p)
}

#Synovium
######################
pcs = read.table("/media/d1/KG_Outputs/Syn_out_KG/matqtl/inputs/PCA10.PEER10.txt")
rownames(pcs) = pcs$rn
pcs = pcs[, colnames(pcs) != "rn"]
pcs = t(pcs)


clin = read.table("/home/kgoldmann/Documents/PEAC_eqtl/Data/PEAC/PEAC_eth_syn.txt")
clin = clin[match(rownames(pcs), clin$vcf_id), ]
identical(as.character(clin$vcf_id), rownames(pcs))

clin = clin[, which(apply(clin, 2, function(x) length(unique(x))) < nrow(clin))]
clin = clin[, colSums(is.na(clin)) < nrow(clin)]

syn_drivers = plot_drivers(pcs, clin, label=F, title="Synovium")
syn_drivers4 = plot_drivers(pcs[, c(paste0("EV", 1:4), paste0("PEER", 1:4))], clin, label=F, title="Synovium")



# Blood
#######################
pcs = read.table("/media/d1/KG_Outputs/Bld_out_KG/matqtl/inputs/PCA10.PEER10.txt")
rownames(pcs) = pcs$rn
pcs = pcs[, colnames(pcs) != "rn"]
pcs = t(pcs)


clin = read.table("/home/kgoldmann/Documents/PEAC_eqtl/Data/PEAC/PEAC_eth_bld.txt")
clin = clin[match(rownames(pcs), clin$vcf_id), ]
identical(as.character(clin$vcf_id), rownames(pcs))

clin = clin[, which(apply(clin, 2, function(x) length(unique(x))) < nrow(clin))]
clin = clin[, which(apply(clin, 2, function(x) length(unique(x))) > 1)]
clin = clin[, colSums(is.na(clin)) < nrow(clin)]

bld_drivers = plot_drivers(pcs, clin, label=F, title="Blood")
bld_drivers4 = plot_drivers(pcs[, c(paste0("EV", 1:4), paste0("PEER", 1:4))], clin, label=F, title="Blood")


# Output 
#########################
pdf("/home/kgoldmann/Documents/PEAC_eqtl/Results/Covariate_drivers_all.pdf", width=12)
ggpubr::ggarrange(syn_drivers, bld_drivers, ncol=2, nrow=1, common.legend = T)
ggpubr::ggarrange(syn_drivers4, bld_drivers4, ncol=2, nrow=1, common.legend = T)
dev.off()
