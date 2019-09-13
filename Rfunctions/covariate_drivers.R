covar_drivers<-function(pca, clin, index='Sample', block=NULL, unblock=NULL, kernel=NULL, 
                        kpar=NULL, top=NULL, n.pc=5L, label=FALSE, 
                        alpha=NULL,
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

covar_drivers_specific<-function(pca, clin, index='Sample', block=NULL, unblock=NULL, kernel=NULL, kpar=NULL, top=NULL, n.pc=5L, label=FALSE, alpha=NULL,
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
