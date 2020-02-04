#######################
#  EQTL
#######################


# Covariate chooser, select which set of covariates provides the most signficant covariates

cvs = list.files("/home/kgoldmann/Documents/PEAC_eqtl/Results", full.names = TRUE, recursive = T)
cvs = cvs[! grepl("fantom", cvs)]
cvs = cvs[grepl("_summary.csv", cvs)]


# Output the summary data
################################
summary = data.frame()
for(cov in cvs){
  temp = cbind(read.table(cov, nrows=4, skip=1, header=F, sep=","), "cvs"=unlist(lapply(strsplit(cov, split="/"), "[[", 8)), "tissue"=unlist(lapply(strsplit(cov, split="/"), "[[", 7) ) )
  colnames(temp) = c("cutoff", "no", "cvs", "tissue")
  summary = rbind(summary, temp)
}

summary$no = as.numeric(as.character(trimws(summary$no, "both")))

write.table(summary, "/home/kgoldmann/Documents/PEAC_eqtl/Results/Covariate_significance_all.csv", sep=",", row.names = FALSE)

max_sig_combination = summary[summary$cutoff == "unique genes with p < 5e-8 and maf > 0.05", ]


# Pick covariates for synovium
################################
max_sig_combination_syn = max_sig_combination[max_sig_combination$tissue == "Synovium", ]
max_sig_combination_syn = as.character(max_sig_combination_syn$cvs[which(max_sig_combination_syn$no == max(max_sig_combination_syn$no))])

files = unlist(lapply(strsplit(cvs[grepl("Synovium", cvs)], split="/"), "[[", 8))
files = files[files != max_sig_combination_syn]
syn.remove = paste0("/home/kgoldmann/Documents/PEAC_eqtl/Results/Synovium/", files, "/Synovium_", files, ".csv")

# Important file: 
syn.keep = paste0("/home/kgoldmann/Documents/PEAC_eqtl/Results/Synovium/", max_sig_combination_syn, "/Synovium_", max_sig_combination_syn, ".csv")


file.exists(syn.remove)
file.remove(syn.remove)

# Output a file for significant snps only!
mat.df = fread(syn.keep)
mat.df = mat.df[mat.df$pvalue < 5e-8, ]
write.table(mat.df, gsub(".csv", "_significant_snps", syn.keep), sep=",", quote=F, row.names = F)


# Pick covariates for blood
################################
max_sig_combination_bld = max_sig_combination[max_sig_combination$tissue == "Blood", ]
max_sig_combination_bld = as.character(max_sig_combination_bld$cvs[which(max_sig_combination_bld$no == max(max_sig_combination_bld$no))])


files = unlist(lapply(strsplit(cvs[grepl("Blood", cvs)], split="/"), "[[", 8))
files = files[files != max_sig_combination_bld]
bld.remove = paste0("/home/kgoldmann/Documents/PEAC_eqtl/Results/Blood/", files, "/Blood_", files, ".csv")

# Important file: 
bld.keep = paste0("/home/kgoldmann/Documents/PEAC_eqtl/Results/Blood/", max_sig_combination_bld, "/Blood_", max_sig_combination_bld, ".csv")


file.exists(bld.remove)
file.remove(bld.remove)

# Output a file for significant snps only!
mat.df = fread(bld.keep)
mat.df = mat.df[mat.df$pvalue < 5e-8, ]
write.table(mat.df, gsub(".csv", "_significant_snps", bld.keep), sep=",", quote=F, row.names = F)



#######################
#  trait
#######################


# Covariate chooser, select which set of covariates provides the most signficant covariates

cvs = list.files("/home/kgoldmann/Documents/PEAC_eqtl/Results", full.names = TRUE, recursive = T)
cvs = cvs[grepl("fantom", cvs)]
cvs = cvs[grepl("_summary.csv", cvs)]


# Output the summary data
################################
summary = data.frame()
for(cov in cvs){
  temp = cbind(read.table(cov, nrows=4, skip=1, header=F, sep=","), "cvs"=unlist(lapply(strsplit(cov, split="/"), "[[", 8)), "tissue"=unlist(lapply(strsplit(cov, split="/"), "[[", 7) ) )
  colnames(temp) = c("cutoff", "no", "cvs", "tissue")
  summary = rbind(summary, temp)
}

summary$no = as.numeric(as.character(trimws(summary$no, "both")))

write.table(summary, "/home/kgoldmann/Documents/PEAC_eqtl/Results/Covariate_significance_all_trait.csv", sep=",", row.names = FALSE)

max_sig_combination = summary[summary$cutoff == "unique modules with p < 5e-8 and maf > 0.05", ]


# Pick covariates for synovium
################################
max_sig_combination_syn = max_sig_combination[max_sig_combination$tissue == "Synovium", ]
if(all(max_sig_combination_syn$no == 0)){max_sig_combination_syn = "EV1234_PEER1234"} else{
  max_sig_combination_syn = as.character(max_sig_combination_syn$cvs[which(max_sig_combination_syn$no == max(max_sig_combination_syn$no))])[1]
}

files = unlist(lapply(strsplit(cvs[grepl("Synovium", cvs)], split="/"), "[[", 8))
files = files[files != max_sig_combination_syn]
syn.remove = paste0("/home/kgoldmann/Documents/PEAC_eqtl/Results/Synovium/", files, "/Synovium_", files, "_p1_fantom5mods.csv")


# Important file: 
syn.keep = paste0("/home/kgoldmann/Documents/PEAC_eqtl/Results/Blood/", max_sig_combination_syn, "/Blood_", max_sig_combination_syn, "_p1_fantom5mods.csv")

file.exists(syn.remove)
file.remove(syn.remove)

# Output a file for significant snps only!
mat.df = fread(paste0("/home/kgoldmann/Documents/PEAC_eqtl/Results/Synovium/", max_sig_combination_syn, "/Synovium_", max_sig_combination_syn, ".csv"))
mat.df = mat.df[mat.df$pvalue < 5e-5, ]
write.table(mat.df, paste0("/home/kgoldmann/Documents/PEAC_eqtl/Results/Synovium/", max_sig_combination_syn, "/Synovium_trait_", max_sig_combination_syn, "_significant_snps.csv"), sep=",", quote=F, row.names = F)


# Pick covariates for blood
################################
max_sig_combination_bld = max_sig_combination[max_sig_combination$tissue == "Blood", ]
max_sig_combination_bld = as.character(max_sig_combination_bld$cvs[which(max_sig_combination_bld$no == max(max_sig_combination_bld$no))])[1]


files = unlist(lapply(strsplit(cvs[grepl("Blood", cvs)], split="/"), "[[", 8))
files = files[files != max_sig_combination_bld]
bld.remove = paste0("/home/kgoldmann/Documents/PEAC_eqtl/Results/Blood/", files, "/Blood_", files, "_p1_fantom5mods.csv")

# Important file: 
bld.keep = paste0("/home/kgoldmann/Documents/PEAC_eqtl/Results/Blood/", max_sig_combination_bld, "/Blood_", max_sig_combination_bld, "_p1_fantom5mods.csv")

file.exists(bld.remove)
file.remove(bld.remove)

# Output a file for significant snps only!
mat.df = fread(paste0("/home/kgoldmann/Documents/PEAC_eqtl/Results/Blood/", max_sig_combination_bld, "/Blood_", max_sig_combination_bld, "_p1_fantom5mods.csv"))
mat.df = mat.df[mat.df$pvalue < 5e-5, ]
write.table(mat.df, paste0("/home/kgoldmann/Documents/PEAC_eqtl/Results/Blood/", max_sig_combination_bld, "/Blood_trait_", max_sig_combination_bld, "_significant_snps.csv"), sep=",", quote=F, row.names = F)

