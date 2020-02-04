

#' cis-eqtl analysis using Matrixeqtl
#' http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/runit.html
#' calling Matrix_eQTL main function

dir.create("/home/kgoldmann/Documents/PEAC_eqtl/Results", showWarnings = F)

source("/home/kgoldmann/Documents/peac/Rfunctions/matrixEQTL.R")

###############
# Synovium
###############
# Covariates file name
args = commandArgs()

i = args[6]
tiss = args[7]
if(tiss == "Synovium"){
  covariates_mat =  read.table(paste0("/media/d1/KG_Outputs/Syn_out_KG/matqtl/inputs/PCA4.PEER4.txt"))
  rownames(covariates_mat) = covariates_mat[, "rn"]
  covariates_mat = as.matrix(covariates_mat[, colnames(covariates_mat) != "rn"])

  meta = read.table("/home/kgoldmann/Documents/PEAC_eqtl/Data/PEAC/PEAC_eth.txt")
  covariates_mat2 = t(apply(covariates_mat, 1, as.numeric))
  dimnames(covariates_mat2) = dimnames(covariates_mat)
  covariates_mat = rbind(covariates_mat2, "none"=1)

  dir = "/media/d1/KG_Outputs/Syn_out_KG/"
  opts = list(c(paste0("EV", 1:4), "PEER1"), c(paste0("EV", 1:4), paste0("PEER", c(1, 3, 4))), c(paste0("EV", 1:4), paste0("PEER", 1:4)), c("none"), paste0("EV", 1:4))
} else{
  covariates_mat =  read.table(paste0("/media/d1/KG_Outputs/Bld_out_KG/matqtl/inputs/PCA4.PEER4.txt"))
  rownames(covariates_mat) = covariates_mat[, "rn"]
  covariates_mat = as.matrix(covariates_mat[, colnames(covariates_mat) != "rn"])

  meta = read.table("/home/kgoldmann/Documents/PEAC_eqtl/Data/PEAC/PEAC_eth_blood.txt")
  meta = meta[meta$vcf_id %in% colnames(covariates_mat), ]
  covariates_mat2 = t(apply(covariates_mat, 1, as.numeric))
  dimnames(covariates_mat2) = dimnames(covariates_mat)
  covariates_mat = rbind(covariates_mat2, "none"=1)

  dir = "/media/d1/KG_Outputs/Bld_out_KG/"
  opts = list(c(paste0("EV", 1:4), paste0("PEER", 1:3)), c(paste0("EV", 1:4), paste0("PEER", 1:4)), c("none"), paste0("EV", 1:4))
}


# runner = function(i) {
#   keep(opts, covariates_mat, output.creater, runner, sure=T, all=T)
#   gc()
cv = as.character(opts[[as.numeric(i)]])
print("###")
print(c("###", cv))
cm = covariates_mat
if (length(cv) == 1) cm = matrix(ncol = ncol(covariates_mat), nrow=0)
output.creater(cv.df=cm, results.dir=dir, tissue=tiss, cv.names=cv, rs.plotter=c())
#   keep(opts, covariates_mat, output.creater, runner, sure=T, all=T)
#   
