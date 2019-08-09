library(MatrixEQTL)

#' cis-eqtl analysis using Matrixeqtl
#' http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/runit.html
#' calling Matrix_eQTL main function

# Select model: Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR

# Genotype file name
SNP_file_name = "/home/kgoldmann/Documents/PEAC_eqtl/Outputs/matqtl/inputs/genotype.txt"
snps_location_file_name =  "/home/kgoldmann/Documents/PEAC_eqtl/Outputs/matqtl/inputs/snp_location.txt"

# Gene expression file name
expression_file_name = "/home/kgoldmann/Documents/PEAC_eqtl/Outputs/matqtl/inputs/gene_expression_cqn.txt"
gene_location_file_name = "/home/kgoldmann/Documents/PEAC_eqtl/Outputs/matqtl/inputs/gene_location.txt"

# Covariates file name
# Set to character() for no covariates
covariates_file_name = character() #"/home/kgoldmann/Documents/PEAC_eqtl/Outputs/matqtl/inputs/{pcs}.{peerBSex}.txt"
covariates_mat =  read.table(paste0("/home/kgoldmann/Documents/PEAC_eqtl/Outputs/matqtl/inputs/PCA4.PEER4.txt"))
rownames(covariates_mat) = covariates_mat[, 1]
covariates_mat = covariates_mat[paste0("PEER", 1:4), 2:ncol(covariates_mat)]
covariates_mat = as.matrix(covariates_mat)

# Output file name
output_file_name_cis = "/home/kgoldmann/Documents/PEAC_eqtl/Outputs/matqtl/Output.txt"


# Only associations significant at this level will be saved
pvOutputThreshold_cis = 0.02
#pvOutputThreshold_tra = 0; ## only cis-eqtl

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();

# Distance for local gene-SNP pairs
cisDist = 5e5;

## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = " ";      # the ' ' character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = " ";      # the  character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

## Load covariates
cvrt = SlicedData$new();
cvrt$fileDelimiter = " ";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 0;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
    #cvrt$CreateFromMatrix(mat)
    cvrt$LoadFile(covariates_file_name);
} else{cvrt$CreateFromMatrix(covariates_mat)}

## Run the analysis
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);


me = Matrix_eQTL_main(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name     = NULL, ## for trans-eqtl
    pvOutputThreshold     = 0, ## only cis-eqtl
    useModel = useModel,
    errorCovariance = errorCovariance,
    verbose = TRUE,
    output_file_name.cis = NULL,
    pvOutputThreshold.cis = pvOutputThreshold_cis,
    snpspos = snpspos,
    genepos = genepos,
    cisDist = cisDist,
    pvalue.hist = "qqplot",
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE);

out <- me$cis$eqtls

write.table(out, row.names=F, file=output_file_name_cis)

#unlink(output_file_name_cis);
