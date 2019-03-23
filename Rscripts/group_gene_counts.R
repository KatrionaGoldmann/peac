library(data.table)

## get args from snakemake

in.files <- unlist(snakemake@input)

out.file <- snakemake@output[[1]]

lib.size.file <-  snakemake@output[[2]]

#in.files = paste("/home/kgoldmann/Documents/PEAC_eqtl/Outputs/RNA_counts/", output.df$SampleID..QMUL.or.Genentech., ".txt", sep="")

##filter <- as.numeric(snakemake@params[['filter']])

filter =0

## process inputs
lcounts <- lapply(in.files, fread)

df <- data.frame(rep(NA, 57905))
for (i in lcounts){
  df = cbind(df, i)
}
df = df[, 2:ncol(df)]
rownames(df) = df$V1
df = df[, c("V1", colnames(df)[grepl("SAM", colnames(df))])]
counts=data.table(df)

lib_size= log(colSums(as.matrix(counts[,2:ncol(counts),with=F])))

lib_size <- matrix(lib_size, ncol=1, dimnames=list(names(counts)[2:ncol(counts)], "lib.s"))

## filter counts, remove genes with 0 counts in all samples
counts <- counts[which(rowSums(counts[, 2:ncol(counts), with=F])>filter),]

## save files

write.table(counts, out.file, row.names=F)

saveRDS(lib_size, lib.size.file)

