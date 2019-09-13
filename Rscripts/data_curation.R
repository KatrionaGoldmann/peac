######################
# For Synovium
######################


# List the PEAC vcf files
vcf_list = list.files("/home/kgoldmann/tmp", full.names=T, pattern=".vcf.gz")
vcf_list = vcf_list[! grepl("tbi", vcf_list)]
write.table(vcf_list, "/home/kgoldmann/Documents/PEAC_eqtl/Data/vcf_list.txt", row.names = F, col.names = FALSE, quote=FALSE)

# Create a file for the sample names and path to the corresponding fastq file
samples = list.files("/home/kgoldmann/Documents/rnaseq/PEAC_RawData", full.names = T, recursive = T, pattern = "fastq.gz")

samples.df = data.frame("samples"=as.character(samples), "path"=as.character(samples))
samples.df$samples = lapply(samples.df$samples, function(x) unlist(lapply(  strsplit(as.character(x), split="/"), "[[", 7)) )
samples.df$samples = gsub("\\_.*", "", samples.df$samples)

load("/home/kgoldmann/Documents/gcpeac/PEAC/PEAC_main_160517.rdata")
meta.use.peac = metadata[metadata$SampleID..QMUL.or.Genentech. %in% samples.df$samples | metadata$SampleID..QMUL.ID.only. %in% samples.df$samples, ]
meta.use.peac = meta.use.peac[meta.use.peac$Tissue == "Synovium", ]
meta.use.peac = meta.use.peac[meta.use.peac$Timepoint == "Baseline", ]
meta.use.peac = meta.use.peac[meta.use.peac$Reads == "Paired", ]
meta.use.peac = meta.use.peac[! is.na(meta.use.peac$Timepoint), ]
meta.use.peac$R1 = NA
meta.use.peac$R2 = NA

for(i in 1:nrow(meta.use.peac)){
  sam = meta.use.peac$SampleID..QMUL.or.Genentech.[i]
  if(any(grepl("R1|_1", samples.df[samples.df$samples == sam, ]$path))) {
  r1 = as.character(samples.df[samples.df$samples == sam, ]$path[grepl("R1|_1", samples.df[samples.df$samples == sam, ]$path)])
  r2 = as.character(samples.df[samples.df$samples == sam, ]$path[grepl("R2|_2", samples.df[samples.df$samples == sam, ]$path)])
  meta.use.peac$R1[meta.use.peac$SampleID..QMUL.or.Genentech. == sam] = r1
  meta.use.peac$R2[meta.use.peac$SampleID..QMUL.or.Genentech. == sam] = r2
  } else{
    r1 = as.character(samples.df[samples.df$samples == sam, ]$path[grepl("R1|_1", samples.df[samples.df$samples == sam, ]$path)])
    meta.use.peac$R1[meta.use.peac$SampleID..QMUL.or.Genentech. == sam] = r1
  }
}

meta.use.peac$Diagnosis = PEAC2viii$DIAGNOSIS[match(meta.use.peac$SampleID..QMUL.ID.only., PEAC2viii$Sample.ID)]
meta.use.peac$Dir = unlist(lapply(meta.use.peac$R1, function(x) gsub("/[^/]+$", "", x)))

# Should be <SampleID..QMUL.or.Genentech.> <unique_individual_id=HospitalNumber><ID in vcf=SampleID..QMUL.ID.only.><Batch><Reads><Timepoint><Tissue><Diagnosis><GenotypeDir><fastq1.R1_path> <fastq1.R2_path> <fastq2.R1_path> <fastq2.R2_path>. This function produces a dictionary of sample_id keys and values(individualID, Batch, Reads, Timepoint, Tissue, Diagnosis, Dir with genotype, and a sub-list of (fq1, fq2, ...,)

meta.use.peac = meta.use.peac[, c("SampleID..QMUL.or.Genentech.", "HospitalNumber", "SampleID..QMUL.ID.only.", "Batch", "Reads", "Timepoint", "Tissue", "Diagnosis", "Dir", "R1", "R2" )]
write.table(meta.use.peac, "/home/kgoldmann/Documents/PEAC_eqtl/Data/RA.csv", col.names = F, row.names = F, quote=F, sep=",")




# Create a file for the reference panel vcf files
vcf  = list.files("/home/kgoldmann/Documents/PEAC_eqtl/Data/Reference/bcf", full.names = T)
vcf = vcf[! grepl("tbi|chrX|wgs", vcf)]
write.table(vcf, "/home/kgoldmann/Documents/PEAC_eqtl/Data/Reference/refPanel_list.txt", quote=F, row.names=F, col.names=F)





# Create a file for the reference panel annotations
ann = list.files("/home/kgoldmann/Documents/PEAC_eqtl/Data/Reference/annotated_legends", full.names = T)
ann = ann[! grepl("chrX|wgs", ann)]
write.table(ann, "/home/kgoldmann/Documents/PEAC_eqtl/Data/Reference/refPanel_leg_list.txt", quote=F, row.names = F, col.names = F)


# Create the Genentech file
temp = list.files("/media/d1/KG_Outputs/Syn_out_KG/RNA_counts/", pattern=".txt")
samps = gsub(".txt", "", temp)

#/home/kgoldmann/Documents/PEAC_eqtl/Outputs/RNA_counts/QMUL2008043.txt
df= read.table(paste0("/media/d1/KG_Outputs/Syn_out_KG/RNA_counts/", samps[1], ".txt"))
for(i in samps[2:length(samps)]){
  tt = read.table(paste0("/media/d1/KG_Outputs/Syn_out_KG/RNA_counts/", i, ".txt"))

  df = cbind(df, tt[match(rownames(df), rownames(tt)), ])
  colnames(df)[ncol(df)] = i
}
dir.create("/media/d1/KG_Outputs/Syn_out_KG/RNA_counts/groups/")
write.table(df, "/media/d1/KG_Outputs/Syn_out_KG/RNA_counts/groups/Genentech.txt", row.names=T, col.names=T, quote = F)

# Create the sample_meta file
df = meta.use.peac[, 1:4]
load("/home/kgoldmann/Documents/gcpeac/Chris/PEAC2_data/PEAC2_alldata.RData")
m2 = meta
df$Gender = m2$GENDER[match(df$SampleID..QMUL.ID.only., m2$Sample.ID)]

load("/home/kgoldmann/Documents/gcpeac/Chris/PEAC_data/PEACmetadataV12.RData")
m1 = metadata$baseline
df$Gender[is.na(df$Gender)] = m1$GENDER[match(df$SampleID..QMUL.ID.only.[is.na(df$Gender)], m1$SampleID)]
df$Ethnicity = m1$Ethnicity[match(df$SampleID..QMUL.ID.only., m1$SampleID)]

df = df[, c("SampleID..QMUL.or.Genentech.", "HospitalNumber" , "Ethnicity", "Gender", "SampleID..QMUL.ID.only.", "Batch" )]
colnames(df) = c("samp.et", "HostpitalNumber", "Ethnicity", "Gender", "vcf_id", "Batch")
write.table(df, "/home/kgoldmann/Documents/PEAC_eqtl/Data/PEAC/PEAC_eth_syn.txt", row.names=T, col.names=T)




# 94 for GenentechID Peac_data
# 104 for SampleID Peac_data
# 137 for sampleID in PEAC2


######################
# For Blood
######################

# List the PEAC vcf files
vcf_list = list.files("/home/kgoldmann/tmp", full.names=T, pattern=".vcf.gz")
vcf_list = vcf_list[! grepl("tbi", vcf_list)]
write.table(vcf_list, "/home/kgoldmann/Documents/PEAC_eqtl/Data/vcf_list.txt", row.names = F, col.names = FALSE, quote=FALSE)

# Create a file for the sample names and path to the corresponding fastq file
samples = list.files("/home/kgoldmann/Documents/rnaseq/PEAC_RawData", full.names = T, recursive = T, pattern = "fastq.gz")

samples.df = data.frame("samples"=as.character(samples), "path"=as.character(samples))
samples.df$samples = lapply(samples.df$samples, function(x) unlist(lapply(  strsplit(as.character(x), split="/"), "[[", 7)) )
samples.df$samples = gsub("\\_.*", "", samples.df$samples)

load("/home/kgoldmann/Documents/gcpeac/PEAC/PEAC_main_160517.rdata")
meta.use.peac = metadata[metadata$SampleID..QMUL.or.Genentech. %in% samples.df$samples | metadata$SampleID..QMUL.ID.only. %in% samples.df$samples, ]
meta.use.peac = meta.use.peac[meta.use.peac$Tissue == "Blood", ]
meta.use.peac = meta.use.peac[meta.use.peac$Timepoint == "Baseline", ]
#meta.use.peac = meta.use.peac[meta.use.peac$Reads == "Paired", ]
meta.use.peac = meta.use.peac[! is.na(meta.use.peac$Timepoint), ]
meta.use.peac$R1 = NA
meta.use.peac$R2 = NA

for(i in 1:nrow(meta.use.peac)){
  sam = meta.use.peac$SampleID..QMUL.or.Genentech.[i]
  r1 = as.character(samples.df[samples.df$samples == sam, ]$path[grepl("R1|_1", samples.df[samples.df$samples == sam, ]$path)])
  meta.use.peac$R1[meta.use.peac$SampleID..QMUL.or.Genentech. == sam] = r1
  if(any(grepl("R2|_2", samples.df[samples.df$samples == sam, ]$path))) {
    r2 = as.character(samples.df[samples.df$samples == sam, ]$path[grepl("R2|_2", samples.df[samples.df$samples == sam, ]$path)])
    meta.use.peac$R2[meta.use.peac$SampleID..QMUL.or.Genentech. == sam] = r2
  } 
}

meta.use.peac$Diagnosis = PEAC2viii$DIAGNOSIS[match(meta.use.peac$SampleID..QMUL.ID.only., PEAC2viii$Sample.ID)]
meta.use.peac$Dir = unlist(lapply(meta.use.peac$R1, function(x) gsub("/[^/]+$", "", x)))

# Should be <SampleID..QMUL.or.Genentech.> <unique_individual_id=HospitalNumber><ID in vcf=SampleID..QMUL.ID.only.><Batch><Reads><Timepoint><Tissue><Diagnosis><GenotypeDir><fastq1.R1_path> <fastq1.R2_path> <fastq2.R1_path> <fastq2.R2_path>. This function produces a dictionary of sample_id keys and values(individualID, Batch, Reads, Timepoint, Tissue, Diagnosis, Dir with genotype, and a sub-list of (fq1, fq2, ...,)

meta.use.peac = meta.use.peac[, c("SampleID..QMUL.or.Genentech.", "HospitalNumber", "SampleID..QMUL.ID.only.", "Batch", "Reads", "Timepoint", "Tissue", "Diagnosis", "Dir", "R1", "R2" )]
meta.use.peac$SampleID..QMUL.ID.only. = gsub("b", "", meta.use.peac$SampleID..QMUL.ID.only.)

write.table(meta.use.peac, "/home/kgoldmann/Documents/PEAC_eqtl/Data/RA_blood.csv", col.names = F, row.names = F, quote=F, sep=",")

# Create the Genentech file
temp = list.files("/home/kgoldmann/Documents/PEAC_eqtl/Outputs_Blood/RNA_counts/", pattern=".txt")
samps = gsub(".txt", "", temp)

#/home/kgoldmann/Documents/PEAC_eqtl/Outputs_Blood_Blood/RNA_counts/QMUL2008043.txt
df= read.table(paste0("/home/kgoldmann/Documents/PEAC_eqtl/Outputs_Blood/RNA_counts/", samps[1], ".txt"))
for(i in samps[2:length(samps)]){
  tt = read.table(paste0("/home/kgoldmann/Documents/PEAC_eqtl/Outputs_Blood/RNA_counts/", i, ".txt"))
  # if(length(rownames(tt)[! rownames(tt) %in% rownames(df)]) > 0){
  #   df = rbind(df, NA)
  #   rownames(df)[-length(rownames(tt)[! rownames(tt) %in% rownames(df)])] = rownames(tt)[! rownames(tt) %in% rownames(df)]
  # }
  df = cbind(df, tt[match(rownames(df), rownames(tt)), ])
  colnames(df)[ncol(df)] = i
}
dir.create("/home/kgoldmann/Documents/PEAC_eqtl/Outputs_Blood/RNA_counts/groups/", showWarnings = F)
write.table(df, "/home/kgoldmann/Documents/PEAC_eqtl/Outputs_Blood/RNA_counts/groups/Genentech.txt")

# Create the sample_meta file
df = meta.use.peac[, 1:4]
load("/home/kgoldmann/Documents/gcpeac/Chris/PEAC2_data/PEAC2_alldata.RData")
m2 = meta
df$Gender = m2$GENDER[match(df$SampleID..QMUL.ID.only., m2$Sample.ID)]

load("/home/kgoldmann/Documents/gcpeac/Chris/PEAC_data/PEACmetadataV12.RData")
m1 = metadata$baseline
df$Gender[is.na(df$Gender)] = m1$GENDER[match(df$SampleID..QMUL.ID.only.[is.na(df$Gender)], m1$SampleID)]
df$Ethnicity = m1$Ethnicity[match(df$SampleID..QMUL.ID.only., m1$SampleID)]


df = df[, c("SampleID..QMUL.or.Genentech.", "HospitalNumber" , "Ethnicity", "Gender", "SampleID..QMUL.ID.only.", "Batch" )]
colnames(df) = c("samp.et", "HostpitalNumber", "Ethnicity", "Gender", "vcf_id", "Batch")
write.table(df, "/home/kgoldmann/Documents/PEAC_eqtl/Data/PEAC/PEAC_eth_blood.txt", row.names=T, col.names=T)

# 94 for GenentechID Peac_data
# 104 for SampleID Peac_data
# 137 for sampleID in PEAC2
