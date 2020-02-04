##############################################
##############################################
# First lets look at the count data
##############################################
##############################################


# Create a file for the sample names and path to the corresponding fastq file (on NAS07)
samples = list.files("/home/kgoldmann/Documents/rnaseq/PEAC_RawData", full.names = T, recursive = T, pattern = "fastq.gz")

samples.df = data.frame("samples"=as.character(samples), "path"=as.character(samples))
samples.df$samples = lapply(samples.df$samples, function(x) unlist(lapply(  strsplit(as.character(x), split="/"), "[[", 7)) )
samples.df$samples = gsub("\\_.*", "", samples.df$samples)

samples.df = samples.df[! duplicated(samples.df$samples), ] # 404 unique samples



load("/home/kgoldmann/Documents/gcpeac/PEAC/PEAC_main_160517.rdata")
samples.df$QMUL.ID = gsub("b", "", metadata$SampleID..QMUL.ID.only.[match(samples.df$samples, metadata$SampleID..QMUL.or.Genentech.)])

# There are some BHAM* and QMUL* SAM IDs which do not match to a QMULID 
# - I am going to carry the BHAM and QMUL over to QMULID column iff 
# - there is only one sample registered to that ID to avoid confusion between tissues
# - Then I will use the new QMUL to get the SAMID
for(x in samples.df$samples[is.na(samples.df$QMUL.ID)]){
  all.sampls = metadata$SampleID..QMUL.or.Genentech.[gsub("b", "", metadata$SampleID..QMUL.ID.only.) %in% x]
  if(length(all.sampls) ==1){ # Then there is only one sample, no risk of getting syn/bld mixed up
    samples.df$QMUL.ID[samples.df$samples == x] = x
    samples.df$samples[samples.df$samples == x] = paste(all.sampls, "fixed?")
  } else{ # non-unique qmul id for that sam id, lets remove these from analysis to avoi confusion.
    samples.df$QMUL.ID[samples.df$samples == x] = "remove"
  }
}

samples.df = samples.df[! grepl("remove", samples.df$QMUL.ID), ] # 358 unique samples
samples.df$samples = gsub(" fixed\\?", "", samples.df$samples)


#########
# Align to timepoint
#########
# Match the timepoint accross
samples.df$Timepoint_lookup = lookuptable$timeline[match(samples.df$samples, lookuptable$Sample_ID)]
samples.df$Timepoint_metadata = metadata$Timepoint[match(samples.df$samples, metadata$SampleID..QMUL.or.Genentech.)]
samples.df$Timepoint_PEAC2viii = PEAC2viii$Timepoint[match(samples.df$QMUL.ID, PEAC2viii$Sample.ID)]

load("/home/kgoldmann/Documents/gcpeac/Chris/PEAC2_data/PEAC2_alldata.RData")
samples.df$Timepoint_Chris2 = meta$Timepoint[match(samples.df$samples, meta$Sample.ID)]

load("/home/kgoldmann/Documents/gcpeac/Chris/PEAC_data/PEACmetadataV12.RData")
samples.df$Timepoint_Chris = metadata$both$Timepoint[match(samples.df$samples, metadata$both$GenentechID)]

# Fix the time names
for(i in colnames(samples.df)[grepl("Timepoint", colnames(samples.df))]){
  samples.df[, i] = gsub("0", "Baseline", samples.df[, i])
  samples.df[, i] = gsub("6|6mth|9mth \\(3mths late\\)", "Six-month", samples.df[, i])
}

# run through each row to check no clashes in Timepoint
samples.df$TP = NA
for(i in 1:nrow(samples.df)){
  vals = samples.df[i, colnames(samples.df)[grepl("Timepoint", colnames(samples.df))]]
  vals = vals[! is.na(vals)]
  if(length(vals) > 0){
    if(all(vals == vals[1]) == FALSE){print(samples.df$samples[i] )} # print error if differences
    samples.df$TP[i] = vals[1]
  } 
}

samples.df = samples.df[samples.df$TP == "Baseline", ] # 232 unique samples

#########
# Align to timepoint
#########
# Match the timepoint accross
load("/home/kgoldmann/Documents/gcpeac/PEAC/PEAC_main_160517.rdata")
samples.df$Tissue_metadata = metadata$Tissue[match(samples.df$samples, metadata$SampleID..QMUL.or.Genentech.)]

load("/home/kgoldmann/Documents/gcpeac/Chris/PEAC_data/PEACmetadataV12.RData")
samples.df$Tissue_Chris = metadata$both$Tissue[match(samples.df$samples, metadata$both$GenentechID)]

# Check no mismatches
samples.df$samples[samples.df$Tissue_metadata != samples.df$Tissue_Chris & 
                     ! is.na(samples.df$Tissue_metadata) & ! is.na(samples.df$Tissue_Chris)]

samples.df$Tissue = as.character(samples.df$Tissue_metadata)
samples.df$Tissue[is.na(samples.df$Tissue)] = samples.df$Tissue_Chris[is.na(samples.df$Tissue)]

samples.df = samples.df[! is.na(samples.df$Tissue), ]

# 153 synovium samples
# 65 blood samples


##############################################
##############################################
# Now lets look at the genotype data
##############################################
##############################################

library(vcfR)


vcf <- read.vcfR("~/Documents/PEAC_eqtl/Outputs/DNA/PEAC_chr22_4PCA_all.vcf")
gt = extract.gt(vcf, element = 'DP', as.numeric = TRUE)


genotype.samps = colnames(gt)

samples.df$match = NA
samples.df$match[samples.df$QMUL.ID %in% genotype.samps] = "*"


matched.df = samples.df[! is.na(samples.df$match), ]
#matched.df = matched.df[! duplicated(matched.df$QMUL.ID), ]


table(matched.df$Tissue)
