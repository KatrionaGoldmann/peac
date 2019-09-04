library(gdsfmt)
library(SNPRelate)
library(data.table)
library(peer)
library(DESeq2)
library(biomaRt)
library(cqn)
library(biomaRt)

pcs.peer <- function(gds.in, snpL.in, n=10,  ld, eaf, counts, meta.in, prefix, gene.coord, out) {

    ## get meta data file to relate ids and covariates:
    meta <- read.table(meta.in)
    
    ## get genotypes
    #    snpgdsClose(sampgds)
    sampgds <- snpgdsOpen(gds.in)
                      
    ## get eigenvectors
    snpL <- readRDS(snpL.in)  
    SL <- snpgdsPCASampLoading(snpL, sampgds)  
    ev <- data.table(sample.id = SL$sample.id)
    for(i in 1:n){ ev[ , paste0("EV",i) := SL$eigenvect[,i] ]  }
    # if blood convert to the genentech id
    # if (nrow(meta) < 50) ev$sample.id = meta$GenentechID[match(ev$sample.id, meta$Sample.final)]

    ## select samples with expression data and make sure both datasets are in the same order
    expr <- read.table(counts, header = T)
    matExp <- as.matrix(expr)
    # remove columns all NA
    matExp = matExp[, colSums(is.na(matExp)) < nrow(matExp)]
    
    #meta$Batch = "GenentechBatch1"
    colnames(meta) = c("GenentechID", "HostpitalNumber", "Ethnicity", "Gender", "Sample.final", "Batch")
    
    vcf_in <- meta[meta$GenentechID %in% colnames(matExp) | meta$Sample.final %in% colnames(matExp), c("GenentechID", "Sample.final", "Batch", "Gender")]

    ## select and order ev 
    vcf_in = vcf_in[vcf_in$Sample.final %in% ev$sample.id, ]
    vcf_in = vcf_in[vcf_in$GenentechID %in% colnames(matExp), ]
    ev  <- ev[ev$sample.id %in% vcf_in$Sample.final, ]
    vcf_in = vcf_in[match(ev$sample.id, vcf_in$Sample.final), ]
    
    ## order matExp as in vcf_in
    matExp <- matExp[, as.character(vcf_in$GenentechID)]
    
    ## for eqtl only use genes in chrom 1-22
    ## get gene coordinates and subset
    gene.c <- fread(gene.coord)
    gene.c <- gene.c[chrom %in% 1:22,]

    matExp <- matExp[rownames(matExp) %in% gene.c$gene_id, ]

    ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh = 37)
    chr_genes <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol','chromosome_name','start_position','end_position'), filters =
                            'ensembl_gene_id', values =unique(gene.c$gene_id), mart = ensembl)
    gene.c$symbol = chr_genes$hgnc_symbol[match(gene.c$gene_id, chr_genes$ensembl_gene_id)]
    
    
    ## save gene.c as for matrixqtl gene_location file
    setcolorder(gene.c, c("gene_id", "chrom", "start", "end", "symbol"))
    write.table(gene.c, row.names=F, file=out[['geneLoc']])
                
    ## get cqn normalised expression and  peer factors    
    ## need gene length and GC content from biomart library, use longest transcript

    ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")

    gc <- data.table(getBM(attributes = c("ensembl_gene_id", "transcript_length", "percentage_gene_gc_content"),
                           filters = "ensembl_gene_id",
                           values = rownames(matExp),
                           mart = ensembl))

    setkey(gc, ensembl_gene_id, transcript_length)

    ## get gc and length for the longest transcript
    gc <- gc[, .SD[.N], ensembl_gene_id]
    setnames(gc, names(gc), c("gene_id", "length", "gc"))

    ## cqn normalise
    cqn.peac <- cqn(counts = matExp, x = gc$gc/100, lengths = gc$length,   verbose=T)

    ## get normalised values
    rpkm.peac <- cqn.peac$y + cqn.peac$offset

    ## save normalised values
    write.table(rpkm.peac, file=out[['expressionCqn']])

    ## run model
    model = PEER()
    PEER_setPhenoMean(model,t(rpkm.peac))
    PEER_setNk(model,n)
    PEER_update(model)
    factors = PEER_getX(model)
    colnames(factors) <- paste0("PEER",1:10)
    rownames(factors) <- colnames(rpkm.peac)

    ## need to transpose ev and factors and then cbind and save
    ev=data.frame(ev)
    rownames(ev) = ev[, 1]
    ev <- t(ev[,colnames(ev) != "sample.id"])
    factors <- t(factors)
    
    ## I need to save output covars as in rule: 1 pcs with 1:10 peer factors, then 2 pcs with peer 1:10, etc, rows covariate name, cols samples
    indx=0
    for(i in 1:nrow(ev)){
        for(j in 1:nrow(factors)){
            indx= indx + 1 ## to count files 1:100
            tmp <- data.table(rbind(ev[1:i, ,drop=F], factors[1:j, , drop=F]), keep.rownames=T)
            write.table(tmp, row.names=T, col.names=T, file=out[['covs']][indx])
        }
    }

    ## Same with pcs and batch and sex, add pcs 1:10
    Bsex <- vcf_in[,  c("Batch", "Gender")]
    Bsex = data.table(Bsex)
    ## recode as numeric
    Bsex[Gender == "M", Gender:="0"][Gender == "F", Gender:="1"]
    b <- unique(Bsex$Batch)
    names(b) <- 0:(length(b)-1)
    for (i in seq_along(b)){
        Bsex[Batch == b[i], Batch:= names(b)[i]]
    }
    ## transpose and make numeric
    mBsex <- t(Bsex)
    mBsex <- apply(mBsex, 2, as.numeric)
    rownames(mBsex) <- names(Bsex)
    
    for( i in 1:nrow(ev)){
        tmp <- data.table(rbind(ev[1:i,, drop=F], mBsex), keep.rownames=T)
        write.table(tmp, row.names=T, col.names=T, file=out[['covfix']][i])

    }
          
    ## Genotypes: LD prune
    gprune <- snpgdsLDpruning(sampgds, maf=eaf, ld.threshold=ld)
    SNPprune <- unlist(gprune)

    ## get genotype matrix for selected snps
    gmat <- snpgdsGetGeno(sampgds, snp.id=SNPprune)

    ## add rownames (sample id) and colnames (snps) to gmat for snp.id
    ## get chrom:pos:ref:alt. snp.id in gds is a an inernal id integer
    ## from 1:N
    
    sampId <- read.gdsn(index.gdsn(sampgds, "sample.id"))
    rs <- readex.gdsn(index.gdsn(sampgds, "snp.rs.id"),sel=SNPprune)
    chrom <- readex.gdsn(index.gdsn(sampgds, "snp.chromosome"),sel=SNPprune)
    pos <- readex.gdsn(index.gdsn(sampgds, "snp.position"), sel=SNPprune)
    allele <- readex.gdsn(index.gdsn(sampgds, "snp.allele"), sel=SNPprune)
    ref <- gsub("/.*", "", allele)
    alt <- gsub(".*/", "", allele)

    rownames(gmat) <- sampId
    colnames(gmat) <- paste(chrom, pos,ref,alt,sep=":")

    
   
    
    ## SNPrelate/gds counts the number of "A" alleles, need to change 2 to 0 and 0 to 2
    g <- copy(gmat)
    g[gmat==2] <- 0
    g[gmat==0] <- 2

    ## transpose g, select samples and order as vcf_in and save
    g <- t(g)
    g <- g[, colnames(ev)]
    
    write.table(g, file=out[['geno']], col.names=T, row.names=T, quote=F)

    ## make snp location file as per matrixqtl
    snpDT <- data.table(snp=rownames(g), chr=chrom, pos=pos, rs=rs)
    snpDT$rs.id = gsub("\\:.*", "", snpDT$rs)
    write.table(snpDT, row.names=F, col.names=T, file=out[['snpLoc']], quote=F)
       
    snpgdsClose(sampgds)
    
}

###################
# Synovium
###################

counts="/media/d1/KG_Outputs/Syn_out_KG/RNA_counts/groups/Genentech.txt"
meta.in="/home/kgoldmann/Documents/PEAC_eqtl/Data/PEAC/PEAC_eth_syn.txt"
gds.in='/media/d1/KG_Outputs/Syn_out_KG/DNA/PEAC_PCA.gds'
snpL.in= "/media/d1/KG_Outputs/Syn_out_KG/DNA/RP_loads.rds"
gene.coord="/media/d1/KG_Outputs/Syn_out_KG/gene_coord.txt"
ld=1 # dont use a snp threshold
eaf=NaN # dont set a threshold 0.05
n=10
prefix=c("pcs", "peerCqn")

out = list(covs=paste("/media/d1/KG_Outputs/Syn_out_KG/matqtl/inputs/PCA", rep(1:10, 10), ".PEER", rep(1:10, each=10), ".txt", sep=""),
           geneLoc="/media/d1/KG_Outputs/Syn_out_KG/matqtl/inputs/gene_location.txt",
           expressionCqn="/media/d1/KG_Outputs/Syn_out_KG/matqtl/inputs/gene_expression_cqn.txt",
           covfix=paste("/media/d1/KG_Outputs/Syn_out_KG/matqtl/inputs/PCA", 1:10, ".covSexBatch.txt", sep=""),
           geno="/media/d1/KG_Outputs/Syn_out_KG/matqtl/inputs/genotype.txt",
           snpLoc="/media/d1/KG_Outputs/Syn_out_KG/matqtl/inputs/snp_location.txt")

pcs.peer(gds.in, snpL.in, n=n,  ld=ld, eaf=eaf, counts=counts, meta.in=meta.in, prefix="", gene.coord=gene.coord, out=out)

###################
# Blood
###################

bld.counts="/home/kgoldmann/Documents/PEAC_eqtl/Outputs_Blood/RNA_counts/groups/Genentech.txt"
bld.meta.in="/home/kgoldmann/Documents/PEAC_eqtl/Data/PEAC/PEAC_eth_blood.txt"
bld.gds.in='/home/kgoldmann/Documents/PEAC_eqtl/Outputs_Blood/DNA/PEAC_PCA.gds'
bld.snpL.in= "/home/kgoldmann/Documents/PEAC_eqtl/Outputs_Blood/DNA/RP_loads.rds"
bld.gene.coord="/home/kgoldmann/Documents/PEAC_eqtl/Outputs_Blood/gene_coord.txt"
bld.ld=1
bld.eaf=NaN
bld.n=10
bld.prefix=c("pcs", "peerCqn")

bld.out = list(covs=paste("/home/kgoldmann/Documents/PEAC_eqtl/Outputs_Blood/matqtl/inputs/PCA", rep(1:10, 10), ".PEER", rep(1:10, each=10), ".txt", sep=""),
           geneLoc="/home/kgoldmann/Documents/PEAC_eqtl/Outputs_Blood/matqtl/inputs/gene_location.txt",
           expressionCqn="/home/kgoldmann/Documents/PEAC_eqtl/Outputs_Blood/matqtl/inputs/gene_expression_cqn.txt",
           covfix=paste("/home/kgoldmann/Documents/PEAC_eqtl/Outputs_Blood/matqtl/inputs/PCA", 1:10, ".covSexBatch.txt", sep=""),
           geno="/home/kgoldmann/Documents/PEAC_eqtl/Outputs_Blood/matqtl/inputs/genotype.txt",
           snpLoc="/home/kgoldmann/Documents/PEAC_eqtl/Outputs_Blood/matqtl/inputs/snp_location.txt")


pcs.peer(gds.in = bld.gds.in, snpL.in=bld.snpL.in, n=bld.n,  ld=bld.ld, eaf=bld.eaf, counts=bld.counts, 
         meta.in=bld.meta.in, prefix="", gene.coord=bld.gene.coord, out=bld.out)



