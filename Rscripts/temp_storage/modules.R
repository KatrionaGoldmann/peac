# Look at the modules
f5_genes = read.csv2("/home/kgoldmann/NAS/GCPEAC/Sharmila/Fantom_5/fantom5 modules v4.csv", sep=",")
B_Cells = as.character(f5_genes$gene[f5_genes$celltype == "CD19+ B Cells"])
B_Cells

B_Cells %in% output$gene.name
bcs <- getBM(attributes=c('ensembl_gene_id', 'ensembl_transcript_id','hgnc_symbol', "chromosome_name", "start_position", "end_position"), values =B_Cells, mart = ensembl)
bcs <- bcs[match(B_Cells, bcs$hgnc_symbol), ]

# Get the list of chromosomes these genes belong to
