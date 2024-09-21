
# create new gene list for FIT

d <- read.table("230728_500_genes_id_conversion.tsv",
                header=T, sep="\t")

genes <- d$rnaseq
genes <- genes[!is.na(genes)]
genes <- genes[-475]

genes <- as.character(genes)

save(genes, file="genes")
