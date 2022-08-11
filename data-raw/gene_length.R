BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
library("TxDb.Mmusculus.UCSC.mm10.knownGene")
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
ls(txdb)
exon_txdb <- exons(txdb)
head(exon_txdb)
genes_txdb <- genes(txdb)
head(genes_txdb)
o  <-  findOverlaps(exon_txdb,genes_txdb)
o
t1 <- exon_txdb[queryHits(o)]
t2 <- genes_txdb[subjectHits(o)]
t1 <- as.data.frame(t1)
t1$geneid <- mcols(t2)[,1]
g_l <- lapply(split(t1,t1$geneid),function(x){
  head(x)
  tmp <- apply(x,1,function(y){
    y[2]:y[3]
  })
  length(unique(unlist(tmp)))
})
head(g_l)
g_l <- data.frame(gene_id=names(g_l), length=as.numeric(g_l))
library(org.Mm.eg.db)
s2g <- toTable(org.Mm.egSYMBOL)
s3g <- toTable(org.Mm.egENSEMBL)
head(s2g)
head(s3g)
g_l <- merge(g_l, s2g, by = 'gene_id')
g_l <- merge(g_l, s3g, by = 'gene_id')
head(g_l)
mm_gene_length <- g_l
usethis::use_data(mm_gene_length, overwrite = TRUE)




########################################### Human
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
ls('package:TxDb.Hsapiens.UCSC.hg38.knownGene')
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
ls(txdb)
exon_txdb <- exons(txdb)
head(exon_txdb)
genes_txdb <- genes(txdb)
head(genes_txdb)
o  <-  findOverlaps(exon_txdb,genes_txdb)
o
t1 <- exon_txdb[queryHits(o)]
t2 <- genes_txdb[subjectHits(o)]
t1 <- as.data.frame(t1)
t1$geneid <- mcols(t2)[,1]
head(t1)
g_l <- lapply(split(t1,t1$geneid),function(x){
  # x=split(t1,t1$geneid)[[1]]
  head(x)
  tmp <- apply(x,1,function(y){
    y[2]:y[3]
  })
  length(unique(unlist(tmp)))
})
head(g_l)
g_l <- data.frame(gene_id=names(g_l),length=as.numeric(g_l))
library(org.Hs.eg.db)
s2g <- toTable(org.Hs.egSYMBOL)
s3g <- toTable(org.Hs.egENSEMBL)
head(s2g)
head(s3g)
g_l <- merge(g_l, s2g, by = 'gene_id')
g_l <- merge(g_l, s3g, by = 'gene_id')
head(g_l)
hs_gene_length <- g_l
usethis::use_data(hs_gene_length, overwrite = TRUE)

