############
library(SummarizedExperiment)
library(jaffelab)

## make directory
dir.create("count_data")

## load counts
load("preprocessed_data/rpkmCounts_Clonetech_CST_HiSeq_n6.rda")
load("preprocessed_data/rawCounts_Clonetech_CST_HiSeq_n6.rda")

# read in pheno
pd = read.csv("clonetech_cst_ip_vs_input_pheno.csv")
rownames(pd) = pd$SampleID

# merge with metrics
pd = cbind(pd[rownames(metrics),], metrics[,-1])

## gene
gr_genes <- GRanges(seqnames = geneMap$Chr,
    IRanges(geneMap$Start, geneMap$End), strand = geneMap$Strand)
names(gr_genes) <- rownames(geneMap)
mcols(gr_genes) <- DataFrame(geneMap[, - which(colnames(geneMap) %in%
    c('Chr', 'Start', 'End', 'Strand'))])

rse_gene <- SummarizedExperiment(assays = list('counts' = geneCounts),
    rowRanges = gr_genes, colData = pd)

## exon
gr_exons <- GRanges(seqnames = exonMap$Chr,
    IRanges(exonMap$Start, exonMap$End), strand = exonMap$Strand)
names(gr_exons) <- rownames(exonMap)
mcols(gr_exons) <- DataFrame(exonMap[, - which(colnames(exonMap) %in%
    c('Chr', 'Start', 'End', 'Strand'))])

rse_exon <- SummarizedExperiment(assays = list('counts' = exonCounts),
    rowRanges = gr_exons, colData = pd)

## jxn
jIndex = which(rowSums(as.data.frame(jCounts) > 0) > 2)
rse_jx <- SummarizedExperiment(
	assays = list('counts' = jCounts[jIndex,]),
    rowRanges = jMap[jIndex,], colData = pd)

## tx
txMap_coord = geneMap[txMap$gencodeID,]
gr_txs <- GRanges(seqnames = txMap_coord$Chr,
    IRanges(txMap_coord$Start, txMap_coord$End), 
	strand = txMap_coord$Strand)
names(gr_txs) <- rownames(txMap)
mcols(gr_txs) <- DataFrame(txMap)

rse_tx <- SummarizedExperiment(assays = list('tpm' = txTpm),
    rowRanges = gr_txs, colData = pd)

	
## save
save(rse_gene, file="count_data/rse_gene_Clonetech_CST_IpVsInput_n6_annotated.Rdata")
save(rse_exon, file="count_data/rse_exon_Clonetech_CST_IpVsInput_n6_annotated.Rdata")
save(rse_jx, file="count_data/rse_jx_Clonetech_CST_IpVsInput_n6_annotated.Rdata")
save(rse_tx, file="count_data/rse_tx_Clonetech_CST_IpVsInput_n6_annotated.Rdata")
