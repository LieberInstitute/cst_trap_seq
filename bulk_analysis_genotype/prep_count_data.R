############
library(SummarizedExperiment)
library(jaffelab)

## make directory
dir.create("count_data")

## load rses
load("preprocessed_data/rse_gene_TrkB_CST_Cortex_n10.Rdata")
load("preprocessed_data/rse_exon_TrkB_CST_Cortex_n10.Rdata")
load("preprocessed_data/rse_jx_TrkB_CST_Cortex_n10.Rdata")
load("preprocessed_data/rse_tx_TrkB_CST_Cortex_n10.Rdata")

# read in pheno
pd = data.frame(SampleID = paste0("CST_Mouse", 1:10),
	Genotype = rep(c("WT", "CST"), each=5),
	stringsAsFactors=FALSE)
rownames(pd) = pd$SampleID
pd$Genotype =factor(pd$Genotype, levels = c("WT", "CST"))	
pheno = cbind(pd[colnames(rse_gene),], colData(rse_gene))
pheno$SampleID = NULL

colData(rse_gene) = pheno
colData(rse_exon) = pheno
colData(rse_jx) = pheno
colData(rse_tx) = pheno

## save
save(rse_gene, file="count_data/rse_gene_TrkB_CST_Cortex_n10_annotated.Rdata")
save(rse_exon, file="count_data/rse_exon_TrkB_CST_Cortex_n10_annotated.Rdata")
save(rse_jx, file="count_data/rse_jx_TrkB_CST_Cortex_n10_annotated.Rdata")
save(rse_tx, file="count_data/rse_tx_TrkB_CST_Cortex_n10_annotated.Rdata")
