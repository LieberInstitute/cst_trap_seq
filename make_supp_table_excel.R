######
library(readxl)
library(readr)
library(xlsx)

dir.create("supp_tables/")

## main tables
files = c(Bulk_DE = "bulk_analysis_genotype/tables/all_genes_voom_sva_trkB_CST.csv",
	Bulk_GO = "bulk_analysis_genotype/tables/bulk_genotype_GOsets_hypergeo_Gene-p005.csv",
	IPvsInput_DE = "ip_vs_input_analysis/tables/all_genes_voom_CST_IPvsInput_lmer.csv",
	IPvsInput_GO = "ip_vs_input_analysis/tables/all_genes_voom_CST_IPvsInput_lmer.csv",
	IP_Geno_DE = "ip_analysis_genotype/tables/mutantVsWT_statistics_all.csv.gz",
	IP_Geno_GO = "ip_analysis_genotype/tables/mutVsWt_GO_FDR05.csv")
	
fileOut = c("Extended_Data_Figure1-1", "Extended_Data_Figure1-2", 
	"Extended_Data_Figure2-1", "Extended_Data_Figure2-3",
	"Extended_Data_Figure3-1", "Extended_Data_Figure3-2")

## de lists
deColnames = c("Symbol",  "logFC",  "t", "P.Value",
	"adj.P.Val", "B", "gene_type","EntrezID",
	"AveExpr","Length", "ensemblID")
for(i in c(1,3,5)) {
	x = read.csv(files[i], row.names=1,as.is=TRUE)
	x = x[,deColnames]
	write.csv(x, paste0("supp_tables/", fileOut[i], ".csv"),row.names=FALSE)
}

## GO lists
for(i in c(2,4,6)) {
	x = read.csv(files[i], as.is=TRUE)
	write.csv(x, paste0("supp_tables/", fileOut[i], ".csv"),row.names=FALSE)
}

#####################
## sfari stuff
sfariRdas = c(Bulk = "bulk_analysis_genotype/tables/SFARI_annotated_results.rda",
	IPvsInput = "ip_vs_input_analysis/tables/SFARI_annotated_results.rda",
	IP_Geno = "ip_analysis_genotype/tables/SFARI_annotated_results.rda")
	
## load in
humanSfariList = lapply(sfariRdas, function(x) {
	load(x) 
	humanSFARI
})

mouseSfariList = lapply(sfariRdas, function(x) {
	load(x) 
	mouseSFARI
})

## filter by project
humanSfariList$Bulk$Sig = humanSfariList$Bulk$adj.P.Val < 0.1
mouseSfariList$Bulk$Sig = mouseSfariList$Bulk$adj.P.Val < 0.1
humanSfariList$IPvsInput$Sig = humanSfariList$IPvsInput$Bonf < 0.05 &
	humanSfariList$IPvsInput$logFC > 0
mouseSfariList$IPvsInput$Sig = mouseSfariList$IPvsInput$Bonf < 0.05 &
	mouseSfariList$IPvsInput$logFC > 0
humanSfariList$IP_Geno$Sig = humanSfariList$IP_Geno$adj.P.Val < 0.05
mouseSfariList$IP_Geno$Sig = mouseSfariList$IP_Geno$adj.P.Val < 0.05

## merge
names(humanSfariList) = paste0(names(humanSfariList), "_human")
names(mouseSfariList) = paste0(names(mouseSfariList), "_mouse")
sfariList = c(humanSfariList, mouseSfariList)
sfariList = sfariList[c(4,1,5,2,6,3)] ## reorder

## get unique genes
sigList = lapply(sfariList, function(x) x[x$Sig,])
sapply(sigList,nrow)
uGenes = unique(unlist(sapply(sigList, function(x) x$gencodeID)))
length(uGenes)

## write otu
sfariMat = sapply(sigList, function(x) uGenes %in% x$gencodeID)
rownames(sfariMat) = uGenes
sfariDat = as.data.frame(sfariMat)
geneSym = do.call("rbind", lapply(sigList, function(x) x[,c("gencodeID", "Symbol")]))
sfariDat$Symbol = geneSym$Symbol[match(rownames(sfariDat), geneSym$gencodeID)]

sfariDat = sfariDat[,c(7,1:6)]
write.csv(sfariDat, "supp_tables/Extended_Data_Figure3-3.csv")

#######################
#### Harmonizome data #
#######################

harmFiles = c(Bulk = "bulk_analysis_genotype/tables/Harmonizome_CTD-Dx_CST_bulk_effects.csv",
	IPvsInput = "ip_vs_input_analysis/tables/Harmonizome_CTD-Dx_CST_effects.csv",
	IP_Geno = "ip_analysis_genotype/tables/Harmonizome_CST_IP-Genotype_effects.csv")
harmList = lapply(harmFiles, read.csv,as.is=TRUE, row.names=1)
for(i in seq(along=harmList)) {
	write.xlsx(harmList[[i]], file = "supp_tables/Extended_Data_Figure3-4.xlsx",
		sheetName = names(harmList)[i],append=TRUE)
}

###############
## zip everything
zip("maynardKardian_CST_Extended-Data.zip", 
	files = list.files("supp_tables", full.names=TRUE))