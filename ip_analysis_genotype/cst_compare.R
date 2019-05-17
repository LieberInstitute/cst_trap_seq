
library(SummarizedExperiment)
library(readxl)
library(jaffelab)
library(edgeR)
library(limma)
library(recount)
library(clusterProfiler)
library(org.Mm.eg.db)
library(VariantAnnotation)

######################
### NEW DATA #########
######################

### load in data
load("count_data/rse_gene_SoLo_Cst_MiSeq_n12.Rdata")

## add phenotype data
dat = read_excel("CstTrkBSamples_For Martha.xlsx")
colnames(dat) = gsub(" ", "_", colnames(dat))
colnames(dat) = gsub("'s", "", colnames(dat))
dat$Sample_ID = gsub(" ", "-", dat$Sample_ID)

## match up
mm = match(colnames(rse_gene), dat$Sample_ID)
colData(rse_gene) = cbind(colData(rse_gene), as.data.frame(dat[mm,3:8]))
colnames(colData(rse_gene))[16:17] = c("Qubit_Conc", "KAPA_Conc")
rse_gene$Genotype = factor(ifelse(rse_gene$Sample_Genotype == "Cst-Cre", "Control", "Mutant"))

# subset controls
rse_gene_new = rse_gene[,rse_gene$Genotype=="Control"]

######################
### OLD DATA #########
######################

load("/dcl01/lieber/ajaffe/lab/oxt_trap_seq/preprocessed_data/Clonetech_CST/rse_gene_Clonetech_CST_HiSeq_n6.Rdata")
rse_gene_old = rse_gene[,grep("IP", colnames(rse_gene))]

######################
## COMBINE ###########
######################

n = intersect(colnames(colData(rse_gene_new)), colnames(colData(rse_gene_old)))
colData(rse_gene_new) = colData(rse_gene_new)[,n]
colData(rse_gene_old) = colData(rse_gene_old)[,n]
rowData(rse_gene_new)$meanExprs = NULL

mm = match(rownames(rse_gene_new), rownames(rse_gene_old))
table(diff(mm)) # same
rowData(rse_gene_old)$Symbol = rowData(rse_gene_new)$Symbol
rowData(rse_gene_old)$EntrezID = rowData(rse_gene_new)$EntrezID

rse_gene = cbind(rse_gene_new, rse_gene_old)
rse_gene$Experiment = rep(c("New","Old"),times=c(6,3))

################
### QC #######

## PCA
exprsIndex = rowMeans(getRPKM(rse_gene,"Length")) > 0.1
pca = prcomp(t(log2(getRPKM(rse_gene,"Length")[exprsIndex,]+1)))
pcaVars = getPcaVars(pca)

palette(brewer.pal(4,"Set1"))
plot(pca$x, xlab=paste0("PC1: ", pcaVars[1], "% Var Expl"), 
	ylab=paste0("PC2: ", pcaVars[2], "% Var Expl"), pch = 21, 
	bg = factor(rse_gene$Experiment))
	
plot(pca$x[,1] ~ rse_gene$totalAssignedGene) # definitely PC1

################
rse_gene = rse_gene[exprsIndex,]

## modeling
mod = model.matrix(~Experiment + totalAssignedGene, 
	data=colData(rse_gene))
dge = DGEList(counts = assays(rse_gene)$counts, 
	genes = rowData(rse_gene))
dge = calcNormFactors(dge)

## mean-variance
vGene = voom(dge,mod,plot=TRUE)
fitGene = lmFit(vGene)
eBGene = eBayes(fitGene)
outGene = topTable(eBGene,coef=2,
	p.value = 1,number=nrow(rse_gene),sort="none")
sum(outGene$adj.P.Val < 0.05)
sigGene = outGene[outGene$adj.P.Val < 0.05,]
sigGene = sigGene[order(sigGene$P.Value),]
vars = c("gencodeID","Symbol","logFC","t", "P.Value", "adj.P.Val",
		"B", "ensemblID", "Length", "gene_type", "EntrezID", "AveExpr")

write.csv(sigGene[,vars], file="tables/CST_compare_CloneVsSolo.csv",
	row.names=FALSE)