######

library(SummarizedExperiment)
library(recount)
library(org.Mm.eg.db)
library(clusterProfiler)
library(edgeR)
library(limma)
library(jaffelab)
library(sva)
library(RColorBrewer)
library(TeachingDemos) # for shadow text
library(biomaRt)

## load counts
load("count_data/rse_gene_TrkB_CST_Cortex_n10_annotated.Rdata")

## outputs
dir.create("rdas")
dir.create("tables")
dir.create("plots")

#### explore
geneExprs = log2(getRPKM(rse_gene,"Length")+1)
pca = prcomp(t(geneExprs))
pcaVars = getPcaVars(pca)
plot(pca$x, pch = 21, bg=factor(rse_gene$Genotype))

plot(pca$x[,1] ~ rse_gene$Genotype)
plot(pca$x[,1] ~ rse_gene$mitoRate) # kind of
plot(pca$x[,2] ~ rse_gene$totalAssignedGene) # kind of

## check sex
pdf("plots/sex_check.pdf")
par(mar=c(5,6,2,2), cex.axis=2, cex.lab=2)
plot(geneExprs[rownames(rse_gene)[which(rowData(rse_gene)$Symbol == "Xist")],],
	colSums(geneExprs[as.character(seqnames(rse_gene)) == "chrY",]),
	xlab="Xist", ylab="chrY", pch=21,bg="grey",cex=1.6)
abline(h=15,lty=2)
text(x=1, y=15.5, "Male",cex=1.6)
text(x=1, y=12.5, "Female",cex=1.6)
dev.off()

# use 15 as a cutoff
rse_gene$predSex = 	ifelse(colSums(geneExprs[as.character(seqnames(rse_gene)) == "chrY",]) > 15, "Male", "Female")

## check Bdnf
bdnf = geneExprs[rownames(rse_gene)[which(rowData(rse_gene)$Symbol == "Bdnf")],]
pdf("plots/bdnf_vs_genotype.pdf")
par(mar=c(5,6,2,2), cex.axis=2, cex.lab=2)
boxplot(bdnf ~ rse_gene$Genotype,las=3, ylab = "Bdnf: log2(RPKM+1)")
dev.off()

## hc
dd = dist(t(geneExprs))
myplclust(hclust(dd),lab.col = as.numeric(rse_gene$Genotype))

#########
## DE ###
#########

## GENE
rse_gene_filter = rse_gene[rowMeans(getRPKM(rse_gene, "Length")) > 0.1,]
mod = model.matrix(~Genotype+totalAssignedGene+mitoRate,
	data = colData(rse_gene_filter))
	
dge = DGEList(counts = assays(rse_gene_filter)$counts, 
	genes = rowData(rse_gene_filter))
dge = calcNormFactors(dge)

## mean-variance
vGene = voom(dge,mod,plot=TRUE)
fitGene = lmFit(vGene, mod)
ebGene = eBayes(fitGene)
outGene = topTable(ebGene, coef=2, sort="none", n = nrow(dge))
sum(outGene$adj.P.Val < 0.1)
sum(outGene$P.Value < 0.001)
sum(outGene$P.Value < 0.005)
hist(outGene$P.Value)
plot(outGene$logFC, -log10(outGene$P.Value))
sigGene = topTable(ebGene , coef=2,
	p.value = 0.2, n = nrow(dge))

## with sva
svaobj = sva(vGene$E, mod, mod0 = model.matrix(~totalAssignedGene+mitoRate, 
						data=colData(rse_gene_filter)))
vGeneSva = voom(dge,cbind(mod, svaobj$sv),plot=TRUE)
fitGeneSva =  lmFit(vGeneSva, cbind(mod, svaobj$sv))
ebGeneSva = eBayes(fitGeneSva)
outGeneSva = topTable(ebGeneSva, coef=2, sort="none", n = nrow(dge))
sum(outGeneSva$adj.P.Val < 0.1)
sum(outGeneSva$adj.P.Val < 0.1)
sum(outGeneSva$adj.P.Val < 0.2)
hist(outGeneSva$P.Value)
plot(outGeneSva$logFC, -log10(outGeneSva$P.Value))
sum(outGeneSva$P.Value < 0.001)
sum(outGeneSva$P.Value < 0.005)

sigGeneSva = topTable(ebGeneSva, coef=2,
	p.value = 0.2, n = nrow(dge))

write.csv(outGeneSva, file = "tables/all_genes_voom_sva_trkB_CST.csv")
write.csv(sigGeneSva, file = "tables/de_genes_voom_sva_trkB_CST_fdr20.csv")

####################
## volcano	plot ###
outGeneSva$sigColor = as.numeric(outGeneSva$adj.P.Val < 0.1)+1

## volanco plot
g = c("Scg2", "Tnmd", "Matn2","Npy", "Syt12","Nptx2", "Chrna4")
m = match(g, outGeneSva$Symbol)


pdf("plots/figure1b_bulk_genotype_volcano_plot.pdf",useDingbats=FALSE)
palette(brewer.pal(5, "Dark2"))
par(mar = c(5,6,2,2), cex.axis=2,cex.lab=2)
plot(-log10(P.Value) ~ logFC, pch = 21, bg=sigColor, 
	data = outGeneSva, xlab = "CST vs WT log2FC")
abline(h = -log10(max(outGeneSva$P.Value[outGeneSva$sigColor==2])), lty=2,lwd=2)
shadowtext(outGeneSva$logFC[m]+0.35, -log10(outGeneSva$P.Value[m]),
	letters[18:24],font=2,cex=1.25,col="grey")
dev.off()


#########################
##### gene ontology #####
sIndex=  which(outGeneSva$P.Value< 0.005)
sigGene = split(outGeneSva$EntrezID[sIndex], sign(outGeneSva$logFC[sIndex]))
sigGene = lapply(sigGene, function(x) x[!is.na(x)])

geneUniverse = as.character(outGeneSva$EntrezID)
geneUniverse = geneUniverse[!is.na(geneUniverse)]

go <- compareCluster(sigGene, fun = "enrichGO",
                universe = geneUniverse, OrgDb = org.Mm.eg.db,
                ont = "ALL", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1,
				readable= TRUE)
	
##########
## gsea ## 
geneStat = outGeneSva$t
names(geneStat) = outGeneSva$EntrezID
geneStat = sort(geneStat, decreasing=TRUE)
geneStat = geneStat[!is.na(names(geneStat))]

goGse = gseGO(geneStat, OrgDb = org.Mm.eg.db, ont = "ALL",
	nPerm = 100000, minGSSize   = 20, maxGSSize    = 500,
    pvalueCutoff = 0.05, verbose=TRUE)
goGseDf = DataFrame(as.data.frame(goGse))
goGseDf$core_enrichment = CharacterList(strsplit(goGseDf$core_enrichment, "/"))
goGseDf$core_enrichment = endoapply(goGseDf$core_enrichment, function(x) {
	outGeneSva$Symbol[match(x, outGeneSva$EntrezID)]
})
goGseDf$leading_edge = CharacterList(strsplit(goGseDf$leading_edge, ", "))

keggGse = gseKEGG(geneStat, organism = "mmu",
	nPerm = 100000, minGSSize    = 20, maxGSSize  = 500,
    pvalueCutoff = 0.05, verbose=TRUE)
keggGseDf = DataFrame(as.data.frame(keggGse))
keggGseDf$core_enrichment = CharacterList(strsplit(keggGseDf$core_enrichment, "/"))
keggGseDf$core_enrichment = endoapply(keggGseDf$core_enrichment, function(x) {
	outGeneSva$Symbol[match(x, outGeneSva$EntrezID)]
})
keggGseDf$leading_edge = CharacterList(strsplit(keggGseDf$leading_edge, ", "))
save(go, goGse, goGseDf, keggGse, keggGseDf, file = "rdas/geneSet_runs_geneLevel_limma_sva.rda")

#########################3
## check sets ########
goDf = as.data.frame(go)
goDf = goDf[order(goDf$pvalue),]
goDf$Cluster = ifelse(goDf$Cluster == "-1", "CST<WT", "WT<CST")
write.csv(goDf, row.names=FALSE,
	file="tables/bulk_genotype_GOsets_hypergeo_Gene-p005.csv")

keggGseDf$Description[keggGseDf$NES > 0]
goGseDf$Description[goGseDf$NES > 0]

## write out
keggGseDf$core_enrichment = sapply(keggGseDf$core_enrichment,paste,collapse=";")
goGseDf$core_enrichment = sapply(goGseDf$core_enrichment,paste,collapse=";")
write.table(keggGseDf, file ="tables/KEGG_gsea.tsv",
	sep="\t",row.names=FALSE)
write.table(goGseDf, file ="tables/GO_gsea.tsv",sep="\t",row.names=FALSE)

###########################
### autism gene sets ######
###########################
## via badoi
lookfor= function(this,inThat) { #gene name matching
  this = toupper(this); inThat = toupper(inThat);
  tmp = sapply(this,function(x) grep(paste0('^',x,'$'),inThat))
  return(sapply(tmp,function(x) ifelse(length(x)==0,NA,x[1])))}

ensembl = useMart("ensembl",dataset = 'mmusculus_gene_ensembl')
MMtoHG = getBM(attributes = c('ensembl_gene_id','hsapiens_homolog_ensembl_gene'),
               mart = ensembl)

## add mouse homologs
outGeneSva$hsapien_homolog = MMtoHG$hsapiens_homolog_ensembl_gene[
	match(outGeneSva$ensemblID,MMtoHG$ensembl_gene_id)]
outGeneSva$sigColor = NULL 
sigGene = outGeneSva[outGeneSva$adj.P.Val < 0.1,]

########################
# load SFARI human genes
humanSFARI = read.csv('../gene_sets/SFARI-Gene_genes_05-06-2019release_05-16-2019export.csv')
humanSFARI = cbind(humanSFARI,outGeneSva[lookfor(humanSFARI$gene.symbol,outGeneSva$Symbol),])
humanSFARI= humanSFARI[!is.na(humanSFARI$Symbol),]
humanSFARI = humanSFARI[!duplicated(humanSFARI),]
nrow(humanSFARI) # 939 expressed in mouse cst bulk dataset

#########################
# load SFARI mouse models
mouseSFARI = read.csv('../gene_sets/SFARI-Gene_animal-genes_05-06-2019release_05-16-2019export.csv')
mouseSFARI = with(mouseSFARI,mouseSFARI[model.species=='Mus musculus',])
mouseSFARI = cbind(mouseSFARI,outGeneSva[lookfor(mouseSFARI$gene.symbol,outGeneSva$Symbol),])
mouseSFARI= mouseSFARI[!is.na(mouseSFARI$Symbol),]
mouseSFARI = mouseSFARI[order(match(mouseSFARI$Symbol,humanSFARI$Symbol)),]
mouseSFARI = mouseSFARI[!duplicated(mouseSFARI),]
rownames(mouseSFARI) = mouseSFARI$Symbol
nrow(mouseSFARI) # 238 expressed in mouse oxt dataset

###########################
# list of DEG in mouse SFARI
outGeneSva$inMouseSFARI = outGeneSva$Symbol %in% mouseSFARI$Symbol
(t1 = with(outGeneSva,table(inMouseSFARI,inDEG = adj.P.Val < 0.1)))
fisher.test(t1) # OR = 5.861999 p-value = 0.05056
ind1 = which(mouseSFARI$adj.P.Val < 0.1)
mouseSFARI[ind1,]
#######################################
# list of DEG in scored human SFARI list
# 131 human SFARI ASD-linked genes in our list of DEGs
outGeneSva$inHumanSFARI =  toupper(outGeneSva$Symbol) %in% humanSFARI$gene.symbol
outGeneSva_hs = outGeneSva[grep("^ENSG", outGeneSva$hsapien_homolog),]
(t2 = with(outGeneSva_hs,table(inHumanSFARI,inDEG = adj.P.Val < 0.1)))
fisher.test(t2) # OR = 1.768907, p-value = 0.4207
ind2 = which(humanSFARI$adj.P.Val < 0.1)
humanSFARI[ind2,]
	
# Bdnf, Grin2b, Ntrk3


