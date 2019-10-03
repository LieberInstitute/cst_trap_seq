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
library(readxl)

## load counts
load("count_data/rse_gene_Clonetech_CST_IpVsInput_n6_annotated.Rdata")

## outputs
dir.create("rdas")
dir.create("tables")
dir.create("plots")

#### explore
geneExprs = log2(getRPKM(rse_gene,"Length")+1)
pca = prcomp(t(geneExprs))
pcaVars = getPcaVars(pca)
plot(pca$x, pch = 21, bg=factor(rse_gene$Fraction))

plot(pca$x[,1] ~ rse_gene$Fraction)
plot(pca$x[,2] ~ rse_gene$overallMapRate) # yes
boxplot(rse_gene$overallMapRate ~rse_gene$Fraction ) # no assoc

#########
## DE ###
#########

## GENE
rse_gene_filter = rse_gene[rowMeans(getRPKM(rse_gene, "Length")) > 0.1,]
nrow(rse_gene_filter) # 21776
mod = model.matrix(~Fraction+overallMapRate,
	data = colData(rse_gene_filter))
	
dge = DGEList(counts = assays(rse_gene_filter)$counts, 
	genes = rowData(rse_gene_filter))
dge = calcNormFactors(dge)

## mean-variance
vGene = voom(dge,mod,plot=TRUE)

## duplicate correlation?
rse_gene$Mouse = ss(colnames(rse_gene), "_", 2)
geneCor = duplicateCorrelation(vGene$E, design = mod, block = rse_gene$Mouse)
fitGene = lmFit(vGene, design = mod,
	correlation =  geneCor$consensus.correlation,block = rse_gene$Mouse)
## or regular?
# fitGene = lmFit(vGene)
ebGene = eBayes(fitGene)
outGene = topTable(ebGene, coef=2, sort="none", n = nrow(dge))
sum(outGene$adj.P.Val < 0.05) ## 5362
outGene$Bonf = p.adjust(outGene$P.Value, "bonf")
sum(outGene$Bonf < 0.05) ## 868

table(abs(outGene$logFC[outGene$Bonf < 0.05]) > 1)
write.csv(outGene, file = "tables/all_genes_voom_CST_IPvsInput_lmer.csv")

####################
## volcano	plot ###
outGene$sigColor = as.numeric(outGene$Bonf < 0.05)+1

## volanco plot
g = c("Syt2", "Nxph1", "Sptan1", "Shank2", "Mal", "Apoe")
m = match(g, outGene$Symbol)

pdf("plots/figure2c_IPvsInput_volcano_plot.pdf",useDingbats=FALSE)
palette(brewer.pal(5, "Dark2"))
par(mar = c(5,6,2,2), cex.axis=2,cex.lab=2)
plot(-log10(P.Value) ~ logFC, pch = 21, bg=sigColor, 
	data = outGene, xlab = "CST IP vs Input log2FC")
abline(h = -log10(max(outGene$P.Value[outGene$sigColor==2])), lty=2,lwd=2)
shadowtext(outGene$logFC[m]+0.35, -log10(outGene$P.Value[m]),
	letters[12:17],font=2,cex=1.25,col="grey")
dev.off()

#####################
## doyle compare ####
#####################

cst_pub = read_excel("../gene_sets/1-s2.0-S0092867408013664-mmc2.xls", 
	sheet = "Cortical Cortistatin Neurons")
colnames(cst_pub)[1] = "IP_Ratio"
cst_pub = as.data.frame(cst_pub)

## line back up
matchCst = match(cst_pub$Symbol, outGene$Symbol)
cst_pub$logFC_cst = outGene$logFC[matchCst]
cst_pub$P.Value_cst = outGene$P.Value[matchCst]
cst_pub$log2_IP_ratio = log2(cst_pub$IP_Ratio)

plot(cst_pub$log2_IP_ratio, cst_pub$logFC_cst)
## looks bad
table(!is.na(matchCst))

#########################
##### gene ontology #####
sIndex=  which(outGene$Bonf < 0.05)
sigGene = split(outGene$EntrezID[sIndex], sign(outGene$logFC[sIndex]))
sigGene = lapply(sigGene, function(x) x[!is.na(x)])

geneUniverse = as.character(outGene$EntrezID)
geneUniverse = geneUniverse[!is.na(geneUniverse)]

go <- compareCluster(sigGene, fun = "enrichGO",
                universe = geneUniverse, OrgDb = org.Mm.eg.db,
                ont = "ALL", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1,
				readable= TRUE)
save(go, file = "rdas/GO_voomBonf_IPvsInput_signed.rda")

## check sets ########
goDf = as.data.frame(go)
goDf$Cluster = ifelse(goDf$Cluster == "-1", "Input", "IP")
goDf = goDf[order(goDf$pvalue),]
write.csv(goDf, row.names=FALSE,
	file="tables/GO_voomBonf_IPvsInput_signed.csv")

#### plots #####
goDf = read.csv("tables/GO_voomBonf_IPvsInput_signed.csv")

setsToPlot = c("GO:0033267", "GO:0044306",
	"GO:0042063", "GO:0050808","GO:0098918",
	"GO:0051015", "GO:0021800", "GO:0005509",
	"GO:0043209")
	
goExample = goDf[goDf$Cluster == "IP",]
goExample = goExample[match(setsToPlot, goExample$ID),]	
goDown = goDf[goDf$Cluster == "Input",]
goExample$input_pAdj = goDown$p.adjust[match(goExample$ID, goDown$ID)]
goExample$input_pAdj[is.na(goExample$input_pAdj)] = 1
goExample = goExample[order(goExample$pvalue),]

negFdrMat = -log10(as.matrix(goExample[,c("p.adjust", "input_pAdj")]))
rownames(negFdrMat) = goExample$Description

pdf("plots/go_figure2d_barplot_ipVsInput.pdf",h=4,w=7)
par(mar=c(5,18,2,2),cex.axis=1.2,cex.lab=1.5)
barplot(t(negFdrMat),width=0.5,# names = goExample$Description,
	horiz=TRUE, beside=TRUE, col = 2:1,
	xlab="-log10(FDR)",las=1)
abline(v=-log10(0.05), col="blue")
dev.off()


###############
### CSEA ######
###############
library(pSI)
library(pSI.data)

upGenes = outGene$Symbol[outGene$logFC > 0 & outGene$Bonf < 0.05]
csea_IP = fisher.iteration(mouse$psi.out, upGenes)
write.csv(csea_IP, file = "tables/csea_enrichment_IPenriched.csv")

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
outGene$hsapien_homolog = MMtoHG$hsapiens_homolog_ensembl_gene[
	match(outGene$ensemblID,MMtoHG$ensembl_gene_id)]
outGene$sigColor = NULL 

########################
# load SFARI human genes
humanSFARI = read.csv('../gene_sets/SFARI-Gene_genes_05-06-2019release_05-16-2019export.csv')
humanSFARI = cbind(humanSFARI,outGene[lookfor(humanSFARI$gene.symbol,outGene$Symbol),])
humanSFARI= humanSFARI[!is.na(humanSFARI$Symbol),]
humanSFARI = humanSFARI[!duplicated(humanSFARI),]
nrow(humanSFARI) # 937 expressed in mouse cst ip dataset

#########################
# load SFARI mouse models
mouseSFARI = read.csv('../gene_sets/SFARI-Gene_animal-genes_05-06-2019release_05-16-2019export.csv')
mouseSFARI = with(mouseSFARI,mouseSFARI[model.species=='Mus musculus',])
mouseSFARI = cbind(mouseSFARI,outGene[lookfor(mouseSFARI$gene.symbol,outGene$Symbol),])
mouseSFARI= mouseSFARI[!is.na(mouseSFARI$Symbol),]
mouseSFARI = mouseSFARI[order(match(mouseSFARI$Symbol,humanSFARI$Symbol)),]
mouseSFARI = mouseSFARI[!duplicated(mouseSFARI),]
rownames(mouseSFARI) = mouseSFARI$Symbol
nrow(mouseSFARI) # 239 expressed in mouse oxt dataset

###########################
# list of DEG in mouse SFARI
outGene$inMouseSFARI = outGene$Symbol %in% mouseSFARI$Symbol
(t1 = with(outGene,table(inMouseSFARI,inDEG = Bonf < 0.05 & logFC > 0)))
fisher.test(t1) # OR = 6.340408 p-value = 9.247e-13
ind1 = which(mouseSFARI$Bonf < 0.05 & mouseSFARI$logFC > 0)
nrow(mouseSFARI[ind1,])

#######################################
# list of DEG in scored human SFARI list
outGene$inHumanSFARI =  toupper(outGene$Symbol) %in% humanSFARI$gene.symbol
outGene_hs = outGene[grep("^ENSG", outGene$hsapien_homolog),]
(t2 = with(outGene_hs,table(inHumanSFARI,inDEG = Bonf < 0.05 & logFC > 0)))
fisher.test(t2) # OR = 4.175872, p-value = 9.026209e-25
ind2 = which(humanSFARI$Bonf < 0.05 & humanSFARI$logFC > 0)
humanSFARI[ind2,]
nrow(humanSFARI[ind2,])
save(mouseSFARI, humanSFARI, file = "tables/SFARI_annotated_results.rda")

#######################
### other DX ##########
#######################
library(biomaRt)
ensembl = useMart("ENSEMBL_MART_ENSEMBL",  
	dataset="hsapiens_gene_ensembl", host="jul2016.archive.ensembl.org")
sym = getBM(attributes = c("ensembl_gene_id","hgnc_symbol","entrezgene"), 
		values=outGene$hsapien_homolog, mart=ensembl)
outGene$hsapien_EntrezID = sym$entrezgene[match(outGene$hsapien_homolog, sym$ensembl_gene_id)]

dxSets = read.delim("../gene_sets/Harmonizome_CTD Gene-Disease Associations Dataset.txt",
	as.is=TRUE,skip=1)
dxSets$isExpMouse = dxSets$GeneID %in% outGene$hsapien_EntrezID
dxSets$CST_enrich_Mouse = dxSets$GeneID %in% outGene$hsapien_EntrezID[
	outGene$Bonf <0.05 & outGene$logFC > 0]
## split
dxSetsList = split(dxSets, dxSets$Disease)
dxSetsList = dxSetsList[sapply(dxSetsList, function(x) sum(x$CST_enrich_Mouse) > 0)]
length(dxSetsList)
univ = outGene[order(outGene$P.Value),]
univ = univ[!is.na(univ$hsapien_EntrezID),]
univ = univ[!duplicated(univ$hsapien_EntrezID),]
univ$CST = univ$Bonf < 0.05 & univ$logFC > 0

dxStats = do.call("rbind", mclapply(dxSetsList, function(x) {
	x = x[x$isExpMouse,]
	inSet = univ$hsapien_EntrezID %in% x$GeneID
	tt = table(inSet, univ$CST)
	data.frame(OR = getOR(tt), p.value = chisq.test(tt)$p.value)
},mc.cores=4))
	
dxStats = dxStats[order(dxStats$p.value),]

## negative control
univ$nonCST = univ$Bonf < 0.05 & univ$logFC < 0

dxStatsDepl = do.call("rbind", mclapply(dxSetsList, function(x) {
	x = x[x$isExpMouse,]
	inSet = univ$hsapien_EntrezID %in% x$GeneID
	tt = table(inSet, univ$nonCST)
	data.frame(OR = getOR(tt), p.value = chisq.test(tt)$p.value)
},mc.cores=4))
dxStatsDepl = dxStatsDepl[rownames(dxStats),]

dxStats = cbind(dxStats, dxStatsDepl)
colnames(dxStats) = c("Enrich_OR" ,"Enrich_Pval", "Deplete_OR", "Deplete_Pval")
dxStats$setSize = sapply(dxSetsList,nrow)[rownames(dxStats)]

write.csv(dxStats, "tables/Harmonizome_CTD-Dx_CST_effects.csv")

## compare
