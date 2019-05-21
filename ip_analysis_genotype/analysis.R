###
library(SummarizedExperiment)
library(readxl)
library(jaffelab)
library(edgeR)
library(limma)
library(recount)
library(TeachingDemos) # for shadow text
library(clusterProfiler)
library(org.Mm.eg.db)
library(VariantAnnotation)
library(RColorBrewer)

dir.create("plots")
dir.create("tables")

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

#########
## QC ###
#########

## sex
rse_gene$totalY = colSums(getRPKM(rse_gene,"Length")[which(seqnames(rowRanges(rse_gene)) == "chrY"),])
rse_gene$xist = log2(getRPKM(rse_gene,"Length")[which(rowData(rse_gene)$Symbol == "Xist"),]+1)
rse_gene$Sex = ifelse(rse_gene$xist > 1, "F", "M")

## spot check
as.data.frame((colData(rse_gene)[,c("totalY", "xist", "Julie_Tube_Code")]))

plot(rse_gene$totalY, rse_gene$xist, pch =21 ,bg=rse_gene$Genotype,
	xlab="total chrY exprs", ylab="Xist log2 exprs")
legend("top", levels(rse_gene$Genotype), pch = 15, col=1:2)

## PCA
exprsIndex = rowMeans(getRPKM(rse_gene,"Length")) > 0.1
pca = prcomp(t(log2(getRPKM(rse_gene,"Length")[exprsIndex,]+1)))
pcaVars = getPcaVars(pca)

palette(brewer.pal(4,"Dark2"))
plot(pca$x, xlab=paste0("PC1: ", pcaVars[1], "% Var Expl"), 
	ylab=paste0("PC2: ", pcaVars[2], "% Var Expl"), pch = 21, 
	bg = rse_gene$Genotype)
	
plot(pca$x[,1] ~ rse_gene$totalAssignedGene) # definitely PC1
plot(pca$x[,2] ~ rse_gene$mitoRate) #  PC2?
plot(pca$x[,2] ~ rse_gene$KAPA_Conc) #  PC2? sample had highest chrM and highest KAPA

################
### genotypes ##
################
vcf = readVcf("preprocessed_data/Genotypes/mergedVariants.Ntrk2.vcf.gz", "mm10")
colnames(vcf) = ss(ss(colnames(vcf), "/",9),"_")
snpsGeno = geno(vcf)$GT
snpsGeno[snpsGeno == "."] = NA
snpsGeno[snpsGeno == "0/0"] = 0
snpsGeno[snpsGeno == "0/1"] = 1
snpsGeno[snpsGeno == "1/1"] = 2
class(snpsGeno) = "numeric"

table(rowSums(is.na(snpsGeno))) # revisit

###################
## analysis #######
###################

exprsIndex = rowMeans(getRPKM(rse_gene,"Length")) > 0.1
rse_gene = rse_gene[exprsIndex,]

## modeling
mod = model.matrix(~Genotype + totalAssignedGene, 
	data=colData(rse_gene))

##### DGE ######
dge = DGEList(counts = assays(rse_gene)$counts, 
	genes = rowData(rse_gene))
dge = calcNormFactors(dge)

## mean-variance
vGene = voom(dge,mod,plot=TRUE)

## do analysis
fitGene = lmFit(vGene)

## top table
eBGene = eBayes(fitGene)
outGene = topTable(eBGene,coef=2,
	p.value = 1,number=nrow(rse_gene), sort="none")
outGene$sigColor = as.numeric(outGene$adj.P.Val < 0.05)+1
sigGene = outGene[outGene$adj.P.Val < 0.05,]
sigGene = sigGene[order(sigGene$P.Value),]
statsOut = sigGene[,c(2,5, 8,10:13,3,1,4,6,9)]
write.csv(statsOut, file = "tables/mutantVsWT_statistics_FDR05.csv",row.names=FALSE)
write.csv(outGene[,c(2,5, 8,10:13,3,1,4,6,9)], 
	file = gzfile("tables/mutantVsWT_statistics_all.csv.gz"),row.names=FALSE)
save(outGene, file="tables/mutantVsWT_statistics_all.rda")

dim(sigGene)

#######
## volanco plot
g = c("Cort", "Ntrk2")
m = match(g, outGene$Symbol)


pdf("plots/volcano_plot_mutVsWt.pdf",useDingbats=FALSE)
palette(brewer.pal(5, "Dark2"))
par(mar = c(5,6,2,2), cex.axis=2,cex.lab=2)
plot(-log10(P.Value) ~ logFC, pch = 21, bg=sigColor, 
	data = outGene, xlab = "MUT vs WT log2FC")
shadowtext(outGene$logFC[m], -log10(outGene$P.Value[m]),
	letters[21:22],font=2,cex=2,col="grey")
abline(v=c(-1,1), lty=2,lwd=2)
dev.off()


## volcano	again
g2 = c("Wt1", "Calb1", "Lgals1", "Trpc6", "Syt6", "Gng4")
m2 = match(g2, outGene$Symbol)

pdf("plots/volcano_plot_validatedGenes.pdf",useDingbats=FALSE)
palette(brewer.pal(5, "Dark2"))
par(mar = c(5,6,2,2), cex.axis=2,cex.lab=2)
plot(-log10(P.Value) ~ logFC, pch = 21, bg=sigColor, 
	data = outGene, xlab = "MUT vs WT log2FC")
shadowtext(outGene$logFC[m2], -log10(outGene$P.Value[m2]),
	letters[21:26],font=2,cex=2,col="grey")
abline(v=c(-1,1), lty=2,lwd=2)
dev.off()

####################
## gene ontology ###
####################


################
## gene set ####
sigGeneList = split(sigGene$EntrezID, sign(sigGene$logFC))
sigGeneList = lapply(sigGeneList, function(x) x[!is.na(x)])
geneUniverse = as.character(outGene$EntrezID)
geneUniverse = geneUniverse[!is.na(geneUniverse)]

go<- compareCluster(sigGeneList, fun = "enrichGO",
                universe = geneUniverse, OrgDb = org.Mm.eg.db,
                ont = "ALL", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1,
				readable= TRUE)
	
## write out				
goOut = as.data.frame(go)
goSig = goOut[goOut$p.adjust < 0.05,]
colnames(goSig)[1] = "Direction"
goSig = goSig[order(goSig$p.adjust),]
save(go, goOut,goSig, file = "tables/mutVsWt_GO_FDR05.rda")
write.csv(goSig, file = "tables/mutVsWt_GO_FDR05.csv",row.names=FALSE)

## plots
goIDs = c("GO:0033267", "GO:0007409", "GO:0061564",	
	"GO:0010976", "GO:0006816", "GO:0006874", "GO:0030003",
	"GO:1904062", "GO:0005509", "GO:0034702")
goPlot = goOut[goOut$ID %in% goIDs,]

goExample = goPlot[!duplicated(goPlot[,c("ID", "Description")]),3:4]
upGo = goPlot[goPlot$Cluster == 1,]
goExample$Up = upGo$pvalue[match(goExample$ID, upGo$ID)]
downGo = goPlot[goPlot$Cluster == -1,]
goExample$Down = downGo$pvalue[match(goExample$ID, downGo$ID)]

goExample$Description[goExample$Description == "regulation of cation transmembrane transport"] = "regulation of cation\ntransmembrane transport"
goExample$Description[goExample$Description == "positive regulation of neuron projection development"] = "positive regulation of\nneuron projection development"

goExample$Label = paste0(goExample$ID, ": ", goExample$Description)
pdf("plots/go_barplot.pdf",h=6,w=8)
par(mar=c(5,21.5,1,1),cex.axis=1.2,cex.lab=1.5)
barplot(t(-log10(as.matrix(goExample[,c("Up", "Down")]))),
	width=0.75, names = goExample$Label,
	horiz=TRUE,xlim=c(0,8),ylim=c(0.5,23.5),
	xlab="-log10(P-Value)",las=1,beside=TRUE, col=c("blue","red"))
abline(v=-log10(max(goOut$pvalue[goOut$p.adj < 0.05])), col="blue")
legend("topright", c("Mut>Wt", "Mut<Wt"), col=c("blue","red"),
	cex=1.2,nc=1,pch=15)
dev.off()
	
goSub = go
goSub@compareClusterResult = goSub@compareClusterResult[goSub@compareClusterResult$ID %in% goIDs,]
goSub@compareClusterResult$ONTOLOGY = NULL

pdf("plots/go_dotplot.pdf")
dotplot(goSub)
dev.off()