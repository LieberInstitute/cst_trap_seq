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
library(biomaRt)
library(readxl)

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

## label
outGene$sigColor = as.numeric(outGene$adj.P.Val < 0.05)+1

sigGene = outGene[outGene$adj.P.Val < 0.05,]
sigGene = sigGene[order(sigGene$P.Value),]

## write out
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
pdf("plots/figure3d_GO_barplot_ipGenotype.pdf",h=6,w=8)
par(mar=c(5,21.5,1,1),cex.axis=1.2,cex.lab=1.5)
barplot(t(-log10(as.matrix(goExample[,c("Up", "Down")]))),
	width=0.75, names = goExample$Label,
	horiz=TRUE,xlim=c(0,8),ylim=c(0.5,23.5),
	xlab="-log10(P-Value)",las=1,beside=TRUE, col=c("blue","red"))
abline(v=-log10(max(goOut$pvalue[goOut$p.adj < 0.05])), col="blue")
legend("topright", c("Mut>Wt", "Mut<Wt"), col=c("blue","red"),
	cex=1.2,nc=1,pch=15)
dev.off()

###############
### CSEA ######
###############
library(pSI)
library(pSI.data)

## weird i have to force this
load("~/R/x86_64-pc-linux-gnu-library/3.5/pSI.data/data/mouse.rda")

## run enrichment
csea_IP = fisher.iteration(mouse$psi.out, outGene$Symbol[outGene$adj.P.Val < 0.05])
write.csv(csea_IP, file = "tables/csea_enrichment_genotypeEffect_IPdata.csv")


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
sigGene = outGene[outGene$adj.P.Val < 0.05,]

########################
# load SFARI human genes
humanSFARI = read.csv('../gene_sets/SFARI-Gene_genes_05-06-2019release_05-16-2019export.csv')
humanSFARI = cbind(humanSFARI,outGene[lookfor(humanSFARI$gene.symbol,outGene$Symbol),])
humanSFARI= humanSFARI[!is.na(humanSFARI$Symbol),]
humanSFARI = humanSFARI[!duplicated(humanSFARI),]
nrow(humanSFARI) # 917 expressed in mouse cst ip dataset

#########################
# load SFARI mouse models
mouseSFARI = read.csv('../gene_sets/SFARI-Gene_animal-genes_05-06-2019release_05-16-2019export.csv')
mouseSFARI = with(mouseSFARI,mouseSFARI[model.species=='Mus musculus',])
mouseSFARI = cbind(mouseSFARI,outGene[lookfor(mouseSFARI$gene.symbol,outGene$Symbol),])
mouseSFARI= mouseSFARI[!is.na(mouseSFARI$Symbol),]
mouseSFARI = mouseSFARI[order(match(mouseSFARI$Symbol,humanSFARI$Symbol)),]
mouseSFARI = mouseSFARI[!duplicated(mouseSFARI),]
rownames(mouseSFARI) = mouseSFARI$Symbol
nrow(mouseSFARI) # 237 expressed in mouse oxt dataset

###########################
# list of DEG in mouse SFARI
outGene$inMouseSFARI = outGene$Symbol %in% mouseSFARI$Symbol
(t1 = with(outGene,table(inMouseSFARI,inDEG = adj.P.Val < 0.05)))
fisher.test(t1) # OR = 6.022539,  p-value = 6.845e-12
ind1 = which(mouseSFARI$adj.P.Val < 0.05)
rownames(mouseSFARI[ind1,])
 # [1] "App"     "Atp1a3"  "Bdnf"    "Clstn3"  "Cntnap2" "Crhr2"   "Dlgap1"
 # [8] "Dpp6"    "Dscam"   "Erbb4"   "Gad1"    "Grin2b"  "Grpr"    "Kirrel3"
# [15] "Mef2c"   "Nrg1"    "Pcdh19"  "Reln"    "Rims1"   "Pvalb"   "Slc6a1"
# [22] "Slc9a6"  "Srrm4"   "St8sia2" "Tcf4"    "Trpc6"

#######################################
# list of DEG in scored human SFARI list
outGene$inHumanSFARI =  toupper(outGene$Symbol) %in% humanSFARI$gene.symbol
outGene_hs = outGene[grep("^ENSG", outGene$hsapien_homolog),]
(t2 = with(outGene_hs,table(inHumanSFARI,inDEG = adj.P.Val < 0.05)))
fisher.test(t2) # OR = 2.784823, p-value = 3.057e-11
ind2 = which(humanSFARI$adj.P.Val < 0.05)
humanSFARI[ind2,"Symbol"]
 # [1] "Ace"      "Ache"     "Adarb1"   "Arhgap15" "Arhgap5"  "Abat"
 # [7] "App"      "Atp1a3"   "Bdnf"     "Clstn3"   "Cnr1"     "Cntnap2"
# [13] "Cntnap3"  "Crhr2"    "Dip2a"    "Dlgap1"   "Dlgap2"   "Dpp10"
# [19] "Dpp6"     "Dpysl2"   "Dscam"    "Erbb4"    "Fat1"     "Fam135b"
# [25] "Gabbr2"   "Gad1"     "Gria1"    "Grin2b"   "Grpr"     "Il1rapl2"
# [31] "Kcnc1"    "Kirrel3"  "Lamb1"    "Mapk3"    "Lpl"      "Mef2c"
# [37] "Mcc"      "Nos1ap"   "Ntrk2"    "Ntrk3"    "Nxph1"    "Nrg1"
# [43] "Pcdh19"   "Pcdh8"    "Plxna4"   "Ptprt"    "Reln"     "Rims1"
# [49] "Robo1"    "Pcdhac1"  "Pcdhac2"  "Pvalb"    "Rit2"     "Sez6l2"
# [55] "Sik1"     "Slc6a1"   "Slc9a6"   "Smarca2"  "Snap25"   "Snx14"
# [61] "Sparcl1"  "Srrm4"    "St8sia2"  "Syt1"     "Tcf4"     "Trpc6"

#######################
### other DX ##########
#######################

ensembl = useMart("ENSEMBL_MART_ENSEMBL",  
	dataset="hsapiens_gene_ensembl", host="jul2016.archive.ensembl.org")
sym = getBM(attributes = c("ensembl_gene_id","hgnc_symbol","entrezgene"), 
		values=outGene$hsapien_homolog, mart=ensembl)
outGene$hsapien_EntrezID = sym$entrezgene[match(outGene$hsapien_homolog, sym$ensembl_gene_id)]

dxSets = read.delim("../gene_sets/Harmonizome_CTD Gene-Disease Associations Dataset.txt",
	as.is=TRUE,skip=1)
dxSets$isExpMouse = dxSets$GeneID %in% outGene$hsapien_EntrezID
dxSets$CST_enrich_Mouse = dxSets$GeneID %in% outGene$hsapien_EntrezID[
	outGene$adj.P.Val <0.05]
	
## split
dxSetsList = split(dxSets, dxSets$Disease)
dxSetsList = dxSetsList[sapply(dxSetsList, function(x) sum(x$CST_enrich_Mouse) > 0)]
length(dxSetsList)
univ = outGene[order(outGene$P.Value),]
univ = univ[!is.na(univ$hsapien_EntrezID),]
univ = univ[!duplicated(univ$hsapien_EntrezID),]
univ$CST_Geno = univ$adj.P.Val < 0.05 

dxStats = do.call("rbind", mclapply(dxSetsList, function(x) {
	x = x[x$isExpMouse,]
	inSet = univ$hsapien_EntrezID %in% x$GeneID
	tt = table(inSet, univ$CST_Geno)
	data.frame(OR = getOR(tt), p.value = chisq.test(tt)$p.value)
},mc.cores=4))

dxStats$adj.P.Val = p.adjust(dxStats$p.value)
dxStats = dxStats[order(dxStats$p.value),]
dxStats$setSize = sapply(dxSetsList,nrow)[rownames(dxStats)]
dxStats$numSig= sapply(dxSetsList,function(x) sum(x$CST_enrich_Mouse))[rownames(dxStats)]
dxStats$ID = dxSets$Mesh.or.Omim.ID[match(rownames(dxStats), dxSets$Disease)]

dxStats$sigGenes = sapply(dxSetsList,
	function(x) paste0(x$GeneSym[x$CST_enrich_Mouse], collapse=";"))[rownames(dxStats)]

write.csv(dxStats, "tables/Harmonizome_CST_IP-Genotype_effects.csv")
