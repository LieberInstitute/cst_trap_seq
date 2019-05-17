## Check Oxt dataset DEGs with human asd and animal model genes from SFARI

## via badoi
lookfor= function(this,inThat) { #gene name matching
  this = toupper(this); inThat = toupper(inThat);
  tmp = sapply(this,function(x) grep(paste0('^',x,'$'),inThat))
  return(sapply(tmp,function(x) ifelse(length(x)==0,NA,x[1])))}

library(jaffelab)
library(biomaRt)
library(WriteXLS)

options(stringsAsFactors = FALSE)

##################################
# get Ensembl mouse to human genes
ensembl = useMart("ensembl",dataset = 'mmusculus_gene_ensembl')
MMtoHG = getBM(attributes = c('ensembl_gene_id','hsapiens_homolog_ensembl_gene'),
               mart = ensembl)


##############################
# load the Oxt results ####
load('tables/mutantVsWT_statistics_all.rda')

#############################

outGene$hsapien_homolog = MMtoHG$hsapiens_homolog_ensembl_gene[
	match(outGene$ensemblID,MMtoHG$ensembl_gene_id)]
outGene$sigColor = NULL 
sigGene = outGene[outGene$adj.P.Val < 0.05,]

########################
# load SFARI human genes
humanSFARI = read.csv('SFARI-Gene_genes_export03-11-2017.csv')
humanSFARI = with(humanSFARI, humanSFARI[order(-grepl('S',gene.score),ss(gene.score,'S')),])
humanSFARI = cbind(humanSFARI,outGene[lookfor(humanSFARI$gene.symbol,outGene$Symbol),])
humanSFARI= humanSFARI[!is.na(humanSFARI$Symbol),]
humanSFARI = humanSFARI[!duplicated(humanSFARI),]
nrow(humanSFARI) # 813 expressed in mouse dataset

#########################
# load SFARI mouse models
mouseSFARI = read.csv('SFARI-Gene_animal-genes_export03-11-2017.csv')
mouseSFARI = with(mouseSFARI,mouseSFARI[model.species=='Mus musculus',])
mouseSFARI = cbind(mouseSFARI,outGene[lookfor(mouseSFARI$gene.symbol,outGene$Symbol),])
mouseSFARI= mouseSFARI[!is.na(mouseSFARI$Symbol),]
mouseSFARI = mouseSFARI[order(match(mouseSFARI$Symbol,humanSFARI$Symbol)),]
mouseSFARI = mouseSFARI[!duplicated(mouseSFARI),]
nrow(mouseSFARI) # 218 expressed in mouse dataset

###########################
# list of DEG in mouse SFARI
# 24 genes in Oxt expressed and in SFARI 
outGene$inMouseSFARI = outGene$Symbol %in% mouseSFARI$Symbol
(t1 = with(outGene,table(inMouseSFARI,inDEG = adj.P.Val < 0.05)))
fisher.test(t1) # OR = 5.73, p-value = 2.52e-10
ind1 = which(mouseSFARI$Symbol %in% mouseSFARI$Symbol & mouseSFARI$adj.P.Val < 0.05)

#######################################
# list of DEG in scored human SFARI list
# 131 human SFARI ASD-linked genes in our list of DEGs
outGene$inHumanSFARI =  outGene$Symbol %in% humanSFARI$Symbol
outGene_hs = outGene[grep("^ENSG", outGene$hsapien_homolog),]
(t2 = with(outGene_hs,table(inHumanSFARI,inDEG = adj.P.Val < 0.05)))
fisher.test(t2) # OR = 2.80, p-value <1.75e-10
ind2 = which(humanSFARI$Symbol %in% humanSFARI$Symbol & humanSFARI$adj.P.Val < 0.05)
	
## venn diagrams
library(limma)
vennCount = vennCounts(data.frame(OxtEnr = outGene$adj.P.Val < 0.05,
	`Hs SFARI` = outGene$inHumanSFARI,
	`Mm SFARI`= outGene$inMouseSFARI))
vennCounts_hs = with(outGene[grep("^ENSG", outGene$hsapien_homolog),],
	vennCounts(data.frame(OxtEnr = adj.P.Val < 0.05,
	`Hs SFARI` = inHumanSFARI,
	`Mm SFARI`= inMouseSFARI)))

pdf("plots/vennDiagram_asd_enrichment.pdf")
vennDiagram(vennCount,cex=1.4,main="Mouse Expressed")
vennDiagram(vennCounts_hs,cex=1.4, main = "Mouse Exprs (Hom)")
dev.off()


####################
# save the two lists
sfariList = list(Scored_Human_SFARI_Genes = humanSFARI[ind2,],
                 Mouse_Model_Genes = mouseSFARI[ind1,])
WriteXLS(sfariList,ExcelFileName = 'tables/SFARI_asd_gene_list.xls')

