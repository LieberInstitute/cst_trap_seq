##

library(jaffelab)


theSamples = read.delim("/dcl01/lieber/ajaffe/Nina/Keri/rids_H7V3FBBXX.txt",
	as.is=TRUE,header=FALSE)$V1

Dir = "/dcl01/lieber/ajaffe/Nina/Keri/data"
reads = list.files(Dir, pattern="fastq.gz$", recur=TRUE, full.names=TRUE)
reads = reads[grep("Keri[0-9]", reads)]

sampleID = ss(ss(reads, "/", 8),"_")
bySample = split(reads,sampleID)
bySample = t(sapply(bySample, function(x) {
	left = paste(x[grep("_R1_",x)],collapse=",")
	right = paste(x[grep("_R2_",x)],collapse=",")
	c(left,right)}))

# new sample IDs	
rownames(bySample) = gsub("Keri", "CST_Mouse", rownames(bySample))

## reorder
bySample = bySample[order(as.numeric(ss(rownames(bySample), "e", 2))),]

## make manifest
man = data.frame(leftRead = bySample[,1], leftMd5 = 0,
	rightRead = bySample[,2], rightMd5 = 0, 
	SampleID = rownames(bySample), stringsAsFactors=FALSE)

## write out
dir.create("preprocessed_data")
write.table(man, file = "preprocessed_data/samples.manifest", 
	sep = "\t", quote=FALSE, row.names=FALSE,col.names=FALSE)