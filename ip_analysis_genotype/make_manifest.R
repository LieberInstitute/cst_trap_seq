##
library(jaffelab)
library(readxl)

fq= list.files("/dcl01/lieber/ajaffe/lab/cst_trap_seq/preprocessed_data/FASTQ",
	pattern = "R1", full.names=TRUE, recur=TRUE)
fq = fq[!grepl("Undetermined", fq)]
names(fq) = ss(ss(fq, "/", 10), "_")

## combine
manifest = data.frame(read = fq, md5=0, SampleID = ss(names(fq), "_"))
manifest = manifest[order(manifest$SampleID),]

write.table(manifest, "preprocessed_data/samples.manifest", sep="\t",
	row.names=FALSE, col.names=FALSE, quote=FALSE)


