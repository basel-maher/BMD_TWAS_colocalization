library("DESeq2")

ids=read.csv("data/RNA_seq/psomagen_key.csv")

countData = as.matrix(read.csv("data/RNA_seq/gene_count_matrix.csv", row.names="gene_id"))

#convert Psomagen IDs to mouse IDs
newColNames = gsub(pattern="Farber_",replacement = "",x=colnames(countData))
newColNames = gsub(pattern="_assembled.gtf",replacement = "",x=newColNames)

colnames(countData) = newColNames
countData=countData[,c(11,22,26:32,1:10,12:21,23:25)]
#
counts_femur = countData[,(1:16)]
counts_spine = countData[,(17:32)]

colData_femur = ids[c(1:16),]
colData_spine = ids[c(17:32),]
#



counts_femur = countData[,which(colnames(countData) %in% c(1:16))]
#counts_femur=counts_femur[,c(8:16,1:7)]
counts_spine = countData[,which(colnames(countData) %in% c(17:32))]

colnames(counts_femur) = ids$Mouse[match(colnames(counts_femur), ids$Psomagen.ID)]
colnames(counts_spine) = ids$Mouse[match(colnames(counts_spine), ids$Psomagen.ID)]


colData_femur = ids[c(1:16),]
rownames(colData_femur) = colData_femur$Mouse

colData_spine = ids[c(17:32),]
rownames(colData_spine) = colData_spine$Mouse

all(colData_femur$Mouse %in% colnames(counts_femur))
counts_femur <- counts_femur[, rownames(colData_femur)]
all(colData_femur$Mouse == colnames(counts_femur))

all(colData_spine$Mouse %in% colnames(counts_spine))
counts_spine <- counts_spine[, rownames(colData_spine)]
all(colData_spine$Mouse == colnames(counts_spine))

##SAMPLE MIXUP FOR FEMURS. EVEN XIST NOT DIFFERENTIALLY EXPRESSED BTWN SEXES

#femur
dds <- DESeqDataSetFromMatrix(countData = counts_femur,
                                colData = colData_femur, design = ~ Sex + Genotype)

dds <- DESeq(dds)
res <- results(dds)
(resOrdered <- res[order(res$padj), ])

x = (as.data.frame(res[order(res$padj), ]))
grep(pattern = "lrp5",ignore.case = T,x = rownames(x))


#spine
dds <- DESeqDataSetFromMatrix(countData = counts_spine,
                              colData = colData_spine, design = ~ Sex + Genotype)

dds <- DESeq(dds)
res <- results(dds)
(resOrdered <- res[order(res$padj), ])

x = (as.data.frame(res[order(res$padj), ]))
grep(pattern = "lrp5",ignore.case = T,x = rownames(x))

