library(biomaRt)


lsbmd = read.table("data/estrada/GEFOS2_LSBMD_POOLED_GC.txt", header = T)
fnbmd = read.table("data/estrada/GEFOS2_FNBMD_POOLED_GC.txt", header = T)

chr = c(1:22, "X")

ensembl = useMart("ENSEMBL_MART_SNP",dataset="hsapiens_snp") #uses human ensembl annotations

##lsbmd
lsbmd_pos <- unname(getBM(attributes=c('refsnp_id', 'chr_name', 'chrom_start'),
                          filters = 'snp_filter', values = lsbmd$MarkerName, mart = ensembl))

lsbmd_pos_clean = lsbmd_pos[which(lsbmd_pos[,2] %in% chr), ]

colnames(lsbmd_pos_clean) = c("rsid","chr","pos38")

x = merge(lsbmd, lsbmd_pos_clean, by.x = "MarkerName", by.y= "rsid")

lsbmd_38 = x

write.table(lsbmd_38, "results/lsbmd_38.txt", row.names = F, quote = F)

##fnbmd
fnbmd_pos <- unname(getBM(attributes=c('refsnp_id', 'chr_name', 'chrom_start'),
                          filters = 'snp_filter', values = fnbmd$MarkerName, mart = ensembl))

fnbmd_pos_clean = fnbmd_pos[which(fnbmd_pos[,2] %in% chr), ]

colnames(fnbmd_pos_clean) = c("rsid","chr","pos38")
x = merge(fnbmd, fnbmd_pos_clean, by.x = "MarkerName", by.y= "rsid")

fnbmd_38 = x

write.table(fnbmd_38, "results/fnbmd_38.txt", row.names = F, quote = F)
