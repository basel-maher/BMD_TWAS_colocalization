#prepare Morris eBMD GWAS data for ENLOC
options(scipen=999,stringsAsFactors = FALSE,digits = 15)


BMD = read.table("data/Morrisetal2018.NatGen.SumStats/Biobank2-British-Bmd-As-C-Gwas-SumStats.txt", header=T,stringsAsFactors = F) #read all gwas snps


BMD_sub = BMD[,c(1:6,9,10,13)] #take relevant cols

BMD_sub$zscore = BMD_sub$BETA/BMD_sub$SE #calculate z-score


BMD_sub$CHR = as.numeric(BMD_sub$CHR) #numeric
BMD_sub$BP = as.numeric(BMD_sub$BP) #numeric

BMD_sub[which(BMD_sub$CHR == 23),"CHR"] = "X"

BMD_sub$chr_str = paste0("chr",BMD_sub$CHR) #make a chrN column

BMD_sub$ID = paste0(BMD_sub$chr_str,"_",BMD_sub$BP,"_",BMD_sub$EA,"_",BMD_sub$NEA,"_","b37") #GTEx format, hg19




#Used for UCSC liftOver, accessed online (sept 27 2020)
# Minimum ratio of bases that must remap: 1
##################
#lift = BMD_sub[,c(3,4)]
#lift$end = lift$BP
#lift$CHR = paste0("chr",lift$CHR)

#lift$format = paste0(lift$CHR, ":",lift$BP,"-",lift$BP)

#write.table(lift$format, file="~/Documents/projects/GWAS_project/results/lift_input_BMD", quote = F, row.names = F, col.names = F)

#this was the output. lifted contains lifted over coords, fail are snps that couldnt be lifted over. ?deleted in hg38?
lifted = read.table("~/Documents/projects/GWAS_project/results/liftOver_pass_BMD.bed",stringsAsFactors = F)
fail = read.table("~/Documents/projects/GWAS_project/results/liftOver_fails_BMD.txt",stringsAsFactors = F)
################
#

BMD_sub$format = paste0(BMD_sub$chr_str, ":",BMD_sub$BP,"-",BMD_sub$BP)

BMD_sub = BMD_sub[-which(BMD_sub$format %in% fail[,1]),]

BMD_sub$new = lifted$V1



BMD_sub$chr_hg38 = sapply(strsplit(BMD_sub$new,":"),"[",1)
BMD_sub$chr2_hg38 = sapply(strsplit(BMD_sub$chr_hg38,"chr"),"[",2)

#BMD_sub = BMD_sub[-which(BMD_sub$chr2_hg38 %in% c(1:22) == FALSE),] #remove chrX because its not in the LD annotation file below.

BMD_sub$pos_hg38 = sapply(strsplit(BMD_sub$new,"-"),"[",2)
BMD_sub$ID_new = paste0(BMD_sub$chr_hg38,"_",BMD_sub$pos_hg38,"_",BMD_sub$EA,"_",BMD_sub$NEA,"_","b38")

write.table(BMD_sub, file = "~/Documents/projects/GWAS_project/results/BMD_hg38",quote = F,row.names = F,col.names = T, sep="\t")

BMD_sub = BMD_sub[-which(BMD_sub$chr2_hg38 %in% c(1:22) == FALSE),] #remove chrX because its not in the LD annotation file below.

BMD_sub = BMD_sub[order(as.numeric(BMD_sub$chr2_hg38), as.numeric(BMD_sub$pos_hg38 )),] #order by chrom and BP

####

loc = read.table("~/Documents/projects/GWAS_project/data/eur_ld.hg38.bed", stringsAsFactors = F, header = T)
loc$loc = paste0("Loc",rownames(loc))


BMD_sub$pos_hg38 = as.numeric(BMD_sub$pos_hg38)

BMD_sub$loc = NA

val = nrow(BMD_sub)


output = c(rep(NA,val))

for(i in 1:nrow(loc)){
  x = which(BMD_sub$chr_hg38 == loc$chr[i] & BMD_sub$pos_hg38 >= loc$start[i] & BMD_sub$pos_hg38 <= loc$stop[i])
  output[x] = loc$loc[i]
  print(i)
}

BMD_sub$loc = output

BMD_sub = BMD_sub[-which(is.na(BMD_sub$loc)),]


write.table(BMD_sub[,c("ID_new","loc","zscore")], file = "~/Documents/projects/GWAS_project/results/BMD_zscore",quote = F,row.names = F,col.names = F, sep="\t")

#######################################################################################
###########################DO THE SAME BUT FOR FRACTURE GWAS###########################
#######################################################################################
options(scipen=999,stringsAsFactors = FALSE,digits = 15)

frax = as.data.frame(data.table::fread("data/Morrisetal2018.NatGen.SumStats/Biobank2-British-FracA-As-C-Gwas-SumStats.txt", header=T,stringsAsFactors = F)) #read all gwas snps


frax_sub = frax[,c(1:6,9,10,15)] #take relevant cols

frax_sub$zscore = frax_sub$logOR/frax_sub$logOR.SE #calculate z-score


frax_sub$CHR = as.numeric(frax_sub$CHR) #numeric
frax_sub$BP = as.numeric(frax_sub$BP) #numeric

frax_sub[which(frax_sub$CHR == 23),"CHR"] = "X"

frax_sub$chr_str = paste0("chr",frax_sub$CHR) #make a chrN column

frax_sub$ID = paste0(frax_sub$chr_str,"_",frax_sub$BP,"_",frax_sub$ALLELE1,"_",frax_sub$ALLELE0,"_","b37") #GTEx format, hg19




#Used for UCSC liftOver, accessed online (sept 27 2020)
# Minimum ratio of bases that must remap: 1
##################
#lift = frax_sub[,c(3,4)]
#lift$end = lift$BP
#lift$CHR = paste0("chr",lift$CHR)

#lift$format = paste0(lift$CHR, ":",lift$BP,"-",lift$BP)

#write.table(lift$format, file="~/Documents/projects/GWAS_project/results/lift_input_frax", quote = F, row.names = F, col.names = F)

#this was the output. lifted contains lifted over coords, fail are snps that couldnt be lifted over. ?deleted in hg38?
lifted = read.table("~/Documents/projects/GWAS_project/results/liftOver_pass_frax.bed",stringsAsFactors = F)
fail = read.table("~/Documents/projects/GWAS_project/results/liftOver_fails_frax.txt",stringsAsFactors = F)
################
#

frax_sub$format = paste0(frax_sub$chr_str, ":",frax_sub$BP,"-",frax_sub$BP)

frax_sub = frax_sub[-which(frax_sub$format %in% fail[,1]),]

frax_sub$new = lifted$V1



frax_sub$chr_hg38 = sapply(strsplit(frax_sub$new,":"),"[",1)
frax_sub$chr2_hg38 = sapply(strsplit(frax_sub$chr_hg38,"chr"),"[",2)


frax_sub$pos_hg38 = sapply(strsplit(frax_sub$new,"-"),"[",2)
frax_sub$ID_new = paste0(frax_sub$chr_hg38,"_",frax_sub$pos_hg38,"_",frax_sub$ALLELE1,"_",frax_sub$ALLELE0,"_","b38")


write.table(frax_sub, file = "~/Documents/projects/GWAS_project/results/frax_hg38",quote = F,row.names = F,col.names = T, sep="\t")
####

frax_sub = frax_sub[-which(frax_sub$chr2_hg38 %in% c(1:22) == FALSE),] #remove chrX because its not in the LD annotation file below.
frax_sub = frax_sub[order(as.numeric(frax_sub$chr2_hg38), as.numeric(frax_sub$pos_hg38 )),] #order by chrom and BP

loc = read.table("~/Documents/projects/GWAS_project/data/eur_ld.hg38.bed", stringsAsFactors = F, header = T)
loc$loc = paste0("Loc",rownames(loc))


frax_sub$pos_hg38 = as.numeric(frax_sub$pos_hg38)

frax_sub$loc = NA

val = nrow(frax_sub)


output = c(rep(NA,val))

for(i in 1:nrow(loc)){
  x = which(frax_sub$chr_hg38 == loc$chr[i] & frax_sub$pos_hg38 >= loc$start[i] & frax_sub$pos_hg38 <= loc$stop[i])
  output[x] = loc$loc[i]
  print(i)
}

frax_sub$loc = output

frax_sub = frax_sub[-which(is.na(frax_sub$loc)),]


write.table(frax_sub[,c("ID_new","loc","zscore")], file = "~/Documents/projects/GWAS_project/results/frax_zscore",quote = F,row.names = F,col.names = F, sep="\t")


