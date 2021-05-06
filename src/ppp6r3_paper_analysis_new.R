library(ggplot2)
#library(patchwork)
library(biomaRt)
#library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(topGO)
#library(IMPCdata)
#library(PhenStat)
library(RACER)
library(ggrepel)
library(LDlinkR)
library(GenomicRanges)
library(qqman)

set.seed(8675309)

setwd("../../../../scratch/bma8ne/GTEx_v8/src")


# read_in -----------------------------------------------------------------

#read in hg38 version of eBMD GWAS

#bmd = as.data.frame(data.table::fread("../BMD_hg38", header=T))

#saveRDS(bmd, "assets/bmd.rds")


#read in FastENLOC sig files for BMD

enloc_sig = read.table("../enloc_out/all_enloc_bmd.sig.out")
colnames(enloc_sig) = c("cluster_name","nSNPs", "cluster_PIP_eqtl", "cluster_PIP_gwas_noeqtl","cluster_PIP_gwas_w_eqtl", "RCP","tissue")
enloc_sig$tissue = sapply(strsplit(enloc_sig$tissue,"[.]"),"[",1)

#saveRDS(enloc_sig, "assets/enloc_sig.rds")

#split ensembl name
enloc_sig$ens = sapply(strsplit(enloc_sig$cluster_name, ":"),"[",1)

#get annotation, remove everything but coding genes
genes = unique(enloc_sig$ens)
mart = useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
ovr = getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol","gene_biotype"),values=genes,mart= mart)

keep = ovr[which(ovr$gene_biotype == "protein_coding"),"ensembl_gene_id"]
enloc_sig_prot = enloc_sig[which(enloc_sig$ens %in% keep),]

#saveRDS(enloc_sig_prot, "assets/enloc_sig_prot.rds")

#RCP >=0.1
enloc_sig_prot_0.1 = enloc_sig_prot[which(enloc_sig_prot[,"RCP"]>=0.1),]
enloc_sig_prot_0.5 = enloc_sig_prot[which(enloc_sig_prot[,"RCP"]>=0.5),]

#saveRDS(enloc_sig_prot_0.1, "assets/enloc_sig_prot_0.1.rds")
#saveRDS(enloc_sig_prot_0.5, "assets/enloc_sig_prot_0.5.rds")


#read in multixcan bmd TWAS

multixcan = read.delim("../bmd_morris_12_1_20_multixcan.txt")

multixcan$gene = sapply(strsplit(multixcan$gene, "[.]"),"[",1) #strsplit to get gene name


length(unique(multixcan$gene)) #22337
length(unique(multixcan$gene_name)) #22331

#LYNX1 and MAL2 each have two different ENS, for mal2 one is mal2-as1, for lynx1, one is SLURP2
multixcan[which(multixcan$gene == "ENSG00000283992"),"gene_name"] = "SLURP2"
multixcan[which(multixcan$gene == "ENSG00000253972"),"gene_name"] = "MAL2-AS1" #lnc, check is removed



#get annotation, remove everything but coding genes
genes = unique(multixcan$gene)
mart = useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
ovr = getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol","gene_biotype"),values=genes,mart= mart)

#FDR correction
multixcan$fdr = p.adjust(multixcan$pvalue,method = "fdr")

#Bonferroni correction
multixcan$bonf = p.adjust(multixcan$pvalue,method = "bonf")

#saveRDS(multixcan, "assets/multixcan.rds")

#remove non-protein coding
keep = ovr[which(ovr$gene_biotype == "protein_coding"),"ensembl_gene_id"]
multixcan_prot = multixcan[which(multixcan$gene %in% keep),]
#saveRDS(multixcan_prot, "assets/multixcan_prot.rds")

length(unique(multixcan_prot$gene)) #17310
length(unique(multixcan_prot$gene_name)) #17310

multixcan_prot_fdr_sig = multixcan_prot[which(multixcan_prot$fdr <=0.05),]
multixcan_prot_bonf_sig = multixcan_prot[which(multixcan_prot$bonf <=0.05),]

#saveRDS(multixcan_prot_fdr_sig, "assets/multixcan_prot_fdr_sig.rds")
#saveRDS(multixcan_prot_bonf_sig, "assets/multixcan_prot_bonf_sig.rds")

length(unique(multixcan_prot_fdr_sig$gene_name)) #6483
length(unique(multixcan_prot_bonf_sig$gene_name)) #2156

### enloc snps for ppp6r3 in thyroid ###

#thyroid eQTL
thyroid = as.data.frame(data.table::fread("../GTEx_Analysis_v8_eQTL_all_associations/Thyroid.allpairs.txt",stringsAsFactors = FALSE,header = TRUE, sep="\t"))#readtissue file
thyroid = thyroid[grep(pattern = "chr11_",x = thyroid$variant_id),]

#get RSIDs from GTEx lookup tables
chr11_lookup = read.table("../GTEx_v8_lookup/chr11_GTEx_V8_lookup.txt", header=T)

thyroid = merge(thyroid, chr11_lookup, by.x = "variant_id", by.y = "variant_id")

ppp6r3_thyroid_eqtl = thyroid[grep("ENSG00000110075",thyroid$gene_id),]


saveRDS(ppp6r3_thyroid_eqtl, "assets/ppp6r3_thyroid_eqtl.rds")
saveRDS(thyroid, "assets/all_thyroid_chr11_gtex_eqtl.rds")


#### read from RDS ####
bmd = readRDS("assets/bmd.rds") #hg38 eBMD GWAS
enloc_sig = readRDS("assets/enloc_sig.rds") #fastENLOC sig files
enloc_sig_prot = readRDS("assets/enloc_sig_prot.rds") #fastENLOC sig files that are protein coding
enloc_sig_prot_0.1 = readRDS("assets/enloc_sig_prot_0.1.rds") #fastENLOC sig files that are protein coding, RCP >= 0.1
enloc_sig_prot_0.5 = readRDS("assets/enloc_sig_prot_0.5.rds") #fastENLOC sig files that are protein coding, RCP >= 0.5

multixcan = readRDS("assets/multixcan.rds") #multixcan results after removal of dup ensembls
multixcan_prot = readRDS("assets/multixcan_prot.rds") #protein-coding multixcan results
multixcan_prot_fdr_sig = readRDS("assets/multixcan_prot_fdr_sig.rds") #protein coding, FDR <= 0.05
multixcan_prot_bonf_sig = readRDS("assets/multixcan_prot_bonf_sig.rds") #protein coding, Bonferroni <= 0.05

t.all = readRDS("assets/t.all") # topGO gene ontology for the 512 prioritized genes

ppp6r3_thyroid_eqtl = readRDS("assets/ppp6r3_thyroid_eqtl.rds")

length(unique(enloc_sig$ens)) #38534
length(unique(enloc_sig_prot$ens)) #17884
length(unique(enloc_sig_prot_0.1$ens)) #1182
length(unique(enloc_sig_prot_0.5$ens)) #278
length(unique(multixcan_prot$gene)) #17310
length(unique(multixcan_prot_fdr_sig$gene_name)) #6483
length(unique(multixcan_prot_bonf_sig$gene_name)) #2156




#### prioritized genes (merge TWAS and ENLOC), superset pruning and enrichment of prioritized genes in bone genes ####

#merge enloc and TWAS
prioritized_bonf_0.1 = merge(enloc_sig_prot_0.1,multixcan_prot_bonf_sig, by.x="ens", by.y="gene")

length(unique(prioritized_bonf_0.1$ens)) #512

#superset
superset = read.table("../superset_humanized_proteinCoding.txt", stringsAsFactors = F)
colnames(superset) = c("SYMBOL", "ENSEMBL")

#number prioritized genes in bone superset
length(unique(prioritized_bonf_0.1[which(prioritized_bonf_0.1$ens %in% superset$ENSEMBL),"ens"])) #66


#phyper, enrichment of prioritized genes in bone genes
#total protein coding genes
total_prot = unique(c(multixcan_prot$gene, enloc_sig_prot$ens))

length(unique(superset[which(superset$ENSEMBL %in% total_prot), "ENSEMBL"])) #1,425/1,484 bone genes are 

superset_pruned = superset[which(superset$ENSEMBL %in% total_prot),]


# m = length(which((superset_pruned$ENSEMBL) %in% total_prot)) #num bone genes in total 
# n = length(total_prot) - m #number non-bone genes in total
# k = length(unique(prioritized_bonf_0.1$gene_name)) #length prioritized genes, thresh 0.1
# q = length(unique(prioritized_bonf_0.1[which(prioritized_bonf_0.1$ens %in% superset_pruned$ENSEMBL),"ens"])) #number prioritized genes in bone superset
# phyper(q=q-1,m=m,n=n,k=k,lower.tail = F) #6.768268e-05 if using bonf

a = 66 #prioritized genes in known bone list
b = length(unique(prioritized_bonf_0.1$ens)) - a #priortized genes not in known bone list
c = length(unique(total_prot[which( (total_prot %in% prioritized_bonf_0.1$ens == F) & (total_prot %in% superset_pruned$ENSEMBL == TRUE) )]))#not prioritized and IN bone
d = length(unique(total_prot[which( (total_prot %in% prioritized_bonf_0.1$ens == F) & (total_prot %in% superset_pruned$ENSEMBL == F) )]))#not prioritized and not bone

mat = matrix(nrow=2, ncol=2)
mat[1,] = c(a,c)
mat[2,] = c(b,d)

fisher.test(mat,alternative = "g") #OR=1.74
fisher.test(mat,alternative = "g")$p.value #pval = 6.768268e-05

#### prioritized gene localization ####

#where are prioritized genes located, locus wise? how are they distributed across the genome?

#get gwas lead snps location +/- 1 mbp
bmd_lead_hg19 = read.csv("../BMD_morris_lead_snps_hg19.csv", header=F) #hg19

#convert to hg38
ensembl <- useEnsembl("snp",dataset = "hsapiens_snp")

#get genomic position
bmd_lead_hg38 <- getBM(attributes=c("refsnp_id",
                           "chr_name",
                           "chrom_start",
                           "chrom_end"),
              filters ="snp_filter", values =bmd_lead_hg19$V2, mart = ensembl, uniqueRows=TRUE)

table(bmd_lead_hg38$chr_name)
bmd_lead_hg38 = bmd_lead_hg38[-which(bmd_lead_hg38$chr_name %in% c(1:22, "X") == FALSE),] #remove chromosome patches

#add in the one SNP without an RSID, 7:120812727_G_C
#using ucsc liftover: #chr7:121172673-121172673 , overlaps rs6150319 but this a del/ins that is 42 bases long
bmd_lead_hg38 = rbind(bmd_lead_hg38, c("7:120812727_G_C", "7", "121172673", "121172673"))

length(unique(bmd_lead_hg38$refsnp_id))

#get prioritized gene start and end from ensembl
prioritized_unique = unique(prioritized_bonf_0.1$ens)
ensembl=useMart("ensembl", dataset = "hsapiens_gene_ensembl")



prioritized_gene_pos <- getBM(attributes=c('ensembl_gene_id',
                                           'chromosome_name',
                                            'start_position',
                                            'end_position'),
                       filters ="ensembl_gene_id", values =prioritized_unique, mart = ensembl, uniqueRows=TRUE)

#why no X chromosome? multixcan did not include sex chromosomes
prioritized_bonf_0.1 = merge(prioritized_bonf_0.1, prioritized_gene_pos, by.x = "ens",by.y="ensembl_gene_id") #merge in positions

#overlap, count number of overlaps per lead snp
#first convert snps to loci, where a locus is a snps pos +/- 1 Mbp (1000000 bp)
bmd_lead_hg38$chr_name = paste0("chr",bmd_lead_hg38$chr_name)
bmd_lead_hg38$start_locus = as.numeric(bmd_lead_hg38$chrom_start) - 1000000
bmd_lead_hg38$end_locus = as.numeric(bmd_lead_hg38$chrom_end) + 1000000

bmd_lead_granges = makeGRangesFromDataFrame(bmd_lead_hg38, seqnames.field = "chr_name", start.field = "start_locus",end.field = "end_locus", keep.extra.columns = T)
#for above, make sure 0 based vs 1 based

#now convert genes granges

prioritized_gene_pos$chr_name = paste0("chr",prioritized_gene_pos$chromosome_name)

prioritized_gene_granges = makeGRangesFromDataFrame(prioritized_gene_pos, seqnames.field = "chr_name", start.field = "start_position",end.field = "end_position", keep.extra.columns = T)

overlaps = findOverlaps(bmd_lead_granges, prioritized_gene_granges)
overlaps = as.data.frame(overlaps)
# 456 genes overlap 542 snps (919 total overlaps, so a gene can overlap multiple snps)

x = bmd_lead_hg38[overlaps$queryHits,]
y = prioritized_gene_pos[overlaps$subjectHits,]
z = cbind(x,y)

#match gene names from prioritized_bonf to z via ensembl
x = prioritized_bonf_0.1[,c("ens","gene_name")]
x = unique(x)

z = merge(x, z, by.x="ens",by.y="ensembl_gene_id")
##### gene ontology for prioritized (merged) genes ####

allGenes = unique(c(multixcan_prot$gene, enloc_sig_prot$ens))

interesting.genes = unique(prioritized_bonf_0.1$ens)

geneList<-factor(as.integer(allGenes %in% interesting.genes)) #If TRUE returns 1 as factor, otherwise 0
names(geneList)<-allGenes
###MF###
GOdata <- new("topGOdata", ontology = "MF", allGenes =geneList,
              annot = annFUN.org, mapping='org.Hs.eg.db', ID='ENSEMBL')
test.stat<-new("classicCount", testStatistic = GOFisherTest, name='Fisher test')
result<-getSigGroups(GOdata,test.stat)
t1<-GenTable(GOdata, classic=result, topNodes=length(result@score))
head(t1)
###CC###
GOdata <- new("topGOdata", ontology = "CC", allGenes = geneList,
              annot = annFUN.org, mapping='org.Hs.eg.db', ID='ENSEMBL')
test.stat<-new("classicCount", testStatistic = GOFisherTest, name='Fisher test')
result<-getSigGroups(GOdata,test.stat)
t2<-GenTable(GOdata, classic=result, topNodes=length(result@score))
head(t2)
###BP###
GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,
              annot = annFUN.org, mapping='org.Hs.eg.db', ID='ENSEMBL')
test.stat<-new("classicCount", testStatistic = GOFisherTest, name='Fisher test')
result<-getSigGroups(GOdata,test.stat)
t3<-GenTable(GOdata, classic=result, topNodes=length(result@score))
head(t3)
####
t.all = NULL
t.all<-rbind(t1,t2,t3)
t.all$classic<-as.numeric(as.character(t.all$classic))

#saveRDS(t.all, "assets/t.all")


#### IMPC cross-ref ####
impc = readRDS("assets/impc_out_5e2.rds") #all genes with at least one pvalue <=0.05 (male, female, genotype)
impc = impc[,"marker_symbol"]

impc = unique(impc) #2,508 genes

#convert mouse to human gene names

human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")


impc_hum = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = impc , mart = mouse, attributesL = c("hgnc_symbol", "ensembl_gene_id"), martL = human, uniqueRows=T)

impc_hum = unique(impc_hum) 
impc_hum = impc_hum[-which(impc_hum$HGNC.symbol == ""),]

length(unique(impc_hum$MGI.symbol)) #2,240 unique mouse genes in impc data have homologues in humans. some, like Sirpa, have multiple human homologs.
length(unique(impc_hum$Gene.stable.ID)) #2,322 total human homologues

length(unique(prioritized_bonf_0.1[which(prioritized_bonf_0.1$ens %in% impc_hum$Gene.stable.ID),"ens"])) #64

x = (unique(prioritized_bonf_0.1[which(prioritized_bonf_0.1$ens %in% impc_hum$Gene.stable.ID),"ens"]))
length(which(x %in% superset_pruned$ENSEMBL == FALSE)) # 48 of these are unique to IMPC, not in superset


a = 64 #prioritized genes in IMPC
b = length(unique(prioritized_bonf_0.1$ens)) - a #priortized genes not in IMPC
c = length(unique(total_prot[which( (total_prot %in% prioritized_bonf_0.1$ens == F) & (total_prot %in% impc_hum$Gene.stable.ID == TRUE) )]))#not prioritized and IN IMPC
d = length(unique(total_prot[which( (total_prot %in% prioritized_bonf_0.1$ens == F) & (total_prot %in% impc_hum$Gene.stable.ID == F) )]))#not prioritized and not bone

mat = matrix(nrow=2, ncol=2)
mat[1,] = c(a,c)
mat[2,] = c(b,d)

fisher.test(mat,alternative = "g") #OR=1.03
fisher.test(mat,alternative = "g")$p.value #pval = 0.3


#### PPP6R3 locus analyses ####

#prioritized genes with RCP >=0.5
prioritized_bonf_0.5 = prioritized_bonf_0.1[which(prioritized_bonf_0.1$RCP >= 0.5),]
 
#remove known bone genes
prioritized_bonf_0.5 = prioritized_bonf_0.5[-which(prioritized_bonf_0.5$ens %in% superset_pruned$ENSEMBL),]
#remove IMPC
prioritized_bonf_0.5 = prioritized_bonf_0.5[-which(prioritized_bonf_0.5$ens %in% impc_hum$Gene.stable.ID),]

length(unique(prioritized_bonf_0.5$ens))

#ppp6r3_enloc = enloc_sig_prot[which(enloc_sig_prot$ens == "ENSG00000110075"),]

#LD of rs10047483 with all other BMD lead snps in locus

#lead snps in locus
lead_rsid = bmd_lead_hg38[which(bmd_lead_hg38$chr_name == "chr11" & bmd_lead_hg38$start_locus <= (68228199) & bmd_lead_hg38$end_locus >= (68228199)),] #check locations

LDmatrix(c(lead_rsid$refsnp_id, "rs10047483"), pop = "EUR", r2d = "r2", token = "2ebcd9d51762") #LD between snp and bmd lead snps in ppp6r3 locus

# rs10047483 ref G alt A in gtex, with positive slope. effect allele is G in GWAS, with positive slope. However, in GTEx, eQTL allele is the ALT allele. Therefore, the G allele effect is negative in GTEx, and positive in GWAS

##LRP5 = "ENSG00000162337
LDmatrix(c(lead_rsid$refsnp_id, "rs3736228"), pop = "EUR", r2d = "r2", token = "2ebcd9d51762") #LD between snp and bmd lead snps in ppp6r3 locus
LDpair("rs10047483","rs3736228", pop = "EUR", token = "2ebcd9d51762")
####

#### Frax GWAS analysis ####
## Frax GWAS ##
frax = as.data.frame(data.table::fread("../frax_hg38",header=TRUE,stringsAsFactors = F))
#get RSIDs from GTEx lookup tables
chr11_lookup = read.table("../GTEx_v8_lookup/chr11_GTEx_V8_lookup.txt", header=T)

frax_merged_rsid = merge(frax, chr11_lookup, by.x = "ID_new", by.y = "variant_id")

###
#hg38 version of eBMD GWAS
frax = as.data.frame(data.table::fread("../frax_hg38", header=T))


#read in enloc sig files for BMD

enloc_frax_sig = read.table("../enloc_frax_out/all_enloc_frax.sig.out")
colnames(enloc_frax_sig) = c("cluster_name","nSNPs", "cluster_PIP_eqtl", "cluster_PIP_gwas_noeqtl","cluster_PIP_gwas_w_eqtl", "RCP","tissue")
enloc_frax_sig$tissue = sapply(strsplit(enloc_frax_sig$tissue,"[.]"),"[",1)

#split ensembl name
enloc_frax_sig$ens = sapply(strsplit(enloc_frax_sig$cluster_name, ":"),"[",1)
saveRDS
ppp6r3_frax_enloc = enloc_frax_sig[which(enloc_frax_sig$ens == "ENSG00000110075"),]


#get annotation, remove everything but coding genes
genes = unique(enloc_frax_sig$ens)
mart = useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
ovr = getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol","gene_biotype"),values=genes,mart= mart)

keep = ovr[which(ovr$gene_biotype == "protein_coding"),"ensembl_gene_id"]
enloc_frax_sig_prot = enloc_frax_sig[which(enloc_frax_sig$ens %in% keep),]

saveRDS


#RCP >=0.1
enloc_frax_sig_prot_0.1 = enloc_frax_sig_prot[which(enloc_frax_sig_prot[,"RCP"]>=0.1),]
enloc_frax_sig_prot_0.5 = enloc_frax_sig_prot[which(enloc_frax_sig_prot[,"RCP"]>=0.5),]

length(unique(enloc_frax_sig_prot_0.1$ens)) #

#multixcan bmd TWAS
multixcan_frax = read.delim("../frac_morris_8_14_20_multixcan.txt")

multixcan_frax$gene = sapply(strsplit(multixcan_frax$gene, "[.]"),"[",1)

length(unique(multixcan_frax$gene))
length(unique(multixcan_frax$gene_name))

which(duplicated(multixcan_frax$gene_name))


multixcan_frax[which(multixcan_frax$gene == "ENSG00000283992"),"gene_name"] = "SLURP2"
multixcan_frax[which(multixcan_frax$gene == "ENSG00000253972"),"gene_name"] = "MAL2-AS1" #lnc, check is removed

#get annotation, remove everything but coding genes
genes = unique(multixcan_frax$gene)
mart = useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
ovr = getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol","gene_biotype"),values=genes,mart= mart)

#FDR correction
multixcan_frax$fdr = p.adjust(multixcan_frax$pvalue,method = "fdr")
multixcan_frax$bonf = p.adjust(multixcan_frax$pvalue,method = "bonf")
#remove non-protein coding
keep = ovr[which(ovr$gene_biotype == "protein_coding"),"ensembl_gene_id"]
multixcan_frax_prot = multixcan_frax[which(multixcan_frax$gene %in% keep),]

multixcan_frax_prot_fdr_sig = multixcan_frax_prot[which(multixcan_frax_prot$fdr <=0.05),]
multixcan_frax_prot_bonf_sig = multixcan_frax_prot[which(multixcan_frax_prot$bonf <=0.05),]

length(unique(multixcan_frax_prot_fdr_sig$gene_name))
length(unique(multixcan_frax_prot_bonf_sig$gene_name))

##
##

#merge enloc and TWAS
merged_enloc_multixcan_frax_fdr = merge(enloc_frax_sig_prot_0.1,multixcan_frax_prot_fdr_sig, by.x="ens", by.y="gene")
merged_enloc_multixcan_frax_bonf = merge(enloc_frax_sig_prot_0.1,multixcan_frax_prot_bonf_sig, by.x="ens", by.y="gene")

length(unique(merged_enloc_multixcan$ens))
length(unique(merged_enloc_multixcan$bonf))

frax_ppp6r3 = frax[which(frax$chr_hg38 == "chr11" & frax$BP <= (68228199+1000000) & frax$BP >= (68228199-1000000)),] #check locations
#ppp6r3_thyroid_eqtl[which(ppp6r3_thyroid_eqtl$rs_id_dbSNP151_GRCh38p7 == "rs35989399"),]
LDpair("rs10047483","rs35989399", pop = "EUR", token = "2ebcd9d51762")


####FIGURE 1####
require(qqman)

par(mfrow=c(2,1))
par(mar=c(0,5,3,3))

twas = multixcan_prot[,c("gene","gene_name","bonf")]

ensembl=useMart("ensembl", dataset = "hsapiens_gene_ensembl")

twas.pos <- getBM(attributes=c('ensembl_gene_id','chromosome_name',
                                           'start_position'),
                              filters ="ensembl_gene_id", values =twas$gene, mart = ensembl, uniqueRows=TRUE)

twas = merge(twas, twas.pos, by.x="gene", by.y="ensembl_gene_id")

twas[which(twas$bonf == 0),"bonf"] = 1e-323 #0 is less than 1e-324

colnames(twas) = c("ens","gene","P","CHR","BP")
  
manhattan(twas,ylim=c(0,350),cex=1,cex.lab=1,font.lab=1,font.axis=1,cex.axis=1,las=2,font=4, suggestiveline = F, genomewideline = 1.3)

par(mar=c(5,5,3,3))

# colnames(bmd)[which(colnames(bmd) == "P.NI")] = "P"
# bmd = bmd[-which(bmd$CHR == "X"),]
# bmd$CHR = as.numeric(bmd$CHR)
# manhattan(bmd,ylim=c(325,0),cex=2.2,cex.lab=2.5,font.lab=2,font.axis=2,cex.axis=1.6,las=2,font=4,xlab="",xaxt="n", suggestiveline = F)

colocs = enloc_sig_prot
colocs.agg = aggregate(RCP~ens, colocs, max)

colocs.agg.pos <- getBM(attributes=c('ensembl_gene_id','chromosome_name',
                               'start_position'),
                  filters ="ensembl_gene_id", values =colocs.agg$ens, mart = ensembl, uniqueRows=TRUE)

colocs.agg.pos = merge(colocs.agg, colocs.agg.pos, by.x="ens", by.y="ensembl_gene_id")
#colocs.agg.pos$P = 10^(-colocs.agg.pos$RCP)
#colocs.agg.pos$P[which(colocs.agg.pos$RCP == 0)] = 1
colnames(colocs.agg.pos) = c("ens","P","CHR","BP")
manhattan(colocs.agg.pos,ylim=c(1.1,0),cex=1,cex.lab=1,font.lab=1,font.axis=1,cex.axis=1,las=2,font=4,xlab="",xaxt="n", suggestiveline = F, logp=F, genomewideline = 0.1)



###Do the same but for significant TWAS/ENLOCs. Maybe also for the shared genes between the two (intersection)

par(mfrow=c(2,1))

  
twas = multixcan_prot_bonf_sig[,c("gene","gene_name","bonf")]

ensembl=useMart("ensembl", dataset = "hsapiens_gene_ensembl")

twas.pos <- getBM(attributes=c('ensembl_gene_id','chromosome_name',
                               'start_position'),
                  filters ="ensembl_gene_id", values =twas$gene, mart = ensembl, uniqueRows=TRUE)

twas = merge(twas, twas.pos, by.x="gene", by.y="ensembl_gene_id")

twas[which(twas$bonf == 0),"bonf"] = 1e-323 #0 is less than 1e-324

colnames(twas) = c("ens","gene","P","CHR","BP")

manhattan(twas,ylim=c(0,350),cex=1,cex.lab=1,font.lab=1,font.axis=1,cex.axis=1,las=2,font=4, suggestiveline = F, genomewideline = 1.3)




# colnames(bmd)[which(colnames(bmd) == "P.NI")] = "P"
# bmd = bmd[-which(bmd$CHR == "X"),]
# bmd$CHR = as.numeric(bmd$CHR)
# manhattan(bmd,ylim=c(325,0),cex=2.2,cex.lab=2.5,font.lab=2,font.axis=2,cex.axis=1.6,las=2,font=4,xlab="",xaxt="n", suggestiveline = F)

colocs = enloc_sig_prot_0.1
colocs.agg = aggregate(RCP~ens, colocs, max)

colocs.agg.pos <- getBM(attributes=c('ensembl_gene_id','chromosome_name',
                                     'start_position'),
                        filters ="ensembl_gene_id", values =colocs.agg$ens, mart = ensembl, uniqueRows=TRUE)

colocs.agg.pos = merge(colocs.agg, colocs.agg.pos, by.x="ens", by.y="ensembl_gene_id")
#colocs.agg.pos$P = 10^(-colocs.agg.pos$RCP)
#colocs.agg.pos$P[which(colocs.agg.pos$RCP == 0)] = 1
colnames(colocs.agg.pos) = c("ens","P","CHR","BP")
manhattan(colocs.agg.pos,ylim=c(1.1,0),cex=1,cex.lab=1,font.lab=1,font.axis=1,cex.axis=1,las=2,font=4,xlab="",xaxt="n", suggestiveline = F, logp=F, genomewideline = 0.1)

twas_shared = twas[which(twas$ens %in% colocs.agg.pos$ens),]
colocs_shared = colocs.agg.pos[which(colocs.agg.pos$ens %in% twas$ens),]

manhattan(twas_shared,ylim=c(0,350),cex=1,cex.lab=1,font.lab=1,font.axis=1,cex.axis=1,las=2,font=4, suggestiveline = F, genomewideline = F)
manhattan(colocs_shared,ylim=c(1.1,0),cex=1,cex.lab=1,font.lab=1,font.axis=1,cex.axis=1,las=2,font=4,xlab="",xaxt="n", suggestiveline = F, logp=F, genomewideline = F)

#sort both by ens then do spearman cor
twas_shared = twas_shared[order(twas_shared$ens),]
colocs_shared = colocs_shared[order(colocs_shared$ens),]
cor(twas_shared$P, colocs_shared$P, method = "s")
#dev.off()


####Locus Fig, lead snps in ppp6r3 locus####
#get BMD lead SNPS fpr PPP6R3
lead_snps = read.csv("../BMD_morris_lead_snps_hg19.csv", header=F)
#ppp6r3 starts at chr11:68,228,199 on hg19

#lead snp rsids within +/- a megabase of the "TSS" or ppp6r3 (ENCODE beginning of gene)
lead_rsid = lead_snps[which(lead_snps$V3 == 11 & lead_snps$V4 <= (68228199+1000000) & lead_snps$V4 >= (68228199-1000000)),]


#plot
ppp6r3_gwas_min = bmd[,c("pos_hg38","P.NI","CHR","RSID")]
ppp6r3_gwas_min = ppp6r3_gwas_min[which(ppp6r3_gwas_min$CHR == 11),]#remove non chromosome 11 snps otherwise theres repeats

ppp6r3_gwas_min$P.NI = as.numeric(ppp6r3_gwas_min$P.NI)
ppp6r3_gwas_min$pos_hg38 = as.numeric(ppp6r3_gwas_min$pos_hg38)
ppp6r3_gwas_min$col = "notleadsnp"
ppp6r3_gwas_min[which(ppp6r3_gwas_min$RSID %in% lead_rsid$V2),"col"]="leadsnp"

#PLOT ALL 7 LEAD SNPS WITHIN 1MBP OF PPP6R3. TAKE GENE LEGEND FROM RACER SINGLE PLOT BELOW
ggplot(ppp6r3_gwas_min, aes(pos_hg38,-log10(P.NI),color=col, label=RSID)) + geom_point(size=1, alpha=0.5) +xlim(68460724-1000000,68460724+1000000) + ggtitle("eBMD GWAS independent SNPs within 1 Mbp of PPP6R3 start site")+
  geom_text_repel(data = subset(ppp6r3_gwas_min, col == "leadsnp"), point.padding = 0.5, box.padding = 0.7, force = 50, nudge_y = 150, direction="x", segment.size = 0.2, size=3) + theme_bw() + 
  annotate("point",ppp6r3_gwas_min$pos_hg38[which(ppp6r3_gwas_min$col=="leadsnp")],-log10(ppp6r3_gwas_min$P.NI[which(ppp6r3_gwas_min$col=="leadsnp")]),size=1)
#add ppp6r3 pos


gwas_racer = RACER::formatRACER(assoc_data = ppp6r3_gwas_min, chr_col = 3, pos_col = 1, p_col = 2)
gwas_racer = RACER::ldRACER(assoc_data = gwas_racer, rs_col = 4, pops = "EUR", lead_snp = "rs3736228")
#gwas_racer$LD = 0
#gwas_racer$LD_BIN = "0.2-0.0"
#gwas_racer$newlabel = NA
#gwas_racer[which(gwas_racer$RS_ID %in% lead_rsid),"newlabel"]= gwas_racer[which(gwas_racer$RS_ID %in% lead_rsid),"RS_ID"]

# GENES FOR ABOVE FROM HERE
RACER::singlePlotRACER(assoc_data = gwas_racer, chr = 11, build = "hg38", plotby = "coord", label_lead = F,start_plot =68460724-1000000, end_plot = 68460724+1000000)
#### mirrorplot of ppp6r3 eQTL and eBMD ####
ppp6r3_gwas_min = bmd[,c("pos_hg38","P.NI","CHR","RSID")]
ppp6r3_gwas_min = ppp6r3_gwas_min[which(ppp6r3_gwas_min$CHR == 11),]#remove non chromosome 11 snps otherwise theres repeats

ppp6r3_gwas_min$P.NI = as.numeric(ppp6r3_gwas_min$P.NI)
ppp6r3_gwas_min$pos_hg38 = as.numeric(ppp6r3_gwas_min$pos_hg38)
ppp6r3_gwas_min$col = "notleadsnp"
ppp6r3_gwas_min[which(ppp6r3_gwas_min$RSID %in% lead_rsid$V2),"col"]="leadsnp"


gwas_racer = RACER::formatRACER(assoc_data = ppp6r3_gwas_min, chr_col = 3, pos_col = 1, p_col = 2)

gwas_racer = RACER::ldRACER(assoc_data = gwas_racer, rs_col = 4, pops = "EUR", lead_snp = "rs10047483")

#
ppp6r3_thyroid_eqtl$chr = sapply(strsplit(ppp6r3_thyroid_eqtl$chr, "chr"),"[",2)
ppp6r3_thyroid_eqtl$pval_nominal = as.numeric(ppp6r3_thyroid_eqtl$pval_nominal)
#ppp6r3_thyroid_eqtl = ppp6r3_thyroid_eqtl[,c(1,7,17,18,23)]
eqtl_racer = RACER::formatRACER(assoc_data = ppp6r3_thyroid_eqtl, chr_col = 10, pos_col = 11, p_col = 7)


eqtl_racer = RACER::ldRACER(assoc_data = eqtl_racer, rs_col = 15, pops = "EUR", lead_snp = "rs10047483")

eqtl_racer = unique(eqtl_racer)

mirrorPlotRACER(assoc_data1 = gwas_racer, assoc_data2 = eqtl_racer, chr = 11, plotby = "coord",build = "hg38",label_lead = T,name1 = "eBMD SNPs", name2 = "PPP6R3 eQTL - Thyroid", start_plot =68460724-1000000, end_plot = 68460724+1000000)

scatterPlotRACER(assoc_data1 = gwas_racer, assoc_data2 = eqtl_racer, chr = 11, name1 = "BMD GWAS p_values", name2 = "PPP6R3 eQTL p_values", region_start = 68460724-1000000, region_end = 68460724+1000000, ld_df = 1)


#### mirrorplot of ppp6r3 eQTL and frax ####
frax = as.data.frame(data.table::fread("../frax_hg38",header=TRUE,stringsAsFactors = F))
#get RSIDs from GTEx lookup tables
chr11_lookup = read.table("../GTEx_v8_lookup/chr11_GTEx_V8_lookup.txt", header=T)

frax_merged_rsid = merge(frax, chr11_lookup, by.x = "ID_new", by.y = "variant_id")

ppp6r3_frax_min = frax_merged_rsid[,c("pos_hg38","P.I","chr2_hg38","rs_id_dbSNP151_GRCh38p7")]
ppp6r3_frax_min = ppp6r3_frax_min[which(ppp6r3_frax_min$chr2_hg38 == 11),]#remove non chromosome 11 snps otherwise theres repeats

ppp6r3_frax_min$P.I = as.numeric(ppp6r3_frax_min$P.I)
ppp6r3_frax_min$pos_hg38 = as.numeric(ppp6r3_frax_min$pos_hg38)


frax_racer = RACER::formatRACER(assoc_data = ppp6r3_frax_min, chr_col = 3, pos_col = 1, p_col = 2)

frax_racer = RACER::ldRACER(assoc_data = frax_racer, rs_col = 4, pops = "EUR", lead_snp = "rs10047483")
#

#RACER::singlePlotRACER(assoc_data = frax_racer, chr = 11, build = "hg38", plotby = "coord", label_lead = T)
mirrorPlotRACER(assoc_data1 = frax_racer, assoc_data2 = eqtl_racer, chr = 11, plotby = "coord",build = "hg38",label_lead = T,name1 = "Fracture SNPs", name2 = "PPP6R3 eQTL - Thyroid", start_plot =68460724-1000000, end_plot = 68460724+1000000)

scatterPlotRACER(assoc_data1 = frax_racer, assoc_data2 = eqtl_racer, chr = 11, name1 = "BMD GWAS p_values", name2 = "PPP6R3 eQTL p_values", region_start = 68460724-1000000, region_end = 68460724+1000000, ld_df = 1)


