library(ggplot2)
library(biomaRt)
library(topGO)
library(RACER)
library(LDlinkR)
library(GenomicRanges)

set.seed(8675309)

setwd("../../../../scratch/bma8ne/GTEx_v8/src")


# read_in -----------------------------------------------------------------

#read in hg38 version of eBMD GWAS

bmd = as.data.frame(data.table::fread("../BMD_hg38", header=T))

#saveRDS(bmd, "assets/bmd.rds")


#read in FastENLOC sig files for BMD

enloc_sig = read.table("../enloc_out/all_enloc_bmd.sig.out")
colnames(enloc_sig) = c("cluster_name","nSNPs", "cluster_PIP_eqtl", "cluster_PIP_gwas_noeqtl","cluster_PIP_gwas_w_eqtl", "RCP","tissue")
enloc_sig$tissue = sapply(strsplit(enloc_sig$tissue,"[.]"),"[",1) #format tissue name

#saveRDS(enloc_sig, "assets/enloc_sig.rds")

#split ensembl name
enloc_sig$ens = sapply(strsplit(enloc_sig$cluster_name, ":"),"[",1)

#get annotation, remove everything but coding genes
genes = unique(enloc_sig$ens)
mart = useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
ovr = getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol","gene_biotype"),values=genes,mart= mart)

keep = ovr[which(ovr$gene_biotype == "protein_coding"),"ensembl_gene_id"]
enloc_sig_prot = enloc_sig[which(enloc_sig$ens %in% keep),] #protein coding genes

#saveRDS(enloc_sig_prot, "assets/enloc_sig_prot.rds")

#RCP >=0.1
enloc_sig_prot_0.1 = enloc_sig_prot[which(enloc_sig_prot[,"RCP"]>=0.1),]
#RCP >=0.5
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
#from https://www.gtexportal.org/home/datasets
thyroid = as.data.frame(data.table::fread("../GTEx_Analysis_v8_eQTL_all_associations/Thyroid.allpairs.txt",stringsAsFactors = FALSE,header = TRUE, sep="\t"))#readtissue file
thyroid = thyroid[grep(pattern = "chr11_",x = thyroid$variant_id),]

#get RSIDs from GTEx lookup tables
#from https://www.gtexportal.org/home/datasets
chr11_lookup = read.table("../GTEx_v8_lookup/chr11_GTEx_V8_lookup.txt", header=T)

thyroid = merge(thyroid, chr11_lookup, by.x = "variant_id", by.y = "variant_id")

ppp6r3_thyroid_eqtl = thyroid[grep("ENSG00000110075",thyroid$gene_id),]

#saveRDS(ppp6r3_thyroid_eqtl, "assets/ppp6r3_thyroid_eqtl.rds")
#saveRDS(thyroid, "assets/all_thyroid_chr11_gtex_eqtl.rds")

#whole blood eqtl
#from https://www.gtexportal.org/home/datasets
blood = as.data.frame(data.table::fread("../GTEx_Analysis_v8_eQTL_all_associations/Whole_Blood.allpairs.txt",stringsAsFactors = FALSE,header = TRUE, sep="\t"))#readtissue file
blood = blood[grep(pattern = "chr19_",x = blood$variant_id),]

#get RSIDs from GTEx lookup tables
#from https://www.gtexportal.org/home/datasets
chr19_lookup = read.table("../GTEx_v8_lookup/chr19_GTEx_V8_lookup.txt", header=T)

blood = merge(blood, chr19_lookup, by.x = "variant_id", by.y = "variant_id")

gpatch1_blood_eqtl = blood[grep("ENSG00000076650",blood$gene_id),]

#saveRDS(gpatch1_blood_eqtl, "assets/gpatch1_blood_eqtl.rds")
#saveRDS(blood, "assets/all_blood_chr19_gtex_eqtl.rds")


#### prioritized genes (merge TWAS and ENLOC), superset pruning and enrichment of prioritized genes in bone genes ####

#merge enloc and TWAS
prioritized_bonf_0.1 = merge(enloc_sig_prot_0.1,multixcan_prot_bonf_sig, by.x="ens", by.y="gene")

length(unique(prioritized_bonf_0.1$ens)) #512

#superset
superset = read.table("../superset_humanized_proteinCoding.txt", stringsAsFactors = F)
colnames(superset) = c("SYMBOL", "ENSEMBL")

#number prioritized genes in bone superset
length(unique(prioritized_bonf_0.1[which(prioritized_bonf_0.1$ens %in% superset$ENSEMBL),"ens"])) #66


#enrichment of prioritized genes in bone genes
#total protein coding genes
total_prot = unique(intersect(multixcan_prot$gene, enloc_sig_prot$ens))

length(unique(superset[which(superset$ENSEMBL %in% total_prot), "ENSEMBL"])) #1,399/1,484 bone genes are 

superset_pruned = superset[which(superset$ENSEMBL %in% total_prot),]

#number of prioritized genes 
length(which(unique(prioritized_bonf_0.1$ens) %in% superset_pruned$ENSEMBL))

a = 66 #prioritized genes in known bone list
b = length(unique(prioritized_bonf_0.1$ens)) - a #priortized genes not in known bone list
c = length(unique(total_prot[which( (total_prot %in% prioritized_bonf_0.1$ens == F) & (total_prot %in% superset_pruned$ENSEMBL == TRUE) )]))#not prioritized and IN bone
d = length(unique(total_prot[which( (total_prot %in% prioritized_bonf_0.1$ens == F) & (total_prot %in% superset_pruned$ENSEMBL == F) )]))#not prioritized and not bone

mat = matrix(nrow=2, ncol=2)
mat[1,] = c(a,c)
mat[2,] = c(b,d)

fisher.test(mat,alternative = "g") #OR=1.72
fisher.test(mat,alternative = "g")$p.value #pval = 0.0001020339


#### prioritized gene localization ####

#where are prioritized genes located, locus wise? how are they distributed across the genome?

#get gwas lead snps location +/- 1 mbp
#from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6358485/
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

length(unique(bmd_lead_hg38$refsnp_id))#1103

#get prioritized gene start and end from ensembl
prioritized_unique = unique(prioritized_bonf_0.1$ens)
ensembl=useMart("ensembl", dataset = "hsapiens_gene_ensembl")



prioritized_gene_pos <- getBM(attributes=c('ensembl_gene_id',
                                           'chromosome_name',
                                           'start_position',
                                           'end_position'),
                              filters ="ensembl_gene_id", values =prioritized_unique, mart = ensembl, uniqueRows=TRUE)

#why no X chromosome? because multixcan did not include sex chromosomes
prioritized_bonf_0.1 = merge(prioritized_bonf_0.1, prioritized_gene_pos, by.x = "ens",by.y="ensembl_gene_id") #merge in positions

#overlap, count number of overlaps per lead snp
#first convert snps to loci, where a locus is a snps pos +/- 1 Mbp (1000000 bp)
bmd_lead_hg38$chr_name = paste0("chr",bmd_lead_hg38$chr_name)
bmd_lead_hg38$start_locus = as.numeric(bmd_lead_hg38$chrom_start) - 1000000
bmd_lead_hg38$end_locus = as.numeric(bmd_lead_hg38$chrom_end) + 1000000

bmd_lead_granges = makeGRangesFromDataFrame(bmd_lead_hg38, seqnames.field = "chr_name", start.field = "start_locus",end.field = "end_locus", keep.extra.columns = T)


#now convert genes granges

prioritized_gene_pos$chr_name = paste0("chr",prioritized_gene_pos$chromosome_name)

prioritized_gene_granges = makeGRangesFromDataFrame(prioritized_gene_pos, seqnames.field = "chr_name", start.field = "start_position",end.field = "end_position", keep.extra.columns = T)

overlaps = findOverlaps(bmd_lead_granges, prioritized_gene_granges)
overlaps = as.data.frame(overlaps)
# 456 genes overlap 542 snps (919 total overlaps, so a gene can overlap multiple snps)

x1 = bmd_lead_hg38[overlaps$queryHits,]
y = prioritized_gene_pos[overlaps$subjectHits,]
z = cbind(x1,y)

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
t1<-GenTable(GOdata, classic=result, topNodes=length(result@score), numChar=1000)
t1$ontology = "MF"
head(t1)
###CC###
GOdata <- new("topGOdata", ontology = "CC", allGenes = geneList,
              annot = annFUN.org, mapping='org.Hs.eg.db', ID='ENSEMBL')
test.stat<-new("classicCount", testStatistic = GOFisherTest, name='Fisher test')
result<-getSigGroups(GOdata,test.stat)
t2<-GenTable(GOdata, classic=result, topNodes=length(result@score), numChar=1000)
t2$ontology = "CC"
head(t2)
###BP###
GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,
              annot = annFUN.org, mapping='org.Hs.eg.db', ID='ENSEMBL')
test.stat<-new("classicCount", testStatistic = GOFisherTest, name='Fisher test')
result<-getSigGroups(GOdata,test.stat)
t3<-GenTable(GOdata, classic=result, topNodes=length(result@score), numChar=1000)
t3$ontology = "BP"
head(t3)
####
t.all = NULL
t.all<-rbind(t1,t2,t3)
t.all$classic<-as.numeric(as.character(t.all$classic))

#saveRDS(t.all, "assets/t.all")


#### IMPC cross-ref ####

#from impc_SOLR.R
##REMOVED SKULL.
impc = readRDS("assets/impc_out_5e2.rds") #all genes with at least one pvalue <=0.05 (male, female, genotype)
impc = impc[,"marker_symbol"]

impc = unique(impc) #2,360 genes

#convert mouse to human gene names

human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")


impc_hum = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = impc , mart = mouse, attributesL = c("hgnc_symbol", "ensembl_gene_id"), martL = human, uniqueRows=T)

impc_hum = unique(impc_hum) 

length(unique(impc_hum$MGI.symbol)) #2,245 unique mouse genes in impc data have homologues in humans. some, like Sirpa, have multiple human homologs.
length(unique(impc_hum$Gene.stable.ID)) #2,367 total human homologues

length(unique(prioritized_bonf_0.1[which(prioritized_bonf_0.1$ens %in% impc_hum$Gene.stable.ID),"ens"])) #64 genes with nominally significant alteration in BMD

x = (unique(prioritized_bonf_0.1[which(prioritized_bonf_0.1$ens %in% impc_hum$Gene.stable.ID),"ens"]))
length(which(x %in% superset_pruned$ENSEMBL == FALSE)) # 49 of these are unique to IMPC, not in superset



#how many of the 512 prioritized genes were tested by the IMPC?
#get all impc genes first. These have "successful" status and dont include skull
load("assets/impc_results_raw_processed_noskull") #all impc BMD genes. These have "successful" status and dont include skull
impc_all = out
impc_all = impc_all[,"marker_symbol"]
impc_all = unique(impc_all)
#convert mouse to human gene names

human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")


impc_all_hum = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = impc_all , mart = mouse, attributesL = c("hgnc_symbol", "ensembl_gene_id"), martL = human, uniqueRows=T)

impc_all_hum = unique(impc_all_hum) 

length(unique(impc_all_hum$MGI.symbol)) #5,154 unique mouse genes in impc data have homologues in humans. some, like Sirpa, have multiple human homologs.
length(unique(impc_all_hum$Gene.stable.ID)) #5,375 total human homologues

x = unique(prioritized_bonf_0.1$ens)
length(which(x %in% impc_all_hum$Gene.stable.ID)) #142 of the 512 (27.7%) prioritized genes are assayed by IMPC

#novel prioritized genes
novel_0.5_ens = prioritized_bonf_0.1[which(prioritized_bonf_0.1$RCP >= 0.5),] #RCP >= 0.5
novel_0.5_ens = novel_0.5_ens[-which(novel_0.5_ens$ens %in% superset_pruned$ENSEMBL),] #novel not in known bone gene list
novel_0.5_ens = novel_0.5_ens[-which(novel_0.5_ens$ens %in% impc_hum$Gene.stable.ID),] # novel also not significant by IMPC

##### phenomexcan ####
#data from https://advances.sciencemag.org/content/6/37/eaba2083
phx = read.csv("phenomexcan_results.csv", stringsAsFactors = F)
phx = phx[which(phx$trait == "3148_raw-Heel_bone_mineral_density_BMD"),]
phx = phx[which(phx$gene_type == "protein_coding"),]
phx = unique(phx)

length(which(phx$gene_id %in% prioritized_bonf_0.1$ens)) #55


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
lead_rsid = bmd_lead_hg38[which(bmd_lead_hg38$chr_name == "chr11" & bmd_lead_hg38$start_locus <= (68460724) & bmd_lead_hg38$end_locus >= (68460724)),]

token = "token_here" #insert your token here, by registereing on the LDLink wehsite
LDmatrix(c(lead_rsid$refsnp_id, "rs10047483"), pop = "EUR", r2d = "r2", token = token) #LD between snp and bmd lead snps in ppp6r3 locus

# rs10047483 ref G alt A in gtex, with positive slope. effect allele is G in GWAS, with positive slope. However, in GTEx, eQTL allele is the ALT allele. Therefore, the G allele effect is negative in GTEx, and positive in GWAS

##LRP5 = "ENSG00000162337
LDmatrix(c(lead_rsid$refsnp_id, "rs3736228"), pop = "EUR", r2d = "r2", token = token) #LD between snp and bmd lead snps in ppp6r3 locus
LDpair("rs10047483","rs3736228", pop = "EUR", token = token)
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

ppp6r3_frax_enloc = enloc_frax_sig[which(enloc_frax_sig$ens == "ENSG00000110075"),]


#get annotation, remove everything but coding genes
genes = unique(enloc_frax_sig$ens)
mart = useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
ovr = getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol","gene_biotype"),values=genes,mart= mart)

keep = ovr[which(ovr$gene_biotype == "protein_coding"),"ensembl_gene_id"]
enloc_frax_sig_prot = enloc_frax_sig[which(enloc_frax_sig$ens %in% keep),]



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


frax_ppp6r3 = frax[which(frax$chr_hg38 == "chr11" & frax$BP <= (68228199+1000000) & frax$BP >= (68228199-1000000)),] #check locations
#ppp6r3_thyroid_eqtl[which(ppp6r3_thyroid_eqtl$rs_id_dbSNP151_GRCh38p7 == "rs35989399"),]
LDpair("rs10047483","rs35989399", pop = "EUR", token = token)#0.458


####FIGURE 1####

###1B###
#get all start sites for genes in prioritized list
ensembl=useMart("ensembl", dataset = "hsapiens_gene_ensembl")

total.pos = merge(multixcan_prot, enloc_sig_prot, by.y="ens", by.x="gene")


total.pos.bm <- getBM(attributes=c('ensembl_gene_id','chromosome_name',
                                   'start_position'),
                      filters ="ensembl_gene_id", values =total.pos$gene, mart = ensembl, uniqueRows=TRUE)

total.pos = merge(total.pos, total.pos.bm, by.y="ensembl_gene_id", by.x="gene")

#get genes whose start sites are within 1 MBP of RUNX2 start site
#chr6, pos 45328157

x = total.pos[which((total.pos$chromosome_name == "6") & (total.pos$start_position >= 45328157-1500000) & (total.pos$start_position <= 45328157+1500000)),]

xx=data.frame()

for(genes in unique(x$gene)){
  s = subset(x, x$gene==genes)
  s = s[which(as.numeric(s$RCP) == max(as.numeric(s$RCP))),]
  xx = rbind(xx,s)
}

xx=xx[order(xx$start_position,decreasing = F),]

xx$color = "black"

xx$color[seq(2, nrow(xx), 2)] = "red"

xx$shape=1
xx[which(xx$gene_name == "RUNX2"),"shape"] = 17


plot = ggplot(xx, aes(y = reorder(gene_name,-start_position), x=-log10(bonf))) + geom_point(color=xx$color, shape=xx$shape) + theme_classic() + theme(axis.text.y = element_text(color=xx$color), title=element_text(size=6),axis.title=element_text(size=6),axis.text = element_text(size=6))
ggsave(filename="runx2_twas.pdf", plot=plot, device="pdf",
       path="plots/", height=6, width=1.5, units="in", dpi=500)


plot = ggplot(xx, aes(y = reorder(gene_name,-start_position), x=RCP)) + geom_point(color=xx$color, shape=xx$shape) + theme_classic() + theme(axis.text.y = element_text(color=xx$color), title=element_text(size=6),axis.title=element_text(size=6),axis.text = element_text(size=6))+ scale_x_reverse()

ggsave(filename="runx2_rcp.pdf", plot=plot, device="pdf",
       path="plots/", height=6, width=1.5, units="in", dpi=500)


plot = ggplot(xx, aes(y = reorder(start_position,-start_position), x=-log10(bonf))) + geom_point(color=xx$color, shape=xx$shape) + theme_classic() + theme(title=element_text(size=6),axis.title=element_text(size=6),axis.text = element_text(size=6))
ggsave(filename="runx2_pos.pdf", plot=plot, device="pdf",
       path="plots/", height=6, width=1.5, units="in", dpi=500)

##same for RERE, 1C

#chr1, pos 8352397
x = total.pos[which((total.pos$chromosome_name == "1") & (total.pos$start_position >= 8352397-1500000) & (total.pos$start_position <= 8352397+1500000)),]

xx=data.frame()

for(genes in unique(x$gene)){
  s = subset(x, x$gene==genes)
  s = s[which(as.numeric(s$RCP) == max(as.numeric(s$RCP))),]
  xx = rbind(xx,s)
}

xx=xx[order(xx$start_position,decreasing = F),]

xx$color = "black"

xx$color[seq(2, nrow(xx), 2)] = "red"

xx$shape=1
xx[which(xx$gene_name == "RERE"),"shape"] = 17

plot = ggplot(xx, aes(y = reorder(gene_name,-start_position), x=-log10(bonf))) + geom_point(color=xx$color, shape=xx$shape) + theme_classic() + theme(axis.text.y = element_text(color=xx$color), title=element_text(size=6),axis.title=element_text(size=6),axis.text = element_text(size=6))

ggsave(filename="rere_twas.pdf", plot=plot, device="pdf",
       path="plots/", height=6, width=1.5, units="in", dpi=500)


plot = ggplot(xx, aes(y = reorder(gene_name,-start_position), x=RCP)) + geom_point(color=xx$color, shape=xx$shape) + theme_classic() + theme(axis.text.y = element_text(color=xx$color), title=element_text(size=6),axis.title=element_text(size=6),axis.text = element_text(size=6))+ scale_x_reverse()

ggsave(filename="rere_rcp.pdf", plot=plot, device="pdf",
       path="plots/", height=6, width=1.5, units="in", dpi=500)


plot = ggplot(xx, aes(y = reorder(start_position,-start_position), x=-log10(bonf))) + geom_point(color=xx$color, shape=xx$shape) + theme_classic() + theme(title=element_text(size=6),axis.title=element_text(size=6),axis.text = element_text(size=6))
ggsave(filename="RERE_pos.pdf", plot=plot, device="pdf",
       path="plots/", height=6, width=1.5, units="in", dpi=500)

####TRY DOING THIS WITH MIRRORPLOTS INSTEAD###



#1D
#"z" and "bmd_lead_hg38" from "prioritized gene localization" section above

bmd_lead_hg38 = bmd_lead_hg38[-which(bmd_lead_hg38$chr_name == "chrX"),]

x=table(z$refsnp_id)

x=as.data.frame(x)

xx = bmd_lead_hg38[which(bmd_lead_hg38$refsnp_id %in% z$refsnp_id == FALSE),"refsnp_id"]
xx = as.data.frame(xx)
xx$Freq = 0
colnames(xx)[1] = "Var1"

xxx = rbind(x,xx)


plot = ggplot(xxx, aes(Freq)) +
  geom_histogram(binwidth = 0.5, fill=c("black")) + scale_x_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9)) + theme_classic()

ggsave(filename="hist.pdf", plot=plot, device="pdf",
       path="plots/", height=2.5, width=6, units="in", dpi=500)




#### FIGURE 2 ####
#GO plot
#get t3 (BP) from GO analysis section above
t3$classic = as.numeric(t3$classic)
t3.s = t3[which(t3$classic <= 0.05),]

t3.ss = t3.s[c(1,2,3,5,6,9,10,12,13,15,17,18,19,27,28,35,37,38,43:49,52,56,58,60,64,69,71,72,84,85,87,95,97,111,113),]

t3.ss$classic = -log10(t3.ss$classic)

plot = ggplot(t3.ss, aes(x = reorder(Term,-classic), y=classic)) + geom_point(size=.6) + theme_classic() + theme(axis.text.x = element_text(angle = 90), title=element_text(size=6),axis.title=element_text(size=6),axis.text = element_text(size=6))


ggsave(filename="GO.pdf", plot=plot, device="pdf",
       path="plots/", height=3.5, width=7, units="in", dpi=500)

###

#gene significant as both known bone gene and IMPC gene, plot RCP/TWAS as above, and also IMPC data
#get data from IMPC section above

# which prioritized genes are in known bone gene and also IMPC
x = (prioritized_bonf_0.1[which(prioritized_bonf_0.1$ens %in% impc_hum$Gene.stable.ID),])
x = (x[which(x$ens %in% superset_pruned$ENSEMBL),])

#GPATCH1, chr19, 33080899
x = total.pos[which((total.pos$chromosome_name == "19") & (total.pos$start_position >= 33080899-1500000) & (total.pos$start_position <= 33080899+1500000)),]

xx=data.frame()

for(genes in unique(x$gene)){
  s = subset(x, x$gene==genes)
  s = s[which(as.numeric(s$RCP) == max(as.numeric(s$RCP))),]
  xx = rbind(xx,s)
}
xx = xx[,c("gene_name", "start_position", "bonf", "RCP")]
xx = unique(xx)

xx=xx[order(xx$start_position,decreasing = F),]

xx$color = "black"

xx$color[seq(2, nrow(xx), 2)] = "red"

xx$shape=1
xx[which(xx$gene_name == "GPATCH1"),"shape"] = 17


plot = ggplot(xx, aes(y = reorder(gene_name,-start_position), x=-log10(bonf))) + geom_point(color=xx$color, shape=xx$shape) + theme_classic() + theme(axis.text.y = element_text(color=xx$color), title=element_text(size=6),axis.title=element_text(size=6),axis.text = element_text(size=6))

ggsave(filename="gpatch1_twas.pdf", plot=plot, device="pdf",
       path="plots/", height=6, width=1.5, units="in", dpi=500)


plot = ggplot(xx, aes(y = reorder(gene_name,-start_position), x=RCP)) + geom_point(color=xx$color, shape=xx$shape) + theme_classic() + theme(axis.text.y = element_text(color=xx$color), title=element_text(size=6),axis.title=element_text(size=6),axis.text = element_text(size=6))+ scale_x_reverse()

ggsave(filename="gpatch1_rcp.pdf", plot=plot, device="pdf",
       path="plots/", height=6, width=1.5, units="in", dpi=500)

plot = ggplot(xx, aes(y = reorder(start_position,-start_position), x=-log10(bonf))) + geom_point(color=xx$color, shape=xx$shape) + theme_classic() + theme(title=element_text(size=6),axis.title=element_text(size=6),axis.text = element_text(size=6))

ggsave(filename="gpatch1_pos.pdf", plot=plot, device="pdf",
       path="plots/", height=6, width=1.5, units="in", dpi=500)




#mirrorplot gpatch1?

gpatch1_gwas_min = bmd[,c("pos_hg38","P.NI","CHR","RSID")]
gpatch1_gwas_min = gpatch1_gwas_min[which(gpatch1_gwas_min$CHR == 19),]#remove non chromosome 19 snps otherwise theres repeats

gpatch1_gwas_min$P.NI = as.numeric(gpatch1_gwas_min$P.NI)
gpatch1_gwas_min$pos_hg38 = as.numeric(gpatch1_gwas_min$pos_hg38)
gpatch1_gwas_min$col = "notleadsnp"

#get BMD lead SNPS fpr GPATCH1
lead_snps = read.csv("../BMD_morris_lead_snps_hg19.csv", header=F)
#gpatch1 starts at chr19:33,571,786 on hg19

#lead snp rsids within +/- a megabase of the "TSS" of gpatch1 (ENCODE beginning of gene)
lead_rsid = lead_snps[which(lead_snps$V3 == 19 & lead_snps$V4 <= (33571786+1000000) & lead_snps$V4 >= (33571786-1000000)),]

gpatch1_gwas_min[which(gpatch1_gwas_min$RSID %in% lead_rsid$V2),"col"]="leadsnp"

gpatch1_gwas_min = unique(gpatch1_gwas_min)

gwas_racer = RACER::formatRACER(assoc_data = gpatch1_gwas_min, chr_col = 3, pos_col = 1, p_col = 2)

gwas_racer = RACER::ldRACER(assoc_data = gwas_racer, rs_col = 4, pops = "EUR", lead_snp = "rs11881367")

#
gpatch1_blood_eqtl$chr = sapply(strsplit(gpatch1_blood_eqtl$chr, "chr"),"[",2)
gpatch1_blood_eqtl$pval_nominal = as.numeric(gpatch1_blood_eqtl$pval_nominal)

eqtl_racer = RACER::formatRACER(assoc_data = gpatch1_blood_eqtl, chr_col = 10, pos_col = 11, p_col = 7)


eqtl_racer = RACER::ldRACER(assoc_data = eqtl_racer, rs_col = 15, pops = "EUR", lead_snp = "rs11881367")

eqtl_racer = unique(eqtl_racer)

plot = mirrorPlotRACER(assoc_data1 = gwas_racer, assoc_data2 = eqtl_racer, chr = 19, plotby = "coord",build = "hg38",label_lead = T,name1 = "eBMD SNPs", name2 = "GPATCH1 eQTL - Whole Blood", start_plot =33080899-500000, end_plot = 33080899+500000)

plot = plot + theme(axis.text.y = element_text(color=xx$color), title=element_text(size=6),axis.title=element_text(size=6),axis.text = element_text(size=6))


ggsave(filename="gpatch1_mirror.pdf", plot=plot, device="pdf",
       path="plots/", height=5, width=5, units="in", dpi=500)


plot = scatterPlotRACER(assoc_data1 = gwas_racer, assoc_data2 = eqtl_racer, chr = 19, name1 = "BMD GWAS p_values", name2 = "GPATCH1 eQTL p_values", region_start = 33080899-500000, region_end = 33080899+500000, ld_df = 1)

ggsave(filename="gpatch1_scatter.pdf", plot=plot, device="pdf",
       path="plots/", height=5, width=8.5, units="in", dpi=500)



#IMPC data
impc = readRDS("assets/impc_out_5e2.rds") #all genes with at least one pvalue <=0.05 (male, female, genotype)

impc_gpatch1 = impc[which(impc$marker_symbol == "Gpatch1"),]




####FIGURE 3 ####
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
plot = ggplot(ppp6r3_gwas_min, aes(pos_hg38,-log10(P.NI),color=col, label=RSID)) + geom_point(size=1, alpha=0.5) +xlim(68460724-1000000,68460724+1000000) +ylim(0,102) + ggtitle("eBMD GWAS independent SNPs within 1 Mbp of PPP6R3 start site")+
  geom_text_repel(data = subset(ppp6r3_gwas_min, col == "leadsnp"), point.padding = 0.5, box.padding = 0.7, force = 50, nudge_y = 150, direction="x", segment.size = 0.2, size=3) + theme_classic() + 
  annotate("point",ppp6r3_gwas_min$pos_hg38[which(ppp6r3_gwas_min$col=="leadsnp")],-log10(ppp6r3_gwas_min$P.NI[which(ppp6r3_gwas_min$col=="leadsnp")]),size=1) +
  theme(axis.text.y = element_text(color=xx$color), title=element_text(size=6),axis.title=element_text(size=6),axis.text = element_text(size=6), legend.position = "none")

ggsave(filename="ppp6r3_locus.pdf", plot=plot, device="pdf",
       path="plots/", height=2, width=6.5, units="in", dpi=500)


gwas_racer = RACER::formatRACER(assoc_data = ppp6r3_gwas_min, chr_col = 3, pos_col = 1, p_col = 2)
gwas_racer = RACER::ldRACER(assoc_data = gwas_racer, rs_col = 4, pops = "EUR", lead_snp = "rs3736228")
#gwas_racer$LD = 0
#gwas_racer$LD_BIN = "0.2-0.0"
#gwas_racer$newlabel = NA
#gwas_racer[which(gwas_racer$RS_ID %in% lead_rsid),"newlabel"]= gwas_racer[which(gwas_racer$RS_ID %in% lead_rsid),"RS_ID"]

# GENES FOR ABOVE FROM HERE
RACER::singlePlotRACER(assoc_data = gwas_racer, chr = 11, build = "hg38", plotby = "coord", label_lead = F,start_plot =68460724-1000000, end_plot = 68460724+1000000)





#3B

#PPP6R3, chr11, 68228199
x = total.pos[which((total.pos$chromosome_name == "11") & (total.pos$start_position >= 68460724-1000000) & (total.pos$start_position <= 68460724+1000000)),]

xx=data.frame()

for(genes in unique(x$gene)){
  s = subset(x, x$gene==genes)
  s = s[which(as.numeric(s$RCP) == max(as.numeric(s$RCP))),]
  xx = rbind(xx,s)
}

xx = xx[,c("gene_name", "start_position", "bonf", "RCP")]
xx = unique(xx)
xx=xx[order(xx$start_position,decreasing = F),]

xx$color = "black"

xx$color[seq(2, nrow(xx), 2)] = "red"

xx$shape=1
xx[which(xx$gene_name == "PPP6R3"),"shape"] = 17


plot = ggplot(xx, aes(y = reorder(gene_name,-start_position), x=-log10(bonf))) + geom_point(color=xx$color, shape=xx$shape) + theme_classic() + theme(axis.text.y = element_text(color=xx$color), title=element_text(size=6),axis.title=element_text(size=6),axis.text = element_text(size=6))

ggsave(filename="ppp6r3_twas.pdf", plot=plot, device="pdf",
       path="plots/", height=6, width=1.5, units="in", dpi=500)


plot = ggplot(xx, aes(y = reorder(gene_name,-start_position), x=RCP)) + geom_point(color=xx$color, shape=xx$shape) + theme_classic() + theme(axis.text.y = element_text(color=xx$color), title=element_text(size=6),axis.title=element_text(size=6),axis.text = element_text(size=6))+ scale_x_reverse()

ggsave(filename="ppp6r3_rcp.pdf", plot=plot, device="pdf",
       path="plots/", height=6, width=1.5, units="in", dpi=500)

plot = ggplot(xx, aes(y = reorder(start_position,-start_position), x=-log10(bonf))) + geom_point(color=xx$color, shape=xx$shape) + theme_classic() + theme(title=element_text(size=6),axis.title=element_text(size=6),axis.text = element_text(size=6))

ggsave(filename="ppp6r3_pos.pdf", plot=plot, device="pdf",
       path="plots/", height=6, width=1.5, units="in", dpi=500)







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




####supp tables####
#1, 2,156 protein coding twas genes
write.csv(x=multixcan_prot_bonf_sig,file = "tables/supp_1.csv", quote = F,row.names = F)

#2, 1,182 genes with RCP >= 0.1
ensembl <- useEnsembl("snp",dataset = "hsapiens_snp")

#get genomic position
ensembl=useMart("ensembl", dataset = "hsapiens_gene_ensembl")
x <- getBM(attributes=c("ensembl_gene_id",
                        "external_gene_name"),
           filters ="ensembl_gene_id", values =enloc_sig_prot_0.1$ens, mart = ensembl, uniqueRows=TRUE)

x = unique(x)
enloc_sig_prot_0.1 = merge(enloc_sig_prot_0.1, x, by.x="ens", by.y="ensembl_gene_id")

write.csv(x=enloc_sig_prot_0.1,file = "tables/supp_2.csv", quote = F,row.names = F)

#3, 512 protein coding genes significant by both TWAS and colocalization @ 0.1
write.csv(x=prioritized_bonf_0.1,file = "tables/supp_3.csv", quote = F,row.names = F)

#4, number of significantly colocalizing genes (in above table) per tissue

temp = as.data.frame(matrix(nrow = 49, ncol=2))
colnames(temp) = c("tissue","number_unique_coloc_genes")
counter=1
for(i in unique(prioritized_bonf_0.1$tissue)){
  print(i)
  x = subset(prioritized_bonf_0.1, tissue == i)
  num = length(unique(x$gene_name))
  temp[counter,] = c(i, as.numeric(num))
  counter = counter +1
}
temp$number_unique_coloc_genes = as.numeric(temp$number_unique_coloc_genes)

write.csv(x=temp,file = "tables/supp_4.csv", quote = F,row.names = F)

#5, known bone gene list
write.csv(x=superset_pruned,file = "tables/supp_5.csv", quote = F,row.names = F)

#6, which of the 512 prioritized genes are in known bone gene list (66)
#length(unique(prioritized_bonf_0.1[which(prioritized_bonf_0.1$ens %in% superset_pruned$ENSEMBL),"ens"])) #66
x = (unique(prioritized_bonf_0.1[which(prioritized_bonf_0.1$ens %in% superset_pruned$ENSEMBL),])) #66
x = x[,c("ens","gene_name")]
x = unique(x)

write.csv(x=x,file = "tables/supp_6.csv", quote = F,row.names = F)

#7, GO enrichments for the 512 prioritized genes, subset pval <= 0.05
t.all = readRDS("assets/t.all") # topGO gene ontology for the 512 prioritized genes
t.all = t.all[t.all$classic <=0.05,]

write.table(x=t.all,file = "tables/supp_7.csv", quote = F,row.names = F,sep = "\t")

#8, Novel prioritized genes (137)
x = novel_0.5_ens[,c("gene_name","ens","RCP","tissue","bonf")]
xx = as.data.frame(matrix(nrow=1, ncol=5))
colnames(xx) = colnames(x)

for(i in unique(x$gene_name)){
  
  s = as.data.frame(subset(x, gene_name == i))
  s = s[which(s$RCP == max(s$RCP)),]
  xx = rbind(xx,s)
  
}

write.table(x=xx,file = "tables/supp_8.csv", quote = F,row.names = F,sep = "\t") #TLN2 has two entries, same RCP in nerve_tibial and esophagus_muscularis


#### supp figure 1 ####
#plot estrada lsbmd and fnbmd in ppp6r3 region, highlight significant eQTL SNPs

#read in lsbmd. converted to hg38 coordinates using get_estrada_gwas_hg38_pos.R
lsbmd_38 = read.table("lsbmd_38.txt", header=T, stringsAsFactors = F)

#get snps in PPP6R3 region (+/- 1mbp)
lsbmd_38 = lsbmd_38[which(lsbmd_38$chr == 11),]
lsbmd_38 = lsbmd_38[which(lsbmd_38$pos38 <= (68460724+1000000)),]
lsbmd_38 = lsbmd_38[which(lsbmd_38$pos38 >= (68460724-1000000)),]

eqtl_sig = ppp6r3_thyroid_eqtl[which(as.numeric(ppp6r3_thyroid_eqtl$pval_nominal) <=0.05),]

lsbmd_38$eqtl = "noteqtl"
lsbmd_38[which(lsbmd_38$MarkerName %in% eqtl_sig$rs_id_dbSNP151_GRCh38p7),"eqtl"] = "eqtl"

lsbmd_38$pval_log = -log10(lsbmd_38$P.value)

colors = c("red", "black")
plot = ggplot(lsbmd_38, aes(pos38,pval_log,color=eqtl)) + geom_point(size=1, alpha=0.5) +xlim(68460724-1000000,68460724+1000000) +
  ylim(0,11) + scale_fill_manual(values=colors) + scale_colour_manual(values=colors) + geom_hline(yintercept=-log10(5e-8), linetype="dashed", color = "black", size=0.5)

ggsave(filename="lsbmd_eqtl_highlighted.pdf", plot=plot, device="pdf",
       path="plots/", height=5, width=5, units="in", dpi=500)
#same for fnbmd
fnbmd_38 = read.table("fnbmd_38.txt", header=T, stringsAsFactors = F)

#get snps in PPP6R3 region (+/- 1mbp)
fnbmd_38 = fnbmd_38[which(fnbmd_38$chr == 11),]
fnbmd_38 = fnbmd_38[which(fnbmd_38$pos38 <= (68460724+1000000)),]
fnbmd_38 = fnbmd_38[which(fnbmd_38$pos38 >= (68460724-1000000)),]

eqtl_sig = ppp6r3_thyroid_eqtl[which(as.numeric(ppp6r3_thyroid_eqtl$pval_nominal) <=0.05),]

fnbmd_38$eqtl = "noteqtl"
fnbmd_38[which(fnbmd_38$MarkerName %in% eqtl_sig$rs_id_dbSNP151_GRCh38p7),"eqtl"] = "eqtl"

fnbmd_38$pval_log = -log10(fnbmd_38$P.value)

colors = c("red", "black")
plot = ggplot(fnbmd_38, aes(pos38,pval_log,color=eqtl)) + geom_point(size=2, alpha=0.5) +xlim(68460724-1000000,68460724+1000000) +
  ylim(0,11) + scale_fill_manual(values=colors) + scale_colour_manual(values=colors) + geom_hline(yintercept=-log10(5e-8), linetype="dashed", color = "black", size=0.5)

ggsave(filename="fnbmd_eqtl_highlighted.pdf", plot=plot, device="pdf",
       path="plots/", height=5, width=5, units="in", dpi=500)

#also plot fnbmd and lsbmd mirrorplots with rs10047483 as lead snp
fnbmd_38 = read.table("fnbmd_38.txt", header=T, stringsAsFactors = F)
fnbmd_racer = RACER::formatRACER(assoc_data = fnbmd_38, chr_col = 13, pos_col = 14, p_col = 8)
fnbmd_racer = RACER::ldRACER(assoc_data = fnbmd_racer, rs_col = 1, pops = "EUR", lead_snp = "rs10047483") 
RACER::singlePlotRACER(assoc_data = fnbmd_racer, chr = 11, build = "hg38", plotby = "coord", gene_plot = "PPP6R3", label_lead = T,start_plot =68460724-1000000, end_plot = 68460724+1000000)

ppp6r3_thyroid_eqtl = readRDS("assets/ppp6r3_thyroid_eqtl.rds")
ppp6r3_thyroid_eqtl$chr = sapply(strsplit(ppp6r3_thyroid_eqtl$chr, "chr"),"[",2)
ppp6r3_thyroid_eqtl$pval_nominal = as.numeric(ppp6r3_thyroid_eqtl$pval_nominal)

eqtl_racer = RACER::formatRACER(assoc_data = ppp6r3_thyroid_eqtl, chr_col = 10, pos_col = 11, p_col = 7)


eqtl_racer = RACER::ldRACER(assoc_data = eqtl_racer, rs_col = 15, pops = "EUR", lead_snp = "rs10047483")

eqtl_racer = unique(eqtl_racer)

plot = mirrorPlotRACER(assoc_data1 = fnbmd_racer, assoc_data2 = eqtl_racer, chr = 11, plotby = "coord",build = "hg38",label_lead = T,name1 = "eBMD SNPs", name2 = "PPP6R3 eQTL - Thyroid", start_plot =68460724-1000000, end_plot = 68460724+1000000)

ggsave(filename="fnbmd_mirrorplot.pdf", plot=plot, device="pdf",
       path="plots/", height=5, width=5, units="in", dpi=500)

LDpair("rs10047483","rs608343", pop = "EUR", token = "2ebcd9d51762")# r2=0.54 with the lead assayed snp in the region


#ggplot(fnbmd_38, aes(pos38,-log10(P.value))) + geom_point() +xlim(68460724-1000000,68460724+1000000) + ggtitle("eBMD GWAS independent SNPs within 1 Mbp of PPP6R3 start site")+ theme_bw()

lsbmd_38 = read.table("lsbmd_38.txt", header=T, stringsAsFactors = F)
lsbmd_racer = RACER::formatRACER(assoc_data = lsbmd_38, chr_col = 10, pos_col = 11, p_col = 8)
lsbmd_racer = RACER::ldRACER(assoc_data = lsbmd_racer, rs_col = 1, pops = "EUR", lead_snp = "rs10047483")
RACER::singlePlotRACER(assoc_data = lsbmd_racer, chr = 11, build = "hg38", plotby = "coord", gene_plot = "PPP6R3", label_lead = T,start_plot =68460724-1000000, end_plot = 68460724+1000000)

LDpair("rs10047483","rs3736228", pop = "EUR", token = "2ebcd9d51762") #r2=0.36 with the lead assayed snp in the region

plot = mirrorPlotRACER(assoc_data1 = lsbmd_racer, assoc_data2 = eqtl_racer, chr = 11, plotby = "coord",build = "hg38",label_lead = T,name1 = "eBMD SNPs", name2 = "PPP6R3 eQTL - Thyroid", start_plot =68460724-1000000, end_plot = 68460724+1000000)

ggsave(filename="lsbmd_mirrorplot.pdf", plot=plot, device="pdf",
       path="plots/", height=5, width=5, units="in", dpi=500)
#genes for above
#RACER::singlePlotRACER(assoc_data = gwas_racer, chr = 11, build = "hg38", plotby = "coord", label_lead = F,start_plot =68460724-1000000, end_plot = 68460724+1000000)
#ggsave(filename="ppp6r3_pos.pdf", plot=plot, device="pdf",
#       path="plots/", height=6, width=1.5, units="in", dpi=500)