library(biomaRt)
library(tidyverse)
library("org.Hs.eg.db")
###

setwd("C:/projects/GWAS_project/src")

#convert mouse to human
#use mgi homolog list
#http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt (Nov-22-2020)

homology = read.table("../data/mgi_homologs.txt",sep = "\t", header = T)

convertMousetoHuman = function(x){
  human=c()
  for(i in 1:length(x)){
    id = homology[which((homology$Common.Organism.Name == "mouse, laboratory") & tolower(homology$Symbol) == tolower(x[i])),"HomoloGene.ID"]
    hum = homology[which((homology$Common.Organism.Name == "human") & homology$HomoloGene.ID == id),"Symbol"]
    human=append(human,hum)
  }
  human=unique(human)
  return(human)
  
}


###


#Using AmiGO2, downloaded GO terms for the following terms (osteo, bone, ossif). Accessed 7/28/19
#Used filters: is_obsolete:False and idspace: GO
bone_terms = read.delim("../data/GO_term_bone.txt", stringsAsFactors = FALSE, header = FALSE)
#trim to exclude some terms

ex = c("monocyte","pyridine","backbone","megakaryocyte","hair","kidney","neuro","ureter","B cell","tolerance","tendon","muscle","heart","cardio","beak","nephric","tooth","chemotaxis","hemopoiesis","amniotic","wishful")
bone_terms = bone_terms[-(grep(pattern = paste(ex,collapse = "|"), x=bone_terms$V2,ignore.case = TRUE)),]

bone_terms = bone_terms$V1

#
osteo_terms = read.delim("../data/GO_term_osteo.txt", stringsAsFactors = FALSE, header = FALSE)
osteo_terms = osteo_terms$V1
#
ossif_terms = read.delim("../data/GO_term_ossif.txt", stringsAsFactors = FALSE, header = FALSE)
ossif_terms = ossif_terms$V1
#
terms = c(ossif_terms, bone_terms, osteo_terms)
terms = unique(terms)


#gets gene symbol, transcript_id and go_id for all genes annotated with terms, from Mus Ensembl
ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl") #uses mus ensembl annotations


genes_mus = as.data.frame(matrix(ncol=4))
no_annot=c()
for(i in terms){
  print(i)
  gene.data <- unname(getBM(attributes=c('external_gene_name', 'ensembl_gene_id', 'go_id', "gene_biotype"),
                            filters = 'go', values = i, mart = ensembl))
  
  if(nrow(gene.data) >0){
    colnames(gene.data) = colnames(genes_mus)
    genes_mus = rbind(genes_mus,gene.data)
  }else{
    no_annot = append(no_annot, i)
  }
  
}
genes_mus = genes_mus[which(genes_mus$V4 == "protein_coding"),]
genes_mus = genes_mus[-1,1]
genes_mus = unique(genes_mus)

humanized_mus = convertMousetoHuman(genes_mus)
#same for human
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl") #uses human ensembl annotations

genes_hum = as.data.frame(matrix(ncol=4))
no_annot=c()
for(i in terms){
  print(i)
  gene.data <- unname(getBM(attributes=c('external_gene_name', 'ensembl_gene_id', 'go_id', "gene_biotype"),
                            filters = 'go', values = i, mart = ensembl))
  
  if(nrow(gene.data) >0){
    colnames(gene.data) = colnames(genes_hum)
    genes_hum = rbind(genes_hum,gene.data)
  }else{
    no_annot = append(no_annot, i)
  }
  
}
genes_hum = genes_hum[which(genes_hum$V4 == "protein_coding"),]
genes_hum = genes_hum[-1,1]
genes_hum = unique(genes_hum)



genes = c(humanized_mus, genes_hum)

genes=tolower(genes)
unq_genes = unique(genes)


#add MGI genes from "human-mouse disease connection". manually downloaded osteoporosis, bone mineral density, osteoblast clast and cyte. human and mouse genes
mgi = read.delim("../data/MGIhdpQuery_markers_20201123_124412.txt",stringsAsFactors = FALSE) #Nov-23-2020

#uses mus ensembl annotations

ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl") #uses mus ensembl annotations

mgi_mouse_id = mgi[which(mgi$Organism=="mouse"),"ID"]

gene.data <- unname(getBM(attributes=c('mgi_id', 'mgi_symbol', 'gene_biotype'),
                            filters = 'mgi_id', values = mgi_mouse_id, mart = ensembl))

mgi_mouse = gene.data[which(gene.data[,3] == "protein_coding"), 2]  
mgi_mouse = mgi_mouse[,1]


humanized_mus = convertMousetoHuman(mgi_mouse)

#
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl") #uses human ensembl annotations

mgi_hum_id = mgi[which(mgi$Organism=="human"),"ID"]


gene.data <- unname(getBM(attributes=c('entrezgene_id', 'hgnc_symbol', 'gene_biotype'),
                          filters = 'entrezgene_id', values = mgi_hum_id, mart = ensembl))

mgi_human = gene.data[which(gene.data[,3] == "protein_coding"), 2]  
mgi_human = mgi_human[,1]


mgi = c(humanized_mus, mgi_human)

mgi = tolower(mgi)
mgi = unique(mgi)

superduperset = append(mgi, unq_genes)

superduperset = na.omit(superduperset)

superduperset = unique(superduperset)
superduperset = superduperset[-which(superduperset == "")] #clean up

#add ENSEMBL ID
superduperset = as.data.frame(superduperset)
colnames(superduperset) = "SYMBOL"
superduperset$SYMBOL = toupper(superduperset$SYMBOL)
superduperset$ens = "NA"

for(i in 1:nrow(superduperset)){
  tryCatch({
    superduperset$ens[i] = mapIds(org.Hs.eg.db, keys = superduperset$SYMBOL[i], keytype = "SYMBOL", column = "ENSEMBL")
  }, error=function(e){})
}

#manual annotation for some
superduperset$ens[which(superduperset$SYMBOL == "C1ORF43")] = 'ENSG00000143612'
superduperset$SYMBOL[which(superduperset$SYMBOL == "C1ORF43")] = 'C1orf43'

superduperset$ens[which(superduperset$SYMBOL == "CIBAR1")] = 'ENSG00000188343'

superduperset$ens[which(superduperset$SYMBOL == "C19ORF81")] = 'ENSG00000235034'
superduperset$SYMBOL[which(superduperset$SYMBOL == "C19ORF81")] = 'C19orf81'

superduperset$ens[which(superduperset$SYMBOL == "C3ORF14")] = 'ENSG00000114405'
superduperset$SYMBOL[which(superduperset$SYMBOL == "C3ORF14")] = 'C3orf14'

superduperset$ens[which(superduperset$SYMBOL == "PEDS1")] = 'ENSG00000240849'

superduperset$ens[which(superduperset$SYMBOL == "C1ORF112")] = 'ENSG00000000460'
superduperset$SYMBOL[which(superduperset$SYMBOL == "C1ORF112")] = 'C1orf112'

superduperset$ens[which(superduperset$SYMBOL == "IHO1")] = 'ENSG00000173421'
superduperset$ens[which(superduperset$SYMBOL == "BPNT2")] = 'ENSG00000104331'
superduperset$ens[which(superduperset$SYMBOL == "ELAPOR2")] = 'ENSG00000164659'
superduperset$ens[which(superduperset$SYMBOL == "DNAI3")] = 'ENSG00000162643'

superduperset$ens[which(superduperset$SYMBOL == "C12ORF4")] = 'ENSG00000047621'
superduperset$SYMBOL[which(superduperset$SYMBOL == "C12ORF4")] = 'C12orf4'

write.table(superduperset, "../results/superset_humanized_proteinCoding.txt", sep = "\t", col.names = FALSE, row.names=FALSE, quote=FALSE)
