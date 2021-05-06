library(solrium)
library(jsonlite)
library(OpenStats)
library(ggplot2)
#packageVersion("Solrium") version: 1.1.4
setwd("C:/projects/GWAS_project/src")





#open conn to IMPC "experiment" core
(impc = SolrClient$new(host = "www.ebi.ac.uk",path="/mi/impc/solr/experiment/select", scheme = "https",errors = "complete",port=NULL))




#query core for all "experimental" mice with BMD parameter name, and fact by gene symbol.  faceting gives counts of occurence of gene symbol, with a mincount of 14, which is the minimum (7 males + 7 females) for DXA phenotyping protocol
(res = impc$facet(params = list(q="parameter_name:Bone*Mineral*Density*", facet.field='gene_symbol',facet.limit=-1, facet.sort="count", facet.mincount=14,rows=0),progress = httr::progress()))

#(res = impc$search(params = list(q="parameter_name:Bone*Mineral*Density* AND biological_sample_group:experimental", rows=-1),progress = httr::progress()))

#(ctrl = impc$search(params = list(q="parameter_name:Bone*Mineral*Density* AND biological_sample_group:control"),progress = httr::progress()))


#these are the genes with experimental mice for BMD query term
impc_genes = (res$facet_fields$gene_symbol)
#vector of the genes
impc_genes_vec = unique(as.character(impc_genes$term))

#remove gene names with parentheses. There are two of these and were causing issues with getting the results below.
#"Gt(ROSA)26Sor" "Tg(Thy1-MAPT*P301S)2541Godt"

impc_genes_vec = impc_genes_vec[-grep(pattern = "\\(",impc_genes_vec)] 


#For each of the genes above, get stats report.
#conn to "statistical-result" core
(stat = SolrClient$new(host = "www.ebi.ac.uk",path="/mi/impc/solr/statistical-result/select", scheme = "https",errors = "complete",port=NULL))
#

#example for search function. I chose this because it has 85 columns, which is the largest number.
#Apparently the genes vary in columns available, so i chose this to encompass (hopefully) the largest colnames for matching in next loop
q = paste0("parameter_name:Bone*Mineral*Density* AND marker_symbol:",as.character(impc_genes_vec[2]))
res_stat = stat$search(params = list(q=q, rows=-1), progress = httr::progress())
#

#get stats output for each gene and put in a dataframe
out = as.data.frame(matrix(ncol=84,nrow=0))
colnames(out) = colnames(res_stat)

#j=j2=1
for(i in 1:length(impc_genes_vec)){
 geneName =  impc_genes_vec[i]
 q = paste0("parameter_name:Bone*Mineral*Density* AND marker_symbol:",as.character(geneName))
 res_stat = as.data.frame(stat$search(params = list(q=q, rows=-1)))
 
 out = merge(out, res_stat, all.x = T, all.y = T)
 #print(ncol(res_stat))
 #j2 = j+(nrow(res_stat)-1)
 #out[j:j2,] = res_stat
 #j = j2+1 
 
 print(i)
}

#last accessed Mar. 8, 2021
#save(out, file = "../results/impc_results_raw")
load(file = "../results/impc_results_raw")

#keep only "successful" analyses
out = out[which(out$status == "Successful"),]
#some are just "not processed".

#remove BMD excluding skull
table(out$parameter_name)
out = out[-grep(pattern = "including skull", x = out$parameter_name),]

#NOTE: SPTBN1 is NOT PROCESSED. WHEN LOOKING, THERE ARE NO CONTROL MICE AND NO WEIGHTS FOR THE EXPERIMENTAL MICE
save(out, file = "../results/impc_results_raw_processed_noskull")
#genes with significant terms at 5e-2 pval
out_5e2 = out[which(as.numeric(out$male_ko_effect_p_value)<=0.05 | as.numeric(out$female_ko_effect_p_value)<=0.05 | as.numeric(out$genotype_effect_p_value)<= 0.05),]
saveRDS(out, file="../results/impc_out_5e2.rds")

### get gpatch1 raw data for plotting ###
#Gpatch1 doc id : 64d6ee80a46520cf580c5588b588b0ba
(stat = SolrClient$new(host = "www.ebi.ac.uk",path="/mi/impc/solr/statistical-raw-data/select", scheme = "https",errors = "complete",port=NULL))




gpatch1 = stat$search(params = list(q="doc_id:64d6ee80a46520cf580c5588b588b0ba", rows=-1),progress = httr::progress())

decode_conn = rawConnection(base64enc::base64decode(as.character(gpatch1[,"raw_data"])))

#decoded=stream_in(gzcon(decode_conn))
decoded[is.na(decoded)] = "NA"
#decoded = decoded[-which(decoded$body_weight == "NA"),]
closeAllConnections()

obj = OpenStatsList(decoded,dataset.colname.sex = "specimen_sex", dataset.colname.weight = "body_weight", debug = F ,dataset.values.missingValue = "NA")

#analyze
analysis_out = OpenStatsAnalysis(OpenStatsListObject = obj, method = "MM",MM_BodyWeightIncluded = TRUE,debug = F)

#wrap
analysis_out = OpenStatsReport(analysis_out)
#get output from wrapper

geno_pval = analysis_out$`Genotype p-value`
male_pval = analysis_out$`Sex MvKO p-value`
female_pval = analysis_out$`Sex FvKO p-value`
analysis_method = analysis_out$`Applied method`






analysis_out$`Additional information`$Data$`Summary statistics`

t1 = obj@datasetPL

ggplot(t1, aes(x=Sex, y=data_point, fill=factor(Genotype))) +
  geom_boxplot() + 
  labs(fill="Genotype") +
  geom_point(position = position_jitterdodge(), alpha=0.3) +
  theme_bw(base_size = 16)

