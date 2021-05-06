# BMD TWAS/colocalization project
This repository contains analysis code for the BMD TWAS/colocalization project. Preprint coming soon.

### 1) MetaXcan/S-Multixcan
TWAS analyses for eBMD and fracture GWAS were performed by Will Rosenow. These can be found [here](https://github.com/Farber-Lab/BMD-MetaXcan).

### 2) fastENLOC colocalization
The following steps were perfomed for fastENLOC colocalization:
  * eBMD GWAS and fracture GWAS summary statistics were obtained [here](http://www.gefos.org/?q=content/data-release-2018).
  * GWAS coordinates were converted to hg38 coordinates using [UCSC's online liftOver tool](https://genome.ucsc.edu/cgi-bin/hgLiftOver). The [prep_enloc.R](./src/prep_enloc.R) script was used to prepare input and parse the output. This script was also used to prepare the input for fastENLOC, where z-scores were calculated (beta/SE) and loci were defined based on the results of Berisa and Pickrell, 2015 ([PMC4731402](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4731402/)). 
  * The resulting files were gzipped and z-scores were then converted to PIPs using torus, using the following commands in Bash:
  ```bash
  #BMD
  
  gzip BMD_zscore
  
  torus -d BMD_zscore.gz --load_zval -dump_pip BMD.gwas.pip
  
  gzip BMD.gwas.pip
  
  #Fracture
  
  gzip frax_zscore
  
  torus -d frax_zscore.gz --load_zval -dump_pip frax.gwas.pip
  
  gzip frax.gwas.pip
```

  * fastENLOC colocalization was performed using pre-computed GTEx multi-tissue eQTL annotations with hg38 position ID, obtained from [here](https://github.com/xqwen/fastenloc). The following commands (in Bash) were used to perform the colocalizations:
  ```bash
  
  #insert relevant tissue name here, do for all 49 GTEx tissues. Alternatively, run a script to do this for you.
  
  tissue=tissue_name
  
  #colocalization
  
  #BMD
  
  fastenloc -eqtl ../gtex_v8.eqtl_annot.vcf.gz -gwas ../BMD.gwas.pip.gz -total_variants 14000000 -t $tissue -prefix ../enloc_out/${tissue}
  
  #Fracture
  
  fastenloc -eqtl ../gtex_v8.eqtl_annot.vcf.gz -gwas ../frax.gwas.pip.gz -total_variants 14000000 -t $tissue -prefix ../enloc_frax_out/${tissue}
```
  
### 3) Construction of "known bone gene" superset
  * The "known bone gene" list was constructed using the [make_superset_humanized.R](./src/make_superset_humanized.R) script.
  
### 4) IMPC results
  * Obtaining IMPCs results for BMD, as well as performing (and plotting) the statistical analyses of *Gpatch1* knockdowns was performed in the [impc_SOLR.R](./src/impc_SOLR.R) script.
  
### 5) **Downstream analyses**
  * Downstream analyses integrating TWAS and colocalization, as well as many plots in our publication, were preformed in the [ppp6r3_paper_analysis_new.R](./src/ppp6r3_paper_analysis_new.R) script. 
  * The Estrada *et al.* GWAS results were used here. Summary statistics can be obtained [here](http://www.gefos.org/?q=content/data-release-2012), and the [get_estrada_gwas_hg38_pos.R](./src/get_estrada_gwas_hg38_pos.R) was used to convert GWAS SNPs to hg38 coordinates.

### 6)***Ppp6r3* functional validation***
  * Statistical analyses and plotting of *Ppp6r3* experimental mice was performed in the [ppp6r3_analysis.R](./src/ppp6r3_analysis.R) script.
  
  
---  
  
For more information, please contact Basel Al-Barghouthi (bma8ne AT virginia DOT edu) or Charles Farber (crf2s AT virginia DOT edu).
