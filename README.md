# Enrichment-test
Enrichment test based on PCA loadings

Required packages:

allez; in bash run R CMD INSTALL allez.tar.gz

.db pkgs; in R run

source("http://bioconductor.org/biocLite.R")

biocLite("org.Hs.eg.db”)

biocLite("GO.db")


Example commands to run the script:
- Rscript Enrich.R ExampleLoading.csv MarkerLists.csv
- Rscript Enrich.R ExampleLoading.txt MarkerLists.txt
- Rscript Enrich.R ExampleLoading.csv 

The 3rd term indicates the name of the loadings (or other scores). 
The file should contain only one column. Row names should be gene names. Each entry represents loading (or score) 
of that particular gene.
Currently the program takes csv files or tab delimited file.
The input file will be treated as a tab delimited file if the suffix is not '.csv'.

The 4th term defines the marker lists of interest. It could be csv or tab delimited file. Each row represents a marker list. 
Row names are the list names. The number of columns should be N, in which N is the number of genes in the longest list. 
For lists with length M - N, the M+1, ..., N's column in that row should be filled with "" or " ". If the 4th term
is not specified, only GO terms will be considered.

Outputs:

XX_allsets.txt: enrichment results of all gene sets (GO sets + local sets). 

XX_localsets.txt: enrichment results of only local sets. The summary statistics of these sets are the same as those in the XX_allsets.txt

The output files contain GO term (NA for local sets), set p value, adjusted p value, z score, set size, set mean and set sd. Sets are sorted by p value. Sets with large absoloute z scores are expected to have small p values.


## EACI test
Install EACI package:

in bash

run R CMD INSTALL EACI_0.0.1.tar.gz



Example commands to run the script:
- Rscript Enrich_eaci.R ExampleLoading.csv MarkerLists.csv
- Rscript Enrich_eaci.R ExampleLoading.txt MarkerLists.txt
- Rscript Enrich_eaci.R ExampleLoading.csv 
 
Outputs

XX_EACIenrichment_allsets.txt

XX_EACIenrichment_localsets.txt

