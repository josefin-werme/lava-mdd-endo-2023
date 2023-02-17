This repository contains the results and analysis scripts from the article [paper title], [preprint link]. In this article, we used LAVA to perform local genetic correlation analyses between depression and endophenotypes across different levels of function (gene expression & splicing, regional brain volume, and brain network connectivity). 


### Results
The output from all analyses can be found in the **'results'** directory, with the following subdirectories:
- **local_rg**: further subsetted by **'eqtl'**, **'sqtl'**, **'regions'**, and **'network'** directories, containing the reults from all univariate (**\*.univ**) and bivariate (**\*.bivar**) lava results for the respective class of endophenotype
- **enrichment**: endophenotype specific enrichment, based on the bivariate LAVA results



### LAVA Analysis Scripts
The main LAVA analyses were performed on a slurm cluster. The scripts can be found in the **'scripts'** directory, and are divided into three main categories: 

1. Job submission scripts (**\*.sh**) - iterate through all endophenotypes (listed in the **'scripts/\*.txt'** files) and submit the relevant job scripts
2. Job scripts (***.job**) - for any given endophenotype, copy all the data to temp dir and run the local rg analyses
3. Analysis scripts (***.R**) - perform the local rg analyses

Since each category of endophenotype requires slightly different analysis setup, there are separate submission, job, and analysis scripts for these, indicated by different prefixes:
- **gtex.*** - Tissue specific gene expression & splicing
- **regions.*** - Regional brain volume
- **networks.*** - Functional and structural brain network connectivity

There is also a **settings.sh** file which specifies useful paths and files required by *sh and *job scripts, as well as the text files that list the relevant tissues, regions, and networks (note that the setup here is in reference to how I organised things on my account, and this would need to be adjusted for anyone wishing to run these scripts themselves) 



### Enrichment Analysis Script
Based on the results from the main LAVA analyses, endophenotype specific enrichment analyses were performed locally using this script: **'functional_enrichment.R'**. NOTE: to run this script, you'll need the output from the LAVA bivariate analyses in **'results/local_rg'**



### Analysis Files
Support files needed to run the analyses can be found in the **'data/analysis_files'** directory, and include the:

- Input info files (input.info.*)
- Sample overlap files (sample.overlap.*)
- Gene definition and locus file (gencode.v26.GRCh38.protein_coding.*)
- MSigDB gene set definition file (magma_GS_v6.2.txt)

