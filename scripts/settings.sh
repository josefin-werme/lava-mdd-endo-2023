# This file defines a bunch of bash variables that point to relevant directories and files used across the analyses. Note that the directories here are in reference to my specific setup on my cluster home dir, and would require adjustment for anyone wishing to use these scripts.

# NOTE that all files used for the analysis (except the sumstats and reference data) can be found in the data/support_files directory

# script and data dirs
LAVA="$HOME/LAVA/depression"
SCRIPTS=$LAVA/scripts
DATA=$HOME/LAVA/2020.LAVA1/Data # dir with depression sumstats
GTEX=$LAVA/gtex
GTEXv8=$GTEX/v8
LDSC=$HOME/Programs/ldsc

# reference data
REFDAT=g1000_eur.maf005
REFDIR=$HOME/Programs/g1000_eur

# locus file
LOCFILE=gencode.v26.GRCh38.protein_coding.1Mb-cis.SNP-IDs.dbsnp-g1000-subset.loci
# this locus file was created using the gencode gene locations with 1Mb windows +/- the TSS. dbSNP rsids (subsetted to g1000_eur.maf005.bim SNPs) were then added to each locus based on location
# Note that its necessary to use SNP IDs since location file was created using a GRCh38 coordinates, while g1000 are using GRCh37.

# result directories
OUT=$LAVA/out
MRI=$LAVA/mri/ukb.mri
NET=$LAVA/networks

