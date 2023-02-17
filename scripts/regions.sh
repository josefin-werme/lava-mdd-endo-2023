# this script will submit the regional brain volume analyses, separately for each region
# note that homologous regions are analysed within the same job since they will be meta-analysed in loci where the rg between them is high

source settings.sh
WT=13:05:00	# takes about 1h ; was set to 3, but now increasing load x 4 as we're not filtering on univ significance in depression a priori
if [[ $WT == 00:05:00 ]]; then NODE=short; else NODE=normal; fi

for REGION in $(cat regions.txt); do
	echo $REGION
	sbatch -J $REGION.mri -o slurm-$REGION.%J -t $WT -p $NODE regions.job $REGION $SCRIPTS/settings.sh
done

