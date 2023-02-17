source settings.sh
WT=13:05:00 # walltime (max for sQTLs was 7:30h and for eQTLs it was 9h)

if [[ $WT == 00:05:00 ]]; then NODE=short; else NODE=normal; fi

QTL=sqtl # whether to analyse splicing or expression 
for TISSUE in $(cat tissues.txt); do
	echo $TISSUE

	## Define the chrs to analyse
	# splitting for the sqtls cause they are more memory intensive, while analysing the eqtls all in one go
	if [[ $QTL == eqtl ]]; then
		CHROMS='all'
	else
		CHROMS=($(seq 1 22))
	fi

	## Submit jobs
	for CHR in ${CHROMS[@]}; do
		sbatch -J $TISSUE.$QTL.$CHR -o slurm-$TISSUE-$QTL-$CHR.%J -t $WT -p $NODE gtex.job $SCRIPTS/settings.sh $TISSUE $QTL $CHR
	done
done
