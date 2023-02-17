source settings.sh
WT=8:00:00 # about 2h with univ filtering, about 5.5h with all genes
if [[ $WT == 00:05:00 ]]; then NODE=short; else NODE=normal; fi

for NETWORK in $(cat networks.txt); do 
	sbatch -J $NETWORK.net -o slurm-$NETWORK-%J -t $WT -p $NODE networks.job $NETWORK $SCRIPTS/settings.sh
done

