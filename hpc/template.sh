#!/bin/sh
#
# Annotator script for Ginsburg 
#
#SBATCH --account=mjlab        
#SBATCH --job-name=annotator-single-run
#SBATCH --cpus-per-task 4
#SBATCH --time=01:00:00            
#SBATCH --mem-per-cpu=20G
 
module load anaconda/2-2019.10
source /burg/opt/anaconda2-2019.10/anaconda2/etc/profile.d/conda.sh
conda activate annotator

annotator --input ${INPUT} --output ${OUTPUT} --gtfdb ${GTFDB} --species hg38
