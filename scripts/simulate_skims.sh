#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -t 96:00:00
#SBATCH -J art-skims
#SBATCH --output=art-skims_%a.out
#SBATCH --error=art-skims_%a.err
#SBATCH --export=ALL
#SBATCH --array=0-5

source $HOME/.bashrc
source ${HOME}/programs/miniconda3/etc/profile.d/conda.sh

conda activate msprime-env
set -x
sample_arr=("nematode" "rotifer" "leech")
pop_size_arr=("62_500" "125_000" "250_000" "500_000" "1_000_000" "2_000_000")
inds=(0 1 2 3)
pop_size=${pop_size_arr[$SLURM_ARRAY_TASK_ID]}
samp=${sample_arr[$1]}


ART_PATH="/calab_data/mirarab/home/echarvel3/programs/art_bin_MountRainier/"

mkdir "./genome_skims/"
pushd "./genome_skims/"

#for distance in $pop_size; o
	mkdir -p ./"${samp}/${pop_size}"
	for ind in 0 1 2 3; do
		genome=$(realpath ../genomes/${samp}/${pop_size}/samp_${ind}*)
		mv ${genome%".fna"} ${genome%".fna"}".fna" 
		phred_quality=25
		read_length=150
	
		for coverage in 0.25 0.5 1 2 4; do 
			${ART_PATH}/art_illumina -q \
				-i ${genome%".fna"}".fna" \
				-l $read_length \
				-f $coverage \
				-qL ${phred_quality} -qU ${phred_quality} \
				-o ./${samp}/${pop_size}/samp-${ind}_${coverage}x &
		done
		wait	
		#find ${samp}/${pop_size}/samp-${ind}_*x.fq | xargs -I {} -P 4 gzip {}
		rm ./${samp}/${pop_size}/samp-${ind}_*x.aln
	done
#done
