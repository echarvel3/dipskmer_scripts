#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=7
#SBATCH -t 96:00:00
#SBATCH -J simulate-fst
#SBATCH --output=simulate_%a.out
#SBATCH --error=simulate_%a.err
#SBATCH --export=ALL
#SBATCH --array=0-5

source $HOME/.bashrc
source ${HOME}/programs/miniconda3/etc/profile.d/conda.sh

conda activate msprime-env
set -x
sample_arr=("nematode" "rotifer" "leech")
pop_size_arr=("62_500" "125_000" "250_000" "500_000" "1_000_000" "2_000_000")
pop_size=${pop_size_arr[$SLURM_ARRAY_TASK_ID]}
samp=${sample_arr[$1]}

mkdir genomes
genome=$(realpath /scratch03/echarvel/dipskmer/genomes/${samp}/*.fna)
mkdir -p ./genomes/${samp}/${pop_size}
x=$(ls ./genomes/${samp}/${pop_size} | wc -l)
if [ $x -eq 0 ]; then
	python /scratch03/echarvel/dipskmer/simulate_diploids.py -i ${genome} -o ./genomes/${samp}/${pop_size}/ -p 2 -s 2 -z 500_000 -t ${pop_size}
fi

