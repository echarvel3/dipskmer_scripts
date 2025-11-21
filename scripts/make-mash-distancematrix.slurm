#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH -t 96:00:00
#SBATCH -J mash
#SBATCH --output=mash_%a.out
#SBATCH --error=mash_%a.err
#SBATCH --export=ALL
#SBATCH --array=0-5
##SBATCH --begin=now+3hours

source $HOME/.bashrc
source ${HOME}/programs/miniconda3/etc/profile.d/conda.sh

conda activate raw_skmer
set -x
sample_arr=("nematode" "rotifer" "leech")
pop_size_arr=("62_500" "125_000" "250_000" "500_000" "1_000_000" "2_000_000")
pop_size=${pop_size_arr[$SLURM_ARRAY_TASK_ID]}
samp=${sample_arr[$1]}

mkdir genomes
genome=$(realpath /scratch03/echarvel/dipskmer/genomes/${samp}/*.hist)

mkdir -p ./libraries/${samp}/${pop_size}
mkdir -p ./distance_matrices/${samp}/${pop_size}
mkdir -p ./stats/${samp}/${pop_size}


#######  MASH   ###############
echo -en "sample\t" > ./distance_matrices/${samp}/${pop_size}-mash.txt
for x in $(ls ./libraries/${samp}/${pop_size}/*/*x.msh); do
                echo -en "${x##*/}\t"
done >> ./distance_matrices/${samp}/${pop_size}-mash.txt
echo "" >> ./distance_matrices/${samp}/${pop_size}-mash.txt

for x in $(ls ./libraries/${samp}/${pop_size}/*/*x.msh); do
        echo -en "${x##*/}\t" >> ./distance_matrices/${samp}/${pop_size}-mash.txt
        for y in $(ls ./libraries/${samp}/${pop_size}/*/*x.msh); do
                echo -en $(mash dist ${x} ${y} | cut -f3)"\t" >> ./distance_matrices/${samp}/${pop_size}-mash.txt
        done
        echo "" >> ./distance_matrices/${samp}/${pop_size}-mash.txt
done
exit
########## TRUE #############
set +x 
conda activate msprime-env
set -x
#true_count=$(python /scratch03/echarvel/dipskmer/count_diploid_snps.py)

echo -en "sample\t" > ./distance_matrices/${samp}/${pop_size}-true.txt
for x in $(ls ./genomes/${samp}/${pop_size}/*fna); do
                echo -en "${x##*/}\t"
done >> ./distance_matrices/${samp}/${pop_size}-true.txt
echo "" >> ./distance_matrices/${samp}/${pop_size}-true.txt

for x in $(ls ./genomes/${samp}/${pop_size}/*fna); do
        echo -en "${x##*/}\t" >> ./distance_matrices/${samp}/${pop_size}-true.txt
        for y in $(ls ./genomes/${samp}/${pop_size}/*fna); do
                echo -en $(python /scratch03/echarvel/dipskmer/count_diploid_snps.py ${x} ${y})"\t" >> ./distance_matrices/${samp}/${pop_size}-true.txt
        done
        echo "" >> ./distance_matrices/${samp}/${pop_size}-true.txt
done



