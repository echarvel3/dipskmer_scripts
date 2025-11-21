# README

## Overview

This repository contains a set of scripts for **simulating diploid genomes**, **simulating population-level variation**, and **producing genome skims** at different coverage levels.  
The pipeline is designed around:

- **msprime** for coalescent simulation of diploid and population genomes  
- **ART MountRainier** for generating Illumina-style sequencing reads  
- **Slurm** batch scheduling for large-scale population simulations  

These tools allow controlled evaluation of distance-based methods (e.g., Skmer, Mash, ReSkmer, DipSkmer, Respect) under different genome sizes, ploidy levels, and sequencing regimes.

---

## Main Script Summaries

### 1. `simulate_diploids.py`

A Python script that uses **msprime** to simulate diploid genomes from a reference genome and outputs genomes for each of the desired individuals.

#### Arguments

| Flag | Description |
|------|-------------|
| `-i, --input_genome` | Input reference genome (FASTA) |
| `-o, --out_dir` | Output directory for simulated diploids |
| `-p, --ploidy` | Ploidy level |
| `-s, --sample_num` | Number of individuals to simulate |
| `-z, --initial_pop_size` | Effective population size (N₀) for populations |
| `-t, --time` | Divergence time in generations |

---

### 2. `simulate_populations.slurm`

Given population size or time, simulated genomes. Calls `simulate_diploids.py` script.

Usage:

```
sbatch simulate_populations.slurm <species_index>
```

---

### 3. `simulate_skims.sh`

A Slurm job that creates **Illumina genome skims** using ART.

- Produces FASTQ reads using ART MountRainier.

Usage:

```
sbatch simulate_skims.sh <species_index>
```

---

### 4. 'make-true-distancematrix.slurm'

Calls on count_diploid_snps.py to iterate over diploid genomes and counts SNPs to obtain true distance.

---

### 5. 'make-mash-distancematrix.slurm

Iterates over samples and calls `mash dist` command to make mash distance matrix.
