
import msprime
import tskit
import numpy as np
from screed import fasta
import argparse
import sys
#import pandas as pd

def simulate_ancestry(initial_size=500_000, 
                   mutation_rate=.3*10**-8, 
                   n_samples=2,
                   recombination_rate=.3*10**-8,
                   sequence_length=100_000,
                   ploidy=2,
                   random_seed=None
                   ):
    '''Divergence is a product of population size and mutation rate product. 
    Since mutation rate is published data, change pop size to change divergence. '''

    demography = msprime.Demography()
    demography.add_population(name="A", initial_size=initial_size) 
    ts = msprime.sim_ancestry(
    samples=[msprime.SampleSet(n_samples, population="A", ploidy=ploidy)], 
    demography=demography,
    recombination_rate=recombination_rate, 
    sequence_length=sequence_length,
    ploidy=ploidy,
    random_seed=random_seed)

    mts = msprime.sim_mutations(ts,
                             rate=mutation_rate,
                             model="JC69",
                             random_seed=random_seed)
    return(mts)

def apply_variants_to_ref(reference, mts):

    #reference = reference.upper()
    if len([m for m in mts.variants()]) == 0:
        return [reference for x in range(0, len(mts.samples()))]
    
    samples = []

    for s in range(0, len(mts.samples())):
        # print("sample:", s)
        vars = [v for v in mts.variants(samples=[s])]

        var_count=0
        mutated_genome=""
        for nt, pos in zip(reference, range(0, len(reference))):
            #print(nt, pos, var_count, vars[var_count].site.position)
            nucleotide = reference[pos]

            if vars[var_count].site.position == pos:
                v = vars[var_count]
                nucleotide = v.alleles[v.genotypes[0]]

                if (reference[pos] == nucleotide):
                    nucleotide = v.site.ancestral_state

                if (reference[pos] == v.alleles[1]):
                    #if reference is the same as the mutated...
                    if v.genotypes[0] == 1:
                    # and there is a mutation then change to ancestral
                        nucleotide = v.site.ancestral_state
                    else:
                    # if there is no mutation, keep the same character (in this case the mutated one)
                        nucleotide = v.alleles[1]

                print("is mutated?", v.genotypes[0], "original", reference[pos],  
                      "alleles", v.alleles, "mutated", v.alleles[v.genotypes[0]], "ancestral", v.site.ancestral_state, 
                      "final nuc", nucleotide)
                var_count += 1
                if var_count == len(vars):
                    mutated_genome = mutated_genome + nucleotide + reference[pos+1:]
                    break
            mutated_genome=mutated_genome+nucleotide

        samples.append(mutated_genome)
    return(samples)

def simulate_fixation_index(initial_size=500_000, 
                   mutation_rate=0.3*10**-8, 
                   n_samples=2,
                   recombination_rate=.3*10**-8,
                   sequence_length=100_000,
                   ploidy=2,
                   random_seed=None,
                   time=None
                   ):
    '''Divergence is a product of population size and mutation rate product. 
    Since mutation rate is published data, change pop size to change divergence. '''
    if time is None:
        time = mutation_rate

    demography = msprime.Demography()
    demography.add_population(name="A", initial_size=initial_size)
    #print("A")
    demography.add_population(name="B", initial_size=initial_size)
    #print("B")
    demography.add_population(name="C", initial_size=initial_size)
    #print("C")
    demography.add_population_split(time=time, derived=["B", "C"], ancestral="A")
    #print("SPLIT")
    ts = msprime.sim_ancestry(
    samples=[msprime.SampleSet(n_samples, population="B", ploidy=ploidy),
             msprime.SampleSet(n_samples, population="C", ploidy=ploidy)], 
    demography=demography,
    recombination_rate=recombination_rate, 
    sequence_length=sequence_length,
    ploidy=ploidy,
    random_seed=random_seed)
    #print("ANCESTRY")
    # print([x for x in ts.populations()])
    # print(tuple(int(x) for x in ts.samples(population = 1)))

    # print(ts.diversity(sample_sets = [0,1]))

    sample_set = []
    for x in ts.individuals():
        sample_set.append([int(n) for n in x.nodes])
    #print(sample_set)

    mts = msprime.sim_mutations(ts,
                             rate=mutation_rate,
                             model="JC69",
                             random_seed=random_seed)
    
    pop1_nodes = [s for s in range(0,n_samples*ploidy)]
    pop2_nodes = [s for s in range(n_samples*ploidy, 2*n_samples*ploidy)]

    print("POP1_DIVERSITY:", mts.diversity(sample_sets =  pop1_nodes))
    print("POP2_DIVERSITY:", mts.diversity(sample_sets =  pop2_nodes))
    print("DIVERGENCE:", mts.divergence(sample_sets = [pop1_nodes,pop2_nodes]))
    #print("FST", mts.Fst(sample_sets = [pop1_nodes,pop2_nodes]))
    #print()

    return(mts)

def msprime_mutate_existing(args):
    print("hoi")
    in_path=args.input_genome
    out_path=args.out_dir
    ploidy=args.ploidy
    n_samples=args.sample_num
    initial_size = args.initial_pop_size
    time = args.time
    print(time)

    with open(in_path, "r") as fin:
        for record in fasta.fasta_iter(fin):
            if time is None:
                mts = simulate_ancestry(sequence_length=len(record.sequence.replace("N", "")), 
                                        n_samples=n_samples, 
                                        ploidy=ploidy, 
                                        initial_size=initial_size)
            else:
                mts = simulate_fixation_index(sequence_length=len(record.sequence.replace("N", "")), 
                        n_samples=n_samples, 
                        ploidy=ploidy, 
                        initial_size=initial_size,
                        time = time)

            samples = apply_variants_to_ref(record.sequence.upper().replace("N", ""), mts)

            for i in range(0, len(samples), ploidy):
                sample_file=out_path + "samp_" + str(i//ploidy)
                with open(sample_file, "a") as fout:
                    for j in range(0,ploidy):
                        fout.write(f">{record.name}.{str(j)}_mutated\n")
                        fout.write(f"{samples[i+j]}\n")

    return

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
                    prog='Genome Simulator',
                    description='Simulates Population Genomes based on Input Genome')
    
    parser.add_argument('-i', '--input_genome', type = str)   
    parser.add_argument('-o', '--out_dir', type = str) 
    parser.add_argument('-p', '--ploidy', type = int) 
    parser.add_argument('-s', '--sample_num', type = int) 
    parser.add_argument('-z', '--initial_pop_size', type = int) 
    parser.add_argument('-t', '--time', type = int, default = None) 
    
    parser.set_defaults(func=msprime_mutate_existing)
    args = parser.parse_args()
    args.func(args)

#print(apply_variants_to_ref(reference = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", 
#                            mts = simulate_fixation_index(time = 500_000, 
#                                                          sequence_length = len("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"))))
#print(apply_variants_to_ref(reference = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", 
#                            mts = simulate_ancestry(n_samples=2, 
#                                                    sequence_length=len("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
#                                                    initial_size=1_000_000)))
