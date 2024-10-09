
# %% PBS
# PBS uses FST to identify genomic regions showing greater evolutionary change in one group 
# (here, the resistant samples) 
# relative to a closely related group (susceptible samples) and an outgroup. While originally designed to detect
# positive selection, it has also been used to detect phenotypic association (Grau-Bové et al., 2021).
# Note, For both H12 and PBS, phenotype permutations were performed as for FST to filter out false positives
# caused by the presence of extended swept haplotypes.
# calculate in 1000 bp windows and plot against the genome. What is classed as significant?
# can I use the control samples as ac3 here?

# %% create ac3 from control samples - use a separate population
# select samples
#con_samples = df_samples[df_samples['phenotype'] == 'control'].index.values
#con_samples

# select genotypes for samples
#gt_con_samples = gt.take(con_samples, axis=1)
#gt_con_samples

# create allele counts array
#ac_con_samples = gt_con_samples.count_alleles()
#ac_con_samples

# %% compute PBS

# allel.pbs(ac_sus_samples, ac_res_samples, ac_con_samples, 1000)
