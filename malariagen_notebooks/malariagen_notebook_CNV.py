## Example MalariaGEN analysis: CNV https://malariagen.github.io/vector-data/ag3/examples/cnv-explore.html
## WORKSHOP 2 - MODULE 3
## activate conda environment in the terminal below with installed python packages
## then use here for import

# %%

import malariagen_data
from bisect import bisect_left, bisect_right
import numpy as np
import dask.array as da
from dask.diagnostics import ProgressBar
import bokeh.io as bkio
import pandas as pd
import matplotlib.pyplot as plt
%matplotlib inline
import seaborn as sns

# %% to make interactive plots with bokeh we need to set up plotting with boken 
bkio.output_notebook()

# %% setup access to Ag3 data in Google Cloud
ag3 = malariagen_data.Ag3()
ag3

# %%
df_sample_sets = ag3.sample_sets(release="3.0")
df_sample_sets

# %% load the CNV HMM data
# this includes
# The raw coverage calls (“call_RawCov” - the number of reads aligning in each 300bp window along the chromosome).
# Normalised coverage calls (“call_NormCov” - the raw coverage divided by expected coverage and multiplied by 2, such unduplicated regions have the expected diploid copy number value of 2).
# Copy number calls, which is the output of the HMM (“call_CN” - the estimated copy number in each 300 bp window).
# Each variant is a 300bp window

hmm = ag3.cnv_hmm(region="2R", sample_sets="AG1000G-GM-A").set_index(samples="sample_id")
hmm

# %% Example  - look at coverage data in these selected samples
# AG1000G-GM-C
sample_names = ['AG0004-C','AG0018-C','AG0019-C','AG0043-C','AG0049-C',
                'AG0058-C','AG0060-C','AG0064-C','AG0074-C','AG0284-C',
                'AG0291-C','AG0303-C','AG0318-C','AG0328-C','AG0340-C',
                'AG0354-C']

# %% AG1000G-GM-A
sample_names = ['AG0135-C','AG0225-C','AG0226-C','AG0234-C']

# %% AG1000G-GM-B
sample_names = ['AG0442-CW']

# %% AG1000G-GW
sample_names = ['AJ0023-C', 'AJ0037-C', 'AJ0038-C', 'AJ0039-C', 'AJ0049-C',
                'AJ0056-C','AJ0059-C','AJ0061-C','AJ0071-C','AJ0072-C','AJ0085-C',
                'AJ0097-C','AJ0107-C','AJ0113-C','AJ0117-C','AJ0120-C','AJ0121-C',
                'AJ0122-C','AJ0123-C','AJ0124-C','AJ0125-C','AJ0126-C','AJ0127-C',
                'AJ0131-C','AJ0133-C','AJ0138-C','AJ0140-C','AJ0145-C','AJ0155-C']

# %% We can plot coverage in the Cyp6aa / Cyp6p cluster of genes, which have been shown to be involved in 
# insecticide resistance. The plots below show this region of the genome, which is on chromosome 2R between
#  positions 28,460,000–28,580,000, with genomic position on the X axis. The black markers show normalised coverage,
#  while the blue line shows the estimated copy number (i.e., the output of the HMM) at each 300 bp window

# %% for Cyp6aa - Cyp6p region

for s in sample_names:
    ag3.plot_cnv_hmm_coverage(s, sample_set='AG1000G-GW', region='2R:28,460,000-28,580,000')

# %% for Cyp6m - Cyp6z region

for s in sample_names:
    ag3.plot_cnv_hmm_coverage(s, sample_set='AG1000G-GW', region='3R:6,900,000-7,000,000')

# %% there is a function in the malariagen package to visualise th HMM output for many samples at once as a heatmap

ag3.plot_cnv_hmm_heatmap(
    region='3R:6,900,000-7,000,000',
    sample_sets='AG1000G-GM-C', 
    row_height=5
);

# %% Load the CNV discordant read calls. These provide allel-specific CNV genotypes at a selection of genomic loci
# known to be involved in insecticide resistance
discordant_read_calls = (
    ag3
    .cnv_discordant_read_calls(contig="2R", sample_sets="AG1000G-GW")
    .set_index(samples="sample_id", variants="variant_id")
)
discordant_read_calls
# %% Pull out discordant read calls in the Cyp6aa/Cyp6p region 
# Locate the CNV alleles from the Cyp6aa / Cyp6p region
loc_cyp6aap = discordant_read_calls["variant_Region"].values == 'Cyp6aa/Cyp6p'
# Select data for our genes and samples of interest
discordant_read_calls_cyp6aap = (
    discordant_read_calls
    .isel(variants=loc_cyp6aap)
    .sel(samples=sample_names)
)
# Convert genotypes to a pandas DataFrame
df_cyp6aap = discordant_read_calls_cyp6aap["call_genotype"].to_pandas()
# Remove CNVs that are absent in all samples
df_cyp6aap = df_cyp6aap.loc[(df_cyp6aap != 0).any(axis=1)]
df_cyp6aap
df_cyp6aap.to_csv('AG1000G-GW_cyp6aap_discordant_reads.csv')

########################################################################
# %% Load the CNV discordant read calls for Cyp6m/Cyp6z region
discordant_read_calls = (
    ag3
    .cnv_discordant_read_calls(contig="3R", sample_sets="AG1000G-GW")
    .set_index(samples="sample_id", variants="variant_id")
)
discordant_read_calls
# %% Pull out discordant read calls in the Cyp6aa/Cyp6p region 
# Locate the CNV alleles
loc_cyp6mz = discordant_read_calls["variant_Region"].values == 'Cyp6m/Cyp6z'
# Select data for our genes and samples of interest
discordant_read_calls_cyp6mz = (
    discordant_read_calls
    .isel(variants=loc_cyp6mz)
    .sel(samples=sample_names)
)
# Convert genotypes to a pandas DataFrame
df_cyp6mz = discordant_read_calls_cyp6mz["call_genotype"].to_pandas()
# Remove CNVs that are absent in all samples
df_cyp6mz = df_cyp6mz.loc[(df_cyp6mz != 0).any(axis=1)]
df_cyp6mz
df_cyp6mz.to_csv('AG1000G-GW_cyp6mz_discordant_reads.csv')

########################################################################
# %% Load the CNV discordant read calls for Gstue region
discordant_read_calls = (
    ag3
    .cnv_discordant_read_calls(contig="3R", sample_sets="AG1000G-GW")
    .set_index(samples="sample_id", variants="variant_id")
)
discordant_read_calls
# %% Pull out discordant read calls in the Cyp6aa/Cyp6p region 
# Locate the CNV alleles
loc_gstue = discordant_read_calls["variant_Region"].values == 'Gstu/Gste'
# Select data for our genes and samples of interest
discordant_read_calls_gstue = (
    discordant_read_calls
    .isel(variants=loc_gstue)
    .sel(samples=sample_names)
)
# Convert genotypes to a pandas DataFrame
df_gstue = discordant_read_calls_gstue["call_genotype"].to_pandas()
# Remove CNVs that are absent in all samples
df_gstue = df_gstue.loc[(df_gstue != 0).any(axis=1)]
df_gstue

df_gstue.to_csv('AG1000G-GW_gstue_discordant_reads.csv')

########################################################################
# %% Load the CNV discordant read calls for Cyp9k1 region
discordant_read_calls = (
    ag3
    .cnv_discordant_read_calls(contig="X", sample_sets="AG1000G-GW")
    .set_index(samples="sample_id", variants="variant_id")
)
discordant_read_calls
# %% Pull out discordant read calls in the Cyp6aa/Cyp6p region 
# Locate the CNV alleles
loc_cyp9k1 = discordant_read_calls["variant_Region"].values == 'Cyp9k1'
# Select data for our genes and samples of interest
discordant_read_calls_cyp9k1 = (
    discordant_read_calls
    .isel(variants=loc_cyp9k1)
    .sel(samples=sample_names)
)
# Convert genotypes to a pandas DataFrame
df_cyp9k1 = discordant_read_calls_cyp9k1["call_genotype"].to_pandas()
# Remove CNVs that are absent in all samples
df_cyp9k1 = df_cyp9k1.loc[(df_cyp9k1 != 0).any(axis=1)]
df_cyp9k1
df_cyp9k1.to_csv('AG1000G-GW_cyp9k1_discordant_reads.csv')

# %% Check where each of these CNVs start and end, and see if this corresponds to the coverage increases 
# that we saw in the HMM data

dup_positions = (
    discordant_read_calls
    .reset_coords()
    [["variant_position", "variant_end"]]
    .sel(variants=["Cyp6aap_Dup16"])
    .to_dataframe()
)
dup_positions

# %% Load the per-gene modal copy number data
# We can also look at individual genes and ask whether copy number was increased in each of those genes. This table shows the copy number observed for each gene in the Cyp6aa and Cyp6p cluster, with 2 being the normal diploid copy number (no CNV)

# %% Access per-gene modal copy number data for all genes in chromosome arm 2R 
gene_copy_number = (
    ag3.gene_cnv(region="2R", sample_sets="AG1000G-GM-C")
    .set_index(genes="gene_id", samples="sample_id")
)

# %% Define genes of interest
cyp6aap_genes = dict(AGAP002862='Cyp6aa1',
                     AGAP013128='Cyp6aa2',
                     AGAP002868='Cyp6p1',
                     AGAP002869='Cyp6p2',
                     AGAP002865='Cyp6p3',
                     AGAP002867='Cyp6p4',
                     AGAP002866='Cyp6p5')

gste_genes = dict(AGAP009195='Gste1',
                  AGAP009194='Gste2',
                  AGAP009197='Gste3',
                  AGAP009193='Gste4',
                  AGAP009192='Gste5',
                  AGAP009191='Gste6',
                  AGAP009196='Gste7')

# %% Select data for genes of interest, excluding samples with poor quality HMM data
cyp6aap_gene_copy_number = (
    gene_copy_number["CN_mode"]
    .sel(genes=list(cyp6aap_genes.keys()))
    .where(~gene_copy_number.sample_is_high_variance, drop=True)
    .transpose()
    .to_pandas()
    .rename(cyp6aap_genes, axis='columns')
)

cyp6aap_gene_copy_number

# %% Select data for genes of interest, excluding samples with poor quality HMM data
gste_gene_copy_number = (
    gene_copy_number["CN_mode"]
    .sel(genes=list(gste_genes.keys()))
    .where(~gene_copy_number.sample_is_high_variance, drop=True)
    .transpose()
    .to_pandas()
    .rename(cyp6aap_genes, axis='columns')
)

gste_gene_copy_number

# %% Plot this as a heatmap for easier viewing
def plot_gene_cnv_heatmap(data):

    # create a figure
    figsize = data.shape[1], data.shape[0]*0.7
    fig, ax = plt.subplots(figsize=figsize)

    # plot a heatmap
    sns.heatmap(
        data, 
        annot=True, 
        cbar=False, 
        cmap='bwr', 
        center=2, 
        linewidths=0.1,
        ax=ax
    )

    # tidy plot
    ax.tick_params(
        axis='both', 
        which='major', 
        labelsize=10, 
        labelbottom=False, 
        bottom=False, 
        top=False, 
        labeltop=True, 
        left=False
    )
    ax.xaxis.set_label_position('top')

# %% 
plot_gene_cnv_heatmap(cyp6aap_gene_copy_number.iloc[:25])

# %%
plot_gene_cnv_heatmap(gste_gene_copy_number.iloc[:25])

# %% WORKSHOP 2 - MODULE 4
# Summarising copy number variation by gene, by calculating the modal copy number for windows overlapping a gene
# compute gene CNV frequencies using the gene_cnv_frequencies() function

# %%
ag3.gene_cnv_frequencies?

#%% analyse all genes in te Cyp6aa/p cluster on chromosome arm 2R
cyp6aap_region = "2R:28,480,000-28,510,000"
cyp9k1_region = "X:15,240,000-15,250,000"

# %%
cohorts = "admin1_year"
# %% have a look at the sample sets in Ag3
ag3.sample_sets()

# %% select mosquitoes from Burkina Faso
sample_sets = ["AG1000G-BF-A", "AG1000G-BF-B", "AG1000G-BF-C"]

# %% compute gene CNV frequencies. This outputs a pandas DF, where rows represent CNVs. 
# One row represents amplification and the other represents a deletion in copy number
burkina_cyp6aap_cnv_freqs_df = ag3.gene_cnv_frequencies(
    region=cyp6aap_region,
    cohorts=cohorts,
    sample_sets=sample_sets,
    drop_invariant=False,
)
burkina_cyp6aap_cnv_freqs_df

# %%
burkina_cyp9k1_cnv_freqs_df = ag3.gene_cnv_frequencies(
    region=cyp9k1_region,
    cohorts=cohorts,
    sample_sets=sample_sets,
    drop_invariant=False,
)
burkina_cyp9k1_cnv_freqs_df

# %% there are 11 annotated genes on the chr arm that we selected, so there are 2 times 11 = 22 rows in this dataframe
len(burkina_cyp6aap_cnv_freqs_df)

# %% look at the columns that provide info about the genes
burkina_cyp6aap_cnv_freqs_df[["gene_strand", "gene_description", "contig", "start", "end"]]

# %%
burkina_cyp9k1_cnv_freqs_df[["gene_strand", "gene_description", "contig", "start", "end"]]

# %% select cohort frequency columns
frequency_columns = [
    col for col in burkina_cyp6aap_cnv_freqs_df.columns 
    if col.startswith("frq_")
]
frequency_columns

# %% select cohort frequency columns
frequency_columns = [
    col for col in burkina_cyp9k1_cnv_freqs_df.columns 
    if col.startswith("frq_")
]
frequency_columns

# %% look at the values of these frequency columsn in our DF
# frequency values range between 0 and 1. 0 means no individuals have the variant. 
# 1 means all individuals have the variant. 0.5 means 50% of individuals have the variant.
# The max_af column simply shows the highest frequency found in the cohorts and is useful for filtering the output 
# down to just CNVs at appreciable frequency in at least one cohort, 
# and the windows columns shows how many 300 bp windows the frequency is calculated over
burkina_cyp6aap_cnv_freqs_df[frequency_columns + ["max_af", "windows"]]

# %% cyp9k1
burkina_cyp9k1_cnv_freqs_df[frequency_columns + ["max_af", "windows"]]

# %% Plotting frequency heatmaps
# also filter dataframe by max_af to just keep CNS that are present in at least one cohort over 5% using a pandas query
ag3.plot_frequencies_heatmap(
    burkina_cyp6aap_cnv_freqs_df.query("max_af > 0.05"), 
    title="Gene CNV frequencies, Burkina Faso, Cyp6aa/p locus"
)


# %% plot for cyp9k1
ag3.plot_frequencies_heatmap(
    burkina_cyp9k1_cnv_freqs_df.query("max_af > 0.05"), 
    title="Gene CNV frequencies, Burkina Faso, Cyp6aa/p locus"
)
