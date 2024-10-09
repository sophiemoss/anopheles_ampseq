## Example MalariaGEN analysis: SNPs https://malariagen.github.io/vector-data/ag3/examples/snps.html
## activate conda environment in the terminal below with installed python packages
## then use here for import

# %%
import malariagen_data

# %%
ag3 = malariagen_data.Ag3()
ag3

seq = ag3.genome_sequence("2R").compute()
seq[:100]

#%%
from collections import Counter
from functools import lru_cache
import numpy as np
import pandas as pd
import dask
import dask.array as da
# silence some dask warnings
dask.config.set(**{'array.slicing.split_large_chunks': True})
import matplotlib as mpl
import matplotlib.pyplot as plt
%matplotlib inline
import seaborn as sns
import allel
from dask.diagnostics.progress import ProgressBar
import malariagen_data

ag3 = malariagen_data.Ag3()
ag3

def plot_sequence_composition(seq_id, window_size=100_000, ax=None):
    
    # load reference sequence
    seq = ag3.genome_sequence(seq_id).compute()
    
    if ax is None:
        # make the figure size relative to largest contig
        figw = 7 * len(seq) / len(ag3.genome_sequence("2R"))
        fig, ax = plt.subplots(figsize=(figw, 1))
 
    # convert to upper-case
    seq = np.char.upper(seq)

    # locate nucleotides
    is_a = seq == b'A'
    is_c = seq == b'C'
    is_g = seq == b'G'
    is_t = seq == b'T'
    is_n = seq == b'N'
    # check there's nothing unexpected
    is_other = ~is_a & ~is_c & ~is_g & ~is_t & ~is_n
    assert np.sum(is_other) == 0

    # construct windows
    bins = np.arange(0, len(seq), window_size)

    # count nucleotides in each window
    h_a, _ = np.histogram(np.nonzero(is_a)[0] + 1, bins=bins)
    h_c, _ = np.histogram(np.nonzero(is_c)[0] + 1, bins=bins)
    h_g, _ = np.histogram(np.nonzero(is_g)[0] + 1, bins=bins)
    h_t, _ = np.histogram(np.nonzero(is_t)[0] + 1, bins=bins)
    h_n, _ = np.histogram(np.nonzero(is_n)[0] + 1, bins=bins)

    # plot frequence of nucleotides within each bin
    left = bins[:-1]
    bottom = 0
    width = np.diff(bins)
    palette = sns.color_palette('colorblind')
    colors = [palette[i] for i in [2, 0, 3, 8]] + ['k']
    for h, c, l in zip([h_a, h_t, h_g, h_c, h_n], colors, 'ATGCN'):
        ax.bar(left, h, width=width, bottom=bottom, color=c, align='edge', label=l)
        bottom += h
        
    # tidy up plot
    ax.set_xlim(0, len(seq))
    ax.set_yticks(ax.get_ylim())
    ax.set_yticklabels(['0%', '100%'])
    ax.set_title(seq_id)
    # convert X axis to Mbp
    ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda x, p: x//1e6))
    ax.set_xlabel("Position (Mbp)")
    
    # add centromere/telomere annotation
    if seq_id in {'2L', '3L'}:
        ltxt = "centromere"
        rtxt = "telomere"
    else:
        ltxt = "telomere"
        rtxt = "centromere"
    ax.annotate(ltxt, xy=(0, 1), xycoords="axes fraction", 
                xytext=(0, 2), textcoords='offset points', va='bottom', ha='left')
    ax.annotate(rtxt, xy=(1, 1), xycoords="axes fraction", 
                xytext=(0, 2), textcoords='offset points', va='bottom', ha='right')
    
    # add legend - reverse order so matches the plot
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(reversed(handles), reversed(labels), 
              loc='center left', bbox_to_anchor=(1, .5), 
              prop=dict(family='monospace'), ncol=1,
              frameon=False)
        
plot_sequence_composition('2R')

# %%
plot_sequence_composition('2L')
# %%
plot_sequence_composition('X')

# %%
ds_snps = ag3.snp_calls(region="2R", sample_sets="3.0")
ds_snps


# %%
# variant positions here holds the 1 based positions (coordinates) where the SNPs were called
# eg. first ten positions

pos = ds_snps["variant_position"].data
pos[:10].compute()

# %%
# these positions do not increase by 1 the whole way through the genome, because
# there are gaps in the reference genome. Visualise that by plotting distance between adjacent sites:

def plot_gaps(contig):
    ds_snps = ag3.snp_calls(region=contig)
    pos = ds_snps["variant_position"].data.compute()
    d = np.diff(pos)
    loc_gap = d > 1
    x = pos[1:][loc_gap]
    y = d[loc_gap]
    fig, ax = plt.subplots(figsize=(7, 2))
    ax.plot(x, y, color='black', linewidth=.5, marker='o', linestyle=' ', markersize=2, mfc='none')
    ax.set_title("2R")
    # convert X axis to Mbp
    ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda x, p: x//1e6))
    ax.set_xlabel("Position (Mbp)")
    ax.set_ylabel('Gap size (bp)')

# %%
plot_gaps(contig="2R")

#%%
len(seq) - len(pos)

# %%
alleles = ds_snps["variant_allele"].data
alleles[:10].compute()

# %%
seq[pos[:10].compute() - 1]


# access SNP genotypes for chromosome arm 2R for all wild samples in Ag3

# %%
ds_snps = ag3.snp_calls(region="2R", sample_sets="3.0")
gt = ds_snps["call_genotype"].data
gt.shape

# %%
gt[9, 0].compute()
gt[9, 2392].compute()
gt[9, 248].compute()

#%%
# wrap the gt variable with a genotype dask array
gt = allel.GenotypeDaskArray(ds_snps["call_genotype"].data)
gt

# %% allele count for first 10 sites on chromosome arm 2
ac = gt[:10].count_alleles(max_allele=3).compute()
ac.displayall()

# %%
# how many alleles were observed
ac.allelism()

# %%
# whether any alternate alleles were observed
ac.is_variant()
# %%
# whether more than one allele was observed
ac.is_segregating()
# %%
# whether exactly two alleles were observed
ac.is_biallelic()
# %%
# total number of alleles called
an_called = ac.sum(axis=1)
an_called
# %%
# total number of missing allele calls
an_missing = (gt.shape[1] * 2) - an_called
an_missing
# %%
### SITE FILTERS
# gamb_colu_arab site filter for chromosome 2
ds_snps = ag3.snp_calls(region="2R")
loc_pass = ds_snps["variant_filter_pass_gamb_colu_arab"].data
loc_pass
# %%
# examine the first ten values
loc_pass[:10].compute()
# %%
# how many sites in total pass the filter?
loc_pass.sum().compute()
# %%
# how many sites pass in each chr arm?
def tabulate_filter_counts():
    rows = []
    for contig in ag3.contigs:
        contig_length = len(ag3.genome_sequence(contig))
        ds_snps = ag3.snp_calls(region=contig)
        n_sites = ds_snps.dims["variants"]
        row = [contig, contig_length, n_sites]
        for mask in "gamb_colu_arab", "gamb_colu", "arab":
            loc_pass = ds_snps[f"variant_filter_pass_{mask}"].data
            n_pass = loc_pass.sum().compute()
            pc_pass = "{:.1%}".format(n_pass / n_sites)
            row.append(pc_pass)
        rows.append(row)
    df = pd.DataFrame.from_records(
        rows, 
        columns=["contig", "contig length", "no. sites", "pass gamb_colu_arab", "pass gamb_colu", "pass arab"]
    )
    return df

df_filter_counts = tabulate_filter_counts()
df_filter_counts
# %% plot the fraction of genome sites that pass the site filters in moving windows over each of the chromosome arms.
def plot_site_filters(contig, window_size=200_000, ax=None):

    # setup figure and axes
    if ax is None:
        contig_length = len(ag3.genome_sequence(contig))
        figw = 7 * contig_length / len(ag3.genome_sequence("2R"))
        fig, ax = plt.subplots(figsize=(figw, 2))
        sns.despine(ax=ax, offset=5)

    # setup X variable
    ds_snps = ag3.snp_calls(region=contig)
    pos = ds_snps["variant_position"].values
    x = allel.moving_statistic(pos, statistic=np.mean, size=window_size)

    # setup Y variables and plot
    for mask in "gamb_colu", "arab", "gamb_colu_arab":
        loc_pass = ds_snps[f"variant_filter_pass_{mask}"].values
        y = allel.moving_statistic(loc_pass, statistic=np.mean, size=window_size)
        ax.plot(x, y * 100, label=mask)

    # tidy plot
    ax.set_xlim(0, contig_length)
    ax.set_ylim(0, 100)
    ax.set_title(contig)
    # convert X axis to Mbp
    ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda x, p: x//1e6))
    ax.set_xlabel("Position (Mbp)")
    ax.set_ylabel('% pass')
    ax.legend(loc='upper left', bbox_to_anchor=(1, 1), title="Filter",
              frameon=False)

# %%
plot_site_filters(contig="2R")
# %%
plot_site_filters(contig="2L")
# %%
plot_site_filters(contig="3L")
# %%
plot_site_filters(contig="X")
# %% access SNP calls with gamb_colu_arab filter already called
ds_snps_filtered = ag3.snp_calls(region="2R", site_mask="gamb_colu_arab", sample_sets="3.0")
ds_snps_filtered

# %% locating a single SNP
ds_snps = ag3.snp_calls(region="2L", sample_sets="3.0")
ds_snps = ds_snps.set_index(variants="variant_position", samples="sample_id")
ds_snps

# %%
# access the SNP data at that site
ds_site = ds_snps.sel(variants=2422652)
ds_site
# %%
# access the alternate alleles
alleles_site = ds_site["variant_allele"].values
alleles_site
# %%
# access the genotypes
gt_site = ds_site["call_genotype"].values
gt_site
# %%
# locate data for the Vgsc gene, which spans the region 2L:2358158-2431617.
ds_region = ds_snps.sel(variants=slice(2358158, 2431617))
ds_region
# %%
# obtain SNP positions within selected region
pos_region = ds_region["variants"].values
pos_region
# %%
# alleles within selected region
alleles_region = ds_region["variant_allele"].values
alleles_region
# %%
# genotypes within selected region
gt_region = allel.GenotypeDaskArray(ds_region["call_genotype"].data)
gt_region
# %%
# allele counts
ac_region = gt_region.count_alleles(max_allele=3).compute()
ac_region
# %%
# how many segregating SNPs within our region of interest?
ac_region.is_segregating().sum()

# %% LOCATING SAMPLES
# load metadata
df_samples = ag3.sample_metadata(sample_sets="3.0").set_index("sample_id")
df_samples

# %% access SNPs for whatever chromosome arm you are analysing
ds_snps = ag3.snp_calls(region="2L", sample_sets="3.0")
ds_snps = ds_snps.set_index(variants="variant_position", samples="sample_id")
ds_snps

# %% the number of rows in the sample metadata dataframe and the size of the “samples” dimension in the SNP calls dataset are the same
len(df_samples) == ds_snps.dims["samples"]

# %% Let’s extract both sample metadata and genotypes for a single sample, e.g., sample AB0110-C.
s = df_samples.loc['AB0110-C']
s

# %%
ds_sample = ds_snps.sel(samples='AB0110-C')
ds_sample
# %%
# compute some genotype counts for our sample of interest
gt_sample = allel.GenotypeVector(ds_sample["call_genotype"].values)
gt_sample
# %%
gt_sample.count_missing()
# %%
gt_sample.count_hom_ref()
# %%
gt_sample.count_het()
# %%
gt_sample.count_hom_alt()

# %% locating multiple samples
# select An. gambiae samples from Burkina Faso collected in 2012
query = "(country == 'Burkina Faso') and (year == 2012) and (taxon == 'gambiae')"
loc_samples = df_samples.eval(query).values
loc_samples

# %%
# examing their metadata
df_samples_selected = df_samples.loc[loc_samples]
df_samples_selected

# %%
# select SNP calls
ds_snps_selected = ds_snps.isel(samples=loc_samples)
ds_snps_selected

# %% combine site selections with sample selections
# obtain SNP calls at 2L:2422652 in sample AB0110-C
ds_snps.sel(variants=2422652, samples=loc_samples)

# %%
# access SNP calls in the Vgsc gene in An. gambiae samples from Burkina Faso collected in 2012
ds_snps_vgsc = ds_snps.sel(variants=slice(2358158, 2431617), samples=loc_samples)
ds_snps_vgsc

# %%
# how many segregating sites now?
gt_vgsc = allel.GenotypeDaskArray(ds_snps_vgsc["call_genotype"].data)
ac_vgsc = gt_vgsc.count_alleles(max_allele=3).compute()
ac_vgsc.is_segregating().sum()

#############
# %% scikit allel may require a genome accessibility mask
is_accessible = ag3.is_accessible(region="2R", site_mask="gamb_colu_arab")
is_accessible.shape
