## Example MalariaGEN analysis: Investigating population structure with PCA https://malariagen.github.io/vector-data/ag3/examples/pca.html
## activate conda environment in the terminal below with installed python packages
## then use here for import

# %%
# import packages
import os
import bisect
import hashlib
import json
import allel
import numpy as np
import dask
import dask.array as da
from dask.diagnostics import ProgressBar
# quieten dask warnings about large chunks
dask.config.set(**{'array.slicing.split_large_chunks': True})
import pandas as pd
import malariagen_data
import plotly.express as px

# %%
# setup access to malariagen data in google cloud
ag3 = malariagen_data.Ag3()
ag3
# %%

results_dir = "/mnt/storage11/sophie/vo_agam_release/ag3-pca-results"
os.makedirs(results_dir, exist_ok=True)

# %% define some functions for running a PCA
def hash_params(*args, **kwargs):
    """Helper function to hash analysis parameters."""
    o = {
        'args': args,
        'kwargs': kwargs
    }
    s = json.dumps(o, sort_keys=True).encode()
    h = hashlib.md5(s).hexdigest()
    return h


def run_pca(
    region, 
    sample_sets="3.0",
    sample_query=None,
    site_mask="gamb_colu_arab",
    min_minor_ac=3,
    max_an_missing=0,
    n_snps=100_000,
    snp_offset=0,
    n_components=10):
    """Main function to run a PCA.
    
    Parameters
    ----------
    region : str
        Chromosome arm, e.g., '3L', or region, e.g., '3L:1000000-37000000'.
    sample_sets : str or list of str, optional
        Sample sets to analyse.
    sample_query : str, optional
        A pandas query string to select specific samples.
    site_mask : {'gamb_colu_arab', 'gamb_colu', 'arab'}
        Which site mask to apply.
    min_minor_ac : int
        Minimum minor allele count.
    max_an_missing : int
        Maximum number of missing allele calls.
    n_snps : int
        Approximate number of SNPs to use.
    snp_offset : int
        Offset when thinning SNPs.
    n_components : int
        Number of PCA components to retain.

    Returns
    -------
    data : pandas DataFrame
        Data frame with one row per sample, including columns "PC1", "PC2", etc.
    evr : numpy array
        Explained variance ratio per principal component.
    
    """
    
    # construct a key to save the results under
    results_key = hash_params(
        region=region,
        sample_sets=sample_sets,
        sample_query=sample_query,
        site_mask=site_mask,
        min_minor_ac=min_minor_ac,
        max_an_missing=max_an_missing,
        n_snps=n_snps,
        snp_offset=snp_offset,
        n_components=n_components
    )

    # define paths for results files
    data_path = f'{results_dir}/{results_key}-data.csv'
    evr_path = f'{results_dir}/{results_key}-evr.npy'

    try:
        # try to load previously generated results
        data = pd.read_csv(data_path)
        evr = np.load(evr_path)
        return data, evr
    except FileNotFoundError:
        # no previous results available, need to run analysis
        print(f'running analysis: {results_key}')
    
    print('setting up inputs')

    # load sample metadata
    df_samples = ag3.sample_metadata(sample_sets=sample_sets)
    
    # access SNP calls
    ds_snps = ag3.snp_calls(region=region, sample_sets=sample_sets, site_mask=site_mask)

    if sample_query:
        # locate selected samples
        loc_samples = df_samples.eval(sample_query).values
        df_samples = df_samples.loc[loc_samples, :]
        ds_snps = ds_snps.isel(samples=loc_samples)

    # access SNP genotypes
    gt = ds_snps["call_genotype"].data

    print('locating segregating sites within desired frequency range')

    # perform allele count
    ac = allel.GenotypeDaskArray(gt).count_alleles(max_allele=3).compute()
    
    # calculate some convenience variables
    n_chroms = gt.shape[1] * 2
    an_called = ac.sum(axis=1)
    an_missing = n_chroms - an_called
    min_ref_ac = min_minor_ac
    max_ref_ac = n_chroms - min_minor_ac

    # here we choose biallelic sites involving the reference allele
    loc_seg = np.nonzero(ac.is_biallelic() & 
                         (ac[:, 0] >= min_ref_ac) & 
                         (ac[:, 0] <= max_ref_ac) & 
                         (an_missing <= max_an_missing))[0]
    
    print('preparing PCA input data')

    # thin SNPs to approximately the desired number
    snp_step = loc_seg.shape[0] // n_snps
    loc_seg_ds = loc_seg[snp_offset::snp_step]

    # subset genotypes to selected sites
    gt_seg = da.take(gt, loc_seg_ds, axis=0)
    
    # convert to genotype alt counts
    gn_seg = allel.GenotypeDaskArray(gt_seg).to_n_alt().compute()
    
    # remove any edge-case variants where all genotypes are identical
    loc_var = np.any(gn_seg != gn_seg[:, 0, np.newaxis], axis=1)
    gn_var = np.compress(loc_var, gn_seg, axis=0)

    print('running PCA')

    # run the PCA
    coords, model = allel.pca(gn_var, n_components=n_components)
    
    # add PCs to dataframe
    data = df_samples.copy()
    for i in range(n_components):
        data[f'PC{i+1}'] = coords[:, i]
    
    # save results
    evr = model.explained_variance_ratio_
    data.to_csv(data_path, index=False)
    np.save(evr_path, evr)
    print(f'saved results: {results_key}')
    
    return data, evr
    

# %% define plotting functions to help visualise the results

def plot_variance(evr, **kwargs):
    """Plot a bar chart showing variance explained by each principal
    component."""
    
    # prepare variables
    y = evr * 100
    x = [str(i+1) for i in range(len(y))]
    
    # setup plotting options
    plot_kwargs = dict(
        labels={
            'x': 'Principal component',
            'y': 'Explained variance (%)',
        },
        template='simple_white',
        width=600,
        height=400
    )
    # apply any user overrides
    plot_kwargs.update(kwargs)

    # make a bar plot
    fig = px.bar(x=x, y=y, **plot_kwargs)
    fig.show()
    


# %%
def jitter(a, f):
    r = a.max() - a.min()
    return a + f * np.random.uniform(-r, r, a.shape)


def plot_coords(
    data,
    x='PC1',
    y='PC2',
    jitter_frac=0.02,
    random_seed=42,
    **kwargs,
    ):

    # setup data
    data = data.copy()
    
    # apply jitter if desired - helps spread out points when tightly clustered
    if jitter_frac:
        np.random.seed(random_seed)
        data[x] = jitter(data[x], jitter_frac)
        data[y] = jitter(data[y], jitter_frac)
            
    # convenience variables
    data['country_location'] = data['country'] + ' - ' + data['location']
    data['size'] = 1  # hack to allow us to control marker size
    
    # setup plotting options
    plot_kwargs = dict(
        width=700,
        height=500,
        template='simple_white',
        hover_name='sample_id',
        hover_data=[
            'partner_sample_id',
            'sample_set',
            'aim_species', 
            'country', 
            'location', 
            'year', 
        ],
        size='size',
        size_max=8,
        opacity=0.9,
        render_mode='svg',
    )
    # apply any user overrides
    plot_kwargs.update(kwargs)

    # 2D scatter plot
    fig = px.scatter(data, x=x, y=y, **plot_kwargs)
    fig.show()


def plot_coords_3d(
    data,
    x='PC1',
    y='PC2',
    z='PC3',
    jitter_frac=0.02,
    random_seed=42,
    **kwargs,
    ):

    # setup data
    data = data.copy()
    
    # apply jitter if desired - helps spread out points when tightly clustered
    if jitter_frac:
        np.random.seed(random_seed)
        data[x] = jitter(data[x], jitter_frac)
        data[y] = jitter(data[y], jitter_frac)
        data[z] = jitter(data[z], jitter_frac)
            
    # convenience variables
    data['country_location'] = data['country'] + ' - ' + data['location']
    
    # setup plotting options
    plot_kwargs = dict(
        width=700,
        height=500,
        hover_name='sample_id',
        hover_data=[
            'partner_sample_id',
            'sample_set',
            'aim_species', 
            'country', 
            'location', 
            'year', 
        ],
    )
    # apply any user overrides
    plot_kwargs.update(kwargs)

    # 3D scatter plot
    fig = px.scatter_3d(data, x=x, y=y, z=z, **plot_kwargs)
    fig.show()
    

# %%
# choose colours for species
species_palette = px.colors.qualitative.Plotly
species_color_map = {
    'gambiae': species_palette[0],
    'coluzzii': species_palette[1],
    'arabiensis': species_palette[2],
    'intermediate_gambiae_coluzzii': species_palette[3],
    'intermediate_arabiensis_gambiae': species_palette[4],
}
# %% EXAMPLE. Run PCA using samples from Central African Republic and SNPs from chromosome arm 3L.

data, evr = run_pca(
    region='3L', 
    sample_sets="AG1000G-CF",
)

# %%
data.head()
# %%
evr

# %%
# make a scatterplot for the PCA results for the first two components:
title = 'Central African Republic (3L)'
plot_coords(data, x='PC1', y='PC2',
            color='aim_species', 
            color_discrete_map=species_color_map, 
            title=title)
# %% how many principal components should I look at? Gives boxplot of explained variance for each principal component
plot_variance(evr, title=title)

#%% Analysis: Uganda

title = 'Uganda (3L)'
data, evr = run_pca(
    region='3L', 
    sample_sets="AG1000G-UG",
)
plot_variance(evr, title=title)
plot_coords(data, x='PC1', y='PC2',
            color='aim_species', 
            color_discrete_map=species_color_map, 
            title=title)


# %% re-run with just An. gambiae and exclude two outliers that were on the previous PCA
title = 'Uganda <i>An. gambiae</i> (3L)'
data, evr = run_pca(
    region='3L', 
    sample_sets="AG1000G-UG",
    sample_query=(
        "aim_species == 'gambiae' and "
        "sample_id not in ['AC0223-C', 'AC0240-C']"
    ),
)
plot_variance(evr, title=title)
plot_coords(data, x='PC1', y='PC2',
            color='location', 
            title=title)

# %% Analysis: Tanzania
title = 'Tanzania (3L)'
data, evr = run_pca(
    region='3L', 
    sample_sets="AG1000G-TZ",
)
plot_variance(evr, title=title)

# %%
plot_coords_3d(data, x='PC1', y='PC2', z='PC3', 
               color='aim_species', 
               color_discrete_map=species_color_map, 
               title=title,
               jitter_frac=0.05)

# %% Analysis: East African An. arabiensis

title = 'East African <i>An. arabiensis</i> (3L)'
data, evr = run_pca(
    region='3L', 
    sample_sets="3.0",
    sample_query=(
        "aim_species == 'arabiensis' and "
        "country in ['Uganda', 'Tanzania', 'Kenya', 'Malawi']"
    )
)
plot_variance(evr, title=title)
plot_coords(data, x='PC1', y='PC2',
            color='country', 
            title=title)

# %%
