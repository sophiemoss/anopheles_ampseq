## Example MalariaGEN analysis: SNP and Gene Effects https://malariagen.github.io/vector-data/ag3/examples/snp-effects.html
## activate conda environment in the terminal below with installed python packages
## then use here for import

# %%
import numpy as np
import bokeh.plotting as bkplt
import bokeh.models as bkmod
import bokeh.layouts as bklay
import bokeh.io as bkio
import bokeh.palettes as bkpal
import pandas as pd
import plotly.express as px
import malariagen_data

# %%
bkio.output_notebook()
# %% set up access to malariagen in google cloud
ag3 = malariagen_data.Ag3()
ag3

######### GENE ANNOTATIONS ##########

# %%
df_geneset = ag3.geneset().set_index("ID")
df_geneset
# %% nice interactive plot
ag3.plot_genes("2L");
# %% looking up the transcripts for the doublesex gene, AGAP004050.
df_geneset.query("Parent == 'AGAP004050'")

## visualise transcripts for the doublesex gene
## the boxes are exons and the connecting lines are introns.
## CDS (coding sequences) are blue, five prime UTRs are green, 3 prime UTRs are red.

# %%
ag3.plot_transcript("AGAP004050-RA");

# %%
ag3.plot_transcript("AGAP004050-RB");

# %% Different example, RDL gene

df_geneset.query("Parent == 'AGAP006028'")

# %%
ag3.plot_transcript("AGAP006028-RA");

# %%
ag3.plot_transcript("AGAP006028-RB");

# %%
ag3.plot_transcript("AGAP006028-RC");


######### SNP EFFECTS ##########
# %%
# Rdl gene, first transcript - change this to investigate a different gene or transcript
transcript = "AGAP006028-RA"

# compute effects for all SNPs in chosen transcript
df_effects = ag3.snp_effects(
    transcript=transcript, 
)
df_effects

# %% 
df_effects.groupby(['impact', 'effect']).size()

## note that not all of these SNPs will have been observed in the Ag3.0 datset, so we also need information on allele frequencies

### SNP ALLELE FREQUENCIES ###

# %% compute allele frequencies for our gene of interest

df_af = ag3.snp_allele_frequencies(
    transcript=transcript, 
    cohorts="admin1_year", 
    sample_sets="3.0",
    sample_query="country in ['Burkina Faso', 'Uganda']"
)
df_af

## The allele frequencies are given in columns named after the cohorts. The allele frequency is given as a fraction, 
# e.g., 1.0 means all samples carry the alternate allele, 
#  0.5 means half the samples carry the alternate allele, and 0 means all samples carry the reference allele.

###### Anlysis of genetic variation present in our gene and populations of interest ######
# %% filter dataframe of SNPs, keep only non-synonymous snps and remove low fq SNPs where AF is less than 2% across all populations of interest
df_snps_filtered = df_af.query("effect == 'NON_SYNONYMOUS_CODING' and max_af > 0.02")
df_snps_filtered

# %% visualise these frequencies as a heatmap
import nbformat

#%%
ag3.plot_frequencies_heatmap(df_snps_filtered)

# %% visualise data in genomic context

def plot_snps_transcript(transcript, data, width=750, height=300, palette='Category10'):
    data = data.reset_index()

    # hover tooltips
    tooltips = [
        ("position", '@contig:@position{,}'),
        ("alleles", '@ref_allele>@alt_allele'),
        ("pass", "@pass_gamb_colu_arab, @pass_gamb_colu, @pass_arab"),
        ("impact", '@impact'),
        ("effect", '@effect'),
        ("aa_change", '@aa_change'),
        ("frequency", '@frequency{%f} (@cohort)'),
    ]

    fig1 = bkplt.figure(
        title=f'Transcript - {transcript}',
        tools='xpan,xzoom_in,xzoom_out,xwheel_zoom,reset,hover',
        active_scroll='xwheel_zoom',
        active_drag='xpan',
        width=width, 
        height=height, 
        tooltips=tooltips,
        toolbar_location="above")

    # set up colors
    frq_cols = [c for c in data.columns if c.startswith("frq_")]
    cohorts = [c.split("frq_")[1] for c in frq_cols]
    palette = bkpal.all_palettes[palette]
    colors = palette[len(cohorts)]

    # plot allele frequencies
    for coh, color in zip(cohorts, colors):
        df = data.copy()
        df['frequency'] = df[f"frq_{coh}"]
        df['cohort'] = coh
        fig1.triangle("position", "frequency", 
                      size=8, 
                      color=color,
                      source=df,
                      legend_label=coh)

    # tidy up the plot
    fig1.y_range = bkmod.Range1d(0, 1)
    fig1.yaxis.axis_label = f'Alt allele frequency'
    fig1.xaxis.visible = False
    fig1.add_layout(fig1.legend[0], 'right')
    fig1.legend.click_policy="hide"

    # plot transcript
    fig2 = ag3.plot_transcript(
        transcript, 
        width=width, 
        height=80, 
        show=False, 
        x_range=fig1.x_range
    )
    fig2.toolbar.logo = None 
    fig2.toolbar_location = None
    fig2.title = None

    bkplt.show(bklay.column(fig1, fig2))


# %%
plot_snps_transcript(transcript, data=df_snps_filtered)
# %%
