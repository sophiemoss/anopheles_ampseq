### ______ Run amplicon sequencing pipeline on David samples ______ ###

## You need to make an index file with the samples to process and no others, and use random barcode eg CCCCC because you have already demultiplexed these files

ls *R1_001.fastq.gz | sed 's/_R1_001.fastq.gz//' > amp_index.txt

####_____ make samples.txt file ––––####

ls *R1_001.fastq.gz | sed 's/_R1_001.fastq.gz//' > samples.txt

####_____ Rename the fastq files to be the correct format for the python script ––––####

for f in *_R1_001.fastq.gz; do
    base=${f%%_R1_001.fastq.gz}
    mv "$f" "${base}_1.fastq.gz"
done

for f in *_R2_001.fastq.gz; do
    base=${f%%_R2_001.fastq.gz}
    mv "$f" "${base}_2.fastq.gz"
done

### Edit python script to remove demux step ###

###################

echo "🔬 Starting amplicon analysis for all samples in amp_index.csv"

# Using dummy values for --read1 and --read2 because they are required by argparse but not actually used in the script
python amplicon_script_nodemux.py \
--read1 dummy_1.fastq.gz \
--read2 dummy_2.fastq.gz \
--ref Pfalciparum.genome.fasta \
--gff Pfalciparum.genome.modified.new.gff3 \
--bed desired_regions_fixed.bed \
--position-info pf3d7.positions.txt \
--index-file amp_index.csv \
> ampseq_log.txt 2>&1

echo "✅ Analysis complete. See ampseq_log.txt for details."

####____WHAT IS THE AMPLICON PIPELINE PRODUCING?____####

## File | Description
## *coverage.txt | Raw per-base depth (from sambamba depth base)
## *mosdepth.summary.txt | Summary stats per region or whole genome
## *region_coverage.txt | Mean coverage per region (via bedtools coverage -mean)
## *per-base.bed.gz | Depth at every base (from mosdepth)
## *thresholds.bed.gz | BED with threshold coverage (e.g., how many bases >1x, >10x, etc.)
## *regions.bed.gz | Per-region coverage (from mosdepth)

### Instead of using naive variant caller, I used just the gatk variant caller and did the merging of the gatk vcfs instead, followed by the rest of the script.

####____COMBINE THE AMPLICON SNPS WITH THE WGS SNPS AND FILTER TOGETHER____####

# output from WGS
pfalciparum.2025_05_07.genotyped.vcf.gz 

# keep only SNPs
snps_pfalciparum.2025_05_07.genotyped.vcf.gz

# slide to be just the bedfile used for amplicons
bcftools view -R desired_regions_fixed.bed -Oz -o sliced_snps_pfalciparum.2025_05_07.genotyped.vcf.gz snps_pfalciparum.2025_05_07.genotyped.vcf.gz


### __________________________________________________________________ ###
### _______MERGE THE AMPLICON AND THE WGS DATA____ ###
### __________________________________________________________________ ###

bcftools merge sliced_snps_pfalciparum.2025_05_07.genotyped.vcf.gz amplicon_snps.vcf.gz -Oz -o wgs_and_amp_merged.vcf.gz


### __________________________________________________________________ ###
### _______FILTER THE AMPLICON AND WGS SNP DATA ____ ###
### __________________________________________________________________ ###


# STEP ONE: Do GATK variant filtration

gatk VariantFiltration \
-R Pfalciparum.genome.fasta \
-V wgs_and_amp_merged.vcf.gz \
-filter "QD < 5.0" --filter-name "QD5" \
-filter "QUAL < 20.0" --filter-name "QUAL30" \
-filter "SOR > 3.0" --filter-name "SOR3" \
-filter "FS > 60.0" --filter-name "FS60" \
-filter "MQ < 40.0" --filter-name "MQ40" \
-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
-O gatk_tagged_wgs_and_amp_merged.vcf.gz

bcftools view -f 'PASS' gatk_tagged_wgs_and_amp_merged.vcf.gz -Oz -o gatk_filtered_wgs_and_amp_merged.vcf.gz

# STEP TWO: Filter out variants with read depth below 5

bcftools filter -i 'FMT/DP>5' -S . gatk_filtered_wgs_and_amp_merged.vcf.gz | bcftools sort -Oz -o DP_gatk_filtered_wgs_and_amp_merged.vcf.gz


# STEP THREE: Filter for Minor Allele Frequency or Allele Count to make sure all variants remaining are true variants.

bcftools filter -i 'MAF>0.01' DP_gatk_filtered_wgs_and_amp_merged.vcf.gz -Oz -o MAF_DP_gatk_filtered_wgs_and_amp_merged.vcf.gz

### Double check this has defo got rid of non variants by doing an additional filter for AC > 0

bcftools filter -i 'AC>0' MAF_DP_gatk_filtered_wgs_and_amp_merged.vcf.gz -Oz -o AC_MAF_DP_gatk_filtered_wgs_and_amp_merged.vcf.gz


### ___________________________________________ ###
### _________________ Analysis –––––––––––––––– ###
### ___________________________________________ ###


### Make a PCA of the merged vcf file to see if WGS and amplicon are sitting together ###

# Create distance matrix

plink --vcf AC_MAF_DP_gatk_filtered_wgs_and_amp_merged.vcf.gz --distance square --double-id --out pfal_wgs_and_amp --threads 6 --allow-extra-chr

# Make PCA

Rscript /mnt/storage11/sophie/gitrepos/david_analysis/PCA_wgs_vs_amplicon.R


### TESTED TO SEE WHETHER SPLIT IS DUE TO MISSINGNESS, IT IS NOT ###

# It looks like the samples do separate based on amplicon vs WGS data. 
# Checking if this is due to coverage (missingness), by changing all missing values to 0 and seeing if the groups still separate

bcftools +fixploidy DP_MAF_AC_wgs_and_amp_merged.vcf.gz -- -f 2 | \
bcftools +setGT -- -t . -n 0/0 | \
bcftools view -Oz -o test_set_missing_to_ref.vcf.gz





### ______ PREPROCESSING OF WGS DATA TO OBTAIN VCF FILE USED ____ ###
### ______ Create Genomics DB for WGS Data ____ ###

### To make a genomics database of sample VCFs, use the following:

ls *.g.vcf.gz | sed 's/.g.vcf.gz//' > database_samples.txt #removed samples that do not pass the cut off

## make a list of the sample names that pass the cutoff to be included in the database
## Use merge_vcfs.py import

/mnt/storage11/sophie/fastq2matrix/scripts/merge_vcfs.py import --sample-file database_samples.txt --ref Pfalciparum.genome.fasta --prefix pfalciparum --vcf-dir . --threads 25

## now merge VCF files
## use merge_vcfs.py genotype

/mnt/storage11/sophie/fastq2matrix/scripts/merge_vcfs.py genotype --ref Pfalciparum.genome.fasta --prefix pfalciparum --threads 25 > david_plasmodium_mergevcf_log.txt 2>&1

# resulting vcf should be called yourprefix.genotyped.vcf.gz


