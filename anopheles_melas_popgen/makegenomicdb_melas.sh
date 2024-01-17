## download additional samples from ncbi using fasterq-dump
#conda activate fasterq-dump
#fasterq-dump SRR000001

### Put all of the samples in one directory for making database
### Put reference genome in the same directory

## to make a genomics database of sample VCFs, use the following
ls *.g.vcf.gz | sed 's/.g.vcf.gz//' > fastq2vcfsamples.txt

/mnt/storage11/sophie/fastq2matrix/scripts/merge_vcfs.py import --sample-file fastq2vcfsamples.txt --ref genomic.fasta --prefix melas_global_melas_bijagos_2019 --vcf-dir .

## now merge VCF files

merge_vcfs.py genotype --ref genomic.fasta --prefix gambiae_bijagos_2022 > mergevcf_log.txt 2>&1

# resulting vcf is called gambiae_bijagos_2022.2023_07_25.genotyped.vcf.gz

# concerned that this genotyping does not use a set of validated variants like I had for plasmodium
# this has a genotyping pipeline:
# https://github.com/malariagen/pipelines/blob/v0.0.4/docs/specs/snp-genotyping-vector.md
# pipeline being sent from MalariaGEN.