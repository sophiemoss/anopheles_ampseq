########## MALARIAGEN DATA - phase 3 ##################

## download vcfs from malariagen

# view species calls for malariagen samples

head ~/vo_agam_release/v3/metadata/species_calls_aim_20220528/AG1000G-GM-C/samples.species_aim.csv

awk -F ',' '$7 == "gambiae" && $7 != "intermediate_gambiae_coluzzii" {print $1, $7}' ~/vo_agam_release/v3/metadata/species_calls_aim_20220528/AG1000G-GM-C/samples.species_aim.csv

# view sample ID and snp_genotypes_vcf for Guinea-Bissau dataset

head ~/vo_agam_release/v3/metadata/general/AG1000G-GW/wgs_snp_data.csv | cut -f1,4 -d,

# made file with all sample IDs and now downloading using python script
# python /mnt/storage11/sophie/gitrepos/anophelesmelas_popgen/download_malariagen_samples.py

# merge the malariagen VCF files
bcftools merge --output-type z --output malariagen_gambiae_GB_GM-ABC_merged.vcf.gz *.vcf.gz

# combine with my samples

## bcftools concatenate is used for when you want to merge 2 or more VCFs with the exact same samples, 
## but for example, each VCF is for a single chromosome. 
## bcftools merge is when you want to merge 2 or more VCFs with non-overlapping samples 

bcftools merge /mnt/storage11/sophie/bijagos_mosq_wgs/2022_gambiae_fq2vcf_agamP4/gambiae_nov2022_genomicdb/gambiae_nov2022_genotypedvcf/gambiae_nov2022.2023_07_05.genotyped.vcf.gz malariagen_gambiae_GB_GM-ABC_merged.vcf.gz -O z -o gambiae_malariagen_GB_GM-ABC_Bijagos_merged.vcf.gz

# This vcf is not filtered yet