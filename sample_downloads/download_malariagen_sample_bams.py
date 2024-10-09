## this python script is designed to download a set of vcf and vcf.tbi files
## where each srr number is on a different line
## to run the script, go to the dir where you want to download the files and use the command
## python /mnt/storage11/sophie/gitrepos/anophelesmelas_popgen/download_malariagen_sample_bams.py

import subprocess

# Initialize an empty list to store the sample IDs

sampleIDs = []

# Read the sample ID numbers from the file 'malariagen_samples_for_comparison.csv'
with open('malariagen_sampleIDs.csv', 'r') as f:
    lines = f.readlines()  # Read all lines from the file
    
    # Iterate through the lines, starting from the second line (index 1) to exclude the header
    for line in lines[1:]:
        # Split each line by comma and extract the value in the first column (index 0)
        values = line.strip().split(',')
        sampleID = values[0]
        sampleIDs.append(sampleID)

# Loop through each SRR number and run the fasterq-dump command
for sampleID in sampleIDs:
    bam_command = f'wget --no-clobber https://vo_agam_output.cog.sanger.ac.uk/{sampleID}.bam'
    bai_command = f'wget --no-clobber https://vo_agam_output.cog.sanger.ac.uk/{sampleID}.bam.bai'

    try:
        # Run the vcf command using subprocess.run()
        subprocess.run(bam_command, shell=True, check=True)
        print(f"Successfully downloaded bam file: {sampleID}")
        # Run the tbi command using subprocess.run()
        subprocess.run(bai_command, shell=True, check=True)
    except subprocess.CalledProcessError:
        print(f"Error downloading bam files: {sampleID}")


#
