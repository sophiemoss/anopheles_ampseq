## this python script is designed to download a set of vcf and vcf.tbi files
## where each srr number is on a different line
## to run the script, go to the dir where you want to download the files and use the command
## python /mnt/storage11/sophie/gitrepos/anophelesmelas_popgen/download_malariagen_samples.py

import subprocess

# Initialize an empty list to store the sample IDs

sampleIDs = []

# Read the sample ID numbers from the file 'malariagen_samples_for_comparison.csv'
with open('malariagen_samples_for_comparison.csv', 'r') as f:
    lines = f.readlines()  # Read all lines from the file
    
    # Iterate through the lines, starting from the second line (index 1) to exclude the header
    for line in lines[1:]:
        # Split each line by comma and extract the value in the first column (index 0)
        values = line.strip().split(',')
        sampleID = values[0]
        sampleIDs.append(sampleID)

# Loop through each SRR number and run the fasterq-dump command
for sampleID in sampleIDs:
    vcf_command = f'wget --no-clobber https://vo_agam_output.cog.sanger.ac.uk/{sampleID}.vcf.gz'
    tbi_command = f'wget --no-clobber https://vo_agam_output.cog.sanger.ac.uk/{sampleID}.vcf.gz.tbi'

    try:
        # Run the vcf command using subprocess.run()
        subprocess.run(vcf_command, shell=True, check=True)
        print(f"Successfully downloaded vcf file: {sampleID}")
        # Run the tbi command using subprocess.run()
        subprocess.run(tbi_command, shell=True, check=True)
    except subprocess.CalledProcessError:
        print(f"Error downloading vcf files: {sampleID}")


#
