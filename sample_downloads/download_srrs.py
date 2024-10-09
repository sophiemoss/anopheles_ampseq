## this python script is designed to download a set of SRR numbers in a file called srr.txt
## where each srr number is on a different line
## to run the script, go to the dir where you want to download the files and use the command
## python /mnt/storage11/sophie/gitrepos/anophelesmelas_popgen/download_srrs.py

import subprocess

# Read the SRR numbers from the file 'srrs.txt' and store them in a list
with open('srrs.txt', 'r') as f:
    srr_numbers = f.read().splitlines()

# Loop through each SRR number and run the fasterq-dump command
for srr_number in srr_numbers:
    command = f'fasterq-dump {srr_number}'
    try:
        # Run the command using subprocess.run()
        subprocess.run(command, shell=True, check=True)
        print(f"Successfully processed SRR number: {srr_number}")
    except subprocess.CalledProcessError:
        print(f"Error processing SRR number: {srr_number}")
