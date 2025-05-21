#!/bin/bash

# Master output file
output="md5s.txt"
> "$output"  # Clear existing content

# Loop through all .fastq.gz files
for file in *.fastq.gz; do
    [ -e "$file" ] || continue

    # Generate checksum
    checksum=$(md5sum "$file")

    # Save to master file
    echo "$checksum" >> "$output"

    # Save to individual .md5 file
    echo "$checksum" > "${file}.md5"

    echo "Checksums written for: $file"
done

echo "All checksums saved to md5s.txt and individual .md5 files."


### Upload to ENA

lftp -u Webin-65104,R9v@0aLBId -e "mput *.fastq.gz; bye" ftp://webin2.ebi.ac.uk

lftp -u Webin-65104,R9v@0aLBId -e "mput *.bam *.md5; bye" ftp://webin2.ebi.ac.uk

