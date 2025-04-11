import os
import pandas as pd
from glob import glob

# Path to your folder of coverage files
coverage_folder = "/mnt/storage11/sophie/env_dna/double_demux_fastqs/"
coverage_files = glob(os.path.join(coverage_folder, "*_coverage_mean.txt"))
print("Found coverage files:", len(coverage_files))
print(coverage_files)

all_data = []

for filepath in coverage_files:
    sample_name = os.path.basename(filepath).replace("_coverage_mean.txt", "")
    df = pd.read_csv(filepath, sep="\t", header=None, names=["CHROM", "START", "END", "AMPLICON", "DEPTH"])
    df["SAMPLE"] = sample_name
    all_data.append(df)

# Combine all rows
combined_df = pd.concat(all_data)

# Pivot into a matrix: rows = amplicons, columns = samples
matrix_df = combined_df.pivot(index="AMPLICON", columns="SAMPLE", values="DEPTH")

# Save matrix
matrix_df.to_csv("amplicon_coverage_matrix.tsv", sep="\t")

# Binary success matrix: coverage > threshold (e.g., 10x)
success_matrix = (matrix_df > 10).astype(int)

# % success per amplicon and per sample
amplicon_success_rate = success_matrix.mean(axis=1) * 100
sample_success_rate = success_matrix.mean(axis=0) * 100

# Print to terminal
print("Amplicon success %: % of samples with coverage ≥ 10x for each amplicon")
print(amplicon_success_rate)

print("\nSample success %: % of amplicons within each sample which have coverage ≥ 10x")
print(sample_success_rate)

# Save summary to file
with open("amplicon_coverage_matrix_summary.txt", "w") as f:
    f.write("Amplicon success (% of samples with coverage ≥ 10x for each amplicon):\n")
    f.write(amplicon_success_rate.sort_index().to_string())
    f.write("\n\nSample success (% of amplicons within each sample which have coverage ≥ 10x):\n")
    f.write(sample_success_rate.sort_index().to_string())
