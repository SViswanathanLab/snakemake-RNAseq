import sys
# log file
sys.stderr = open(snakemake.log[0], "w")
import pandas as pd
import os
import subprocess

try:
    import pyreadr
except ImportError:
    subprocess.check_call([sys.executable, "-m", "pip", "install", "pyreadr"])
    import pyreadr

# load the annotation file
annotation = pd.read_table(snakemake.input[-1], header=None)
annotation.columns = ["Gene_ID", "Transcript_ID", "Gene_type", "Gene_symbol", "gene_type"]

# Initialize an empty DataFrame
merged_count = pd.DataFrame(columns=["Name", "Length", "EffectiveLength", "TPM", "NumReads", "Sample"])

# load the quant.sf for each sample and merge them together
for quant_path in snakemake.input[:-1]:
    quant = pd.read_table(quant_path)

    sample_name = quant_path.split('/')[2] # Extract sample name from path
    # sample_count.columns = ["Gene_ID", "Unstranded_count", "First_strand_count", "Second_strand_count"] # Assign columns
    quant["Sample"] = sample_name # Add a column of sample names
    # sample_count["Strandedness_ratio"] = sample_count["First_strand_count"].sum() / sample_count["Second_strand_count"].sum() # help decide strandedness

    merged_count = pd.concat([merged_count, quant], ignore_index=True) # merge quant matrices of all samples

merged_count = merged_count.merge(annotation, left_on="Name", right_on="Transcript_ID", how="left")

# Pivot the dataframe to create a wide format
TPM_wide = merged_count.pivot(index=["Transcript_ID", "Gene_symbol", "Gene_type"], columns="Sample", values="TPM")
# Reset the index to make transcript_ID, gene_symbol, and gene_type regular columns
TPM_wide.reset_index(inplace=True)


# save the final output countMatrix
TPM_wide.to_csv(snakemake.output[0], sep="\t", index=False, header=True, index_label=False)
pyreadr.write_rds(snakemake.output[1], TPM_wide)
