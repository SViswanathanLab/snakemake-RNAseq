import sys
# log file
sys.stderr = open(snakemake.log[0], "w")
import pandas as pd
import os
import subprocess

try:
    import pyreadr
except ImportError:
    subprocess.check_all([sys.executable, "-m", "pip", "install", "pyreadr"])
    import pyreadr

# load the annotation file
annotation = pd.read_table(snakemake.input[-1], header=None)
annotation.columns = ["Gene_ID", "Transcript_ID", "Gene_type", "Gene_symbol", "gene_type"]

# Initialize an empty DataFrame
merged_count = pd.DataFrame(columns=["gene_id", "transcript_id(s)", "length", "effective_length", "expected_count", "TPM", "FPKM", "Sample"])

# load the quant.sf for each sample and merge them together
for quant_path in snakemake.input[:-1]:
    quant = pd.read_table(quant_path)
    sample_name = os.path.basename(quant_path).replace(".genes.results", "") # extract sample name from path
    quant["Sample"] = sample_name # Add a column of sample names

    merged_count = pd.concat([merged_count, quant], ignore_index=True) # merge quant matrices of all samples

merged_count = merged_count.merge(annotation, left_on="gene_id", right_on="Gene_ID", how="left")


# remove duplicates by taking the mean of the count
merged_count = merged_count[["Gene_ID", "Gene_symbol", "Gene_type", "Sample", "TPM"]]
TPM_rm_duplicates = merged_count.groupby(["Gene_ID", "Gene_symbol", "Gene_type", "Sample"]).mean().reset_index()
# Pivot the dataframe to create a wide format
TPM_wide = TPM_rm_duplicates.pivot(index=["Gene_ID", "Gene_symbol", "Gene_type"], columns="Sample", values="TPM")
# Reset the index to make transcript_ID, gene_symbol, and gene_type regular columns
TPM_wide.reset_index(inplace=True)

# save the final output countMatrix
TPM_wide.to_csv(snakemake.output[0], sep="\t", index=False, header=True, index_label=False)
pyreadr.write_rds(snakemake.output[1], TPM_wide)
