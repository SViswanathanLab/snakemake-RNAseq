import sys
# log file
sys.stderr = open(snakemake.log[0], "w")
import pandas as pd
import os
import subprocess

# install pyreadr if not previously installed
try:
    import pyreadr
except ImportError:
    subprocess.check_call([sys.executable, "-m", "pip", "install", "pyreadr"])
    import pyreadr


# load the annotation file
annotation = pd.read_table(snakemake.input[-1], header=None)
annotation.columns = ["Gene_ID", "Transcript_ID", "Gene_type", "Gene_symbol", "gene_type"]

# Initialize an empty DataFrame
merged_count = pd.DataFrame(columns=["Gene_ID", "Unstranded_count", "First_strand_count", "Second_strand_count", "Sample", "Strandedness_ratio"])

# load the .ReadsPerGene.out.tab for each sample and merge them together
for count_path in snakemake.input[:-1]:
    count = pd.read_table(count_path, header=None)
    sample_count = count.iloc[4:] # remove the first 4 QC metrics rows

    sample_name = os.path.basename(count_path).replace("_ReadsPerGene.out.tab", "") # Extract sample name from path
    sample_count.columns = ["Gene_ID", "Unstranded_count", "First_strand_count", "Second_strand_count"] # Assign columns
    sample_count["Sample"] = sample_name # Add a column of sample names
    sample_count["Strandedness_ratio"] = sample_count["First_strand_count"].sum() / sample_count["Second_strand_count"].sum() # help decide strandedness

    merged_count = pd.concat([merged_count, sample_count], ignore_index=True) # merge count matrices of all samples

merged_count = merged_count.merge(annotation, left_on="Gene_ID", right_on="Gene_ID", how="left")

# Check the strandedness of the library used
if all(merged_count["Strandedness_ratio"] < 0.8):
    star_Count = merged_count[["Gene_ID", "Gene_symbol", "Gene_type", "Sample", "Second_strand_count"]]
    star_Count = star_Count.rename({"Second_strand_count" : "Count"}, axis=1) # rename the last column as Count
elif all(merged_count["Strandedness_ratio"] > 1.2):
    star_Count = merged_count[["Gene_ID", "Gene_symbol", "Gene_type", "Sample", "First_strand_count"]]
    star_Count = star_Count.rename({"First_strand_count" : "Count"}, axis=1) # rename the last column as Count
elif all(merged_count["Strandedness_ratio"].between(0.8, 1.2)):
    star_Count = merged_count[["Gene_ID", "Gene_symbol", "Gene_type", "Sample", "Unstranded_count"]]
    star_Count = star_Count.rename({"Unstranded_count" : "Count"}, axis=1) # rename the last column as Count
else:
    print("The strandedness cannot be determined.")

# remove duplicates by taking the mean of the count
star_Count_rm_duplicates = star_Count.groupby(["Gene_ID", "Gene_symbol", "Gene_type", "Sample"]).mean().reset_index()
# Pivot the dataframe to create a wide format
Count_wide = star_Count_rm_duplicates.pivot(index=["Gene_ID", "Gene_symbol", "Gene_type"], columns="Sample", values="Count")

# Reset the index to make Gene_ID, gene_symbol, and gene_type regular columns
Count_wide.reset_index(inplace=True)

# save the final output countMatrix
Count_wide.to_csv(snakemake.output[0], sep="\t", index=False, columns=list(Count_wide.columns), header=True, index_label=False)
pyreadr.write_rds(snakemake.output[1], Count_wide)
