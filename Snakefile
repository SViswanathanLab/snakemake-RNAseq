configfile: "./config/config.yaml"

# reference files paths
annotation_path = config["annotation_path"]
STAR_index_path = config["STAR_index_path"]
gtf_path = config["gtf_path"]
salmon_index_path = config["salmon_index_path"]
fa_path = config["fa_path"]

# Read samples.tsv to extract Samples, fq1, fq2
Samples = []
fq1 = []
fq2 = []
with open(config["samples"]) as file:
    for line in file:
        l = line.strip().split(',')
        if len(l) == 3:
            Samples.append(l[0])
            fq1.append(l[1])
            fq2.append(l[2])

# Assume the full sample name is within the fq1 and fq2 file name
# Get the suffix of the seq file name excluding the sample name, used in starAlign.smk and salmonQuant.smk
fq1_suffix = fq1[0].replace(Samples[0], "")
fq2_suffix = fq2[0].replace(Samples[0], "")

###### Get wildcards for Reads (accomodate 3 patterns: .fastq.gz, .fq.gz, .fastq)
# Define the patterns with possible suffixes
patterns = ["data/{Read}.fastq.gz", "data/{Read}.fastq", "data/{Read}.fq.gz", "data/{Read}.fq"]
# Use glob_wildcards to find matching files
Reads = set()
for pattern in patterns:
    Reads.update(glob_wildcards(pattern).Read)
# Convert set to list
Reads = list(Reads)

###### Perform fastQC based on the file type
if ".fastq.gz" in fq1_suffix:
    include: "rules/preAlignQC_fastq_gz.smk",
elif ".fq.gz" in fq1_suffix:
    include: "rules/preAlignQC_fq_gz.smk",
elif ".fastq" in fq1_suffix:
    include: "rules/preAlignQC_fastq.smk",
else:
    include: "rules/preAlignQC_fq.smk",

include: "rules/starAlign.smk",
include: "rules/salmonQuant.smk",
include: "rules/countMatrix.smk",
include: "rules/rsem.smk",

rule all:
    input:
        expand("results/fastqc_results/{Read}_fastqc.zip", Read=Reads),
        expand("results/fastqc_results/{Read}_fastqc.html", Read=Reads),
        expand("results/STAR_results/{Sample}/{Sample}_Log.final.out", Sample=Samples),
        "results/star_wide_countMatrix.csv",
        "results/star_wide_countMatrix.Rds",
        expand("results/salmon_results/{Sample}/quant.sf", Sample=Samples),
        "results/salmon_wide_TPM_Matrix.csv",
        "results/salmon_wide_TPM_Matrix.Rds",
        multiext("rsem_ref/human_gencode", ".chrlist", ".n2g.idx.fa", ".transcripts.fa", ".grp", ".seq", ".idx.fa", ".ti"),
        expand("results/rsem_results/{Sample}/{Sample}.genes.results", Sample=Samples),
        expand("results/rsem_results/{Sample}/{Sample}.isoforms.results", Sample=Samples),
        "results/rsem_geneLevel_wide_TPM_Matrix.csv",
        "results/rsem_geneLevel_wide_TPM_Matrix.Rds",
        "results/rsem_isoformLevel_wide_TPM_Matrix.csv",
        "results/rsem_isoformLevel_wide_TPM_Matrix.Rds",
