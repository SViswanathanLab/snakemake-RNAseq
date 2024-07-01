# snakemake-RNAseq
This pipeline performs a standard RNAseq analysis, including fastQC, STAR alignment, RSEM & salmon quantification.

## Directory structure
```
.
├── config            # Contains sample sheet (samples.tsv) and config file (config.yaml)
├── rules             # Snakemake rules
├── scripts           # Scripts to run each step of RNAseq
├── README.md
├── Snakefile         # Snakemake workflow
└── job_template.sh   # A template for argos job submission, modified by users as needed

```
## Installation
Clone the pipeline using the following command
```
git clone https://github.com/SViswanathanLab/snakemake-RNAseq.git
```
There should be a folder named ```snakemake-RNAseq``` in users' working directory. If users modify the name of this folder, they should also modify the contents of ```job_template.sh``` to make the name of the directory containing ```Snakefile``` consistent with the modified folder name. 

## Usage 
### Instructions for preparing sample sheet
* Paired-end data is assumed.
* 3 types of RNAseq data formats are accommodated: **.fastq.gz, .fq.gz, .fastq**
* The ```samples.tsv``` file is an example sample sheet.
* Users should modify ```samples.tsv``` to have the first column consisting of sample names, the second column consisting of fq1 file names, and the third column consisting of fq2 file names. Each column is separated by **one space**. 
* The fq1 & fq2 file names must contain the full sample names.
  
  For example: 
```
293T-TFE3-1 293T-TFE3-1_R1_001.fastq.gz 293T-TFE3-1_R2_001.fastq.gz
293T-TFE3-2 293T-TFE3-2_R1_001.fastq.gz 293T-TFE3-2_R2_001.fastq.gz
```
### Input files
* Users should **create a folder named ```data``` in the directory of ```snakemake-RNAseq```**.
* The fq1 & fq2 files for analysis should be copied to ```data```.

### Run snakemake
* Since the analysis takes time, we would recommend submitting it as a job to avoid being interrupted. The template of the job submission script is provided as ```job_template.sh```.
* Modify the name of the directory containing the Snakefile as needed.

### Output
* The results are saved in the folder ```snakemake-RNAseq/results```.
    * ```snakemake-RNAseq/results/fastqc_results``` contains the fastqc results.
    * ```snakemake-RNAseq/results/STAR_results``` contains the STAR results, and each subfolder is named by the sample name.
    * ```snakemake-RNAseq/results/salmon_results``` contains the salmon results, and each subfolder is named by the sample name.
    * ```snakemake-RNAseq/results/star_wide_countMatrix.csv``` and ```snakemake-RNAseq/results/star_wide_countMatrix.Rds``` contain the star count matrix, with genes as rows and samples as columns, in both .csv and .Rds format.
    * ```snakemake-RNAseq/results/salmon_wide_TPM_Matrix.csv``` and ```snakemake-RNAseq/results/salmon_wide_TPM_Matrix.Rds``` contain the salmon TPM matrix, with genes as rows and samples as columns, in both .csv and .Rds format.
    * ```snakemake-RNAseq/results/rsem_geneLevel_wide_TPM_Matrix.csv``` and ```snakemake-RNAseq/results/rsem_geneLevel_wide_TPM_Matrix.Rds``` contain the rsem gene-level TPM matrix, with genes as rows and samples as columns, in both .csv and .Rds format.
    * ```snakemake-RNAseq/results/rsem_isoformLevel_wide_TPM_Matrix.csv``` and ```snakemake-RNAseq/results/rsem_isoformLevel_wide_TPM_Matrix.Rds``` contain the rsem transcript-level TPM matrix, with genes as rows and samples as columns, in both .csv and .Rds format.
* ```snakemake-RNAseq/rsem_ref``` contains the reference files generated for rsem quantification.
* ```snakemake-RNAseq/logs``` contains the log files for running each step of this analysis, for debugging.
* ```snakemake-RNAseq/joblogs``` contains the log files for job submission, for debugging.

### config
```config.yaml``` contains the information about versions of each tool used, reference file paths
* Module versions (latest ones globally installed on argos):
    * fastqc: 0.11.7
    * star: 2.7.10a
    * rsem: 1.3.1
    * salmon: 1.10.1
    * snakemake: 7.25.0
* Reference file directory: ```/mnt/storage/labs/sviswanathan/snakemake_RNAseq_2024/Human_genome_2024/```
    * Can be modified as needed
