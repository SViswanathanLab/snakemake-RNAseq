# snakemake-RNAseq
This pipeline performs a standard RNAseq analysis, including fastQC, STAR alignment, RSEM & salmon quantification.

## Directory structure
```
.
├── config            # Contains sample sheet (samples.tsv) and config file (config.yaml)
├── rules             # Snakemake rules
├── scripts           # Scripts to run each step of RNAseq
├── README.md
└── Snakefile         # Snakemake workflow

```
## Installation
### Option 1: Download the package

* Choose "Download ZIP"
<img width="935" alt="Screenshot 2024-07-09 at 3 37 12 PM" src="https://github.com/SViswanathanLab/snakemake-RNAseq/assets/143852554/394a5529-ffca-4222-acf1-6936989d65a8">

* The folder named ```snakemake-RNAseq-main``` is downloaded.
* Transfer the folder to users' working directory on argos.
  ```
  scp -r path/to/snakemake-RNAseq-main <USER_ID>@argos-stgw2.dfci.harvard.edu:/mnt/storage/home/<USER_ID>/
  ```
* Log onto argos:
  ```
  ssh USER_ID@argos.dfci.harvard.edu
  ```
* Change the name of the folder to ```snakemake-RNAseq```
  ```
  mv $HOME/snakemake-RNAseq-main $HOME/snakemake-RNAseq
  ```

### Option 2: git clone
Clone the pipeline using the following command
```
git clone https://github.com/SViswanathanLab/snakemake-RNAseq.git
```


**Make sure there is a folder named ```snakemake-RNAseq``` in users' working directory.**


## Usage 
### Instructions for preparing sample sheet
* Paired-end data is assumed.
* 4 types of RNAseq data formats are accommodated: **.fastq.gz, .fq.gz, .fastq, .fq**
* The ```config/samples.tsv``` file is an example sample sheet.
* Users should modify ```config/samples.tsv``` to have the first column consisting of sample names, the second column consisting of fq1 file names, and the third column consisting of fq2 file names. Each column is separated by **one space**. 
* The fq1 & fq2 file names must contain the full sample names.
  
  For example: 
```
293T-TFE3-1 293T-TFE3-1_R1_001.fastq.gz 293T-TFE3-1_R2_001.fastq.gz
293T-TFE3-2 293T-TFE3-2_R1_001.fastq.gz 293T-TFE3-2_R2_001.fastq.gz
```
### Input files
* The fq1 & fq2 files for analysis should be copied to ```data```.
  ```
  cp -r path/to/<fq_files_folder> $HOME/snakemake-RNAseq/
  ```
* Users should **change the name of the folder containing fq files into ```data```**.
  ```
  mv $HOME/snakemake-RNAseq/<fq_files_folder> $HOME/snakemake-RNAseq/data
  ```
  
### Run snakemake
* **Step 1: Change into the directory ```snakemake-RNAseq```**
  
  ```
  cd $HOME/snakemake-RNAseq
  ```
  
* **Step 2: Activate the environment with snakemake installed & install plugin for cluster submission**
  
  ```
  source /mnt/storage/apps/Mambaforge-23.1.0-1/etc/profile.d/conda.sh
  conda activate snakemake

  pip install snakemake-executor-plugin-cluster-generic
  ```
  
* **Step 3: Run snakemake pipeline**
  
  ```
  snakemake --unlock
  snakemake --executor cluster-generic --jobs 50 --latency-wait 60 --cluster-generic-submit-cmd "qsub -l h_vmem=256G, -pe pvm 32 -o $HOME/snakemake-RNAseq/joblogs/ -e $HOME/snakemake-RNAseq/joblogs/"
  ```
  * This step might take long, depending on the sample sizes.
  * If the command execution is interrupted, users need to rerun Step 3 to generate all results expected.
  
* **Step 4: Deactivate the environment as needed**
  
  ```
  conda deactivate
  ```
  
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
    * snakemake: 8.15.2
    * snakemake-executor-plugin-cluster-generic: 1.0.9
* Reference file directory: ```/mnt/storage/labs/sviswanathan/snakemake_RNAseq_2024/Human_genome_2024/```
    * Can be modified as needed
