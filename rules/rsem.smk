rule rsem_reference:
    input:
        fa=fa_path,
        gtf=gtf_path,
    output:
        multiext("rsem_ref/human_gencode", ".chrlist", ".n2g.idx.fa", ".transcripts.fa", ".grp", ".seq", ".idx.fa", ".ti")
    log:
        "logs/rsem_ref.log"
    threads: 16
    params:
        threads=16,
    shell:
        """
        source /etc/profile.d/modules.sh

        exec 2> {log}

        module load /mnt/storage/spack/share/spack/modules/linux-ubuntu16.04-x86_64/rsem-1.3.1-gcc-5.4.0-ml3p6ok

        mkdir -p rsem_ref

        rsem-prepare-reference -p {params.threads} --gtf {input.gtf} {input.fa} rsem_ref/human_gencode
        """

rule rsem_quant:
    input:
        ref=multiext("rsem_ref/human_gencode", ".chrlist", ".n2g.idx.fa", ".transcripts.fa", ".grp", ".seq", ".idx.fa", ".ti"),
        counts="results/STAR_results/{Sample}/{Sample}_Aligned.toTranscriptome.out.bam",
    output:
        "results/rsem_results/{Sample}/{Sample}.genes.results",
        "results/rsem_results/{Sample}/{Sample}.isoforms.results",
    log:
        "logs/rsem_quant/{Sample}.log"
    threads: 16
    params:
        threads=16,
    script:
        "../scripts/rsemQuant.sh"

rule rsem_gene_TPM:
    input:
        gene_quants=expand("results/rsem_results/{Sample}/{Sample}.genes.results", Sample=Samples),
        annotation=annotation_path,
    output:
        "results/rsem_geneLevel_wide_TPM_Matrix.csv",
        "results/rsem_geneLevel_wide_TPM_Matrix.Rds",
    log:
        "logs/rsem_geneLevel_TPM_Matrix.log"
    script:
        "../scripts/rsem_gene_TPM.py"

rule rsem_isoform_TPM:
    input:
        isoform_quants=expand("results/rsem_results/{Sample}/{Sample}.isoforms.results", Sample=Samples),
        annotation=annotation_path,
    output:
        "results/rsem_isoformLevel_wide_TPM_Matrix.csv",
        "results/rsem_isoformLevel_wide_TPM_Matrix.Rds",
    log:
        "logs/rsem_isoformLevel_TPM_Matrix.log"
    script:
        "../scripts/rsem_isoform_TPM.py"
