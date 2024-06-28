rule salmon_quant:
    input:
        index=salmon_index_path,
        reads=["data/{Sample}"+fq1_suffix, "data/{Sample}"+fq2_suffix],
    output:
        "results/salmon_results/{Sample}/quant.sf"
    log:
        "logs/salmon_quant/{Sample}.log"
    threads: 30
    script:
        "../scripts/salmonQuant.sh"

rule salmon_TPM:
    input:
        quants=expand("results/salmon_results/{Sample}/quant.sf", Sample=Samples),
        annotation=annotation_path,
    output:
        "results/salmon_wide_TPM_Matrix.csv",
        "results/salmon_wide_TPM_Matrix.Rds",
    log:
        "logs/salmonTPM_Matrix.log"
    threads: 8
    script:
        "../scripts/salmonTPM.py"
