rule STARcount:
    input:
        counts=expand("results/STAR_results/{Sample}/{Sample}_ReadsPerGene.out.tab", Sample=Samples),
        annotation=annotation_path,
    output:
        "results/star_wide_countMatrix.csv",
        "results/star_wide_countMatrix.Rds",
    log:
        "logs/countMatrix.log"
    threads: 8
    script:
        "../scripts/countMatrix.py"
