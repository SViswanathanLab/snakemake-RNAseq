rule Align:
    input:
        index=STAR_index_path,
        gtf=gtf_path,
        reads=["data/{Sample}"+fq1_suffix, "data/{Sample}"+fq2_suffix],
    output:
        "results/STAR_results/{Sample}/{Sample}_Aligned.toTranscriptome.out.bam",
        "results/STAR_results/{Sample}/{Sample}_ReadsPerGene.out.tab",
        "results/STAR_results/{Sample}/{Sample}_Log.final.out",
    log:
        "logs/star/{Sample}.log"
    threads: 16
    params:
        threads=16,
    script:
        "../scripts/starAlign.sh"
