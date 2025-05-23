rule fastqc:
    input:
        "data/{Read}.fastq.gz"
    output:
        "results/fastqc_results/{Read}_fastqc.zip",
        "results/fastqc_results/{Read}_fastqc.html",
    log:
        "logs/fastqc/{Read}.log"
    threads: 16
    params:
        threads=16,
    shell:
        """
        source /etc/profile.d/modules.sh

        exec 2> {log}

        module load /mnt/storage/spack/share/spack/modules/linux-ubuntu16.04-x86_64/fastqc-0.11.7-gcc-5.4.0-hjgfsw2

        mkdir -p results/fastqc_results

        fastqc -o results/fastqc_results/ -t 6 {input} --threads {params.threads}
        """
