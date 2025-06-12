rule bwa_index:
    input:
        "../config/linear_genomes/sequence/{ref}.fa",
    output:
        amb="../config/linear_genomes/sequence/{ref}.fa.amb",
        ann="../config/linear_genomes/sequence/{ref}.fa.ann",
        bwt="../config/linear_genomes/sequence/{ref}.fa.bwt",
        pac="../config/linear_genomes/sequence/{ref}.fa.pac",
        sa="../config/linear_genomes/sequence/{ref}.fa.sa",
    conda:
        "../envs/bwa.yaml"
    log:
        "logs/bwa_index/{ref}.log",
    shell:
        """
        # index reference
        bwa index {input} &>> {log}
        """
