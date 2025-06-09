rule bwa_index:
    input:
        "../config/linear_genomes/sequence/{ref}.fa",
    output:
        amb=temp("../config/linear_genomes/sequence/{ref}.fa.amb"),
        ann=temp("../config/linear_genomes/sequence/{ref}.fa.ann"),
        bwt=temp("../config/linear_genomes/sequence/{ref}.fa.bwt"),
        pac=temp("../config/linear_genomes/sequence/{ref}.fa.pac"),
        sa=temp("../config/linear_genomes/sequence/{ref}.fa.sa"),
    conda:
        "../envs/bwa.yaml"
    log:
        "logs/bwa_index/{ref}.log",
    shell:
        """
        # index reference
        bwa index {input} &>> {log}
        """
