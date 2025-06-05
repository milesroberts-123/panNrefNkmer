rule bwa_index:
    input:
        "../config/linear_genomes/{ref}.fa",
    output:
        amb=temp("../config/linear_genomes/{ref}.fa.amb"),
        ann=temp("../config/linear_genomes/{ref}.fa.ann"),
        bwt=temp("../config/linear_genomes/{ref}.fa.bwt"),
        pac=temp("../config/linear_genomes/{ref}.fa.pac"),
        sa=temp("../config/linear_genomes/{ref}.fa.sa"),
    conda:
        "../envs/bwa.yaml"
    log:
        "logs/bwa_index/{ref}.log",
    shell:
        """
        # index reference
        bwa index {input} &>> {log}
        """
