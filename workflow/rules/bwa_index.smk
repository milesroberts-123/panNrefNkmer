rule bwa_index:
    input:
        "reference_genomes/{ref}.fasta",
    output:
        amb=temp("reference_genomes/{ref}.fasta.amb"),
        ann=temp("reference_genomes/{ref}.fasta.ann"),
        bwt=temp("reference_genomes/{ref}.fasta.bwt"),
        pac=temp("reference_genomes/{ref}.fasta.pac"),
        sa=temp("reference_genomes/{ref}.fasta.sa"),
    conda:
        "../envs/bwa.yaml"
    log:
        "logs/bwa_index/{ID}.log",
    shell:
        """
        # index reference
        bwa index {input} &>> {log}
        """
