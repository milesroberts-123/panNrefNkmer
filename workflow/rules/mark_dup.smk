# source:
# https://www.htslib.org/algorithms/duplicate.html

rule mark_dup:
    input:
        "bwa_results/{ID}_{ref}.bam"
    output:
        temp("mark_dup_results/{ID}_{ref}.bam")
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools collate -@ {threads} -O -u {input} | samtools fixmate -@ {threads} -m -u - - | samtools sort -@ {threads} -u - | samtools markdup -@ {threads} - {output}"
