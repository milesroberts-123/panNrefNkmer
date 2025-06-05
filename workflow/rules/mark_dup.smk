# source:
# https://www.htslib.org/algorithms/duplicate.html

rule mark_dup:
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools collate -@ 4 -O -u example.bam | samtools fixmate -@ 4 -m -u - - | samtools sort -@ 4 -u - | samtools markdup -@ 4 - markdup.bam"
