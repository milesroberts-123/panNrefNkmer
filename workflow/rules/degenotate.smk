rule degenotate:
    input:
        fasta = ,
        gff = 
    output:
    conda:
        "../envs/degenotate.yaml"
    shell:
        "degenotate.py -a {input.gff} -g {input.fasta} -o degen_results/{wildcards.ref}"