rule degenotate:
    input:
        fasta = "../config/linear_genomes/sequence/{ref}.fa",
        gff = "../config/linear_genomes/annotation/{ref}.gff",
    output:
        "degenotate_results/{ref}/degeneracy-all-sites.bed"
    conda:
        "../envs/degenotate.yaml"
    shell:
        "degenotate.py -a {input.gff} -g {input.fasta} -l -o degenotate_results/{wildcards.ref}"
