rule degenotate:
    input:
        fasta = "../config/linear_genomes/sequence/{ref}.fa",
        gff = "../config/linear_genomes/annotation/{ref}.gff",
    output:
        "degenotate_results/{ref}/degeneracy-all-sites.bed",
        #"degenotate_results/{ref}/cds-nt-longest.fa"
    conda:
        "../envs/degenotate.yaml"
    shell:
        "degenotate.py --overwrite -a {input.gff} -g {input.fasta} -o degenotate_results/{wildcards.ref}"
