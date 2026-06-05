rule degenotate:
    input:
        fasta = "../config/linear_genomes/sequence/{ref}.fa",
        gff = "../config/linear_genomes/annotation/{ref}.gff",
    output:
        temp("degenotate_results/{ref}/degeneracy-all-sites.bed"),
        "degenotate_results/{ref}/cds-nt-longest.fa"
    conda:
        "../envs/degenotate.yaml"
    shell:
        """
        degenotate.py --overwrite -a {input.gff} -g {input.fasta} -o degenotate_results/{wildcards.ref}

        degenotate.py --overwrite -a {input.gff} -g {input.fasta} -l -o degenotate_results/{wildcards.ref}
        """

rule mk_test:
    input:
        gff = "../config/linear_genomes/annotation/{ref}.gff",
        fa = "../config/linear_genomes/sequence/{ref}.fa",
        vcf = "bcftools_concat_results/{ref}_sorted.vcf.gz"
    output:
        "mk_test_{ref}/mk.tsv"
    params:
        outgroup = config["outgroup"]
    conda:
        "../envs/degenotate.yaml"
    shell:
        "degenotate.py -a {input.gff} -g {input.fa} -v {input.vcf} -u $(echo {params.outgroup} | sed 's: :,:g' ) -o mk_test_{wildcards.ref} -sfs --overwrite"
