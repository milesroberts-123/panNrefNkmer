rule mk_test:
    input:
        gff = "../config/linear_genomes/annotation/{ref}.gff",
        fa = "../config/linear_genomes/sequence/{ref}.fa",
        vcf = "bcftools_concat_results/{ref}.vcf.gz"
    output:
        "mk_test_{ref}/mk.tsv"
    params:
        outgroup = config["outgroup"]
    conda:
        "../envs/degenotate.yaml"
    shell:
        "degenotate.py -a {input.gff} -g {input.fa} -v {input.vcf} -u {params.outgroup} -o mk_test_{wildcards.ref} -sfs"
