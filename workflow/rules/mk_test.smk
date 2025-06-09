rule mk_test:
    input:
        gff =,
        fa =,
        vcf =
    output:
    conda:
        "../envs/degenotate.yaml"
    shell:
        "degenotate.py -a [annotation file] -g [genome fasta file] -v [vcf file] -u [sample ID(s) of outgroup in VCF file] -o [output directory] -sfs"
