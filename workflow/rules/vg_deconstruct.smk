rule vg_deconstruct:
    input:
        "{chr}_{ref}.xg"
    output:
        "{chr}_{ref}.vcf"
    conda:
        "../envs/vg.yaml"
    shell:
        "vg deconstruct -t {threads} -p {wildcards.ref} {input} > {output}"
