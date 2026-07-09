rule vg_index:
    input:
        "../config/pangenomes/{ref}.gfa"
    output:
        dist = "../config/pangenomes/{ref}.dist",
        gbz = "../config/pangenomes/{ref}.giraffe.gbz",
        min = "../config/pangenomes/{ref}.shortread.withzip.min",
        zip = "../config/pangenomes/{ref}.shortread.zipcodes"
    log:
        "logs/vg_index/{ref}.log"
    conda:
        "../envs/vg.yaml"
    shell:
        "vg autoindex -w giraffe -g {input} -p {wildcards.ref}"

rule vg_filter:
    input:
        gbz = "../config/pangenomes/{ref}.giraffe.gbz",
        gam = "results/vg_giraffe/{ID}_{ref}.gam"
    output:
        "results/vg_filter/{ID}_{ref}.gam"
    conda:
        "../envs/vg.yaml"
    shell:
        "vg filter {input.gam} -r 0.90 -fu -m 1 -q 15 -D 999 -x {input.gbz} > {output}"

rule vg_stats:
    input:
        "results/vg_filter/{ID}_{ref}.gam"
    output:
        "results/vg_stats/{ID}_{ref}.txt"
    conda:
        "../envs/vg.yaml"
    shell:
        """
        vg stats -a {input}  {output}
        """

rule vg_giraffe:
    input:
        read1 = "results/fastp/trimmed_paired_R1_{ID}.fastq.gz",
        read2 = "results/fastp/trimmed_paired_R1_{ID}.fastq.gz",
        dist = "../config/pangenomes/{ref}.dist",
        gbz = "../config/pangenomes/{ref}.giraffe.gbz",
        min = "../config/pangenomes/{ref}.shortread.withzip.min",
        zip = "../config/pangenomes/{ref}.shortread.zipcodes"
    output:
        "results/vg_giraffe/{ID}_{ref}.gam"
    conda:
        "../envs/vg.yaml"
    benchmark:
        "benchmarks/vg_giraffe/{ref}/{ID}.log"
    shell:
        """
        vg giraffe -p -t {threads} -Z {input.gbz} -m {input.min} -z {input.zip} -d {input.dist} -f {input.read1} -f {input.read2} > {output}
        """

rule vg_surject:
    input:
        gbz = "../config/pangenomes/{panref}.giraffe.gbz",
        gam = "results/vg_filter/{ID}_{panref}.gam",
        paths = "../config/pangenomes/paths/{linref}.txt"
    output:
        temp("results/vg_surject/{ID}_{panref}_{linref}.bam")
    conda:
        "../envs/vg.yaml"
    shell:
        "vg surject -x {input.gbz} -F {input.paths} -b {input.gam} > {output}"

rule vg_deconstruct:
    input:
        "results/vg_deconstruct/{chr}_{ref}.xg"
    output:
        "results/vg_deconstruct/{chr}_{ref}.vcf"
    conda:
        "../envs/vg.yaml"
    shell:
        "vg deconstruct -t {threads} -p {wildcards.ref} {input} > {output}"

rule vg_paths:
    input:
        "../config/pangenomes/{ref}.giraffe.gbz",
    output:
        "results/vg_paths/{ref}.fa"
    conda:
        "../envs/vg.yaml"
    shell:
        "vg paths -F -x {input} > {output}"
