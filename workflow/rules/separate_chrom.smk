rule separate_chrom:
    input:
        expand("../config/linear_genomes/sequence/{ref}.fa", ref = config["linrefs"])
    output:
        gz="{chr}.fa.gz",
        tbi="{chr}.fa.gz.tbi"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        echo Grab all instances of {wildcards.chr} 
        grep "{wildcards.chr}$" {input} > $(basename {output} .gz)

        echo Compress fasta file
        bgzip $(basename {output} .gz)

        echo Index fasta file
        tabix {output.gz}
        """
