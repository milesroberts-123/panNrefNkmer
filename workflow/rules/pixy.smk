rule pixy:
    input:
        vcf = "results/bcftools_concat/{ref}_sorted.vcf.gz",
        tbi = "results/bcftools_concat/{ref}_sorted.vcf.gz.tbi",
        gff = "../config/linear_genomes/annotation/{ref}.gff",
    output:
        bed = temp("results/pixy/{ref}_genes.bed"),
        populations = temp("results/populations_{ref}.txt"),
        pixy = expand("results/pixy/{{ref}}_{stat}.txt",stat = ["pi", "watterson_theta", "tajima_d", "dxy", "fst"])
    params:
        ingroup = config["samples"],
        outgroup = config["outgroup"]
    conda:
        "../envs/pixy.yaml"
    shell:
        r"""
        # get bed intervals by gene
        awk '(($3 == "gene"))' {input.gff} | cut -f1,4,5 > {output.bed}
        
        echo {params.ingroup} | sed 's: :\n:g' | sed 's:$:\tingroup:g' > {output.populations}
        echo {params.outgroup} | sed 's: :\n:g' | sed 's:$:\toutgroup:g' >> {output.populations}

        # calculate statistics by gene
        pixy --populations {output.populations} --vcf {input.vcf} --bed_file {output.bed} --stats pi fst dxy watterson_theta tajima_d --output_folder results/pixy --output_prefix {wildcards.ref}
        """
