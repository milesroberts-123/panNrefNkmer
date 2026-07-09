rule multiqc:
    input:
        expand("results/fastp/{ID}_R1R2.json", ID = config["samples"] + config["outgroup"] + config["rna"]),
        expand("results/fastp/{ID}_U1.json", ID = config["samples"] + config["outgroup"] + config["rna"]),
        expand("results/fastp/{ID}_U2.json", ID = config["samples"] + config["outgroup"] + config["rna"]),
        expand("results/fastp/no_contam_{ID}.json", ID = config["samples"])
    output:
        "results/multiqc/multiqc_report.html"
    conda:
        "../envs/multiqc.yaml"
    shell:
        "multiqc results/fastp/ -o results/multiqc"
