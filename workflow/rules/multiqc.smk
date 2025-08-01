rule multiqc:
    input:
        expand("fastp_results/{ID}_R1R2.json", ID = config["samples"] + config["outgroup"] + config["rna"]),
        expand("fastp_results/{ID}_U1.json", ID = config["samples"] + config["outgroup"] + config["rna"]),
        expand("fastp_results/{ID}_U2.json", ID = config["samples"] + config["outgroup"] + config["rna"]),
        expand("fastp_results/no_contam_{ID}.json", ID = config["samples"])
    output:
        "multiqc_report.html"
    conda:
        "../envs/multiqc.yaml"
    shell:
        "multiqc fastp_results/"
