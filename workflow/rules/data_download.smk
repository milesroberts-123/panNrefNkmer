rule sra:
    output:
        r1=temp("results/raw_reads/{ID}_1.fastq"),
        r2=temp("results/raw_reads/{ID}_2.fastq"),
        sra=temp("{ID}/{ID}.sra")
    conda:
        "../envs/sra.yaml"
    shell:
        """
        # create directory
        if [ ! -d "results/raw_reads" ]; then
            mkdir -p results/raw_reads
        fi

        prefetch --max-size 500G {wildcards.ID}

        # download data
        fasterq-dump --progress --temp /tmp --threads {threads} --outdir ./results/raw_reads --split-files --skip-technical ./{wildcards.ID} 
        """

rule datasets:
    output:
        "results/contam.fa"
    params:
        contams = config["ncbi_contams"]
    conda:
        "../envs/datasets.yaml"
    shell:
        r"""
        if [ -d "contam_genomes" ]; then
            rm -r contam_genomes/
        fi

        # get all contaminating genomes
        datasets download genome taxon {params.contams} --reference --dehydrated --filename contam.zip

        # unpack metadata
        unzip contam.zip -d contam_genomes

        # download genomes based on metadata
        datasets rehydrate --directory contam_genomes/

        # search directory for all genomes and copy them into one file
        find contam_genomes/ncbi_dataset/data/ -type f -name '*.fna' -exec cat {{}} \; > {output}
        """

rule biosample:
    input:
        pread1=expand("results/fastp/trimmed_paired_R1_{run}.fastq.gz", run = lookup(query="BioSample == '{bio}'", within=reads, cols="Run")),
        pread2=expand("results/fastp/trimmed_paired_R2_{run}.fastq.gz", run = lookup(query="BioSample == '{bio}'", within=reads, cols="Run")),
        uread1=expand("results/fastp/trimmed_unpaired_R1_{run}.fastq.gz", run = lookup(query="BioSample == '{bio}'", within=reads, cols="Run")),
        uread2=expand("results/fastp/trimmed_unpaired_R2_{run}.fastq.gz", run = lookup(query="BioSample == '{bio}'", within=reads, cols="Run")),
    output:
        pread1="results/biosample/{bio}_paired_R1.fastq.gz",
        pread2="results/biosample/{bio}_paired_R2.fastq.gz",
        uread1="results/biosample/{bio}_unpaired_R1.fastq.gz",
        uread2="results/biosample/{bio}_unpaired_R2.fastq.gz",
    params:
        no_inputs=lambda wildcards, input: len(input)
    shell:
        """
        if [ {params.no_inputs} -gt 4 ]; then
            zcat {input.pread1} | gzip > {output.pread1}
            zcat {input.pread2} | gzip > {output.pread2}
            zcat {input.uread1} | gzip > {output.uread1}
            zcat {input.uread2} | gzip > {output.uread2}
        else
            mv {input.pread1} {output.pread1}
            mv {input.pread2} {output.pread2}
            mv {input.uread1} {output.uread1}
            mv {input.uread2} {output.uread2}
        fi
        """
