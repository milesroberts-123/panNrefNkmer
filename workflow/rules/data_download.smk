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
        reads=expand("results/fastp/trimmed_{{pairing}}_{{read}}_{run}.fastq.gz", run=lookup(query="BioSample == '{bio}'", within=reads, cols="Run")),
    output:
        "results/biosample/{bio}_{pairing}_{read}.fastq.gz",
    params:
        no_inputs=lambda wildcards, input: len(input.reads)
    shell:
        """
        if [ {params.no_inputs} -gt 1 ]; then
            zcat {input.reads} | gzip > {output}
        else
            mv {input.reads} {output}
        fi
        """
