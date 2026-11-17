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

        # check if prefetch directory exists already
        if [ -d "{wildcards.ID}" ]; then
            rm -rf "{wildcards.ID}"
            echo "Folder deleted successfully."
        fi

        # drop stale split files from earlier attempts: they are not declared
        # outputs, but leftover _3.fastq / partial _1.fastq would both inflate
        # the file count below and waste tens of GB per sample
        rm -f results/raw_reads/{wildcards.ID}_*.fastq

        prefetch --max-size 500G {wildcards.ID}

        # download data
        fasterq-dump --progress --temp /tmp --threads {threads} --outdir ./results/raw_reads --split-files --skip-technical ./{wildcards.ID}

        # fasterq-dump numbers the split files by SRA read-group index, not
        # contiguously. Runs whose first division is technical (barcode read)
        # therefore come out as _2/_3 instead of _1/_2 once --skip-technical
        # drops the leading division. Renumber the emitted biological mates to
        # _1/_2 so the fastp rule's paired-end inputs exist.
        n=$(ls -1 results/raw_reads/{wildcards.ID}_*.fastq 2>/dev/null | wc -l)
        if [ "$n" -ne 2 ]; then
            echo "ERROR: expected 2 read divisions for {wildcards.ID} after --skip-technical, found $n:" >&2
            ls -l results/raw_reads/{wildcards.ID}_*.fastq 2>/dev/null >&2
            exit 1
        fi
        first=$(ls -1 results/raw_reads/{wildcards.ID}_*.fastq | sort -V | sed -n 1p)
        second=$(ls -1 results/raw_reads/{wildcards.ID}_*.fastq | sort -V | sed -n 2p)
        if [ "$first" != "results/raw_reads/{wildcards.ID}_1.fastq" ]; then
            mv -f "$first" results/raw_reads/{wildcards.ID}_1.fastq
        fi
        if [ "$second" != "results/raw_reads/{wildcards.ID}_2.fastq" ]; then
            mv -f "$second" results/raw_reads/{wildcards.ID}_2.fastq
        fi
        echo "Split files for {wildcards.ID}: $(basename "$first") $(basename "$second") -> _1.fastq _2.fastq"

        # fail here, with a clear message, rather than as a MissingOutputException
        for f in results/raw_reads/{wildcards.ID}_1.fastq results/raw_reads/{wildcards.ID}_2.fastq; do
            if [ ! -s "$f" ]; then
                echo "ERROR: $f is missing or empty after fasterq-dump" >&2
                exit 1
            fi
        done
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
    #group: "trim"
    input:
        reads=expand("results/fastp/trimmed_{{pairing}}_{{read}}_{run}.fastq.gz", run=lookup(query="BioSample == '{bio}'", within=reads, cols="Run")),
    output:
        temp("results/biosample/{bio}_{pairing}_{read}.fastq.gz"),
    params:
        no_inputs=lambda wildcards, input: len(input.reads)
    shell:
        """
        if [ {params.no_inputs} -gt 1 ]; then
            cat {input.reads} > {output}
        else
            mv {input.reads} {output}
        fi
        """
