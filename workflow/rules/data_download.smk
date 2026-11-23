rule sra:
    output:
        r1=temp("results/raw_reads/{ID}_1.fastq"),
        r2=temp("results/raw_reads/{ID}_2.fastq"),
        sra=temp(directory("{ID}"))
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

        # NCBI serves the SRA Lite variant ({wildcards.ID}.sralite) for runs whose
        # original base qualities were never submitted, and the normalized
        # {wildcards.ID}.sra otherwise, so the cache is declared as the directory rather
        # than one specific filename. Report which one arrived: Lite files
        # carry reduced-precision qualities, which weakens quality-based
        # trimming downstream even though the sequences are identical.
        if [ -f "{wildcards.ID}/{wildcards.ID}.sra" ]; then
            sra_file="{wildcards.ID}/{wildcards.ID}.sra"
            echo "prefetch: normalized SRA (full qualities)"
        elif [ -f "{wildcards.ID}/{wildcards.ID}.sralite" ]; then
            sra_file="{wildcards.ID}/{wildcards.ID}.sralite"
            echo "prefetch: SRA Lite (reduced-precision qualities for {wildcards.ID})" >&2
        else
            echo "ERROR: prefetch produced no .sra/.sralite under {wildcards.ID}" >&2
            ls -l "{wildcards.ID}" >&2
            exit 1
        fi

        # download data
        #
        # The extraction chain, cheapest route first:
        #   1. fasterq-dump, multi-threaded; works for ordinary unaligned runs.
        #   2. fastq-dump; slower, handles some layouts fasterq-dump rejects.
        #   3. ENA FASTQ mirror; required for runs stored in alignment context
        #      (the ones whose prefetch pulls an extra
        #      ncbi-acc:...?vdb-ctx=refseq dependency). Those objects abort
        #      fasterq-dump with "double free or corruption" (4 threads) or
        #      "sorter.c producer_thread_func: rec.read.len = 0" (1 thread), and
        #      fastq-dump cannot open them at all ("cannot be opened as database
        #      or table"), but ENA serves the same reads as plain FASTQ.
        dump_ok=0
        if fasterq-dump --progress --temp /tmp --threads {threads} --outdir ./results/raw_reads --split-files --skip-technical "$sra_file"; then
            dump_ok=1
        else
            echo "fasterq-dump failed for {wildcards.ID}; retrying with fastq-dump (alignment-context run?)" >&2
            rm -f results/raw_reads/{wildcards.ID}_*.fastq
            if fastq-dump --outdir ./results/raw_reads --split-3 --skip-technical "$sra_file"; then
                dump_ok=1
            else
                # Neither SRA tool can read this object. Runs submitted in
                # alignment context are stored against a reference and are not
                # extractable by fasterq-dump/fastq-dump at all; the ENA mirror
                # serves the same reads as plain FASTQ, so take them from there.
                echo "both SRA tools failed for {wildcards.ID}; falling back to the ENA FASTQ mirror" >&2
                rm -f results/raw_reads/{wildcards.ID}_*.fastq
                ena_urls=$(curl -sf "https://www.ebi.ac.uk/ena/portal/api/filereport?accession={wildcards.ID}&result=read_run&fields=fastq_ftp&format=tsv" \
                    | tail -n +2 | cut -f1)
                if [ -z "$ena_urls" ]; then
                    echo "ERROR: ENA has no FASTQ for {wildcards.ID}; no extraction route left" >&2
                    exit 1
                fi
                i=1
                for u in $(echo "$ena_urls" | tr ';' ' '); do
                    echo "downloading https://$u -> results/raw_reads/{wildcards.ID}_$i.fastq" >&2
                    curl -sfL "https://$u" | gzip -dc > results/raw_reads/{wildcards.ID}_$i.fastq
                    i=$((i + 1))
                done
                # the cache directory is a declared output; create it so the
                # rule's outputs are complete even though SRA was never readable
                mkdir -p "{wildcards.ID}"
                dump_ok=1
            fi
        fi
        if [ "$dump_ok" -ne 1 ]; then
            exit 1
        fi

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
                echo "ERROR: $f is missing or empty after read extraction" >&2
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

        # fastp emits a 0-byte file for --unpaired1/--unpaired2 when a library
        # has no unpaired reads. Such a file satisfies Snakemake's existence
        # check but is not a valid gzip stream, so downstream KMC fails with
        # "Unexpected end of gzip file". Substitute a valid empty gzip.
        if [ ! -s {output} ]; then
            echo "no reads for {output}; writing an empty gzip" >&2
            gzip -c /dev/null > {output}
        fi
        """
