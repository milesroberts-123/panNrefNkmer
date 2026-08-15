rule kmc_nofilt:
    input:
        pread1="results/biosample/{ID}_paired_R1.fastq.gz",
        pread2="results/biosample/{ID}_paired_R2.fastq.gz",
        uread1="results/biosample/{ID}_unpaired_R1.fastq.gz",
        uread2="results/biosample/{ID}_unpaired_R2.fastq.gz"
    output:
        counts=temp("results/nofilt/kmc/{ID}.txt"),
        tmp_pre=temp("results/nofilt/kmc_db_{ID}.kmc_pre"),
        tmp_suf=temp("results/nofilt/kmc_db_{ID}.kmc_suf"),
        tmp_sort=temp(expand("results/nofilt/sorted_kmc_db_{{ID}}.{ext}", ext=["kmc_pre", "kmc_suf"]))
    conda:
        "../envs/kmc.yaml"
    params:
        mincount=config["mincount"],
        maxcount=config["maxcount"],
        k=config["k"],
    shell:
        """
        if [ -d "tmp_kmc_nofilt_{wildcards.ID}" ]; then
            rm -r tmp_kmc_nofilt_{wildcards.ID}
        fi

        mkdir -p results/nofilt tmp_kmc_nofilt_{wildcards.ID}

        kmc -sm -m25 -t{threads} -ci{params.mincount} -cs{params.maxcount} -k{params.k} \
            {input.pread1} {input.pread2} {input.uread1} {input.uread2} \
            results/nofilt/kmc_db_{wildcards.ID} tmp_kmc_nofilt_{wildcards.ID}

        kmc_tools -t{threads} transform results/nofilt/kmc_db_{wildcards.ID} sort results/nofilt/sorted_kmc_db_{wildcards.ID}

        kmc_tools -t{threads} transform results/nofilt/sorted_kmc_db_{wildcards.ID} dump {output.counts}

        rm -r tmp_kmc_nofilt_{wildcards.ID}
        """


def get_unique_biosample_from_species_nofilt(wildcards):
    matched_rows = reads[reads["Species"] == wildcards.species]
    unique_outputs = matched_rows["BioSample"].unique().tolist()
    return unique_outputs


rule kmc_combine_dbs_nofilt:
    input:
        pre=expand("results/nofilt/kmc_db_{ID}.kmc_pre", ID=get_unique_biosample_from_species_nofilt),
        suf=expand("results/nofilt/kmc_db_{ID}.kmc_suf", ID=get_unique_biosample_from_species_nofilt)
    output:
        db=temp(expand("results/nofilt/kmc_combine_dbs/{{species}}.{suffix}", suffix=["kmc_pre", "kmc_suf"])),
        complex=temp("results/nofilt/kmc_combine_dbs/{species}.complex")
    conda:
        "../envs/kmc.yaml"
    params:
        prefix="results/nofilt/kmc_combine_dbs/{species}"
    shell:
        """
        mkdir -p results/nofilt/kmc_combine_dbs

        {{
            echo "INPUT:"
            printf '%s\\n' {input.pre} | grep '\\.kmc_pre$' | sed 's/\\.kmc_pre$//' | awk '{{print "set" NR " = " $0 " -ci1"}}'
            echo "OUTPUT:"
            printf "results/nofilt/kmc_combine_dbs/{wildcards.species} = "
            printf '%s\\n' {input.pre} | grep '\\.kmc_pre$' | sed 's/\\.kmc_pre$//' | awk '{{printf "%sset%d", (NR>1?" + ":""), NR}} END{{print ""}}'
            echo "OUTPUT_PARAMS:"
            echo "-cs10000000000"
        }} > {output.complex}

        kmc_tools -t{threads} complex {output.complex}
        """


rule dump_combined_kmers_nofilt:
    input:
        expand("results/nofilt/kmc_combine_dbs/{{species}}.{suffix}", suffix=["kmc_pre", "kmc_suf"]),
    output:
        temp("results/nofilt/dump_combined_kmers/{species}.txt"),
    conda:
        "../envs/kmc.yaml"
    shell:
        """
        kmc_tools transform results/nofilt/kmc_combine_dbs/{wildcards.species} dump {output}
        """


def get_species_from_biosample_nofilt(wildcards):
    matched_rows = reads[reads["BioSample"] == wildcards.ID]
    unique_outputs = matched_rows["Species"].unique().tolist()
    return unique_outputs


rule prejoin_nofilt:
    input:
        comb=expand("results/nofilt/dump_combined_kmers/{species}.txt", species=get_species_from_biosample_nofilt),
        sample="results/nofilt/kmc/{ID}.txt"
    output:
        temp("results/nofilt/prejoin/{ID}.txt")
    shell:
        "join -t $'\t' -a1 -a2 -e '0' -o auto {input.comb} {input.sample} | cut -f 3 > {output}"


rule paste_nofilt:
    input:
        kmer_list="results/nofilt/dump_combined_kmers/{species}.txt",
        kmer_dumps=expand("results/nofilt/prejoin/{ID}.txt", ID=get_unique_biosample_from_species_nofilt)
    output:
        "results/nofilt/paste/{species}.txt"
    shell:
        "paste <(cut -f 1 {input.kmer_list}) {input.kmer_dumps} > {output}"


rule pigz_nofilt:
    input:
        "results/nofilt/paste/{species}.txt"
    output:
        "results/nofilt/paste/{species}.txt.gz"
    shell:
        "pigz -p {threads} {input}"


rule counting_bloom_filter_nofilt:
    input:
        "results/nofilt/kmc/{ID}.txt"
    output:
        temp("results/nofilt/cbf/{ID}.txt")
    params:
        array_size=config["array_size"],
        num_hash=config["num_hash"]
    conda:
        "../envs/cbf.yaml"
    shell:
        """
        if [ ! -d "results/nofilt/cbf" ]; then
            mkdir -p results/nofilt/cbf
        fi

        python scripts/counting_bloom_filter.py --input {input} --output {output} --array-size {params.array_size} --num-hash {params.num_hash}
        """


rule cbind_nofilt:
    input:
        expand("results/nofilt/cbf/{ID}.txt", ID=lookup(query="Species == '{species}'", within=reads, cols="BioSample"))
    output:
        "results/nofilt/cbf_table/{species}.txt"
    shell:
        r"""
        paste -d' ' {input} > {output}
        """


rule kmer_distances_nofilt:
    input:
        "results/nofilt/cbf_table/{species}.txt"
    output:
        "results/nofilt/kmer_distances/{species}.txt"
    conda:
        "../envs/cbf.yaml"
    shell:
        "python scripts/kmer_distances.py --input {input} --output {output}"


rule subset_kmer_table_nofilt:
    input:
        "results/nofilt/paste/{species}.txt.gz"
    output:
        lines=temp("results/nofilt/subset/kmer_lines/{species}.txt"),
        table="results/nofilt/subset/kmer_tables/{species}.txt"
    params:
        num_kmers_kept=config["num_kmers_kept"]
    shell:
        """
        TOTAL_LINES=$(zcat {input} | wc -l)

        shuf -i 1-"$TOTAL_LINES" -n {params.num_kmers_kept} | sort -n > {output.lines}

        zcat {input} | awk 'NR==FNR{{a[$1];next}} FNR in a' {output.lines} - | tr '\t' ' ' | cut --delimiter=" " -f 2- > {output.table}
        """


rule subset_kmer_distances_nofilt:
    input:
        "results/nofilt/subset/kmer_tables/{species}.txt"
    output:
        "results/nofilt/subset/kmer_distances/{species}.txt"
    conda:
        "../envs/cbf.yaml"
    shell:
        "python scripts/kmer_distances.py --input {input} --output {output}"


rule batch_per_species_nofilt:
    input:
        "results/nofilt/kmer_distances/{species}.txt",
        "results/nofilt/paste/{species}.txt.gz",
        "results/nofilt/subset/kmer_distances/{species}.txt",
        expand("results/fastp/post_trim/{ID}.json", ID=lookup(query="Species == '{species}'", within=reads, cols="Run"))
    output:
        touch("results/nofilt/batch_tracker/{species}.done")
