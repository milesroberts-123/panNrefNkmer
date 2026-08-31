rule kmc_kraken:
    input:
        unclass1="results/kraken2_unclassified/{ID}_unclassified_1.fastq",
        unclass2="results/kraken2_unclassified/{ID}_unclassified_2.fastq"
    output:
        counts=temp("results/kraken/kmc/{ID}.txt"),
        tmp_pre=temp("results/kraken/kmc_db_{ID}.kmc_pre"),
        tmp_suf=temp("results/kraken/kmc_db_{ID}.kmc_suf"),
        tmp_sort=temp(expand("results/kraken/sorted_kmc_db_{{ID}}.{ext}", ext=["kmc_pre", "kmc_suf"]))
    conda:
        "../envs/kmc.yaml"
    params:
        mincount=config["mincount"],
        maxcount=config["maxcount"],
        k=config["k"],
    shell:
        """
        if [ -d "tmp_kmc_kraken_{wildcards.ID}" ]; then
            rm -r tmp_kmc_kraken_{wildcards.ID}
        fi

        mkdir -p results/kraken tmp_kmc_kraken_{wildcards.ID}

        kmc -sm -m25 -t{threads} -ci{params.mincount} -cs{params.maxcount} -k{params.k} \
            {input.unclass1} {input.unclass2} \
            results/kraken/kmc_db_{wildcards.ID} tmp_kmc_kraken_{wildcards.ID}

        kmc_tools -t{threads} transform results/kraken/kmc_db_{wildcards.ID} sort results/kraken/sorted_kmc_db_{wildcards.ID}

        kmc_tools -t{threads} transform results/kraken/sorted_kmc_db_{wildcards.ID} dump {output.counts}

        rm -r tmp_kmc_kraken_{wildcards.ID}
        """


rule kmc_combine_dbs_kraken:
    input:
        pre=expand("results/kraken/kmc_db_{ID}.kmc_pre", ID=get_unique_biosample_from_species),
        suf=expand("results/kraken/kmc_db_{ID}.kmc_suf", ID=get_unique_biosample_from_species)
    output:
        db=temp(expand("results/kraken/kmc_combine_dbs/{{species}}.{suffix}", suffix=["kmc_pre", "kmc_suf"])),
        complex=temp("results/kraken/kmc_combine_dbs/{species}.complex")
    conda:
        "../envs/kmc.yaml"
    shell:
        """
        mkdir -p results/kraken/kmc_combine_dbs

        {{
            echo "INPUT:"
            printf '%s\\n' {input.pre} | grep '\\.kmc_pre$' | sed 's/\\.kmc_pre$//' | awk '{{print "set" NR " = " $0 " -ci1"}}'
            echo "OUTPUT:"
            printf "results/kraken/kmc_combine_dbs/{wildcards.species} = "
            printf '%s\\n' {input.pre} | grep '\\.kmc_pre$' | sed 's/\\.kmc_pre$//' | awk '{{printf "%sset%d", (NR>1?" + ":""), NR}} END{{print ""}}'
            echo "OUTPUT_PARAMS:"
            echo "-cs10000000000"
        }} > {output.complex}

        kmc_tools -t{threads} complex {output.complex}
        """


rule dump_combined_kmers_kraken:
    input:
        expand("results/kraken/kmc_combine_dbs/{{species}}.{suffix}", suffix=["kmc_pre", "kmc_suf"]),
    output:
        temp("results/kraken/dump_combined_kmers/{species}.txt"),
    conda:
        "../envs/kmc.yaml"
    shell:
        """
        kmc_tools transform results/kraken/kmc_combine_dbs/{wildcards.species} dump {output}
        """


rule prejoin_kraken:
    input:
        comb=expand("results/kraken/dump_combined_kmers/{species}.txt", species=get_species_from_biosample),
        sample="results/kraken/kmc/{ID}.txt"
    output:
        temp("results/kraken/prejoin/{ID}.txt")
    shell:
        "join -t $'\t' -a1 -a2 -e '0' -o auto {input.comb} {input.sample} | cut -f 3 > {output}"


rule paste_kraken:
    input:
        kmer_list="results/kraken/dump_combined_kmers/{species}.txt",
        kmer_dumps=expand("results/kraken/prejoin/{ID}.txt", ID=get_unique_biosample_from_species)
    output:
        "results/kraken/paste/{species}.txt"
    shell:
        "paste <(cut -f 1 {input.kmer_list}) {input.kmer_dumps} > {output}"


rule pigz_kraken:
    input:
        "results/kraken/paste/{species}.txt"
    output:
        "results/kraken/paste/{species}.txt.gz"
    shell:
        "pigz -p {threads} {input}"


rule subset_kmer_table_kraken:
    input:
        "results/kraken/paste/{species}.txt.gz"
    output:
        lines=temp("results/kraken/subset/kmer_lines/{species}.txt"),
        table="results/kraken/subset/kmer_tables/{species}.txt"
    params:
        num_kmers_kept=config["num_kmers_kept"]
    shell:
        """
        TOTAL_LINES=$(zcat {input} | wc -l)

        shuf -i 1-"$TOTAL_LINES" -n {params.num_kmers_kept} | sort -n > {output.lines}

        zcat {input} | awk 'NR==FNR{{a[$1];next}} FNR in a' {output.lines} - | tr '\t' ' ' | cut --delimiter=" " -f 2- > {output.table}
        """


rule subset_kmer_distances_kraken:
    input:
        "results/kraken/subset/kmer_tables/{species}.txt"
    output:
        "results/kraken/subset/kmer_distances/{species}.txt"
    conda:
        "../envs/cbf.yaml"
    shell:
        "python scripts/kmer_distances.py --input {input} --output {output} --threads {threads}"


rule batch_per_species_kraken:
    input:
        "results/kraken/paste/{species}.txt.gz",
        "results/kraken/subset/kmer_distances/{species}.txt",
        expand("results/kraken2_reports/{ID}.txt", ID=lookup(query="Species == '{species}'", within=reads, cols="BioSample")),
        expand("results/fastp/post_trim/{ID}.json", ID=lookup(query="Species == '{species}'", within=reads, cols="Run"))
    output:
        touch("results/kraken/batch_tracker/{species}.done")
