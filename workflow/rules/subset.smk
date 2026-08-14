rule subset_kmer_table:
    input:
        "results/paste/{species}.txt.gz" 
    output:
        lines=temp("results/subset/kmer_lines/{species}.txt"),
        table="results/subset/kmer_tables/{species}.txt"
    params:
        num_kmers_kept = config["num_kmers_kept"]
    shell:
        """
        # 1. Get the total line count of the gzipped file
        TOTAL_LINES=$(zcat {input} | wc -l)

        # 2. Extract exactly 10000 random line numbers, sorted numerically
        shuf -i 1-"$TOTAL_LINES" -n {params.num_kmers_kept} | sort -n > {output.lines}

        # 3. Stream the file once and pull those exact matching rows
        zcat {input} | awk 'NR==FNR{{a[$1];next}} FNR in a' {output.lines} - | tr '\t' ' ' | cut --delimiter=" " -f 2- > {output.table}
        """

rule subset_kmer_distances:
    input:
        "results/subset/kmer_tables/{species}.txt"
    output:
        "results/subset/kmer_distances/{species}.txt"
    conda:
        "../envs/cbf.yaml"
    shell:
        "python scripts/kmer_distances.py --input {input} --output {output}"
