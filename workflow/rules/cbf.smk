# https://unix.stackexchange.com/questions/16443/combine-text-files-column-wise

rule cbind:
    input:
        expand("cbf_results/{ID}.txt", ID=lookup(query="Species == '{species}'", within = reads, cols="BioSample"))
    output:
        "cbf_table_{species}.txt"
    shell:
        r"""
        paste -d' ' {input} > {output}
        """

rule counting_bloom_filter:
    input:
        "kmc_results/{ID}.txt"
    output:
        temp("cbf_results/{ID}.txt")
    params:
        array_size = config["array_size"],
        num_hash = config["num_hash"]
    conda:
        "../envs/cbf.yaml"
    shell:
        """
        if [ ! -d "cbf_results" ]; then
            mkdir cbf_results
        fi
        
        python scripts/counting_bloom_filter.py --input {input} --output {output} --array-size {params.array_size} --num-hash {params.num_hash}
        """

rule kmer_distances:
    input:
        "cbf_table_{species}.txt"
    output:
        "kmer_distances_{species}.txt"
    conda:
        "../envs/cbf.yaml"
    shell:
        "python scripts/kmer_distances.py --input {input} --output {output}"
