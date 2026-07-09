# https://unix.stackexchange.com/questions/16443/combine-text-files-column-wise

rule cbind:
    input:
        expand("results/cbf/{ID}.txt", ID=lookup(query="Species == '{species}'", within = reads, cols="BioSample"))
    output:
        "results/cbf_table_{species}.txt"
    shell:
        r"""
        paste -d' ' {input} > {output}
        """

rule counting_bloom_filter:
    input:
        "results/kmc/{ID}.txt"
    output:
        temp("results/cbf/{ID}.txt")
    params:
        array_size = config["array_size"],
        num_hash = config["num_hash"]
    conda:
        "../envs/cbf.yaml"
    shell:
        """
        if [ ! -d "results/cbf" ]; then
            mkdir -p results/cbf
        fi
        
        python scripts/counting_bloom_filter.py --input {input} --output {output} --array-size {params.array_size} --num-hash {params.num_hash}
        """

rule kmer_distances:
    input:
        "results/cbf_table_{species}.txt"
    output:
        "results/kmer_distances_{species}.txt"
    conda:
        "../envs/cbf.yaml"
    shell:
        "python scripts/kmer_distances.py --input {input} --output {output}"
