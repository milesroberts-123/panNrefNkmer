rule counting_bloom_filter:
    input:
        "kmc_results/{ID}.txt"
    output:
        "cbf_results/{ID}.txt"
    params:
        array_size = config["array_size"],
        num_hash = config["num_hash"]
    conda:
        "../envs/cbf.smk"
    shell:
        "scripts/counting_bloom_filter.py --input {input} --output {output} --array-size {params.array_size} --num-hash {params.num_hash}"

