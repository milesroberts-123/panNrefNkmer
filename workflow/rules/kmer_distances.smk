rule kmer_distances:
    input:
        "cbf_table_{species}.txt"
    output:
        "kmer_distances_{species}.txt"
    conda:
        "../envs/cbf.yaml"
    shell:
        "python scripts/kmer_distances.py --input {input} --output {output}"
