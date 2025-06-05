rule bcftools_call:
    shell:
        "bcftools mpileup --no-reference -Ou {input} | bcftools call -mv > {output}"
