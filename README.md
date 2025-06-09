# vg-kmer-ref-snakes

A snakemake workflow to analyze short reads with a either pangenome reference, linear reference, and/or no reference

# To-do

- [x] add mitochondrial and chloroplast genomes as contaminants

- [x] figure out vg surject

- [x] add kmer counts

- [x] add genotype calling for linear references

- [ ] sort bam files

- [ ] find read duplicates with samtools

- [ ] add MK test in degenotate?

- [ ] counting bloom filter

- [ ] k-mer distances

- [ ] remove k-mers from cds sequences?

- [ ] create linear references for pangenome snp calling depending on path: Use `vg paths -F` to get fasta sequences for paths, then use that as the linear reference for genotype calling
