# AGENTS.md

## Must run from `workflow/`

All paths in the Snakefile and rules are relative to `workflow/` (e.g. `../config/config.yaml`). Always `cd workflow` first.

## Environment

```bash
mamba env create --name snakemake --file setup.yaml
```

Per-rule conda envs live in `workflow/envs/`. The `--sdm conda` flag is required to use them.

## Common commands

```bash
# full pipeline on slurm
snakemake --sdm conda --rerun-incomplete --rerun-triggers mtime --scheduler greedy --retries 1 --keep-going

# local execution (testing)
snakemake --sdm conda --workflow-profile local

# jules_only (legacy per-sample track) — batched, NOT part of `all`
num_batch=50
for i in $(seq 1 $num_batch); do
  snakemake --sdm conda --rerun-triggers mtime --scheduler greedy --rerun-incomplete --batch jules_only=$i/$num_batch jules_only
done
```

Profiles in `workflow/profiles/`: `default` (slurm), `icer` (MSU ICER), `local` (no slurm).

## Target rules

- `all` (default) — samtools/bcftools stats, MK test, salmon, pixy, k-mer distances, multiqc
- `jules_only` — legacy PSMC/ROH track; **not included in `all`**, must be run separately
- `kmers_only` — k-mer distance outputs only (no reference needed)
- `batch_per_sample` — per-BioSample VCF splits, counting bloom filter, salmon quant

## Linear references are auto-detected

The Snakefile discovers linear references from the filesystem:
```python
linrefs = [os.path.splitext(os.path.basename(f))[0]
           for f in glob.glob("../config/linear_genomes/sequence/*.fa")]
```

There is **no `linrefs` key in config.yaml**. To add a reference, place `{ref}.fa` in `config/linear_genomes/sequence/` and `{ref}.gff` in `config/linear_genomes/annotation/`.

## Split numbering starts at 10

The `{split}` wildcard uses `--numeric-suffixes=10` in `split.smk` and `range(10, 10 + config["splits"])` in the Snakefile. Split indices are 10, 11, 12, ... — **not 0-based or 1-based**.

## Three independent tracks

The pipeline has three largely non-overlapping tracks:
1. **Linear-reference** (`bwa` → `samtools`/`bcftools`) — keyed by `{ref}`
2. **Pangenome** (`vg`/`pggb`) — keyed by `{panref}`/`{chr}`
3. **Jules** (legacy) — `jules.smk`, uses its own `reference_genome_path` config key

Rules from different tracks generally don't share intermediates.

## Wildcard conventions

- `{ID}` — SRA Run accession (per-run rules) or BioSample (per-sample rules); context-dependent
- `{ref}` — linear reference name (auto-detected from `config/linear_genomes/sequence/`)
- `{panref}` — pangenome name
- `{chr}` — chromosome
- `{split}` — genome partition index (starts at 10)
- `{species}` — species name (for k-mer/contaminant DBs)

## Configuration

- `config/config.yaml` — pipeline parameters; schema at `config/config.schema.yaml`
- `config/samples.tsv` — columns: `Run`, `BioSample`, `Species`, `LibraryType`, `Group`
  - `Group` column: used for ingroup/outgroup in MK test (`degenotate.smk`)
  - `LibraryType`: `rna` values are used for salmon quantification
- `reference_genome_path` — only used by the jules track, points at `{species}/{species}.fasta` directories
- Contaminant genomes: NCBI taxa in `config["ncbi_contams"]` + custom fastas in `config/contaminants/`

## No lint or test suite

`tests/` contains small fixture files (fastq/fasta/vg) for manual use. There is no automated test runner, linter, or typechecker.

## `localrules` bypass slurm

`all` and `batch_per_sample` are declared as `localrules` in the Snakefile — they run on the submit node regardless of profile.
