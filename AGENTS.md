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

- `jules_only` — legacy PSMC/ROH track; **not included in `all`**, must be run separately
- `jules_roh_only` — ROH outputs only, plus raw-read gzip verification of every sample
- `jules_verify_downloads_only` — audit-only: re-verify existing raw reads, produce no analysis outputs; run with `--keep-going`
- `jules_msmc2_only` — per-sample MSMC2 outputs; drafted but not yet proven end-to-end, run separately
- `kmers_only` — k-mer distance outputs only (no reference needed)
- `kmers_subset_only` — subset k-mer distance outputs only
- `kmers_nofilt_only` — k-mer branch with no reference contamination removal (outputs under `results/nofilt/`)
- `kmers_kraken_only` — k-mer branch with kraken2 microbial screening (outputs under `results/kraken/`)
- `batch_per_species` — per-species aggregation (k-mer distances, paste table, kmc histo, fastp JSONs)

There is no `all` target currently — the linear-reference, pangenome, and population-genetics tracks are disabled (see below).

## Disabled tracks are archived

Rules for the linear-reference (`bwa`, `samtools`), pangenome (`vg`, `pggb`), and population-genetics (`salmon`, `degenotate`, `pixy`, `split`, `mash`, `multiqc`, `change_headers`) tracks are disabled. They live in `workflow/archives/rules/` (with their envs in `workflow/archives/envs/` and the unused `call_fastdfe.py` in `workflow/archives/scripts/`). To re-enable a track, move the file back to `workflow/rules/` and add its `include:` line to the Snakefile.

`workflow/rules/msmc2.smk` (pangenome-keyed MSMC2) also sits in `rules/` but is **not included** — it reads config keys `vg_ref_backbone`/`vg_ref_variant` that do not exist in `config.yaml`. Enabling it requires both the `include:` line and those config keys.

## Split numbering starts at 10

The `{split}` wildcard uses `--numeric-suffixes=10` in `archives/rules/split.smk` and `range(10, 10 + config["splits"])` in the Snakefile. Split indices are 10, 11, 12, ... — **not 0-based or 1-based**.

## Three independent tracks

The pipeline has three largely non-overlapping tracks:
1. **Linear-reference** (`bwa` → `samtools`/`bcftools`) — keyed by `{ref}` (archived)
2. **Pangenome** (`vg`/`pggb`) — keyed by `{panref}`/`{chr}` (archived)
3. **Jules** (legacy) — `jules_v2.smk` + `jules_msmc2.smk`, uses its own `reference_genome_path` config key and `envs/msmc2.yaml`

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
- `config/samples_large.tsv` — full sample sheet driving the k-mer branches (read as `reads`); columns: `Run`, `BioSample`, `Species`
- `config/samples_medium.tsv` — chromosome-level-assembly subset driving the jules_* targets (read as `jules_reads`); same columns
  - The `Group` (ingroup/outgroup for MK test) and `LibraryType` (salmon) columns referenced by the archived rules are not present in either samples file
- `reference_genome_path` — used by the jules track only; points at `{species}/{species}.fasta` directories
- The k-mer reference contamination rules (`kmc_ref_db`, `kmc_rm_contam`) read from `../config/reference_genomes/{species}/{species}.fasta` instead
- Contaminant genomes: NCBI taxa in `config["ncbi_contams"]` + custom fastas in `config/contaminants/`

## No lint or test suite

`tests/` contains small fixture files (fastq/fasta/vg) for manual use. There is no automated test runner, linter, or typechecker.

## `localrules` bypass slurm

`all` and `batch_per_sample` are declared as `localrules` in the Snakefile — they run on the submit node regardless of profile.

## Adding rules: update the slurm profile

Every new cluster rule needs `set-resources` and `set-threads` entries in `workflow/profiles/default/config.yaml`, mirroring the analogous existing rule. Rules without entries silently fall back to `default-resources` (30m runtime, 1 cpu), which will kill long jobs. Also check `workflow/profiles/icer/config.yaml` if the rule will run on ICER.
