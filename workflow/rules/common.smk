# Shared helpers for the k-mer branches (kmc.smk, kmc_nofilt.smk, kmc_kraken.smk).

def get_unique_biosample_from_species(wildcards):
    """BioSample IDs for every sample of the given species."""
    return reads[reads["Species"] == wildcards.species]["BioSample"].unique().tolist()


def get_species_from_biosample(wildcards):
    """Species name for the given BioSample ID."""
    return reads[reads["BioSample"] == wildcards.ID]["Species"].unique().tolist()


def get_paste_groups(wildcards):
    """Group indices for the grouped-paste rules."""
    n = len(get_unique_biosample_from_species(wildcards))
    return list(range((n + config["paste_group_size"] - 1) // config["paste_group_size"]))
