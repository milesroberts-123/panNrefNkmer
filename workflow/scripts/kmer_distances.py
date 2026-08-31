"""Mean pairwise k-mer distance metrics (Bray-Curtis, cosine) over a count matrix."""

import os

# Single-threaded BLAS inside workers; parallelism comes from joblib threads.
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

import click
import numpy as np
import pandas as pd
from joblib import Parallel, delayed

CHUNK_ROWS = 100000


def process_chunk(values, col_sums):
    """Partial pairwise sums for one row-chunk of the count matrix.

    Returns (bc_sum, dot, sq): the sum of Bray-Curtis distances over all
    sample pairs, the pairwise dot-product matrix, and per-sample squared
    norms, each restricted to this chunk's rows.
    """
    Y = values.T  # n_samples x n_rows
    with np.errstate(divide="ignore", invalid="ignore"):
        Yn = Y / col_sums[:, None]
    # For L1-normalized columns, Bray-Curtis equals half the L1 distance.
    # Sum of pairwise L1 distances via per-column sort: each sorted value
    # s_k contributes (2k - n + 1) * s_k.
    n = Yn.shape[0]
    k = np.arange(n)
    bc = (np.sort(Yn, axis=0) * (2 * k - n + 1)[:, None]).sum() / 2
    dot = Y @ Y.T
    sq = (Y * Y).sum(axis=1)
    return bc, dot, sq


@click.command(context_settings={"show_default": True})
@click.option("-i", "--input", required=True, help="Path to k-mer count matrix")
@click.option("-o", "--output", required=True, help="Path to output file")
@click.option("-t", "--threads", default=1, help="Number of threads")
def main(input, output, threads):
    """Average Bray-Curtis and cosine distances over all column pairs."""
    print("Loading k-mer count matrix...")
    reader = pd.read_csv(input, sep=" ", header=None, chunksize=CHUNK_ROWS, dtype=np.float64)

    print("Pass 1: column sums...")
    col_sums = None
    for chunk in reader:
        s = chunk.sum(axis=0).values
        col_sums = s if col_sums is None else col_sums + s

    if col_sums is None:
        raise ValueError(f"Input file {input} is empty")

    n = len(col_sums)
    num_pairs = n * (n - 1) // 2
    print(f"Number of columns: {n}")
    print(f"Number of column pairs: {num_pairs}")

    print("Pass 2: pairwise distances...")
    reader = pd.read_csv(input, sep=" ", header=None, chunksize=CHUNK_ROWS, dtype=np.float64)
    results = Parallel(n_jobs=threads, prefer="threads", pre_dispatch=threads)(
        delayed(process_chunk)(chunk.values, col_sums) for chunk in reader
    )

    bc_total = 0.0
    dot_total = None
    sq_total = None
    for bc, dot, sq in results:
        bc_total += bc
        dot_total = dot if dot_total is None else dot_total + dot
        sq_total = sq if sq_total is None else sq_total + sq

    norms = np.sqrt(sq_total)
    with np.errstate(divide="ignore", invalid="ignore"):
        cos_matrix = 1 - dot_total / np.outer(norms, norms)
    zero = norms == 0
    cos_matrix[zero, :] = 1.0
    cos_matrix[:, zero] = 1.0
    cos_total = (cos_matrix.sum() - np.trace(cos_matrix)) / 2

    final_result = [str(bc_total / num_pairs), str(cos_total / num_pairs)]

    print(f"Write result to {output}...")
    with open(output, "w") as file:
        file.write(",".join(final_result))

    print("Done!")


if __name__ == "__main__":
    main()
