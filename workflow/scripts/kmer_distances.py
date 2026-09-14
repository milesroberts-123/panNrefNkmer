"""Mean pairwise k-mer distance metrics (Bray-Curtis, cosine) over a count matrix."""

import os

# Single-threaded BLAS inside workers; parallelism comes from joblib threads.
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

import click
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
import gzip
from pathlib import Path
import random
import mmh3

#CHUNK_ROWS = 100000

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


import heapq

def add_heap(heap, hash_value, index_value, n):
    """
    Maintain a max-heap of size n containing the n smallest values seen.
    heap is a list used as a max-heap (store negated values).
    """
    if len(heap) < n:
        heapq.heappush(heap, -index_value)  # negate for max-heap behavior
    elif hash_value < -heap[0]:  # -heap[0] is the current max
        heapq.heapreplace(heap, -index_value)  # pop max, push new value
    return heap

@click.command(context_settings={"show_default": True})
@click.option("-i", "--input", required=True, help="Path to k-mer count matrix")
@click.option("-p", "--prefix", required=True, help="Path to prefix file")
@click.option("-t", "--threads", default=1, help="Number of threads")
@click.option("-c", "--chunk-size", default=100000, help="Number of rows per chunk (lower means less memory is required)")
@click.option("-v", "--print-freq", default=100, help="Print a message every v chunks to indicate progress.")
@click.option("--ignore-first/--no-ignore-first", default=True, help="Toggle whether to ignore first column of text file (e.g. k-mer sequence column)")
@click.option("--subset-columns", type=click.IntRange(min=2), default=None, help="Randomly subset to this many columns before calculating distances. Must be an integer.")
@click.option("--subset-rows", type=click.IntRange(min=2), default=None, help="Randomly subset to this many rows before calculating distances. --no-ignore-first must be set to allow for subsetting to by minhashing k-mers.")
@click.option("-s", "--seperator", default=" ", help="Separator between columns")


def main(input, prefix, threads, chunk_size, print_freq, ignore_first, subset_columns, subset_rows, seperator):
    """Average Bray-Curtis and cosine distances over all column pairs."""
    # Interpret escape sequences like '\t' (single quotes in bash pass the
    # literal two characters backslash+t) so both --seperator '\t' and
    # --seperator $'\t' work.
    seperator = seperator.encode().decode("unicode_escape")

    print("Determining number of columns in k-mer table...")

    if Path(input).suffix.lower() == ".gz":
        print("Detected .gz at end of file name. Assuming file is gzip compressed...")
        with gzip.open(input, 'rt') as x:
            ncols = len(x.readline().split(seperator))
    else:
        print("No .gz detected at end of file name. Assuming file is not gzip compressed...")
        with open(input, 'rt') as x:
            ncols = len(x.readline().split(seperator))

    if ncols < 2:
        raise ValueError(
            f"Detected {ncols} column(s) in {input} using separator "
            f"{seperator!r}. Check that the separator matches the file "
            f"(e.g. pass --seperator $'\\t' for tab-separated files)."
        )

    print(f"Number of columns is: {ncols}")

    if ignore_first:
        print("Ignoring first column...")
        cols_kept=range(1,ncols)
    else:
        print("Not ignoring first column...")
        cols_kept=range(0,ncols)

    if subset_columns is None:
        print("Not subsetting columns...")
    else:
        print(f"Subsetting to {subset_columns} columns...")
        cols_kept = random.sample(cols_kept, subset_columns)
        print(f"Kept columns are: {cols_kept}")

    # Min-hash k-mer subsetting if needed
    print("Loading k-mer count matrix by chunk...")
    if subset_rows is None:
        print("Not subsetting rows by min-hash...")
        reader = pd.read_csv(input, sep=seperator, header=None, chunksize=chunk_size, dtype=np.float64, usecols=cols_kept)
    elif subset_rows is not None and ignore_first:
        raise click.ClickException("If subset_rows is set, then ignore_first should not be set because the first column needs to be a k-mer string.")
    elif subset_rows is not None and not ignore_first:
        print("Pass 0: subseting columns based on min-hash...")
        rows_kept = []
        reader = pd.read_csv(input, sep=seperator, header=None, chunksize=1, usecols=[0])
        index_value = 0
        for chunk in reader:
            #print(chunk)
            kmer = str(chunk.items())
            #kmer = chunk.astype(str)
            hash_value = mmh3.hash(kmer, 1, signed=False)
            add_heap(rows_kept, hash_value, index_value, subset_rows)
            index_value += 1
        reader = pd.read_csv(input, sep=seperator, header=None, chunksize=chunk_size, dtype=np.float64, usecols=cols_kept, skiprows=lambda i: i != 0 and (i - 1) not in rows_kept)
        
    #print("Loading k-mer count matrix by chunk...")
    #reader = pd.read_csv(input, sep=seperator, header=None, chunksize=chunk_size, dtype=np.float64, usecols=cols_kept)

    print("Pass 1: column sums...")
    col_sums = None
    col_non_zero = None
    chunk_tracker = 0
    for chunk in reader:
        s = chunk.sum(axis=0).values
        b = chunk.astype(bool).sum(axis=0) 
        col_sums = s if col_sums is None else col_sums + s
        col_non_zero = b if col_non_zero is None else col_non_zero + b
        chunk_tracker += 1
        if chunk_tracker % print_freq == 0:
            print(f"{chunk_tracker*chunk_size} rows processed.")

    if col_sums is None:
        raise ValueError(f"Input file {input} is empty")
    
    print(f"Write summed k-mer counts per column to {prefix}" + "_sums.txt" + "...")
    with open(prefix + "_sums.txt", "w") as file:
        file.write(",".join(map(str, col_sums)))

    print(f"Write total unique k-mers per column to {prefix}" + "_non_zero.txt" + "...")
    with open(prefix + "_non_zero.txt", "w") as file:
        file.write(",".join(map(str, col_non_zero)))

    n = len(col_sums)
    num_pairs = n * (n - 1) // 2
    print(f"Number of columns: {n}")
    print(f"Number of column pairs: {num_pairs}")

    print("Pass 2: pairwise distances...")
    reader = pd.read_csv(input, sep=seperator, header=None, chunksize=chunk_size, dtype=np.float64, usecols=cols_kept)
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

    print(f"Write k-mer distances to {prefix}" + "_dist.txt" + "...")
    with open(prefix + "_dist.txt", "w") as file:
        file.write(",".join(final_result))

    print("Done!")


if __name__ == "__main__":
    main()
