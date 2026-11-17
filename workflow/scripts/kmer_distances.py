"""k-mer count matrix utilities.

Two modes over a K x N count matrix (rows = k-mers, columns = samples):

1. Distance mode (default): streams the matrix and reports average pairwise
   Bray-Curtis and cosine distances over the sample columns, with optional
   min-hash row subsetting and cumulative per-column sum series.

2. Spectrum mode (--spectrum-out): streams the matrix once and computes
   incidence spectra h(k) (number of k-mers present in exactly k samples)
   at each depth standardization point and count threshold, plus the
   empirical rarefaction curve and neutral-theory (Ewens) comparison
   curve derived from the unthinned spectrum.
"""

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
    Maintain a max-heap of size n containing the n smallest packed keys seen.

    Key packing: (hash << 32) | row_index, so heap order breaks hash ties by
    row order and each entry knows which row it came from. Stored negated to
    emulate a max-heap with heapq (a min-heap).
    """
    key = (hash_value << 32) | index_value
    if len(heap) < n:
        heapq.heappush(heap, -key)
    elif key < -heap[0]:  # -heap[0] is the current max packed key
        heapq.heapreplace(heap, -key)
    return heap

def open_text(path):
    """Open path for reading text lines, transparently handling gzip."""
    if Path(path).suffix.lower() == ".gz":
        return gzip.open(path, "rt")
    return open(path, "rt")

M_GRID = (1, 2, 5, 10, 25, 50, 100, 150, 200)


def count_lines(path):
    """Number of rows in a text or gzip table file."""
    n = 0
    with open_text(path) as x:
        for _ in x:
            n += 1
    return n


def parse_grid(text, cast):
    """Parse a comma-separated CLI grid like '0,0.5,1'."""
    values = [cast(v) for v in str(text).split(",") if v.strip()]
    if not values:
        raise ValueError(f"empty grid: {text!r}")
    return values


def fmt_num(x):
    """Filename-safe label for a numeric grid value (0 not 0.0)."""
    xi = int(x)
    return str(xi) if x == xi else str(x)


def fit_power_law(ms, us):
    """Least-squares log-log fit u = c*m^z; returns (z, r2)."""
    if len(ms) < 2 or any(u <= 0 for u in us):
        return float("nan"), float("nan")
    lm = np.log(np.asarray(ms, dtype=float))
    lu = np.log(np.asarray(us, dtype=float))
    if lu.var() == 0:
        return float("nan"), float("nan")
    z, logc = np.polyfit(lm, lu, 1)
    resid = lu - (z * lm + logc)
    return float(z), float(1 - resid.var() / lu.var())


def rarefy_disjointprob(N, m_grid, h):
    """Sampling-based rarefaction of an incidence histogram.

    h: array of length N+1, h[k] = number of k-mers present in exactly k
    samples. Returns a list of (m, U(m)) where
    U(m) = sum_k h[k] * (1 - C(N-m, k)/C(N, k)) is the expected number of
    distinct k-mers present in a uniformly random subset of m samples
    (Heck et al. 1975 sample-based rarefaction).
    """
    ks = np.arange(1, N + 1, dtype=float)
    hv = h[1:].astype(float)
    points = []
    seen = set()
    for m in m_grid:
        m = int(min(m, N))
        if m in seen:
            continue
        seen.add(m)
        valid = ks <= (N - m)
        # P(disjoint) = prod_{t=0}^{m-1} (N-k-t)/(N-t), stable log form;
        # k > N-m is impossible to miss, so those terms are exactly 1.
        logp = np.zeros(len(ks))
        with np.errstate(divide="ignore"):
            for t in range(m):
                logp[valid] += np.log((N - ks[valid] - t) / (N - t))
        Pd = np.exp(logp)
        Pd[~valid] = 0.0
        points.append((m, float((hv * (1.0 - Pd)).sum())))
    if len(seen) == 0 or max(int(min(m, N)) for m in m_grid) < N:
        # endpoint: the full union (u(N) is the primary cross-species scalar)
        points.append((N, float(hv.sum())))
    return points


def ewens_theta_from_total(N, k_total):
    """Neutral infinite-sites theta hat matching the observed union size.

    Under Ewens with theta, the expected number of segregating sites in a
    sample of N is theta * H_N, so theta_hat = k_total / H_N.
    """
    H = sum(1.0 / i for i in range(1, N + 1))
    return k_total / H


def run_spectrum(input, spectrum_out, depth_totals, genome_size, depth_grid,
                 count_grid, seed, subset_rows_count, union_rows, chunk_size,
                 print_freq):
    """Spectrum mode: incidence spectra with depth standardization and count
    thresholds, plus rarefaction and Ewens comparison curves.

    Per-sample depth is estimated from the per-sample k-mer counts in
    depth_totals (computed over this same table), scaled by the ratio of
    union-table rows to subset-table rows; the thinning probability for a
    target depth d (in units of genome_size) is p_j = min(1, d*G/D_j).
    """
    print("Spectrum mode: computing k-mer incidence spectra...")

    if depth_totals is None:
        raise click.ClickException("--spectrum-out requires --depth-totals (per-sample sums from a distance run over the same table).")

    rows_subset = subset_rows_count if subset_rows_count is not None else count_lines(input)
    union = union_rows if union_rows is not None else rows_subset
    print(f"Subset table rows: {rows_subset}; union rows: {union}")

    col_sums = np.array([float(v) for v in open(depth_totals).read().strip().split(",")])
    depth = col_sums * (union / float(rows_subset))

    # sanity + layout: match the sums length against the table's column count
    # to auto-detect whether column 0 is a k-mer string column (paste tables)
    # or a sample column (subset tables, which have the k-mer column stripped)
    with open_text(input) as x:
        header_cols = len(x.readline().rstrip("\n").split())
    N = len(col_sums)
    if N == header_cols - 1:
        first_sample_col = 1
        print("First column looks like a k-mer string column; excluding it.")
    elif N == header_cols:
        first_sample_col = 0
        print("Table has no k-mer column (counts only); all columns are samples.")
    else:
        raise ValueError(
            f"depth_totals has {N} columns but the table has {header_cols} columns; "
            "the sums must come from a distance run over this same table."
        )
    if N < 2:
        raise ValueError(f"Spectrum mode needs at least 2 sample columns; found {N}.")

    print(f"Samples: {N}; genome size: {genome_size}")
    print(f"Per-sample depth estimates (in k-mer counts): min {depth.min():.3g}, max {depth.max():.3g}")
    print(f"Depth grid (in G-units): {depth_grid}; count grid: {count_grid}")
    if 0.0 not in depth_grid:
        raise click.ClickException("--depth-grid must include 0 (the unthinned spectrum, which rarefaction and the Ewens comparison are computed from).")

    # one histogram per (depth, count-threshold) pair; hists[d] is the c=1
    # (plain presence) spectrum, written to its own file
    ks_max = N + 1
    hists = {}
    for d in depth_grid:
        key = fmt_num(d)
        hists[key] = np.zeros(ks_max, dtype=np.int64)
        for c in count_grid:
            if c == 1:
                continue
            hists[(key, c)] = np.zeros(ks_max, dtype=np.int64)

    rng = np.random.default_rng(seed)
    chunk_tracker = 0
    reader = pd.read_csv(input, sep=" ", header=None, chunksize=chunk_size,
                         usecols=list(range(first_sample_col, first_sample_col + N)), dtype=np.float64)
    h0_total = 0
    rows_seen = 0
    for chunk in reader:
        vals = np.nan_to_num(chunk.values, nan=0.0)
        # only draw thinning variates for d > 0, in sorted order so the
        # seeded stream is deterministic for a given grid
        thinned = {}
        for d in sorted(depth_grid):
            if d == 0:
                continue
            p = np.clip(d * genome_size / depth, 0.0, 1.0)
            if np.all(p >= 1.0):
                thinned[d] = vals.astype(np.int64)
            else:
                thinned[d] = rng.binomial(vals.astype(np.int64), p[None, :] * np.ones((vals.shape[0], 1)))
        for d in depth_grid:
            base = vals if d == 0 else thinned[d]
            key = fmt_num(d)
            k = (base > 0).sum(axis=1)
            h0_total += int((k == 0).sum())
            hists[key] += np.bincount(k, minlength=ks_max)[:ks_max]
            for c in count_grid:
                if c == 1:
                    continue
                kc = (base >= c).sum(axis=1)
                hists[(key, c)] += np.bincount(kc, minlength=ks_max)[:ks_max]
        rows_seen += vals.shape[0]
        chunk_tracker += 1
        if chunk_tracker % print_freq == 0:
            print(f"{rows_seen} rows processed.")

    if rows_seen == 0:
        raise ValueError(f"Input file {input} is empty")

    for d in depth_grid:
        key = fmt_num(d)
        out_path = spectrum_out + "_incidence_d" + key + ".txt"
        with open(out_path, "w") as dst:
            for k in range(1, ks_max):
                dst.write(f"{k},{hists[key][k]}\n")
        print(f"Wrote incidence spectrum to {out_path}")
        for c in count_grid:
            if c == 1:
                continue
            out_path = spectrum_out + "_incidence_d" + key + "_c" + str(c) + ".txt"
            with open(out_path, "w") as dst:
                for k in range(1, ks_max):
                    dst.write(f"{k},{hists[(key, c)][k]}\n")
            print(f"Wrote c-thresholded incidence spectrum to {out_path}")

    k_total = int(hists[fmt_num(0.0)][1:].sum())
    if k_total == 0:
        raise ValueError("No k-mer is present in any sample; check the input table.")

    # rarefaction (empirical) + neutral Ewens comparison from the unthinned
    # presence spectrum
    h0 = hists[fmt_num(0.0)][0]
    U_points = rarefy_disjointprob(N, M_GRID, hists[fmt_num(0.0)])
    U1 = U_points[0][1]
    if U1 <= 0:
        raise ValueError("U(1) is zero; all sample columns are empty?")
    u_emp = [(m, U / U1) for m, U in U_points]
    z_emp, r2_emp = fit_power_law([m for m, _ in u_emp if m >= 2], [u for m, u in u_emp if m >= 2])

    theta = ewens_theta_from_total(N, k_total)
    theo_h = np.zeros(ks_max, dtype=np.float64)
    for i in range(1, N + 1):
        theo_h[i] = theta / i
    theo_points = rarefy_disjointprob(N, M_GRID, theo_h)
    theo1 = theo_points[0][1]
    u_theo = [(m, U / theo1) for m, U in theo_points]
    z_theo, r2_theo = fit_power_law([m for m, _ in u_theo if m >= 2], [u for m, u in u_theo if m >= 2])

    rare_path = spectrum_out + "_rarefaction.txt"
    with open(rare_path, "w") as dst:
        for (m, ue), (_, ut) in zip(u_emp, u_theo):
            dst.write(f"{m},{ue},{ut}\n")
    print(f"Wrote rarefaction and Ewens comparison curves to {rare_path}")

    info_path = spectrum_out + "_spectrum_info.txt"
    with open(info_path, "w") as dst:
        dst.write(f"N {N}\n")
        dst.write(f"table_rows {rows_seen}\n")
        dst.write("union_rows " + str(union_rows) + "\n")
        dst.write(f"genome_size {genome_size}\n")
        dst.write("k_total " + str(k_total) + "\n")
        dst.write("h0_frac " + str(h0 / float(rows_seen)) + "\n")
        dst.write("depth_grid " + ",".join(fmt_num(d) for d in depth_grid) + "\n")
        dst.write("count_grid " + ",".join(str(c) for c in count_grid) + "\n")
        dst.write("seed " + str(seed) + "\n")
        dst.write("theta " + str(theta) + "\n")
        dst.write("z_empirical " + str(z_emp) + "\n")
        dst.write("z_empirical_r2 " + str(r2_emp) + "\n")
        dst.write("u_empirical_last " + str(u_emp[-1][1]) + "\n")
        dst.write("u_empirical_last_m " + str(u_emp[-1][0]) + "\n")
        dst.write("z_theoretical " + str(z_theo) + "\n")
        dst.write("z_theoretical_r2 " + str(r2_theo) + "\n")
        dst.write("chunk_size " + str(chunk_size) + "\n")
        dst.write("per_sample_depth " + ",".join(repr(float(v)) for v in depth) + "\n")
    print(f"Wrote spectrum summary to {info_path}")
    print("Done!")

@click.command(context_settings={"show_default": True})
@click.option("-i", "--input", required=True, help="Path to k-mer count matrix")
@click.option("-p", "--prefix", required=True, help="Path to prefix file")
@click.option("-t", "--threads", default=1, help="Number of threads")
@click.option("-c", "--chunk-size", default=100000, help="Number of rows per chunk (lower means less memory is required)")
@click.option("-v", "--print-freq", default=100, help="Print a message every v chunks to indicate progress.")
@click.option("--ignore-first/--no-ignore-first", default=True, help="Toggle whether to ignore first column of text file (e.g. k-mer sequence column). Has no effect when --subset-rows is set, since subset mode always hashes and excludes the first column.")
@click.option("--subset-columns", type=click.IntRange(min=2), default=None, help="Randomly subset to this many columns before calculating distances. Must be an integer.")
@click.option("--subset-rows", type=click.IntRange(min=2), default=None, help="Randomly subset to this many rows by bottom-k min-hash of the first column (k-mer string) before calculating distances. The k-mer column is always excluded from the distance matrix in this mode. Requires --no-ignore-first. Writes the kept rows to {prefix}_subset_rows.txt.")
@click.option("--seed", default=1, help="Seed for the min-hash function used by --subset-rows and for the binomial thinning draws in spectrum mode.")
@click.option("--cumsum-chunk", type=click.IntRange(min=1), default=None, help="Opt-in: track the whole cumulative series every C rows of k-mer counts through one snapshot per chunk, not just the final total. Writes {prefix}_cumsum.txt (per-column sums, one row per snapshot) and {prefix}_cumsum_total.txt (two columns: cumulative unique k-mers = rows processed, and cumulative total k-mer count summed across all processed samples). Unset, the script behaves as before and writes neither file.")
@click.option("--normalize/--no-normalize", default=True, help="Divide cumulative sums so each snapshot is a fraction of full scale: per-column series by each column's total k-mer count from pass 1, and the two-column series by total rows (K) and the grand total of all processed counts (final row is 1.0 in every column; zero-total columns become 0.0). Applies only to the cumulative sums outputs, not to the distance metrics.")
@click.option("--spectrum-out", type=click.Path(), default=None, help="Enable spectrum mode: compute incidence spectra h(k) (k-mers present in exactly k samples) at each depth standardization point and count threshold, plus rarefaction and Ewens comparison curves; PREFIX for the outputs ({prefix}_incidence_d{d}.txt, {prefix}_rarefaction.txt, {prefix}_spectrum_info.txt). Requires --depth-totals. Distance outputs are not produced in this mode. Works on both tables with a k-mer string first column (excluded) and counts-only subset tables (auto-detected against the --depth-totals length).")
@click.option("--depth-totals", type=click.Path(exists=True), default=None, help="Per-sample summed k-mer counts as a comma-separated line, from a distance-mode run (--no-ignore-first) over the same table ({prefix}_sums.txt). Used to estimate per-sample depth for depth standardization.")
@click.option("--genome-size", type=float, default=None, help="Genome size G (assembly bp); depth standardization d is expressed in units of G.")
@click.option("--depth-grid", default=None, help="Comma-separated depth standardization points in G-units (default: 0,0.5,1,2,5,10).")
@click.option("--count-grid", default=None, help="Comma-separated count thresholds c for the c-thresholded incidence spectra (default: 1,2,5,10,30,100; c=1 is the plain presence spectrum).")
@click.option("--subset-rows-count", type=int, default=None, help="Rows in the input table (skips the row-counting pass).")
@click.option("--union-rows", type=int, default=None, help="Rows in the full union table the subset was drawn from; scales per-sample depth estimates. Defaults to the subset row count.")
@click.option("-s", "--seperator", default=" ", help="Separator between columns")


def main(input, prefix, threads, chunk_size, print_freq, ignore_first, subset_columns, subset_rows, seed, cumsum_chunk, normalize, spectrum_out, depth_totals, genome_size, depth_grid, count_grid, subset_rows_count, union_rows, seperator):
    """Average Bray-Curtis and cosine distances over all column pairs, or k-mer incidence spectra if --spectrum-out is set."""
    # Interpret escape sequences like '\t' (single quotes in bash pass the
    # literal two characters backslash+t) so both --seperator '\t' and
    # --seperator $'\t' work.
    seperator = seperator.encode().decode("unicode_escape")

    if spectrum_out is not None:
        run_spectrum(
            input=input,
            spectrum_out=spectrum_out,
            depth_totals=depth_totals,
            genome_size=genome_size,
            depth_grid=parse_grid(depth_grid, float) if depth_grid is not None else [0.0, 0.5, 1.0, 2.0, 5.0, 10.0],
            count_grid=parse_grid(count_grid, int) if count_grid is not None else [1, 2, 5, 10, 30, 100],
            seed=seed,
            subset_rows_count=subset_rows_count,
            union_rows=union_rows,
            chunk_size=chunk_size,
            print_freq=print_freq,
        )
        return

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

    if subset_rows is not None:
        # Subset-rows mode always hashes and excludes column 0 (the k-mer
        # string) from the distance matrix.
        print("Subset-rows mode: excluding first column (k-mer strings) from distance matrix...")
        cols_kept=range(1,ncols)
    elif ignore_first:
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
    elif ignore_first:
        raise click.ClickException(
            "If --subset-rows is set, then --no-ignore-first should also be set, "
            "because the first column needs to be a k-mer string to min-hash."
        )
    else:
        print(f"Pass 0: bottom-k min-hash subsetting to {subset_rows} rows (seed={seed})...")

        # Single pass over the raw lines: hash the k-mer (first field only,
        # split() avoids hashing the rest of the line), keep the n smallest
        # hashes in a bounded max-heap of packed keys.
        heap = []
        index_value = 0
        with open_text(input) as x:
            for line in x:
                kmer = line.split(seperator, 1)[0]
                hash_value = mmh3.hash(kmer, seed, signed=False)
                add_heap(heap, hash_value, index_value, subset_rows)
                index_value += 1

        total_rows = index_value
        rows_kept = sorted(-k & 0xFFFFFFFF for k in heap)
        print(f"Hashed {total_rows} rows; kept {len(rows_kept)} rows.")
        if len(rows_kept) < subset_rows:
            raise ValueError(
                f"Input file {input} has {total_rows} rows, fewer than "
                f"--subset-rows {subset_rows}."
            )

        # Second linear scan: write the kept rows, in original file order,
        # byte-for-byte, to {prefix}_subset_rows.txt. This file is both the
        # provenance record and the input for passes 1 and 2.
        subset_path = prefix + "_subset_rows.txt"
        print(f"Writing min-hash subset rows to {subset_path}...")
        with open_text(input) as src, open(subset_path, "w") as dst:
            keep_iter = iter(rows_kept)
            next_keep = next(keep_iter, None)
            for index_value, line in enumerate(src):
                if index_value == next_keep:
                    dst.write(line)
                    next_keep = next(keep_iter, None)
                    if next_keep is None:
                        break

        print("Loading k-mer count matrix from subset file by chunk...")
        # Subset file has the same columns as the input; column 0 was already
        # dropped from cols_kept above, --subset-columns already applied.
        reader = pd.read_csv(subset_path, sep=seperator, header=None, chunksize=chunk_size, dtype=np.float64, usecols=cols_kept)
    #print("Loading k-mer count matrix by chunk...")
    #reader = pd.read_csv(input, sep=seperator, header=None, chunksize=chunk_size, dtype=np.float64, usecols=cols_kept)

    print("Pass 1: column sums...")
    col_sums = None
    col_non_zero = None
    chunk_tracker = 0

    cumsum_path = prefix + "_cumsum.txt"
    totals_path = prefix + "_cumsum_total.txt"
    cumsum_file = None
    totals_file = None
    cum_runs = None
    total_runs = 0.0
    rows_since_snapshot = 0
    rows_seen = 0
    if cumsum_chunk is not None:
        print(f"Writing cumulative per-column k-mer sums every {cumsum_chunk} rows to {cumsum_path}...")
        print(f"Writing cumulative k-mers-seen and total-k-mer-count series every {cumsum_chunk} rows to {totals_path}...")
        cumsum_file = open(cumsum_path, "w")
        totals_file = open(totals_path, "w")

    for chunk in reader:
        s = chunk.sum(axis=0).values
        b = chunk.astype(bool).sum(axis=0) 
        col_sums = s if col_sums is None else col_sums + s
        col_non_zero = b if col_non_zero is None else col_non_zero + b

        if cumsum_file is not None:
            # Both series snapshot at the same C-row boundaries, which can
            # fall inside a read chunk when C != chunk_size: split the chunk
            # at the nearest C multiple each time. cum_runs (per-column) and
            # total_runs (scalar, summed across samples) keep accumulating
            # across snapshots; only the boundary counter resets.
            vals = chunk.values
            pos = 0
            while pos < vals.shape[0]:
                take = min(cumsum_chunk - rows_since_snapshot, vals.shape[0] - pos)
                seg_sum = vals[pos:pos + take].sum(axis=0)
                cum_runs = seg_sum if cum_runs is None else cum_runs + seg_sum
                total_runs += seg_sum.sum()
                rows_since_snapshot += take
                pos += take
                if rows_since_snapshot == cumsum_chunk:
                    rows_processed = str(rows_seen + pos)
                    cumsum_file.write(rows_processed + "," + ",".join(map(str, cum_runs)) + "\n")
                    totals_file.write(rows_processed + "," + str(total_runs) + "\n")
                    rows_since_snapshot = 0
            rows_seen += vals.shape[0]

        chunk_tracker += 1
        if chunk_tracker % print_freq == 0:
            print(f"{chunk_tracker*chunk_size} rows processed.")

    if col_sums is None:
        raise ValueError(f"Input file {input} is empty")

    if cumsum_file is not None:
        # Counts are integers, so the C-segment sums and the whole-chunk
        # col_sums sums are both exact float64 values: each final total equals
        # its pass-1 counterpart bitwise, making the normalized final rows
        # exactly 1.0 and the series deterministic across any --chunk-size
        # grouping.
        if rows_since_snapshot > 0:
            cumsum_file.write(str(rows_seen) + "," + ",".join(map(str, cum_runs)) + "\n")
            totals_file.write(str(rows_seen) + "," + str(total_runs) + "\n")
        cumsum_file.close()
        totals_file.close()
        if normalize:
            print(f"Normalizing cumulative sums in {cumsum_path} by per-column totals...")
            tmp_path = cumsum_path + ".tmp"
            with open(cumsum_path) as src, open(tmp_path, "w") as dst:
                for line in src:
                    head, tail = line.rstrip("\n").split(",", 1)
                    snap = np.array(tail.split(","), dtype=np.float64)
                    with np.errstate(divide="ignore", invalid="ignore"):
                        snap = snap / col_sums
                    snap[col_sums == 0] = 0.0
                    dst.write(head + "," + ",".join(map(str, snap)) + "\n")
            os.replace(tmp_path, cumsum_path)
            print(f"Normalizing cumulative series in {totals_path} by total rows and grand total...")
            k_total = rows_seen
            grand_total = col_sums.sum()
            tmp_path = totals_path + ".tmp"
            with open(totals_path) as src, open(tmp_path, "w") as dst:
                for line in src:
                    head, tail = line.rstrip("\n").split(",", 1)
                    k_frac = float(head) / k_total
                    count_frac = 0.0 if grand_total == 0 else float(tail) / grand_total
                    dst.write(str(k_frac) + "," + str(count_frac) + "\n")
            os.replace(tmp_path, totals_path)
    
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
    if subset_rows is None:
        reader = pd.read_csv(input, sep=seperator, header=None, chunksize=chunk_size, dtype=np.float64, usecols=cols_kept)
    else:
        # Read the same subset rows as pass 1 so both passes are consistent.
        reader = pd.read_csv(subset_path, sep=seperator, header=None, chunksize=chunk_size, dtype=np.float64, usecols=cols_kept)
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
