# K-mer rarefaction across individuals: math reference

This note documents the mathematics behind the k-mer incidence-spectrum
outputs ("tier 2"), produced by `kmer_distances.py --spectrum-out` (see
`workflow/rules/kmc_nofilt.smk`, rule `kmer_spectrum_nofilt`) and analyzed in
`workflow/notebooks/red_list_genomics.Rmd` ("k-mer incidence spectra" section).

The goal: a k-mer accumulation curve that can be compared **across species**,
where the variation reflects within-population genomic diversity rather than
genome size or sequencing depth.

The math below is written for the **full union table**: the complete
count-filtered k-mer matrix for one species (all distinct k-mers of the
species' samples), not a sampled subset of its rows. All quantities are then
exact properties of that table — nothing in sections 1-5 involves sampling
from the union. (Section 8's pilot was run on minhash subsets of the union;
the same formulas apply unchanged, with a sampling error noted there.)

---

## 1. Setup and notation

Input: the full union count matrix for one species, rows = k-mers (every row
of the species' paste table), columns = samples (individuals), cells = counts
of that k-mer's reads in that sample. Let

| symbol | meaning |
|---|---|
| $N$ | number of samples (individuals) |
| $K$ | number of rows of the table = number of union k-mers above the count filter |
| $G$ | genome size in bp (assembly length, from the reference `.fai`) |
| $S_j$ | sample $j$'s summed k-mer count over the table (column sum) |
| $c_{ij}$ | count of k-mer $i$ in sample $j$ (an element of the matrix) |
| $h(k)$ | incidence spectrum: number of k-mers present (count $\ge 1$) in exactly $k$ of the $N$ samples |
| $h_c(k)$ | count-thresholded spectrum: k-mers with count $\ge c$ in exactly $k$ samples |

By construction of the union dump, every row has at least one sample with a
nonzero count, so $h(0) = 0$ and

$$
K \;=\; \sum_{k=1}^{N} h(k).
$$

Because the table holds all union k-mers exactly once each per copy... (rows
are distinct k-mer sequences), no rescaling of table statistics is needed
anywhere: sums, spectra, and curves below are computed directly.

---

## 2. Depth estimation and thinning (the sampling model)

### 2.1 Per-sample depth

Sample $j$'s depth in k-mer reads is exactly its column sum of the table:

$$
D_j \;=\; S_j.
$$

Two subtleties worth stating because the downstream target assumes them:

- $D_j$ counts **all** k-mer copies, so it includes repetitive and organellar
  content: it is "reads per genome-equivalent," not per-base nuclear coverage.
- $D_j$ is defined over the union table (post count filter). If a per-sample
  k-mer count file from an earlier stage of the pipeline is used instead
  (pre-filter), the thinning probabilities shift by the filter's truncation —
  the depth totals and the table must come from the **same filtering step**.

### 2.2 Thinning probability

To evaluate every species at the same standardized depth, a target expressed
in **genome units** is used:

$$
\mathrm{target}(d) = d \cdot G, \qquad d \in \{0,\ 0.5,\ 1,\ 2,\ 5,\ 10\}
$$

The per-sample thinning probability is

$$
p_j(d) \;=\; \mathrm{clip}\!\left(\frac{d \cdot G}{D_j},\ 0,\ 1\right)
\qquad\text{with } \mathrm{clip}(x, 0, 1) = \min\{1,\max\{0, x\}\}.
$$

Samples already shallower than the target get $p_j = 1$ (the *clamp*: data
cannot be thinned up, so those points reflect native depth). Because $G$ is
an assembly length (§4.3), it is a proxy; the $d$-axis is comparable across
species only as "coverage-equivalent reads kept," not as absolute truth.

### 2.3 The binomial draw

Depth standardization never re-counts reads. Downsample a sample to a
fraction $p$ of its reads and each observed count becomes a binomial
resample:

$$
c'_{ij} \sim \mathrm{Binomial}(c_{ij},\ p_j(d)).
$$

This is statistically equivalent to read-level subsampling for every marginal
statistic the analysis uses (per-k-mer counts, presence probabilities), the
same equivalence underlying rarefaction in ecology and `samtools view -s`.
The equivalence assumes reads are independent draws (i.i.d.); PCR duplicates
make counts clump, so real thinning is slightly *more* informative than
read-subsampling with the same seed — a second-order caveat recorded in §7.
Draws are seeded (`--seed`, default 1) for reproducibility; the $d = 0$ point
makes **no draws at all** (raw presence) and its file is what the curves in
section 3 are computed from.

Expectation and presence after thinning:

$$
\mathbb{E}[c'_{ij}] = c_{ij}\, p_j(d), \qquad
\Pr(c'_{ij} \ge 1) = 1 - (1-p_j(d))^{c_{ij}}.
$$

Note the semantics: thinning equalizes the expected **count** ($d \cdot G$
per sample), not the detection probability of rare k-mers; a single-copy
k-mer observed once ($c = 1$) in a $5{\times}G$-deep sample survives a
$d=1$ thinning with probability $1/5$, whereas hypothetical re-sequencing
at $1\times$ would detect it with probability $\approx 0.63$. The $d$-axis
means "reads kept" (rarefaction semantics), not re-sequencing coverage.

**Where depth actually matters for presence.** At standard WGS coverage
(counts per single-copy k-mer $\gtrsim 10$), raw presence is essentially the
true presence spectrum: the miss probability for a single-copy k-mer at
$10\times$ is $e^{-10} \approx 5\times10^{-5}$, at $40\times$ ~$10^{-18}$.
So the $d = 0$ spectrum is the depth-truth-leaning point, and thinning to
$d < 1$ mainly serves (a) count-resolution standardization for the
$c$-thresholded spectra, (b) equalization for genuinely shallow samples
(the clamp's insurance role), and (c) robustness probing of the
accumulation-curve summaries across detection regimes — not correction of a
large presence bias. If some samples sit near $\approx 1\times$ effective
nuclear coverage (organellar-heavy or otherwise poor libraries), their
$d=0$ presence is noisy; the remedy is the per-species report of
$d_{\max} = \min_j (D_j / G)$ and comparison at $d \le d_{\max}$.

### 2.4 Edge cases

- Zero-depth columns: $p_j = 1$ (the clip's ceiling), and every cell is 0,
  so the sample contributes nothing regardless.
- $p_j = 1$ (clamp or $d$ large): the thinned equals the unthinned spectrum
  bit-for-bit (asserted in tests).
- Reproducibility: same seed → byte-identical thinned spectra; the $d=0$
  spectrum is seed-independent.

---

## 3. Rarefaction: the accumulation curve

### 3.1 Definition

For a subset of $m$ of the $N$ individuals, the expected number of **distinct**
k-mers present in at least one of them is the sample-based rarefaction
estimator (Heck et al. 1975; incidence-interpolation curve of iNEXT at hill
number $q=0$):

$$
U(m) \;=\; \sum_{k=1}^{N} h(k)\;\Bigl(1 - \frac{\binom{N-m}{k}}{\binom{N}{k}}\Bigr).
$$

### 3.2 Where the formula comes from

Fix one k-mer $i$ and let $k = $ the number of samples that observed it. A
random $m$-subset $M$ *misses* $i$ exactly when all $m$ picks fall in the
$N-k$ samples that lack it:

$$
\Pr(S \cap M = \emptyset) = \frac{\binom{N-m}{k}}{\binom{N}{k}}.
$$

(two symmetric countings of the same product $N!/(m!\,k!\,(N-m-k)!)$: choose
$S$ first then $M$, or $M$ first then $S$). The k-mer contributes 1 to the
union with the complementary probability. Summing over all K k-mers is exact
by **linearity of expectation** — the joint dependence of overlapping
k-mers, LD, and population structure never enters; only the one-dimensional
histogram $h(k)$ is needed. That is the compression that makes the
K-dimensional multinomial tractable: $U(m)$ is *the exact expectation*,
no per-k-mer parameters appear.

Sanity checks (all hold on the full table by the identities in §1 and §3.1):

- $U(1) = \sum_k h(k)\, k / N$ = the average per-sample k-mer load $S$.
- $U(N) = \sum_{k\ge 1} h(k) = K$ = the full union (the m-grid always
  includes this endpoint).
- k-mers present in every sample ($k = N$) contribute the constant
  $h(N)$ — the flat baseline; all *growth* of the curve is diversity.
- Monotonically increasing in $m$.

### 3.3 Per-step novelty (the original motivating quantity)

The expected number of so-far unseen k-mers contributed by individual
$m{+}1$ is the discrete derivative; summing over the histogram,

$$
U(m{+}1) - U(m) \;=\; \sum_k h(k)\;\cdot\;\Pr(\text{disjoint after } m)\;\cdot\;\frac{k}{N-m},
$$

— survive $m$ draws unseen (weight $\binom{N-m}{k}/\binom{N}{k}$), then
appear with probability $k/(N{-}m)$. The per-sample-conditional form is the
"probability of sampling something new given everything sampled so far"
quantity; the curve is its cumulative version.

### 3.4 Scaling behavior (why a power law)

When $h(k)$ is smoothly declining (the usual shape for k-mer data), $U(m)$
is well-approximated by a power law over the fitted range:

$$
u(m) \approx c\, m^{\,z}.
$$

The fitted $z$ (least squares on $(\log m, \log u)$ over $m \ge 2$ including
the $m = N$ endpoint, reported with $R^2$) is the one-number-per-species
summary, directly analogous to the species–area / MAR exponent. There is no
universal $z$: the exponent depends on the shape of $h(k)$ — how much of the
union is shared by most individuals (flat baseline, pulls $z$ down) versus
private to few (pulls $z$ up) — and on where the fitted range ends, which is
why the fitted *range* and the endpoint are reported rather than assumed.
The neutral comparison curve (§5) has its own reference exponent,
$z_{\text{theo}}$, so departures can be read species by species.

Two cautions carried over from the pilot (§8): near-perfect $R^2$ is
expected for ANY smooth saturating curve on these grids, so it does not by
itself validate the power law against competitors (exponential saturation,
two-pool mixtures); and because different species have different $N$, the
fitted windows differ — report both the common-window fit ($m \le \min N$)
and the endpoint-inclusive fit, plus the model-free $u(2), u(5), u(10)$
points, before reading $z$ as mechanism.

---

## 4. Genome-size anchoring

### 4.1 The anchor: $u(1) \equiv 1$

$U(1)$ is the per-individual k-mer repertoire: every single-copy genomic
k-mer of the reference haploid genome $G$, plus distinct repeat sequences,
plus anything above the count filter. It is dominated by $G$ — compare two
species with $G$ = 700 Mb and 8 Gb: raw $U(m)$ separates by genome size
almost as much as by anything else. The normalized curve

$$
u(m) = \frac{U(m)}{U(1)}
$$

starts at exactly 1 for *every* species; the unit of the y-axis becomes
"one individual's k-mer repertoire". The variation along the curve,
$U(m) - U(1)$, counts only k-mers seen in some individuals but not in a
given first one — within-population genomic variation — and to first order
in $G$ the ratios are scale-free: a larger genome inflates $U(1)$ and the
single-copy share of every $U(m)$ roughly proportionally.

The normalization divides by the **measured** $U(1)$, not by the assembly
size: assembly incompleteness, organelle contigs, and haploid-vs-diploid
assembly artifacts would otherwise enter the denominator. $G$ itself is used
only where a *target* is needed (§2.2). The ratio $U(1)/G$ ("k-mer repertoire
per assembly bp") is reported per species as the diagnostic for when this
proportionality fails: values far above 1 indicate repeat-expanded or
organelle-inflated repertoires, i.e. exactly the conditions under which the
anchoring argument weakens.

### 4.2 The second $G$: depth targets in genome units

The depth grid of §2.2 expresses effort as $d$ genome-equivalents of reads.
Comparing species at equal read *totals* would over-sample a small genome
(and saturate it) and under-sample a large one; equal *coverage-equivalent*
effort puts a single-copy k-mer at expected count $\approx d$ in every
species. In expectation the two normalizations divide out "how big is one
genome" and "how much did we sequence", and what survives is what the
cross-species comparison measures.

### 4.3 What the anchoring does not fix

- $U(1) \ne G$: repeat-rich genomes have $U(1)/G > 1$; above-threshold
  residual errors also contribute. The anchor cancels $G$ only to first
  order — and the leading violation is repeat architecture, which varies
  *between* species, not just within. Watch the $U(1)/G$ column.
- Copy-number novelty (structural variation, recent WGD, ploidy §7) raises
  $U(m)$ alongside sequence novelty — part of the signal, but not
  allele-frequency diversity.
- Clamped depths (samples below target $d$) make those $d$-points reflect
  native depth; the standardization is only as good as the depth estimates
  (which inherit any column-sum bias or pre-filter mismatch, §2.1).
- $G$ from the reference `.fai` is the assembly bp, which *understates* the
  true genome size systematically for repeat-rich species; the k-mer-derived
  alternative (a GenomeScope2-style fit on this species' `results/kmc_histo/
  {species}.txt`) estimates $G$ and repeat content from the data and is the
  preferred target denominator when available.

---

## 5. The comparison curve: neutral theory (Ewens) baseline

To interpret an empirical $z$, the script also computes the curve a neutral
model would produce. Under the infinite-sites Ewens model with parameter
$\theta$, the expected number of segregating sites in a sample of $N$ is
$\theta H_N$ ($H_N = \sum_{i \le N} 1/i$), matching $\theta$ to the table's
union size:

$$
\hat\theta = \frac{K}{H_N}, \qquad
\hat h(k) \;\propto\; \frac{\theta}{k}, \qquad k = 1 \dots N.
$$

Running the same rarefaction (§3.1) on $\hat h$ yields the *theoretical* curve

$$
u_{\text{theo}}(m) = \frac{U_{\text{theo}}(m)}{U_{\text{theo}}(1)}, \qquad
U_{\text{theo}}(m) = \sum_{k=1}^{N} \frac{\theta}{k}
\Bigl(1 - \frac{\binom{N-m}{k}}{\binom{N}{k}}\Bigr),
$$

with two exact endpoints worth knowing: $u_{\text{theo}}(1) = 1$ by
construction, and

$$
u_{\text{theo}}(N) = H_N
$$

independent of $\theta$ — because the neutral union size is exactly
$\sum_k \theta/k = \theta H_N$ while the single-sample load is exactly
$\theta$ ($\sum_k (\theta/k)\cdot k/N = \theta$). The theoretical curve is
therefore a pure function of $N$ (it saturates like $H_m \sim \log m$), and
over a fixed fitted range its exponent $z_{\text{theo}}$ is a known
benchmark.

Both curves are written side by side (`{prefix}_rarefaction.txt`:
`m, u_emp, u_theo`). This is the same construction as the `mar` package's
theoretical mutations-area curve (`MARtheory` on an incidence spectrum),
computed independently of that package's pipeline.

**Diagnostic use:** the gap $z_{\text{emp}} - z_{\text{theo}}$ (both reported
with $R^2$ in the info file) measures how far a species' k-mer incidence
structure deviates from neutral expectations — selection, population
structure, and error structure all move it. A large gap is itself a
comparable per-species parameter, not just a nuisance.

**Interpretation limits for $\hat\theta$.** One real SNP contributes ~31
overlapping k-mer rows with nearly identical incidence patterns, so $K$ —
and hence $\hat\theta$ — is an SFS *overcounted by the variant-cluster
factor* (Khaya: $\hat\theta = 4.1$M "segregating sites" against ~130k
plausible SNPs, factor ~31). The global rescaling cancels in
$u_{\text{theo}}$ (and in $z_{\text{theo}}$) but not in $\hat\theta$ itself;
the factor is also mildly k-dependent (different variant types yield
different cluster sizes). Report $\hat\theta$ only alongside the
expected-redundancy correction $\hat\theta/\langle \text{cluster size} \rangle$
(cross-checkable against the jules track's SNP counts where available).

---

## 6. Outputs produced per species

All under `results/nofilt/spectrum/` (rule `kmer_spectrum_nofilt`):

| file | content |
|---|---|
| `{sp}_incidence_d{d}.txt` | $h(k)$ after depth-$d$ thinning; rows `k,h(k)`, $k \ge 1$ |
| `{sp}_incidence_d{d}_c{c}.txt` | $h_c(k)$ (count $\ge c$ after thinning) |
| `{sp}_rarefaction.txt` | rows `m, u_emp(m), u_theo(m)` |
| `{sp}_spectrum_info.txt` | key-value summary: $N$, rows ($K$), $G$, `h0_frac`, $\hat\theta$, $z_{\text{emp}}/z_{\text{theo}}$ with $R^2$, per-sample depths, grids, seed |

(QA note: on the full union table $h(0) = 0$ structurally, so `h0_frac`
should be exactly 0; if it is not, the inputs did not come from the same
paste/dump generation — re-run the paste before trusting the spectra. The
same reading applies verbatim to subset tables: nonzero `h0_frac` there
flags a sample-set mismatch between table versions, a condition the count
filter by itself cannot produce.)

---

## 7. Assumptions, caveats, and validity

1. **Exchangeability of samples.** The rarefaction design is "random subset
   of samples"; if the species' samples are structured (geography, subpopulations)
   the curve measures structure + diversity jointly — intended here, but say so.
   Note spatial/sampling design (e.g. `MARsampling`) is *not* used: subsets are
   uniform.
2. **Presence and depth.** At standard WGS coverage raw presence ($d = 0$)
   is essentially true presence (§2.3); presence bias grows only as samples
   approach $\approx 1\times$ effective coverage. The $d$-grid's primary
   function is count-resolution standardization and shallow-sample insurance,
   not correction of presence bias.
3. **Count thresholds vs. errors.** Small-$c$ spectra are contaminated by
   sequencing errors surviving the count filter; fit and interpret from the
   mid-$c$ tail (the count grid $\{1,2,5,10,30,100\}$ exposes the region).
4. **Small $N$.** $z$ from $N \le 4$ is not meaningful (the m-grid has too
   few points); the pipeline still writes their spectra, and notebook
   analyses filter accordingly. Bootstrap over individuals for CIs in
   downstream use.
5. **Row dependence.** The rarefaction needs only marginal counts, so k-mer
   correlations (overlapping windows) are harmless for $U(m)$, $u(m)$, and
   $z$ — they are NOT harmless for $\hat\theta$ (§5, cluster-size
   overcount).
6. **Ploidy/copy number.** $h(k)$ mixes sequence novelty with copy-number
   variation; polyploid species carry structurally different per-locus count
   distributions at the same $d$ (each homeolog sampled at roughly
   $d$/ploidy effective coverage). Record and test ploidy before reading
   $z$ as diversity alone.
7. **Model competition.** The power law is the default summary, not an
   established mechanism; near-perfect $R^2$ on these grids is expected for
   competitors too (§3.4). Model-free curve values ($u(2), u(5), u(10)$,
   $u(N)$) are the robust cross-species vectors.

---

## 8. Pilot validation (2026-10-30, Savio)

- 241 minhash-subset tables (2M-10M rows of the union) processed; 167 with
  $N \ge 2$; 93 with a valid $z$ fit ($N \ge 5$). Sampling SE per histogram
  bin $\sim \sqrt{h(k)}/R$-scale, negligible for the fitted exponents.
- $z$ range $[0.021, 0.737]$, median 0.31; power-law fits median $R^2 =
  0.9991$; 83% of fits $R^2 > 0.99$. (Interpret with §3.4 and §7.7 in mind.)
- $N$-confound weak overall (corr$(z, \log N) \approx -0.03$ among all
  z-fitted species; $\approx +0.26$ restricted to the five $N \ge 200$
  species — mostly window artifact candidates, §3.4).
- Biologically coherent extremes: lowest $z$ — *Kokia* (×2), *Nicotiana*,
  *Kandelia*, *Ochetophila*, *Juglans*, *Lupinus*, *Pistacia*
  (restricted endemics / bottlenecked taxa); highest — *Erythroxylum*,
  *Pinus*, *Isoetes*, *Juncus*, *Mentha*, *Psidium*, *Typha*, *Brassica*
  (widespread / weedy / large-population species).
- End-to-end rule validation on `Khaya_senegalensis` (10M rows × 6 samples ×
  6 depths × 6 count thresholds) in 31 s; `u(2) = 1.651` vs 1.640 from the
  independent pilot sweep (difference explained by the extra sample and
  endpoint), $z$ 0.513 vs 0.510.

---

## References

- Heck, van Belle & Simberloff (1975) — exact sample-based rarefaction.
- Hurlbert (1971) — the "unbiased" form of expected species.
- Sanders (1968) — original rarefaction.
- Good (1953); Good & Toulmin (1956) — sample-coverage/new-species probability
  ($f_1/n$; the depth-axis analogue of §3.3).
- Ewens (1972) — the neutral SFS used in §5.
- Chao et al. / iNEXT — interpolation/extrapolation and coverage standardization
  (the framework the $u(m)$ curve sits in).
- Clauset, Shalizi & Newman (2009) — power-law model competition (the basis
  of the §3.4/§7.7 cautions).
- Serrano et al. (2014); Meixilin `mar` package — mutations-area relationship;
  the $z$-fit parallels.