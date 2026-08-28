# Next analysis — spatial XCI, 9w vs 78w Visium HD

Task list for Claude Code, working in `/Users/lachlang/Postdoc`.

Ordered so that the cheap, decisive tests come first. Tasks 1–3 determine whether
the headline result is real; do not invest in 5–9 until they pass.

---

## Status, 2026-08-28

Tasks **1–4**, **9** and **10** are implemented. Tasks **5–8** are deliberately
untouched, per the gate above.

Task 9 was done out of order because it protects 1–4 rather than following them:
this round adds a second SNP mask and a second set of PDFs per sample, so
unlabelled figures and un-fingerprinted beds would have made the new results
harder to attribute than the old ones.

Nothing has been run. Run order on the cluster:

```bash
# 0. Build the SNP and interval beds. Once, on the reference tree.
cd $BASE/GRCm39 && OUT_DIR=. Rscript ~/Postdoc/OCM_heart/core_escape_SNPs.R

# 1. Recount. 'force' is REQUIRED: the existing tables have neither the window
#    table nor the subset columns, and the job now refuses to reuse them.
sbatch ~/Postdoc/slurm/spatial_ase_sweep.slurm 64 force

# 2. Task 1 - where on chrX the excess sits. The decisive one.
sbatch ~/Postdoc/slurm/spatial_window_profile.slurm 64

# 3. Task 2 - how much of the CAST difference survives Xic masking.
sbatch ~/Postdoc/slurm/spatial_mask_comparison.slurm 64

# Tasks 3 and 4 come out of step 1's own log and PDF - no separate job.
# Only if step 3 says the Xic carries the effect, rebuild the tile map on the
# masked bed (~14 h/sample; the tile BAMs are reused, only scoring re-runs):
#   sbatch ~/Postdoc/slurm/spatial_sinto_tiles.slurm 64 0 no_Xic500kb
#   sbatch ~/Postdoc/slurm/spatial_tile_map.slurm 64 9w,78w no_Xic500kb
```

Design decisions worth knowing before reading the output:

- **One counting pass, many masks.** Rather than re-running the counter against a
  second SNP bed, `ase_bin_allele_counts.py` gained `--subset-bed NAME=PATH`
  (and `:complement`), which adds a `NAME_ref`/`NAME_alt` column pair per
  interval set within the single pass. Escape vs non-escape, Xic-masked vs not,
  and the imprinted controls are therefore all measured on **the same UMI
  calls**, which is the only way the partitions are comparable. Two separate
  runs would have differed by whatever else changed between them.
- **The imprinted control is split by imprinting direction.** H19/Igf2,
  Airn/Igf2r and Meg3/Dlk1 are *reciprocally* imprinted pairs. Pooled into one
  set they contribute UMIs carrying opposite parental alleles, and C(d) over that
  mixture sits near 0.5 on perfectly clean data — indistinguishable from a failed
  control. Hence `imppat` and `impmat` as separate sets, with the reciprocal
  overlaps subtracted from both (Airn 118.6→90.2 kb against Igf2r, Igf2r
  87.3→58.9 kb against Airn, Rtl1 and Meg3 likewise in the Dlk1–Dio3 cluster).
  The Kcnq1 domain is excluded: its maternal genes carried 40% of the maternal
  SNPs — Osbpl5 alone 29% — and their imprinting in mouse is largely
  placenta-specific, so in adult heart they would have been biallelic and would
  have dragged the control toward 0.5 on their own. Both sets are unavoidably
  concentrated (Snrpn is 82% of the paternal SNPs), so `impsnrpn` and `impmeg3`
  carry one locus each as a per-locus check on the aggregates — free in the same
  pass, where getting one afterwards would cost a second counting pass.
- **The imprinted control tests the per-UMI error floor, not spatial
  resolution.** Task 3 as written expects `C_imprinted(d)` to decay with
  distance; it should not. Imprinting is uniform across the tissue, so the
  expected curve is flat and high, and the pass condition is only the `C(4um)`
  near 1 part. No monoallelic-*and*-patchy locus exists to calibrate resolution
  against, so the chrX excess over the autosomal null remains the only handle on
  patch structure. Task 3's decision rule survives this; its stated mechanism
  does not.
- **The concentration test in task 1 is exact, not a fold-change threshold.**
  Ranking windows by excess and reporting the top few overstates concentration by
  construction: on a synthetic control of 200 pure-Poisson windows with no spike
  at all, the top 10 carried 13% of the "excess", 3.1x their share of normal
  output. The verdict is instead a per-window binomial test conditioned on each
  window's total across the two tile groups, which absorbs the global rate factor
  the Xi-wide hypothesis is about and cannot be fooled by a window merely being
  busy. Validated both ways on synthetic data: 0 significant windows on the null,
  1 (in the Xic, carrying 76% of the excess) on a planted spike.

## Context: what prompted this list

Working from `tile_ratio_map_64um.csv` and the two
`tile_size_sweep.pdf` outputs:

- The **tile allelic ratio** shows almost nothing between ages. Unweighted mean
  0.736 (9w) vs 0.720 (78w); in the shallowest autosomal-depth tercile 78w is
  *higher* (0.753 vs 0.737). The apparent "more tiles > 0.9 at 78w" is explained
  by depth (reweighting 9w to the 78w depth distribution moves its fraction > 0.9
  from 0.226 to 0.296 [0.281–0.311] against 78w's 0.352) and by 78w being the
  noisier library (autosomal null SD 0.048 vs 0.028).
- A **per-allele rate** — informative chrX UMIs of one allele per 1000 informative
  autosomal UMIs — looked much more promising: BL6 (Xa) 12.28 → 12.54
  (FC 1.021 [0.982, 1.060]), CAST (Xi) 4.79 → 6.94 (FC 1.449 [1.371, 1.531]),
  flat across six autosomal-depth bins.
- **But that result is carried by ~10% of tiles.** Dropping the top decile of
  tiles by X-share (`x_n/(x_n+a_n)`) collapses the CAST fold-change to
  1.092 [1.036, 1.154] and turns the BL6 fold-change into 0.897 [0.877, 0.918].
  In 78w those tiles have `x_n` = 198 against `a_n` = 3879 — 3.4x the chrX UMIs
  of other tiles at *lower* autosomal depth, which is not a plausible
  transcriptional state. They are full tiles (`n_bins` = 1024), their autosomal
  ratio is normal (0.514), and they are **not** spatially clustered
  (Moran's I on the indicator = 0.019), so this is not a fold or an edge artefact.
  9w has the same spikes (top-decile `x_n` 168 vs 77) but they are not
  CAST-biased there (0.722 vs 0.738); in 78w they are (0.561 vs 0.738).
- The **tile size sweep does not identify a patch scale.** `C(d)` is flat at 0.785
  from 4 um to 2 mm in both samples, the chrX−autosome excess is flat at 0.28 and
  identical between ages, and `half_decay()` returns NA. rho(s) declines
  monotonically from 8 um with no plateau. 64 um is defensible as the coverage
  knee, not as a patch-matched scale.
- `C(4um) = 0.785`, not ~1. At 4 um separation two UMIs are usually inside one
  cardiomyocyte and should nearly always share an allele. This implies ~12% of
  chrX UMI allele calls disagree with their own cell's XCI state, identically at
  4 um and 2 mm. Whether that 12% is escape/leak or per-UMI assignment error
  decides whether the +45% means anything.

Everything below is aimed at those four problems. Note throughout: **n = 1 animal
per age**, so nothing here can be reported as an age effect without replication.

---

## 1. Locate the excess chrX UMIs on the chromosome — DECISIVE

The single analysis that settles the headline result. The chrX aggregate is built
against whole-chromosome intervals with no gene requirement
(`slurm/spatial_ase_sweep.slurm:97-101`, `chr_annotation_mm39.bed`), so intronic
and intergenic chrX reads are included and gene identity is discarded before any
downstream script sees it.

- [x] Modify `spatial/ase_bin_allele_counts.py` to emit, in addition to its
      current per-bin output, a **per-bin genomic-window count table**: chrX
      partitioned into 100 kb windows, counts split by allele. Keep it optional
      behind a flag so the existing output path is unchanged.
- [x] Aggregate to 64 um tiles, stratify tiles by X-share decile, and produce a
      positional profile of chrX UMIs (allele-split) for the top decile vs the
      rest, per sample.
- [x] Acceptance criterion: we can state whether the excess CAST UMIs in
      high-X-share tiles are concentrated in one or a few windows — the Xic
      around chrX:100–104 Mb (mm39) being the prime suspect — or spread across
      the chromosome.
- [x] If concentrated: the +45% is a locus artefact, and tasks 5–8 are not worth
      running until task 2 removes it. If spread across many windows in
      proportion to their normal expression: it is a genuine Xi-wide effect and
      becomes the paper's central result.

## 2. Mask the Xic properly, not just the Xist gene span

`core_escape_SNPs.R:11-16` drops SNP positions inside the single annotated Xist
interval. Ftx, Jpx, Tsix outside that span and all intergenic Xic reads survive,
and in this design they are effectively CAST-only.

- [x] Build a second SNP bed masking Xist ± 500 kb. Name it explicitly
      (e.g. `SNPfile_..._mm39_no_Xic500kb.bed`) and record which analyses use it.
- [x] Fix the latent bug at `core_escape_SNPs.R:13`: `seq(Xist_annotation$V2,
      Xist_annotation$V3)` assumes exactly one Xist row in `annotation_us.bed`
      and errors otherwise. Use the `unlist(Map(seq, ...))` pattern already used
      eleven lines later at `core_escape_SNPs.R:31`.
- [x] Note the coordinate convention while you are there: `core_escape_SNPs.R:14`
      matches SNP column `V2` against a span built from BED `V2..V3`, but
      `slurm/spatial_ase_sweep.slurm:118-121` establishes that column 2 of this
      SNP bed is a 1-based VCF position. Off-by-one at interval edges only —
      document it or fix it, don't leave it ambiguous.
- [x] Re-run the per-allele rates and the tile ratio map with the new mask.
      Report both masks side by side; do not silently switch.
- [x] Acceptance criterion: we know how much of the CAST rate difference
      (4.79 vs 6.94) survives Xic masking.

## 3. Imprinted-locus positive control for spatial resolution

Nothing currently establishes that the 2 um binning resolves within-cell allele
sharing at all. The autosomal `C(d)` control sits at 0.50, which is exactly what
biallelic expression predicts whether or not the spatial signal is intact, so it
cannot detect a resolution failure. An imprinted locus is monoallelic in every
cell and gives a real positive control, in the canonical BL6xCAST mm39 cross.

- [x] Build a SNP bed restricted to well-characterised imprinted loci: H19/Igf2,
      Airn/Igf2r, Snrpn, Peg3, Meg3/Dlk1. Check per-locus informative UMI depth
      first and keep only loci with usable coverage.
- [x] Run `pair_correlation()` from `spatial/ase_tile_sweep.R:263` on that set.
- [x] Acceptance criterion: if `C_imprinted(4um)` is near 1 and decays with
      distance, the method resolves within-cell allele sharing and chrX's flat
      0.785 is biology. If `C_imprinted(4um)` is also ~0.78, spatial resolution
      or per-UMI allele assignment is the ceiling, the tile analysis cannot go
      further, and the project should pivot to the snRNA-seq.
- [x] Report the realised separations, not the nominal ones —
      `ase_tile_sweep.R:281-283` already notes the integer-rounding drift, which
      matters most at exactly the short range this test depends on.

## 4. Split C(d) into escape vs non-escape chrX

Cheap, and it partitions the 12% floor. The core-escape SNP bed already exists
(`core_escape_SNPs.R:23-25`, the 11-gene `new_core_escape_genes` set).

- [x] Run `pair_correlation()` three ways: all chrX, core-escape SNPs only,
      chrX minus core-escape SNPs.
- [x] Acceptance criterion: if removing escapees lifts short-range `C(d)` toward
      1, the floor is escape and the global skew is ~0.88. If `C(d)` stays at
      ~0.78, the floor is per-UMI error and the +45% is measuring noise rate.

---

Only proceed past this line once 1–4 have answers.

## 5. Make the per-allele rate the primary statistic

- [ ] New script `spatial/ase_allele_rates.R`. For each sample compute
      `1000 * sum(x_a1)/sum(a_n)` and `1000 * sum(x_a2)/sum(a_n)`, with
      tile-level bootstrap CIs (10k draws).
- [ ] Report BL6 and CAST rates separately in every table and figure. Never
      report the ratio alone — it conflates a numerator change with a
      denominator change, which is exactly what went wrong here.
- [ ] Built-in sensitivity panel, always printed: X-share deciles, autosomal
      depth sextiles, `x_n` thresholds (all / >=50 / >=100), full tiles only.
      A result that does not survive all four is not reported as a result.
- [ ] Leave-one-region-out: recompute dropping each spatial quadrant in turn, so
      a localised anomaly cannot carry the estimate.

## 6. Depth-match the two sections before comparing anything

The two sections do not have matched coverage at any tile size (9w median `x_n`
82, 78w 65; autosomal 5374 vs 4239), so no fixed tile size makes the two maps
comparable.

- [ ] Downsample the deeper section's UMIs to the shallower section's
      per-bin depth distribution, at the UMI level in
      `ase_bin_allele_counts.py` output, not by subsetting tiles.
- [ ] Re-derive tasks 5, 7 and 8 on the downsampled data.
- [ ] Acceptance criterion: the autosomal null SDs match between samples
      (currently 0.028 vs 0.048). Until they do, any dispersion comparison
      between ages is confounded.

## 7. Fix the tile map's null and calling

- [ ] `SINTO_MIN_UMI <- 10` (`ase_tile_sweep.R:105`) puts 10-UMI tiles in the
      published map. Raise it, and state the coverage cost: at 64 um only ~2%
      (9w) and ~0% (78w) of tiles reach 50 UMIs, so the honest options are a
      higher threshold with far fewer tiles, or a larger tile.
- [ ] Replace the null used by `spatial/tile_ratio_map.R` for the `z` column.
      Fitting `var(a_ratio) = phi * p(1-p)/n + tau^2` on the autosomal control
      gives phi = 10.6, tau = 0.0145 (9w) and phi = 14.4, tau = 0.0371 (78w),
      i.e. the autosomal ratio is 15–35x overdispersed relative to binomial.
      Propagated to chrX depth the null SD is ~0.18–0.22 at n = 65, so only
      1.2% (9w) / 3.6% (78w) of tiles differ significantly from their own sample
      mean — against ~62% currently called skewed.
- [ ] Use a beta-binomial null centred on the sample's own global ratio, not on
      0.5, and not the autosomal SD directly (the autosomal SD is measured at
      n ~ 5000 and is being applied at n ~ 65).
- [ ] Acceptance criterion: the recalled map should be almost entirely "mixed",
      and the figure legend should say so. If it isn't, the null is still wrong.

## 8. Re-derive the tile size choice honestly

- [ ] Rewrite the sweep's interpretation. The subtitle currently reads "chrX
      above the autosomal null is real structure; the fall-off marks where tiles
      start averaging over more than one patch". A *scale-invariant* excess is
      not patchiness — `C(d)` excess flat at 0.28 across three orders of
      magnitude is what a global skew of p ~ 0.88 with no spatial structure
      predicts (0.88^2 + 0.12^2 - 0.5 = 0.27).
- [ ] Check whether rho(s)'s slower-than-1/s^2 decline is explained by the
      spatial smoothness of *depth* rather than of ratio: Moran's I on `x_n` is
      0.19–0.21 and on `a_n` is 0.67–0.80, while on `x_ratio` it is 0.003 (9w)
      and 0.023 (78w). Regress rho(s) on tile-level depth variance to test this.
- [ ] Restate the 64 um justification as coverage-based, with the numbers:
      92% (9w) / 73% (78w) of tiles clear 10 UMIs at 64 um.

## 9. Reproducibility fixes

- [x] Put `SAMPLE` in the plot titles in `ase_tile_sweep.R:409-421`. Both
      `tile_size_sweep.pdf` files are currently byte-distinct but visually
      unidentifiable, distinguishable only by parent directory — the two files
      have already been confused once.
- [x] Record SNP-bed provenance in the output. Three "no-Xist" builds are named
      in the repo (`..._mm39_no_Xist.bed`, `..._withoutXist.clean.bed`, and the
      unfiltered original) and only one is live. Write the bed path and its
      md5 into a sidecar file next to every result.
- [x] Emit informative-SNP and informative-gene counts per sample. Nothing
      currently reports how many distinct chrX SNPs or genes contribute to the
      aggregate, which makes the 12% floor impossible to attribute.

## 10. Bugs found in the snRNA-seq pipeline (independent of the above)

- [x] `MIN_READS` is referenced at `03_all_genes.R:309,325,340` and defined
      nowhere; line 325 will throw at runtime. The actual filter is a hardcoded
      `10` at `03_all_genes.R:298`.
- [x] `04_core_escape.R:21-22` hardcodes barcode field indices 7/8 in a
      `strsplit(dirname(barcodes),'_')` parse, under the tree name
      `Allelome.PRO2_core_escape_new`. `00_functions.R:213-216` documents that
      this index is 4/5 under `Allelome.PRO2/` and 6/7 under
      `Allelome.PRO2_all_genes/` and is "silently wrong under any new tree
      name". A wrong index yields barcodes that fail to match Seurat rather
      than erroring. **04 is the LOX analysis** — verify before trusting it.
      Use `parse_allelome_cell()` instead.
- [x] `core_escape_genes_gene_df.txt` is read at `03_all_genes.R:211-215` but
      never written anywhere in the repo, so `is_core_escape` is silently
      all-FALSE and the escape overlay in
      `VCM_subcluster_per_gene_delta_scatter.pdf` is empty.
- [x] `allelic_ratio` means different things in different files:
      `10_build_ratio_table.R:63` sets it to `ar_dom` (bounded [0.5,1]),
      `04_core_escape.R:31,41` passes through the directional Allelome.PRO2
      column, and `05_lox_sensitivity.R:28-30` joins across the two conventions
      without documenting it. Rename one of them.
- [x] `01_setup.R:98-99` uses mm10 Xist coordinates in an otherwise mm39
      pipeline. Affects only the approximate-TPM plot.

---

## Do not do

- Do not report the tile allelic ratio as an age comparison. It is dominated by
  depth and by library noise, and the two sections are not depth-matched.
- Do not read a patch size off rho(s). It is monotone decreasing with no
  plateau; there is no scale to read.
- Do not treat n = 1 per age as a comparison of ages. Every statement needs
  either replication or explicit framing as a single-animal observation.
  `00_functions.R:41-53` already makes this point about the snRNA-seq.
