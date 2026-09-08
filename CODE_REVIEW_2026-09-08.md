# Code review: Postdoc repo, 2026-09-08

Reviewed by Claude (Fable 5.1). Written to be handed to another Claude session
(Opus 5) that will implement the fixes.

## Instructions for the implementing session

- Repo root is `/Users/lachlang/Postdoc`. Nothing in it can be run locally: every
  data path is on the LRZ cluster (`/dss/dssfs03/...`), and there are no tests.
  Verify changes by reading code and, where a check needs data, write the check
  as a command the user can run on the cluster.
- Read `spatial/NEXT_ANALYSIS.md` and `multiome/ANALYSIS_PLAN.md` first. They
  carry the biology that decides what the numbers mean: B6 mother x CAST father,
  Xist deleted on the B6 allele, so the CAST X is inactive in every cell and
  CAST fraction on chrX is escape.
- Section 5 lists things that look wrong but are deliberate. Do not "fix" them.
- Findings are ordered by priority within each section. Each has a file and
  line, what is wrong, why it matters, and a suggested fix.
- The most recent work (commits of 2026-09-04 to 09-07) is the cellranger-arc
  multiome pipeline, the `spatial_sinto_tiles` tile chain, and the
  `OCM_heart/allelic_ratio` series. The review concentrated there.

## 1. Correctness bugs

### 1.1 `spatial/tile_ratio_map.R:587` — `cut()` on quantile breaks will error at shallow tile sizes

`depth_confound()` does:

```r
q <- cut(s$x_n, breaks = quantile(s$x_n, probs = seq(0, 1, 0.25), na.rm = TRUE),
         include.lowest = TRUE, labels = FALSE)
```

`x_n` is a small integer at 32 um (the file's own comments say median 4 to 7).
When two quartiles coincide, `cut()` stops with "'breaks' are not unique" and
the whole `panel_call()` fails, taking the PDF with it. The same file already
solves this for `clustering_test()` with `depth_strata()` (line 623), which
rank-splits and cannot collide.

Fix: replace the quantile cut with `q <- depth_strata(s$x_n, 4L)`.

### 1.2 `slurm/spatial_sinto_tiles.slurm:573` — locus-table validator splits on whitespace, not tabs

`locus_table_state()` uses `awk` with the default field separator and treats any
row whose `NF` differs from the header as a truncation. With whitespace
splitting, a row with an empty tab-delimited field has one fewer field than the
header, so a complete table is classed as truncated, deleted (line 707), and the
tile is re-scored at about 10 min of CPU. This repeats on every chain link.

Fix: `awk -F'\t' -v min="$MIN_LOCI" '...'`. Also confirm on the cluster whether
`locus_table.txt` ever has an empty field:

```bash
f=$(ls $SCRATCH/spatial_tiles_9w_64um/allelome/*/locus_table.txt | head -1)
awk -F'\t' 'NR>1 { for (i=1;i<=NF;i++) if ($i=="") { print FILENAME": empty field row "NR" col "i; exit } }' "$f"
```

### 1.3 `spatial/tile_ratio_map.R:198` — tile ratio uses Allelome.PRO2's `total_reads` as the denominator

`read_locus()` sets `x_n = sum(x$total_reads)` and the ratio is `x_a1 / x_n`.
The OCM side deliberately does not trust that column: `00_functions.R:348`
recomputes `total_reads <- A1_reads + A2_reads`. `spatial/compare_ap2_pysam.R:293`
(section 2b) was written to test whether Allelome's `total_reads` equals
`a1 + a2`, but its verdict is not recorded anywhere in the repo. If the column
counts all reads over the locus rather than informative reads, every tile ratio
in the spatial figures has the wrong denominator.

Fix: run section 2b on the cluster (or check one table by hand), record the
answer in `spatial/NEXT_ANALYSIS.md`, and make `read_locus()` use
`a1_reads + a2_reads` regardless so the two sides of the project agree by
construction. `spatial/ase_tile_locus_counts.py:279` writes a `total_reads`
column too; state in its header comment what that column holds.

```bash
f=$(ls $SCRATCH/spatial_tiles_9w_64um/allelome/*/locus_table.txt | head -1)
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next} $h["A1_reads"]+$h["A2_reads"]!=$h["total_reads"]{print; bad++} END{print (bad+0)" rows where A1+A2 != total"}' "$f"
```

### 1.4 `OCM_heart/allelic_ratio/02_whole_chrX.R:50` — subset without the barcode check

`subset(heart, cells = chr_allelic_ratio$cell_barcode)` runs without
`check_barcode_match()`, which exists in `00_functions.R:284` precisely for this
step and which `04_core_escape.R:56` does call. On Seurat v5 a wrong barcode
parse errors with "Cannot find cells"; on v4 it silently drops them and every
downstream number is computed on the remainder.

Fix: insert before line 50

```r
check_barcode_match(chr_allelic_ratio$cell_barcode, colnames(heart),
                    what = paste0("02 whole chrX (", ALLELIC_RATIOS_FILE, ")"))
chr_allelic_ratio <- chr_allelic_ratio[chr_allelic_ratio$cell_barcode %in% colnames(heart), ]
```

### 1.5 `OCM_heart/allelic_ratio/02_whole_chrX.R:291 and :375` — monoallelic boundary is inconsistent

`frac_escaping` uses `AR <= 0.9` (escaping), `prop_tbl` uses `allelic_ratio < 0.9`
(biallelic), and 04/08 use `>= 0.90` for LOX. A cell at exactly 0.90 (9/10,
18/20, 27/30 reads, all common at `MIN_TOTAL_READS = 30`) is "escaping" in one
table and "not biallelic" in the other.

Fix: define `MONO_AR <- 0.90` in `00_functions.R` next to `MIN_TOTAL_READS`, and
use `AR < MONO_AR` for escaping/biallelic and `AR >= MONO_AR` for monoallelic
everywhere (02, 04, 08, and the `tile_ratio_map.R` OCM bins are already
`right = TRUE` so they agree with `>=`).

### 1.6 `spatial/ase_tile_sweep.R:618` — one SE is not deflated by the duplication factor

The script's rule (lines 122 to 139) is that every SE uses
`n_eff = n / EFF_DIVISOR`. The escape-by-set table breaks it:

```r
esc[, escape_se := sqrt(escape * (1 - escape) / umis)]
```

In a `dup` run `umis` are reads, so this SE is optimistic by `sqrt(EFF_DIVISOR)`
(about 3.4x at the measured factor of 11.8). The imprinted-control bounds at
lines 570 to 571 (`3 / (r + a)` and `qbeta(0.95, wrong + 1, n - wrong)`) have the
same problem.

Fix: divide `umis` and `r + a` by `EFF_DIVISOR` in those three places, and say
in the printed note that the imprinted `err` is per independent molecule.

### 1.7 `OCM_heart/allelic_ratio/08_lox_calling.R:135` — `optim()` convergence never checked

`fit_bb()` returns `plogis(o$par)` without looking at `o$convergence`. A
non-converged fit feeds `ab_all` and `ab_ref`, and every expected count in steps
3 and 4 inherits it silently.

Fix: `if (o$convergence != 0) warning("fit_bb did not converge: ", o$message)`
and carry `converged = o$convergence == 0` in the returned list; print it in
the `say()` lines that quote `mu` and `rho`.

## 2. Cluster pipeline robustness

### 2.1 `slurm/spatial_sinto_tiles.slurm:356` — raw tile BAMs are reused on count alone

The split is skipped when `N_BAM >= N_TILES`. The filtered set gets
`samtools quickcheck` (lines 477 to 484); the raw set in `$BAM_DIR` never does.
If sinto is killed while flushing its last files, the count is complete but the
tail files are truncated. `samtools view` then fails on them, the `&&` in the
filter command means no index is written, and the tile lands on the filter list
of every subsequent link forever without ever scoring.

Fix: when reusing, run `samtools quickcheck -v "$BAM_DIR"/*.bam`; if any fail,
delete the whole `$BAM_DIR` and let sinto re-split (sinto has no resume, so
say in the log that this costs a full split). Cheap to check, header-only.

### 2.2 Thirteen slurm scripts request 18 CPUs on `serial_std`

`spatial_sinto_tiles.slurm:57` records that `serial_std` caps `cpus-per-task`
at 16 and that sbatch rejects rather than trims. These files have
`--clusters=serial --partition=serial_std --cpus-per-task=18`:

- `slurm/Allelome.LINK.bodymap.slurm`
- `slurm/Allelome.LINK.pseudobulk.slurm`
- `slurm/Allelome.LINK.snRNAseq.slurm`
- `slurm/ML.analysis.slurm`
- `slurm/create_bigwig.slurm`
- `slurm/map_F1_TAC.slurm`
- `slurm/map_TAC_FACS.slurm`
- `slurm/map_adult_FACS.slurm`
- `slurm/map_adult_bodymap.slurm`
- `slurm/map_aged_FACS.slurm`
- `slurm/map_aged_bodymap.slurm`
- `slurm/mm10_analysis.slurm`
- `slurm/visualise_snRNAseq.slurm`

Fix: confirm the cap first with the command in `spatial_sinto_tiles.slurm:109`,
then set them to 16 (and check any `--threads 18` or `-p 18` inside the body).

### 2.3 `slurm/submit_sinto_tiles_chain.sh:28` — stale watch command

The comment says to watch with `squeue --clusters=cm4`, but the job header now
targets `--clusters=serial`. Change to `--clusters=serial`, or better,
`squeue -M serial,cm4`.

### 2.4 `slurm/cellranger_arc_multiome.slurm:106` — `--localmem=90` against `--mem=96G`

cellranger's `--localmem` is the amount its stages may request at once, not a
hard cap on the process tree. Martian plus the OS need headroom, and the cgroup
kills at 96G. 90 of 96 is tighter than the usual 10 to 15 percent margin. This
is a judgment call, not a bug, but an OOM kill 20 hours in costs a day.

Suggestion: `--localmem=84`. Leave `--localcores=16`.

### 2.5 `slurm/cellranger_arc_preflight.slurm:59,63` — `2>/dev/null` hides a missing input

A missing or unreadable fastq produces an empty subsample and the job proceeds
to run cellranger-arc for two hours on nothing. The read-count printout at
line 67 makes this visible only to a reader. The preflight already passed for
9w (job 5471258), so this matters for re-use on 78w.

Fix: test `[ -r "$f" ]` before each `zcat`, and after subsampling exit 1 if any
count is not `$NREADS`.

### 2.6 `slurm/sinto.slurm` and `slurm/sc_Allelome.PRO2.slurm` are edit-and-rerun scripts

Both are driven by commenting blocks in and out. `multiome/ANALYSIS_PLAN.md`
step 4 says the multiome allelic split will reuse them. Before the multiome
BAMs arrive, rewrite them in the style of `slurm/scAllelome_dedup.slurm`:
array task maps to sample, inputs come from positional args (the header uses
`--export=NONE`, so environment variables do not reach the job), and the task
list is built inside the job. Two multiome-specific points for that rewrite:

- The ATAC BAM from cellranger-arc has `CB` but no `UB`, so the UMI-based
  filters in `spatial/ase_bin_allele_counts.py` and the `-F 0x400` dedup logic
  need an ATAC branch that dedups on the duplicate flag only.
- The ATAC cell barcode sits at R2 positions 9 to 24 (plan, line 85). Nothing
  parsing the raw fastq may assume position 1.

## 3. Statistics and consistency

### 3.1 `spatial/tile_ratio_map.R:350` — z ignores the uncertainty in the autosomal ratio

`se <- sqrt(a_ratio * (1 - a_ratio) / x_n + auto_sd^2)` treats `a_ratio` as
known. It is estimated from `a_n` reads in the same tile. Add
`+ a_ratio * (1 - a_ratio) / a_n`. The effect is small (autosomes are about 7x
deeper) but it is free and makes the call panel slightly more honest.

### 3.2 `OCM_heart/allelic_ratio/11_autosomal_control.R:58` — the excess does not use the null floors it computes

`ar_excess <- ar_dom_x - ar_dom_a`, but the script then computes the depth-
dependent floor of `ar_dom` for each arm (`ar_dom_a_null`, `ar_dom_x_null`,
lines 69 to 70) and never uses them. The two arms differ in depth by roughly
7x, so their floors differ, and the raw difference is biased toward chrX.

Fix: add `ar_excess_adj <- (ar_dom_x - ar_dom_x_null) - (ar_dom_a - ar_dom_a_null)`,
report it beside `ar_excess`, and run the paired Wilcoxon on both.

### 3.3 `OCM_heart/allelic_ratio/02_whole_chrX.R:509` — the tricycle block always runs

Lines 425 to 459 mark it "NEGATIVE RESULT - DO NOT INTERPRET", yet it loads
`tricycle` and `SingleCellExperiment` (which mask dplyr verbs for the rest of
the script) and fits several hundred glmmTMB models every run.

Fix: wrap it in `if (Sys.getenv("RUN_TRICYCLE", "0") == "1") { ... }`. Keep
the depth block after it (line 823 on) running, since that result is live.
While there, consider `FPR_N_REPS` defaulting to 0 or 200 rather than 1000;
the comment at line 217 puts the default at hours.

### 3.4 `OCM_heart/allelic_ratio/00_functions.R:97` — warnings count as failures in the FPR simulation

`tryCatch(..., warning = function(w) NULL)` drops any replicate that warns
(glmmTMB convergence warnings are common on small groups) and reports it in
`n_failed`. That is conservative, and fine, but the FPR is then computed only
over clean replicates. Print `n_failed / n_reps` next to the FPR in
`calibrate_fpr()` so a 40 percent failure rate is visible.

## 4. Repository hygiene

### 4.1 `.claude/settings.json` — hook binary path is for a different user

Both `PreToolUse` hooks call `/Users/graylachlan/.local/bin/graphify`. The
current user is `lachlang` and neither `/Users/graylachlan/...` nor
`~/.local/bin/graphify` exists. Every Bash, Grep, Read and Glob call in Claude
Code fires a failing hook. Move machine-specific hooks to
`.claude/settings.local.json` (gitignored) or drop them.

### 4.2 `graphify-out/` is tracked (117 files, 1.6 MB, includes a cache)

Generated output. Add `graphify-out/` to `.gitignore` and
`git rm -r --cached graphify-out`. `graphify-out/cache/last_query_stamp` alone
appears in 19 recent commits as noise.

### 4.3 `.gitignore` ignores `CLAUDE.md`

There is therefore no committed project guidance for any Claude session. The
material already exists in `spatial/NEXT_ANALYSIS.md` (genotype, the three
meanings of `allelic_ratio`, the no_Xist bed rule) and `multiome/ANALYSIS_PLAN.md`.
Write a short `CLAUDE.md` pointing at those two files, stating the conventions
in section 5 below, and remove it from `.gitignore`. Also delete the duplicated
`Xist_TPM*.pdf` lines in `.gitignore`.

### 4.4 Absolute `source()` paths in every `OCM_heart/allelic_ratio/*.R`

`source("/dss/dssfs03/.../Postdoc/OCM_heart/allelic_ratio/00_functions.R")` is
repeated in 01 to 11. Replace with

```r
source(file.path(Sys.getenv("POSTDOC_ROOT", "/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Postdoc"),
                 "OCM_heart/allelic_ratio/00_functions.R"))
```

so the scripts can be syntax-checked off the cluster.

### 4.5 Smaller items

- `README.md` is one line. A paragraph per top-level directory would help
  any future reader (and any future model).
- `ASE/GTEx_LINK.R:16` has `setwd("/Users/graylachlan/LRZ Sync+Share/...")`.
- `T`/`F` used for `TRUE`/`FALSE` in `pSS_preprocessing.R:70`, `merge_strands.R:31-32`,
  `spatial/visium_hd_test.R` (several). Low risk, but `T` can be reassigned.
- Two `.r` files (`Xist_TPM.r`, `OCM_heart/visualising_data.r`) beside `.R`.
- `slurm/spatial_ase_sweep.slurm:250-266`: changing arg 3 (`WINDOW_TILE_UM`)
  invalidates the completeness check and triggers a full recount (hours). The
  window table is derived from the same pass, so this is by design, but the
  usage block at line 113 should say so.

## 5. Verified correct. Do not change these

These were checked and are either deliberate or right, even where they look odd.

- The sinto map has no header. `ase_tile_sweep.R:1045` writes it with
  `col.names = FALSE`; `tile_ratio_map.R:229` reads with `header = FALSE`;
  `spatial_sinto_tiles.slurm:323` counts tiles with `cut -f2 | sort -u | wc -l`.
  All three agree.
- `spatial_sinto_tiles.slurm:417`: the `sort` shim resolves the real binary
  with the shim directory removed from `PATH`, so a re-run cannot build a shim
  that execs itself.
- `submit_sinto_tiles_chain.sh:117`: `${res[@]+"${res[@]}"}` is the bash < 4.4
  `set -u` idiom, and `afterany` rather than `afterok` is deliberate (a
  TIMEOUT exit is the normal case).
- `spatial_sinto_tiles.slurm:686-708`: the chrX-free redo counter increments
  once per pass because the old table is deleted at line 707; `MAX_REDO = 2`
  means two real redos.
- `--export=NONE` in every slurm header means arguments must be positional.
  Do not "simplify" to environment variables.
- `ase_bin_allele_counts.py:1200`: the offset probe loads an unshifted SNP
  table on purpose. `--snp-offset -1` is measured (100.0 percent match), not
  assumed.
- `ase_bin_allele_counts.py:1460`: the loop variable is `sp`, not `pos`,
  because `pos` is the barcode-to-coordinate dict for all of `main()`.
- `04_core_escape.R:408-413`: `dbetabinom` matches glmmTMB's `Beta(mu*phi,
  (1-mu)*phi)` parameterisation and has a pmf sanity check.
- `08_lox_calling.R:461-465`: `mantelhaen.test` is fed a double array to avoid
  32-bit integer overflow in its variance term.
- `allelic_ratio` means `ar_dom` (0.5 to 1) in 02, 06 and 10, and `ar_a1`
  (0 to 1, directional) in 04. `10_build_ratio_table.R:64-68` documents this.
  Do not unify the column name without renaming both.
- `02_whole_chrX.R:503-508`: the `dplyr::` prefixes after `library(tricycle)`
  are load-bearing (namespace masking). Keep them.
- `ase_tile_sweep.R`: `rho_bb` and `C(d)` are fitted on raw read counts in
  dup mode and the script says so at lines 997 to 999. That is a documented
  limitation, not a bug.
- `multiome/libraries_9w.csv` and `libraries_78w.csv` list the four ATAC
  `iset` directories as four fastq rows of one sample. That is correct; the
  plan (lines 41 to 49) shows the isets are one library.

## 6. Suggested order of work

1. Section 1.1, 1.2, 1.4, 1.5, 1.6, 1.7 (small, local, no data needed).
2. Section 4.1 and 4.2 (five minutes, stops the hook noise and the cache churn).
3. Section 1.3: write the two cluster checks into `NEXT_ANALYSIS.md` for the
   user to run, then make `read_locus()` use `a1 + a2`.
4. Section 2.1, 2.3, 2.5, 3.1, 3.2, 3.3.
5. Section 2.2 after the user confirms the partition cap.
6. Section 2.6 before the multiome BAMs exist (the count job is queued now).
7. Section 4.3 and 4.4 last; they touch many files and are easiest to review
   as one commit.
