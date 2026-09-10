# OCM hearts: per-cell-type pseudobulk allelic ratios

Rerun of Allelome.PRO2 on the OCM hearts pseudobulked **per cell type per
condition** instead of per nucleus, with scDblFinder doublets excluded, over
gene bodies on **all chromosomes** — so the question "are genes biallelic at
chromosome ends?" can actually be asked.

Modelled on the `adult_aged_heart_snRNAseq` pseudobulk run
(`slurm/Allelome.LINK.pseudobulk.slurm`): one BAM in, one `locus_table.txt`
out, per group.

## Why not just reuse `Allelome.PRO2_all_genes/`

That tree is **one Allelome.PRO2 run per nucleus**, and despite the directory
name it was scored against `annotation_us_mm39_chrX.bed` — chrX transcripts
only. It cannot answer a positional question about chromosome ends, and a
per-gene ratio from a single nucleus is mostly 0 or 1 by construction.

The new tree is `Allelome.PRO2_pseudobulk_celltype/`, 42 pseudobulks
(4 cohorts × 9–11 cell types), scored against
`annotation_us_mm39_gene_level.bed` — 39 881 gene bodies over 21 chromosomes.

## Run order

One command, from the repo root. Four stages chained with `--dependency`, so
it is fire and forget:

```bash
cd ~/Postdoc && git pull
./slurm/submit_pseudobulk_celltype_chain.sh
```

| stage | job | roughly |
| --- | --- | --- |
| `ids` | `pseudobulk_celltype_R.slurm ids` | minutes |
| `bams` | `pseudobulk_celltype_bams.slurm` | array 1-4, 1-3 h |
| `allelome` | `pseudobulk_celltype_allelome.slurm` | array 1-4, the long one |
| `analysis` | `pseudobulk_celltype_R.slurm analysis` | minutes |

Watch it:

```bash
squeue -M cm4 -u $USER -o '%.10i %.28j %.8T %.10M %.10l %.6D %R'
```

Arguments are positional (`--export=NONE` everywhere):
`submit_pseudobulk_celltype_chain.sh [merge|sinto] [RESULTS_ROOT] [TERMINAL_MB]`.

Re-run a single stage without redoing the rest — the normal way to widen the
terminal window after reading the coverage table:

```bash
STAGES=analysis ./slurm/submit_pseudobulk_celltype_chain.sh merge Allelic_ratio_results 10
```

`afterok`, not `afterany` — the opposite choice to
`submit_sinto_tiles_chain.sh` and for the opposite reason. Nothing here
resumes: a timed-out BAM build would leave the Allelome.PRO2 stage scoring a
partial cell type and reporting a number that looks fine and is wrong.

All four jobs are on `cm4`, including the two small single-threaded R stages.
An LRZ `--dependency` does not cross clusters, so putting the R stages on
`serial_std` would break the chain.

Two `cm4_tiny` limits the headers are pinned to, from the
[LRZ partition table](https://doku.lrz.de/job-processing-on-the-linux-cluster-10745970.html):

- **CPU range 17–112 physical cores.** 17 is a *floor*, and sbatch rejects
  anything below it with `QOSMinCpuNotSatisfied` rather than rounding up. So
  `pseudobulk_celltype_R.slurm` asks for 17 cores for a single-threaded job —
  that is deliberate, not slack to be trimmed. Do not take a core count from a
  neighbouring script: `index_igv_bams.slurm` asks for 4 and is itself below
  the floor.
- **Max 4 running jobs per user.** Both array stages are `--array=1-4`, which
  sits exactly at the cap, so if you have other `cm4_tiny` jobs running some
  array tasks will queue rather than start. That is throttling, not failure.

The stages can still be run by hand if you want to inspect between them:

```bash
cd /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/OCM
Rscript $POSTDOC_ROOT/OCM_heart/pseudobulk_celltype_ids.R
cd ~/Postdoc && sbatch slurm/pseudobulk_celltype_bams.slurm merge
                sbatch slurm/pseudobulk_celltype_allelome.slurm
cd /dss/.../OCM
Rscript $POSTDOC_ROOT/OCM_heart/allelic_ratio/13_pseudobulk_celltype_ar.R
```

Both R stages need the `seurat_env` R, not the repo-wide `RNAseq` env, because
stage 3 sources `allelic_ratio/00_functions.R`. The chain's `analysis` stage
pre-flights `tidyr` before claiming its allocation.

## Doublets

`pseudobulk_celltype_ids.R` reads
`Allelic_ratio_results/scDblFinder_per_cell.txt` and keeps only
`class == "singlet"`, so the doublets are gone **before** any read is assigned
to a cell type — no downstream filter is doing the work. Measured on the
current table: 1807 of 19 885 nuclei dropped (9.1%), 118/133/650/906 in
9w/78w/Sham/TAC. The script refuses to run if that table contains no doublets
at all, rather than quietly producing a "no doublet" result identical to the
unfiltered one.

## The two routes in step 1

Both produce the same reads.

- **`merge` (default).** Merges the per-cell BAMs already in `sinto/<cohort>/`,
  taking only singlet barcodes. Those per-cell sets are an exact 1:1 match to
  the QC'd cell list — 2376 / 2489 / 6927 / 8093 BAMs against the same counts
  in `scDblFinder_per_cell.txt`, verified with zero missing on both sides — and
  sinto assigns each read to exactly one barcode, so this is the same read set
  `filterbarcodes` would emit from a cell type map. Batched
  `samtools merge` (400 at a time) because the largest group is ~2000 files,
  past the open-file limit; inputs are already coordinate-sorted so the output
  needs no sort pass.
- **`sinto`.** `sinto filterbarcodes` on the cellranger per-sample BAM with the
  cell type map. Blocked right now: `Hearts_OCM/` exists only as
  `Hearts_OCM.tar.gz` (79 GB) and nothing under `OCM/` has
  `sample_alignments.bam` unpacked. The script checks and says so.

Cost of `merge`: ~30 GB of new BAMs (≈1.5 MB per nucleus). They are worth
keeping for IGV, but can be deleted once `Allelome.PRO2_pseudobulk_celltype/`
is written.

## Deliberate parameter choices

- **`annotation_us_mm39_gene_level.bed`**, all 21 chromosomes. The autosomes
  are not padding — they are the null. They carry the same B6-reference
  mapping bias and the same subtelomeric repeat content as chrX with no
  monoallelic biology, so a positional trend on chrX only means something if
  it beats the autosomal one.
- **`SNPfile_..._mm39_no_Xist.bed`**, per the house rule: Xist is deleted on
  the B6 X in this cross, so SNPs inside the Xist span cannot be read as
  allelic imbalance.
- **`-t 1`** (gene enters the table on one read). The depth cutoff lives in
  one place, `MIN_GENE_READS` in the R stage, so the same 42 tables can be
  re-thresholded without re-running anything. Read depth is the dominant
  confounder of an allelic ratio and must stay tunable.
- **The chr-prefix pre-flight in step 2 is load-bearing.** `Allelome.PRO2.sh`
  rewrites the chr prefix of the SNP and annotation files *in place*. With 11
  parallel workers on shared reference files that is a data race on a 630 MB
  bed. The pre-flight asserts BAM, SNP bed and annotation already agree, which
  makes those steps read-only no-ops. Do not remove it.

## Reading the output

Everything lands in `$RESULTS_ROOT/pseudobulk_celltype/`. Read
`pseudobulk_celltype_qc.txt` first: three of the 42 pseudobulks have fewer
than 20 nuclei (78w Epicardial-Mesothelial n=6, 78w stressed CM n=1, Sham
stressed CM n=12) and cannot support a per-gene ratio. `MIN_PB_CELLS`
(default 20) drops them from the models.

Column names follow the CLAUDE.md rule. Allelome.PRO2's own `allelic_ratio`
is renamed:

| column | meaning |
| --- | --- |
| `ar_a1` | A1/total, **directional**, 0–1. A1 is B6. |
| `ar_dom` | max(A1,A2)/total, **undirected**, 0.5–1. |
| `escape_frac` | `1 - ar_a1`, **chrX only** — the CAST X is the Xi, so this is escape. |
| `biallelic` | `ar_dom < MONO_AR` (0.90), the repo-wide boundary. |

`pseudobulk_direction_check.txt` is the sanity check to look at before
anything else: chrX must come out B6-dominant and the autosomes ~0.5. If it
does not, the allele assignment is wrong and nothing else means anything.

## The positional question, and what can confound it

"Biallelic at chromosome ends" is two different questions:

- **chrX** — every gene is expected monoallelic except escapers, so
  "biallelic" here *is* escape, and the question is whether escape
  concentrates towards the ends.
- **autosomes** — everything is biallelic already, imprinting aside, so any
  positional trend there is the technical baseline, not biology.

Three confounders are held explicitly, because each one can manufacture the
result on its own:

1. **Depth.** At `total_reads = 10` the smallest non-zero minor fraction is
   0.1, so shallow genes are pushed towards `ar_dom = 1` whether or not they
   are monoallelic. Depth is stratified (quintiles, `..._depth_matched.txt`)
   *and* one of the two matching variables in the terminal-vs-interior
   contrast.
2. **SNP density.** B6/CAST divergence is not uniform along a chromosome, and
   subtelomeric and pericentromeric sequence is repeat-rich. Counted once per
   gene (`gene_level_snp_count.bed`) and the second matching variable.
3. **Mouse chromosomes are telocentric.** Low coordinates are the centromere,
   not a telomere, and the annotation itself only starts ~3 Mb in. "Distance
   to nearest end" pools two very different places, so `dist_prox` and
   `dist_dist` are kept apart and reported separately in
   `pseudobulk_end_type_summary.txt`.

The test (`pseudobulk_end_vs_matched_interior.txt`) compares genes within
`TERMINAL_MB` (default 5) Mb of an end against **interior genes matched on
depth quintile and SNP-count tertile**, on the continuous ratio — escape
fraction on chrX, `ar_dom` on the autosomes. The reference distribution comes
from shuffling the terminal/interior label *within* strata, which is exactly
the null being claimed and assumes no distribution. `diff_matched > 0` means
higher at the ends.

There is deliberately **no regression of ratio on position**. A logistic fit of
the thresholded `biallelic` call on distance-to-end was written first and
removed: it throws away the continuous measurement in favour of a MONO_AR
indicator, and its p-value assumes genes are independent draws when terminal
genes are spatially clustered *by definition* — the hypothesis cannot also be
the null. On null synthetic data that GLM returned p = 2e-4; the matched
permutation returns a false-positive rate of 0.045 at α = 0.05 with finite
z ~ mean 0.05 / sd 0.99.

Its one real limitation, stated once: neighbouring genes share mappability and
reads, and shuffling within strata breaks that spatial correlation, so the
permutation p is somewhat optimistic. It is a reference, not a certificate.
The check that actually matters is `pseudobulk_end_sign_consistency.txt` —
whether the same sign shows up in independent cell types. There is no FDR
across the contrasts on purpose: they share genes and positions, so they are
not independent tests and a q-value would be arithmetic without meaning.

`pseudobulk_end_contrast_coverage.txt` says which contrasts were possible at
all, and it matters here: a 5 Mb window at the end of the 169 Mb X holds few
genes, and fewer pass the depth cutoff, so `matched_contrast()` declines any
pseudobulk with under 10 terminal or 20 interior genes. An absent row is a
pseudobulk that could not be measured, not one with no effect. If chrX comes
back mostly skipped, widen the window — `TERMINAL_MB=10` — rather than reading
the gaps as null results.

**n = 1 animal per condition.** Nothing here compares cohorts inferentially.
The positional models are within-pseudobulk; cohort differences are
descriptive only.
