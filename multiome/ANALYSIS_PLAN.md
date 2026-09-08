# Multiome heart (Eckstein/Engelhardt, Project 1729/1730)

## The data

Two 10x Multiome samples. Raw data lives in the shared lab area, not in
`go93qiw2`:

```
.../andergassen_lab/00_raw_data/Project_1729_..._snATAC   (56G, 48 fastq)
.../andergassen_lab/00_raw_data/Project_1729_..._snRNA   (127G,  8 fastq)
```

| | 9w (adult) | 78w (aged) |
|---|---|---|
| animal | `DAN-4106` | `DAN-0561` |
| GEX library | `25L008363` (S4) | `25L008364` (S5) |
| ATAC library | `25L008352` (S13-S16) | `25L008356` (S17-S20) |
| fastq prefix | `He_9w_multiome_XBxC_XX` | `25L008364` / `25L008356` |

Both pairings are **confirmed from the sample sheets**, not inferred. In
`Project_1730_lims.csv` and `Project_1729_lims.csv` the GEX and ATAC libraries
of each animal carry an identical `Sample_NameLIMS`:
`DAN-4106_He_9w_multiome_XBxC_-plus_XX` and
`DAN-0561_He_78w_multiome_XBxC_-plus_XX`. Protocol is recorded as
"Chromium Next GEM Single Cell Multiome ATAC + Gene Expression" for all ten
libraries. The `-plus` in those names is almost certainly a filesystem-safe
rendering of `-/+`, i.e. heterozygous for the deletion, matching the `het` in
the OCM sample names (`He_9w_snRNA_XBxC_het_XX`).

**Three other samples share this submission and are not ours.** The sample
sheets list five animals per project: `E74_REV_WT`, `E56_DMT_KO` and
`E68_KMA_KO` alongside the two `DAN-*` ones. Only the two `DAN-*` samples were
delivered into `andergassen_lab/00_raw_data` - the file counts confirm it (48
ATAC fastq = 2 libraries x 4 isets x 2 lanes x 3 reads; 8 GEX fastq = 2 x 2 x
2). The `E*` samples belong to the submitting group, which is why the project
directory is named for them; `DAN-*` is this lab's. So there is no missing
data, but do not be surprised by the extra rows in the sheet.

Two things about the directory layout are easy to get wrong:

- The ATAC `iset1..iset4` directories are the **four i7 index oligos of one
  sample index set**, demultiplexed separately. They are a single library, and
  treating them as four samples would quadruple the apparent n. This is not an
  inference from the directory names: the top 1000 cell barcodes of iset1
  overlap isets 2/3/4 at 966/967/964, i.e. all four sample the same nuclei.
  Four separate GEM wells would each draw ~8k barcodes from a 737k whitelist
  and overlap almost not at all. The four isets share one LIMS id and one
  24bp i5 spacer, and differ only on i7, which bcl2fastq consumed during
  demultiplexing and never wrote to a read file.
- The 78w sample was never given a descriptive rename, so its fastq prefix is
  the bare LIMS id. `libraries_78w.csv` accounts for this.

## Read geometry (verified from the fastqs, not assumed)

| | R1 | R2 | R3 |
|---|---|---|---|
| ATAC | 50 (genomic) | 24 (barcode) | 49 (genomic) |
| GEX | 101 | 101 | - |

All reads are uniform length - no trimming was applied.

Note a **documentation discrepancy**: the sample sheets record the intended
recipe as `50-8-24-49` for ATAC (R1-i7-i5-R2) and `28-10-10-90` for GEX. ATAC
was delivered exactly as recorded, and the `i7=8` there independently confirms
that the four isets differ on i7. GEX was **not** - the sheet says 28/90 but
the flowcell delivered 101/101. The recorded value is the protocol default, not
what was run, so trust the measurements below over the sheet.

**GEX used a symmetric 101/101 recipe, not the 10x multiome 28/90 recipe.**
The library itself is nonetheless genuine 10x 3' architecture, confirmed by
per-position base composition of R1:

```
pos  1-28   all four bases ~0.25   -> 16bp barcode + 12bp UMI
pos 29-52   T ~ 0.91               -> poly-dT
pos 53-60   T decays 0.86 -> 0.45  -> end of the 30nt dT tail, into cDNA
```

A poly-T block can only sit at position 29 if positions 1-28 are barcode+UMI,
so cellranger-arc reads the correct 28bp and the extra ~73 cycles on R1 are
simply wasted. Note the consequence: R1 carries no usable cDNA, so *all*
transcript sequence - and therefore all SNP overlap for ASE - comes from R2.
At 101bp that beats the 90bp spec, which is a real gain for allelic power.

**The ATAC cell barcode is at R2 positions 9-24, not 1-16.** Positions 1-8 are
a constant spacer (`CAGACGCG`, identical in all four isets). This is why 10x
specifies a 24-cycle i5 read to capture a 16bp barcode. cellranger-arc handles
the offset natively, but anything that parses these barcodes directly - which
the allelic ATAC work may well need to - must skip the first 8 bases or it
will read pure garbage.

Depth, exact, from `Project_1730_SequencingReport.csv`:

| GEX library | lane 1 | lane 2 | total pairs |
|---|---|---|---|
| 9w  `25L008363` | 339,538,551 | 347,544,986 | **687,083,537** |
| 78w `25L008364` | 364,239,443 | 369,170,677 | **733,410,120** |

ATAC is still estimated at ~330M pairs per sample; `Project_1729`'s report has
not been read.

The over-sequencing is unambiguously deliberate. On the same flowcell the three
samples belonging to the submitting group got 77-127M pairs each, while these
two got 687M and 733M - five to nine times the depth. At ~8k nuclei that is
roughly 86k read pairs per nucleus against the 20-25k 10x recommends, which is
the depth ASE needs and not an accident.

**Read 1 quality is much lower than Read 2, and this is expected rather than a
problem.** The report gives R1 at 70-73% bases >= Q30 (mean Q ~33.6) against
R2 at 96-97% (mean Q ~39.2). That asymmetry follows directly from the geometry
above: R1 spends positions 29-58 in a poly-T homopolymer and the rest reading
through into low-complexity sequence, both of which degrade base calling, and
it is also why R1 compresses worse than R2 despite equal length. What matters
for the pipeline is only R1 positions 1-28, which are read before phasing
degrades. Confirm that rather than assume it: check the "Q30 bases in barcode"
and "Q30 bases in UMI" metrics in the cellranger-arc summary, since poor
barcode quality would cost cells at the calling stage.

That depth sets the compute problem. The count job runs on `serial_std`, not
cm4_tiny, because cm4 was refusing to place jobs on *fairshare* rather than on
resources - see the rationale in `slurm/spatial_sinto_tiles.slurm`, which
established that neither waiting nor a smaller shape helps in that situation.
The cost is 16 cpus / 96G instead of 112 / 256G. Whether that fits inside the
24h limit is genuinely uncertain: a generic throughput estimate says it is
marginal, but the closest local precedent argues it is fine -
`slurm/cellranger_multi.slurm` counted four OCM snRNA samples on 18 cpus
within 24h. Local calibration beats the generic estimate, so 24h on
`serial_std` is the working assumption.

The timeout path therefore exists as insurance, not as an expectation.
cellranger-arc checkpoints per stage, so if a run is killed, resubmitting the
same array resumes it and the sample that already finished exits in seconds.
`serial_long` (168h, `--qos=cm4_serial_long`) is the fallback if 24h turns out
to be short, at the cost of running the samples sequentially - it caps the
user at 100G across all running jobs. And if `sshare -U -u $USER` later shows
EffectvUsage decayed back toward NormShares, cm4_tiny becomes much the better
home; the script has that swap ready in comments.


## Preflight QC (job 5471258, 1M-read subsample of 9w)

`Pipestance completed successfully` on cellranger-arc 2.2.0 against GRCm39 -
not merely past validation but through the whole pipeline including BAMs.
Recorded here because the run lived on scratch and will be purged; a copy of
`summary.csv` / `web_summary.html` is in
`adult_aged_multiome/preflight_outs/`.

**Depth-independent metrics, i.e. the ones that carry over to the full run:**

| metric | value | reading |
|---|---|---|
| GEX Q30 in barcode | 0.9719 | resolves the R1 quality worry |
| GEX Q30 in UMI | 0.9720 | ditto |
| GEX valid barcodes | 0.9114 | genuine `737K-arc-v1` multiome library |
| ATAC valid barcodes | 0.9808 | barcode read at the right 8bp offset |
| ATAC Q30 in barcode | 0.9211 | fine |
| ATAC TSS enrichment | 12.05 | good chromatin signal (>10 is healthy) |
| ATAC confidently mapped | 0.8246 | healthy |
| ATAC non-nuclear (chrM) | 0.0227 | low - nuclei prep worked |
| GEX mapped to genome | 0.9015 | healthy |
| GEX mapped to transcriptome | 0.6260 | good for single-nucleus |
| GEX intronic / exonic | 0.4565 / 0.3923 | intronic > exonic, as nuclei should be |
| GEX total genes detected | 13346 | good complexity from only 1M reads |

The barcode result is the one that mattered. The sequencing report had R1 at
only 70-73% bases >= Q30 across its full 101bp, which raised the prospect of
losing cells at the calling stage. At 97.2% Q30 over positions 1-28, the
degradation is confined to the poly-T and read-through region exactly as the
base-composition profile predicted, and the informative part of R1 is clean.

The two valid-barcode figures together close the last chemistry question. A 3'
library mismatched against the multiome whitelist would score near zero rather
than 91%; and had cellranger-arc taken ATAC positions 1-16 (spacer plus 8bp of
barcode) instead of 9-24, its 98% would likewise have been near zero.

One number to keep an eye on: **GEX reads mapped antisense to gene = 0.2203**.
Typical is nearer 5-15%. Elevated antisense is not unusual for single-nucleus
data, but if it persists at full depth it is worth understanding before
interpreting expression, and it interacts with allelic counting because
antisense reads over a SNP still carry allele information.

**Ignore these from the preflight** - all depth-limited at 1M reads and
meaningless here: estimated cells (363), feature linkages (3), reads/UMIs/genes
per cell, number of peaks (1176), fraction of genome in peaks (0.0004),
fraction of fragments overlapping peaks, and both duplicate rates, which will
rise steeply at 687-733M reads.

Disk is not a constraint: dssfs03 had 18T free against the ~250G the two runs
need. Note `df` reports the filesystem rather than the DSS container quota, so
that is the number to check if writes ever fail.

## The genetic system

B6 mother x CAST father, female, **Xist deleted on the B6 allele** -
confirmed directly, not carried over from the spatial/snRNA work, because the
direction of this one fact determines the sign of every result below. Because
Xist is required in *cis* to inactivate its own chromosome, the B6 X cannot be
silenced. So the CAST X is the inactive X in **100% of nuclei** - fully
skewed, non-mosaic, no patch structure. That is what makes the readouts clean:

- **RNA**: CAST fraction on chrX *is* the escape level. Prior snRNA/spatial
  work in this project puts it near 12.7%, and flat between 9w and 78w.
- **ATAC**: CAST fraction over chrX peaks is accessibility *of the inactive
  X*. This is the axis multiome adds and RNA alone cannot give. Expect a low
  genome-wide CAST fraction on chrX with focal accessibility at escape genes.

Two corrections that must be applied, both of which already have assets in
`GRCm39/`:

1. **Mask the Xic.** The Xist deletion removes B6 sequence, so that interval
   reads as 100% CAST in both modalities - a deletion artefact, not escape.
   Use the existing `SNPfile_..._mm39_no_Xist.bed` /
   `SNPfile_..._mm39_no_Xic500kb.bed` builds rather than the unmasked SNP set.
2. **Reference mapping bias.** Reads are mapped to a B6 reference, so CAST
   reads align slightly less well and the CAST fraction is biased low. The
   autosomal control in `OCM_heart/allelic_ratio/11_autosomal_control.R`
   already estimates this; the same correction applies here, and needs
   re-estimating separately for ATAC because the read length differs.

## Count results, and what they change

Both samples completed on cellranger-arc 2.2.0. BAMs are 23G/21G (9w
gex/atac) and 26G/26G (78w) under `adult_aged_multiome/<id>/outs/`.

| | 9w | 78w |
|---|---|---|
| GEX read pairs | 687,083,537 | 733,410,120 |
| **GEX duplicates** | **98.58%** | **97.79%** |
| non-duplicate reads | ~9.8M | ~16.2M |
| estimated nuclei | 1,758 | 3,822 |
| median UMI/nucleus | 1,353 | 863 |
| GEX valid barcodes | 0.9128 | 0.9056 |
| GEX Q30 barcode / UMI | .973 / .976 | .971 / .975 |
| GEX antisense | 0.2199 | 0.1850 |
| ATAC read pairs | 412,392,071 | 444,372,323 |
| ATAC duplicates | 77.74% | 66.93% |
| ATAC median HQ fragments/nucleus | 9,994 | 7,820 |
| ATAC TSS enrichment | 5.96 | 7.42 |
| ATAC FRiP (fragments in peaks) | 0.1413 | 0.2047 |
| ATAC peaks | 52,231 | 71,694 |
| feature linkages | 8,675 | 40,966 |

**Deduplication dominates everything.** `scAllelome_dedup.slurm` measured 39.4%
duplicates in the OCM data and expected `total_reads` to "roughly halve". Here
it is 98.58% - a ~70x reduction. The depth bought amplification, not
molecules: ~390k reads per nucleus collapsed to ~1,350 UMIs, because library
complexity was the ceiling rather than sequencing. Consequences:

- Read-level Allelome.PRO2 is not a variant worth keeping for this dataset. A
  single amplified molecule would vote ~70 times, which is exactly the failure
  `09_dedup_comparison.R` describes ("three real molecules and heavy
  amplification on one reads as near-monoallelic").
- **Every threshold calibrated on the read-level tree must be re-derived**:
  `MIN_TOTAL_READS` in `02_whole_chrX.R`, the sweep grids in `06`/`07`, the
  exact-count stratification in `08`. They were tuned where dedup halved
  counts.

**MAPQ encodings differ between the two modalities.** GEX is STAR (unique =
255); ATAC is BWA (caps at 60). Reusing `-q 255` on the ATAC BAM returns zero
reads *without erroring* - use `-q 30`. Verified on the BAMs themselves.

**Per-nucleus allelic ratios will be thin; pseudobulk is the level to work
at.** chrX carries 626,029 SNPs over 169.5 Mb in the `_no_Xist` build, i.e. one
per 271 bp, so a 101 bp read is informative ~37% of the time. At 1,353 median
UMIs with chrX order 5% of expression, that is roughly 20-25 informative chrX
molecules per nucleus - a proportion whose standard error is near 0.10, which
is the noise regime that made the AR >= 0.90 population in `08` hard to read.
Per cell type, an abundant population gives ~11k informative molecules, which
is ample.

**The age contrast has gained a technical confound on top of n=1.** 78w is the
better sample on nearly every axis - 2.2x the nuclei, 1.7x the unique
molecules, better ATAC on both TSS enrichment and FRiP, and 4.7x the feature
linkages. So a 9w-vs-78w difference is confounded with data quality, and one
animal per group cannot separate them. It was already descriptive-only; it now
has a named alternative explanation.

**ATAC is the weaker modality.** TSS enrichment of 5.96/7.42 is modest (the
preflight's 12.05 was a selection effect - only the 363 best nuclei were
called at 1M reads), and 14-20% FRiP is low against a typical 40-70%. The
allelic-ATAC aim is the most novel and the least well supported by the data,
so RNA escape should be made solid first.

## Statistical ceiling - read this before designing any test

**n = 1 animal per age.** One 9w, one 78w. The 9w-vs-78w contrast therefore
has no biological replication and supports **no inferential test of the age
effect** - the same limitation already documented for the OCM_heart data,
where dispersion and monoallelic LRTs on n=1 gave false positive rates of
8-81%. Report the age comparison descriptively.

What *is* testable: contrasts *within* a sample, where nuclei provide
replication - escape gene vs non-escape gene, cell type vs cell type, RNA vs
ATAC concordance. Those are the claims to build the analysis around. Anything
framed as "escape changes with age" is a two-animal anecdote and must be
worded as one.

A further power constraint specific to allelic ATAC: fragments are only
assignable to an allele if they overlap a B6/CAST SNP. With 50/49 bp reads
that is a minority of fragments, so effective allelic depth is well below
nominal. Measure the assignable fraction from the mm39 SNP bed before
interpreting any per-peak ratio - `GRCm39/imprinted_snp_density.tsv` shows the
density calculation already exists in the project.

## Steps

0. Install `cellranger-arc` + `refdata-cellranger-arc-GRCm39-2024-A`. DONE:
   cellranger-arc 2.2.0 with the 2024-A GRCm39 reference.
   GRCm39 matters: it keeps coordinates compatible with every existing
   `*_mm39_*.bed` in `GRCm39/`. An mm10 reference would orphan all of them.
   Verified: the reference uses `chr`-prefixed contigs (`chr1`, `chrX`,
   `chrM`), matching `SNPfile_..._mm39_*.bed`, `chr_annotation_mm39.bed` and
   `mm39.chrom.sizes`. Ensembl-style bare names would have produced zero
   allelic counts with no error at all.
   One mismatch to reconcile later: the ARC reference is annotated with
   GENCODE **vM33**, whereas the existing GEX work uses **vM37**
   (`refdata-gex-GRCm39-2024-A`, `gencode.vM37...gtf` in `ASE/DESeq2.R`).
   Coordinates are identical so SNP beds are unaffected, but gene IDs and
   names drift between releases - expect some when comparing multiome cell
   types against the existing snRNA object.
1. `slurm/cellranger_arc_preflight.slurm` - validate GEX read geometry.
2. `slurm/cellranger_arc_multiome.slurm` - count both samples (array 1-2).
3. Joint QC + WNN cell typing, against the cardiac labels already used in
   `OCM_heart/Seurat_preprocessing.R` so cell types are comparable.
4. **Pseudobulk allelic ratio, no cell typing needed** -
   `slurm/multiome_allelome_pseudobulk.slurm`. Filters the GEX BAM to
   `-F 0x400 -q 255`, then runs Allelome.PRO2 against
   `chr_annotation_mm39.bed` and the `_no_Xist` SNP build. That annotation is
   one row per chromosome, so a single run returns both halves of the
   question: chrX gives the escape level, and chr1-chr19 give the
   mapping-bias baseline under identical filters and SNP set. Autosomes should
   sit at 0.5; the offset is the correction chrX needs, which supersedes
   borrowing OCM's estimate from `11_autosomal_control.R`.
   **Gate:** does chrX reproduce the ~12.7% from the snRNA/spatial work? If
   not, stop and find out why before building anything on top.
5. Joint WNN cell typing, against the cardiac labels in
   `OCM_heart/Seurat_preprocessing.R`. Everything per-celltype blocks on this,
   because `sinto filterbarcodes` needs a barcode->group `cell_index.txt`.
6. Per-celltype Allelome.PRO2 via `sinto`, reusing the `scAllelome_*.slurm`
   GNU-parallel pattern - but those target cm4_tiny at 110 cpus and will need
   the serial_std 16-cpu treatment.
7. Allelic ATAC. New work: peaks as the annotation rather than genes, `-q 30`
   not `-q 255`, and no UMIs so deduplication is positional only.
8. Peak-to-gene linkage, then ask whether escape genes differ from
   non-escape chrX genes in Xi accessibility and in enhancer usage.
