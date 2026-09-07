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
| GEX library | `25L008363` (S4) | `25L008364` (S5) |
| ATAC library | `25L008352` (S13-S16) | `25L008356` (S17-S20) |
| fastq prefix | `He_9w_multiome_XBxC_XX` | `25L008364` / `25L008356` |

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

Depth: roughly 650-700M GEX read pairs and ~330M ATAC read pairs per sample
(estimated from file sizes; the sequencing report has never been read). The
GEX is about 5x the 10x recommendation - deliberate over-sequencing for ASE.

That depth sets the compute problem. The count job runs on `serial_std`, not
cm4_tiny, because cm4 was refusing to place jobs on *fairshare* rather than on
resources - see the rationale in `slurm/spatial_sinto_tiles.slurm`, which
established that neither waiting nor a smaller shape helps in that situation.
The cost is 16 cpus / 96G instead of 112 / 256G, roughly a 6x wall-clock
penalty, which puts one sample past the 24h limit. That is handled rather than
avoided: cellranger-arc checkpoints per stage, so a timeout is resolved by
resubmitting the same array, and the sample that already finished exits in
seconds. Budget two rounds. If `sshare -U -u $USER` later shows EffectvUsage
decayed back toward NormShares, cm4_tiny becomes much the better home and the
script has the swap ready in comments.


## The genetic system

B6 mother x CAST father, female, **Xist deleted on the B6 allele**. Because
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
4. Allelic split of both modalities over the Xic-masked mm39 SNP set.
5. Peak-to-gene linkage, then ask whether escape genes differ from
   non-escape chrX genes in Xi accessibility and in enhancer usage.
