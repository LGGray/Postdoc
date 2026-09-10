# IPT seminar / lab meeting, 15 Sep 2026

Figures and slides for the results part of the talk. Everything reads summary
tables that already exist on the cluster (mounted at `/Users/graylachlan/cluster`)
and redraws them in one style; nothing here recomputes an analysis.

| File | What it does |
|---|---|
| `consolidate_gene_tables.py` | Concatenates the per-nucleus `Allelome.PRO2_all_genes` locus tables (chrX, all four samples) into one TSV, excluding scDblFinder doublets (`--keep-doublets` to keep them). ~4 min over the mount. |
| `make_figures.R` | Writes `F00`-`F20` PNGs (300 dpi), adult vs aged only; Sham/TAC versions of the biology figures go to `figures/sham_tac/`. `F09*` need the consolidated TSV. |
| `dump_seurat_meta.R` | Cluster-side one-off: dumps QC metrics, cell type and UMAP for every nucleus. Optional; slide 6 falls back to the existing tables without it (see below). |
| `crop_xist_lox_panel.py` | Re-crops `existing_figures/snRNA_Xist_by_LOX_VCM_adult_aged.png` out of panel b of the cluster PDF `core_escape_block_new_AR_Xist_umap_panel_VCM.pdf`. |
| `build_deck.py` | Appends the results slides (with speaker notes) to `IPT Seminar 15.9.26.pptx` using the deck's own layouts, and writes `IPT Seminar 15.9.26 - v2.pptx`. |

Outputs go to `~/LRZ Sync+Share/LGray/Presentations/IPT_Seminar_2026-09-15/`:
`figures/` (new), `existing_figures/` (cluster PDFs converted to PNG),
`data/` (consolidated table + logs), and the finished deck.

Run order, on the laptop:

```bash
python3 presentations/ipt_seminar_2026-09-15/consolidate_gene_tables.py "$HOME/LRZ Sync+Share/LGray/Presentations/IPT_Seminar_2026-09-15/data/all_genes_per_cell.tsv"
Rscript presentations/ipt_seminar_2026-09-15/make_figures.R
python3 presentations/ipt_seminar_2026-09-15/build_deck.py
```

Slide 6 (cell-type UMAP + QC panels, adult vs aged) is drawn from
`cutoff_sweep_cell_table.txt` (UMAP, UMI counts, all nuclei) and the cutoff_30
metadata (features, mito %, 86-89% of nuclei). To use every nucleus for all four
panels, run `dump_seurat_meta.R` on the cluster and copy the resulting
`OCM/seurat_metadata_umap.tsv` into `~/LRZ Sync+Share/LGray/Presentations/IPT_Seminar_2026-09-15/data/`,
then rerun `make_figures.R` and `build_deck.py`; the figure caption updates itself.

## Doublets

The snRNA-seq section is drawn from the doublet-free tree,
`OCM/Allelic_ratio_results_nodoublet/` (scDblFinder, `slurm/allelic_ratio_nodoublet.slurm`).
`make_figures.R` reads it by default; set `RESULTS_ROOT=Allelic_ratio_results`
to redraw the older doublet-containing figures for comparison. Removing
doublets drops 5.9% of nuclei at the 30-read chrX cutoff (worst in ventricular
CM, 11-15%) and lowers every escape estimate by 1-6 points; the age direction
is unchanged.

Two figures cannot be redrawn doublet-free. `F07` and `F08` read
`core_escape_genes_bayes_posterior_by_gene.txt` and
`core_escape_genes_pseudobulk_*`, which no script in this repo writes - they
date from 2026-07-07/08 and came from code that was never committed to the
modular pipeline. `make_figures.R` skips them with a message rather than
silently falling back to the doublet-containing tree, so those two slides keep
the 2026-09-08 PNGs. `F00` is also not doublet-filtered (its
`cell_counts_per_celltype_and_condition.txt` lives outside the results tree),
but it is not used in the deck.

Numbers quoted in slide text that are not drawn from a table are taken from
`spatial/NEXT_ANALYSIS.md` (2026-09-03 status: gene-body split of the escape
estimate, imprinted-locus error rates, C(d) expectations) and
`multiome/ANALYSIS_PLAN.md` (sequencing metrics).
