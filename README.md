# Postdoc

Allele-specific expression and X-chromosome inactivation escape in mouse heart.
B6 mother x CAST father, Xist deleted on the B6 allele, so the CAST X is the
inactive X and CAST fraction on chrX measures escape.

All computation runs on the LRZ cm4/serial cluster via SLURM. Nothing here is
meant to run on a laptop: every data path is under `/dss/dssfs03/...`.

| Directory | What is in it |
|---|---|
| `OCM_heart/` | snRNA-seq of mouse heart. `allelic_ratio/` is the numbered 00-11 analysis series: `00_functions.R` holds the shared constants and helpers, and each later script is a stage that reads the previous one's handoff. |
| `spatial/` | Visium HD spatial ASE. Tile-level allelic ratio maps, the tile-size sweep, and the per-bin allele counters. Start at `NEXT_ANALYSIS.md`. |
| `multiome/` | cellranger-arc joint RNA+ATAC. Start at `ANALYSIS_PLAN.md`. |
| `slurm/` | Every job script. Submit from the repo root: `sbatch slurm/<name>.slurm`. All use `--export=NONE`, so arguments are positional. |
| `ASE/` | Bulk and GTEx allele-specific expression. |
| `miRNA/`, `lncRNA_prediction/` | Earlier side projects: cardiac host-miRNA tables, and lncRNA classification from sequence. |
| `fasta_synth/` | A small standalone browser tool (plain HTML/JS) for synthesising FASTA sequences. |
| `SLE_paper/` | Machine-learning classification of systemic lupus from pseudobulk scRNA-seq, using chrX and SLE-associated gene features. Has its own README. |
| `graphify-out/` | Generated knowledge graph of this repo. See CLAUDE.md. |

See `CLAUDE.md` for the conventions that look wrong but are deliberate, and for
the cluster partition limits.
