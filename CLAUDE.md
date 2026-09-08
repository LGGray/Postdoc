# Environment

- Files in this repo are edited locally on my Mac and pushed to GitHub. Editing files here is fine.
- All computation runs on the LRZ cm4 HPC cluster via SLURM, NOT on this machine
- Do NOT run the pipeline, R scripts, or job submissions in the local terminal, and don't treat local runs as representative of the cluster.
- When a step needs to execute, give me the command to run on the HPC rather than running it here. Local git is fine.

## graphify

This project has a knowledge graph at graphify-out/ with god nodes, community structure, and cross-file relationships.

Rules:
- For codebase questions, first run `graphify query "<question>"` when graphify-out/graph.json exists. Use `graphify path "<A>" "<B>"` for relationships and `graphify explain "<concept>"` for focused concepts. These return a scoped subgraph, usually much smaller than GRAPH_REPORT.md or raw grep output.
- If graphify-out/wiki/index.md exists, use it for broad navigation instead of raw source browsing.
- Read graphify-out/GRAPH_REPORT.md only for broad architecture review or when query/path/explain do not surface enough context.
- After modifying code, run `graphify update .` to keep the graph current (AST-only, no API cost).

# Read these first

Two documents carry the biology that decides what the numbers mean. Read them
before changing anything in `spatial/` or `OCM_heart/`:

- `spatial/NEXT_ANALYSIS.md` - genotype, the three different meanings of
  `allelic_ratio`, the no_Xist bed rule, and the open task list.
- `multiome/ANALYSIS_PLAN.md` - the cellranger-arc multiome pipeline.

The genotype, because every allelic number is read through it: B6 mother x CAST
father, with Xist deleted on the B6 allele. The CAST X is therefore the inactive
X in every cell, and CAST fraction on chrX is **escape** (~12.7%). Never infer
the model from the ratios - ask.

# Conventions that look wrong and are not

Do not "simplify" these. Each was checked and is deliberate:

- **`--export=NONE` in every slurm header**, so arguments must be positional.
  Do not convert them to environment variables.
- **The sinto tile map has no header.** `ase_tile_sweep.R` writes it with
  `col.names = FALSE`, `tile_ratio_map.R` reads with `header = FALSE`, and
  `spatial_sinto_tiles.slurm` counts tiles with `cut -f2 | sort -u | wc -l`.
  All three agree.
- **`allelic_ratio` means two different things.** It is `ar_dom` (0.5-1,
  undirected) in 02, 06 and 10, and `ar_a1` (0-1, directional) in 04.
  `10_build_ratio_table.R:64-68` documents this. Do not unify the column name
  without renaming both.
- **`dplyr::` prefixes in the tricycle block of `02_whole_chrX.R`** are
  load-bearing: SingleCellExperiment and tricycle mask several dplyr verbs.
- **`${res[@]+"${res[@]}"}`** in `submit_sinto_tiles_chain.sh` is the bash < 4.4
  `set -u` idiom, and `afterany` rather than `afterok` is deliberate - a TIMEOUT
  exit is the normal case for that chain.
- **`--snp-offset -1`** in `ase_bin_allele_counts.py` is measured (100% match),
  not assumed.
- **`MONO_AR` (0.90) in `OCM_heart/allelic_ratio/00_functions.R`** is the single
  monoallelic/LOX boundary. `AR < MONO_AR` is escaping/biallelic, `AR >=
  MONO_AR` is monoallelic. Do not reintroduce a per-script literal.

# Cluster

Read the LRZ partition table before writing or resizing a SLURM script:
https://doku.lrz.de/job-processing-on-the-linux-cluster-10745970.html
Do not copy core counts from neighbouring scripts. `serial_std` caps
`cpus-per-task` at 16 and sbatch **rejects** an over-request rather than
trimming it. `serial_long` gives 168h but caps the user at 100 GB across all
running jobs, so throttle arrays with `--array=1-N%1`.
