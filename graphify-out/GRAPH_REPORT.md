# Graph Report - .  (2026-07-20)

## Corpus Check
- Large corpus: 77 files · ~663,354 words. Semantic extraction will be expensive (many Claude tokens). Consider running on a subfolder.

## Summary
- 176 nodes · 235 edges · 18 communities (14 shown, 4 thin omitted)
- Extraction: 91% EXTRACTED · 9% INFERRED · 0% AMBIGUOUS · INFERRED: 21 edges (avg confidence: 0.78)
- Token cost: 79,117 input · 0 output

## Community Hubs (Navigation)
- Gene Synth UI State & Config
- SLE Classifier Model Scripts
- Gene Synth App + miRNA Datasets
- SLE Feature Selection & SHAP
- SLE Ensemble & Validation Pipeline
- Gene Synth Playback Engine
- Gene Synth FASTA Loading UI
- Allelome scRNA Error Parsing
- Gene Synth Sequence Sonification
- Gene Synth Helix Visualization
- PyTorch Learning Exercise
- SLURM Allelome Job Submission
- Cardiac miRNA Target Data
- Gene Synth Sequence Window Controls
- RNA-seq Experiment Inference
- ASE Conservation Visualization
- Repository Overview README

## God Nodes (most connected - your core abstractions)
1. `SLE Classification Pipeline README` - 15 edges
2. `handlePlay()` - 11 edges
3. `reportStatus()` - 10 edges
4. `playNotes()` - 9 edges
5. `SLE Paper Python Requirements` - 8 edges
6. `Gene Synth — FASTA to Monophonic Synth (index.html)` - 7 edges
7. `SLE Classification Pipeline (multi-model ML framework)` - 6 edges
8. `pushHelixEvent()` - 5 edges
9. `buildGenomeBrowser()` - 5 edges
10. `getOctaveRangeLabel()` - 5 edges

## Surprising Connections (you probably didn't know these)
- `Gene Synth — FASTA to Monophonic Synth (index.html)` --semantically_similar_to--> `SLE Classification Pipeline README`  [INFERRED] [semantically similar]
  fasta_synth/index.html → SLE_paper/README.md
- `Xist gene (X-inactive specific transcript)` --conceptually_related_to--> `Adult chrX Host Genes with miRNA Overlap (adult_host_miRNAs.txt)`  [INFERRED]
  fasta_synth/index.html → miRNA/adult_host_miRNAs.txt
- `Xist gene (X-inactive specific transcript)` --conceptually_related_to--> `Aged chrX Host Genes with miRNA Overlap (aged_host_miRNAs.txt)`  [INFERRED]
  fasta_synth/index.html → miRNA/aged_host_miRNAs.txt
- `Xist gene (X-inactive specific transcript)` --conceptually_related_to--> `Embryo chrX Host Genes with miRNA Overlap (embryo_host_miRNAs.txt)`  [INFERRED]
  fasta_synth/index.html → miRNA/embryo_host_miRNAs.txt
- `Xist gene (X-inactive specific transcript)` --conceptually_related_to--> `Young chrX Host Genes with miRNA Overlap (young_host_miRNAs.txt)`  [INFERRED]
  fasta_synth/index.html → miRNA/young_host_miRNAs.txt

## Import Cycles
- None detected.

## Hyperedges (group relationships)
- **chrX allelic-ratio & miRNA-target dataset series across mouse developmental stages** — mirna_adult_organ, mirna_adult_host_mirnas, mirna_aged_organ, mirna_aged_host_mirnas, mirna_embryo_organ, mirna_embryo_host_mirnas, mirna_young_organ, mirna_young_host_mirnas, mirna_host_mirnas [INFERRED 0.85]
- **Individual classifiers combined into the SLE voting ensemble** — sle_paper_logit_regression, sle_paper_random_forest, sle_paper_support_vector_machine, sle_paper_gradient_boosting_machine, sle_paper_multilayer_perceptron, sle_paper_voting_classifier [EXTRACTED 1.00]
- **Cardiac cell-type resolution allelic-ratio and miRNA-target dataset group** — mirna_cardiac_celltypes, mirna_cardiac_cell_types_host_mirnas, mirna_host_mirnas [INFERRED 0.75]

## Communities (18 total, 4 thin omitted)

### Community 0 - "Gene Synth UI State & Config"
Cohesion: 0.05
Nodes (35): activeVoices, BASE_COLORS, baseMap, browserBaseCells, browserCells, browserEl, COMPLEMENT_BASES, controls (+27 more)

### Community 1 - "SLE Classifier Model Scripts"
Cohesion: 0.10
Nodes (13): 01_model_comparison.R, 02_feature_concordance.R, 03_independent_validation.R, 04_chrX_feature_analysis.R, Gradient Boosting Machine classifier for SLE classification.  Trains a GBM model, Logistic Regression classifier for SLE classification.  Trains an Elastic Net-re, Multilayer Perceptron classifier for SLE classification.  Trains an MLP neural n, Plot overlaid Precision-Recall curves for all trained ML models.  Loads all indi (+5 more)

### Community 2 - "Gene Synth App + miRNA Datasets"
Cohesion: 0.12
Nodes (17): Gene Synth — FASTA to Monophonic Synth (index.html), Gene Sonification (FASTA-to-audio translation), Mobile silent-mode Web Audio warning rationale, bindRangeKnob(), bindSelectKnob(), setKnobAngle(), styles.css (Ableton-style rack UI), XIST.fasta (default reference sequence) (+9 more)

### Community 3 - "SLE Feature Selection & SHAP"
Cohesion: 0.14
Nodes (11): Boruta and Elastic Net feature selection pipeline.  Performs Boruta (Random Fore, boruta (Python package), matplotlib (Python package), numpy (Python package), pandas (Python package), pyreadr (Python package), scikit-learn (Python package), shap (Python package) (+3 more)

### Community 4 - "SLE Ensemble & Validation Pipeline"
Cohesion: 0.15
Nodes (12): analysis_pipeline.R (cleaned analysis script), Analysis Pipeline README, chrX_biomaRt.txt (chrX gene annotation), Nehar-Belaid Validation Cohort (nehar_belaid_update/), Perez Cohort Dataset (perez_update/), SLE_DisGeNet.tsv (SLE-associated gene set), Predict SLE status on an independent validation cohort.  Loads a pre-trained ens, Soft Voting Ensemble (+4 more)

### Community 5 - "Gene Synth Playback Engine"
Cohesion: 0.24
Nodes (11): advanceGenome(), buildGenomeBrowser(), clampRelease(), clampToOctaveWindow(), gatherParams(), getOctaveLimits(), midiToFreq(), playheadOffsetPx() (+3 more)

### Community 6 - "Gene Synth FASTA Loading UI"
Cohesion: 0.36
Nodes (9): dismissMobileReminder(), handleCustomFasta(), handleFileSelection(), loadDefaultFasta(), reportStatus(), resetGenomeBrowser(), setFastaStatus(), setPlayButtonState() (+1 more)

### Community 7 - "Allelome scRNA Error Parsing"
Cohesion: 0.50
Nodes (7): extract_first(), extract_sample(), main(), mean_ignore_nan(), parse_elapsed_seconds(), parse_file(), percentile_ignore_nan()

### Community 8 - "Gene Synth Sequence Sonification"
Cohesion: 0.33
Nodes (7): getRootMidi(), handlePlay(), resetHelix(), sanitizeSequence(), sequenceToNotes(), stopActiveVoices(), stopPlaybackNow()

### Community 9 - "Gene Synth Helix Visualization"
Cohesion: 0.33
Nodes (6): formatMidiNote(), getComplementBase(), getOctaveRangeLabel(), handleOctaveChange(), pushHelixEvent(), renderHelix()

### Community 11 - "SLURM Allelome Job Submission"
Cohesion: 0.53
Nodes (4): chunk_has_failures(), submit_sc_allelome_chunks.sh script, summarise_job(), wait_for_job()

### Community 12 - "Cardiac miRNA Target Data"
Cohesion: 0.50
Nodes (5): Cardiac Cell-Type Host Genes with miRNA Overlap, Cardiac Cell-Type Allelic Ratio Table (cardiac_celltypes.txt), Histogram: Number of Targets per miRNA (all genes), Histogram: Number of Targets per chrX miRNA, Aggregate Host Genes with miRNA Overlap (host_miRNAs.txt)

### Community 13 - "Gene Synth Sequence Window Controls"
Cohesion: 0.50
Nodes (4): getSequenceWindow(), handleSequenceWindowChange(), normalizeSequenceWindowInputs(), parsePositiveInt()

## Knowledge Gaps
- **58 isolated node(s):** `visualise_conservation.sh script`, `setup_directories.sh script`, `rangeOutputs`, `statusEl`, `playBtn` (+53 more)
  These have ≤1 connection - possible missing edges or undocumented components.
- **4 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `Gene Synth — FASTA to Monophonic Synth (index.html)` connect `Gene Synth App + miRNA Datasets` to `Gene Synth UI State & Config`, `SLE Classifier Model Scripts`?**
  _High betweenness centrality (0.390) - this node is a cross-community bridge._
- **Why does `SLE Classification Pipeline README` connect `SLE Classifier Model Scripts` to `Gene Synth App + miRNA Datasets`, `SLE Feature Selection & SHAP`, `SLE Ensemble & Validation Pipeline`?**
  _High betweenness centrality (0.350) - this node is a cross-community bridge._
- **Why does `SLE Classification Pipeline (multi-model ML framework)` connect `SLE Ensemble & Validation Pipeline` to `SLE Classifier Model Scripts`, `SLE Feature Selection & SHAP`?**
  _High betweenness centrality (0.143) - this node is a cross-community bridge._
- **What connects `visualise_conservation.sh script`, `setup_directories.sh script`, `rangeOutputs` to the rest of the system?**
  _58 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `Gene Synth UI State & Config` be split into smaller, more focused modules?**
  _Cohesion score 0.05128205128205128 - nodes in this community are weakly interconnected._
- **Should `SLE Classifier Model Scripts` be split into smaller, more focused modules?**
  _Cohesion score 0.1 - nodes in this community are weakly interconnected._
- **Should `Gene Synth App + miRNA Datasets` be split into smaller, more focused modules?**
  _Cohesion score 0.12418300653594772 - nodes in this community are weakly interconnected._