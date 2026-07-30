# Graph Report - /Users/graylachlan/Postdoc  (2026-07-30)

## Corpus Check
- Large corpus: 56 files · ~660,384 words. Semantic extraction will be expensive (many Claude tokens). Consider running on a subfolder.

## Summary
- 504 nodes · 744 edges · 52 communities (40 shown, 12 thin omitted)
- Extraction: 81% EXTRACTED · 14% INFERRED · 5% AMBIGUOUS · INFERRED: 104 edges (avg confidence: 0.89)
- Token cost: 669,685 input · 0 output

## Community Hubs (Navigation)
- Allelome.LINK Bulk Mapping Jobs
- miRNA Target & DGE Analysis
- ChIP-seq, miRNA & Coverage Jobs
- Single-Cell Allelome HPC Jobs
- pSS Preprocessing & IL2RG Power
- OCM Allelic Ratio Pipeline
- Link-DE Overlap & IGV Tracks
- CellRanger & CellBender Jobs
- Bayesian Allelic Ratio Priors
- DESeq2 & NMF Pathway Analysis
- SLE ML Job Array
- Female chrX Escape Posteriors
- Bodymap vs TAC Link Comparison
- GTEx Aging Link Analysis
- Cardiac TPM Normalisation
- Seurat to scVI Export
- Replicate Merge Consensus
- Pseudobulk Cell-Type DE
- Allelome.LINK Aging Transitions
- Link Gain/Loss Classification
- Adult/Aged snRNA-seq Processing
- Gene Correlation Rewiring
- Dispersion LRT & FPR Simulation
- scANVI Celltype Label Transfer
- Depth & Ambient AR Confounders
- Xist TPM Developmental Series
- Delta Allelic Ratio Overlaps
- BAM Strand Separation
- scAllelome Error Log Parsing
- Allelome.PRO2 Bayesian Update
- LOX Call & Xist repA Model
- FACS Link-DEG Overlap
- RefSeq Annotation Filtering
- Cardiac Xist Expression Matrices
- Beta-Binomial AR Model Helpers
- Core Escape Block Analysis
- PyTorch Learning Exercise
- SLURM Chunk Submission Driver
- Mechanism Permutation Nulls
- Cardiac Celltype miRNA Tables
- Conservation Visualisation Script
- Adult Organ miRNA Tables
- Aged Organ miRNA Tables
- Embryo Organ miRNA Tables
- CAST miRNA SNP Atlas
- Young Organ miRNA Tables
- Xist-AR Spearman Correlation
- Directory Archiving Job
- GTEx Link Presence Helper
- Xist Per-Celltype TPM Plot
- Repository Overview README

## God Nodes (most connected - your core abstractions)
1. `Bodymap Adult/Aged vs TAC Link Comparison` - 27 edges
2. `conda env RNAseq (LRZ miniconda3 24.7.1)` - 24 edges
3. `Allelome.LINK Aging Analysis Script` - 19 edges
4. `Allelome.PRO2.sh (Allelome.LINK/00_src, external toolkit)` - 19 edges
5. `SLURM job TAC_FACS - TAC cardiac RNAseq allelic mapping (clusters=serial, partition=serial_std, 18 cpus, 100G, 8h, array 1-32)` - 17 edges
6. `SLURM job adult_FACSs - adult cardiac RNAseq Allelome.PRO2/LINK (serial_std, 18 cpus, 100G, 8h, array 1-24)` - 15 edges
7. `SLURM job adult_Allelome_LINK - adult bodymap STAR/htseq-count (serial_std, 18 cpus, 100G, 8h, array 1-24)` - 15 edges
8. `DESeq2 Differential Expression Pipeline` - 14 edges
9. `SLURM job mm10_analysis - mm10 STAR mapping + Allelome.PRO2 + merge_strands (serial_std, 18 cpus, 50G, 8h, array 1-27)` - 14 edges
10. `SLURM job sc_Allelome.PRO2.chrX_gene - all cohorts (9w/78w/Sham/TAC) single-cell Allelome.PRO2 (cm4_tiny, 110 cpus, 150G, 24h)` - 14 edges

## Surprising Connections (you probably didn't know these)
- `Cardiac Cell Type Xist TPM Stage` --semantically_similar_to--> `Combined Xist Boxplot Stage`  [INFERRED] [semantically similar]
  Xist_TPM.r → cardiac_RNAseq_Xist.R
- `Median-AR / Min-Allelic-Score Replicate Consensus` --semantically_similar_to--> `Depth-Weighted Posterior AR Aggregation`  [INFERRED] [semantically similar]
  merge_replicates.R → Bulk_cardiac_bayesian_update.R
- `Escape definition: median allelic ratio <= 0.9` --semantically_similar_to--> `Fraction of cells escaping XCI (AR <= 0.9)`  [INFERRED] [semantically similar]
  miRNA/map_escape_to_miRNA.R → OCM_heart/allelic_ratio/02_whole_chrX.R
- `scANVI semi-supervised celltype classifier` --semantically_similar_to--> `Manual cluster-to-celltype annotation (11 celltypes)`  [INFERRED] [semantically similar]
  OCM_heart/scVI_cardiac_finetune.ipynb → OCM_heart/Seurat_preprocessing.R
- `CAST/BL6 SNP overlap with chrX mature miRNAs` --semantically_similar_to--> `core_escape_SNPs.R (SNP/annotation BED preparation)`  [INFERRED] [semantically similar]
  miRNA/miRNA_SNP.R → OCM_heart/core_escape_SNPs.R

## Import Cycles
- None detected.

## Hyperedges (group relationships)
- **Cardiac cell-type resolution allelic-ratio and miRNA-target dataset group** — mirna_cardiac_celltypes, mirna_cardiac_cell_types_host_mirnas, mirna_host_mirnas [INFERRED 0.75]
- **chrX allelic-ratio & miRNA-target dataset series across mouse developmental stages** — mirna_adult_organ, mirna_adult_host_mirnas, mirna_aged_organ, mirna_aged_host_mirnas, mirna_embryo_organ, mirna_embryo_host_mirnas, mirna_young_organ, mirna_young_host_mirnas, mirna_host_mirnas [INFERRED 0.85]
- **Allelic ratio pipeline: 00 helpers through 05 sensitivity** — ocm_heart_allelic_ratio_00_functions_module, ocm_heart_allelic_ratio_01_setup_stage, ocm_heart_allelic_ratio_02_whole_chrx_stage, ocm_heart_allelic_ratio_03_all_genes_stage, ocm_heart_allelic_ratio_04_core_escape_stage, ocm_heart_allelic_ratio_05_lox_sensitivity_stage [EXTRACTED 1.00]
- **Ambient RNA / technical artefact diagnostics for allelic ratio** — ocm_heart_allelic_ratio_03_all_genes_xist_ambient_probe, ocm_heart_allelic_ratio_03_all_genes_per_gene_discriminator, ocm_heart_allelic_ratio_02_whole_chrx_sequencing_depth_confounder, ocm_heart_allelic_ratio_04_core_escape_extreme_ar_calibration, ocm_heart_allelic_ratio_00_functions_n1_animal_caveat [INFERRED 0.85]
- **Seurat -> h5ad -> scVI/scANVI -> back to Seurat label transfer flow** — ocm_heart_seurat_preprocessing_module, ocm_heart_export_seurat_to_h5ad_module, ocm_heart_scvi_cardiac_finetune_scvi_reference_training, ocm_heart_scvi_cardiac_finetune_scanvi_classifier, ocm_heart_scvi_cardiac_finetune_scarches_surgery, ocm_heart_scvi_cardiac_finetune_reimport_to_seurat [EXTRACTED 1.00]
- **Aging Link Gain/Loss Concordance Pipeline** — ase_allelome_link_analysis_test_concordance_expr, ase_allelome_link_analysis_score_gained, ase_allelome_link_analysis_score_lost, ase_allelome_link_analysis_permute_mechanisms, ase_adult_aged_snrnaseq_adult_aged_log2fc_results_rds, ase_allelome_link_analysis_mechanism_match_perm_results_txt [EXTRACTED 1.00]
- **Aging vs TAC Injury Parallel Comparison Across Modalities** — ase_deseq2, ase_nmf, ase_gene_correlation_rewiring, ase_adult_age_bodymap_tac, ase_delta_ar, ase_overlapping_links_de, ase_deseq2_aging_vs_injury_parallel_concept [INFERRED 0.85]
- **Shared Noncoding-to-Coding Link Filtering Contract** — ase_adult_age_bodymap_tac_refseq_coding_noncoding_txt, ase_filter_bedpe, ase_compare_refseq, ase_adult_age_bodymap_tac, ase_deseq2, ase_adult_age_bodymap_tac_noncoding_coding_filter_concept [EXTRACTED 1.00]
- **Bayesian Allelic-Ratio Escape Gene Detection Flow** — allelome_pro2_bayes_bayes_function, bulk_cardiac_bayesian_update_female_posterior_stage, bulk_cardiac_bayesian_update_male_control_stage, bulk_cardiac_bayesian_update_escape_filter_stage, bulk_cardiac_bayesian_update_fisher_combination, merge_strands_locus_table_txt [EXTRACTED 1.00]
- **Strand-Aware Allelome.PRO Pipeline (BAM split to replicate consensus)** — separate_bam_strand_read_splitting_loop, separate_bam_strand_fwd_rev_bam_output, merge_strands_merge_locus_stage, merge_strands_merged_locus_table_txt, merge_replicates_median_consensus_summarisation, merge_replicates_locus_table_reps_txt [INFERRED 0.85]
- **Cross-Cohort Xist TPM Quantification** — cardiac_rnaseq_xist_get_gene_lengths, cardiac_rnaseq_xist_calculate_tpm, cardiac_rnaseq_xist_adult_mcl_stage, cardiac_rnaseq_xist_aged_mcl_stage, cardiac_rnaseq_xist_tac_mcl_stage, cardiac_rnaseq_xist_combined_boxplot_stage [EXTRACTED 1.00]
- **Allelome.PRO2 -> locus_table -> Allelome.LINK.R allelic-imbalance pipeline across HPC jobs** — slurm_allelome_pro, slurm_allelome_link_snrnaseq, slurm_allelome_link_pseudobulk, slurm_allelome_link_bodymap, allelome_pro2_sh, allelome_pro_locus_table, allelome_link_r, mouse_snp_snpfile_c57bl6_nj_cast_eij [INFERRED 0.85]
- **snRNA-seq preprocessing chain: cellranger -> cellbender -> ptrepack/Seurat -> per-cell BAM indexing** — slurm_cellranger_multi, slurm_cellranger_snrnaseq, slurm_cellbender, slurm_cellbender_ocm, slurm_index_igv_bams, ocm_hearts_ocm_per_sample_outs, adult_aged_heart_snrnaseq_dir [INFERRED 0.85]
- **SLE ML model zoo: one feature-selection step feeding six classifiers plus curve plotting, run over four feature sets** — slurm_ml_analysis, sle_paper_boruta_enet, sle_paper_logit_regression, sle_paper_random_forest, sle_paper_gradient_boosting_machine, sle_paper_support_vector_machine, sle_paper_multilayer_perceptron, sle_paper_voting_classifier, sle_paper_plot_all_curves [EXTRACTED 1.00]
- **Bulk RNA-seq allele-specific expression pipeline: STAR align, samtools index, Allelome.PRO2, Allelome.LINK, filter_bedpe** — star, samtools, allelome_pro2_sh, allelome_link_r, ase_filter_bedpe [EXTRACTED 1.00]
- **Single-cell chrX ASE cohort pipeline: sinto barcode split feeding per-cohort Allelome.PRO2 jobs (9w/78w/Sham/TAC)** — slurm_sinto, data_ocm_sinto_cell_bams, slurm_scallelome_9w, slurm_scallelome_78w, slurm_scallelome_sham, slurm_scallelome_tac, allelome_pro2_sh [INFERRED 0.95]
- **miRNA quantification pipeline: bowtie to miRBase then genome, samtools counting, bedtools tagging into merged count matrix** — bowtie, samtools, bedtools, ref_mirbase_mmu_mature_bed, data_mirna_counts_matrix [EXTRACTED 1.00]

## Communities (52 total, 12 thin omitted)

### Community 0 - "Allelome.LINK Bulk Mapping Jobs"
Cohesion: 0.11
Nodes (43): adult_aged_bodymap (dir_path.txt / dir_names.txt sample dirs), adult_aged_heart_snRNAseq/pseudobulk (adult.bam, aged.bam), Allelome.LINK.R (Allelome.LINK/00_src, external toolkit), Allelome.PRO2.sh (Allelome.LINK/00_src, external toolkit), Allelome.PRO locus_table.txt / locus_table_reps.txt output, BEDPE Noncoding-Coding Filter CLI, Data: AL_65.flt bedpe / links_full_table (Allelome.LINK output at allelic bias 0.65), Data: annotation_us_locus_table_reps.txt (merged-replicate locus table, Allelome.LINK input) (+35 more)

### Community 1 - "miRNA Target & DGE Analysis"
Cohesion: 0.06
Nodes (33): aged_host_miRNAs.txt, CAST_miRNA_SNPs.RData (loaded), combined_counts_gene_symbols.csv, <tissue>.edgeR.csv per-tissue result tables, ENSMUSG to gene symbol mapping via biomaRt, fgsea enrichment of miRDB miRNA target sets, m3.mirdb.v2025.1.Mm.symbols.gmt miRNA target sets, edgeR_DEG.R (HTseq -> edgeR DEG -> fgsea) (+25 more)

### Community 2 - "ChIP-seq, miRNA & Coverage Jobs"
Cohesion: 0.08
Nodes (28): bam2wig.py (RSeQC BAM -> WIG coverage), bamCoverage (deepTools, inverted-strand BigWig), bedtools (intersect), blacklist_GRCm39.bed (ENCODE blacklist regions), bowtie / bowtie-build short-read aligner (miRNA), bwa-mem2 (index + mem short-read aligner), cardiac_RNAseq/GRCm39_mapped *_Aligned.sortedByCoord.out.bam, cutadapt adaptor trimming (TruSeq Small RNA 3' adaptor) (+20 more)

### Community 3 - "Single-Cell Allelome HPC Jobs"
Cohesion: 0.19
Nodes (27): conda env RNAseq (LRZ miniconda3 24.7.1), Data: OCM/Allelome.PRO2_core_escape_new/<cohort> locus_table.txt and BED_files (SCRATCH staged), Data: OCM/Allelome.PRO2/<cohort> single-cell ASE output (task-list + locus tables), Data: Hearts_OCM cellranger per_sample_outs sample_alignments.bam (9w/78w/TAC/Sham snRNA), Data: OCM/sinto/<cohort> per-cell and per-celltype BAM files (sinto output, Allelome.PRO2 input), GNU parallel task-list executor, GRCm38.primary_assembly.genome.fa (mm10 reference FASTA), mm10_STAR_index (GRCm38.p6 + gencode vM25) (+19 more)

### Community 4 - "pSS Preprocessing & IL2RG Power"
Cohesion: 0.13
Nodes (19): Posterior Tail Probability at Null, CD8/CD4 IL2RG Analysis Driver Stage, Cohen's d Pooled Effect Size, merged_data In-Environment Data Frame, Post-hoc Power on Observed Effect, post_hoc_power_ttest, power_analysis_wilcox, sample_size_for_power (+11 more)

### Community 5 - "OCM Allelic Ratio Pipeline"
Cohesion: 0.13
Nodes (19): Allelic ratio shared helpers (00_functions.R), run_monoallelic_lrt, heart_seurat_object_SCT.rds, 01 Setup: Seurat load, VCM subclustering, QC plots, 02 Whole X chromosome allelic ratio analysis, whole_chr_allelic_ratios.txt (Allelome.PRO2 chrX table), whole_chr_cell_metadata.txt (handoff to 03 and 05), read_gene_table (+11 more)

### Community 6 - "Link-DE Overlap & IGV Tracks"
Cohesion: 0.12
Nodes (18): link_key (name_base|name_target|mechanism), get_fc_row, Link Mechanism (enhancing vs repressing), score_gained, score_lost, test_concordance_expr, LINKS_study/DEG/{celltype}_aged_vs_adult.txt and _TAC_vs_sham.txt, plot_volcano (+10 more)

### Community 7 - "CellRanger & CellBender Jobs"
Cohesion: 0.18
Nodes (15): adult_aged_heart_snRNAseq (9w/78w cellranger + cellbender + sinto BAMs), adult_aged_kidney_snRNAseq (fasterq-dump FASTQ output dir), cellranger 9.0.1 (count / multi), refdata-gex-GRCm39-2024-A (cellranger transcriptome), conda env cellbender, OCM/Hearts_OCM/outs/per_sample_outs raw feature bc matrices (9w, 78w, Sham, TAC), ptrepack (PyTables repack of h5 for Seurat), SLURM job cellbender (adult/aged heart snRNAseq; cm4_tiny, 18 cpus, 8G, 8h, array 1-2) (+7 more)

### Community 8 - "Bayesian Allelic Ratio Priors"
Cohesion: 0.18
Nodes (13): Allelic Ratio (theta), bayes_function, Beta-Binomial Conjugate Posterior Update, Beta(10,10) Autosomal Prior, Beta(0.5,0.5) Female chrX Prior, Signed Posterior Score (score_postP), Minimum SNP / Read Depth Restriction, weighted_score (+5 more)

### Community 9 - "DESeq2 & NMF Pathway Analysis"
Cohesion: 0.18
Nodes (13): LINKS_study/adult_aged_TAC_unique_link_overlaps.RData, DESeq2 Differential Expression Pipeline, apeglm (LFC shrinkage), fgsea (fora over-representation), adult_aged_bodymap/DEG/He_GEM.txt, limma (removeBatchEffect), F1_TAC_Sarah/DEG/TAC_GEM.txt, F1_TAC_Sarah/DEG/TAC_vs_sham.txt (+5 more)

### Community 10 - "SLE ML Job Array"
Cohesion: 0.17
Nodes (13): conda env SLE, python3 (ML script runner), SLE_paper/boruta_enet.py (feature selection), SLE_paper/gradient.boosting.machine.py, SLE_paper/logit.regression.py, SLE_paper/multilayer.perceptron.py, SLE_paper/plot.all.curves.py, SLE_paper/predict_independent_SLE.py (commented-out invocation) (+5 more)

### Community 11 - "Female chrX Escape Posteriors"
Cohesion: 0.20
Nodes (12): weighted_mean, Depth-Weighted Posterior AR Aggregation, Female Escape Gene Filter Stage, X-Linked Escape Gene, Female chrX Posterior AR Stage, Fisher's Method p_post Combination, Male Raw AR Control Stage, Male-Sample False Positive Control (+4 more)

### Community 12 - "Bodymap vs TAC Link Comparison"
Cohesion: 0.18
Nodes (12): Bodymap Adult/Aged vs TAC Link Comparison, biomaRt (MGI symbol to Ensembl ID), ComplexHeatmap + circlize, pair_key (name_base|name_target), adult_aged_heart_snRNAseq/adult.Allelome.LINK.RData, LINKS_study/snRNAseq_adult_TAC_links.txt, adult_aged_heart_snRNAseq/aged.Allelome.LINK.RData, LINKS_study/snRNAseq_aged_TAC_links.txt (+4 more)

### Community 13 - "GTEx Aging Link Analysis"
Cohesion: 0.20
Nodes (12): Mouse.CpG.clock.txt, GTEx ASE Link Aging Analysis, age_change_links_logit.RDS, emmeans (estimated marginal means), Epigenetic Clock Gene Overlap (Horvath / PhenoAge / Hannum / Mouse CpG), GTEx_ASE_links.txt, GTEx_Analysis_v10 gene TPM matrix, Hannum_clock.xlsx (+4 more)

### Community 14 - "Cardiac TPM Normalisation"
Cohesion: 0.20
Nodes (12): Aged MCL Counts + TPM Stage, mm39 RefSeq annotation.gtf, calculate_tpm, FID_comment Underscore Field Parsing, get_gene_lengths, Project_1050_lims_fixed.csv (TAC LIMS), Project_699_lims.csv (Aged LIMS), TAC MCL Counts + TPM Stage (+4 more)

### Community 15 - "Seurat to scVI Export"
Cohesion: 0.20
Nodes (12): ensure_rna_data, get_joined_layer_matrix, export_seurat_to_h5ad.R (Seurat -> h5ad for scVI), query_OCM.h5ad, AnnData X must hold raw counts for scVI, read_first_existing_rds, reference_adult_aged.h5ad, write_h5ad_counts_x (+4 more)

### Community 16 - "Replicate Merge Consensus"
Cohesion: 0.24
Nodes (11): Adult Bodymap Replicate Merge Stage, Aged Bodymap Replicate Merge Stage, Allelic Ratio itemRgb Colour Scheme, AR-Coloured ENCODE BED Track Writer, H3K27ac ChIP-seq Replicate Merge Stage, H3K4me3 ChIP-seq Replicate Merge Stage, annotation_locus_table_reps.txt, Median-AR / Min-Allelic-Score Replicate Consensus (+3 more)

### Community 17 - "Pseudobulk Cell-Type DE"
Cohesion: 0.20
Nodes (10): edgeR (cpm normalisation), Per-Cell-Type Pseudobulk Fold Change, pseudobulk_fc, LINKS_study/batch_corrected_VST_adult_aged_TAC_sham.txt, DESeq2 (Bioconductor), FACS cell-type GEM matrices (cardiac_RNAseq / aged_cardiac_RNAseq / TAC_cardiac_RNAseq DEG/*_GEM.txt), run_celltype_deseq, compute_vst (+2 more)

### Community 18 - "Allelome.LINK Aging Transitions"
Cohesion: 0.20
Nodes (10): Allelome.LINK Aging Analysis Script, adult.Allelome.LINK.RData, aged.Allelome.LINK.RData, ggalluvial (link state transition flows), is_lnc_refseq, figures/mechanism_match_perm_results.txt, pair_key (gene_base|gene_target), plot_link_flows (+2 more)

### Community 19 - "Link Gain/Loss Classification"
Cohesion: 0.22
Nodes (9): jaccard_index, Expression Loss vs Mechanism Loss, Gained / Lost / Retained Link Classification, link_key (gene_base|gene_target|mechanism), create_jaccard_matrix, jaccard_index, link_key (ncRNA | pcGene | Mechanism), Per-Link Logistic Age Model (+1 more)

### Community 20 - "Adult/Aged snRNA-seq Processing"
Cohesion: 0.25
Nodes (9): Adult/Aged Heart snRNA-seq Processing, adult_aged_log2FC_results.RDS, adult_aged_merged.RDS (integrated Seurat object), 9w/78w cell_index.txt (sinto barcode-to-celltype export), DoubletFinder, Seurat (snRNA-seq toolkit), stack_one (Allelome list to data.table), adult.Allelome.PRO2.RData (+1 more)

### Community 21 - "Gene Correlation Rewiring"
Cohesion: 0.22
Nodes (9): Aging vs Pressure-Overload Injury Parallel, Gene Correlation Rewiring Aging vs TAC, cor_matrix, LINKS_study/gene_correlation/{ct}/edge_rewiring.tsv, plot_delta_scatter, LINKS_study/gene_correlation/rewiring_summary.tsv, run_one_celltype, select_variable_genes (+1 more)

### Community 22 - "Dispersion LRT & FPR Simulation"
Cohesion: 0.33
Nodes (9): n=1 animal per condition caveat, run_dispersion_lrt, simulate_dispersion_null_fpr, Adult (9w) vs aged (78w) dispersion LRT, Sham vs TAC dispersion LRT, Excluding LOX cells before dispersion testing, 05 LOX sensitivity re-run of the dispersion tests, run_dispersion_lrt (monolith copy) (+1 more)

### Community 23 - "scANVI Celltype Label Transfer"
Cohesion: 0.22
Nodes (9): Ventricular cardiomyocyte subclustering (celltype_sub), cell_type_palette (matched to adult_aged_snRNAseq.R), Prediction uncertainty (entropy) filter, query_OCM_predictions.csv, Re-import scANVI predictions into Seurat, scANVI semi-supervised celltype classifier, scArches surgery fine-tuning on the OCM query, Manual cluster-to-celltype annotation (11 celltypes) (+1 more)

### Community 24 - "Depth & Ambient AR Confounders"
Cohesion: 0.22
Nodes (9): fit_depth, dplyr namespace masking by SingleCellExperiment/tricycle, Replication-timing dosage hypothesis (Xi late, Xa early), Sequencing depth as a confounder of allelic ratio, Tricycle cell-cycle position vs AR (negative result), all_genes_gene_df.txt, Xist allelic purity as ambient-contamination probe, SNPfile_C57BL_6NJxCAST_EiJ_sorted_mm39_no_Xist.bed (+1 more)

### Community 25 - "Xist TPM Developmental Series"
Cohesion: 0.22
Nodes (9): Developmental Age Series (embryo/young/adult/aged), Bodymap Xist TPM Extraction Stage, Cardiac Cell Type Xist TPM Stage, German-Labelled Figure Stage, log2(TPM + 1) Variance Transform, 43587_2025_856_MOESM2_ESM.xlsx Supplementary Data, Organ Colour Palette, Per-Organ Reshape to Long Data Stage (+1 more)

### Community 26 - "Delta Allelic Ratio Overlaps"
Cohesion: 0.25
Nodes (8): LINKS_study/adult_aged_TAC_link_overlaps.RData, annotation_no_predicted_locus_table_reps.txt (Allelome.PRO2 locus table), adult_aged_bodymap/bodymap_links.RData, F1_TAC_Sarah/TAC_SHAM_links.RData, Allelic Ratio (ASE), Delta Allelic Ratio Aged vs TAC, Delta Gene-Gene Correlation Rewiring, matrix_to_edges

### Community 27 - "BAM Strand Separation"
Cohesion: 0.29
Nodes (8): Forward/Reverse Locus Table Merge Stage, fwd/rev BAM Filename Convention, _fwd.bam / _rev.bam Output, BAM Header + samtools Pipe Setup, SAM Flag Read Splitting Loop, SAM Flag Bitmask Convention (0x1/0x10/0x40/0x80), samtools External Dependency, RSeQC Strand Rule Parsing

### Community 28 - "scAllelome Error Log Parsing"
Cohesion: 0.50
Nodes (7): extract_first(), extract_sample(), main(), mean_ignore_nan(), parse_elapsed_seconds(), parse_file(), percentile_ignore_nan()

### Community 29 - "Allelome.PRO2 Bayesian Update"
Cohesion: 0.29
Nodes (7): classify_gene (nested helper), classify_lost_link, Allelome.PRO2 Bayesian Update Script, bayes_function (Beta-Binomial posterior), Bayesian Beta-Binomial Allelic Bias, weighted_mean, weighted_score

### Community 30 - "LOX Call & Xist repA Model"
Cohesion: 0.29
Nodes (7): Bivariate AR x Xist UMAP multi-panel figure, dbetabinom, AR ~ 0 cells as an artefact calibration set, LOX call (loss of the inactive X, AR >= 0.9), Xist repA deletion mouse model (fixed maternal active X), Xist == 0 cross-tabulated against LOX calls, Four-gene confirmation of LOX status

### Community 31 - "FACS Link-DEG Overlap"
Cohesion: 0.33
Nodes (6): cardiac_RNAseq/adult_facs_links.RData, aged_cardiac_RNAseq/aged_facs_links.RData, TAC_cardiac_RNAseq/sham_facs_links.RData, adult_aged_bodymap/DESeq2/He_9w_vs_78w.txt, Link and Differential Expression Overlap Analysis, LINKS_study/DEG/facs_link_overlap/*.tsv

### Community 32 - "RefSeq Annotation Filtering"
Cohesion: 0.33
Nodes (6): Noncoding Base to Coding Target Link Filter, LINKS_study/Refseq_coding_noncoding.txt, GRCm39/GCF_000001635.27_GRCm39_genomic.gtf (RefSeq), RefSeq Annotation and Threshold Comparison, Predicted vs No-Predicted Annotation and 65/70 Ratio Sweep, GRCm39/gencode.vM37.primary_assembly.annotation.gtf

### Community 33 - "Cardiac Xist Expression Matrices"
Cohesion: 0.40
Nodes (6): Adult MCL Counts + TPM Stage, Cardiac Cell Type Codes (MP/CM/CF/EC), Combined Xist Boxplot Stage, ExprMat_TPM.txt Matrix, ExprMat.txt Raw Count Matrix, Gapdh Housekeeping Control Stage

### Community 34 - "Beta-Binomial AR Model Helpers"
Cohesion: 0.33
Nodes (6): Beta-binomial allelic ratio model (glmmTMB), fdr_to_stars, wilson_ci_cc, wilson_ci_depth, Beta-binomial model of AR ~ Xist expression, wilson_ci

### Community 35 - "Core Escape Block Analysis"
Cohesion: 0.33
Nodes (6): Majority-vote active-X allele orientation per cell, Per-gene discriminator: ambient contamination vs genuine escape, VCM_subcluster_per_gene_inactive_fraction.txt, core_escape_block_new_allelic_ratio_table.txt, 04 Core escape gene block analysis, core_escape_block_cell_metadata.txt

### Community 37 - "SLURM Chunk Submission Driver"
Cohesion: 0.53
Nodes (4): chunk_has_failures(), submit_sc_allelome_chunks.sh script, summarise_job(), wait_for_job()

### Community 38 - "Mechanism Permutation Nulls"
Cohesion: 0.50
Nodes (4): Mechanism-Shuffling Permutation Null, permute_mechanisms, acc_curve (link accumulation rarefaction), Link Accumulation (Rarefaction) Curve

### Community 39 - "Cardiac Celltype miRNA Tables"
Cohesion: 0.67
Nodes (3): Cardiac Cell-Type Host Genes with miRNA Overlap, Cardiac Cell-Type Allelic Ratio Table (cardiac_celltypes.txt), Aggregate Host Genes with miRNA Overlap (host_miRNAs.txt)

## Ambiguous Edges - Review These
- `infer_experiment.py` → `SLURM job adult_Allelome_LINK - adult bodymap STAR/htseq-count (serial_std, 18 cpus, 100G, 8h, array 1-24)`  [AMBIGUOUS]
  slurm/map_adult_bodymap.slurm · relation: references
- `infer_experiment.py` → `SLURM job adult_FACSs - adult cardiac RNAseq Allelome.PRO2/LINK (serial_std, 18 cpus, 100G, 8h, array 1-24)`  [AMBIGUOUS]
  slurm/map_adult_FACS.slurm · relation: references
- `infer_experiment.py` → `SLURM job mm10_analysis - mm10 STAR mapping + Allelome.PRO2 + merge_strands (serial_std, 18 cpus, 50G, 8h, array 1-27)`  [AMBIGUOUS]
  slurm/mm10_analysis.slurm · relation: references
- `04 Core escape gene block analysis` → `core_escape_block_cell_metadata.txt`  [AMBIGUOUS]
  OCM_heart/allelic_ratio/04_core_escape.R · relation: shares_data_with
- `CAST_miRNA_SNPs.RData (loaded)` → `CAST_mature_miRNA_SNPs.RData`  [AMBIGUOUS]
  miRNA/edgeR_DEG.R · relation: shares_data_with
- `LINKS_study/snRNAseq_adult_TAC_links.txt` → `gnomAD v4.1.0 chr13:30388771-30404967 variant export`  [AMBIGUOUS]
  ASE/gnomAD.R · relation: conceptually_related_to
- `GRCm39/gencode.vM37.primary_assembly.annotation.gtf` → `GRCm39/GCF_000001635.27_GRCm39_genomic.gtf (RefSeq)`  [AMBIGUOUS]
  ASE/DESeq2.R · relation: conceptually_related_to
- `cellTypist Label Import Stage` → `merged_data In-Environment Data Frame`  [AMBIGUOUS]
  power_calculation_IL2RG.R · relation: shares_data_with
- `SLURM job ML.analysis_complete (serial_std, 18 cpus, 20G, 8h, array 1-48)` → `SLE_paper/predict_independent_SLE.py (commented-out invocation)`  [AMBIGUOUS]
  slurm/ML.analysis.slurm · relation: references
- `SLURM job RNAseq (STAR index build + mapping; cm4_tiny, 18 cpus, 8G, 1h, array 1-24)` → `htseq-count (stranded read counting)`  [AMBIGUOUS]
  slurm/RNAseq.slurm · relation: calls
- `SLURM job RNAseq (STAR index build + mapping; cm4_tiny, 18 cpus, 8G, 1h, array 1-24)` → `samtools (view/sort/index/merge)`  [AMBIGUOUS]
  slurm/RNAseq.slurm · relation: calls
- `SLURM job F1_TAC (STAR map + htseq-count F1 TAC RNAseq; serial_std, 18 cpus, 100G, 8h, array 1-12)` → `Allelome.LINK.R (Allelome.LINK/00_src, external toolkit)`  [AMBIGUOUS]
  slurm/map_F1_TAC.slurm · relation: references
- `SLURM job F1_TAC (STAR map + htseq-count F1 TAC RNAseq; serial_std, 18 cpus, 100G, 8h, array 1-12)` → `Allelome.PRO2.sh (Allelome.LINK/00_src, external toolkit)`  [AMBIGUOUS]
  slurm/map_F1_TAC.slurm · relation: references
- `SLURM job F1_TAC (STAR map + htseq-count F1 TAC RNAseq; serial_std, 18 cpus, 100G, 8h, array 1-12)` → `merge_strands.R (merge forward/reverse locus tables)`  [AMBIGUOUS]
  slurm/map_F1_TAC.slurm · relation: references
- `SLURM job F1_TAC (STAR map + htseq-count F1 TAC RNAseq; serial_std, 18 cpus, 100G, 8h, array 1-12)` → `separate_BAM_strand.pl (split BAM into fwd/rev strand)`  [AMBIGUOUS]
  slurm/map_F1_TAC.slurm · relation: references
- `SLURM job F1_TAC (STAR map + htseq-count F1 TAC RNAseq; serial_std, 18 cpus, 100G, 8h, array 1-12)` → `STAR (RNA-seq aligner / genomeGenerate)`  [AMBIGUOUS]
  slurm/map_F1_TAC.slurm · relation: calls
- `Allelome.LINK.R (Allelome.LINK/00_src, external toolkit)` → `SLURM job testing_code - scratch/debug harness for merge_strands, merge_replicates and Allelome.LINK (cm4_tiny, 18 cpus, 8G, 1h)`  [AMBIGUOUS]
  slurm/testing.code.slurm · relation: references
- `Allelome.PRO2.sh (Allelome.LINK/00_src, external toolkit)` → `SLURM job adult_Allelome_LINK - adult bodymap STAR/htseq-count (serial_std, 18 cpus, 100G, 8h, array 1-24)`  [AMBIGUOUS]
  slurm/map_adult_bodymap.slurm · relation: calls
- `Allelome.PRO2.sh (Allelome.LINK/00_src, external toolkit)` → `SLURM job aged_Allelome_LINK - aged bodymap STAR/htseq-count (serial_std, 18 cpus, 100G, 8h, array 1-24)`  [AMBIGUOUS]
  slurm/map_aged_bodymap.slurm · relation: calls
- `merge_strands.R (merge forward/reverse locus tables)` → `SLURM job adult_Allelome_LINK - adult bodymap STAR/htseq-count (serial_std, 18 cpus, 100G, 8h, array 1-24)`  [AMBIGUOUS]
  slurm/map_adult_bodymap.slurm · relation: references
- `merge_strands.R (merge forward/reverse locus tables)` → `SLURM job aged_Allelome_LINK - aged bodymap STAR/htseq-count (serial_std, 18 cpus, 100G, 8h, array 1-24)`  [AMBIGUOUS]
  slurm/map_aged_bodymap.slurm · relation: references
- `merge_strands.R (merge forward/reverse locus tables)` → `SLURM job testing_code - scratch/debug harness for merge_strands, merge_replicates and Allelome.LINK (cm4_tiny, 18 cpus, 8G, 1h)`  [AMBIGUOUS]
  slurm/testing.code.slurm · relation: references
- `separate_BAM_strand.pl (split BAM into fwd/rev strand)` → `SLURM job adult_Allelome_LINK - adult bodymap STAR/htseq-count (serial_std, 18 cpus, 100G, 8h, array 1-24)`  [AMBIGUOUS]
  slurm/map_adult_bodymap.slurm · relation: references
- `separate_BAM_strand.pl (split BAM into fwd/rev strand)` → `SLURM job adult_FACSs - adult cardiac RNAseq Allelome.PRO2/LINK (serial_std, 18 cpus, 100G, 8h, array 1-24)`  [AMBIGUOUS]
  slurm/map_adult_FACS.slurm · relation: references
- `separate_BAM_strand.pl (split BAM into fwd/rev strand)` → `SLURM job aged_Allelome_LINK - aged bodymap STAR/htseq-count (serial_std, 18 cpus, 100G, 8h, array 1-24)`  [AMBIGUOUS]
  slurm/map_aged_bodymap.slurm · relation: references
- `separate_BAM_strand.pl (split BAM into fwd/rev strand)` → `SLURM job mm10_analysis - mm10 STAR mapping + Allelome.PRO2 + merge_strands (serial_std, 18 cpus, 50G, 8h, array 1-27)`  [AMBIGUOUS]
  slurm/mm10_analysis.slurm · relation: references
- `samtools (view/sort/index/merge)` → `SLURM job TAC_FACS - TAC cardiac RNAseq allelic mapping (clusters=serial, partition=serial_std, 18 cpus, 100G, 8h, array 1-32)`  [AMBIGUOUS]
  slurm/map_TAC_FACS.slurm · relation: calls
- `htseq-count (stranded read counting)` → `SLURM job TAC_FACS - TAC cardiac RNAseq allelic mapping (clusters=serial, partition=serial_std, 18 cpus, 100G, 8h, array 1-32)`  [AMBIGUOUS]
  slurm/map_TAC_FACS.slurm · relation: calls
- `SLURM job TAC_FACS - TAC cardiac RNAseq allelic mapping (clusters=serial, partition=serial_std, 18 cpus, 100G, 8h, array 1-32)` → `STAR RNA-seq aligner`  [AMBIGUOUS]
  slurm/map_TAC_FACS.slurm · relation: calls
- `SLURM job adult_FACSs - adult cardiac RNAseq Allelome.PRO2/LINK (serial_std, 18 cpus, 100G, 8h, array 1-24)` → `STAR RNA-seq aligner`  [AMBIGUOUS]
  slurm/map_adult_FACS.slurm · relation: calls
- `SLURM job adult_Allelome_LINK - adult bodymap STAR/htseq-count (serial_std, 18 cpus, 100G, 8h, array 1-24)` → `STAR RNA-seq aligner`  [AMBIGUOUS]
  slurm/map_adult_bodymap.slurm · relation: calls
- `SLURM job adult_FACS - aged (78w) cardiac RNAseq Allelome.PRO2/LINK (serial_std, 18 cpus, 100G, 8h, array 1-32)` → `STAR RNA-seq aligner`  [AMBIGUOUS]
  slurm/map_aged_FACS.slurm · relation: calls
- `SLURM job aged_Allelome_LINK - aged bodymap STAR/htseq-count (serial_std, 18 cpus, 100G, 8h, array 1-24)` → `STAR RNA-seq aligner`  [AMBIGUOUS]
  slurm/map_aged_bodymap.slurm · relation: calls
- `SLURM job RNAseq - cardiac miRNA bowtie/miRBase quantification (cm4_tiny, 18 cpus, 8G, 1h, array 1-12)` → `cutadapt adaptor trimming (TruSeq Small RNA 3' adaptor)`  [AMBIGUOUS]
  slurm/map_miRNA.slurm · relation: calls
- `SLURM job sinto - split cellranger BAMs into per-cell/per-celltype BAMs (cm4_tiny, 18 cpus, 64G, 8h, array 1-4: 9w/78w/TAC/Sham)` → `SLURM job spaceranger - placeholder batch script (file is currently empty, no directives or commands)`  [AMBIGUOUS]
  slurm/spaceranger.slurm · relation: conceptually_related_to
- `SLURM job testing_code - scratch/debug harness for merge_strands, merge_replicates and Allelome.LINK (cm4_tiny, 18 cpus, 8G, 1h)` → `merge_replicates.R - merge locus tables across replicates`  [AMBIGUOUS]
  slurm/testing.code.slurm · relation: references

## Knowledge Gaps
- **162 isolated node(s):** `visualise_conservation.sh script`, `Postdoc Repository README`, `Adult chrX Host Genes with miRNA Overlap (adult_host_miRNAs.txt)`, `Adult Organ Allelic Ratio Table (adult_organ.txt)`, `Aged chrX Host Genes with miRNA Overlap (aged_host_miRNAs.txt)` (+157 more)
  These have ≤1 connection - possible missing edges or undocumented components.
- **12 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **What is the exact relationship between `infer_experiment.py` and `SLURM job adult_Allelome_LINK - adult bodymap STAR/htseq-count (serial_std, 18 cpus, 100G, 8h, array 1-24)`?**
  _Edge tagged AMBIGUOUS (relation: references) - confidence is low._
- **What is the exact relationship between `infer_experiment.py` and `SLURM job adult_FACSs - adult cardiac RNAseq Allelome.PRO2/LINK (serial_std, 18 cpus, 100G, 8h, array 1-24)`?**
  _Edge tagged AMBIGUOUS (relation: references) - confidence is low._
- **What is the exact relationship between `infer_experiment.py` and `SLURM job mm10_analysis - mm10 STAR mapping + Allelome.PRO2 + merge_strands (serial_std, 18 cpus, 50G, 8h, array 1-27)`?**
  _Edge tagged AMBIGUOUS (relation: references) - confidence is low._
- **What is the exact relationship between `04 Core escape gene block analysis` and `core_escape_block_cell_metadata.txt`?**
  _Edge tagged AMBIGUOUS (relation: shares_data_with) - confidence is low._
- **What is the exact relationship between `CAST_miRNA_SNPs.RData (loaded)` and `CAST_mature_miRNA_SNPs.RData`?**
  _Edge tagged AMBIGUOUS (relation: shares_data_with) - confidence is low._
- **What is the exact relationship between `LINKS_study/snRNAseq_adult_TAC_links.txt` and `gnomAD v4.1.0 chr13:30388771-30404967 variant export`?**
  _Edge tagged AMBIGUOUS (relation: conceptually_related_to) - confidence is low._
- **What is the exact relationship between `GRCm39/gencode.vM37.primary_assembly.annotation.gtf` and `GRCm39/GCF_000001635.27_GRCm39_genomic.gtf (RefSeq)`?**
  _Edge tagged AMBIGUOUS (relation: conceptually_related_to) - confidence is low._