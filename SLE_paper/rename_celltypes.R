###### Rename cellTypist labels based on Perez et al. 2022

library(Seurat)
library(dplyr)
library(edgeR)

perez_celltypes <- list.files('../../SLE/pseudobulk/split_1/combined/ensemble/', pattern='.chrX.sav')
perez_celltypes <- gsub('.chrX.sav', '', perez_celltypes)
perez_celltypes <- gsub('_', ' ', perez_celltypes)
perez_celltypes[6] <- "CD4 T cell Effector-Memory"
perez_celltypes[10] <- "CD8 T cell Cytotoxic-GZMH"
perez_celltypes[11] <- "CD8 T cell Cytotoxic-GZMK"
perez_celltypes[20] <- "Non-classical monocyte"

levels(pbmc) %in% perez_celltypes


pbmc <- readRDS('pbmc.female.RDS')

Perez <- pbmc$cellTypist

Perez <- case_when(
  Perez == "NK cells" ~ "natural killer cell",
  Perez == "Naive B cells" ~ "B cell Naive",
  Perez == "B cells" ~ "B cell",
  Perez == "Classical monocytes" ~ "Classical monocyte",
  Perez == "Non-classical monocytes" ~ "Non-classical monocyte",
  Perez == "Tcm/Naive helper T cells" ~ "CD4 T cell Naive",
  Perez == "Tcm/Naive cytotoxic T cells" ~ "CD8 T cell Naive",
  Perez == "Tem/Trm cytotoxic T cells" ~ "CD8 positive, alpha beta T cell",
  Perez == "Tem/Effector helper T cells" ~ "CD4 positive, alpha beta T cell",
  Perez == "Regulatory T cells" ~ "CD4 T cell Treg",
  Perez == "MAIT cells" ~ "CD8 T cell MAIT",
  Perez == "DC1" ~ "Conventional DC",
  Perez == "DC2" ~ "Conventional DC",
  Perez == "CD16+ NK cells" ~ "NK cell Bright",
  Perez == "pDC" ~ "Plasmacytoid DC",
  Perez == "Age-associated B cells" ~ "B cell Atypical",
  Perez == "Plasmablasts" ~ "Plasmablast",
  Perez == "HSC/MPP" ~ "Progenitor cell",
  TRUE ~ Perez
)

pbmc@meta.data$Perez <- Perez