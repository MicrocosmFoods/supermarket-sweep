library(tidyverse)
library(BiocManager)
library(ComplexHeatmap)
library(colorRamp2)

#################################
# Plotting metabolomics top detected metabolites heatmaps
#################################

# sample metadata
sample_metadata <- read_tsv("metadata/Combined master sheet updated 042126-8ed4f7f4.tsv")

MAJOR_MAP <- c(
  "Dairy" = "Dairy", "Legume - Soy" = "Legume / Soy", "Legume - other" = "Legume / Soy",
  "Legume" = "Legume / Soy", "Soy" = "Legume / Soy",
  "Grain" = "Grain", "Grain/legume" = "Grain", "Grain/Legume" = "Grain",
  "Vegetable - Cabbage" = "Vegetable", "Vegetable - Beet" = "Vegetable",
  "Vegetable - Cucumber" = "Vegetable", "Vegetable - other" = "Vegetable",
  "Vegetable/grain" = "Vegetable", "Vegetable/Grain" = "Vegetable",
  "Vegetable" = "Vegetable", "Fruit" = "Fruit", "Fruit/Grain" = "Grain",
  "Sugar" = "Sugar", "Meat" = "Meat", "Meat - Fish" = "Meat", "Meat - fish" = "Meat"
)

MAJOR_COLORS <- c(
  "Dairy" = "#8FA8E0", "Legume / Soy" = "#E0A3A3", "Grain" = "#D3C285",
  "Vegetable" = "#A3D6A0", "Fruit" = "#F5DA97", "Sugar" = "#D9C6AC", "Meat" = "#B79A8B"
)

MAJOR_ORDER <- c("Vegetable", "Grain", "Legume / Soy", "Dairy", "Meat", "Sugar", "Fruit")

major_lookup <- tibble(substrate = names(MAJOR_MAP), major = unname(MAJOR_MAP))

# reverse phase metabolomics
reverse_phase_data <- read_csv("raw_data/metabolomics/2026_03_17_reverse_phase_lc_ms_metabolomics.csv") %>% 
  select(-1) %>% 
  filter(!is.na(lab_id))

# hilic metabolomics
hilic_sample_metadata <- read_csv("metadata/metabolomics/PTFI_hilic_sample_metadata.csv", col_names = FALSE)

META_FIELDS <- c("Global Unique Sample ID", "Study Sample ID", "Sample Name",
                 "Substrate Category", "Sample Type", "PTFI Batch")
hilic_sample_metadata_wide <- hilic_sample_metadata %>%
  rename(field = 1) %>%
  pivot_longer(-field, names_to = "col_pos", values_to = "value") %>%
  pivot_wider(names_from = field, values_from = value) %>%
  select(-col_pos)  

hilic_data <- read_csv("raw_data/metabolomics/2026_06_12_hilic_lc_ms_polar_metabolomics_modf_cols.csv") %>% 
  filter(!is.na(.[[1]]))
