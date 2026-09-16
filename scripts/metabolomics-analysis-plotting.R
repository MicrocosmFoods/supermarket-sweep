library(tidyverse)
library(BiocManager)
library(ComplexHeatmap)
library(colorRamp2)
library(grid)
library(cowplot)

#################################
# Plotting metabolomics top detected metabolites heatmaps
#################################

#################################
# Data import and prep
#################################

## sample metadata
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

## reverse phase metabolomics
reverse_phase_data <- read_csv("raw_data/metabolomics/2026_03_17_reverse_phase_lc_ms_metabolomics.csv") %>% 
  select(-1) %>% 
  filter(!is.na(lab_id)) %>% 
  pivot_longer(
    cols = -(1:5), 
    names_to = "metabolite_name",
    values_to = "metabolite_intensity"
  ) %>% 
  rename(sample_id = study_sample_identifier, substrate = substrate_category) %>% 
  mutate(substrate = str_to_title(substrate)) %>% 
  left_join(major_lookup) %>% 
  mutate(
    major = replace_na(major, "Other"),
    color = MAJOR_COLORS[major],
    color = replace_na(color, "#AAAAAA"),
    label = if_else(is.na(sample_name), sample_id, sample_name)
  )

## hilic metabolomics
# metadata
hilic_sample_metadata <- read_csv("metadata/metabolomics/PTFI_hilic_sample_metadata.csv", col_names = FALSE)

META_FIELDS <- c("Global Unique Sample ID", "Study Sample ID", "Sample Name",
                 "Substrate Category", "Sample Type", "PTFI Batch")

hilic_sample_metadata_wide <- hilic_sample_metadata %>%
  rename(field = 1) %>%
  pivot_longer(-field, names_to = "col_pos", values_to = "value") %>%
  pivot_wider(names_from = field, values_from = value) %>%
  select(-col_pos)

colnames(hilic_sample_metadata_wide) <- c("global_unique_sample_id", "sample_id", "sample_name", "substrate", "sample_type", "ptfi_batch")

# data
hilic_data <- read_csv("raw_data/metabolomics/2026_06_12_hilic_lc_ms_polar_metabolomics_modf_cols.csv") %>% 
  filter(!is.na(.[[1]])) %>% 
  pivot_longer(
    cols = -(1:17),
    names_to = "sample_id",
    values_to = "metabolite_intensity"
  ) %>% 
  rename(metabolite_name = Annotation) %>% 
  left_join(hilic_sample_metadata_wide) %>% 
  left_join(major_lookup) %>% 
  mutate(
    major = replace_na(major, "Other"),
    color = MAJOR_COLORS[major],
    color = replace_na(color, "#AAAAAA"),
    label = if_else(is.na(sample_name), sample_id, sample_name)
  )


#################################
# normalization and selecting top metabolites to plot
# functions for normalization and selecting top metabolites
# function for building heatmap
#################################
TOP_METABOLITES <- 50
MIN_DETECTION <- 0.10
FREQ_DETECTION_CUTOFF <- 0.60

# normalization - filter by detection, then % relative abundance and log2
normalize_metabolites <- function(df) {
  df <- df %>% 
    mutate(metabolite_intensity = replace_na(metabolite_intensity, 0)) %>% 
    group_by(sample_id) %>% 
    mutate(
      pct_abundance = metabolite_intensity / sum(metabolite_intensity) * 100,
      log2_abundance = log2(pct_abundance + 1)
    ) %>% 
    ungroup()
} 

# selection - mean abundance or by frequency

select_top_by_abundance <- function(df_norm, n=TOP_METABOLITES) {
  df_norm %>% 
    group_by(metabolite_name) %>% 
    summarise(mean_pct = mean(pct_abundance), .groups = "drop") %>% 
    slice_max(mean_pct, n=n, with_ties = FALSE) %>% 
    pull(metabolite_name)
} 

build_metabolomics_heatmap <- function(df_norm, top_metabolites, title) {
  heat_data <- df_norm %>%
    filter(metabolite_name %in% top_metabolites)
  
  heat_mat <- heat_data %>%
    select(sample_id, metabolite_name, log2_abundance) %>%
    pivot_wider(names_from = sample_id, values_from = log2_abundance) %>%
    column_to_rownames("metabolite_name") %>%
    as.matrix()
  
  sample_meta <- df_norm %>% distinct(sample_id, major, label)
  
  dup_ids <- sample_meta %>% count(sample_id) %>% filter(n > 1)
  if (nrow(dup_ids) > 0) {
    warning(sprintf(
      "sample_meta has %d duplicated sample_id(s) — keeping the first row for each.\n%s",
      nrow(dup_ids), paste(dup_ids$sample_id, collapse = ", ")
    ))
    sample_meta <- sample_meta %>% distinct(sample_id, .keep_all = TRUE)
  }
  
  col_meta <- sample_meta[match(colnames(heat_mat), sample_meta$sample_id), ]
  
  missing_ids <- colnames(heat_mat)[is.na(col_meta$sample_id)]
  if (length(missing_ids) > 0) {
    warning(sprintf(
      "%d sample_id(s) in heat_mat have no match in sample_meta: %s",
      length(missing_ids), paste(missing_ids, collapse = ", ")
    ))
  }
  
  top_ann <- HeatmapAnnotation(
    Substrate = factor(col_meta$major, levels = MAJOR_ORDER),
    col = list(Substrate = MAJOR_COLORS),
    annotation_name_side = "right",
    show_annotation_name = FALSE,
    annotation_legend_param = list(
      Substrate = list(
        title_gp = gpar(fontsize = 7.5, fontface = "bold"),
        labels_gp = gpar(fontsize = 6.5),
        grid_height = unit(6, "mm"),
        grid_width = unit(6, "mm"),
        direction = "vertical"
      )
    )
  )
  
  col_fun <- colorRamp2(c(0, max(heat_mat, na.rm = TRUE)), c("#F2F2F2", "#1B4F8A"))
  
  Heatmap(
    heat_mat,
    name = "log2(% rel.\nabundance + 1)",
    col = col_fun,
    top_annotation = top_ann,
    column_labels = col_meta$label,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    clustering_distance_rows = "pearson",
    clustering_distance_columns = "euclidean",
    clustering_method_rows = "average",
    clustering_method_columns = "average",
    row_names_gp = gpar(fontsize = 6.5, fontface = "italic"),
    row_names_side = "left",
    column_names_gp = gpar(fontsize = 5.5),
    column_names_rot = 90,
    row_names_max_width = max_text_width(rownames(heat_mat), gp = gpar(fontsize = 6.5, fontface = "italic")) + unit(2, "mm"),
    column_names_max_height = unit(6, "cm"),
    column_title = title,
    heatmap_legend_param = list(
      title = "log2(% rel.\nabundance + 1)",
      title_gp = gpar(fontsize = 7.5, fontface = "bold"),
      labels_gp = gpar(fontsize = 6.5),
      legend_height = unit(5, "cm"),
      direction = "vertical"
    )
  )
}

save_heatmap <- function(ht, filename, width = 12, height = 10, res = 300) {
  pdf(paste0(filename, ".pdf"), width = width, height = height)
  draw(ht)
  dev.off()
  
  png(paste0(filename, ".png"), width = width, height = height, units = "in", res = res)
  draw(ht)
  dev.off()
}



#################################
# normalized and subset metabolomics dfs
#################################
# reverse phase data
reverse_phase_normalized <- normalize_metabolites(reverse_phase_data)
reverse_phase_top_abundant <- select_top_by_abundance(reverse_phase_normalized, 50)

# hilic data - first filter for high confidence data where confidence score = 2
hilic_data_high_conf <- hilic_data %>% 
  filter(ms_confidence == 2)
hilic_normalized <- normalize_metabolites(hilic_data_high_conf)
hilic_top_abundant <- select_top_by_abundance(hilic_normalized, 50)


#################################
# build heatmap data and plots for both datasets
#################################
reverse_phase_heatmap <- build_metabolomics_heatmap(
  reverse_phase_normalized,
  select_top_by_abundance(reverse_phase_normalized, 50),
  sprintf("Top %d reverse-phase metabolites by mean relative abundance", TOP_METABOLITES)
)

hilic_heatmap <- build_metabolomics_heatmap(
  hilic_normalized,
  select_top_by_abundance(hilic_normalized, 50),
  sprintf("Top %d HILIC metabolites by mean relative abundance", TOP_METABOLITES)
)

save_heatmap(reverse_phase_heatmap, "figures/reverse_phase_metabolomics_top_metabolites_heatmap")
save_heatmap(hilic_heatmap, "figures/hilic_metabolomics_top_metabolites_heatmap")
