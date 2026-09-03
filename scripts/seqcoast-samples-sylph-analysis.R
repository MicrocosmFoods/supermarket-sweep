library(tidyverse)
library(hues)
library(BiocManager)
library(ComplexHeatmap)
library(colorRamp2)

#################################
# Prep metadata and sylph profiles files
#################################

# MAG metadata from curation repo
mag_metadata_url <- "https://raw.githubusercontent.com/MicrocosmFoods/fermentedfood_metadata_curation/refs/heads/main/data/2025-05-21-genome-metadata-food-taxonomy.tsv"

mag_metadata <- read_tsv(mag_metadata_url) %>% 
  mutate(genome_accession = mag_id) %>% 
  select(genome_accession, completeness, contamination, contigs, taxonomy, species, rep_95id, food_name, main_ingredient, ingredient_group, origin, food_type)

rep_mags_metadata <- mag_metadata %>% 
  filter(genome_accession == rep_95id) %>% 
  select(-rep_95id) %>% 
  mutate(species = case_when(
    is.na(species) | str_to_lower(species) == "unknown" ~ str_c(
      str_extract(taxonomy, "[^;]+$"),
      " spp."
    ),
    TRUE ~ species
  )) %>% 
  select(genome_accession, completeness, contamination, contigs, taxonomy, species)

# sylph profiling results
sylph_profiles_preliminary_run <- read_tsv("results/combined_sylph_profiles_preliminary_run.tsv") %>%
  mutate(accession_name = gsub("_trimmed_1.fastq.gz", "", Sample_file)) %>% 
  mutate(genome_accession = gsub(".fa", "", Genome_file)) %>% 
  select(accession_name, genome_accession, Sequence_abundance, Adjusted_ANI, Eff_cov, Contig_name)

sylph_profiles_full_run <- read_tsv("results/combined_sylph_profiles_full_run.tsv") %>% 
  mutate(accession_name = gsub("_trimmed_1.fastq.gz", "", Sample_file)) %>% 
  mutate(genome_accession = gsub(".fa", "", Genome_file)) %>% 
  select(accession_name, genome_accession, Sequence_abundance, Adjusted_ANI, Eff_cov, Contig_name)

all_sylph_profiles <- rbind(sylph_profiles_preliminary_run, sylph_profiles_full_run)

# sample metadata
sample_metadata_preliminary_run <- read.csv("metadata/2025-12-23-seqcoast-preliminary-run-samples-metadata.csv") %>% 
  mutate(accession_name = gsub("_R1.fastq.gz", "", fastq_1)) %>% 
  select(accession_name, sample_name, fermented_food)

sample_metadata_full_run <- read.csv("metadata/2026-01-26-seqcoast-full-run-sample-metadata.csv") %>% 
  mutate(accession_name = gsub("_R1.fastq.gz", "", fastq_1)) %>% 
  select(accession_name, sample_name, fermented_food)

all_metadata <- rbind(sample_metadata_preliminary_run, sample_metadata_full_run)


# merge with genome and sample metadata
sylph_profiles_metadata <- left_join(all_sylph_profiles, rep_mags_metadata) %>% 
  left_join(all_metadata) %>% 
  mutate(genus = str_extract(taxonomy, "[^;]+$"))

write_tsv(sylph_profiles_metadata, "results/2026-01-30-seqcoast-supermarket-sweep-profiles.tsv")

#################################
# Basic summary stats
#################################

# summary stats per sample
sylph_profiles_stats <- sylph_profiles_metadata %>% 
  group_by(sample_name) %>%
  summarise(
    n_genomes = n_distinct(genome_accession),
    percent_mapped = round(sum(Sequence_abundance, na.rm = TRUE), 3),
    percent_unmapped = round(100 - sum(Sequence_abundance, na.rm = TRUE), 3),
    .groups = "drop"
  )

sylph_profile_stats_5p_cutoff <- sylph_profiles_metadata %>% 
  filter(Eff_cov > 0.05) %>% 
  group_by(sample_name) %>%
  summarise(
    n_genomes = n_distinct(genome_accession),
    percent_mapped = round(sum(Sequence_abundance, na.rm = TRUE), 3),
    percent_unmapped = round(100 - sum(Sequence_abundance, na.rm = TRUE), 3),
    .groups = "drop"
  )

#################################
# Stacked bar plot of abundance of top species in samples
#################################
# genus map
genus_map <- tibble::tribble(
  ~pattern,                         ~genus_group,
  "^Bacillus(_.*)?$",               "Bacillus",
  "^Enterococcus(_.*)?$",           "Enterococcus")

sylph_profiles_metadata <- sylph_profiles_metadata %>%
  mutate(
    genus_group = purrr::map_chr(
      genus,
      \(g) {
        hit <- genus_map %>% filter(str_detect(g, pattern))
        if (nrow(hit) > 0) hit$genus_group[[1]] else g
      }
    )
  )


# prep df for showing abundance of top species
abundance_df_labelled <- sylph_profiles_metadata %>% 
  group_by(sample_name) %>% 
  arrange(desc(Sequence_abundance), .by_group = TRUE) %>% 
  mutate(
    rank_in_sample = row_number(),
    genus_label = if_else(rank_in_sample <=4, genus_group, "Other Species")
  ) %>% 
  ungroup() %>% 
  select(sample_name, Sequence_abundance, genus_label)

# leave out other species for specific grey color
other_genera <- "Other Genera"

# Get the levels actually present in the plot
genus_levels <- abundance_df_labelled %>%
  distinct(genus_label) %>%
  pull(genus_label)

main_genera <- setdiff(genus_levels, other_genera)

# palette prep
main_colors <- khroma::colour("smoothrainbow")(length(main_genera))
names(main_colors) <- main_genera

# Add grey for Other Species
fill_colors <- c(main_colors, "Other Genera" = "grey70")

# different palette strategy
genus_colors <- setNames(
  iwanthue(65),
  unique(abundance_df_labelled$genus_label)
)

# plot
supermarket_sweep_samples_abundance_plot <- abundance_df_labelled %>%
  mutate(sample_name = gsub("_", " ", sample_name)) %>% 
  ggplot(aes(x = sample_name, y = Sequence_abundance, fill = genus_label)) +
  geom_col() +
  theme_bw() +
  scale_x_discrete(expand = c(0, 0),
                   labels = function(x) stringr::str_wrap(x, width = 6)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = fill_colors, drop = FALSE) +
  guides(fill = guide_legend(ncol = 2)) +
  theme(
    axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    plot.title = element_text(face = "bold", size = 16),

    legend.text = element_text(size = 14),
    legend.title = element_text(size = 15),
    legend.key.height = unit(0.75, "cm"),
    legend.spacing.y = unit(0.1, "cm"),
    legend.margin = margin(6, 6, 6, 6),

    strip.text = element_text(size = 15),
    strip.background = element_rect(fill = "white", color = "black")
  ) +
  labs(
    x = "Sample Name",
    y = "% Sequence Abundance of Genera",
    fill = "Genus",
    title = "Top Abundant Genera in Supermarket Sweep Preliminary Samples"
  )

supermarket_sweep_samples_abundance_plot

# save supermarket sweep abundance plot
ggsave("figures/supermarket-sweep-abundance-plot.png", supermarket_sweep_samples_abundance_plot, width=45, height=30, units=c("cm"), dpi=300)
 
#################################
# Binary conversion
# For species detected with at least 5X effective coverage, count as "present" in the sample
# Join with the full list of species in the representative database to convert to a binary table of 0s and 1s
#################################

# filter by 5X coverage and create the "detected" column 
profiles_0.05_covg_filtered_binary <- sylph_profiles_metadata %>% 
  filter(Eff_cov >= 0.05) %>% 
  select(sample_name, genome_accession) %>% 
  mutate(sample_name = gsub(" ", "_", sample_name)) %>% 
  distinct() %>% 
  mutate(detected = 1L) %>% 
  pivot_wider(
    names_from  = sample_name,
    values_from = detected,
    values_fill = list(detected = 0L)
  )

# write out presence/absence TSV
write_tsv(profiles_0.05_covg_filtered_binary, "results/all-profiles-0.05x-covg-filtered-binary-matrix.tsv")

#################################
# Plot heatmap of top most abundant species
# top 50 species at least 1% abundance in at least 1 sample
# log-transformed
#################################

# Filter to top 1% abundance in at least 1 sample
# Merge with metadata

supermarket_sample_metadata <- read_tsv("metadata/Combined master sheet updated 042126-8ed4f7f4.tsv") %>% 
  mutate(sample_number = `Sample Number`)

# substrate palette and groupings
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

# clean column names and substrate categories, colors
sylph_profiles_metadata_clean <- sylph_profiles_metadata %>% 
  mutate(sample_name = gsub(" ", "_", sample_name)) %>% 
  separate_wider_delim(cols=sample_name, delim="_", names=c("sample_number", "sample_name"), too_many="merge") %>% 
  left_join(supermarket_sample_metadata) %>% 
  select(sample_number, Sample, `Substrate category`, `Food type`, taxonomy, species, Sequence_abundance, Adjusted_ANI, Eff_cov) %>% 
  rename(substrate = `Substrate category`, food_type = `Food type`, sample = Sample, abundance = Sequence_abundance, ANI = Adjusted_ANI, covg = Eff_cov) %>% 
  left_join(major_lookup) %>% 
  mutate(
    major = replace_na(major, "Other"),
    color = MAJOR_COLORS[major],
    color = replace_na(color, "#AAAAAA"),
    label = if_else(is.na(sample), sample_number, sample)
  )

food_meta <- sylph_profiles_metadata_clean %>% 
  select(sample_number, major, color, label)

# filtered abundance table
abund_long <- sylph_profiles_metadata_clean %>%
  group_by(sample_number, species) %>%
  summarise(abundance = sum(abundance), .groups = "drop") %>%
  complete(sample_number, species, fill = list(abundance = 0)) %>%
  group_by(sample_number) %>%
  mutate(pct_abund = 100 * abundance / sum(abundance)) %>%
  ungroup()

# keep 50 most abundant species by mean relative abundance
N_TOP_SPECIES = 50

keep_species <- abund_long %>%
  group_by(species) %>%
  summarise(max_pct = max(pct_abund), .groups = "drop") %>%
  filter(max_pct > 1) %>%
  pull(species)

top_species <- abund_long %>%
  filter(species %in% keep_species) %>%
  group_by(species) %>%
  summarise(mean_pct = mean(pct_abund), .groups = "drop") %>%
  slice_max(mean_pct, n = N_TOP_SPECIES) %>%
  pull(species)

# log2 transform
heat_data <- abund_long %>%
  filter(species %in% top_species) %>%
  mutate(log2_pct = log2(pct_abund + 1))


# pivot to matrix for ComplexHeatmap
heat_mat <- heat_data %>% 
  select(sample_number, species, log2_pct) %>% 
  pivot_wider(names_from =sample_number, values_from = log2_pct) %>% 
  column_to_rownames("species") %>% 
  as.matrix()

# ─── Column (food) annotation — this is the "grid on top" ─────────────────
# Check food_meta is one row per sample_id before joining — a duplicate
# Sample Number in supermarket_sample_metadata is the usual cause of the
# "Length of column_labels should be the same as ncol of matrix" error,
# since left_join() will silently expand rows for any duplicated key.
dup_ids <- food_meta %>% count(sample_number) %>% filter(n > 1)
if (nrow(dup_ids) > 0) {
  warning(sprintf(
    "food_meta has %d duplicated sample_id(s) — keeping the first row for each.\n%s",
    nrow(dup_ids), paste(dup_ids$sample_number, collapse = ", ")
  ))
  food_meta <- food_meta %>% distinct(sample_number, .keep_all = TRUE)
}

# Index (not join) so col_meta is guaranteed the same length AND same order
# as colnames(z_mat), regardless of how food_meta is sorted.
col_meta <- food_meta[match(colnames(heat_mat), food_meta$sample_number), ]

# Any sample_id present in z_mat but missing from food_meta will show up as
# NA here — check for that too, since it will make Substrate/label blank.
missing_ids <- colnames(heat_mat)[is.na(col_meta$sample_number)]
if (length(missing_ids) > 0) {
  warning(sprintf(
    "%d sample_id(s) in z_mat have no match in food_meta: %s",
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

# color scale

col_fun <- colorRamp2(
  c(0, max(heat_mat)),
  c("#F2F2F2", "#1B4F8A")
)

# heatmap plot
ht <- Heatmap(
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
  row_names_max_width = max_text_width(rownames(z_mat), gp = gpar(fontsize = 6.5, fontface = "italic")) + unit(2, "mm"),
  column_names_max_height = unit(6, "cm"),   # bottom x-axis label space — just increase this number if labels are cut off
  column_title = sprintf("Top %d most abundant species (>1%% in >=1 sample) across fermented foods", N_TOP_SPECIES),
  heatmap_legend_param = list(
    title = "log2(% rel.\nabundance + 1)",
    title_gp = gpar(fontsize = 7.5, fontface = "bold"),
    labels_gp = gpar(fontsize = 6.5),
    legend_height = unit(5, "cm"),
    direction = "vertical"
  )
)

ht

# export plot
OUT <- "figures"

# width and height
FIG_WIDTH  <- 15   # inches
FIG_HEIGHT <- 8   # inches

ht_opt$legend_gap <- unit(4, "mm")

pdf(file.path(OUT, "Metagenomics_Full_Abundance_Heatmap.pdf"), width = FIG_WIDTH, height = FIG_HEIGHT)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right", merge_legend = TRUE,
     align_heatmap_legend = "heatmap_top", align_annotation_legend = "heatmap_top")
dev.off()

png(file.path(OUT, "Metagenomics_Full_Abundance_Heatmap.png"),
    width = FIG_WIDTH, height = FIG_HEIGHT, units = "in", res = 200)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right", merge_legend = TRUE,
     align_heatmap_legend = "heatmap_top", align_annotation_legend = "heatmap_top")
dev.off()

