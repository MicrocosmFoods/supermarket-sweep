library(tidyverse)
library(hues)
library(ComplexUpset)

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
