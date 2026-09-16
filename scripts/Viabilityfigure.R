library(readxl)
library(readr)
library(dplyr)
library(ggplot2)

# -----------------------------------------------------------------------------
# load data
# -----------------------------------------------------------------------------
main_df <- read_excel("raw_data/bioactivity/Fermented_food_data_mastersheet.xlsx", sheet = "full_data")
metadata_df  <- read_tsv("metadata/FF samples - Combined master sheet.tsv")


# -----------------------------------------------------------------------------
# plot constants - category orders and colors
# -----------------------------------------------------------------------------
# Category order placing Kill Control on the reference controls side
category_order <- c("Beet", "Soy", "Cabbage", "Other Vegetable", "Grain", 
                    "Fruit", "Dairy", "Meat", "Sugar", "Kill Control", 
                    "LPS-", "Ruxolitinib", "PDTC", "LPS+")


category_colors <- c(
  "Beet"            = "#E1E5AE",
  "Soy"             = "#BEC559",
  "Cabbage"         = "#3C4228",
  "Other Vegetable" = "#FFAC4D",
  "Grain"           = "#FD966C",
  "Fruit"           = "#DDB9F9",
  "Dairy"           = "#7394E9",
  "Meat"            = "#F6F4EA",
  "Sugar"           = "#E7D4C0",
  "LPS-"            = "#D3D3D3",
  "Ruxolitinib"     = "#A59BCE",
  "PDTC"            = "#9FA8DA", 
  "LPS+"            = "#8C8C8C"
)

# -----------------------------------------------------------------------------
# standardize control samples
# -----------------------------------------------------------------------------
main_df <- main_df %>%
  mutate(
    Sample_type = case_when(
      # Negative control conditions
      Sample_type %in% c("negative control", "negative_control") ~ "negative control",
      Sample_type == "blank" & Insult == "LPS-" ~ "negative control",
      
      # Positive control conditions
      Sample_type %in% c("positive control", "positive_control") ~ "positive control",
      Sample_type == "blank" & Insult == "LPS+" ~ "positive control",
      
      # Kill control variations
      Sample_type %in% c("kill control", "kill_control", "kill _control") ~ "kill control",
      
      TRUE ~ Sample_type
    )
  )

# -----------------------------------------------------------------------------
# subset for controls
# -----------------------------------------------------------------------------
control_types <- c("inhibitor", "positive control", "negative control", "kill control")

subset_df <- main_df %>%
  filter(
    Sample_number %in% metadata_df$`Sample Number` | 
      Sample_type %in% control_types
  )

# -----------------------------------------------------------------------------
# join with metadata and map categories, unfiltered viability
# -----------------------------------------------------------------------------
unfiltered_plot_df <- subset_df %>%
  left_join(metadata_df, by = c("Sample_number" = "Sample Number")) %>%
  mutate(
    # Assign standard group names including Kill Control
    Raw_Group = case_when(
      Sample_type == "negative control" ~ "LPS-",
      Sample_type == "kill control" ~ "Kill Control",
      Sample_type == "inhibitor" & grepl("ruxolitinib", `Pre-treatment`, ignore.case = TRUE) ~ "Ruxolitinib",
      Sample_type == "inhibitor" & grepl("pdtc", `Pre-treatment`, ignore.case = TRUE) ~ "PDTC",
      Sample_type == "positive control" ~ "LPS+",
      Sample_type %in% c("sample", "Sample") ~ coalesce(Substrate, `Substrate category`),
      TRUE ~ NA_character_
    ),
    # Standardize food category names
    Group = case_when(
      Raw_Group %in% c("Beet", "beet") ~ "Beet",
      Raw_Group %in% c("Soy", "soy") ~ "Soy",
      Raw_Group %in% c("Cabbage", "Napa cabbage", "Vegetable - Cabbage") ~ "Cabbage",
      Raw_Group %in% c("Vegetable", "vegetable", "Cucumber", "Carrot", "Other Vegetable") ~ "Other Vegetable",
      Raw_Group %in% c("Grain", "grain") ~ "Grain",
      Raw_Group %in% c("Fruit", "fruit") ~ "Fruit",
      Raw_Group %in% c("Dairy", "dairy", "Milk", "milk") ~ "Dairy",
      Raw_Group %in% c("Meat", "Meat - Fish") ~ "Meat",
      Raw_Group %in% c("Sugar", "sugar") ~ "Sugar",
      TRUE ~ Raw_Group
    )
  ) %>%
  filter(!is.na(Group), !is.na(normalized_viability))

# -----------------------------------------------------------------------------
# factor ordering and groups numbering
# -----------------------------------------------------------------------------
unfiltered_group_counts <- unfiltered_plot_df %>%
  group_by(Group) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(Label = paste0(Group, "\n(n=", n, ")"))

unfiltered_plot_df <- unfiltered_plot_df %>%
  left_join(group_counts, by = "Group") %>%
  mutate(
    Group = factor(Group, levels = category_order),
    Label = factor(Label, levels = unfiltered_group_counts$Label[match(category_order, unfiltered_group_counts$Group)])
  ) %>%
  filter(!is.na(Group))

# -----------------------------------------------------------------------------
# join with metadata and map categories, unfiltered viability
# -----------------------------------------------------------------------------
filtered_plot_df <- subset_df %>%
  left_join(metadata_df, by = c("Sample_number" = "Sample Number")) %>%
  mutate(
    # Assign standard group names including Kill Control
    Raw_Group = case_when(
      Sample_type == "negative control" ~ "LPS-",
      Sample_type == "kill control" ~ "Kill Control",
      Sample_type == "inhibitor" & grepl("ruxolitinib", `Pre-treatment`, ignore.case = TRUE) ~ "Ruxolitinib",
      Sample_type == "inhibitor" & grepl("pdtc", `Pre-treatment`, ignore.case = TRUE) ~ "PDTC",
      Sample_type == "positive control" ~ "LPS+",
      Sample_type %in% c("sample", "Sample") ~ coalesce(Substrate, `Substrate category`),
      TRUE ~ NA_character_
    ),
    # Standardize food category names
    Group = case_when(
      Raw_Group %in% c("Beet", "beet") ~ "Beet",
      Raw_Group %in% c("Soy", "soy") ~ "Soy",
      Raw_Group %in% c("Cabbage", "Napa cabbage", "Vegetable - Cabbage") ~ "Cabbage",
      Raw_Group %in% c("Vegetable", "vegetable", "Cucumber", "Carrot", "Other Vegetable") ~ "Other Vegetable",
      Raw_Group %in% c("Grain", "grain") ~ "Grain",
      Raw_Group %in% c("Fruit", "fruit") ~ "Fruit",
      Raw_Group %in% c("Dairy", "dairy", "Milk", "milk") ~ "Dairy",
      Raw_Group %in% c("Meat", "Meat - Fish") ~ "Meat",
      Raw_Group %in% c("Sugar", "sugar") ~ "Sugar",
      TRUE ~ Raw_Group
    )
  ) %>%
  filter(!is.na(Group), !is.na(normalized_viability)) %>% 
  filter(Sample_type == "kill control" | normalized_viability > 70)

# -----------------------------------------------------------------------------
# factor ordering and groups numbering
# -----------------------------------------------------------------------------
filtered_group_counts <- filtered_plot_df %>%
  group_by(Group) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(Label = paste0(Group, "\n(n=", n, ")"))

filtered_plot_df <- filtered_plot_df %>%
  left_join(group_counts, by = "Group") %>%
  mutate(
    Group = factor(Group, levels = category_order),
    Label = factor(Label, levels = filtered_group_counts$Label[match(category_order, filtered_group_counts$Group)])
  ) %>%
  filter(!is.na(Group))


# -----------------------------------------------------------------------------
# violin plots
# -----------------------------------------------------------------------------
unfiltered_viability_violin_plot <- ggplot(unfiltered_plot_df, aes(x = Group, y = normalized_viability, fill = Group, color = Group)) +
  geom_violin(alpha = 0.5, scale = "width", trim = FALSE, linewidth = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", color = "gray20", outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.7) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "gray50", linewidth = 0.6) +
  geom_vline(xintercept = 9.5, linetype = "dotted", color = "gray60", linewidth = 0.6) +
  annotate(
    "text", x = 9.7, y = 135, label = "Reference\ncontrols →", 
    hjust = 0, size = 3, color = "gray40", fontface = "italic"
  ) +
  scale_fill_manual(values = category_colors, na.translate = FALSE) +
  scale_color_manual(values = category_colors, na.translate = FALSE) +
  scale_x_discrete(
    labels = levels(unfiltered_plot_df$Label), 
    drop = TRUE, 
    na.translate = FALSE
  ) +
  scale_y_continuous(breaks = seq(0, 140, 20), limits = c(-10, 145)) +
  labs(
    x = NULL,
    y = "Viability (%)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "#E5E5E5", linewidth = 0.5),
    axis.text.x = element_text(face = "bold", color = "gray20", size = 9),
    axis.text.y = element_text(color = "gray30"),
    axis.title.y = element_text(face = "bold", size = 11),
    legend.position = "none"
  )

filtered_viability_violin_plot <- ggplot(filtered_plot_df, aes(x = Group, y = normalized_viability, fill = Group, color = Group)) +
  geom_violin(alpha = 0.5, scale = "width", trim = FALSE, linewidth = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", color = "gray20", outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.7) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "gray50", linewidth = 0.6) +
  geom_vline(xintercept = 9.5, linetype = "dotted", color = "gray60", linewidth = 0.6) +
  annotate(
    "text", x = 9.7, y = 135, label = "Reference\ncontrols →", 
    hjust = 0, size = 3, color = "gray40", fontface = "italic"
  ) +
  scale_fill_manual(values = category_colors, na.translate = FALSE) +
  scale_color_manual(values = category_colors, na.translate = FALSE) +
  scale_x_discrete(
    labels = levels(filtered_plot_df$Label), 
    drop = TRUE, 
    na.translate = FALSE
  ) +
  scale_y_continuous(breaks = seq(0, 140, 20), limits = c(-10, 145)) +
  labs(
    x = NULL,
    y = "Viability (%)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "#E5E5E5", linewidth = 0.5),
    axis.text.x = element_text(face = "bold", color = "gray20", size = 9),
    axis.text.y = element_text(color = "gray30"),
    axis.title.y = element_text(face = "bold", size = 11),
    legend.position = "none"
  )


# -----------------------------------------------------------------------------
# export plots
# -----------------------------------------------------------------------------

ggsave("figures/Viability_norm_unfiltered.png", unfiltered_viability_violin_plot, width=30, height=8, units=c("cm"))
ggsave("figures/Viability_norm_filtered.png", filtered_viability_violin_plot, width=30, height=8, units=c("cm"))
