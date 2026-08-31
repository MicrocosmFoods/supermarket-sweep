library(readxl)
library(readr)
library(dplyr)
library(ggplot2)

# -----------------------------------------------------------------------------
# 1. Load Data
# -----------------------------------------------------------------------------
main_df <- read_excel("raw_data/bioactivity/Fermented_food_data_mastersheet.xlsx", sheet = "full_data")
tsv_df  <- read_tsv("metadata/FF samples - Combined master sheet.tsv")

# -----------------------------------------------------------------------------
# 2. Update & Standardize Control Sample Types
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
# 3. Create Subset
# -----------------------------------------------------------------------------
control_types <- c("inhibitor", "positive control", "negative control", "kill control")

subset_df <- main_df %>%
  filter(
    Sample_number %in% tsv_df$`Sample Number` | 
      Sample_type %in% control_types
  )

# -----------------------------------------------------------------------------
# 4. Join TSV Data & Map Categories
# -----------------------------------------------------------------------------
plot_df <- subset_df %>%
  left_join(tsv_df, by = c("Sample_number" = "Sample Number")) %>%
  mutate(
    # Assign standard group names
    Raw_Group = case_when(
      Sample_type == "negative control" ~ "LPS-",
      Sample_type == "inhibitor" & `Pre-treatment` == "100nM ruxolitinib" ~ "Ruxolitinib",
      Sample_type == "positive control" ~ "LPS+",
      Sample_type %in% c("sample", "Sample") ~ coalesce(Substrate, `Substrate category`),
      TRUE ~ NA_character_
    ),
    # Standardize category names
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
  filter(!is.na(Group), !is.na(normalized_IFN)) %>%
  # --- Filter out low positive control outlier values (< 50) ---
  filter(!(Group == "LPS+" & normalized_IFN < 50))

# -----------------------------------------------------------------------------
# 5. Factor Ordering & Group Counts (n=)
# -----------------------------------------------------------------------------
category_order <- c("Beet", "Soy", "Cabbage", "Other Vegetable", "Grain", 
                    "Fruit", "Dairy", "Meat", "Sugar", "LPS-", "Ruxolitinib", "LPS+")

group_counts <- plot_df %>%
  group_by(Group) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(Label = paste0(Group, "\n(n=", n, ")"))

plot_df <- plot_df %>%
  left_join(group_counts, by = "Group") %>%
  mutate(
    Group = factor(Group, levels = category_order),
    Label = factor(Label, levels = group_counts$Label[match(category_order, group_counts$Group)])
  ) %>%
  filter(!is.na(Group))

# -----------------------------------------------------------------------------
# 6. Color Palette 
# -----------------------------------------------------------------------------
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
  "LPS+"            = "#8C8C8C"
)

# -----------------------------------------------------------------------------
# 7. Violin Plot Generation
# -----------------------------------------------------------------------------
ifn_violin_plot <- ggplot(plot_df, aes(x = Group, y = normalized_IFN, fill = Group, color = Group)) +
  # Full Violin outline & fill
  geom_violin(alpha = 0.5, scale = "width", trim = FALSE, linewidth = 0.8) +
  
  # Narrow inner boxplot to display Median and IQR
  geom_boxplot(width = 0.1, fill = "white", color = "gray20", outlier.shape = NA, alpha = 0.8) +
  
  # Individual data points overlaid with jitter
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.7) +
  
  # Reference Lines
  geom_hline(yintercept = 100, linetype = "dashed", color = "gray50", linewidth = 0.6) +
  geom_vline(xintercept = 9.5, linetype = "dotted", color = "gray60", linewidth = 0.6) +
  annotate(
    "text", x = 9.7, y = 168, label = "Reference\ncontrols →", 
    hjust = 0, size = 3, color = "gray40", fontface = "italic"
  ) +
  
  # Aesthetics & Scales
  scale_fill_manual(values = category_colors, na.translate = FALSE) +
  scale_color_manual(values = category_colors, na.translate = FALSE) +
  scale_x_discrete(
    labels = levels(plot_df$Label), 
    drop = TRUE, 
    na.translate = FALSE
  ) +
  scale_y_continuous(breaks = seq(0, 175, 25), limits = c(-10, 180)) +
  labs(
    x = NULL,
    y = expression(atop("IFN activity", paste("(LPS+ = 100, LPS- = 0)")))
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "#E5E5E5", linewidth = 0.5),
    axis.text.x = element_text(face = "bold", color = "gray20", size = 10),
    axis.text.y = element_text(color = "gray30"),
    axis.title.y = element_text(face = "bold", size = 11),
    legend.position = "none"
  )

# -----------------------------------------------------------------------------
# 8. Export plot
# -----------------------------------------------------------------------------

ggsave("figures/IFN_norm.png", ifn_violin_plot, width=30, height=8, units=c("cm"))
