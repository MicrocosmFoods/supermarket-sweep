library(readxl)
library(readr)
library(dplyr)
library(ggplot2)

# -----------------------------------------------------------------------------
# 1. Load Data
# -----------------------------------------------------------------------------
main_df <- read_excel("Fermented_food_data_mastersheet.xlsx", sheet = "full_data")
tsv_df  <- read_tsv("FF samples - Combined master sheet.tsv")

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
# Set your desired threshold here:
VIABILITY_THRESHOLD <- 70 

plot_df <- subset_df %>%
  left_join(tsv_df, by = c("Sample_number" = "Sample Number")) %>%
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
# 5. Factor Ordering & Group Counts (n=)
# -----------------------------------------------------------------------------
# Category order placing Kill Control on the reference controls side
category_order <- c("Beet", "Soy", "Cabbage", "Other Vegetable", "Grain", 
                    "Fruit", "Dairy", "Meat", "Sugar", "Kill Control", 
                    "LPS-", "Ruxolitinib", "PDTC", "LPS+")

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
  "Beet"            = "#D88BB3",
  "Soy"             = "#E8A3A3",
  "Cabbage"         = "#A2D89B",
  "Other Vegetable" = "#A0CBE8",
  "Grain"           = "#E3D5A5",
  "Fruit"           = "#F8D489",
  "Dairy"           = "#A3BFE8",
  "Meat"            = "#C49A8D",
  "Sugar"           = "#E7D4C0",
  "Kill Control"    = "#E57373",
  "LPS-"            = "#D3D3D3",
  "Ruxolitinib"     = "#B39DDB", # Soft purple
  "PDTC"            = "#9FA8DA", # Slate purple/blue
  "LPS+"            = "#8C8C8C"
)

# -----------------------------------------------------------------------------
# 7. Viability Violin Plot Generation
# -----------------------------------------------------------------------------
ggplot(plot_df, aes(x = Group, y = normalized_viability, fill = Group, color = Group)) +
  # Violin outline & fill
  geom_violin(alpha = 0.5, scale = "width", trim = FALSE, linewidth = 0.8) +
  
  # Inner boxplot for Median & IQR
  geom_boxplot(width = 0.1, fill = "white", color = "gray20", outlier.shape = NA, alpha = 0.8) +
  
  # Individual points with jitter
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.7) +
  
  # Baseline reference line at 100%
  geom_hline(yintercept = 100, linetype = "dashed", color = "gray50", linewidth = 0.6) +
  
  # Separator line separating sample groups from controls (x = 9.5)
  geom_vline(xintercept = 9.5, linetype = "dotted", color = "gray60", linewidth = 0.6) +
  annotate(
    "text", x = 9.7, y = 135, label = "Reference\ncontrols →", 
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
