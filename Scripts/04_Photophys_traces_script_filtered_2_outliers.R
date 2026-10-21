rm(list = ls())
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(here)


load(here("Outputs", "photophys_preprocessed_data2.RData"))
top_vectors2 <- readRDS(here("Outputs", "top_vectors2.rds"))
source(here("Scripts/00_visualization_prep.R"))

# Define dendrogram sample IDs from fit_cbb
dend_ids <- fit_cbb2$hclust$labels %>%
  tolower() %>%
  str_replace("\\.processed\\.csv$", "") %>%
  str_replace("\\._?\\d{2,3}r$", "") %>%
  trimws()

# 2. Clean MessyID to match standard ID format (e.g., t1_aabr1 -> 1_aabr_1)
longphoto_shortie <- longphoto_clean %>%
  mutate(
    SampleID_clean = MessyID %>%
      tolower() %>%
      str_replace("\\.processed\\.csv$", "") %>%
      str_replace("\\._?\\d{2,3}r$", "") %>%
      str_replace("^t(\\d+)_([a-z]+)(\\d+)$", "\\1_\\2_\\3") %>%
      trimws(),
    # Create uppercase ID version for joining with meta_filteredK
    SampleID_upper = toupper(SampleID_clean)
  ) %>%
  filter(!str_starts(metric, "m")) %>% 
  filter(SampleID_clean %in% dend_ids | SampleID_upper %in% meta_filteredK$SampleID_clean)

# Join metadata (keeps K_3 and Site)
# Process longphoto and join metadata using roman_map for Clusters
# Process longphoto and join metadata using roman_map for Clusters
df_with_clusters <- longphoto_shortie %>%
  left_join(
    meta_filteredK2 %>% dplyr::select(SampleID_clean, Species, Site, K_3), 
    by = c("SampleID_upper" = "SampleID_clean")
  ) %>%
  mutate(
    stage = as.numeric(as.character(stage)),
    Value = as.numeric(as.character(Value)),
    # Apply cluster mapping from visualization prep
    Cluster_Label = factor(roman_map[as.character(K_3)], levels = target_levels)
  )

top_families <- top_vectors2 %>%
  { if ("Root.Metric" %in% names(.)) pull(., Root.Metric) else pull(., Metric) } %>%
  c(if ("Family" %in% names(top_vectors2)) pull(top_vectors2, Family) else NULL) %>%
  sub(" .*", "", .) %>%       # Drop appended description text if present
  sub("\\..*", "", .) %>%     # Strip channel/stage suffixes (e.g., Quant.1.0 -> Quant)
  unique() %>%
  .[!is.na(.) & . != ""]

# Construct clean word-boundary regex pattern
top_pattern <- paste0("\\b(", paste(top_families, collapse = "|"), ")\\b")

# Exclude list for redundant/raw parameters
exclude_metrics <- c("nCon", "Fm", "Fo", "nTau2", "nTau1")

# 2. Filter sub_photo dynamically
sub_photo <- df_with_clusters %>%
  filter(
    !metric %in% exclude_metrics,
    col %in% c("1", "2", "3", "4"),       # 4 excitation channels
    stage >= 1 & stage <= 27,            # Actinic light steps
    !str_starts(metric, "m|MT"),         # Drop m-prefixed/derived metrics
    metric %in% top_vectors2$Metric | grepl(top_pattern, metric, ignore.case = TRUE)
  )
# Summarize mean values per stage, excitation channel, cluster label, and metric
sum_photo_long <- sub_photo %>%
  group_by(stage, col, Cluster_Label, metric) %>%
  summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop") %>%
  filter(!is.na(Cluster_Label))

# Print summary check
cat("Successfully processed rows in sum_photo_long:", nrow(sum_photo_long), "\n")
print(unique(sum_photo_long$metric))

# Calculate mean and standard error per group/step
sum_photo_stats <- sub_photo %>%
  group_by(stage, col, Cluster_Label, metric) %>%
  summarise(
    mean_val = mean(Value, na.rm = TRUE),
    se_val   = sd(Value, na.rm = TRUE) / sqrt(n()),
    .groups  = "drop"
  ) %>%
  filter(!is.na(Cluster_Label))

# 4. Plot mean trends with shaded confidence intervals
ggplot(sum_photo_stats, aes(x = stage, y = mean_val, color = as.factor(col), fill = as.factor(col))) +
  geom_ribbon(
    aes(ymin = mean_val - se_val, ymax = mean_val + se_val, group = col),
    alpha = 0.9, 
    color = NA
  ) +
  geom_line(aes(group = col), linewidth = 0.8) +
  facet_grid(metric ~ Cluster_Label, scales = "free_y") +
  theme_coral_pub() +
  scale_color_manual(
    values = c("1" = "#8A2BE2", "2" = "violet", "3" = "blue", "4" = "cyan"),
    labels = c("1" = "420-nm excitation", "2" = "442-nm excitation", "3" = "505-nm excitation", "4" = "525-nm excitation")
  ) +
  scale_fill_manual(
    values = c("1" = "#8A2BE2", "2" = "violet", "3" = "blue", "4" = "cyan"),
    guide = "none"
  ) +
  labs(x = "Sample Light Step", y = "Mean Trait Value (± SE)")

# 4. Plot Traces across K3 Clusters (avg)
ggplot() + 
  # Individual sample traces (thin, semi-transparent)
  geom_line(
    data = sub_photo,
    aes(
      x = stage, 
      y = Value, 
      group = interaction(SampleID_clean, col), 
      colour = as.factor(col)
    ), 
    linewidth = 0.3,
    alpha = 0.35
  ) + 
  # Group mean traces (thicker, solid overlay)
  geom_line(
    data = sum_photo_long,
    aes(
      x = stage, 
      y = Value, 
      group = col, 
      colour = as.factor(col)
    ), 
    linewidth = 0.9
  ) + 
  facet_grid(metric ~ Cluster_Label, scales = "free_y") +
  theme_coral_pub() +
  scale_color_manual(
    values = c("1" = "#8A2BE2", "2" = "violet", "3" = "blue", "4" = "cyan"),
    labels = c("1" = "408-nm excitation", "2" = "415-nm excitation", "3" = "470-nm excitation", "4" = "505-nm excitation"),
    name   = "Excitation Channel"
  ) +
  theme(
    legend.position  = "bottom",
    strip.text       = element_text(size = 10, face = "bold"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3)
  ) +
  labs(
    x = "Sample Light Step",
    y = "Trait Value"
  )

# 4. Plot Traces across K3 Clusters (all)
# Plot Traces across mapped Cluster Labels using theme_coral_pub
# Plot Traces across mapped Cluster Labels using theme_coral_pub
ggplot(sum_photo_long) + 
  geom_line(
    aes(x = stage, y = Value, group = interaction(col, Cluster_Label), colour = as.factor(col)), 
    linewidth = 0.8
  ) + 
  facet_grid(metric ~ Cluster_Label, scales = "free_y") +
  theme_coral_pub() +
  scale_color_manual(
    values = c("1" = "#8A2BE2", "2" = "violet", "3" = "blue", "4" = "cyan"),
    labels = c("1" = "408-nm excitation", "2" = "415-nm excitation", "3" = "470-nm excitation", "4" = "505-nm excitation"),
    name   = "Excitation Channel")+
  theme(
    legend.position = "bottom",
    strip.text      = element_text(size = 10, face = "bold"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3)
  ) +
  labs(
    x = "Light Step Stage (1–27)",
    y = "Mean Trait Value",
    #title = "Photophysiological Traces across Functional Clusters"
  )

# 5. Association Tests (Species / Site vs K_3)
table_spp <- table(meta_filteredK$Species, meta_filteredK$K_3)
print(table_spp)
chisq.test(table_spp, simulate.p.value = TRUE)

table_site <- table(meta_filteredK$Site, meta_filteredK$K_3)
print(table_site)
chisq.test(table_site, simulate.p.value = TRUE)