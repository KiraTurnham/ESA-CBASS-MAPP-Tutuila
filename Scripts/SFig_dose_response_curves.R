library(dplyr)
library(ggplot2)
library(here)
library(patchwork)
library(dplyr)
library(tidyr)
library(ggplot2)
library(here)
library(patchwork)
library(CBASSED50)

# -----------------------------------------------------------------------------
# 1. LOAD PREPROCESSED DATA & FIT GLOBAL MODELS
# -----------------------------------------------------------------------------
source(here("Scripts/00_visualization_prep.R"))
load(here("Outputs", "photophys_preprocessed_data.RData"))

# Define file paths to your raw CBASS/PAM run files
run_files <- list.files(
  path       = here("Data"), # Adjust to your actual raw data directory
  pattern    = "\\PAM.csv$",              # Adjust file extension (.csv, .xlsx, etc.)
  full.names = TRUE
)

# Combine raw site files
df_combined <- map_dfr(run_files, function(file_path) {
  df_run <- read_data(file_path)
  df_run %>%
    mutate(
      PAM1 = ifelse(is.na(PAM1), Pam_value, PAM1),
      PAM2 = ifelse(is.na(PAM2), Pam_value, PAM2)
    ) %>%
    preprocess_dataset()
})

grouping_properties <- c("Site", "Species")
drm_formula         <- "Pam_value ~ Temperature"

# Fit global models and predict PAM yield profiles
models_global <- fit_drms(df_combined, grouping_properties, drm_formula, is_curveid = FALSE, LL.4 = FALSE)
temp_ranges   <- define_temperature_ranges(df_combined$Temperature, n = 100)
predictions   <- get_predicted_pam_values(models_global, temp_ranges)

# -----------------------------------------------------------------------------
# 2. JOIN METADATA VIA SampleID_clean
# -----------------------------------------------------------------------------
# Format GroupingProperty to match SampleID_clean (e.g., "1_AABR_1")
df_combined_grouped <- define_grouping_property(df_combined, grouping_properties) %>%
  mutate(SampleID_clean = paste(GroupingProperty, Genotype, sep = "_")) %>%
  dplyr::select(GroupingProperty, SampleID_clean, Site, Species) %>%
  distinct()

# Prepare cluster metadata lookup
cluster_lookup <- meta_filteredK %>%
  dplyr::select(SampleID_clean, K_3) %>%
  mutate(
    # Align cluster level labels with cluster_colors if needed
    K_3 = factor(K_3, labels = c("Cluster I", "Cluster II", "Cluster III"))
  ) %>%
  distinct()

# Combine predictions with meta_filteredK
predictions_df <- predictions %>%
  left_join(df_combined_grouped, by = "GroupingProperty") %>%
  left_join(cluster_lookup, by = "SampleID_clean") %>%
  filter(!is.na(Species))

# -----------------------------------------------------------------------------
# 3. GENERATE FIGURE S1 (DRC BY SPECIES VS. K3 CLUSTER)
# -----------------------------------------------------------------------------

# Panel A: DRC by Species
p_drc_spp <- ggplot(predictions_df, aes(x = Temperature, y = PredictedPAM, color = Species, fill = Species)) +
  geom_line(aes(group = SampleID_clean), alpha = 0.35, linewidth = 0.5) +
  geom_smooth(method = "loess", se = TRUE, alpha = 0.20, linewidth = 1.2) +
  scale_color_manual(name = "Species", values = species_colors, limits = species_order, labels = species_labels) +
  scale_fill_manual(name = "Species", values = species_colors, limits = species_order, labels = species_labels) +
  labs(
    x = "Temperature (°C)",
    y = expression("Photophysiological Yield (" * F[v] / F[m] * ")")
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.title   = element_text(face = "bold", size = 11, color = "black"),
    axis.text    = element_text(size = 10, color = "black"),
    legend.title = element_text(face = "bold", size = 11),
    legend.text  = element_text(size = 10, face = "italic"),
    legend.position = "right"
  )

# Panel B: DRC by K3 Functional Cluster
p_drc_k3 <- ggplot(
  predictions_df %>% filter(!is.na(K_3)), 
  aes(x = Temperature, y = PredictedPAM, color = K_3, fill = K_3)
) +
  geom_line(aes(group = SampleID_clean), alpha = 0.30, linewidth = 0.5) +
  geom_smooth(method = "loess", se = TRUE, alpha = 0.20, linewidth = 1.2) +
  scale_color_manual(name = "Cluster", values = cluster_colors) +
  scale_fill_manual(name = "Cluster", values = cluster_colors) +
  labs(
    x = "Temperature (°C)",
    y = expression("Photophysiological Yield (" * F[v] / F[m] * ")")
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.title   = element_text(face = "bold", size = 11, color = "black"),
    axis.text    = element_text(size = 10, color = "black"),
    legend.title = element_text(face = "bold", size = 11),
    legend.text  = element_text(size = 10),
    legend.position = "right"
  )

# Combine panels
fig_s1_drc <- (p_drc_spp / p_drc_k3) + 
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(face = "bold", size = 14))

print(fig_s1_drc)

# Save output
# ggsave(
#   here("Outputs", "Figure_S1_Dose_Response_Curves.pdf"),
#   fig_s1_drc,
#   width = 9,
#   height = 9,
#   device = "pdf"
# )
# -----------------------------------------------------------------------------
# 1. AGGREGATE SUMMARY METRICS (ED5, ED50, ED95) AND PREDICTIONS
# -----------------------------------------------------------------------------

# Calculate ED metric means by Species
eds_spp_summary <- meta_filteredK %>%
  group_by(Species) %>%
  summarise(
    Mean_ED5  = mean(ED5, na.rm = TRUE),
    Mean_ED50 = mean(ED50, na.rm = TRUE),
    Mean_ED95 = mean(ED95, na.rm = TRUE),
    .groups   = "drop"
  )

# Calculate ED metric means by K_3 Cluster
eds_k3_summary <- meta_filteredK %>%
  filter(!is.na(K_3)) %>%
  mutate(K_3 = factor(K_3, labels = c("Cluster I", "Cluster II", "Cluster III"))) %>%
  group_by(K_3) %>%
  summarise(
    Mean_ED5  = mean(ED5, na.rm = TRUE),
    Mean_ED50 = mean(ED50, na.rm = TRUE),
    Mean_ED95 = mean(ED95, na.rm = TRUE),
    .groups   = "drop"
  )

# Aggregate Mean Yield & 95% Confidence Intervals for Species
drc_spp_avg <- predictions_df %>%
  group_by(Species, Temperature) %>%
  summarise(
    Mean_PAM = mean(PredictedPAM, na.rm = TRUE),
    Lower_CI = quantile(PredictedPAM, 0.025, na.rm = TRUE),
    Upper_CI = quantile(PredictedPAM, 0.975, na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  left_join(eds_spp_summary, by = "Species")

# Aggregate Mean Yield & 95% Confidence Intervals for K3 Cluster
drc_k3_avg <- predictions_df %>%
  filter(!is.na(K_3)) %>%
  group_by(K_3, Temperature) %>%
  summarise(
    Mean_PAM = mean(PredictedPAM, na.rm = TRUE),
    Lower_CI = quantile(PredictedPAM, 0.025, na.rm = TRUE),
    Upper_CI = quantile(PredictedPAM, 0.975, na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  left_join(eds_k3_summary, by = "K_3")

# -----------------------------------------------------------------------------
# 2. BUILD AVERAGE DOSE-RESPONSE PLOTS WITH RIBBONS AND ED LINES
# -----------------------------------------------------------------------------

# Panel A: Average DRC by Species
p_avg_spp <- ggplot(drc_spp_avg, aes(x = Temperature, y = Mean_PAM, color = Species, fill = Species)) +
  geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.2) +
  
  # ED5 Line & Label
  geom_vline(data = eds_spp_summary, aes(xintercept = Mean_ED5, color = Species), linetype = "dashed", alpha = 0.7) +
  geom_text(data = eds_spp_summary, aes(x = Mean_ED5, y = 0.1, label = paste0("ED5: ", round(Mean_ED5, 1))), 
            angle = 90, vjust = -0.4, size = 3, show.legend = FALSE) +
  
  # ED50 Line & Label
  geom_vline(data = eds_spp_summary, aes(xintercept = Mean_ED50, color = Species), linetype = "dashed", alpha = 0.9) +
  geom_text(data = eds_spp_summary, aes(x = Mean_ED50, y = 0.35, label = paste0("ED50: ", round(Mean_ED50, 1))), 
            angle = 90, vjust = -0.4, size = 3, fontface = "bold", show.legend = FALSE) +
  
  # ED95 Line & Label
  geom_vline(data = eds_spp_summary, aes(xintercept = Mean_ED95, color = Species), linetype = "dashed", alpha = 0.7) +
  geom_text(data = eds_spp_summary, aes(x = Mean_ED95, y = 0.1, label = paste0("ED95: ", round(Mean_ED95, 1))), 
            angle = 90, vjust = -0.4, size = 3, show.legend = FALSE) +
  
  facet_wrap(~ Species, labeller = labeller(Species = species_labels), ncol = 1) +
  scale_color_manual(values = species_colors) +
  scale_fill_manual(values = species_colors) +
  labs(
    x = "Temperature (°C)",
    y = expression("Photophysiological Yield (" * F[v] / F[m] * ")")
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.text      = element_text(face = "bold.italic", size = 11),
    axis.title      = element_text(face = "bold", size = 11, color = "black"),
    legend.position = "none"
  )

# Panel B: Average DRC by K3 Cluster
p_avg_k3 <- ggplot(drc_k3_avg, aes(x = Temperature, y = Mean_PAM, color = K_3, fill = K_3)) +
  geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.2) +
  
  # ED5 Line & Label
  geom_vline(data = eds_k3_summary, aes(xintercept = Mean_ED5, color = K_3), linetype = "dashed", alpha = 0.7) +
  geom_text(data = eds_k3_summary, aes(x = Mean_ED5, y = 0.1, label = paste0("ED5: ", round(Mean_ED5, 1))), 
            angle = 90, vjust = -0.4, size = 3, show.legend = FALSE) +
  
  # ED50 Line & Label
  geom_vline(data = eds_k3_summary, aes(xintercept = Mean_ED50, color = K_3), linetype = "dashed", alpha = 0.9) +
  geom_text(data = eds_k3_summary, aes(x = Mean_ED50, y = 0.35, label = paste0("ED50: ", round(Mean_ED50, 1))), 
            angle = 90, vjust = -0.4, size = 3, fontface = "bold", show.legend = FALSE) +
  
  # ED95 Line & Label
  geom_vline(data = eds_k3_summary, aes(xintercept = Mean_ED95, color = K_3), linetype = "dashed", alpha = 0.7) +
  geom_text(data = eds_k3_summary, aes(x = Mean_ED95, y = 0.1, label = paste0("ED95: ", round(Mean_ED95, 1))), 
            angle = 90, vjust = -0.4, size = 3, show.legend = FALSE) +
  
  facet_wrap(~ K_3, ncol = 1) +
  scale_color_manual(values = cluster_colors) +
  scale_fill_manual(values = cluster_colors) +
  labs(
    x = "Temperature (°C)",
    y = expression("Photophysiological Yield (" * F[v] / F[m] * ")")
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.text      = element_text(face = "bold", size = 11),
    axis.title      = element_text(face = "bold", size = 11, color = "black"),
    legend.position = "none"
  )

# Combine panels into Supplemental Figure S1
fig_s1_avg <- (p_avg_spp + p_avg_k3) + 
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(face = "bold", size = 14))

print(fig_s1_avg)

# Export high-resolution output
ggsave(
  here("Outputs", "Figure_S1_Averaged_DRC_Ribbons.pdf"),
  fig_s1_avg,
  width = 10,
  height = 10,
  device = "pdf"
)

library(dplyr)
library(ggplot2)
library(here)
library(patchwork)

# -----------------------------------------------------------------------------
# 1. ATTACH RAW CBASS POINTS TO SPECIES & CLUSTER LOOKUPS
# -----------------------------------------------------------------------------
# Prepare raw data points with matching species shapes & K_3 cluster assignments
df_raw_points <- df_combined %>%
  define_grouping_property(grouping_properties) %>%
  mutate(SampleID_clean = paste(GroupingProperty, Genotype, sep = "_")) %>%
  left_join(
    meta_filteredK %>% 
      dplyr::select(SampleID_clean, K_3) %>% 
      mutate(K_3 = factor(K_3, labels = c("Cluster I", "Cluster II", "Cluster III"))) %>% 
      distinct(),
    by = "SampleID_clean"
  ) %>%
  filter(!is.na(Species))

# -----------------------------------------------------------------------------
# 2. PANEL A: AVERAGED DRC BY SPECIES + RAW POINTS
# -----------------------------------------------------------------------------
p_avg_spp <- ggplot() +
  # Layer 1: Raw Observations (Shaped by Species)
  geom_point(
    data = df_raw_points,
    aes(x = Temperature, y = Pam_value, fill = Species, shape = Species),
    size = 1.8, alpha = 0.35, color = "black", stroke = 0.3
  ) +
  
  # Layer 2: 95% Confidence Ribbons
  geom_ribbon(
    data = drc_spp_avg,
    aes(x = Temperature, ymin = Lower_CI, ymax = Upper_CI, fill = Species),
    alpha = 0.25, color = NA
  ) +
  
  # Layer 3: Average Fitted DRC Lines
  geom_line(
    data = drc_spp_avg,
    aes(x = Temperature, y = Mean_PAM, color = Species),
    linewidth = 1.2
  ) +
  
  # Layer 4: Vertical ED Metric Lines & Labels
  geom_vline(data = eds_spp_summary, aes(xintercept = Mean_ED5, color = Species), linetype = "dashed", alpha = 0.7) +
  geom_text(data = eds_spp_summary, aes(x = Mean_ED5, y = 0.08, label = paste0("ED5: ", round(Mean_ED5, 1))), 
            angle = 90, vjust = -0.4, size = 2.8, color = "black", show.legend = FALSE) +
  
  geom_vline(data = eds_spp_summary, aes(xintercept = Mean_ED50, color = Species), linetype = "dashed", alpha = 0.9) +
  geom_text(data = eds_spp_summary, aes(x = Mean_ED50, y = 0.35, label = paste0("ED50: ", round(Mean_ED50, 1))), 
            angle = 90, vjust = -0.4, size = 2.8, fontface = "bold", color = "black", show.legend = FALSE) +
  
  geom_vline(data = eds_spp_summary, aes(xintercept = Mean_ED95, color = Species), linetype = "dashed", alpha = 0.7) +
  geom_text(data = eds_spp_summary, aes(x = Mean_ED95, y = 0.08, label = paste0("ED95: ", round(Mean_ED95, 1))), 
            angle = 90, vjust = -0.4, size = 2.8, color = "black", show.legend = FALSE) +
  
  facet_wrap(~ Species, labeller = labeller(Species = species_labels), ncol = 4) +
  scale_color_manual(values = species_colors) +
  scale_fill_manual(values = species_colors) +
  scale_shape_manual(values = species_shapes, limits = species_order, labels = species_labels) +
  labs(
    x = "Temperature (°C)",
    y = expression("Photophysiological Yield (" * F[v] / F[m] * ")")
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.text      = element_text(face = "bold.italic", size = 11),
    axis.title      = element_text(face = "bold", size = 11, color = "black"),
    legend.position = "none"
  )

# -----------------------------------------------------------------------------
# 3. PANEL B: AVERAGED DRC BY K3 CLUSTER + RAW POINTS
# -----------------------------------------------------------------------------
p_avg_k3 <- ggplot() +
  # Layer 1: Raw Observations (Shaped by Species, Filled by Cluster Color)
  geom_point(
    data = df_raw_points %>% filter(!is.na(K_3)),
    aes(x = Temperature, y = Pam_value, fill = K_3, shape = Species),
    size = 1.8, alpha = 0.35, color = "black", stroke = 0.3
  ) +
  
  # Layer 2: 95% Confidence Ribbons
  geom_ribbon(
    data = drc_k3_avg,
    aes(x = Temperature, ymin = Lower_CI, ymax = Upper_CI, fill = K_3),
    alpha = 0.25, color = NA
  ) +
  
  # Layer 3: Average Fitted DRC Lines
  geom_line(
    data = drc_k3_avg,
    aes(x = Temperature, y = Mean_PAM, color = K_3),
    linewidth = 1.2
  ) +
  
  # Layer 4: Vertical ED Metric Lines & Labels
  geom_vline(data = eds_k3_summary, aes(xintercept = Mean_ED5, color = K_3), linetype = "dashed", alpha = 0.7) +
  geom_text(data = eds_k3_summary, aes(x = Mean_ED5, y = 0.08, label = paste0("ED5: ", round(Mean_ED5, 1))), 
            angle = 90, vjust = -0.4, size = 2.8, color = "black", show.legend = FALSE) +
  
  geom_vline(data = eds_k3_summary, aes(xintercept = Mean_ED50, color = K_3), linetype = "dashed", alpha = 0.9) +
  geom_text(data = eds_k3_summary, aes(x = Mean_ED50, y = 0.35, label = paste0("ED50: ", round(Mean_ED50, 1))), 
            angle = 90, vjust = -0.4, size = 2.8, fontface = "bold", color = "black", show.legend = FALSE) +
  
  geom_vline(data = eds_k3_summary, aes(xintercept = Mean_ED95, color = K_3), linetype = "dashed", alpha = 0.7) +
  geom_text(data = eds_k3_summary, aes(x = Mean_ED95, y = 0.08, label = paste0("ED95: ", round(Mean_ED95, 1))), 
            angle = 90, vjust = -0.4, size = 2.8, color = "black", show.legend = FALSE) +
  
  facet_wrap(~ K_3, ncol = 4, drop = FALSE) +
  scale_color_manual(values = cluster_colors) +
  scale_fill_manual(values = cluster_colors) +
  scale_shape_manual(values = species_shapes, limits = species_order, labels = species_labels) +
  labs(
    x = "Temperature (°C)",
    y = expression("Photophysiological Yield (" * F[v] / F[m] * ")")
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.text      = element_text(face = "bold", size = 11),
    axis.title      = element_text(face = "bold", size = 11, color = "black"),
    legend.position = "none"
  )

# -----------------------------------------------------------------------------
# 4. COMBINE & EXPORT FIGURE S1
# -----------------------------------------------------------------------------
# Combine panels with explicit column alignment
fig_s1_aligned <- (p_avg_spp / p_avg_k3) + 
  plot_layout(
    ncol = 1,
    heights = c(1, 1),
    guides = "collect",
    axis_titles = "collect"
  ) + 
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(face = "bold", size = 14))

print(fig_s1_aligned)

ggsave(
  here("Outputs", "Figure_S1_DRC_Ribbons_RawPoints.pdf"),
  fig_s1_points,
  width = 10,
  height = 10,
  device = "pdf"
)

#######
p_spp_standalone <- ggplot() +
  geom_point(
    data = df_raw_points,
    aes(x = Temperature, y = Pam_value, fill = Species, shape = Species),
    size = 1.8, alpha = 0.35, color = "black", stroke = 0.3
  ) +
  geom_ribbon(
    data = drc_spp_avg,
    aes(x = Temperature, ymin = Lower_CI, ymax = Upper_CI, fill = Species),
    alpha = 0.25, color = NA
  ) +
  geom_line(
    data = drc_spp_avg,
    aes(x = Temperature, y = Mean_PAM, color = Species),
    linewidth = 1.2
  ) +
  geom_vline(data = eds_spp_summary, aes(xintercept = Mean_ED5, color = Species), linetype = "dashed", alpha = 0.7) +
  geom_text(data = eds_spp_summary, aes(x = Mean_ED5, y = 0.55, label = paste0("ED5: ", round(Mean_ED5, 1))), 
            angle = 90, vjust = -0.4, size = 2.5, color = "black", show.legend = FALSE) +
  geom_vline(data = eds_spp_summary, aes(xintercept = Mean_ED50, color = Species), linetype = "dashed", alpha = 0.9) +
  geom_text(data = eds_spp_summary, aes(x = Mean_ED50, y = 0.35, label = paste0("ED50: ", round(Mean_ED50, 1))), 
            angle = 90, vjust = -0.4, size = 2.5, fontface = "bold", color = "black", show.legend = FALSE) +
  geom_vline(data = eds_spp_summary, aes(xintercept = Mean_ED95, color = Species), linetype = "dashed", alpha = 0.7) +
  geom_text(data = eds_spp_summary, aes(x = Mean_ED95, y = 0.15, label = paste0("ED95: ", round(Mean_ED95, 1))), 
            angle = 90, vjust = -0.4, size = 2.5, color = "black", show.legend = FALSE) +
  
  facet_wrap(~ Species, labeller = labeller(Species = species_labels), ncol = 4) +
  scale_color_manual(values = species_colors) +
  scale_fill_manual(values = species_colors) +
  scale_shape_manual(values = species_shapes, limits = species_order, labels = species_labels) +
  scale_x_continuous(limits = c(29, 41), breaks = c(32, 36, 40)) +
  scale_y_continuous(limits = c(0, 0.75), breaks = seq(0, 0.6, by = 0.2)) +
  labs(
    x = "Temperature (°C)",
    y = expression("Photophysiological Yield (" * F[v] / F[m] * ")")
  ) +
  theme_classic(base_size = 12) +
  theme(
    aspect.ratio     = 1, # LOCKS PANEL SHAPE TO PERFECT SQUARE
    strip.text       = element_text(face = "bold.italic", size = 11),
    strip.background = element_rect(fill = "gray95", color = "black", linewidth = 0.8),
    axis.title       = element_text(face = "bold", size = 11, color = "black"),
    axis.text        = element_text(size = 10, color = "black"),
    legend.position  = "none"
  )

print(p_spp_standalone)

p_k3_standalone <- ggplot() +
  geom_point(
    data = df_raw_points %>% filter(!is.na(K_3)),
    aes(x = Temperature, y = Pam_value, fill = K_3, shape = Species),
    size = 1.8, alpha = 0.35, color = "black", stroke = 0.3
  ) +
  geom_ribbon(
    data = drc_k3_avg,
    aes(x = Temperature, ymin = Lower_CI, ymax = Upper_CI, fill = K_3),
    alpha = 0.25, color = NA
  ) +
  geom_line(
    data = drc_k3_avg,
    aes(x = Temperature, y = Mean_PAM, color = K_3),
    linewidth = 1.2
  ) +
  geom_vline(data = eds_k3_summary, aes(xintercept = Mean_ED5, color = K_3), linetype = "dashed", alpha = 0.7) +
  geom_text(data = eds_k3_summary, aes(x = Mean_ED5, y = 0.55, label = paste0("ED5: ", round(Mean_ED5, 1))), 
            angle = 90, vjust = -0.4, size = 2.5, color = "black", show.legend = FALSE) +
  geom_vline(data = eds_k3_summary, aes(xintercept = Mean_ED50, color = K_3), linetype = "dashed", alpha = 0.9) +
  geom_text(data = eds_k3_summary, aes(x = Mean_ED50, y = 0.35, label = paste0("ED50: ", round(Mean_ED50, 1))), 
            angle = 90, vjust = -0.4, size = 2.5, fontface = "bold", color = "black", show.legend = FALSE) +
  geom_vline(data = eds_k3_summary, aes(xintercept = Mean_ED95, color = K_3), linetype = "dashed", alpha = 0.7) +
  geom_text(data = eds_k3_summary, aes(x = Mean_ED95, y = 0.15, label = paste0("ED95: ", round(Mean_ED95, 1))), 
            angle = 90, vjust = -0.4, size = 2.5, color = "black", show.legend = FALSE) +
  
  facet_wrap(~ K_3, ncol = 3) +
  scale_color_manual(values = cluster_colors) +
  scale_fill_manual(values = cluster_colors) +
  scale_shape_manual(values = species_shapes, limits = species_order, labels = species_labels) +
  scale_x_continuous(limits = c(29, 41), breaks = c(32, 36, 40)) +
  scale_y_continuous(limits = c(0, 0.75), breaks = seq(0, 0.6, by = 0.2)) +
  labs(
    x = "Temperature (°C)",
    y = expression("Photophysiological Yield (" * F[v] / F[m] * ")")
  ) +
  theme_classic(base_size = 12) +
  theme(
    aspect.ratio     = 1, # LOCKS PANEL SHAPE TO PERFECT SQUARE (MATCHES PLOT A)
    strip.text       = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "gray95", color = "black", linewidth = 0.8),
    axis.title       = element_text(face = "bold", size = 11, color = "black"),
    axis.text        = element_text(size = 10, color = "black"),
    legend.position  = "none"
  )

print(p_k3_standalone)

library(dplyr)
library(ggplot2)
library(here)
library(patchwork)

# -----------------------------------------------------------------------------
# 1. PANEL A: SPECIES PROFILES (4 VERTICAL PANELS)
# -----------------------------------------------------------------------------
p_avg_spp_vert <- ggplot() +
  # Raw observations
  geom_point(
    data = df_raw_points,
    aes(x = Temperature, y = Pam_value, fill = Species, shape = Species),
    size = 1.8, alpha = 0.35, color = "black", stroke = 0.3
  ) +
  # 95% Confidence Ribbons
  geom_ribbon(
    data = drc_spp_avg,
    aes(x = Temperature, ymin = Lower_CI, ymax = Upper_CI, fill = Species),
    alpha = 0.25, color = NA
  ) +
  # Fitted DRC Lines
  geom_line(
    data = drc_spp_avg,
    aes(x = Temperature, y = Mean_PAM, color = Species),
    linewidth = 1.2
  ) +
  # Vertical ED Lines
  geom_vline(data = eds_spp_summary, aes(xintercept = Mean_ED5, color = Species), linetype = "dashed", alpha = 0.7) +
  geom_vline(data = eds_spp_summary, aes(xintercept = Mean_ED50, color = Species), linetype = "dashed", alpha = 0.9) +
  geom_vline(data = eds_spp_summary, aes(xintercept = Mean_ED95, color = Species), linetype = "dashed", alpha = 0.7) +
  
  # Floating Text Labels ABOVE 0.6 on the same line (y = 0.72)
  geom_text(data = eds_spp_summary, aes(x = Mean_ED5, y = 0.72, label = paste0("5: ", round(Mean_ED5, 1))), 
            angle = 90, vjust = -0.3, size = 2.7, color = "black", show.legend = FALSE) +
  geom_text(data = eds_spp_summary, aes(x = Mean_ED50, y = 0.72, label = paste0("50: ", round(Mean_ED50, 1))), 
            angle = 90, vjust = -0.3, size = 2.7, fontface = "bold", color = "black", show.legend = FALSE) +
  geom_text(data = eds_spp_summary, aes(x = Mean_ED95, y = 0.72, label = paste0("95: ", round(Mean_ED95, 1))), 
            angle = 90, vjust = -0.3, size = 2.7, color = "black", show.legend = FALSE) +
  
  facet_wrap(~ Species, labeller = labeller(Species = species_labels), ncol = 1) +
  scale_color_manual(values = species_colors) +
  scale_fill_manual(values = species_colors) +
  scale_shape_manual(values = species_shapes, limits = species_order, labels = species_labels) +
  scale_x_continuous(limits = c(29, 41), breaks = c(32, 36, 40)) +
  scale_y_continuous(limits = c(0, 0.85), breaks = seq(0, 0.6, by = 0.2)) +
  labs(
    x = "Temperature (°C)",
    y = expression("Photophysiological Yield (" * F[v] / F[m] * ")")
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.text       = element_text(face = "bold.italic", size = 11),
    strip.background = element_rect(fill = "gray95", color = "black", linewidth = 0.8),
    axis.title       = element_text(face = "bold", size = 11, color = "black"),
    axis.text        = element_text(size = 10, color = "black"),
    legend.position  = "none"
  )

# -----------------------------------------------------------------------------
# 2. PANEL B: K3 CLUSTER PROFILES (3 VERTICAL PANELS)
# -----------------------------------------------------------------------------
p_avg_k3_vert <- ggplot() +
  # Raw observations
  geom_point(
    data = df_raw_points %>% filter(!is.na(K_3)),
    aes(x = Temperature, y = Pam_value, fill = K_3, shape = Species),
    size = 1.8, alpha = 0.35, color = "black", stroke = 0.3
  ) +
  # 95% Confidence Ribbons
  geom_ribbon(
    data = drc_k3_avg,
    aes(x = Temperature, ymin = Lower_CI, ymax = Upper_CI, fill = K_3),
    alpha = 0.25, color = NA
  ) +
  # Fitted DRC Lines
  geom_line(
    data = drc_k3_avg,
    aes(x = Temperature, y = Mean_PAM, color = K_3),
    linewidth = 1.2
  ) +
  # Vertical ED Lines
  geom_vline(data = eds_k3_summary, aes(xintercept = Mean_ED5, color = K_3), linetype = "dashed", alpha = 0.7) +
  geom_vline(data = eds_k3_summary, aes(xintercept = Mean_ED50, color = K_3), linetype = "dashed", alpha = 0.9) +
  geom_vline(data = eds_k3_summary, aes(xintercept = Mean_ED95, color = K_3), linetype = "dashed", alpha = 0.7) +
  
  # Floating Text Labels ABOVE 0.6 on the same line (y = 0.72)
  geom_text(data = eds_k3_summary, aes(x = Mean_ED5, y = 0.72, label = paste0("5: ", round(Mean_ED5, 1))), 
            angle = 90, vjust = -0.3, size = 2.7, color = "black", show.legend = FALSE) +
  geom_text(data = eds_k3_summary, aes(x = Mean_ED50, y = 0.72, label = paste0("50: ", round(Mean_ED50, 1))), 
            angle = 90, vjust = -0.3, size = 2.7, fontface = "bold", color = "black", show.legend = FALSE) +
  geom_text(data = eds_k3_summary, aes(x = Mean_ED95, y = 0.72, label = paste0("95: ", round(Mean_ED95, 1))), 
            angle = 90, vjust = -0.3, size = 2.7, color = "black", show.legend = FALSE) +
  
  facet_wrap(~ K_3, ncol = 1) +
  scale_color_manual(values = cluster_colors) +
  scale_fill_manual(values = cluster_colors) +
  scale_shape_manual(values = species_shapes, limits = species_order, labels = species_labels) +
  scale_x_continuous(limits = c(29, 41), breaks = c(32, 36, 40)) +
  scale_y_continuous(limits = c(0, 0.85), breaks = seq(0, 0.6, by = 0.2)) +
  labs(
    x = "Temperature (°C)",
    y = expression("Photophysiological Yield (" * F[v] / F[m] * ")")
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.text       = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "gray95", color = "black", linewidth = 0.8),
    axis.title       = element_text(face = "bold", size = 11, color = "black"),
    axis.text        = element_text(size = 10, color = "black"),
    legend.position  = "none"
  )

# -----------------------------------------------------------------------------
# 3. FORCE PERFECT 1:1 PANEL DIMENSIONS VIA PATCHWORK DESIGN
# -----------------------------------------------------------------------------
# Map 'A' to 4 units tall, 'B' to 3 units tall, leaving 1 unit empty ('#')
layout_design <- "
  A B
  A B
  A B
  A #
"

fig_s1_perfect_vertical <- p_avg_spp_vert + p_avg_k3_vert + 
  plot_layout(design = layout_design) + 
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(face = "bold", size = 14))

print(fig_s1_perfect_vertical)

# Export PDF
ggsave(
  here("Outputs", "Figure_S1_Vertical_PerfectDimensions.pdf"),
  fig_s1_perfect_vertical,
  width = 8.5,
  height = 11,
  device = "pdf"
)

# ==============================================================================
# DOSE-RESPONSE CURVE PLOTTING PIPELINE
# (Faceted Vertical Layout + Unfaceted Overlay Comparison)
# ==============================================================================
library(dplyr)
library(tidyr)
library(ggplot2)
library(here)
library(patchwork)
library(CBASSED50)

# -----------------------------------------------------------------------------
# 1. LOAD PREPROCESSED DATA & FIT GLOBAL MODELS
# -----------------------------------------------------------------------------
source(here("Scripts/00_visualization_prep.R"))
load(here("Outputs", "photophys_preprocessed_data.RData"))

# Define file paths to raw CBASS/PAM run files
run_files <- list.files(
  path       = here("Data"),
  pattern    = "\\PAM.csv$",
  full.names = TRUE
)

# Combine raw site files safely
if (length(run_files) > 0) {
  df_combined <- map_dfr(run_files, function(file_path) {
    df_run <- read_data(file_path)
    df_run %>%
      mutate(
        PAM1 = ifelse(is.na(PAM1), Pam_value, PAM1),
        PAM2 = ifelse(is.na(PAM2), Pam_value, PAM2)
      ) %>%
      preprocess_dataset()
  })
}

grouping_properties <- c("Site", "Species")
drm_formula         <- "Pam_value ~ Temperature"

# Fit global models and predict PAM yield profiles
models_global <- fit_drms(df_combined, grouping_properties, drm_formula, is_curveid = FALSE, LL.4 = FALSE)
temp_ranges   <- define_temperature_ranges(df_combined$Temperature, n = 100)
predictions   <- get_predicted_pam_values(models_global, temp_ranges)

# -----------------------------------------------------------------------------
# 2. JOIN METADATA VIA SampleID_clean
# -----------------------------------------------------------------------------
df_combined_grouped <- define_grouping_property(df_combined, grouping_properties) %>%
  mutate(SampleID_clean = paste(GroupingProperty, Genotype, sep = "_")) %>%
  dplyr::select(GroupingProperty, SampleID_clean, Site, Species) %>%
  distinct()

# Prepare cluster metadata lookup
cluster_lookup <- meta_filteredK %>%
  dplyr::select(SampleID_clean, K_3) %>%
  mutate(
    K_3 = factor(K_3, labels = c("Cluster I", "Cluster II", "Cluster III"))
  ) %>%
  distinct()

# Combine predictions with meta_filteredK
predictions_df <- predictions %>%
  left_join(df_combined_grouped, by = "GroupingProperty") %>%
  left_join(cluster_lookup, by = "SampleID_clean") %>%
  filter(!is.na(Species))

# Attach raw observations for scatter points
df_raw_points <- df_combined %>%
  define_grouping_property(grouping_properties) %>%
  mutate(SampleID_clean = paste(GroupingProperty, Genotype, sep = "_")) %>%
  left_join(cluster_lookup, by = "SampleID_clean") %>%
  filter(!is.na(Species))

# -----------------------------------------------------------------------------
# 3. AGGREGATE SUMMARY METRICS (ED5, ED50, ED95) AND PREDICTIONS
# -----------------------------------------------------------------------------
# Means by Species
eds_spp_summary <- meta_filteredK %>%
  group_by(Species) %>%
  summarise(
    Mean_ED5  = mean(ED5, na.rm = TRUE),
    Mean_ED50 = mean(ED50, na.rm = TRUE),
    Mean_ED95 = mean(ED95, na.rm = TRUE),
    .groups   = "drop"
  )

# Means by Cluster
# Fix the factor level mapping in the summary table
eds_k3_summary <- meta_filteredK %>%
  filter(!is.na(K_3)) %>%
  mutate(
    # Map factor levels explicitly to target_levels ("Cluster I", "Cluster II", "Cluster III")
    K_3 = roman_map[as.character(K_3)],
    K_3 = factor(K_3, levels = target_levels)
  ) %>%
  group_by(K_3) %>%
  summarise(
    Mean_ED5  = mean(ED5, na.rm = TRUE),
    Mean_ED50 = mean(ED50, na.rm = TRUE),
    Mean_ED95 = mean(ED95, na.rm = TRUE),
    .groups   = "drop"
  )

# Aggregate Mean Yield & CIs for Species
drc_spp_avg <- predictions_df %>%
  group_by(Species, Temperature) %>%
  summarise(
    Mean_PAM = mean(PredictedPAM, na.rm = TRUE),
    Lower_CI = quantile(PredictedPAM, 0.025, na.rm = TRUE),
    Upper_CI = quantile(PredictedPAM, 0.975, na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  left_join(eds_spp_summary, by = "Species")

# Aggregate Mean Yield & CIs for K3 Cluster
drc_k3_avg <- predictions_df %>%
  filter(!is.na(K_3)) %>%
  group_by(K_3, Temperature) %>%
  summarise(
    Mean_PAM = mean(PredictedPAM, na.rm = TRUE),
    Lower_CI = quantile(PredictedPAM, 0.025, na.rm = TRUE),
    Upper_CI = quantile(PredictedPAM, 0.975, na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  left_join(eds_k3_summary, by = "K_3")

# ==============================================================================
# FIGURE 1: FACETED VERTICAL COLUMNS (ALIGNING EDs)
# ==============================================================================

# Column A: Species (4 vertical panels)
p_spp_vert <- ggplot() +
  geom_point(
    data = df_raw_points,
    aes(x = Temperature, y = Pam_value, fill = Species, shape = Species),
    size = 1.6, alpha = 0.30, color = "black", stroke = 0.25
  ) +
  geom_ribbon(
    data = drc_spp_avg,
    aes(x = Temperature, ymin = Lower_CI, ymax = Upper_CI, fill = Species),
    alpha = 0.25, color = NA
  ) +
  geom_line(
    data = drc_spp_avg,
    aes(x = Temperature, y = Mean_PAM, color = Species),
    linewidth = 1.1
  ) +
  # Vertical ED Lines
  geom_vline(data = eds_spp_summary, aes(xintercept = Mean_ED5, color = Species), linetype = "dashed", alpha = 0.7) +
  geom_vline(data = eds_spp_summary, aes(xintercept = Mean_ED50, color = Species), linetype = "dashed", alpha = 0.9) +
  geom_vline(data = eds_spp_summary, aes(xintercept = Mean_ED95, color = Species), linetype = "dashed", alpha = 0.7) +
  # Floating Labels above 0.6
  geom_text(data = eds_spp_summary, aes(x = Mean_ED5, y = 0.72, label = paste0("5: ", round(Mean_ED5, 1))), 
            angle = 90, vjust = -0.3, size = 2.6, color = "black", show.legend = FALSE) +
  geom_text(data = eds_spp_summary, aes(x = Mean_ED50, y = 0.72, label = paste0("50: ", round(Mean_ED50, 1))), 
            angle = 90, vjust = -0.3, size = 2.6, fontface = "bold", color = "black", show.legend = FALSE) +
  geom_text(data = eds_spp_summary, aes(x = Mean_ED95, y = 0.72, label = paste0("95: ", round(Mean_ED95, 1))), 
            angle = 90, vjust = -0.3, size = 2.6, color = "black", show.legend = FALSE) +
  
  facet_wrap(~ Species, labeller = labeller(Species = species_labels), ncol = 1) +
  scale_color_manual(values = species_colors) +
  scale_fill_manual(values = species_colors) +
  scale_shape_manual(values = species_shapes, limits = species_order, labels = species_labels) +
  scale_x_continuous(limits = c(29, 41), breaks = c(32, 36, 40)) +
  scale_y_continuous(limits = c(0, 0.85), breaks = seq(0, 0.6, by = 0.2)) +
  labs(
    x = "Temperature (°C)",
    y = expression("Photophysiological Yield (" * F[v] / F[m] * ")")
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.text       = element_text(face = "bold.italic", size = 11),
    strip.background = element_rect(fill = "gray95", color = "black", linewidth = 0.8),
    axis.title       = element_text(face = "bold", size = 11, color = "black"),
    axis.text        = element_text(size = 10, color = "black"),
    legend.position  = "none"
  )

# Column B: Clusters (3 vertical panels)
p_k3_vert <- ggplot() +
  geom_point(
    data = df_raw_points %>% filter(!is.na(K_3)),
    aes(x = Temperature, y = Pam_value, fill = K_3, shape = Species),
    size = 1.6, alpha = 0.30, color = "black", stroke = 0.25
  ) +
  geom_ribbon(
    data = drc_k3_avg,
    aes(x = Temperature, ymin = Lower_CI, ymax = Upper_CI, fill = K_3),
    alpha = 0.25, color = NA
  ) +
  geom_line(
    data = drc_k3_avg,
    aes(x = Temperature, y = Mean_PAM, color = K_3),
    linewidth = 1.1
  ) +
  # Vertical ED Lines
  geom_vline(data = eds_k3_summary, aes(xintercept = Mean_ED5, color = K_3), linetype = "dashed", alpha = 0.7) +
  geom_vline(data = eds_k3_summary, aes(xintercept = Mean_ED50, color = K_3), linetype = "dashed", alpha = 0.9) +
  geom_vline(data = eds_k3_summary, aes(xintercept = Mean_ED95, color = K_3), linetype = "dashed", alpha = 0.7) +
  # Floating Labels above 0.6
  geom_text(data = eds_k3_summary, aes(x = Mean_ED5, y = 0.72, label = paste0("5: ", round(Mean_ED5, 1))), 
            angle = 90, vjust = -0.3, size = 2.6, color = "black", show.legend = FALSE) +
  geom_text(data = eds_k3_summary, aes(x = Mean_ED50, y = 0.72, label = paste0("50: ", round(Mean_ED50, 1))), 
            angle = 90, vjust = -0.3, size = 2.6, fontface = "bold", color = "black", show.legend = FALSE) +
  geom_text(data = eds_k3_summary, aes(x = Mean_ED95, y = 0.72, label = paste0("95: ", round(Mean_ED95, 1))), 
            angle = 90, vjust = -0.3, size = 2.6, color = "black", show.legend = FALSE) +
  
  facet_wrap(~ K_3, ncol = 1) +
  scale_color_manual(values = cluster_colors) +
  scale_fill_manual(values = cluster_colors) +
  scale_shape_manual(values = species_shapes, limits = species_order, labels = species_labels) +
  scale_x_continuous(limits = c(29, 41), breaks = c(32, 36, 40)) +
  scale_y_continuous(limits = c(0, 0.85), breaks = seq(0, 0.6, by = 0.2)) +
  labs(
    x = "Temperature (°C)",
    y = expression("Photophysiological Yield (" * F[v] / F[m] * ")")
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.text       = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "gray95", color = "black", linewidth = 0.8),
    axis.title       = element_text(face = "bold", size = 11, color = "black"),
    axis.text        = element_text(size = 10, color = "black"),
    legend.position  = "none"
  )

# Combine using patchwork design layout for perfect 1:1 panel alignment
layout_design <- "
  A B
  A B
  A B
  A #
"

fig_s1_faceted_vertical <- p_spp_vert + p_k3_vert + 
  plot_layout(design = layout_design) + 
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(face = "bold", size = 14))

print(fig_s1_faceted_vertical)



# ==============================================================================
# FIGURE 5: UNFACETED OVERLAY PLOT (DIRECT LINE COMPARISONS)
# ==============================================================================
# ==============================================================================
# FIGURE 2: UNFACETED OVERLAY PLOT WITH VERTICAL ED LINES
# ==============================================================================

# Panel A: Unfaceted Species Overlay with Vertical ED Lines
p_spp_overlay <- ggplot(drc_spp_avg, aes(x = Temperature, y = Mean_PAM, color = Species, fill = Species)) +
  # 95% Confidence Ribbons
  geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), alpha = 0.15, color = NA) +
  
  # Dose-Response Curves
  geom_line(linewidth = 1.3) +
  
  # Vertical ED Lines (No Text Labels)
  geom_vline(data = eds_spp_summary, aes(xintercept = Mean_ED5, color = Species), linetype = "dotted", alpha = 0.65, linewidth = 0.7) +
  geom_vline(data = eds_spp_summary, aes(xintercept = Mean_ED50, color = Species), linetype = "dashed", alpha = 0.85, linewidth = 0.9) +
  geom_vline(data = eds_spp_summary, aes(xintercept = Mean_ED95, color = Species), linetype = "dotdash", alpha = 0.65, linewidth = 0.7) +
  
  # Manual Aesthetics & Scales
  scale_color_manual(name = "Species", values = species_colors, limits = species_order, labels = species_labels) +
  scale_fill_manual(name = "Species", values = species_colors, limits = species_order, labels = species_labels) +
  scale_x_continuous(limits = c(29, 41), breaks = c(30, 32, 34, 36, 38, 40)) +
  scale_y_continuous(limits = c(0, 0.70), breaks = seq(0, 0.6, by = 0.2)) +
  labs(
    x = "Temperature (°C)",
    y = expression("Photophysiological Yield (" * F[v] / F[m] * ")")
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.title      = element_text(face = "bold", size = 12, color = "black"),
    axis.text       = element_text(size = 11, color = "black"),
    legend.title    = element_text(face = "bold", size = 11),
    legend.text     = element_text(size = 10, face = "italic"),
    legend.position = "right"
  )

# Panel B: Unfaceted Cluster Overlay with Vertical ED Lines
p_k3_overlay <- ggplot(drc_k3_avg, aes(x = Temperature, y = Mean_PAM, color = K_3, fill = K_3)) +
  # 95% Confidence Ribbons
  geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), alpha = 0.18, color = NA) +
  
  # Dose-Response Curves
  geom_line(linewidth = 1.3) +
  
  # Vertical ED Lines (No Text Labels)
  geom_vline(data = eds_k3_summary, aes(xintercept = Mean_ED5, color = K_3), linetype = "dotted", alpha = 0.65, linewidth = 0.7) +
  geom_vline(data = eds_k3_summary, aes(xintercept = Mean_ED50, color = K_3), linetype = "dashed", alpha = 0.85, linewidth = 0.9) +
  geom_vline(data = eds_k3_summary, aes(xintercept = Mean_ED95, color = K_3), linetype = "dotdash", alpha = 0.65, linewidth = 0.7) +
  
  # Manual Aesthetics & Scales
  scale_color_manual(name = "Cluster", values = cluster_colors) +
  scale_fill_manual(name = "Cluster", values = cluster_colors) +
  scale_x_continuous(limits = c(29, 41), breaks = c(30, 32, 34, 36, 38, 40)) +
  scale_y_continuous(limits = c(0, 0.70), breaks = seq(0, 0.6, by = 0.2)) +
  labs(
    x = "Temperature (°C)",
    y = expression("Photophysiological Yield (" * F[v] / F[m] * ")")
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.title      = element_text(face = "bold", size = 12, color = "black"),
    axis.text       = element_text(size = 11, color = "black"),
    legend.title    = element_text(face = "bold", size = 11),
    legend.text     = element_text(size = 10),
    legend.position = "right"
  )

# Combine into Unfaceted Overlay Figure with Vertical ED Lines
fig_s2_overlay_unfaceted <- (p_spp_overlay / p_k3_overlay) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 14))

print(fig_s2_overlay_unfaceted)

# -----------------------------------------------------------------------------
# AGGREGATE BIN MEANS & SE AT DISCRETE TEMPERATURE STEPS
# -----------------------------------------------------------------------------
# Bin means for Species
raw_spp_bins <- df_raw_points %>%
  mutate(Temp_bin = round(Temperature)) %>% # Rounds to target CBASS temp steps
  group_by(Species, Temp_bin) %>%
  summarise(
    Mean_Yield = mean(Pam_value, na.rm = TRUE),
    SE_Yield   = sd(Pam_value, na.rm = TRUE) / sqrt(n()),
    .groups    = "drop"
  )

# Bin means for K3 Clusters
raw_k3_bins <- df_raw_points %>%
  filter(!is.na(K_3)) %>%
  mutate(Temp_bin = round(Temperature)) %>%
  group_by(K_3, Temp_bin) %>%
  summarise(
    Mean_Yield = mean(Pam_value, na.rm = TRUE),
    SE_Yield   = sd(Pam_value, na.rm = TRUE) / sqrt(n()),
    .groups    = "drop"
  )

# -----------------------------------------------------------------------------
# PANEL A: SPECIES OVERLAY WITH TEMP-BIN AVERAGES
# -----------------------------------------------------------------------------
p_spp_overlay_bins <- ggplot() +
  # 1. Temperature Bin Means & Standard Error Bars
  geom_errorbar(
    data = raw_spp_bins,
    aes(x = Temp_bin, ymin = Mean_Yield - SE_Yield, ymax = Mean_Yield + SE_Yield, color = Species),
    width = 0.2, alpha = 0.6, linewidth = 0.5
  ) +
  geom_point(
    data = raw_spp_bins,
    aes(x = Temp_bin, y = Mean_Yield, fill = Species, shape = Species),
    size = 2.5, alpha = 0.85, color = "black", stroke = 0.4
  ) +
  
  # 2. 95% Confidence Ribbons & Fitted Curves
  geom_ribbon(data = drc_spp_avg, aes(x = Temperature, ymin = Lower_CI, ymax = Upper_CI, fill = Species), alpha = 0.12, color = NA) +
  geom_line(data = drc_spp_avg, aes(x = Temperature, y = Mean_PAM, color = Species), linewidth = 1.2) +
  
  # 3. Vertical ED Lines
  geom_vline(data = eds_spp_summary, aes(xintercept = Mean_ED5, color = Species), linetype = "dotted", alpha = 0.6, linewidth = 0.6) +
  geom_vline(data = eds_spp_summary, aes(xintercept = Mean_ED50, color = Species), linetype = "dashed", alpha = 0.85, linewidth = 0.8) +
  geom_vline(data = eds_spp_summary, aes(xintercept = Mean_ED95, color = Species), linetype = "dotdash", alpha = 0.6, linewidth = 0.6) +
  
  # Scales & Formatting
  scale_color_manual(name = "Species", values = species_colors, limits = species_order, labels = species_labels) +
  scale_fill_manual(name = "Species", values = species_colors, limits = species_order, labels = species_labels) +
  scale_shape_manual(name = "Species", values = species_shapes, limits = species_order, labels = species_labels) +
  scale_x_continuous(limits = c(29, 41), breaks = c(30, 32, 34, 36, 38, 40)) +
  scale_y_continuous(limits = c(0, 0.72), breaks = seq(0, 0.6, by = 0.2)) +
  labs(x = "Temperature (°C)", y = expression("Photophysiological Yield (" * F[v] / F[m] * ")")) +
  theme_classic(base_size = 12) +
  theme(
    axis.title      = element_text(face = "bold", size = 12, color = "black"),
    axis.text       = element_text(size = 11, color = "black"),
    legend.title    = element_text(face = "bold", size = 11),
    legend.text     = element_text(size = 10, face = "italic"),
    legend.position = "right"
  )

# -----------------------------------------------------------------------------
# PANEL B: CLUSTER OVERLAY WITH TEMP-BIN AVERAGES
# -----------------------------------------------------------------------------
p_k3_overlay_bins <- ggplot() +
  # 1. Temperature Bin Means & Standard Error Bars
  geom_errorbar(
    data = raw_k3_bins,
    aes(x = Temp_bin, ymin = Mean_Yield - SE_Yield, ymax = Mean_Yield + SE_Yield, color = K_3),
    width = 0.2, alpha = 0.6, linewidth = 0.5
  ) +
  geom_point(
    data = raw_k3_bins,
    aes(x = Temp_bin, y = Mean_Yield, fill = K_3),
    shape = 21, size = 2.5, alpha = 0.85, color = "black", stroke = 0.4
  ) +
  
  # 2. 95% Confidence Ribbons & Fitted Curves
  geom_ribbon(data = drc_k3_avg, aes(x = Temperature, ymin = Lower_CI, ymax = Upper_CI, fill = K_3), alpha = 0.12, color = NA) +
  geom_line(data = drc_k3_avg, aes(x = Temperature, y = Mean_PAM, color = K_3), linewidth = 1.2) +
  
  # 3. Vertical ED Lines
  geom_vline(data = eds_k3_summary, aes(xintercept = Mean_ED5, color = K_3), linetype = "dotted", alpha = 0.6, linewidth = 0.6) +
  geom_vline(data = eds_k3_summary, aes(xintercept = Mean_ED50, color = K_3), linetype = "dashed", alpha = 0.85, linewidth = 0.8) +
  geom_vline(data = eds_k3_summary, aes(xintercept = Mean_ED95, color = K_3), linetype = "dotdash", alpha = 0.6, linewidth = 0.6) +
  
  # Scales & Formatting
  scale_color_manual(name = "Cluster", values = cluster_colors) +
  scale_fill_manual(name = "Cluster", values = cluster_colors) +
  scale_x_continuous(limits = c(29, 41), breaks = c(30, 32, 34, 36, 38, 40)) +
  scale_y_continuous(limits = c(0, 0.72), breaks = seq(0, 0.6, by = 0.2)) +
  labs(x = "Temperature (°C)", y = expression("Photophysiological Yield (" * F[v] / F[m] * ")")) +
  theme_classic(base_size = 12) +
  theme(
    axis.title      = element_text(face = "bold", size = 12, color = "black"),
    axis.text       = element_text(size = 11, color = "black"),
    legend.title    = element_text(face = "bold", size = 11),
    legend.text     = element_text(size = 10),
    legend.position = "right"
  )

# Combine
fig_overlay_bin_averages <- (p_spp_overlay_bins / p_k3_overlay_bins) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 14))

print(fig_overlay_bin_averages)
