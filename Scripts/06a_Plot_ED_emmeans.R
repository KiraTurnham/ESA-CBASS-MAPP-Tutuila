rm(list=ls())
library(dplyr)
library(lme4)
library(car)
library(performance)
library(stringr)
library(lmerTest)
library(broom.mixed)
library(emmeans)
library(multcomp)
library(here)
library(sjPlot)
source(here("Scripts/00_visualization_prep.R"))
load(here("Outputs", "photophys_preprocessed_data.RData"))

df_modeling <- meta_filteredK %>%
  filter(!is.na(color_slope), !is.na(Site), !is.na(Species)) %>%
  mutate(
    YY      = ED50,
    Species = droplevels(as.factor(Species)),
    Site    = droplevels(as.factor(Site)),
    across(starts_with("K_"), ~ droplevels(as.factor(.)))
  )

#model with site and species as random factor
k3ms <- lmer(YY ~  Species + K_3 + (1 | Site), data = df_modeling, REML = TRUE)

print(anova(k3ms, type = 2))

supp_coefs <- tidy(k3ms, effects = "fixed", conf.int = TRUE) %>%
  transmute(
    Term = term,
    Estimate = round(estimate, 3),
    `Std. Error` = round(std.error, 3),
    `2.5% CI` = round(conf.low, 3),
    `97.5% CI` = round(conf.high, 3),
    `p-value` = round(p.value, 4)
  )
performance::r2(k3ms)

# coefs
coefs <- broom.mixed::tidy(k3ms, effects = "fixed", conf.int = TRUE) %>%
  transmute(
    Term          = term,
    Estimate      = round(estimate, 3),
    `Std. Error`  = round(std.error, 3),
    `2.5% CI`     = round(conf.low, 3),
    `97.5% CI`    = round(conf.high, 3),
    `Statistic`   = round(statistic, 2),
    `p-value`     = ifelse(p.value < 0.001, "< 0.001", sprintf("%.4f", p.value))
  );coefs

# random
random_vars <- broom.mixed::tidy(k3ms, effects = "ran_pars") %>%
  transmute(
    Group    = group,
    Term     = term,
    StdDev   = round(estimate, 3),
    Variance = round(estimate^2, 3)
  );random_vars

plot_model(k3ms, type = "re")
# ==============================================================================
# Diagnostics
# ==============================================================================
# --- Residual Diagnostics ---
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
plot(resid(k3ms) ~ fitted(k3ms), main = "Residuals vs Fitted (k3msr)", xlab = "Fitted", ylab = "Residuals")
qqnorm(resid(k3ms), main = "Normal Q-Q Plot (k3ms)")
qqline(resid(k3ms), col = "red")
par(mfrow = c(1, 1))

# ==============================================================================
# 2. FOR SPP: Extract EMMs & Compact Letter Display (CLD)
# ==============================================================================
k3_emmeans_spp <- emmeans(k3ms, pairwise ~ Species, adjust = "tukey")
k3_emm_df       <- as.data.frame(k3_emmeans_spp$emmeans)
k3_contrasts_df <- as.data.frame(k3_emmeans_spp$contrasts)

cat("\n--- Estimated Marginal Means (ED50) for K_3 Clusters ---\n")
print(k3_emm_df %>% dplyr::select(Species, emmean, SE, df, lower.CL, upper.CL))

cat("\n--- Pairwise Tukey Contrasts Between K_3 Clusters ---\n")
print(k3_contrasts_df %>% dplyr::select(contrast, estimate, SE, p.value))

# Compute Tukey CLD Letters
cld_spp   <- cld(k3_emmeans_spp$emmeans, adjust = "tukey", Letters = letters) %>%
  as.data.frame() %>%
  mutate(.group = trimws(.group))

# # ==============================================================================
# # 3. Align Cluster Labels & Factor Ordering
# # ==============================================================================
# Verify counts in spp
cat("Raw Data Counts per spp Label:\n")
print(table(df_modeling$Species, useNA = "always"))

# ==============================================================================
# 4. Publication-Ready Plot
# ==============================================================================
y_min_val <- min(df_modeling$ED50, na.rm = TRUE) - 0.5
y_max_val <- max(cld_spp$upper.CL, na.rm = TRUE) + 0.8

plot_k3spp <- ggplot() +
  # Layer A: Raw data points
  # geom_jitter(
  #   data = df_modeling, 
  #   aes(x = Species, y = ED50, color = Species),
  #   width = 0.12, 
  #   alpha = 0.5, 
  #   size = 2
  # ) +
  # Layer B: Model Estimated Marginal Means & 95% CIs
  geom_errorbar(
    data = cld_spp,
    aes(x = Species, ymin = lower.CL, ymax = upper.CL),
    width = 0.15, 
    linewidth = 0.8, 
    color = "black"
  ) +
  geom_point(
    data = cld_spp,
    aes(x = Species, y = emmean, fill = Species, shape = Species),
    size = 5,
    #fill = Species,
    color = "black",
    #shape = 21
  ) +
  # Layer C: Significance Letters (CLD)
  geom_text(
    data = cld_spp,
    aes(x = Species, y = upper.CL + 0.5, label = .group),
    #fontface = "bold", 
    hjust = 0.5,
    size = 6,
    color = "black",
    inherit.aes = FALSE
  ) +
  # Custom Aesthetics & Formatting
  coord_cartesian(ylim = c(y_min_val, y_max_val)) +
  labs(
    y = expression("ED50" ~ "(°C)"),
    color = "Species"
  ) +
  
  # FORCES X-AXIS ORDER & ITALIC SPECIES LABELS
  scale_x_discrete(
    limits = species_order, 
    labels = species_labels
  ) +
  
  scale_color_manual(values = species_colors, limits = species_order, labels = species_labels) +
  scale_fill_manual(values = species_colors, limits = species_order, labels = species_labels) +
  scale_shape_manual(values = species_shapes, limits = species_order, labels = species_labels)+
  
  theme_classic(base_size = 12, base_family = "sans") +
  theme(
    axis.title      = element_text(face = "bold", size = 13, color = "black"),
    axis.text       = element_text(size = 11,  color = "black"),
    axis.text.x     = element_text(face = "italic", vjust = 1, hjust = 1, angle = 30,size = 11, color = "black"), # Added italic for species names
    axis.line       = element_line(linewidth = 0.6, color = "black"),
    axis.ticks      = element_line(linewidth = 0.6, color = "black"),
    legend.position = "none",
    # legend.title    = element_text(face = "bold", size = 11),
    # legend.text     = element_text(size = 11, face = "italic"),
    plot.margin     = margin(t = 10, r = 10, b = 5, l = 10, unit = "pt")
  );plot_k3spp

saveRDS(plot_k3spp, "Outputs/plot_ed50K_spp.rds")

# ==============================================================================
# 3. For Cluster: Extract EMMs & Compact Letter Display (CLD)
# ==============================================================================
k3_emmeans_obj <- emmeans(k3ms, pairwise ~ K_3, adjust = "tukey")
k3_emm_df       <- as.data.frame(k3_emmeans_obj$emmeans)
k3_contrasts_df <- as.data.frame(k3_emmeans_obj$contrasts)

cat("\n--- Estimated Marginal Means (ED50) for K_3 Clusters ---\n")
print(k3_emm_df %>% dplyr::select(K_3, emmean, SE, df, lower.CL, upper.CL))

cat("\n--- Pairwise Tukey Contrasts Between K_3 Clusters ---\n")
print(k3_contrasts_df %>% dplyr::select(contrast, estimate, SE, p.value))

# Compute Tukey CLD Letters
cld_clusters   <- cld(k3_emmeans_obj$emmeans, adjust = "tukey", Letters = letters) %>%
  as.data.frame() %>%
  mutate(.group = trimws(.group))

# # ==============================================================================
# # 3. Align Cluster Labels & Factor Ordering
# # ==============================================================================
df_modeling <- df_modeling %>%
  mutate(
    Cluster_Label = roman_map[as.character(K_3)],
    Cluster_Label = factor(Cluster_Label, levels = target_levels)
  )

cld_clusters <- cld_clusters %>%
  mutate(
    Cluster_Label = roman_map[as.character(K_3)],
    Cluster_Label = factor(Cluster_Label, levels = target_levels)
  )

# Verify counts in console
cat("Raw Data Counts per Cluster Label:\n")
print(table(df_modeling$Cluster_Label, useNA = "always"))

cat("\nEMM Summary Counts per Cluster Label:\n")
print(table(cld_clusters$Cluster_Label, useNA = "always"))
# ==============================================================================
# 4. Publication-Ready Plot
# ==============================================================================
y_min_val <- min(df_modeling$ED50, na.rm = TRUE) - 0.5
y_max_val <- max(cld_clusters$upper.CL, na.rm = TRUE) + 0.8

plot_k3msr <- ggplot() +
  # Layer A: Raw data points
  # geom_jitter(
  #   data = df_modeling, 
  #   aes(x = Cluster_Label, y = ED50, color = Species),
  #   width = 0.12, 
  #   alpha = 0.8, 
  #   size = 3
  # ) +
  # Layer B: Model Estimated Marginal Means & 95% CIs
  geom_errorbar(
    data = cld_clusters,
    aes(x = Cluster_Label, ymin = lower.CL, ymax = upper.CL),
    width = 0.15, 
    linewidth = 0.8, 
    color = "black"
  ) +
  geom_point(
    data = cld_clusters,
    aes(x = Cluster_Label, y = emmean, fill = Cluster_Label),
    size = 5,
    shape=21
    # fill = "black",
    # color = "black"
  ) +
  # Layer C: Significance Letters (CLD)
  geom_text(
    data = cld_clusters,
    aes(x = Cluster_Label, y = upper.CL + 0.5, label = .group),
    #fontface = "bold", 
    hjust = 0.5,
    size = 6,
    color = "black",
    inherit.aes = FALSE
  ) +
  # Custom Aesthetics & Formatting
  #scale_color_manual(values = species_colors, limits = species_order, labels = species_labels) +
  # scale_fill_manual(values = species_colors, limits = species_order, labels = species_labels) +
  # scale_color_manual(values = cluster_colors, limits = roman_map) +
  #scale_fill_manual(values = cluster_colors, limits = roman_map) +
  # Discrete Fill Scale using target_levels for limits
  scale_color_manual(name = "Cluster", values = cluster_colors, limits = target_levels) +
  scale_fill_manual(name = "Cluster", values = cluster_colors, limits = target_levels) +
  #scale_shape_manual(values = custom_shapes)+
  coord_cartesian(ylim = c(y_min_val, y_max_val)) +
  labs(
    x = "Photophysiological Cluster (K = 3)",
    y = expression("ED50" ~ "(°C)"),
    color = "Cluster",
    #shape = "Site"
  ) +
  scale_x_discrete(
    limits = target_levels,
    labels = c("Cluster I" = "I", "Cluster II" = "II", "Cluster III" = "III", "Cluster IV" = "IV")
  ) +
  #ggrepel::geom_text_repel()+
  theme_classic(base_size = 12, base_family = "sans") +
  theme(
    axis.title      = element_text(face = "bold", size = 13, color = "black"),
    axis.text       = element_text(size = 11, color = "black"),
    axis.text.x     = element_text(face = "bold", size = 11, color = "black"),
    axis.line       = element_line(linewidth = 0.6, color = "black"),
    axis.ticks      = element_line(linewidth = 0.6, color = "black"),
    legend.position = "right",
    legend.title    = element_text(face = "bold", size = 11),
    legend.text     = element_text(size = 11),
    plot.margin = margin(t = 10, r = 10, b = 5, l = 10, unit = "pt")
  );plot_k3msr

saveRDS(plot_k3msr, "Outputs/plot_ed50K_3.rds")

# plot_k3msr <- ggplot() +
#   # Layer B: Model Estimated Marginal Means & 95% CIs
#   geom_errorbar(
#     data = cld_clusters,
#     aes(x = Cluster_Label, ymin = lower.CL, ymax = upper.CL),
#     width = 0.15, 
#     linewidth = 0.8, 
#     color = "black"
#   ) +
#   geom_point(
#     data = cld_clusters,
#     aes(x = Cluster_Label, y = emmean, fill = Cluster_Label, shape = species_shapes),
#     size = 5,
#     #shape = 21,       # Required so points accept 'fill' colors with a black outline
#     color = "black"
#   ) +
#   # Layer C: Significance Letters (CLD)
#   geom_text(
#     data = cld_clusters,
#     aes(x = Cluster_Label, y = upper.CL + 0.7, label = .group),
#     fontface = "bold", 
#     hjust = 0.5,
#     size = 6,
#     color = "black",
#     inherit.aes = FALSE
#   ) +
#   # Enforce Cluster Order on X-Axis
#   scale_x_discrete(limits = target_levels) +
#   
#   # Discrete Fill Scale using target_levels for limits
#   scale_fill_manual(
#     name = "Cluster",
#     values = cluster_colors, 
#     limits = target_levels
#   ) +
#   scale_x_discrete(
#     limits = target_levels,
#     labels = c("Cluster I" = "I", "Cluster II" = "II", "Cluster III" = "III")
#   ) +
#   
#   coord_cartesian(ylim = c(y_min_val, y_max_val)) +
#   labs(
#     x = "Photophysiological Cluster (K = 3)",
#     y = expression("ED50" ~ "(°C)")
#   ) +
#   theme_classic(base_size = 12, base_family = "sans") +
#   theme(
#     axis.title      = element_text(face = "bold", size = 13, color = "black"),
#     axis.text       = element_text(size = 11, color = "black"),
#     axis.text.x     = element_text( size = 11, color = "black"), # Non-italic for clusters
#     #axis.title.y     = element_text(face = "bold", size = 11, color = "black"),
#     axis.line       = element_line(linewidth = 0.6, color = "black"),
#     axis.ticks      = element_line(linewidth = 0.6, color = "black"),
#     legend.position = "right",
#     legend.title    = element_text(face = "bold", size = 11),
#     legend.text     = element_text(size = 11),
#     plot.margin     = margin(t = 10, r = 10, b = 5, l = 10, unit = "pt")
#   )
# 
# print(plot_k3msr)
# saveRDS(plot_k3msr, "Outputs/plot_ed50K_3.rds")

# ==============================================================================
# FACETED DOSE-RESPONSE TRAIT PLOTTING SCRIPT
# ==============================================================================
library(dplyr)
library(lme4)
library(lmerTest)
library(emmeans)
library(multcomp)
library(ggplot2)
library(purrr)
library(tidyr)
library(here)

# Load helper variables (species_colors, cluster_colors, target_levels, roman_map, etc.)
source(here("Scripts/00_visualization_prep.R"))
load(here("Outputs", "photophys_preprocessed_data.RData"))

# 1. Clean data & ensure factor definitions
df_modeling <- meta_filteredK %>%
  filter(!is.na(Site), !is.na(Species), !is.na(K_3)) %>%
  mutate(
    Species       = droplevels(as.factor(Species)),
    K_3           = droplevels(as.factor(K_3)),
    Site          = droplevels(as.factor(Site)),
    Cluster_Label = roman_map[as.character(K_3)],
    Cluster_Label = factor(Cluster_Label, levels = target_levels)
  )

# Define dose-response metrics to facet
dr_metrics <- c("ED5", "ED50", "ED95", "DW")

# Clean metric display labels for facets
metric_facet_labels <- c(
  "ED5"  = "ED[5]~('°C')",
  "ED50" = "ED[50]~('°C')",
  "ED95" = "ED[95]~('°C')",
  "DW"   = "Decline~Width~('°C')"
)

# ==============================================================================
# 2. ITERATE OVER METRICS: FIT MODELS, EXTRACT EMMEANS & CLDs
# ==============================================================================
all_cld_list <- list()

for (m in dr_metrics) {
  
  # Prepare temporary dataset for target metric
  df_temp <- df_modeling %>%
    filter(!is.na(.data[[m]])) %>%
    mutate(YY = .data[[m]])
  
  # Fit winning additive REML model
  m_reml <- lmer(YY ~ Species + K_3 + (1 | Site), data = df_temp, REML = TRUE)
  
  # Extract EMMs and calculate Tukey CLD letters across Cluster
  emm_obj <- emmeans(m_reml, pairwise ~ K_3, adjust = "tukey")
  
  cld_df <- cld(emm_obj$emmeans, adjust = "tukey", Letters = letters) %>%
    as.data.frame() %>%
    mutate(
      Metric        = m,
      Cluster_Label = roman_map[as.character(K_3)],
      Cluster_Label = factor(Cluster_Label, levels = target_levels),
      .group        = trimws(.group)
    )
  
  all_cld_list[[m]] <- cld_df
}

# Master data frame containing EMMs, CIs, and CLDs for all facets
master_cld_df <- bind_rows(all_cld_list) %>%
  mutate(Metric = factor(Metric, levels = dr_metrics))

# ==============================================================================
# 3. ALL EDS PLOT
# ==============================================================================
plot_faceted_dr <- ggplot(data = master_cld_df, aes(x = Cluster_Label, y = emmean)) +
  
  # Layer A: 95% Confidence Intervals
  geom_errorbar(
    aes(ymin = lower.CL, ymax = upper.CL),
    width = 0.18,
    linewidth = 0.7,
    color = "black"
  ) +
  
  # Layer B: Estimated Marginal Means (Points filled by Cluster)
  geom_point(
    aes(fill = Cluster_Label),
    size = 4,
    shape = 21,
    color = "black"
  ) +
  
  # Layer C: Tukey Significance Letters
  geom_text(
    aes(y = upper.CL + (max(upper.CL) - min(lower.CL)) * 0.08, label = .group),
    fontface = "bold",
    size = 4.5,
    color = "black",
    vjust = 0
  ) +
  
  # Facet Panel Layout (2x2 grid with free y-axes)
  facet_wrap(
    ~ Metric, 
    scales = "free_y", 
    ncol = 2,
    labeller = as_labeller(metric_facet_labels, default = label_parsed)
  ) +
  
  # Manual Aesthetics & Scales
  scale_fill_manual(
    name   = "Cluster",
    values = cluster_colors,
    limits = target_levels
  ) +
  scale_x_discrete(
    limits = target_levels,
    labels = c("Cluster I" = "I", "Cluster II" = "II", "Cluster III" = "III")
  ) +
  
  # Axes & Labels
  labs(
    x = "Photophysiological Cluster (K = 3)",
    y = "Estimated Marginal Mean (°C)"
  ) +
  
  # Publication-Ready Theme
  theme_classic(base_size = 12, base_family = "sans") +
  theme(
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text       = element_text(face = "bold", size = 12, color = "black"),
    axis.title       = element_text(face = "bold", size = 13, color = "black"),
    axis.text        = element_text(size = 11, color = "black"),
    axis.text.x      = element_text(face = "bold", size = 11, color = "black"),
    axis.line        = element_line(linewidth = 0.6, color = "black"),
    axis.ticks       = element_line(linewidth = 0.6, color = "black"),
    panel.spacing    = unit(1.2, "lines"),
    legend.position  = "right",
    legend.title     = element_text(face = "bold", size = 11),
    legend.text      = element_text(size = 11)
  )

# Display plot
print(plot_faceted_dr)

# Save high-resolution outputs
ggsave(
  filename = here("Outputs", "Figure_Faceted_DoseResponse_Clusters.png"),
  plot     = plot_faceted_dr,
  width    = 8.5,
  height   = 7.0,
  dpi      = 300
)

saveRDS(plot_faceted_dr, here("Outputs", "plot_faceted_dr.rds"))