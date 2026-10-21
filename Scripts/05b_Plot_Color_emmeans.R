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
    YY      = color_slope,
    Species = droplevels(as.factor(Species)),
    Site    = droplevels(as.factor(Site)),
    across(starts_with("K_"), ~ droplevels(as.factor(.)))
  )

#model with site and species as random factor
k3ms <- lmer(YY ~  Species + (1 | Site), data = df_modeling, REML = TRUE)

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
y_min_val <- min(cld_spp$lower.CL, na.rm = TRUE) - 0.1
y_max_val <- max(cld_spp$upper.CL, na.rm = TRUE) + 0.1

plot_sppslope<- ggplot() +
  #Layer A: Raw data points
  # geom_jitter(
  #   data = df_modeling,
  #   aes(x = Species, y = color_slope, color = Species),
  #   width = 0.12,
  #   alpha = 0.5,
  #   size = 2
  # ) +
  #Layer B: Model Estimated Marginal Means & 95% CIs
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
    # fill = "black",
    # color = "white",
    #color="black",
    #shape =21
  ) +
  # Layer C: Significance Letters (CLD)
  geom_text(
    data = cld_spp,
    aes(x = Species, y = upper.CL + 0.1, label = .group),
   # fontface = "bold", 
    #hjust = 0.5,
    vjust = 1,
    size = 6,
    color = "black",
    inherit.aes = FALSE
  ) +
  # Custom Aesthetics & Formatting
  coord_cartesian(ylim = c(y_min_val, y_max_val)) +
  labs(
    y = expression("Color slope"),
    color = "Species"
  ) +
  
  # FORCES X-AXIS ORDER & ITALIC SPECIES LABELS
  scale_x_discrete(
    limits = species_order, 
    labels = species_labels
  ) +
  
  scale_color_manual(values = species_colors, limits = species_order, labels = species_labels) +
  scale_fill_manual(values = species_colors, limits = species_order, labels = species_labels) +
  scale_shape_manual(values = species_shapes)+
  
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
  );plot_sppslope
plot_sppslope

saveRDS(plot_sppslope, "Outputs/plot_colorK_spp.rds")
# # ==============================================================================
# # 3. For Cluster: Extract EMMs & Compact Letter Display (CLD)
# # ==============================================================================
# k3_emmeans_obj <- emmeans(k3ms, pairwise ~ K3, adjust = "tukey")
# k3_emm_df       <- as.data.frame(k3_emmeans_obj$emmeans)
# k3_contrasts_df <- as.data.frame(k3_emmeans_obj$contrasts)
# 
# cat("\n--- Estimated Marginal Means (ED50) for K_3 Clusters ---\n")
# print(k3_emm_df %>% dplyr::select(K_3, emmean, SE, df, lower.CL, upper.CL))
# 
# cat("\n--- Pairwise Tukey Contrasts Between K_3 Clusters ---\n")
# print(k3_contrasts_df %>% dplyr::select(contrast, estimate, SE, p.value))
# 
# # Compute Tukey CLD Letters
# cld_clusters   <- cld(k3_emmeans_obj$emmeans, adjust = "tukey", Letters = letters) %>%
#   as.data.frame() %>%
#   mutate(.group = trimws(.group))
# 
# # # ==============================================================================
# # # 3. Align Cluster Labels & Factor Ordering
# # # ==============================================================================
# df_modeling <- df_modeling %>%
#   mutate(
#     Cluster_Label = roman_map[as.character(K_3)],
#     Cluster_Label = factor(Cluster_Label, levels = target_levels)
#   )
# 
# cld_clusters <- cld_clusters %>%
#   mutate(
#     Cluster_Label = roman_map[as.character(K_3)],
#     Cluster_Label = factor(Cluster_Label, levels = target_levels)
#   )
# 
# # Verify counts in console
# cat("Raw Data Counts per Cluster Label:\n")
# print(table(df_modeling$Cluster_Label, useNA = "always"))
# 
# cat("\nEMM Summary Counts per Cluster Label:\n")
# print(table(cld_clusters$Cluster_Label, useNA = "always"))
# # ==============================================================================
# # 4. Publication-Ready Plot
# # ==============================================================================
# y_min_val <- min(df_modeling$color_slope, na.rm = TRUE) 
# y_max_val <- max(cld_clusters$upper.CL, na.rm = TRUE) + 0.1
# 
# plot_k3msr <- ggplot() +
#   # Layer A: Raw data points
#   geom_jitter(
#     data = df_modeling, 
#     aes(x = Cluster_Label, y = color_slope, color = Species),
#     width = 0.12, 
#     alpha = 0.5, 
#     size = 2
#   ) +
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
#     aes(x = Cluster_Label, y = emmean),
#     size = 5,
#     fill = "black",
#     color = "black"
#   ) +
#   # Layer C: Significance Letters (CLD)
#   geom_text(
#     data = cld_clusters,
#     aes(x = Cluster_Label, y = upper.CL + 0.5, label = .group),
#     fontface = "bold", 
#     hjust = 0.5,
#     size = 6,
#     color = "black",
#     inherit.aes = FALSE
#   ) +
#   # Custom Aesthetics & Formatting
#   scale_color_manual(values = species_colors, limits = species_order, labels = species_labels) +
#   scale_fill_manual(values = species_colors, limits = species_order, labels = species_labels) +
#   scale_shape_manual(values = custom_shapes)+
#   coord_cartesian(ylim = c(y_min_val, y_max_val)) +
#   labs(
#     x = "Photophysiological Cluster (K = 3)",
#     y = expression("Color slope"),
#     color = "Species",
#     #shape = "Site"
#   ) +
#   theme_classic(base_size = 12, base_family = "sans") +
#   theme(
#     axis.title      = element_text(face = "bold", size = 13, color = "black"),
#     axis.text       = element_text(size = 11, color = "black"),
#     axis.text.x     = element_text(face = "bold", size = 11, color = "black"),
#     axis.line       = element_line(linewidth = 0.6, color = "black"),
#     axis.ticks      = element_line(linewidth = 0.6, color = "black"),
#     legend.position = "right",
#     legend.title    = element_text(face = "bold", size = 11),
#     legend.text     = element_text(size = 11, face = "italic"),
#     plot.margin = margin(t = 10, r = 10, b = 5, l = 10, unit = "pt")
#   );plot_k3msr
# 
# plot_k3msrcolor <- ggplot() +
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
#     aes(x = Cluster_Label, y = emmean, fill = Cluster_Label),
#     size = 5,
#     shape = 21,       # Required so points accept 'fill' colors with a black outline
#     color = "black"
#   ) +
#   # Layer C: Significance Letters (CLD)
#   geom_text(
#     data = cld_clusters,
#     aes(x = Cluster_Label, y = upper.CL + 0.1, label = .group),
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
#     y = expression("Color slope")
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
# print(plot_k3msrcolor)
# saveRDS(plot_k3msrcolor, "Outputs/plot_colorK_3.rds")
