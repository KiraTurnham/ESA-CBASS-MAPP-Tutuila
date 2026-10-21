rm(list = ls())
library(tidyverse)
library(vegan)
library(here)
library(ggplot2)
library(ggnewscale)
library(pheatmap)
library(viridis)
library(caret)
library(pvclust)
library(patchwork)
library(factoextra)
library(dplyr)
library(pairwiseAdonis)

source(here("Scripts/00_visualization_prep.R"))
load(here("Outputs", "photophys_preprocessed_data2.RData"))
metric_lookup<-read.csv(here("Data","MAPP_Metric_Lookup_Table_Compact.csv"))

##############################################################################
#### MULTIVARIATE ANALYSIS####################################################
##############################################################################
# =========================================================================
# PERMANOVA & Dispersion Tests
# =========================================================================
# Homogeneity test of Multivariate Dispersion (BETADISPER)
dist_filtered <- vegdist(permanova_matrix_filtered2, method = "euclidean")
bd_filtered   <- betadisper(dist_filtered, meta_filteredK2$Species)
bd_filtered_site  <- betadisper(dist_filtered, meta_filteredK2$Site)
anova_bd      <- anova(bd_filtered)
anova_bd_site <- anova(bd_filtered_site)

cat("\n--- FILTERED species DISPERSION TEST p-value:", anova_bd[1, "Pr(>F)"], "---\n")
#pass
cat("\n--- FILTERED site DISPERSION TEST p-value:", anova_bd_site[1, "Pr(>F)"], "---\n")
#does NOT pass

# # spp PERMANOVA
# permanova_filtered <- adonis2(
#   permanova_matrix_filtered ~ Species, 
#   data         = meta_filtered, 
#   method       = "euclidean", 
#   permutations = 999
# )
# cat("\n--- FILTERED PERMANOVA RESULTS ---\n")
# print(permanova_filtered)

meta_filtered2$Site_factor <- as.factor(meta_filtered2$Site)

# type 1 species and site global PERMANOVA
permanova_spp_filtered2 <- adonis2(
  permanova_matrix_filtered2 ~ Site + Species,
  data         = meta_filtered,
  method       = "euclidean",
  by = "terms", #site then spp
  permutations = 999
)
cat("\n--- FILTERED PERMANOVA RESULTS w spp  ---\n")
print(permanova_spp_filtered2)

#BH notes:  vegan can't do random effects AND the order of the factors matter when running by = "terms". So this model design asks:
#how much multivariate variation in the photophys data is explained first by Site?, and then
#how much additional variation is explained by Species, after Site has already be accounted for?
#in other words: "Does a species effect persist after accounting for site-to-site variation?"

# Run pairwise PERMANOVA between all species pairs
spp_pairwise <- pairwise.adonis2(
  permanova_matrix_filtered2 ~ Species, 
  data = meta_filtered2, 
  method = "euclidean",
  p.adjust.m = "bonferroni" # Adjusts p-values for multiple comparisons
)

# View the results
print(spp_pairwise)

# Run pairwise PERMANOVA between all site pairs
site_pairwise <- pairwise.adonis2(
  permanova_matrix_filtered2 ~ Site, 
  data = meta_filtered2, 
  method = "euclidean",
  p.adjust.m = "bonferroni" # Adjusts p-values for multiple comparisons
)

# View the results
print(spp_pairwise)
# =========================================================================
# 3. nMDS & Trait Vector Fitting
# =========================================================================
set.seed(123)
nmds_filtered <- metaMDS(
  permanova_matrix_filtered2, 
  distance      = "euclidean", 
  k             = 2, 
  trymax        = 200, 
  maxit         = 500, 
  autotransform = FALSE,
  trace         = 0
)

cat("\nFiltered nMDS Stress Value:", nmds_filtered$stress, "\n")

# Fit vectors for retained traits
#ef_filtered2 <- envfit(nmds_filtered, permanova_matrix_filtered2, permutations = 999)
#saveRDS(ef_filtered2, here("Outputs", "ef_filtered2.rds"))
ef_filtered2 <- readRDS(here("Outputs", "ef_filtered2.rds"))

vector_scores2 <- as.data.frame(scores(ef_filtered2, display = "vectors"))
vector_scores2$Metric <- rownames(vector_scores2)
vector_scores2$r2     <- ef_filtered2$vectors$r
vector_scores2$p_val  <- ef_filtered2$vectors$pvals

vector_scores2 <- vector_scores2 %>%
  left_join(metric_families, by = "Metric")

# Filter top metric drivers
top_vectors2 <- vector_scores2 %>%
  filter(p_val < 0.002) %>%
  arrange(desc(r2)) %>%
  head(6)

cat("\n--- TOP METRIC DRIVERS OF FILTERED ORDINATION ---\n")
print(top_vectors2)

nmds_scores2 <- as.data.frame(scores(nmds_filtered, display = "sites")) %>%
  mutate(SampleID_clean = rownames(permanova_matrix_filtered2)) %>%
  left_join(meta_filteredK2, by = "SampleID_clean")

max_site_coord   <- max(abs(c(nmds_scores2$NMDS1, nmds_scores2$NMDS2)))
max_vector_coord <- max(sqrt(top_vectors2$NMDS1^2 + top_vectors2$NMDS2^2))
arrow_mult       <- (max_site_coord * 0.8) / max_vector_coord

top_vectors_clean2 <- top_vectors2 %>%
  mutate(
    Metric_Category = Family, 
    Trait_Family    = str_extract(Metric, "(?i)(mQuant|qqP|mqP|ppq|mpq|qP|rqm|qm|rABQ|ABQ|Sigma|Tau1|Tau2|NPQ|nCon)"),
    Phase           = str_extract(Metric, "^[A-Za-z0-9]+"),
    Clean_Label     = case_when(
      !is.na(Trait_Family) & !is.na(Phase) ~ paste0(Trait_Family, " (", Phase, ")"),
      !is.na(Trait_Family)                 ~ Trait_Family,
      TRUE                                 ~ Metric
    ),
    NMDS1_scaled = NMDS1 * arrow_mult,
    NMDS2_scaled = NMDS2 * arrow_mult
  )
# Merge top_vectors with the lookup table
top_vectors_merged2 <- top_vectors2 %>%
  left_join(metric_lookup, by = c("Metric" = "Metric.Name"))

# View the result
print(top_vectors_merged2)
# 
write.csv(
  top_vectors_merged2,
  file = here("Outputs", "Table_Top_Vectors2.csv"),
  row.names = FALSE
)
# In 02_nmds_envfit.R:
saveRDS(top_vectors_merged2, file = here("Outputs", "top_vectors2.rds"))

############################################
#plot w no ellipses
# p <- ggplot() +
#   geom_point(
#     data = nmds_scores, 
#     aes(x = NMDS1, y = NMDS2, color = Species, shape = factor(Site)), 
#     size = 2.8, alpha = 0.75
#   ) + 
#   scale_shape_manual(values = custom_shapes, name = "Site") +
#   scale_color_manual(values = species_colors, name = "Species") +
#   scale_fill_manual(values = species_colors, name = "Species") +
#   new_scale_color() +
#   geom_segment(
#     data = top_vectors_clean, 
#     aes(x = 0, y = 0, xend = NMDS1_scaled, yend = NMDS2_scaled, color = Metric_Category),
#     arrow = arrow(length = unit(0.20, "cm")), linewidth = 0.85
#   ) +
#   geom_text(
#     data = top_vectors_clean,
#     aes(x = NMDS1_scaled * 1.10, y = NMDS2_scaled * 1.10, label = Metric),
#     size = 3
#   ) +
#   scale_color_manual(values = metric_colors, name = "Trait Family") +
#   coord_cartesian(clip = "off") +
#   theme_bw() +
#   theme(plot.margin = unit(c(15, 25, 15, 25), "pt"))+
#   labs( x = "nMDS Dimension 1", y = "nMDS Dimension 2")
# p

centroid_colors <- species_colors
centroid_colors["AGLO"] <-  "#5C8252FF"

# =========================================================================
# plot spp ellipses
# =========================================================================
#calc centroids
spp_centroids <- nmds_scores2 %>%
  filter(!is.na(Species)) %>%
  group_by(Species) %>%
  summarize(
    NMDS1 = mean(NMDS1, na.rm = TRUE),
    NMDS2 = mean(NMDS2, na.rm = TRUE),
    .groups = "drop"
  )
spp <- ggplot() +
  geom_point(
    data = nmds_scores2, 
    aes(x = NMDS1, y = NMDS2, color= Species, fill = Species, shape = Species), 
    size = 3, alpha = 0.7
  ) + 
  stat_ellipse(
    data = nmds_scores2, 
    aes(x = NMDS1, y = NMDS2, color = Species, fill = Species), 
    geom = "polygon", alpha = 0.3, level = 0.70
  ) +
  #centroid :
  geom_point(
    data = spp_centroids,
    aes(x = NMDS1, y = NMDS2, color = Species, fill = Species),
    shape = 4,        
    #color = "black",  
    size = 5,        
    stroke = 4,    
    show.legend = FALSE
  ) +
  scale_fill_manual(values = centroid_colors, limits = species_order) +
  #scale_shape_manual(values = custom_shapes, name = "Site") +
  scale_color_manual(
    name   = "Species",
    values = species_colors, 
    limits = species_order, 
    labels = species_labels
  ) +
  scale_fill_manual(
    name   = "Species",
    values = species_colors, 
    limits = species_order, 
    labels = species_labels
  ) +
  scale_shape_manual(
    name   = "Species",
    values = species_shapes, 
    limits = species_order, 
    labels = species_labels
  ) +
  theme_classic()+
  theme(
    legend.position = "right",
    legend.title    = element_text(face = "bold", size = 11),
    legend.text     = element_text(size = 11, face = "italic"))+
  new_scale_color() +
  
  # geom_segment(
  #   data = top_vectors_clean, 
  #   aes(x = 0, y = 0, xend = NMDS1_scaled, yend = NMDS2_scaled, color = Metric_Category),
  #   arrow = arrow(length = unit(0.20, "cm")), linewidth = 0.85
  # ) +
  # geom_text(
  #   data = top_vectors_clean,
  #   aes(x = NMDS1_scaled * 1.10, y = NMDS2_scaled * 1.10, label = Metric),
  #   size = 3
  # ) +
  #scale_color_manual(values = metric_colors, name = "Trait Family") +
  coord_cartesian(clip = "off") +
  theme(plot.margin = unit(c(15, 25, 15, 25), "pt"))+
  labs( x = "nMDS Dimension 1", y = "nMDS Dimension 2")

spp
# =========================================================================
# plot by K-3 ellipses
# =========================================================================
#prep labels
nmds_scores2 <- nmds_scores2 %>%
  mutate(
    Cluster_Label = roman_map[as.character(K_3)],
    Cluster_Label = factor(Cluster_Label, levels = target_levels),
    Site          = factor(as.character(Site), levels = c("1", "3", "4", "5", "7", "9", "10", "11"))
  )
#calc centroids
k3_centroids <- nmds_scores2 %>%
  filter(!is.na(Cluster_Label)) %>%
  group_by(Cluster_Label) %>%
  summarize(
    NMDS1 = mean(NMDS1, na.rm = TRUE),
    NMDS2 = mean(NMDS2, na.rm = TRUE),
    .groups = "drop"
  )

kplot <- ggplot() +
  # Points: Color by Species, Shape by Site
  geom_point(
    data = nmds_scores2, 
    aes(x = NMDS1, y = NMDS2, color = Cluster_Label, fill=Cluster_Label, shape = Species), 
    size = 3, 
    alpha = 0.5
  ) + 
  # scale_shape_manual(values = species_shapes, limits = species_order, labels = species_labels) +
  scale_color_manual(values = cluster_colors, limits = target_levels, name = "Photophysiological Cluster") +
  scale_fill_manual(values = cluster_colors, limits = target_levels, name = "Photophysiological Cluster") +
  scale_shape_manual(
    name   = "Species",
    values = species_shapes, 
    limits = species_order, 
    labels = species_labels
  ) +
  # 
  # # Scale layer for Ellipses (K_3 functional clusters)
  # new_scale_color() +
  # new_scale_fill() +
  # 
  # 95% Confidence Ellipses
  stat_ellipse(
    data = nmds_scores2 %>% filter(!is.na(Cluster_Label)),
    aes(x = NMDS1, y = NMDS2, color = Cluster_Label, fill = Cluster_Label),
    geom = "polygon", 
    alpha = 0.1, 
    level = 0.7
  ) +
  
  # Centroid Markers
  geom_point(
    data = k3_centroids,
    aes(x = NMDS1, y = NMDS2, color = Cluster_Label),
    shape = 4,        
    size = 5,          
    stroke = 4, 
    show.legend = FALSE
  ) +
  scale_color_manual(values = cluster_colors, limits = target_levels, name = "Photophysiological Cluster") +
  scale_fill_manual(values = cluster_colors, limits = target_levels, name = "Photophysiological Cluster") +
  # Scale layer for envfit vectors
  new_scale_color() +
  geom_segment(
    data = top_vectors_clean2, 
    aes(x = 0, y = 0, xend = NMDS1_scaled, yend = NMDS2_scaled,  color = Metric_Category),
    arrow = arrow(length = unit(0.20, "cm")), 
    linewidth = 1
  ) +
  # Replace geom_text() and ggrepel::geom_text_repel() with this single layer:
  # ggrepel::geom_text_repel(
  #   data = top_vectors_clean2,
  #   aes(x = NMDS1_scaled * 1.01, y = NMDS2_scaled * 1.01, label = Metric),
  #   size = 3,
  #   #fontface = "bold",
  #   box.padding = 0.5,
  #   point.padding = 0.3,
  #   segment.color = "grey50",  # Draws a subtle line connecting label to vector tip if pushed far
  #   segment.size = 0.3,
  #   max.overlaps = Inf        # Ensures all 6 labels get rendered
  # ) +
  scale_color_manual(values = metric_colors, name = "Trait Family", drop = FALSE) +
  # Formatting
  coord_cartesian(clip = "off") +
  theme_classic()+
  theme(
    legend.position = "none",
    legend.title    = element_text(face = "bold", size = 11),
    legend.text     = element_text(size = 11))+
  labs(
    # title = "Hclust Groupings (K = 3)", 
    x     = "nMDS Dimension 1", 
    y     = "nMDS Dimension 2"
  );kplot
##########
# combined_plot <- (spp + kplot) +
#   plot_layout(ncol = 2, guides = "collect") +
#   plot_annotation(
#     title = paste0("Filtered ( nMDS Stress = ", round(nmds_filtered$stress, 3), ")"),
#     tag_levels = 'A'
#   ) &
#   theme(legend.position = "right")
# 
# print(combined_plot)
# 
# ggsave(
#   filename = here("Plots", paste0("Fig4_nMDS_filtered_top10",".png")),
#   plot     = combined_plot,
#   width    = 16.0,
#   height   = 9.0,
#   dpi      = 300,
#   units    = "in"
# )
# 
# #all 3
# # combined_plot3 <- p + spp + kplot + 
# #   plot_layout(
# #     ncol   = 3, 
# #     guides = "collect" # Keeps specific legends per panel; use "collect" to combine shared legends
# #   ) + 
# #   plot_annotation(
# #     title    = "Multivariate Photophysiological Separation across Coral Species and Clusters",
# #     tag_levels = "A" # Automatically labels panels A, B, C
# #   )
# # 
# # print(combined_plot3)
# # 
# # ggsave(
# #   filename = here("Plots", paste0("nMDS_3combined",".png")),
# #   plot     = combined_plot3,
# #   width    = 16.0,
# #   height   = 9.0,
# #   dpi      = 300,
# #   units    = "in"
# # )
# 
# # # 1. Clear stuck graphics devices
# while(!is.null(dev.list())) dev.off()
# 
# # 2. Re-open a clean graphics device
# dev.new()
# 
# # 3. Print your plot object directly
# print(p)
