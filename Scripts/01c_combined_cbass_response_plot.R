library(ggplot2)
library(patchwork)
library(here)

# Load shared theme/scales
source(here("Scripts/00_visualization_prep.R"))

# Load pre-made plots
ed50_emm  <- readRDS(here("Outputs", "plot_ed50.rds"))
color_emm <- readRDS(here("Outputs", "plot_color.rds"))

#removing spp from ed50 since share
ed50_emm <- ed50_emm +  
  theme(
    axis.text.x  = element_blank(),
    axis.title.x = element_blank()
    #axis.ticks.x = element_blank()
  )
# ------------------------------------------------------------------------------
# Combine with Patchwork
# ------------------------------------------------------------------------------
# '+' places plots side-by-side; '/' places them stacked vertically
combined_fig <- (ed50_emm / color_emm) + 
  scale_color_manual(
    values = species_colors, 
    limits = species_order, 
    labels = species_labels
  ) +
  
  # Applies full names to x-axis tick labels
  scale_x_discrete(
    limits = species_order, 
    labels = species_labels
  )+
  plot_layout(
    guides = "collect", # Collects shared species/site legends into a single legend
    axis_titles = "collect") + 
  plot_annotation(
    tag_levels = "A" # Adds publication tags (A, B)
  ) & 
  theme(legend.position = "right")

print(combined_fig)



# ------------------------------------------------------------------------------
# Export Multi-Panel Figure
# ------------------------------------------------------------------------------
# For a 2-panel vertical plot (ED50 on top, Color on bottom):
ggsave(
  filename = here("Plots", "Fig2_ED50_and_Color_Combined.pdf"),
  plot     = combined_fig,
  width    = 6.5,       # Full page width (standard two-column journal width)
  height   = 8.5,       # Proportional height for a 2-panel vertical stack
  units    = "in",
  device   = cairo_pdf
)

ggsave(
  filename = here("Plots", "Fig2_ED50_and_Color_Combined.png"),
  plot     = combined_fig,
  width    = 6.5,
  height   = 8.5,
  units    = "in",
  dpi      = 600
)

#horizontal plot
# 1. Apply discrete scales individually to ensure both plots render identical species axes & legends
library(ggplot2)
library(patchwork)

# 1. Prepare Base Plots
p1 <- ed50_emm +
  theme_coral_pub() +
  scale_x_discrete(limits = species_order, labels = species_labels) +
  #scale_color_manual(name = species_title, values = species_colors, limits = species_order, labels = species_labels) +
  #scale_shape_manual(name = site_title, values = site_shapes)
  aes(shape = NULL, species = NULL)

p2 <- color_emm +
  theme_coral_pub() +
  scale_x_discrete(limits = species_order, labels = species_labels) +
 # scale_color_manual(name = species_title, values = species_colors, limits = species_order, labels = species_labels) +
 # scale_shape_manual(name = site_title, values = site_shapes)
  aes(shape = NULL, species = NULL)

# 2. Assemble and Enforce Formatting Globally via Patchwork '&'
combined_fig <- (p1 + p2) +
  plot_layout(
    guides = "collect",
    axis_titles = "collect"
  ) +
  plot_annotation(tag_levels = "A") &
  # Global guides override: Drop species legend, collapse site to 1 row
 # guides(
    #color = "none",
    #shape = guide_legend(title = site_title, nrow = 1)
  #) &
  # Global theme override: Force horizontal top legend & straight italic X-axis text
  theme(
    # legend.position = "top",
    # legend.direction = "horizontal",
    # legend.box = "horizontal",
    axis.text.x = element_text(angle = 30, hjust = 1, vjust=1, face = "italic", color = "black")
  )

print(combined_fig)
####################################################################################################
ed50K_emm  <- readRDS(here("Outputs", "plot_ed50K_3.rds"))
colorK_emm <- readRDS(here("Outputs", "plot_colorK_3.rds"))
ed50spp_emm  <- readRDS(here("Outputs", "plot_ed50K_spp.rds"))
colorspp_emm <- readRDS(here("Outputs", "plot_colorK_spp.rds"))

ed50spp_emm <- ed50spp_emm +  
  theme(
    #axis.text.x  = element_blank(),
   # axis.title.x = element_blank(),
    #axis.ticks.x = element_blank(),
    legend.position = "none")+
  scale_x_discrete(limits = species_order, labels = species_labels) +
  aes(shape = Species) 
colorspp_emm <- colorspp_emm +  
  theme(
   # axis.text.x  = element_blank(),
   # axis.title.x = element_blank(),
    #axis.ticks.x = element_blank(),
    legend.position = "none")+
  scale_x_discrete(limits = species_order, labels = species_labels) +
  aes(shape = Species) 
  

combined_fig1 <- (ed50spp_emm + colorspp_emm) + 
  # scale_color_manual(
  #   values = species_colors, 
  #   limits = species_order, 
  #   labels = species_labels
  # ) +
  # 
  # Applies full names to x-axis tick labels
  scale_x_discrete(
    limits = species_order, 
    labels = species_labels
  )+
  plot_layout(
    guides = "collect") + # , # Collects shared species/site legends into a single legend
    #axis_titles = "collect") + 
  plot_annotation(tag_levels = list(c("B", "C"))) &
  
  # plot_annotation(
  #   tag_levels = "A" # Adds publication tags (A, B)
  # ) & 
  theme(legend.position = "none", 
        legend.text = element_text(size = 10, face = "italic"),
        axis.title.y  = element_text(size = 11, face = "bold", color = "black"))

print(combined_fig1)

ed50K_emm <- ed50K_emm +  
  theme(
    axis.text.x  = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "right")
    #axis.ticks.x = element_blank()
  

combined_fig2 <- (ed50K_emm / colorK_emm) + 
  # scale_color_manual(
  #   values = species_colors, 
  #   limits = species_order, 
  #   labels = species_labels
  # ) +
  # 
  # Applies full names to x-axis tick labels
  scale_x_discrete(
    limits = target_levels,
    labels = c("Cluster I" = "I", "Cluster II" = "II", "Cluster III" = "III")
  ) +
  plot_layout(
    guides = "collect", # Collects shared species/site legends into a single legend
    axis_titles = "collect") + 
  plot_annotation(tag_levels = list(c("B", "C"))) &
  
  theme(legend.position = "none")

print(combined_fig2)
  ####################################################################
#4 panel graph

colorspp_emm <- colorspp_emm +
  scale_x_discrete(
    limits = species_order, 
    labels = species_labels
  )

colorK_emm <- colorK_emm +
  scale_x_discrete(
    limits = target_levels,
    labels = c("Cluster I" = "I", "Cluster II" = "II", "Cluster III" = "III")
  )

#4 
# 1. Format Top Row (Species: B = ED50, C = Color)
p_B <- ed50spp_emm + 
  scale_x_discrete(limits = species_order, labels = species_labels) +
  theme(
    axis.text.x  = element_text(angle = 35, hjust = 1, vjust =1, face = "italic", color = "black"),
    axis.title.x = element_text(face = "bold", size = 11, color = "black")
  )

p_C <- colorspp_emm + 
  scale_x_discrete(limits = species_order, labels = species_labels) +
  theme(
    axis.text.x  = element_text(angle = 35, hjust = 1, face = "italic", color = "black"),
    axis.title.x = element_text(face = "bold", size = 11, color = "black")
  )

# 2. Format Row 2 (Clusters: D = ED50, E = Color)
p_D <- ed50K_emm + 
  scale_x_discrete(limits = target_levels, labels = c("Cluster I" = "I", "Cluster II" = "II", "Cluster III" = "III")) +
  theme(
    axis.text.x  = element_text(angle = 0, hjust = 0.5, face = "plain", color = "black"),
    axis.title.x = element_text(face = "bold", size = 11, color = "black")
  )

p_E <- colorK_emm + 
  scale_x_discrete(limits = target_levels, labels = c("Cluster I" = "I", "Cluster II" = "II", "Cluster III" = "III")) +
  theme(
    axis.text.x  = element_text(angle = 0, hjust = 0.5, face = "plain", color = "black"),
    axis.title.x = element_text(face = "bold", size = 11, color = "black")
  )

# 3. Assemble 2x2 Grid (Species on Row 1, Clusters on Row 2)
four_panel_fig <- wrap_plots(p_B, p_C, p_D, p_E, ncol = 2) +
  plot_layout(
    guides = "collect",
    axis_titles = "keep"
  ) +
  plot_annotation(tag_levels = list(c("B", "C", "D", "E"))) &
  
  guides(
    color = guide_legend(title = "Species", nrow = 2, byrow = TRUE),
    fill  = "none"
  ) &
  
  theme(
    legend.position   = "top",
    legend.box        = "vertical",
    legend.title      = element_text(face = "bold", size = 11),
    legend.text       = element_text(size = 10, face = "italic"),
    plot.tag          = element_text(face = "bold", size = 14),
    axis.title        = element_text(face = "bold", size = 11, color = "black"),
    axis.text.y       = element_text(size = 10, color = "black"),
    axis.line         = element_line(linewidth = 0.6, color = "black"),
    plot.margin       = margin(t = 5, r = 10, b = 10, l = 5, unit = "pt")
  );four_panel_fig


