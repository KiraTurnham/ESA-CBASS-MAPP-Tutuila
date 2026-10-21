library(ggplot2)
library(here)
library(dplyr)
library(sf)
library(readr)

site_title    <- "Site"
species_title <- "Species"

# 1. Load data
load(here("Outputs", "photophys_preprocessed_data.RData"))

# 2. Read new site density metadata sheet and extract coordinates per original Site
geo_metadata <- read.csv(here("Data", "fixed_sites_sdmtmb_predicted_density_2025.csv")) %>%
  mutate(
    # Extract numeric site code if values look like "ESA-TUT-011" -> "11"
    Site_orig = gsub(".*-", "", as.character(Site)),
    Site_orig = as.character(as.numeric(Site_orig)) # strips leading zeros
  ) %>%
  group_by(Site_orig) %>%
  summarise(
    Latitude  = first(Lat),
    Longitude = first(Long),
    .groups   = "drop"
  )

# 3. Explicitly define original site mapping (Change right-side values if needed)
site_rename_map <- c(
  "1"  = "Site 6",
  "3"  = "Site 1",
  "4"  = "Site 3",
  "5"  = "Site 2",
  "7"  = "Site 8",
  "9"  = "Site 4",
  "10" = "Site 7",
  "11" = "Site 5"
)

# 4. Process meta_filteredK, apply new site names, and attach coordinates
meta_filteredK <- meta_filteredK %>%
  mutate(
    Site_orig = as.character(Site),
    Site      = case_match(
      Site_orig,
      "1"  ~ "Site 6",
      "3"  ~ "Site 1",
      "4"  ~ "Site 3",
      "5"  ~ "Site 2",
      "7"  ~ "Site 8",
      "9"  ~ "Site 4",
      "10" ~ "Site 7",
      "11" ~ "Site 5",
      .default = Site_orig
    )
  ) %>%
  left_join(geo_metadata, by = "Site_orig") %>%
  dplyr::select(-Site_orig)

theme_coral_pub <- function(base_size = 12) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      legend.position = "right",
      legend.title    = element_text(face = "bold", size = 11),
      legend.text     = element_text(size = 10, face = "italic"),
      axis.title      = element_text(face = "bold", size = 12, color = "black"),
      axis.text.y     = element_text(size = 11, color = "black"),
      axis.text.x     = element_text(size = 11, color = "black", face = "bold"),
      axis.line       = element_line(linewidth = 0.6, color = "black"),
      axis.ticks      = element_line(linewidth = 0.6, color = "black"),
      plot.margin     = ggplot2::margin(t = 8, r = 8, b = 5, l = 8, unit = "pt") # Namespaced fix
    )
}
# --- Aesthetic Mappings & Palettes ---

custom_shapes <- c(23, 22, 15, 16, 17, 18, 21, 24)
site_shapes   <- c(23, 22, 15, 16, 17, 18, 21, 24)
all_sites     <- paste("Site", 1:8)

species_colors <- c(
  "AABR" = "#046E8F",
  "AGLO" = "#83BA75FF",
  "AHYA" = "#D44D5C",
  "ICRA" = "#462255"
)

species_order  <- c("ICRA", "AGLO", "AABR", "AHYA")

species_shapes <- c(
  "ICRA" = 23,
  "AGLO" = 24,
  "AABR" = 22,
  "AHYA" = 25
)

species_labels <- c(
  "AGLO" = "A. globiceps",
  "AABR" = "A. abrotanoides",
  "AHYA" = "A. hyacinthus",
  "ICRA" = "I. crateriformis"
)

# Cluster & Interaction Mappings
K_names <- c("2" = "Group I", "1" = "Group II", "3" = "Group III")
K_names_order <- c("Group I", "Group II", "Group III")

roman_map <- c(
  "Cluster_1" = "Cluster II", "Cluster_2" = "Cluster I", "Cluster_3" = "Cluster III",
  "1" = "Cluster II", "2" = "Cluster I", "3" = "Cluster III"
)

target_levels <- c("Cluster I", "Cluster II", "Cluster III")

cluster_colors <- c(
  "Cluster I"   = "violet",
  "Cluster II"  = "orange",
  "Cluster III" = "cyan"
)

# Metric & Annotation Palettes
metric_colors <- c(
  "Quant (Quantum Yield)"             = "green",
  "Sigma (Antenna Size)"              = "#FCF340",
  "Connect (Connectivity)"            = "#1F77B4", 
  "Connect (Normalized Curves)"       = "blue", 
  "Tau1 (Transport Kinetics 1)"       = "#17BECF", 
  "Tau2 (Transport Kinetics 2)"       = "#9467BD", 
  "NPQ (Non-Photochemical Quenching)" = "#8C564B",
  "qP (Photochemical Quenching)"      = "#E377C2", 
  "ABQ (Antenna Bed Quenching)"       = "#FF7F0E",
  "qm (Max Quenching)"                = "#FF13F0", 
  "OJIP Kinetics (Derivatives)"       = "#98DF8A", 
  "OJIP Ratios & Area"                = "#FFBB78", 
  "Other Metric"                      = "#7F7F7F"  
)

site_colors <- c(
  "Site 1" = "#E41A1C", "Site 2" = "#FF7F00", "Site 3" = "#FFD700", "Site 4" = "#4DAF4A", 
  "Site 5" = "green",   "Site 6" = "#377EB8", "Site 7" = "#4B0082", "Site 8" = "#984EA3"
)

ann_colors <- list(
  Family  = metric_colors,
  Species = species_colors,
  K_3     = cluster_colors,
  Site    = site_colors
)

theme_coral_pub <- function(base_size = 12) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      legend.position = "right",
      legend.title    = element_text(face = "bold", size = 11),
      legend.text     = element_text(size = 10, face = "italic"),
      axis.title      = element_text(face = "bold", size = 12, color = "black"),
      axis.text.y     = element_text(size = 11, color = "black"),
      axis.text.x     = element_text(size = 11, color = "black", face = "bold"),
      axis.line       = element_line(linewidth = 0.6, color = "black"),
      axis.ticks      = element_line(linewidth = 0.6, color = "black"),
      plot.margin     = ggplot2::margin(t = 8, r = 8, b = 5, l = 8, unit = "pt") # Namespaced fix
    )
}

