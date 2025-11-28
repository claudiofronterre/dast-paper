# =============================================================================
# Figures produced in the paper
# =============================================================================
# Code used to produce some of the auxiliary figures reported in the paper 
#
# Inputs (expected in ./data/)
#   - KenyaData.csv (survey data with x, y, year, subcounty, noASC, noTT, noHKW)
#   - KenyaCoverage_TotalMDA.shp (subcounty shapefile with MDA data)
#   - simulations.qs (results from the simulation study)
#
# Outputs (written to ./figures/)
#   - sth_kenya_panels.pdf
#   - sim_bias.pdf
#   - sim_RMSE.pdf
#
# =============================================================================

options(stringsAsFactors = FALSE)

# --------------------------- 0. Libraries -----------------------------------
suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(stringr)
  library(scales)
  library(ggspatial)
  library(rgeoboundaries)
})

if (!dir.exists("outputs")) dir.create("outputs")

# 1. IMPORT DATA ---------------------------------------------------------------

# Survey data (point locations in ESRI:102022 as provided)
sth <- read.csv("data/KenyaData.csv")
sth <- st_as_sf(sth, coords = c("x", "y"), crs = "ESRI:102022")

# Subcounty MDA polygons
subcounty <- st_read("data/KenyaCoverage_TotalMDA.shp", quiet = TRUE)

# Kenya boudaries
kenya_bg <- gb_adm0("Kenya")

# --------------------------- 2. Data processing -------------------------------


# Change CRS
sth <- st_transform(sth, crs = 4326)
subcounty <- st_transform(subcounty, crs = 4326)


# Define time panels
panel_labs <- c("2012 (Baseline)", "2013–2014", "2015–2016", "2017–2018", "2021–2022")

sth_pan <- sth |>
  mutate(
    panel = case_when(
      year == 2012 ~ "2012 (Baseline)",
      year %in% 2013:2014 ~ "2013–2014",
      year %in% 2015:2016 ~ "2015–2016",
      year %in% 2017:2018 ~ "2017–2018",
      year %in% 2021:2022 ~ "2021–2022",
      .default = NA_character_
    ),
    panel = factor(panel, levels = panel_labs)
  )

# --------------------- 3. Simulation results (Fig. 4 and Fig. 5) --------------


# --------------------- 4. Prevalence map for STH (Fig. 6) ---------------------

base_theme <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2),
    axis.title = element_blank(),
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    legend.position = "top"
  )

# A function to plot a single year-range panel
plot_panel <- function(lbl) {
  pdat <- filter(sth_pan, panel == lbl)
  ggplot() +
    geom_sf(data = kenya_bg, fill = "grey95", color = "black") +
    geom_sf(
      data = pdat,
      aes(fill = prevSTH),
      stroke = 0.2, shape = 21, size = 1.1, alpha = .95
    ) +
    scale_fill_distiller(
      type = "div",
      palette = "RdYlBu",
      direction = -1,
      limits = c(0, 0.75),
      breaks = seq(0, 0.75, by = 0.1),
      labels = scales::label_percent(accuracy = 1),
      name = "Prevalence (%)",
      guide = guide_colorbar(
        title.position = "top",        
        title.hjust = 0.5,             
        barwidth = unit(6, "in"),      
        barheight = unit(0.4, "cm")
      )) +
    coord_sf() +
    labs(title = lbl) +
    base_theme
}

# Build each panel
p_2012 <- plot_panel("2012 (Baseline)")
p_1314 <- plot_panel("2013–2014")
p_1516 <- plot_panel("2015–2016")
p_1718 <- plot_panel("2017–2018")
p_2122 <- plot_panel("2021–2022")

# Subcounty-only panel
p_subc <- ggplot() +
  geom_sf(data = kenya_bg, fill = "grey95", color = NA) +
  geom_sf(data = subcounty, fill = "grey", color = "grey40", linewidth = 0.15) +
  labs(title = "Subcounties") +
  base_theme +
  theme(legend.position = "none") +
  annotation_scale(location = "bl", width_hint = 0.25) +  
  annotation_north_arrow(location = "tr", which_north = "true", 
                         height = unit(.5, "cm"), width = unit(.5, "cm")) 

# Assemble 3×2 layout and add title/caption
plt <- (guide_area() / (
  p_2012 | p_1314 | p_1516
) /
  (p_1718 | p_2122 | p_subc)) +
  plot_layout(guides = "collect", heights = unit(c(1, 1, 1), c("cm", "null", "null"))) &
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.box.just = "center",
    legend.title = element_text(size = 11, face = "bold", vjust = 0.3),
    legend.text = element_text(size = 9),
  )

plt
ggsave(filename = "figures/sth_kenya_panels.pdf",
       plot = plt, device = cairo_pdf, 
       bg = "white", width = 8, height = 7.5)