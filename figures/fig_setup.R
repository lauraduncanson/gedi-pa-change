# Shared setup for the figure notebooks: paths, packages, data loading, country map, theme.
# Edit the paths here; every notebook sources this file.

INPUTS  <- Sys.getenv("FIG_INPUTS",  "/home/jovyan/GEDI/results/figure_inputs")   # totals.py output
DRIFT   <- Sys.getenv("FIG_DRIFT",   "/home/jovyan/GEDI/results")                  # drift tables
OUTDIR  <- Sys.getenv("FIG_OUTDIR",  "/home/jovyan/GEDI/results/figures")
CACHE   <- Sys.getenv("FIG_CACHE",   "/home/jovyan/GEDI/results/figures/cache")
MAIN    <- "agbd_c"     # headline metric: drift-corrected AGBD (central)
MIN_PCT_COUNTERFACTUAL <- 10   # countries below this % of PAs with own counterfactuals are flagged

pkgs <- c("dplyr", "tidyr", "readr", "ggplot2", "sf", "patchwork", "scales", "forcats")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2); library(sf)
  library(patchwork); library(scales); library(forcats)
})
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)

WWF_BIOMES <- c(`1` = "Tropical moist broadleaf", `2` = "Tropical dry broadleaf",
                `3` = "Tropical conifer", `4` = "Temperate broadleaf & mixed",
                `5` = "Temperate conifer", `6` = "Boreal", `7` = "Tropical grassland & savanna",
                `8` = "Temperate grassland & savanna", `9` = "Flooded grassland",
                `10` = "Montane grassland & shrubland", `11` = "Tundra",
                `12` = "Mediterranean", `13` = "Desert & xeric shrubland", `14` = "Mangroves")
biome_name <- function(b) unname(ifelse(is.na(b), NA, WWF_BIOMES[as.character(b)]))
METRIC_LABEL <- c(agbd = "Uncorrected", agbd_c = "Drift-corrected",
                  agbd_c_lo = "Corrected (low drift)", agbd_c_hi = "Corrected (high drift)")

rd <- function(name) read_csv(file.path(INPUTS, name), show_col_types = FALSE)
load_inputs <- function() {
  list(country = rd("summary_country.csv"), continent = rd("summary_continent.csv"),
       biome = rd("summary_biome.csv") %>% mutate(biome_name = biome_name(biome)),
       cb = rd("summary_continent_biome.csv") %>% mutate(biome_name = biome_name(biome)),
       global = rd("summary_global.csv"), pa = rd("pa_table.csv"),
       boot_country = rd("boot_country.csv"), boot_global = rd("boot_global.csv"))
}

# Country polygons: rnaturalearth if installed, otherwise Natural Earth 1:50m GeoJSON (cached).
get_countries <- function() {
  if (requireNamespace("rnaturalearth", quietly = TRUE) &&
      requireNamespace("rnaturalearthdata", quietly = TRUE)) {
    w <- rnaturalearth::ne_countries(scale = 50, returnclass = "sf")
    w$iso3 <- ifelse(w$iso_a3_eh != "-99", w$iso_a3_eh, w$adm0_a3)
  } else {
    f <- file.path(CACHE, "ne_50m_admin_0_countries.geojson")
    if (!file.exists(f)) download.file(
      "https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/ne_50m_admin_0_countries.geojson",
      f, quiet = TRUE)
    w <- st_read(f, quiet = TRUE)
    w$iso3 <- ifelse(w$ISO_A3_EH != "-99", w$ISO_A3_EH, w$ADM0_A3)
  }
  w <- w[w$iso3 != "ATA", ]
  sf_use_s2(FALSE)
  w <- st_wrap_dateline(st_make_valid(st_transform(w["iso3"], 4326)),
                        options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"), quiet = TRUE)
  st_transform(st_make_valid(w), "+proj=robin")
}

# Outline of the globe in the Robinson projection, drawn behind the countries.
globe_outline <- function() {
  lon <- c(seq(-180, 180, 1), rep(180, 179), seq(180, -180, -1), rep(-180, 179))
  lat <- c(rep(-90, 361), seq(-89, 89, 1), rep(90, 361), seq(89, -89, -1))
  st_transform(st_sfc(st_polygon(list(cbind(c(lon, lon[1]), c(lat, lat[1])))), crs = 4326), "+proj=robin")
}

theme_fig <- function(base = 11) {
  theme_minimal(base_size = base) +
    theme(panel.grid.minor = element_blank(), plot.title.position = "plot",
          plot.title = element_text(face = "bold", size = base + 1),
          legend.title = element_text(size = base - 1), strip.text = element_text(face = "bold"))
}
theme_map <- function(base = 11) {
  theme_void(base_size = base) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold", size = base + 1),
          legend.key.width = unit(1.6, "cm"), legend.key.height = unit(0.3, "cm"),
          plot.title.position = "plot")
}
save_fig <- function(p, name, w = 10, h = 6) {
  ggsave(file.path(OUTDIR, paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path(OUTDIR, paste0(name, ".pdf")), p, width = w, height = h, bg = "white")
  invisible(p)
}
