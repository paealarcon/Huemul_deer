# =============================================================================
# Movement behavior and habitat selection to guide conservation of an
# endangered Neotropical deer
#
# Alarcón, P.A.E., Crego, R.D., Galetto, E., Aprile, G.,
# Morales, J.M., Pauli, J.N. & Donadio, E.
#
# -----------------------------------------------------------------------------
# Script:  03d_suitabilityMaps.R
# Purpose: Predict continuous and 5-class habitat-suitability rasters across
#          the study area from the fitted forest-intercept hRSF
# =============================================================================

# 0. Libraries ----
suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(INLA)
})

# 1. Configuration ----
FIT_RDA       <- "./hrsf/outputs/hrsf_forest_intercept.rda"
RESULTS_RDS   <- "./hrsf/outputs/hrsf_forest_intercept_results.rds"
STATIC_STACK  <- "./data/layers/stack_static.tif"
LC_TIF        <- "./data/layers/others/LCChalten_Pred.tif"

# Merged protected areas polygon (Los Glaciares NP + Reserva Natural Piedra
# del Fraile + Reserva Provincial Lago del Desierto, dissolved into ONE
# multipolygon for cleaner display).
PA_MERGED_SHP <- "./data/layers/others/protectedAreas.shp"
ARG_CH_SHP    <- "./data/layers/others/arg-ch border.shp"

# Vector overlay layers
RIVERS_LAKES_SHP <- "./data/layers/waterCourses/rivers_and_lakes.shp"
ROAD_SHP         <- "./data/layers/others/ruta_23_line.shp"

OUT_DIR_TIF   <- "./hrsf/outputs"
OUT_DIR_FIG   <- "./hrsf/figures"
dir.create(OUT_DIR_TIF, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIR_FIG, recursive = TRUE, showWarnings = FALSE)

CONT_TIF  <- file.path(OUT_DIR_TIF, "group_hRSF_continuous.tif")
CLASS_TIF <- file.path(OUT_DIR_TIF, "group_hRSF_5class.tif")
OUT_JPG   <- file.path(OUT_DIR_FIG, "FigureSD6.jpg")

# IWLR weight used in 03a. Adding log(W_LARGE) to muForest converts the
# weighted-binomial intercept back to the natural log-odds scale.
W_LARGE <- 1000

TARGET_CRS_UTM <- "EPSG:32719"

stopifnot(file.exists(FIT_RDA), file.exists(RESULTS_RDS),
          file.exists(STATIC_STACK), file.exists(LC_TIF))

# 2. Posterior summaries from the INLA fit ----
cat("\n=== Loading hRSF INLA fit + tidied results ===\n")
load(FIT_RDA)                          # provides `hrsf_forest`
results <- readRDS(RESULTS_RDS)

fe          <- hrsf_forest$summary.fixed
indiv_betas <- results$individual_betas
scale_params <- results$scale_params

cat("  fixed-effect rows: ", paste(rownames(fe), collapse = ", "), "\n")

beta_of <- function(name) {
  if (name %in% rownames(fe)) as.numeric(fe[name, "mean"]) else 0
}

# Year-round mean = cold + 0.5 * season-interaction (mathematically equal to
# averaging cold and warm season betas).
year_mean <- function(main, season_int) beta_of(main) + 0.5 * beta_of(season_int)

# Forest baseline — built from the per-individual free intercepts (forest has
# no population-level fixed effect in 03a; it is the free intercept instead).
forest_cold_mean <- mean(indiv_betas$Forest_free$post_mean)
forest_warm_mean <- mean(indiv_betas$Forest_free$post_mean +
                          indiv_betas$ForestS_free$post_mean)
muForest <- (forest_cold_mean + forest_warm_mean) / 2 + log(W_LARGE)

# Habitat dummies (cold + 0.5 * season interaction)
muSteppe <- year_mean("steppe",  "steppe_season")
muMeadow <- year_mean("meadow",  "meadow_season")
muAndean <- year_mean("andean",  "andean_season")
muSoil   <- beta_of("soil")                    # no season interaction in 03a

# Continuous covariates
muElev   <- year_mean("elev_scl",  "elev_season")
muElev2  <- beta_of("elev_scl2")
muBuild  <- year_mean("build_scl", "build_season")

cat("  Population-level year-round means:\n")
print(round(c(muForest = muForest, muSteppe = muSteppe, muMeadow = muMeadow,
              muAndean = muAndean, muSoil   = muSoil,
              muElev   = muElev,   muElev2  = muElev2,
              muBuild  = muBuild), 3))

# 3. Land-cover binary rasters ----
cat("\n=== Loading LCChalten_Pred ===\n")
lc <- rast(LC_TIF)
lc <- project(lc, TARGET_CRS_UTM, method = "near")

# Land-cover codes: 0 forest, 1 steppe, 2 meadow, 3 water, 4 snow,
# 5 highAndean, 6 bareSoil
binary_rasters <- list()
for (i in 0:6) binary_rasters[[i + 1]] <- as.numeric(lc == i)
names(binary_rasters) <- c("forest", "steppe", "meadow", "water",
                           "snow",   "highAndean", "bareSoil")

forest <- binary_rasters$forest
steppe <- binary_rasters$steppe
meadow <- binary_rasters$meadow
andean <- binary_rasters$highAndean
soil   <- binary_rasters$bareSoil

# 4. Scaled elevation + dist-to-build rasters ----
cat("\n=== Loading static stack and scaling continuous covariates ===\n")
stk <- rast(STATIC_STACK)
cat("  layers:", paste(names(stk), collapse = ", "), "\n")

elev_mean  <- scale_params$elevation$mean
elev_sd    <- scale_params$elevation$sd
build_mean <- scale_params$build$mean
build_sd   <- scale_params$build$sd

elevation_scl <- (stk[["elevation"]]     - elev_mean)  / elev_sd
build_scl     <- (stk[["dist_to_build"]] - build_mean) / build_sd

# Align all rasters to the land-cover grid
if (!compareGeom(elevation_scl, forest, stopOnError = FALSE)) {
  cat("  -> resampling continuous rasters onto the LC grid\n")
  elevation_scl <- resample(elevation_scl, forest, method = "bilinear")
  build_scl     <- resample(build_scl,     forest, method = "bilinear")
}

# 5. Group-level (year-round) suitability raster ----
cat("\n=== Computing log(RSF) ===\n")
log_rsf <- muForest * forest +
           muSteppe * steppe +
           muMeadow * meadow +
           muAndean * andean +
           muSoil   * soil   +
           muElev   * elevation_scl +
           muElev2  * elevation_scl^2 +
           muBuild  * build_scl

group_rsf <- exp(log_rsf)
names(group_rsf) <- "hRSF_suitability"
writeRaster(group_rsf, CONT_TIF, overwrite = TRUE,
            gdal = c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=IF_SAFER"))
cat(sprintf("  -> %s\n", CONT_TIF))

# 6. Quintile reclassification (very low ... very high) ----
cat("\n=== Reclassifying into 5 suitability classes ===\n")
q <- global(group_rsf, fun = function(v)
              quantile(v, probs = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                       na.rm = TRUE))
q <- as.numeric(q[1, ])
reclass_mat <- matrix(c(
  q[1], q[2], 1,
  q[2], q[3], 2,
  q[3], q[4], 3,
  q[4], q[5], 4,
  q[5], q[6], 5
), ncol = 3, byrow = TRUE)
group_rsf_cat <- classify(group_rsf, reclass_mat, include.lowest = TRUE)
levels(group_rsf_cat) <- data.frame(
  ID    = 1:5,
  class = c("Very low", "Low", "Moderate", "High", "Very high"))

# 7. Reproject raster + load overlay layers ----
# The raster is shown at full extent as the basemap; the protected-area
# polygon, Argentina-Chile boundary, rivers/lakes and Route 23 are drawn as
# line/polygon overlays on top (no crop or mask of the raster itself).
cat("\n=== Reprojecting raster + loading overlay layers ===\n")

# Reproject the 5-class raster to WGS84 for plotting in lon/lat
rsf_dd <- project(group_rsf_cat, "EPSG:4326")
levels(rsf_dd) <- data.frame(
  ID    = 1:5,
  label = c("Very low", "Low", "Moderate", "High", "Very high"))

writeRaster(rsf_dd, CLASS_TIF, overwrite = TRUE,
            gdal = c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=IF_SAFER"))
cat(sprintf("  -> %s\n", CLASS_TIF))

# Helper — read & reproject any shapefile to WGS84, with a NULL fallback
read_wgs84 <- function(path, label) {
  if (!file.exists(path)) {
    warning(sprintf("Shapefile missing (%s): %s — will be skipped",
                    label, path))
    return(NULL)
  }
  tryCatch({
    sf_obj <- st_read(path, quiet = TRUE)
    if (is.na(st_crs(sf_obj))) st_crs(sf_obj) <- TARGET_CRS_UTM
    st_make_valid(st_transform(sf_obj, 4326))
  }, error = function(e) {
    warning(sprintf("Failed to read %s: %s", path, conditionMessage(e)))
    NULL
  })
}

pa_merged    <- read_wgs84(PA_MERGED_SHP,    "Merged protected areas")
arg_ch_sf    <- read_wgs84(ARG_CH_SHP,       "Argentina-Chile boundary")
rivlakes_sf  <- read_wgs84(RIVERS_LAKES_SHP, "Rivers and lakes")
road_sf      <- read_wgs84(ROAD_SHP,         "Ruta 23")

# Plot extent = raster extent shrunk inward by a small fraction so the PA
# polygon overlay always sits on top of raster pixels (no white sliver
# between the dashed PA boundary and the raster edge).
deg_per_km_lat <- 1 / 111
deg_per_km_lon <- 1 / (111 * cos(49 * pi / 180))

ZOOM_PAD_FRAC <- 0.04          # shrink plot extent by 4 % on each side
Y_SOUTH_CAP   <- -49.45        # hard floor on the southern latitude shown

x_full <- c(xmin(rsf_dd), xmax(rsf_dd))
y_full <- c(ymin(rsf_dd), ymax(rsf_dd))
x_lim  <- c(x_full[1] + ZOOM_PAD_FRAC * diff(x_full),
            x_full[2] - ZOOM_PAD_FRAC * diff(x_full))
y_lim  <- c(y_full[1] + ZOOM_PAD_FRAC * diff(y_full),
            y_full[2] - ZOOM_PAD_FRAC * diff(y_full))
if (y_lim[1] < Y_SOUTH_CAP) {
  cat(sprintf("  Clamping southern y_lim from %.3f to Y_SOUTH_CAP = %.3f\n",
              y_lim[1], Y_SOUTH_CAP))
  y_lim[1] <- Y_SOUTH_CAP
}
cat(sprintf("  Plot extent: lon [%.3f, %.3f], lat [%.3f, %.3f]\n",
            x_lim[1], x_lim[2], y_lim[1], y_lim[2]))

# 8. JPG export ----
cat("\n=== Plotting suitability map ===\n")

# Crop the raster JUST FOR DISPLAY purposes (avoids drawing pixels outside
# the desired extent, which would otherwise force the panel to expand).
# This is a non-destructive crop on a copy — the saved CLASS_TIF above is
# still the full extent.
display_ext <- ext(x_lim[1], x_lim[2], y_lim[1], y_lim[2])
rsf_display <- tryCatch(crop(rsf_dd, display_ext),
                        error = function(e) rsf_dd)

# Aesthetic constants for the overlays (used in both the plot calls and the
# combined legend below)
RIVER_COL  <- "steelblue1"
ROAD_COL   <- "red"
PA_COL     <- "black"
PA_LTY     <- 2          # dashed
ARG_CH_COL <- "orange"

# Default terra::plot colour palette for the 5-class factor raster.
# Pulled directly so the legend swatches match the raster fills exactly.
SUIT_COLS  <- tryCatch(terra::map.pal("viridis", 5),
                       error = function(e) hcl.colors(5, palette = "viridis"))
SUIT_LABS  <- c("Very low", "Low", "Moderate", "High", "Very high")

# Wider JPG + generous right margin to host the legend OUTSIDE the plot.
# At 2400 x 1800 + mar[4] = 12 lines, the right margin holds the full
# combined legend (longest label "Rivers and lakes") without clipping.
jpeg(OUT_JPG, width = 2400, height = 1800, res = 300)
op <- par(mar = c(4, 4, 3, 12))        # generous right margin for legend

# `legend = FALSE` suppresses terra's automatic categorical legend on the
# right — the combined legend at the end of this block replaces it.
plot(rsf_display, main = "Huemul Habitat Suitability",
     axes = TRUE, xlim = x_lim, ylim = y_lim,
     col = SUIT_COLS, legend = FALSE)
grid(col = "grey60", lty = 3)

# North arrow — upper-left of panel
x_n <- x_lim[1] + 0.05 * diff(x_lim)
y_n <- y_lim[1] + 0.85 * diff(y_lim)
arrows(x_n, y_n, x_n, y_n + 0.06 * diff(y_lim),
       length = 0.10, col = "black", lwd = 2)
text(x_n + 0.08 * diff(x_lim),
     y_n + 0.09 * diff(y_lim),
     "N", cex = 1.6, font = 2)

# --- Overlays (drawn ON TOP of the raster, in this order) -------------------
# PolyLine shapefiles (Route 23 + Arg-Ch boundary) use `col = ...` to colour
# the line itself; `border = ...` only applies to polygons.
#
# 1. Rivers and lakes (polygon) — light blue filled
if (!is.null(rivlakes_sf))
  plot(st_geometry(rivlakes_sf), col = RIVER_COL,
       border = RIVER_COL, lwd = 0.5, add = TRUE)
# 2. Route 23 (PolyLine) — red, lwd 2
if (!is.null(road_sf))
  plot(st_geometry(road_sf), col = ROAD_COL,
       lwd = 2, add = TRUE)
# 3. Argentina-Chile boundary (PolyLine) — orange, lwd 3
if (!is.null(arg_ch_sf))
  plot(st_geometry(arg_ch_sf), col = ARG_CH_COL,
       lwd = 3, add = TRUE)
# 4. Merged protected areas (Polygon) — black dashed outline, no fill
if (!is.null(pa_merged))
  plot(st_geometry(pa_merged), col = NA, border = PA_COL,
       lwd = 1.5, lty = PA_LTY, add = TRUE)

# Scale bar — drawn LAST so it stays on top of every overlay
sb_dx <- 10 * deg_per_km_lon            # 10 km wide
sb_x0 <- x_lim[2] - 0.05 * diff(x_lim) - sb_dx
sb_y0 <- y_lim[1] + 0.04 * diff(y_lim)
sb_h  <- 0.01 * diff(y_lim)
rect(sb_x0,           sb_y0, sb_x0 + sb_dx/2, sb_y0 + sb_h,
     col = "black", border = "black")
rect(sb_x0 + sb_dx/2, sb_y0, sb_x0 + sb_dx,   sb_y0 + sb_h,
     col = "white", border = "black")
text(sb_x0,             sb_y0 - 0.012 * diff(y_lim), "0",  cex = 0.8)
text(sb_x0 + sb_dx/2,   sb_y0 - 0.012 * diff(y_lim), "5",  cex = 0.8)
text(sb_x0 + sb_dx,     sb_y0 - 0.012 * diff(y_lim), "10", cex = 0.8)
text(sb_x0 + sb_dx/2,   sb_y0 - 0.028 * diff(y_lim), "km", cex = 0.8)

# --- Single combined legend ------------------------------------------------
# Layout: raster classes first (filled square icons, pch = 22 with pt.bg =
# fill colour), then the overlay layers (line icons with lty/lwd):
#   * raster entries: pch = 22, pt.bg = colour, lty = NA   -> filled square
#   * overlay entries: pch = NA, lty = drawn,   pt.bg = NA -> line only
combined_labels <- c("Suitability",
                     SUIT_LABS,
                     "",                                # blank separator row
                     "Watercourses \nand lakes",
                     "Route 23",
                     " Argentinian \nprotected areas",
                     "Arg-Ch boundary")
combined_pch    <- c(NA, rep(22, length(SUIT_LABS)), NA, NA, NA, NA, NA)
combined_ptbg   <- c(NA, SUIT_COLS,                  NA, NA, NA, NA, NA)
combined_col    <- c(NA, rep("black", length(SUIT_LABS)),
                     NA, RIVER_COL, ROAD_COL, PA_COL, ARG_CH_COL)
combined_lty    <- c(NA, rep(NA, length(SUIT_LABS)), NA, 1, 1, PA_LTY, 1)
combined_lwd    <- c(NA, rep(NA, length(SUIT_LABS)), NA, 4, 2, 1.5, 3)
combined_font   <- c(2, rep(1, length(SUIT_LABS) + 5))   # bold header

# Legend anchored just outside the plot's right edge, extending into the
# generous right margin set by `par(mar = ...)` above. xpd = NA allows
# drawing outside the plot region; text.width pre-allocates space so R
# doesn't clip the longest label.
leg_anchor_x <- x_lim[2] + 0.02 * diff(x_lim)     # just past right plot edge
leg_anchor_y <- y_lim[2]                          # at top of plot
LEG_CEX <- 0.9
leg_text_width <- max(strwidth(combined_labels, cex = LEG_CEX))

legend(x         = leg_anchor_x,
       y         = leg_anchor_y,
       xjust     = 0,                              # left-anchor: extends right
       yjust     = 1,                              # top-anchor: extends down
       legend    = combined_labels,
       pch       = combined_pch,
       pt.bg     = combined_ptbg,
       pt.cex    = 1.4,
       col       = combined_col,
       lty       = combined_lty,
       lwd       = combined_lwd,
       text.font = combined_font,
       bty       = "n",                            # no box at all
       cex       = LEG_CEX,
       text.width = leg_text_width,                # full label width budget
       y.intersp = 1.25,
       seg.len   = 2.3,
       xpd       = NA)                             # allow drawing in margin

par(op)
dev.off()
cat(sprintf("  -> %s\n", OUT_JPG))

cat("\n=== Done ===\n")
