# =============================================================================
# Movement behavior and habitat selection to guide conservation of an
# endangered Neotropical deer
#
# Alarcón, P.A.E., Crego, R.D., Galetto, E., Aprile, G.,
# Morales, J.M., Pauli, J.N. & Donadio, E.
#
# -----------------------------------------------------------------------------
# Script:  00_prepare_used_available.R
# Purpose: Build the used/available dataset for the hRSF.
#          Projects GPS fixes to UTM 19S, computes a two-step 95% MCP as the
#          availability domain, clips it to remove unavailable surfaces
#          (water, snow, elevations above ELEV_CEILING), splits each
#          individual into train/test by day, samples NAV available points
#          per used point, and merges everything into UsedAv_INLA.rda. Also
#          exports the clipped availability polygon as a shapefile and
#          computes the habitat composition inside it.
# =============================================================================

# 0. Libraries ---------------------------------------------------------------
suppressPackageStartupMessages({
  library(adehabitatHR); library(sp); library(sf); library(terra)
  library(lubridate); library(dplyr)
})

# 0.1 CONSTANTS -------------------------------------------------------------
PROJECT      <- "."
ELEV_CEILING <- 1650L    # observed maximum (1637 m) rounded up to 50 m
NAV          <- 10L      # available:used ratio
SEED_SPLIT   <- 123L
SEED_AVAIL   <- 124L

# 0.2 Paths -----------------------------------------------------------------
GPS_FILE          <- file.path(PROJECT, "collarData/data/clean_gpsdata.rda")
STACK_FILE        <- file.path(PROJECT, "data/layers/stack_static.tif")
LC_DIR            <- file.path(PROJECT, "data/layers/habitatTypes")
AVAIL_DIR         <- file.path(PROJECT, "data/layers/availableHabitat")

GPSDATA_OUT_DIR   <- file.path(PROJECT, "hrsf/data")
MODEL_OUT_DIR     <- file.path(PROJECT, "hrsf/outputs")

TRAIN_RDA         <- file.path(GPSDATA_OUT_DIR, "gpsdata_for_hrsf_train.rda")
TEST_RDA          <- file.path(GPSDATA_OUT_DIR, "gpsdata_for_hrsf_test.rda")
USEDAV_RDA        <- file.path(GPSDATA_OUT_DIR, "UsedAv_INLA.rda")
MCP_RDS           <- file.path(MODEL_OUT_DIR,   "mcp_polygons.rds")
COVERAGE_CSV      <- file.path(MODEL_OUT_DIR,   "mcp95_coverage_report.csv")
AVMCP_CLIPPED_DSN <- file.path(AVAIL_DIR,       "AvMCP95_clipped.shp")
HABITAT_PCT_CSV   <- file.path(MODEL_OUT_DIR,   "habitat_pct_within_AvMCP95_clipped.csv")

dir.create(GPSDATA_OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(MODEL_OUT_DIR,   showWarnings = FALSE, recursive = TRUE)
dir.create(AVAIL_DIR,       showWarnings = FALSE, recursive = TRUE)

# 1. Load GPS data -----------------------------------------------------------
cat("=== Loading GPS data ===\n")
load(GPS_FILE)   # provides `gpsdata`
cat("  Rows loaded:", nrow(gpsdata), "\n")

gpsdata$local_datetime <- ymd_hms(gpsdata$local_datetime)

# 1.1 Season (austral)
get_season <- function(date) {
  m <- as.numeric(format(date, "%m")); d <- as.numeric(format(date, "%d"))
  if ((m == 3 & d >= 21) | m %in% c(4,5,6,7,8) | (m == 9 & d < 21)) "cold"
  else                                                              "warm"
}
gpsdata$season <- sapply(gpsdata$local_date, get_season)

# 1.2 Drop rows with missing coordinates
gpsdata <- subset(gpsdata, !is.na(longitude) & !is.na(latitude))
cat("  Rows after dropping NA coords:", nrow(gpsdata), "\n")

# 2. Project to UTM 19S ------------------------------------------------------
cat("\n=== Projecting to UTM zone 19S ===\n")
pts <- SpatialPoints(cbind(gpsdata$longitude, gpsdata$latitude),
                     CRS("+proj=longlat +datum=WGS84"))
pts <- spTransform(pts, CRS("+proj=utm +zone=19 +south +datum=WGS84"))
pts_id <- SpatialPointsDataFrame(
  coords = coordinates(pts),
  data   = data.frame(id = gpsdata$animal_id),
  proj4string = CRS("+proj=utm +zone=19 +south +datum=WGS84"))

# 3. MCPs (100% + two-step 95%) ---------------------------------------------
cat("\n=== Computing MCPs ===\n")
mcp100 <- mcp(pts, percent = 100, unin = "m", unout = "km2")

cat("  Step A: per-individual 95% MCPs (n =",
    length(unique(gpsdata$animal_id)), " huemul)\n")
mcp95_per_id    <- mcp(pts_id, percent = 95, unin = "m", unout = "km2")
mcp95_per_id_sf <- st_as_sf(mcp95_per_id)

gpsdata_sf <- st_transform(
  st_as_sf(gpsdata, coords = c("longitude","latitude"), crs = 4326),
  crs = 32719)

inside_own <- logical(nrow(gpsdata_sf))
for (aid in unique(gpsdata_sf$animal_id)) {
  this_poly <- mcp95_per_id_sf[mcp95_per_id_sf$id == aid, ]
  if (nrow(this_poly) == 0) next
  rows <- gpsdata_sf$animal_id == aid
  inside_own[rows] <- lengths(
    st_intersects(gpsdata_sf$geometry[rows], this_poly$geometry)) > 0
}
cat("  Points kept after per-id 95% trim:",
    sum(inside_own), "/", length(inside_own),
    sprintf("(%.1f%%)\n", 100 * mean(inside_own)))

pts95 <- as(gpsdata_sf[inside_own, ], "Spatial")
mcp95_obj <- mcp(SpatialPoints(coordinates(pts95),
                                proj4string = CRS(proj4string(pts95))),
                  percent = 100, unin = "m", unout = "km2")
mcp95 <- st_as_sf(mcp95_obj); mcp95$id <- "all"

# 3.3 Per-individual coverage sanity check
coverage <- gpsdata_sf %>%
  mutate(inside = lengths(st_intersects(geometry, mcp95$geometry)) > 0) %>%
  st_drop_geometry() %>%
  group_by(animal_id) %>%
  summarise(n_total = n(), n_inside = sum(inside),
            pct_inside = round(100 * mean(inside), 2), .groups = "drop") %>%
  arrange(pct_inside)
print(coverage, n = Inf)
cat("\n  Mean coverage:", round(mean(coverage$pct_inside), 2), "%",
    "  Min:", round(min(coverage$pct_inside),  2), "%\n")

# 4. Available polygon — remove lakes/snow + cap at ELEV_CEILING ------------
cat(sprintf("\n=== Building habitat polygons + clipping (ELEV_CEILING = %d m) ===\n",
            ELEV_CEILING))

lc_list <- list(
  Forest   = st_read(file.path(LC_DIR, "forest.shp"),     quiet = TRUE),
  Steppe   = st_read(file.path(LC_DIR, "steppes.shp"),    quiet = TRUE),
  Meadow   = st_read(file.path(LC_DIR, "meadow.shp"),     quiet = TRUE),
  water    = st_read(file.path(LC_DIR, "water.shp"),      quiet = TRUE),
  snow     = st_read(file.path(LC_DIR, "snow.shp"),       quiet = TRUE),
  Andean   = st_read(file.path(LC_DIR, "highAndean.shp"), quiet = TRUE),
  bareSoil = st_read(file.path(LC_DIR, "bareSoil.shp"),   quiet = TRUE)
)
for (cn in names(lc_list)) lc_list[[cn]]$landCat <- cn
lc <- do.call(rbind, lc_list)
lc <- st_transform(st_transform(lc, crs = 4326), crs = 32719)

unavailable <- subset(lc, lc$landCat %in% c("water", "snow"))
mcp100 <- st_as_sf(mcp100)
st_erase <- function(x, y) st_difference(x, st_union(st_combine(y)))

mcpAv100 <- st_erase(mcp100, unavailable)
mcpAv95  <- st_erase(mcp95,  unavailable)

elev_rast    <- rast(STACK_FILE)$elevation
below_cap    <- ifel(elev_rast < ELEV_CEILING, 1, NA)
below_poly   <- as.polygons(below_cap) |> st_as_sf() |> st_transform(crs = 32719)

mcpAv100 <- st_collection_extract(
  st_intersection(mcpAv100, st_union(below_poly)), "POLYGON")
mcpAv95  <- st_collection_extract(
  st_intersection(mcpAv95,  st_union(below_poly)), "POLYGON")

# 4.1 Export the clipped availability polygon as a shapefile ----------------
suppressWarnings(st_write(mcpAv95, AVMCP_CLIPPED_DSN,
                           driver = "ESRI Shapefile",
                           delete_layer = TRUE, quiet = TRUE))
cat("  Saved clipped availability shapefile:", AVMCP_CLIPPED_DSN, "\n")

# 4.2 Habitat composition within the clipped availability polygon -----------
# Each land-cover polygon is intersected with mcpAv95 and the area summed.
# Percentages are relative to the total area of mcpAv95.
cat("\n=== Computing habitat composition within AvMCP95 (clipped) ===\n")
tot_area_m2 <- as.numeric(sum(st_area(mcpAv95)))

habitat_pct <- do.call(rbind, lapply(names(lc_list), function(nm) {
  lc_i   <- st_transform(lc_list[[nm]], st_crs(mcpAv95))
  inter  <- suppressWarnings(st_intersection(lc_i, mcpAv95))
  area_i <- if (nrow(inter) == 0) 0 else as.numeric(sum(st_area(inter)))
  data.frame(category       = nm,
             area_km2       = round(area_i / 1e6, 3),
             pct_of_AvMCP95 = round(100 * area_i / tot_area_m2, 2),
             stringsAsFactors = FALSE)
}))
habitat_pct <- habitat_pct[order(-habitat_pct$pct_of_AvMCP95), ]
print(habitat_pct, row.names = FALSE)
write.csv(habitat_pct, HABITAT_PCT_CSV, row.names = FALSE)
cat("  Saved habitat composition CSV:", HABITAT_PCT_CSV, "\n")

# 5. Train / test split (random by DAY, 4 fixes/day) ------------------------
cat("\n=== Train / test split ===\n")
set.seed(SEED_SPLIT)

split_and_sample <- function(data) {
  unique_days <- unique(data$date)
  n_days      <- length(unique_days)
  train_days  <- sample(unique_days, size = ceiling(0.7 * n_days))
  test_days   <- setdiff(unique_days, train_days)
  train_data  <- data[data$date %in% train_days, ]
  test_data   <- data[data$date %in% test_days,  ]
  sample_4    <- function(x) x[sample(nrow(x), min(4, nrow(x))), ]
  train_samp  <- do.call(rbind, by(train_data, train_data$date, sample_4))
  test_samp   <- do.call(rbind, by(test_data,  test_data$date,  sample_4))
  list(training = train_samp, test = test_samp)
}

split_list    <- lapply(split(gpsdata, gpsdata$animal_id), split_and_sample)
training_data <- do.call(rbind, lapply(split_list, `[[`, "training"))
test_data     <- do.call(rbind, lapply(split_list, `[[`, "test"))

save(training_data, file = TRAIN_RDA)
save(test_data,     file = TEST_RDA)
cat("  Train:", nrow(training_data), " |  Test:", nrow(test_data), "\n")

# 6. Annotate USED points with covariates ------------------------------------
cat("\n=== Annotating USED points with covariates ===\n")

gpsdataL <- split(training_data, training_data$animal_id)
nid      <- length(gpsdataL)
xyUsedL  <- vector("list", nid); n <- integer(nid)

for (h in seq_len(nid)) {
  xyUsedL[[h]] <- st_as_sf(gpsdataL[[h]],
                           coords = c("x_utm", "y_utm"), crs = 32719)
  n[h] <- nrow(gpsdataL[[h]])
}

stack <- rast(STACK_FILE)
UsedL <- vector("list", nid)

for (h in seq_len(nid)) {
  use         <- rep(1, n[h])
  elevation   <- terra::extract(stack$elevation,     xyUsedL[[h]])$elevation
  habitatCat  <- st_join(xyUsedL[[h]], lc)$landCat
  build       <- terra::extract(stack$dist_to_build, xyUsedL[[h]])$dist_to_build
  animal_id   <- gpsdataL[[h]]$animal_id
  sex         <- gpsdataL[[h]]$sex
  season_h    <- gpsdataL[[h]]$season

  UsedL[[h]] <- data.frame(animal_id, sex, season = season_h,
                            use, elevation, habitatCat, build)
}

for (h in seq_len(nid)) {
  UsedL[[h]] <- subset(UsedL[[h]],
                       !c(UsedL[[h]]$habitatCat == "snow" |
                          UsedL[[h]]$habitatCat == "water"))
  n[h] <- nrow(UsedL[[h]])
}

# 7. Sample AVAILABLE within mcpAv95 -----------------------------------------
cat(sprintf("\n=== Sampling AVAILABLE points within mcpAv95 (nav = %d) ===\n", NAV))
set.seed(SEED_AVAIL)
xyAvL <- vector("list", nid)

for (h in seq_len(nid)) {
  xyAvL[[h]] <- st_sample(mcpAv95, size = n[h] * NAV, type = "random", exact = TRUE)
  xyAvL[[h]] <- st_as_sf(xyAvL[[h]], coords = c("x", "y"), crs = 32719)
}

all_xyAv       <- do.call(rbind, xyAvL)
all_elevation  <- terra::extract(stack$elevation,     all_xyAv)$elevation
all_habitatCat <- st_join(all_xyAv, lc)$landCat
all_build      <- terra::extract(stack$dist_to_build, all_xyAv)$dist_to_build
all_animal_id  <- rep(unlist(lapply(UsedL, function(x) x$animal_id)), each = NAV)
all_sex        <- rep(unlist(lapply(UsedL, function(x) x$sex)),       each = NAV)
all_season     <- rep(unlist(lapply(UsedL, function(x) x$season)),    each = NAV)
all_use        <- rep(0, sum(n) * NAV)

all_av <- data.frame(animal_id  = all_animal_id, sex = all_sex, season = all_season,
                     use = all_use, elevation = round(all_elevation, 0),
                     habitatCat = all_habitatCat, build = all_build)

avL <- split(all_av, rep(seq_len(nid), times = n * NAV))

# 8. Merge used + available -> UsedAv ----------------------------------------
cat("\n=== Merging USED and AVAILABLE ===\n")
UsedAvL <- vector("list", nid)
for (h in seq_len(nid)) UsedAvL[[h]] <- rbind(avL[[h]], UsedL[[h]])
UsedAv <- do.call(rbind.data.frame, UsedAvL)
UsedAv <- UsedAv[!is.na(UsedAv$habitatCat) & !is.na(UsedAv$build), ]

save(UsedAv, file = USEDAV_RDA)
cat("  UsedAv rows:", nrow(UsedAv), " |  used =",
    sum(UsedAv$use == 1), " |  available =", sum(UsedAv$use == 0), "\n")

# 9. Persist MCP outputs + coverage ------------------------------------------
saveRDS(list(mcp100        = mcp100,
             mcp95_per_id  = mcp95_per_id_sf,
             mcp95         = mcp95,
             mcpAv100      = mcpAv100,       # mcp100 - water - snow - >ELEV_CEILING
             mcpAv95       = mcpAv95,        # mcp95  - water - snow - >ELEV_CEILING
             elev_ceiling  = ELEV_CEILING),
        file = MCP_RDS)
write.csv(coverage, COVERAGE_CSV, row.names = FALSE)

cat("\n=== DONE — outputs:\n")
cat("  ", TRAIN_RDA, "\n")
cat("  ", TEST_RDA, "\n")
cat("  ", USEDAV_RDA, "\n")
cat("  ", MCP_RDS, "\n")
cat("  ", COVERAGE_CSV, "\n")
cat("  ", AVMCP_CLIPPED_DSN, "\n", sep = "")
cat("  ", HABITAT_PCT_CSV, "\n", sep = "")
