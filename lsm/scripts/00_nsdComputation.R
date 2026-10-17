# =============================================================================
# Movement behavior and habitat selection to guide conservation of an
# endangered Neotropical deer
#
# Alarcón, P.A.E., Crego, R.D., Galetto, E., Aprile, G.,
# Morales, J.M., Pauli, J.N. & Donadio, E.
#
# -----------------------------------------------------------------------------
# Script:  00_nsdComputation.R
# Purpose: Compute per-individual Net Squared Displacement (NSD) time series
#          from the cleaned GPS data, assemble them in the format expected
#          by the latent-state model (one record per day per animal, chosen
#          at random when multiple fixes exist), and derive the N-to-S
#          ordering lookups used by downstream LSM scripts and figures.
#          Also produces two supplementary inspection figures:
#            lsm/figures/FigS1_NSD_timeseries.png
#            lsm/figures/FigS1_monitoring_periods.png
#          Output: ./collarData/data/gpsdata_for_lsm.rda
# =============================================================================

# Files: ----
# 1. Metadata -----
load(file = "./collarData/data/huemul_metadata.rda")

# Build a quick animal_id -> capture_date lookup from the metadata. We
# need this because `min(gpsdata$Date)` is NOT a reliable proxy for the
# actual capture date — some collars only begin transmitting weeks (or
# months) after deployment, which silently mis-classified the capture
# season as autumn for animals tagged in spring (e.g. Huemul-8/10/11/12
# captured in October 2021 but first GPS fix in April 2022).
# The metadata's capture_date column is in DD/MM/YYYY format.
capture_date_lookup <- setNames(
  as.Date(metadatos$capture_date, format = "%d/%m/%Y"),
  metadatos$animal_id)

# 2. Collar data ----
load("./collarData/data/clean_gpsdata.rda")

# 3. Net Squared Displacement (NSD) computation ----
## 3.1a. Subsetting ----
 gpsdata$Date     <- gpsdata$date
 gpsdata$DateTime <- gpsdata$local_datetime
 # short_id is created in 0. cleaningData/0. scripts/1. cleaningData.R
 # and is already present in 1. clean_gpsdata.rda

 gpsdata <- subset(gpsdata, is.na(gpsdata$latitude) == "FALSE")
 
# retain ONE location per day per animal, chosen at random
library(dplyr)
set.seed(2026)

gpsdata <- gpsdata %>%
  group_by(animal_id, Date) %>%
  slice_sample(n = 1) %>%
  ungroup() %>%
  as.data.frame()

gpsdataL <- split(gpsdata, gpsdata$animal_id)
nid <- length(gpsdataL)

# retain animal ids with > 274 days of monitoring (> 3/4 annual cycle)
## count the days of monitoring for each animal id 
ndays_id <- numeric(nid)

for (g in 1:nid){
  ndays_id[g] <- length(unique(gpsdataL[[g]]$Date))
}

## identify animal ids with > 274 days of monitoring and subset
useful_gth <- which(ndays_id > 274)
useful_ids <- unique(gpsdata$animal_id)[useful_gth] 

gpsdata <- gpsdata[gpsdata$animal_id %in% useful_ids, ]

gpsdata <- gpsdata[,!names(gpsdata) == "dup"]

## 3.1b. Trim spring-captured individuals to following austral autumn ----
# Animals captured in spring/summer (Sep–Mar) start their GPS series in the
# summering range. Data are trimmed to April 1 of the following year so that
# all series begin during austral autumn, making individual NSD trajectories
# directly comparable without modifying the model structure.

gpsdataL <- split(gpsdata, gpsdata$animal_id)
nid      <- length(gpsdataL)

# Record raw monitoring periods before trimming (used in section 3.5b for the
# monitoring-period plot showing excluded vs retained segments per individual).
# short_id and lat_centroid are captured here so dropped animals (those that
# fall below the 274-day threshold after trimming) can still be plotted.
orig_periodL <- lapply(gpsdataL, function(df) {
  # Classify capture season from the metadata's REAL capture_date, not
  # min(df$Date) — collars that started transmitting months after deployment
  # were silently mis-classified as autumn before.
  animal_id_g <- unique(df$animal_id)[1]
  cap_date    <- capture_date_lookup[[animal_id_g]]
  if (is.null(cap_date) || is.na(cap_date))
    cap_date  <- as.Date(min(df$Date, na.rm = TRUE))
  cap_month <- as.integer(format(cap_date, "%m"))

  list(
    start          = min(df$Date, na.rm = TRUE),
    end            = max(df$Date, na.rm = TRUE),
    capture_season = if (cap_month >= 9 | cap_month <= 3) "spring" else "autumn",
    short_id       = unique(na.omit(df$short_id))[1],
    lat_centroid   = mean(df$latitude, na.rm = TRUE)
  )
})

for (g in 1:nid) {
  animal_id_g <- unique(gpsdataL[[g]]$animal_id)[1]
  cap_date    <- capture_date_lookup[[animal_id_g]]

  # If the animal is missing from the metadata fall back to min(GPS), but
  # warn so the user knows. (Should never trigger in practice.)
  if (is.null(cap_date) || is.na(cap_date)) {
    warning("No capture_date in metadatos for ", animal_id_g,
            " — falling back to min(GPS).")
    cap_date <- as.Date(min(gpsdataL[[g]]$Date, na.rm = TRUE))
  }
  cap_month <- as.integer(format(cap_date, "%m"))
  cap_year  <- as.integer(format(cap_date, "%Y"))

  # Spring/summer capture: September (9) through March (3) in the
  # Southern Hemisphere.
  is_spring <- cap_month >= 9 | cap_month <= 3

  if (is_spring) {
    # April 1 of the following austral autumn
    next_year   <- ifelse(cap_month >= 9, cap_year + 1L, cap_year)
    next_april1 <- as.POSIXct(paste0(next_year, "-04-01"),
                               format = "%Y-%m-%d", tz = "GMT")
    gpsdataL[[g]] <- gpsdataL[[g]][gpsdataL[[g]]$Date >= next_april1, ]
  }
}

gpsdata <- do.call(rbind.data.frame, gpsdataL)

# Re-apply >274-day filter in case trimming shortened any series
gpsdataL   <- split(gpsdata, gpsdata$animal_id)
ndays_trim <- sapply(gpsdataL, function(x) length(unique(x$Date)))
useful_ids <- names(ndays_trim)[ndays_trim > 274]
gpsdata    <- gpsdata[gpsdata$animal_id %in% useful_ids, ]
message("Animals retained after spring trim + >274-day filter: ",
        paste(sort(useful_ids), collapse = ", "))

## graphic summary of the data used for model estimation
nid <- length(unique(gpsdata$animal_id))

frec_id <- as.data.frame(table(gpsdata$animal_id))$Freq
tmp <- which(frec_id == 0)

npts_id <- frec_id[-tmp]
col_id<- rainbow(nid)
  
barplot(table(gpsdata$animal_id), 
        col  = col_id,
        names.arg = NULL, 
        axes = T, 
        las =2)

## 3.2. Time series of NSD (reference point: first retained GPS location) ----
# NSD is computed relative to the first location after the spring trim, so the
# reference point is the start of the analysed period for all individuals.

## create an empty column to save nsd values in kilometres
gpsdata$nsd_km <- NA

## create list of individual datasets
gpsdataL <- split(gpsdata, gpsdata$animal_id)
nid <- length(gpsdataL)
nxid <- as.data.frame(table(gpsdata$animal_id))$Freq

# computation of nsd time series
for (g in 1:nid) {
  gpsdataL[[g]]$nsd_km[1] <- 0
  
  for (i in 2:nxid[g]){
    
    # raw nsd
    
    gpsdataL[[g]]$nsd_km[i] <- (
      sqrt(
        (gpsdataL[[g]]$x_utm[i] - gpsdataL[[g]]$x_utm[1]) ^ 2 +
          (gpsdataL[[g]]$y_utm[i] - gpsdataL[[g]]$y_utm[1]) ^ 2
      ) 
      /1000)
  }
}  

# scaled nsd

## function to scale
range01 <- function(x){ 
  (x - min(x))/(max(x) - min(x)) 
}

for (g in 1:nid) {
  gpsdataL[[g]]$nsd_scl[2:nxid[g]] <- range01(gpsdataL[[g]]$nsd_km[2:nxid[g]])
 rownames(gpsdataL[[g]]) <- NULL
  }

gpsdata <- do.call(rbind.data.frame, gpsdataL)

## 3.3. Check out that time series are regular ----
# Insert missing dates within each individual's monitoring period.
library(dplyr)

## split into individual datasets
gpsdataL <- split(gpsdata, gpsdata$animal_id)
nid      <- length(gpsdataL)

## Insert missed dates
nobs  <- numeric(nid)
nmiss <- numeric(nid)

for (g in 1:nid){
  # Regular daily series between first and last observation
  dates_gth <- seq.POSIXt(min(gpsdataL[[g]]$Date, na.rm = TRUE),
                           max(gpsdataL[[g]]$Date, na.rm = TRUE),
                           by = "day")
  dates_gth <- as.POSIXct(dates_gth,
                           format = "%Y-%m-%d",
                           tz = "GMT", usetz = TRUE)

  dates_obs_gth  <- which(dates_gth %in% gpsdataL[[g]]$Date)
  nobs[g]        <- length(dates_obs_gth)

  dates_miss_gth <- which(!dates_gth %in% gpsdataL[[g]]$Date)
  nmiss[g]       <- length(dates_miss_gth)

  ncol_gth  <- ncol(gpsdataL[[g]])
  emptyness <- as.data.frame(matrix(NA, nmiss[g], ncol_gth))
  colnames(emptyness) <- names(gpsdataL[[g]])
  emptyness$Date      <- dates_gth[dates_miss_gth]
  emptyness$site      <- rep(unique(gpsdataL[[g]]$site),                     nmiss[g])
  emptyness$animal_id <- rep(unique(gpsdataL[[g]]$animal_id),                nmiss[g])
  emptyness$sex       <- rep(unique(gpsdataL[[g]]$sex),                      nmiss[g])
  emptyness$short_id  <- rep(unique(na.omit(gpsdataL[[g]]$short_id)),        nmiss[g])

  gpsdataL[[g]] <- rbind(gpsdataL[[g]], emptyness)
  gpsdataL[[g]] <- gpsdataL[[g]][order(gpsdataL[[g]]$Date), ]
}

gpsdata <- do.call(rbind.data.frame, gpsdataL)

## 3.4. Create numerical (integer) day ----
# Assign a sequential integer day index starting at 1 for each individual's
# first retained observation. After the spring trim (section 3.1b) all series
# begin in austral autumn; the index is used as the time axis in the HMM.
nid <- length(gpsdataL)

for (g in 1:nid){
  gpsdataL[[g]] <- subset(gpsdataL[[g]],
                           !is.na(gpsdataL[[g]]$Date))
  gpsdataL[[g]]$day <- seq_len(nrow(gpsdataL[[g]]))
}
 
 gpsdata <- do.call(rbind.data.frame, gpsdataL)
 
## 3.5. Graphics ----

# Latitudinal order: individuals sorted north → south by mean latitude.
# Two versions are maintained:
#   - all_lat_order / all_short_lat_order : ALL animals that passed the first
#     >274-day filter (pre-trim), used in the monitoring-periods plot so that
#     dropped animals still appear.
#   - id_lat_order / short_lat_order      : RETAINED animals only (post-trim),
#     used in NSD panels and all model-output figures.
# Both are saved in the .rda for downstream scripts.

all_lat_order       <- names(sort(sapply(orig_periodL, `[[`, "lat_centroid"),
                                  decreasing = TRUE))
all_short_lat_order <- sapply(all_lat_order,
                               function(id) orig_periodL[[id]]$short_id)

lat_centroid    <- sapply(gpsdataL, function(df) mean(df$latitude, na.rm = TRUE))
id_lat_order    <- names(sort(lat_centroid, decreasing = TRUE))
short_lat_order <- sapply(id_lat_order, function(id)
  unique(na.omit(gpsdata$short_id[gpsdata$animal_id == id]))[1])

### 3.5a. NSD time series (all individuals) ----
# Panels ordered north → south (top-left = northernmost individual).
# Rebuilt with ggplot + facet_wrap so the elevation gradient is accompanied
# by a shared colour-bar legend.
library(ggplot2)

gpsdataL_lat <- gpsdataL[id_lat_order]   # reorder list by latitude
nid          <- length(gpsdataL_lat)

# Flat data frame (all individuals), with short_id as a factor whose levels
# follow the N→S order so facet_wrap draws them in that order.
ns_df <- do.call(rbind, lapply(gpsdataL_lat, function(df) {
  data.frame(short_id = unique(na.omit(df$short_id))[1],
             Date     = df$Date,
             nsd_scl  = df$nsd_scl,
             altitude = df$altitude,
             stringsAsFactors = FALSE)
}))
ns_df$short_id <- factor(ns_df$short_id, levels = short_lat_order)

p_nsd_ts <- ggplot(ns_df, aes(x = Date, y = nsd_scl)) +
  geom_line(colour = "gray70", linewidth = 0.4) +
  geom_point(aes(colour = altitude), size = 0.9) +
  scale_colour_gradientn(
    colours = rev(heat.colors(100)),       # yellow (low) → red (high)
    name    = "Elevation (m a.s.l.)",
    guide   = guide_colourbar(barwidth = 15, barheight = 0.8,
                              title.position = "top")) +
  facet_wrap(~ short_id, ncol = 2, scales = "free_x") +
  labs(x = "Date", y = "NSD (scaled)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        strip.text      = element_text(face = "bold"),
        panel.grid.minor = element_blank())

ggsave("./lsm/figures/FigS1_NSD_timeseries.png",
       p_nsd_ts, width = 10, height = 14, dpi = 300, bg = "white")

### 3.5b. Monitoring periods ----
# Horizontal bar chart: full GPS monitoring period (grey) overlaid with the
# period retained in the analysis (coloured by capture season).
# Autumn captures (Apr–Aug, amber): no excluded pre-period, bar is entirely coloured.
# Spring captures (Sep–Mar, teal): excluded spring/summer segment shown in grey;
# retained segment (from April 1 of the following year onward) shown in teal.
# Open circle marks the capture date (start of the raw GPS record).
# Individuals ordered north → south (top = northernmost), matching 3.5a.

library(ggplot2)

# Build the period data frame from the raw periods (recorded before trimming)
# and the trimmed start/end dates (from the final gpsdata).
# ALL animals that passed the first >274-day filter are included; those that
# were subsequently dropped by the re-filter (spring trim reduced their series
# below 274 days) receive incl_start/incl_end = NA and retained = FALSE.
period_df <- do.call(rbind, lapply(names(orig_periodL), function(id) {
  op  <- orig_periodL[[id]]
  sub <- gpsdata[gpsdata$animal_id == id, ]
  rs  <- as.Date(op$start)
  re  <- as.Date(op$end)
  data.frame(
    animal_id      = id,
    short_id       = op$short_id,
    raw_start      = rs,
    raw_end        = re,
    incl_start     = if (nrow(sub) > 0) as.Date(min(sub$Date, na.rm = TRUE)) else as.Date(NA),
    incl_end       = if (nrow(sub) > 0) as.Date(max(sub$Date, na.rm = TRUE)) else as.Date(NA),
    capture_season = op$capture_season,
    n_days         = as.integer(re - rs),   # total days: latest GPS − capture date
    retained       = nrow(sub) > 0,
    stringsAsFactors = FALSE
  )
}))

# Apply all-animal latitudinal order (north at top).
# ggplot puts the last factor level at the top of a discrete y-axis,
# so levels run south → north (= reverse of north → south order).
period_df$short_id <- factor(period_df$short_id,
                              levels = rev(all_short_lat_order))

season_cols <- c(autumn = "deepskyblue4", spring = "green4")

gg_periods <- ggplot(period_df) +
  # Full raw monitoring bar (light grey background stripe, all individuals)
  geom_segment(aes(x = raw_start, xend = raw_end,
                   y = short_id,  yend = short_id),
               colour = "grey78", linewidth = 4) +
  # Retained period (colour-coded by capture season, retained animals only)
  geom_segment(data = period_df[period_df$retained, ],
               aes(x = incl_start, xend = incl_end,
                   y = short_id,   yend = short_id,
                   colour = capture_season),
               linewidth = 4) +
  # Capture date marker (open circle = start of raw monitoring)
  geom_point(aes(x = raw_start, y = short_id),
             shape = 21, fill = "white", colour = "grey40", size = 3) +
  # Total monitoring days label at the right end of each bar
  geom_text(aes(x = raw_end, y = short_id,
                label = paste0(n_days, "")),
            hjust = -0.2, size = 3.2, colour = "grey30") +
  scale_colour_manual(
    values = season_cols,
    labels = c(autumn = "Autumn capture (Apr\u2013Aug)",
               spring = "Spring capture (Sep\u2013Mar), trimmed to following April"),
    name   = NULL
  ) +
  scale_x_date(date_breaks = "6 months", date_labels = "%b %Y",
               expand = expansion(mult = c(0.02, 0.12))) +
  labs(
    x        = NULL,
    y        = NULL,
    title    = "Monitoring periods",
    subtitle = paste0("Grey = full GPS record\u2003",
                      "Coloured = period retained in analysis\u2003",
                      "Circle = capture date")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    axis.text.x        = element_text(angle = 35, hjust = 1),
    legend.position    = "bottom",
    plot.title         = element_text(face = "bold"),
    plot.subtitle      = element_text(size = 9, colour = "grey40")
  )

ggsave("./lsm/figures/FigS1_monitoring_periods.png",
       plot = gg_periods,
       width = 10, height = 6, units = "in", dpi = 300)

# save a rda file
# id_lat_order and short_lat_order are saved alongside gpsdata so all
# downstream scripts can reuse the same north-to-south individual ordering.
save(gpsdata, id_lat_order, short_lat_order, all_lat_order, all_short_lat_order,
     file = "./collarData/data/gpsdata_for_lsm.rda")

