# =============================================================================
# Movement behavior and habitat selection to guide conservation of an
# endangered Neotropical deer
#
# Alarcón, P.A.E., Crego, R.D., Galetto, E., Aprile, G.,
# Morales, J.M., Pauli, J.N. & Donadio, E.
#
# -----------------------------------------------------------------------------
# Script:  00_cleanRawData.R
# Purpose: Convert the raw GPS collar export + metadata CSV into the single
#          cleaned GPS dataset consumed by every downstream module. Steps:
#          merge collar fixes with individual metadata, parse date/time
#          fields, trim fixes to the strict monitoring period (capture_date
#          < fix_date < death_date), project latlong to UTM 19S, drop
#          unused fields, standardise sex to F/M, build the short animal
#          ID (sex initial + numeric tag), filter by DOP ≤ 5, add
#          day/night (nautical twilight, suncalc) and austral season
#          (cold = 21 Mar – 20 Sep; warm = 21 Sep – 20 Mar) classifiers.
#          Output: ./collarData/data/clean_gpsdata.{csv,rda}
# =============================================================================

# 0. Libraries ----
library(lubridate)
library(sp)
library(suncalc)

# 1. Load raw collar data + metadata ----
gpsdata   <- read.csv("./collarData/data/raw_GPS_042226_124836.CSV", sep = ",")
metadatos <- read.csv("./collarData/data/huemul_metadata.csv",       sep = ",")
save(metadatos, file = "./collarData/data/huemul_metadata.rda")

# 2. Merge collar data with metadata ----
gpsdata <- merge(gpsdata, metadatos,
                 by.x = "Device.ID", by.y = "collar_id")

# 3. Parse date / time fields ----
gpsdata$DateTime <- as.POSIXct(strptime(gpsdata$Date...Time..Local.,
                                        format = "%Y-%m-%d %H:%M:%S",
                                        tz = "America/Argentina/Buenos_Aires"))
gpsdata$Date <- as.POSIXct(strptime(gpsdata$DateTime,
                                    format = "%Y-%m-%d",
                                    tz = "America/Argentina/Buenos_Aires"))
gpsdata$CaptureDate <- as.POSIXct(strptime(gpsdata$capture_date,
                                           format = "%d/%m/%Y",
                                           tz = "America/Argentina/Buenos_Aires"))
gpsdata$StopDate <- as.POSIXct(strptime(gpsdata$death_date,
                                        format = "%d/%m/%Y",
                                        tz = "America/Argentina/Buenos_Aires"))

# 4. Keep fixes strictly between capture date and end of monitoring. The
# comparison is date-level (as.Date) so every fix recorded on the capture
# day and on the last monitoring day is excluded, regardless of time.
gpsdata$fix_date <- as.Date(gpsdata$DateTime)
gpsdata <- gpsdata[gpsdata$fix_date > as.Date(gpsdata$CaptureDate) &
                   gpsdata$fix_date < as.Date(gpsdata$StopDate), ]

# 5. Project latlong to UTM 19S ----
pos  <- which(is.na(gpsdata$Latitude) == "FALSE")
data <- gpsdata[pos, ]

pts <- SpatialPoints(cbind(data$Longitude, data$Latitude),
                     CRS("+proj=longlat +datum=WGS84"))
pts <- spTransform(pts,
                   CRS("+proj=utm +zone=19 +south +datum=WGS84"))
pts <- SpatialPointsDataFrame(pts, data = data)

gpsdata$x_utm <- NA
gpsdata$y_utm <- NA
gpsdata$x_utm[pos] <- pts@coords[, 1]
gpsdata$y_utm[pos] <- pts@coords[, 2]

# 6. Drop unused fields ----
drop <- c("Device.Name", "Date...Time..GMT.", "Date...Time..Local.",
          "Fix.Status", "Temp..C.", "Main..V.", "Back..V.", "capture_date",
          "death_date")
gpsdata <- gpsdata[, !names(gpsdata) %in% drop]

# 7. Reorder and rename fields ----
col_order <- c("animal_id", "Device.ID", "CaptureDate", "StopDate",
               "sex", "weight_kg", "DateTime", "Date",
               "Latitude", "Longitude", "x_utm", "y_utm",
               "Altitude", "DOP")
headings <- c("animal_id", "collar_id", "capture_date", "death_date",
              "sex", "body_weight_kg", "local_datetime", "date",
              "latitude", "longitude", "x_utm", "y_utm", "altitude",
              "dop")
gpsdata <- gpsdata[, col_order]
colnames(gpsdata) <- headings

# Standardise sex to F/M (metadata stores Spanish values).
gpsdata$sex <- ifelse(gpsdata$sex == "hembra", "F", "M")

# Short individual ID: sex initial (F/M) + numeric tag without leading zeros.
# E.g. F + "Huemul-07" -> "F7"; M + "Huemul-5" -> "M5". Used as display label.
gpsdata$short_id <- paste0(
  gpsdata$sex,
  gsub("Huemul-0*([0-9]+)", "\\1", gpsdata$animal_id)
)

# 8. Filter out fixes with DOP > 5 ----
gpsdata <- subset(gpsdata, gpsdata$dop <= 5)

# 9. Day / night classifier (nautical twilight) ----
dats <- as.data.frame(cbind(rep(NA, length(gpsdata$latitude)),
                            gpsdata$latitude, gpsdata$longitude))
dats$V1 <- ymd(gpsdata$date)
colnames(dats) <- c("date", "lat", "lon")

times <- getSunlightTimes(data = dats,
                          tz   = "America/Argentina/Buenos_Aires",
                          keep = c("nauticalDawn", "nauticalDusk"))

gpsdata$night_end   <- times$nauticalDawn
gpsdata$night_start <- times$nauticalDusk
gpsdata$day_night   <- ifelse(gpsdata$local_datetime > gpsdata$night_start |
                               gpsdata$local_datetime <= gpsdata$night_end,
                              "night", "day")

# 10. Austral season classifier ----
gpsdata$local_datetime <- ymd_hms(gpsdata$local_datetime)

get_season <- function(date) {
  month <- as.numeric(format(date, "%m"))
  day   <- as.numeric(format(date, "%d"))
  if ((month == 3 & day >= 21) | month %in% c(4,5,6,7,8) | (month == 9 & day < 21))
    return("cold") else return("warm")
}
gpsdata$season <- sapply(gpsdata$local_datetime, get_season)

# 11. Write cleaned dataset ----
write.csv(gpsdata, file = "./collarData/data/clean_gpsdata.csv")
save(gpsdata,     file = "./collarData/data/clean_gpsdata.rda")
message("Saved: ./collarData/data/clean_gpsdata.{csv,rda}")
