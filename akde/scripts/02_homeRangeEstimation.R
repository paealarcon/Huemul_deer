# =============================================================================
# Movement behavior and habitat selection to guide conservation of an
# endangered Neotropical deer
#
# Alarcón, P.A.E., Crego, R.D., Galetto, E., Aprile, G.,
# Morales, J.M., Pauli, J.N. & Donadio, E.
#
# -----------------------------------------------------------------------------
# Script:  02_homeRangeEstimation.R
# Purpose: Estimate individual home ranges with the continuous-time
#          Autocorrelated Kernel Density Estimator (AKDE), fitted with an
#          anisotropic OUF movement model (Fleming et al. 2015)
# =============================================================================

# 0. Libraries ----
library(ctmm)
library(ggplot2)

# 1. Data ----
load("./collarData/data/clean_gpsdata.rda")

gpsdata <- subset(gpsdata, !is.na(longitude) & !is.na(latitude))

gpsdataL <- split(gpsdata, gpsdata$animal_id)
nid      <- length(gpsdataL)

# 2. AKDE home range estimation ----
# For each individual: format the fixes as a ctmm telemetry object, fit an
# anisotropic OUF movement model, and estimate the 95% AKDE home range.
ctmmL <- list()
akdeL <- list()

home_range_area <- data.frame(matrix(NA, nid, 5))
colnames(home_range_area) <- c("animal_id", "low", "est", "high", "model")

for (h in seq_len(nid)) {

  ctmmL[[h]] <- data.frame(
    animal_id = gpsdataL[[h]]$short_id,
    timestamp = as.character(gpsdataL[[h]]$local_datetime),
    latitude  = gpsdataL[[h]]$latitude,
    longitude = gpsdataL[[h]]$longitude
  )

  huemul <- ctmmL[[h]]
  huemul$timestamp <- as.POSIXct(huemul$timestamp, format = "%Y-%m-%d %H:%M:%S")

  telemetry_data <- as.telemetry(huemul)

  # tau of length 2 + isotropic = FALSE forces an anisotropic OUF model
  guess <- ctmm.guess(telemetry_data,
                      CTMM = ctmm(tau = c(1, 1), isotropic = FALSE),
                      interactive = FALSE)
  fit <- ctmm.fit(telemetry_data, guess)

  akdeL[[h]] <- akde(telemetry_data, fit)

  animal_id <- unique(gpsdataL[[h]]$animal_id)
  home_range_area[h, 1]   <- animal_id
  home_range_area[h, 2:4] <- summary(akdeL[[h]])$CI
  home_range_area[h, 5]   <- summary(guess)$name
}

save(akdeL,           file = "./akde/outputs/akdeL.rda")
save(home_range_area, file = "./akde/outputs/home_range_area.rda")

# 3. Ordering & plotting utilization distributions ----
load("./akde/outputs/akdeL.rda")

short_id_vec <- sapply(gpsdataL, function(d) unique(d$short_id)[1])
names(short_id_vec) <- names(gpsdataL)

centroids_lat  <- sapply(gpsdataL, function(d) mean(d$latitude, na.rm = TRUE))
north_to_south <- order(centroids_lat, decreasing = TRUE)

akdeL_ordered    <- akdeL[north_to_south]
short_id_ordered <- short_id_vec[north_to_south]

jpeg("./akde/figures/FigureS5a.jpg", width = 10, height = 10, units = "in", res = 300)
par(mfrow = c(4, 4))
for (h in seq_along(akdeL_ordered)) {
  plot(akdeL_ordered[[h]], main = short_id_ordered[h])
}
dev.off()

# 4. Home range size vs. latitude ----
centroids_df <- aggregate(latitude ~ animal_id, data = gpsdata,
                          FUN = mean, na.rm = TRUE)
colnames(centroids_df)[2] <- "mean_latitude"

akde_data      <- merge(home_range_area, centroids_df, by = "animal_id")
akde_data$est  <- as.numeric(akde_data$est)
akde_data$low  <- as.numeric(akde_data$low)
akde_data$high <- as.numeric(akde_data$high)

save(akde_data, file = "./akde/outputs/akde_data.rda")

plot_akde_lat <- ggplot(akde_data, aes(x = mean_latitude, y = est)) +
  geom_errorbar(aes(ymin = low, ymax = high), width = 0.05, color = "darkgray") +
  geom_point(color = "black", size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", formula = y ~ x, color = "black",
             linetype = "dashed", se = FALSE) +
  labs(
    title = "Huemul home range size (AKDE) vs. latitude",
    x = "Mean latitude (centroid)",
    y = expression("Home range area (km"^2*")")
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave("./akde/figures/FigureS5b.jpg", plot = plot_akde_lat,
       width = 7, height = 5, dpi = 300)
