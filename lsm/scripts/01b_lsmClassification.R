# =============================================================================
# Movement behavior and habitat selection to guide conservation of an
# endangered Neotropical deer
#
# Alarcón, P.A.E., Crego, R.D., Galetto, E., Aprile, G.,
# Morales, J.M., Pauli, J.N. & Donadio, E.
#
# -----------------------------------------------------------------------------
# Script:  01b_lsmClassification.R
# Purpose: Classify individuals as resident, nomadic, migratory, disperser
#          or unknown, following Bastille-Rousseau et al. (2016)
# =============================================================================

# 0. Libraries ----


# 1. Data ----
load("./collarData/data/gpsdata_for_lsm.rda")
load("./lsm/outputs/jags_samp_lsm.rda")
load("./lsm/outputs/state_prob_lsm.rda")

gpsdataL   <- split(gpsdata, gpsdata$animal_id)
sorted_ids <- names(gpsdataL)
nid        <- length(sorted_ids)

# 2. Short-ID lookup ----
short_id_fix <- c(
  "Huemul-07"="F7",  "Huemul-1"="F1",  "Huemul-10"="F10",
  "Huemul-11"="M11", "Huemul-12"="M12","Huemul-13"="M13",
  "Huemul-14"="M14", "Huemul-15"="F15","Huemul-2"="F2",
  "Huemul-3" ="F3",  "Huemul-4"="M4",  "Huemul-5"="M5",
  "Huemul-8" ="F8"
)

# 3. Classification rule ----

# Bastille-Rousseau et al. (2016, Table 2), based on the posterior median
# transition probabilities q11, q22 and q33:
#   Resident  : q22 <= 0.90 & q33 <= 0.90
#   Nomadic   : q11 >  0.95 & q22 >  0.90 & q33 <= 0.85
#   Migratory : q11 >  0.95 & q22 >  0.95 & q33 >  0.85, and the animal
#               transitioned back from E2 to E1 at least once
#   Disperser : same thresholds as Migratory, but never transitioned back
#               from E2 to E1
#   Unknown   : any other combination

# has the modal state sequence ever gone from E2 back to E1?
ever_E2_to_E1 <- function(s_E1, s_E2, s_Ex) {
  state_seq <- apply(cbind(s_E1, s_E2, s_Ex), 1, which.max)
  any(state_seq[-1] == 1 & state_seq[-length(state_seq)] == 2)
}

classify_br2016 <- function(q11, q22, q33, reversal) {
  if (q11 > 0.95 & q22 > 0.95 & q33 > 0.85) {
    if (reversal) "Migratory" else "Disperser"
  } else if (q11 > 0.95 & q22 > 0.90 & q33 <= 0.85) "Nomadic"
  else if (q22 <= 0.90 & q33 <= 0.90) "Resident"
  else "Unknown"
}

# 4. Per-individual classification ----
classification_df <- do.call(rbind, lapply(seq_len(nid), function(i) {

  id <- sorted_ids[i]
  tm <- jags_sampL[[id]]$sims.list$trans_mat   # [iter, 3, 3]

  q11_med <- median(tm[, 1, 1])
  q22_med <- median(tm[, 2, 2])
  q33_med <- median(tm[, 3, 3])

  reversal <- ever_E2_to_E1(state_probL[[id]]$s_E1,
                            state_probL[[id]]$s_E2,
                            state_probL[[id]]$s_Ex)

  data.frame(
    animal_id      = id,
    short_id       = unname(short_id_fix[id]),
    q11_median     = q11_med,
    q22_median     = q22_med,
    q33_median     = q33_med,
    classification = classify_br2016(q11_med, q22_med, q33_med, reversal),
    stringsAsFactors = FALSE
  )
}))

# 5. Save output ----
save(classification_df, file = "./lsm/outputs/classification_lsm.rda")

# 6. Summary ----
classification_df$q11_median <- round(classification_df$q11_median, 3)
classification_df$q22_median <- round(classification_df$q22_median, 3)
classification_df$q33_median <- round(classification_df$q33_median, 3)

cat("\n=== Classification (Bastille-Rousseau et al. 2016, posterior medians) ===\n")
print(classification_df, row.names = FALSE)

cat("\n=== Summary counts ===\n")
print(table(classification_df$classification))
