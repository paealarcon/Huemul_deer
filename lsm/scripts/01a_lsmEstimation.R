# =============================================================================
# Movement behavior and habitat selection to guide conservation of an
# endangered Neotropical deer
#
# Alarcón, P.A.E., Crego, R.D., Galetto, E., Aprile, G.,
# Morales, J.M., Pauli, J.N. & Donadio, E.
#
# -----------------------------------------------------------------------------
# Script:  01a_lsmEstimation.R
# Purpose: Fit the Bayesian latent-state model (LSM) to individual net
#          squared displacement (NSD) time series, following
#          Bastille-Rousseau et al. (2016) Mov Ecol
# =============================================================================

# 0. Libraries ----
library(jagsUI)

# 1. Data ----
load("./collarData/data/gpsdata_for_lsm.rda")

gpsdataL <- split(gpsdata, gpsdata$animal_id)
nid      <- length(gpsdataL)

# 2. BUG model ----

# Mixture of two normal distributions (encamped states E1, E2) and one
# pseudo-uniform distribution (exploratory state Ex). State switching
# follows a first-order Markov process defined by the 3x3 transition
# matrix Q.

bug_file <- "./lsm/scripts/01_lsm.bug"

cat(
  file = bug_file,

"# =============================================================================
# Movement behavior and habitat selection to guide conservation of an
# endangered Neotropical deer
#
# Alarcón, P.A.E., Crego, R.D., Galetto, E., Aprile, G.,
# Morales, J.M., Pauli, J.N. & Donadio, E.
#
# -----------------------------------------------------------------------------
# Model:   01_lsm.bug
# Purpose: JAGS specification of the three-state Latent-State Markov model
#          (two encamped states E1, E2 and one exploratory state Ex) used to
#          classify huemul movement strategies, following Bastille-Rousseau
#          et al. (2016). Priors are flat Dirichlet for the initial state
#          distribution and each row of the transition matrix.
# =============================================================================

 model {

   # ---- Emission means ----
   # mu[1]: low-NSD encamped state (E1)
   # mu[2] = mu[1] + eps: high-NSD encamped state (E2)
   mu[1] ~ dunif(0, 0.50)
   eps   ~ dunif(0.05, 0.50)
   mu[2] <- mu[1] + eps
   mu[3] <- 0.5

   # ---- Emission standard deviations ----
   sigma[1] ~ dunif(0.01, 0.5)
   sigma[2] ~ dunif(0.01, 0.5)
   sigma[3] <- 10

   tau[1] <- pow(sigma[1], -2)
   tau[2] <- pow(sigma[2], -2)
   tau[3] <- pow(sigma[3], -2)

   # ---- Initial state distribution: flat Dirichlet prior ----
   alpha1[1] <- 1
   alpha1[2] <- 1
   alpha1[3] <- 1
   p_state[1:3] ~ ddirch(alpha1[])
   idx[1] ~ dcat(p_state[])

   # ---- Transition matrix Q: flat Dirichlet prior per row ----
   alpha[1] <- 1
   alpha[2] <- 1
   alpha[3] <- 1
   trans_mat[1,1:3] ~ ddirch(alpha[])
   trans_mat[2,1:3] ~ ddirch(alpha[])
   trans_mat[3,1:3] ~ ddirch(alpha[])

   # ---- State sequence, emissions and state indicators ----
   s_E1[1] <- equals(idx[1], 1)
   s_E2[1] <- equals(idx[1], 2)
   s_Ex[1] <- equals(idx[1], 3)

   for (t in 2:ndays) {
     idx[t]     ~ dcat(trans_mat[idx[t-1], ])
     nsd_scl[t] ~ dnorm(mu[idx[t]], tau[idx[t]])
     s_E1[t]  <- equals(idx[t], 1)
     s_E2[t]  <- equals(idx[t], 2)
     s_Ex[t]  <- equals(idx[t], 3)
   }
  }"
)

# 3. Fit model per individual ----
jags_sampL <- vector("list", nid)
names(jags_sampL) <- names(gpsdataL)

params <- c("mu", "eps", "sigma", "trans_mat", "p_state",
            "s_E1", "s_E2", "s_Ex")

nc <- 3       # chains
ni <- 100000  # iterations
nt <- 4       # thinning rate
nb <- 20000   # burn-in

for (h in seq_len(nid)) {

  nsd_scl <- gpsdataL[[h]]$nsd_scl
  ndays   <- length(nsd_scl)
  data    <- list(nsd_scl = nsd_scl, ndays = ndays)

  inits <- function() {
    list(
      mu    = c(runif(1, 0.08, 0.18), NA, NA),
      eps   = runif(1, 0.25, 0.45),
      sigma = c(runif(1, 0.02, 0.06), runif(1, 0.04, 0.12), NA)
    )
  }

  jags_sampL[[h]] <- jags(data = data,
                          inits = inits,
                          parameters.to.save = params,
                          n.chains = nc,
                          n.thin = nt,
                          n.iter = ni,
                          n.burnin = nb,
                          model.file = bug_file,
                          parallel = FALSE)
  print(h)
}

save(jags_sampL, file = "./lsm/outputs/jags_samp_lsm.rda")

# 4. State probabilities ----
# Posterior mean of each state indicator (s_E1, s_E2, s_Ex), per day.
state_probL <- lapply(seq_len(nid), function(h) {
  list(
    s_E1 = jags_sampL[[h]]$mean$s_E1,
    s_E2 = jags_sampL[[h]]$mean$s_E2,
    s_Ex = jags_sampL[[h]]$mean$s_Ex,
    day  = gpsdataL[[h]]$day,
    Date = gpsdataL[[h]]$Date
  )
})
names(state_probL) <- names(gpsdataL)

save(state_probL, file = "./lsm/outputs/state_prob_lsm.rda")

# 5. MCMC convergence diagnostics ----
# Effective sample size (n.eff) and Gelman-Rubin statistic (Rhat) for the
# key stochastic nodes: mu[1], eps, sigma, and trans_mat
load("./collarData/data/gpsdata_for_lsm.rda")
load("./lsm/outputs/jags_samp_lsm.rda")

gpsdataL <- split(gpsdata, gpsdata$animal_id)
nid      <- length(gpsdataL)

per_mat <- as.data.frame(matrix(NA, nid, 10))
colnames(per_mat) <- c("neff.mu1", "neff.eps", "neff.sigma1",
                       "neff.sigma2", "neff.trans_mat.min",
                       "rhat.mu1", "rhat.eps", "rhat.sigma1",
                       "rhat.sigma2", "rhat.trans_mat.max")
per_mat$animal_id <- names(gpsdataL)

for (h in seq_len(nid)) {

  per_mat[h, 1:5] <- c(
    jags_sampL[[h]]$n.eff$mu[1],
    jags_sampL[[h]]$n.eff$eps,
    jags_sampL[[h]]$n.eff$sigma[1:2],
    min(jags_sampL[[h]]$n.eff$trans_mat)
  )

  per_mat[h, 6:10] <- round(c(
    jags_sampL[[h]]$Rhat$mu[1],
    jags_sampL[[h]]$Rhat$eps,
    jags_sampL[[h]]$Rhat$sigma[1:2],
    max(jags_sampL[[h]]$Rhat$trans_mat, na.rm = TRUE)
  ), 3)
}

save(per_mat, file = "./lsm/outputs/performace_lsm.rda")

# 6. Transition matrices ----
# Posterior median transition matrix (3x3) for each individual.
trans_matL  <- list()
state_names <- c("encamped1", "encamped2", "exploration")

for (h in seq_len(nid)) {
  trans_mat_samp   <- jags_sampL[[h]]$sims.list$trans_mat  # [iter, 3, 3]
  trans_matL[[h]]  <- apply(trans_mat_samp, c(2, 3), median)
  colnames(trans_matL[[h]]) <- state_names
  rownames(trans_matL[[h]]) <- state_names
}
names(trans_matL) <- names(gpsdataL)

save(trans_matL, file = "./lsm/outputs/trans_matL.rda")
