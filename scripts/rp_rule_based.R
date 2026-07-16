# Rule-based reservoir policy

# Based on the decision table structure (Arnold et al., 2018) implemented
# in SWAT+ model (Bieger et al., 2017)

# Structure of the rules is based on Sánchez-Gómez et al., (2025)


library(tidyverse)
library(hydroGOF)
library(DEoptim)
library(lubridate)
library(parallel)
library(gridExtra)
library(patchwork)

#### 1. Input datasets ####

# Reservoir to simulate
studied_reservoirs <- read_table("data/studied_reservoirs.txt")$studied_reservoirs

# Properties (from GDW and adapted considering clean time series)
reservoir_properties <- read_csv("data/Processed_data/reserv_properties.csv")

# Time series for simulations (includes water balance error derived from observations)
reservoirs_time_series <- read_csv("data/Processed_data/time_series_for_simulations.csv")

# Hypsometric curves (used in the evaporation simulation)
hypso_curves <- read_csv("data/Processed_data/hypso_curves.csv")


#### 2. Functions ####

##### 2.1. Common functions ####
###### 2.1.1. Functions used in simulations ######

# Use of hypsometric curve to calculate area from simulated storage
get_area <- function(hypso, storage){
  area <- if(hypso$poly2 == 0){
    hypso$intercept + hypso$poly1 * storage
  } else if(!is.na(hypso$max_stor_rel) && storage > hypso$max_stor_rel){
    hypso$intercept + hypso$poly1 * hypso$max_stor_rel +
      hypso$poly2 * hypso$max_stor_rel^2
  } else{
    hypso$intercept + hypso$poly1 * storage +
      hypso$poly2 * storage^2
  }
  max(area, 0.0005)
}

# Calculating area, mass balance fluxes (P, E, and unknown fluxes extracted from observed data)
mb_fluxes <- function(stor,
                      hypso,
                      pcp,
                      evap,
                      err,
                      pcp_fix,
                      err_fix){
  
  sim_area <- get_area(hypso = hypso, storage = stor)
  
  sim_pcp_Mm3 <- pcp_fix * (sim_area * pcp) / 1e3
  sim_evap_Mm3 <- (sim_area * evap) / 1e3
  
  et_pcp_flux <- sim_pcp_Mm3 - sim_evap_Mm3
  
  wb_err <- if(isTRUE(err_fix)){
    err
  } else {
    0
  }
  
  mb_flux <- et_pcp_flux + wb_err
  
  tibble::new_tibble(
    list(
      area_km = sim_area,
      sim_pcp_Mm3 = sim_pcp_Mm3,
      sim_evap_Mm3 = sim_evap_Mm3,
      mb_flux = mb_flux
    ),
    nrow = 1L
  )
}

# Apply physical constrains to prevent overflow and negative storage
physical_constrains <- function(storage,
                                max_capacity,
                                prov_release){
  spill <- max(0, storage - max_capacity)
  release <- prov_release + spill
  final_storage <- max(0, storage - spill)
  
  tibble::new_tibble(
    list(
      final_release = release,
      final_storage = final_storage
    ),
    nrow = 1L
  )
}

###### 2.1.2. Other functions (helpers) ####

obs_sim_tab <- function(sim_data,
                        obs_data){
  
  sim <- sim_data %>%
    select(date, sim_storage_Mm3, sim_outflow_Mm3)
  
  obs <- obs_data %>%
    select(date, outflow, storage) %>%
    mutate(
      obs_outflow_Mm3 = outflow * 86400 / 1e6,
      obs_storage = storage
    ) %>%
    select(date, obs_outflow_Mm3, obs_storage)
  
  left_join(sim, obs, by = "date")
}

r2_cor <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred)
  if (sum(ok) < 2) return(NA_real_)
  cor(obs[ok], pred[ok])^2
}

##### 2.2. Specific functions ####
###### 2.2.1. Release policy functions

rulebased_prep <- function(grand_id,
                           reservoirs_time_series,
                           reservoir_properties,
                           hypso_curves,
                           full_period = 1990:2014) {
  
  ts_full <- reservoirs_time_series %>%
    filter(GRAND_ID == grand_id, year(date) %in% full_period) %>%
    arrange(date) %>%
    mutate(
      inflow_Mm3_day    = inflow * 86400 / 1e6,
      outflow_Mm3_day   = outflow * 86400 / 1e6,
      evap_mm           = replace_na(ET_mm, 0),
      pcp_mm            = replace_na(PCP_mm, 0),
      error_balance_Mm3 = replace_na(daily_err_Mm3, 0),
      month             = month(date)
    )
  
  if (nrow(ts_full) == 0) {
    stop(paste("No time series found for GRAND_ID", grand_id))
  }
  
  av_inflow <- mean(ts_full$inflow_Mm3_day, na.rm = TRUE)
  
  monthly_median_outflow <- ts_full %>%
    group_by(month) %>%
    summarise(
      median_outflow_Mm3_day = median(outflow_Mm3_day, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    complete(month = 1:12) %>%
    arrange(month)
  
  max_cap <- reservoir_properties %>%
    filter(GRAND_ID == grand_id) %>%
    pull(Cap_mcm)
  
  if (length(max_cap) != 1 || !is.finite(max_cap)) {
    stop(paste("Invalid capacity for GRAND_ID", grand_id))
  }
  
  hypso <- hypso_curves %>%
    filter(GRAND_ID == grand_id)
  
  if (nrow(hypso) != 1) {
    stop(paste("Invalid hypsometric curve for GRAND_ID", grand_id))
  }
  
  prev_year_end_storage <- reservoirs_time_series %>%
    filter(
      GRAND_ID == grand_id,
      year(date) %in% (min(full_period) - 1):(max(full_period) - 1)
    ) %>%
    arrange(date) %>%
    mutate(year = year(date)) %>%
    group_by(year) %>%
    slice_tail(n = 1) %>%
    ungroup() %>%
    select(year, storage)
  
  list(
    grand_id = grand_id,
    ts_full = ts_full,
    av_inflow = av_inflow,
    monthly_median_outflow = monthly_median_outflow$median_outflow_Mm3_day,
    max_cap = max_cap,
    hypso = hypso,
    prev_year_end_storage = prev_year_end_storage
  )
}

rulebased_release <- function(storage_after_inflow,
                              max_cap,
                              month_i,
                              prep,
                              scarcity_threshold,
                              env_flow,
                              emergency_threshold,
                              security_buffer_frac,
                              emergency_lag) {
  
  stor_frac <- max(0, min(1, storage_after_inflow / max_cap))
  
  median_month_outflow <- prep$monthly_median_outflow[month_i]
  if (!is.finite(median_month_outflow)) median_month_outflow <- env_flow
  
  emergency_target <- max(
    0,
    min(max_cap, ((emergency_threshold * max_cap) - (security_buffer_frac * max_cap)))
  )
  
  lag_days <- max(1, round(emergency_lag))
  
  release_Mm3 <- if (stor_frac < scarcity_threshold) {
    env_flow
  } else if (stor_frac < emergency_threshold) {
    median_month_outflow
  } else {
    (storage_after_inflow - emergency_target) / lag_days
  }
  
  release_Mm3 <- max(0, release_Mm3)
  release_Mm3 <- min(release_Mm3, storage_after_inflow)
  
  release_Mm3
}

rulebased_sim <- function(grand_id,
                          res_tseries,
                          sim_period,
                          scarcity_threshold,
                          env_flow_frac,
                          emergency_threshold,
                          security_buffer_frac,
                          emergency_lag,
                          pcp_fix = 0.6,
                          err_fix = TRUE) {
  
  full_series <- res_tseries[res_tseries$GRAND_ID == grand_id, ]
  
  max_cap <- reservoir_properties[reservoir_properties$GRAND_ID == grand_id, "Cap_mcm"] %>%
    pull(.)
  
  hypso_curve <- hypso_curves[hypso_curves$GRAND_ID == grand_id, ]
  
  prep <- rulebased_prep(
    grand_id = grand_id,
    reservoirs_time_series = res_tseries,
    reservoir_properties = reservoir_properties,
    hypso_curves = hypso_curves,
    full_period = 1990:2014
  )
  
  min_year <- min(sim_period)
  
  initial_stor <- prep$prev_year_end_storage %>%
    filter(year == (min_year - 1)) %>%
    pull(storage)
  
  initial_stor <- if (length(initial_stor) == 0 || is.na(initial_stor)) {
    full_series[year(full_series$date) == min_year & !is.na(full_series$storage), ] %>%
      slice_head(n = 1) %>%
      pull(storage)
  } else {
    initial_stor
  }
  
  tseries <- prep$ts_full[year(prep$ts_full$date) %in% sim_period, ]
  tsteps <- nrow(tseries)
  
  sim_storage <- numeric(tsteps)
  sim_outflow <- numeric(tsteps)
  sim_area <- numeric(tsteps)
  sim_pcp <- numeric(tsteps)
  sim_evap <- numeric(tsteps)
  
  sim_series <-
    tibble(
      date = tseries$date,
      inflow_Mm3_day = tseries$inflow_Mm3_day,
      evap_mm = tseries$evap_mm,
      pcp_mm = tseries$pcp_mm,
      error_balance_Mm3 = tseries$error_balance_Mm3
    ) %>%
    mutate(
      sim_storage_Mm3 = NA_real_,
      sim_outflow_Mm3 = NA_real_,
      sim_area_km2 = NA_real_,
      sim_evap_Mm3 = NA_real_,
      sim_pcp_Mm3 = NA_real_,
      GRAND_ID = grand_id
    ) %>%
    mutate(error_balance_Mm3 = if_else(is.na(error_balance_Mm3), 0, error_balance_Mm3))
  
  env_flow <- env_flow_frac * prep$av_inflow
  
  for (tstp in seq_len(tsteps)) {
    
    storage_after_inflow <- if (tstp == 1) {
      initial_stor + sim_series$inflow_Mm3_day[tstp]
    } else {
      sim_storage[tstp - 1] + sim_series$inflow_Mm3_day[tstp]
    }
    
    month_i <- month(sim_series$date[tstp])
    
    release_Mm3 <- rulebased_release(
      storage_after_inflow = storage_after_inflow,
      max_cap = max_cap,
      month_i = month_i,
      prep = prep,
      scarcity_threshold = scarcity_threshold,
      env_flow = env_flow,
      emergency_threshold = emergency_threshold,
      security_buffer_frac = security_buffer_frac,
      emergency_lag = emergency_lag
    )
    
    storage_after_release <- max(0, storage_after_inflow - release_Mm3)
    
    mb_flux <- mb_fluxes(
      stor = storage_after_release,
      hypso = hypso_curve,
      pcp = sim_series$pcp_mm[tstp],
      evap = sim_series$evap_mm[tstp],
      err = sim_series$error_balance_Mm3[tstp],
      pcp_fix = pcp_fix,
      err_fix = err_fix
    )
    
    storage_pre_constr <- storage_after_release + mb_flux$mb_flux
    
    constrains <- physical_constrains(
      storage = storage_pre_constr,
      max_capacity = max_cap,
      prov_release = release_Mm3
    )
    
    sim_area[tstp] <- mb_flux$area_km
    sim_pcp[tstp] <- mb_flux$sim_pcp_Mm3
    sim_evap[tstp] <- mb_flux$sim_evap_Mm3
    sim_storage[tstp] <- constrains$final_storage
    sim_outflow[tstp] <- constrains$final_release
  }
  
  sim_series$sim_area_km2 <- sim_area
  sim_series$sim_pcp_Mm3 <- sim_pcp
  sim_series$sim_evap_Mm3 <- sim_evap
  sim_series$sim_storage_Mm3 <- sim_storage
  sim_series$sim_outflow_Mm3 <- sim_outflow
  
  sim_series
}

###### 2.2.2. Functions used for optimization ######

fitting_rulebased <- function(grand_id, sim_period){
  
  obs_data <- reservoirs_time_series %>%
    filter(GRAND_ID == grand_id,
           year(date) %in% sim_period)
  
  obj_fn_rulebased <- function(params) {
    
    scarcity_threshold   <- params[1]
    env_flow_frac        <- params[2]
    emergency_threshold  <- params[3]
    security_buffer_frac <- params[4]
    emergency_lag        <- params[5]
    
    sim <- tryCatch(
      rulebased_sim(
        grand_id = grand_id,
        res_tseries = obs_data,
        sim_period = sim_period,
        scarcity_threshold = scarcity_threshold,
        env_flow_frac = env_flow_frac,
        emergency_threshold = emergency_threshold,
        security_buffer_frac = security_buffer_frac,
        emergency_lag = emergency_lag
      ),
      error = function(e) NULL
    )
    
    if (is.null(sim)) return(1e6)
    
    tab <- tryCatch(
      obs_sim_tab(sim_data = sim,
                  obs_data = obs_data),
      error = function(e) NULL
    )
    
    if (is.null(tab)) return(1e6)
    
    ok_out <- is.finite(tab$sim_outflow_Mm3) & is.finite(tab$obs_outflow_Mm3)
    ok_sto <- is.finite(tab$sim_storage_Mm3) & is.finite(tab$obs_storage)
    
    if (sum(ok_out) < 2 || sum(ok_sto) < 2) return(1e6)
    
    nrmse_out <- hydroGOF::nrmse(
      sim = tab$sim_outflow_Mm3[ok_out],
      obs = tab$obs_outflow_Mm3[ok_out],
      na.rm = TRUE,
      norm = "sd"
    )
    
    nrmse_sto <- hydroGOF::nrmse(
      sim = tab$sim_storage_Mm3[ok_sto],
      obs = tab$obs_storage[ok_sto],
      na.rm = TRUE,
      norm = "sd"
    )
    
    0.5 * nrmse_out + 0.5 * nrmse_sto
  }
  
  fit <- DEoptim(
    fn = obj_fn_rulebased,
    lower = low_bound,
    upper = upp_bound,
    control = DEoptim_conf
  )
  
  best_pars <- tibble(
    !!parms[1] := fit$optim$bestmem[1],
    !!parms[2] := fit$optim$bestmem[2],
    !!parms[3] := fit$optim$bestmem[3],
    !!parms[4] := fit$optim$bestmem[4],
    !!parms[5] := fit$optim$bestmem[5]
  )
  
  perf_evol <- tibble(
    iteration = seq_along(fit$member$bestvalit),
    obj_f = fit$member$bestvalit
  )
  
  best_ind <- tibble(
    iteration = seq_len(nrow(fit$member$bestmemit)),
    !!parms[1] := fit$member$bestmemit[,1],
    !!parms[2] := fit$member$bestmemit[,2],
    !!parms[3] := fit$member$bestmemit[,3],
    !!parms[4] := fit$member$bestmemit[,4],
    !!parms[5] := fit$member$bestmemit[,5]
  )
  
  all_pop <- purrr::map_dfr(seq_along(fit$member$storepop), function(i) {
    tt <- tibble::as_tibble(fit$member$storepop[[i]])
    colnames(tt) <- parms
    dplyr::mutate(tt, iteration = i, .before = 1)
  })
  
  list(
    best_pars = best_pars,
    perf_evol = perf_evol,
    best_ind = best_ind,
    all_pop = all_pop
  )
}


#### 3. Running simulations for one reservoir ####

##### 3.1. Default simulation #####

gid <- 132

# Observed data
obs_data <- reservoirs_time_series %>%
  filter(GRAND_ID == gid)

# Simulation parameters
parms <- c(
  "scarcity_threshold",
  "env_flow_frac",
  "emergency_threshold",
  "security_buffer_frac",
  "emergency_lag"
)

def_parms <- tibble(
  scarcity_threshold = 0.10,
  env_flow_frac = 0.10,
  emergency_threshold = 0.9,
  security_buffer_frac = 0.05,
  emergency_lag = 10
)

# Run simulation
def_sim <- rulebased_sim(
  grand_id = gid,
  res_tseries = reservoirs_time_series,
  sim_period = 1990:2014,
  scarcity_threshold = def_parms$scarcity_threshold,
  env_flow_frac = def_parms$env_flow_frac,
  emergency_threshold = def_parms$emergency_threshold,
  security_buffer_frac = def_parms$security_buffer_frac,
  emergency_lag = def_parms$emergency_lag
)

# Obs.sim. table
obssim_tab <- obs_sim_tab(
  sim_data = def_sim,
  obs_data = obs_data
)

# Performance
obssim_tab %>%
  summarise(
    KGE_stor = KGE(sim_storage_Mm3, obs_storage),
    KGE_out = KGE(sim_outflow_Mm3, obs_outflow_Mm3)
  )

# Plot
obssim_tab %>%
  pivot_longer(-date) %>%
  mutate(
    var = if_else(str_detect(name, "out"), "Release (Mm3/d)", "Storage (Mm3)"),
    type = if_else(str_detect(name, "sim"), "Simulated", "Observed")
  ) %>%
  ggplot(aes(x = date, y = value, color = type)) +
  geom_line() +
  facet_wrap(~var, scales = "free", ncol = 1) +
  scale_color_manual(values = c("black", "steelblue")) +
  theme_bw() +
  theme(axis.title = element_blank(), legend.title = element_blank())


##### 3.2. Parameters optimization #####

gid <- 132

# Observed data
obs_data <- reservoirs_time_series %>%
  filter(GRAND_ID == gid)

# Simulation parameters
parms <- c(
  "scarcity_threshold",
  "env_flow_frac",
  "emergency_threshold",
  "security_buffer_frac",
  "emergency_lag"
)

def_parms <- tibble(
  scarcity_threshold = 0.10,
  env_flow_frac = 0.10,
  emergency_threshold = 0.75,
  security_buffer_frac = 0.05,
  emergency_lag = 10
)

n_pars <- length(parms)
max_itrs <- 200

# Bounds for the five parameters
low_bound <- c(
  0.00,  # scarcity_threshold
  0.05,  # env_flow_frac
  0.80,  # emergency_threshold
  0.01,  # security_buffer_frac
  1.00   # emergency_lag
)

upp_bound <- c(
  0.30,  # scarcity_threshold
  0.60,  # env_flow_frac
  0.99,  # emergency_threshold
  0.10,  # security_buffer_frac
  15.00  # emergency_lag
)

# DE adjustments --> NP, F, CR default
DEoptim_conf <- DEoptim.control(
  NP = 10 * n_pars,
  itermax = max_itrs,
  F = 0.8,
  CR = 0.9,
  trace = FALSE,
  reltol = 1e-4,
  steptol = 10,
  storepopfrom = 1,
  storepopfreq = 1
)

# Run optimization
deoptim_fit <- fitting_rulebased(
  grand_id = gid,
  sim_period = 1990:2004
)


##### 3.3. Optimized simulation #####

opt_pars <- deoptim_fit$best_pars

# Run simulation
opt_sim <- rulebased_sim(
  grand_id = gid,
  res_tseries = reservoirs_time_series,
  sim_period = 1990:2014,
  scarcity_threshold = opt_pars$scarcity_threshold,
  env_flow_frac = opt_pars$env_flow_frac,
  emergency_threshold = opt_pars$emergency_threshold,
  security_buffer_frac = opt_pars$security_buffer_frac,
  emergency_lag = opt_pars$emergency_lag
)

# Obs.sim. table
obssim_tab <- obs_sim_tab(
  sim_data = opt_sim,
  obs_data = obs_data
)

# Performance
obssim_tab %>%
  summarise(
    KGE_stor = KGE(sim_storage_Mm3, obs_storage),
    KGE_out = KGE(sim_outflow_Mm3, obs_outflow_Mm3)
  )

# Plot
obssim_tab %>%
  pivot_longer(-date) %>%
  mutate(
    var = if_else(str_detect(name, "out"), "Release (Mm3/d)", "Storage (Mm3)"),
    type = if_else(str_detect(name, "sim"), "Simulated", "Observed")
  ) %>%
  ggplot(aes(x = date, y = value, color = type)) +
  geom_line() +
  facet_wrap(~var, scales = "free", ncol = 1) +
  scale_color_manual(values = c("black", "steelblue")) +
  theme_bw() +
  theme(axis.title = element_blank(), legend.title = element_blank())

# Compare fitted vs default
rbind(
  (opt_sim %>%
     select(date, sim_storage_Mm3, sim_outflow_Mm3) %>%
     mutate(pars = "fitted")),
  (def_sim %>%
     select(date, sim_storage_Mm3, sim_outflow_Mm3) %>%
     mutate(pars = "default"))
) %>%
  left_join(., obs_data[, c("date", "outflow", "storage")]) %>%
  mutate(
    obs_storage = storage,
    obs_outflow_Mm3 = outflow * 86400 / 1e6
  ) %>%
  select(date:pars, starts_with("obs")) %>%
  pivot_longer(-c(date, pars)) %>%
  mutate(
    var = if_else(str_detect(name, "out"), "Release (Mm3/d)", "Storage (Mm3)"),
    type = if_else(str_detect(name, "sim"), "Simulated", "Observed")
  ) %>%
  mutate(
    type = case_when(
      type == "Observed" ~ "Observed",
      pars == "default" ~ "Sim. default",
      .default = "Sim. fitted"
    )
  ) %>%
  ggplot(aes(x = date, y = value, color = type)) +
  geom_line() +
  facet_wrap(~var, scales = "free", ncol = 1) +
  scale_color_manual(values = c("black", "steelblue", "red4")) +
  theme_bw() +
  theme(axis.title = element_blank(), legend.title = element_blank())


