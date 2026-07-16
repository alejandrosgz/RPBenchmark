# SDP reservoir policy framework
# Monthly SDP solver + daily forward simulation

# Functions include modifications, but their core have been extracted from 
# reservoir package, specifically reservoir::sdp_multi()
# Reference for original function is Turner and Galelli (2016)

library(tidyverse)
library(hydroGOF)
library(DEoptim)
library(lubridate)
library(zoo)
library(parallel)
library(gridExtra)
library(patchwork)

setwd('C:/ASG/UCDavis/Research/Papers/RPB/')

#### 1. Input datasets ####

# Reservoirs to simulate
studied_reservoirs <- read_table("data/studied_reservoirs.txt")$studied_reservoirs

# Reservoir properties (capacity, area, etc.)
reservoir_properties <- read_csv("data/Processed_data/reserv_properties.csv")

# Daily time series used for simulations and evaluation
reservoirs_time_series <- read_csv("data/Processed_data/time_series_for_simulations.csv")

# Hypsometric curves (used to compute area from storage)
hypso_curves <- read_csv("data/Processed_data/hypso_curves.csv")

# Water demand by reservoir and month (required by the SDP policy)
# Expected columns: GRAND_ID, Month, Total_demand_res_Mm3_month
wat_demand <- read_csv("data/Processed_data/reservoirs_demands_month.csv")

#### 2. Functions ####

##### 2.1. Common functions #####

# Use of hypsometric curve to calculate area from simulated storage
# Storage is in Mm3, area in km2.
get_area <- function(hypso, storage) {
  area <- if (hypso$poly2 == 0) {
    hypso$intercept + hypso$poly1 * storage
  } else if (!is.na(hypso$max_stor_rel) && storage > hypso$max_stor_rel) {
    hypso$intercept + hypso$poly1 * hypso$max_stor_rel +
      hypso$poly2 * hypso$max_stor_rel^2
  } else {
    hypso$intercept + hypso$poly1 * storage +
      hypso$poly2 * storage^2
  }
  max(area, 0.0005) # Avoid zero area
}

# Calculating area, mass balance fluxes (P, E, and unknown fluxes extracted from observed data)
# pcp and evap are in mm for the time step.
# Output fluxes are in Mm3 for the time step.
mb_fluxes <- function(stor,
                      hypso,
                      pcp,
                      evap,
                      err,
                      pcp_fix,
                      err_fix) {
  
  sim_area <- get_area(hypso = hypso, storage = stor)
  
  # mm * km2 -> Mm3
  sim_pcp_Mm3 <- pcp_fix * (sim_area * pcp) / 1e3
  sim_evap_Mm3 <- (sim_area * evap) / 1e3
  
  et_pcp_flux <- sim_pcp_Mm3 - sim_evap_Mm3
  
  wb_err <- if (isTRUE(err_fix)) err else 0
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

# Apply physical constraints to prevent overflow and negative storage
physical_constrains <- function(storage,
                                max_capacity,
                                prov_release) {
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

# Build observed vs simulated table at daily scale
obs_sim_tab <- function(sim_data,
                        obs_data) {
  
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

# Build observed vs simulated table at monthly scale
obs_sim_tab_monthly <- function(sim_data,
                                obs_data) {
  
  sim_m <- sim_data %>%
    mutate(ym = floor_date(date, "month")) %>%
    group_by(ym) %>%
    summarise(
      sim_storage_Mm3 = mean(sim_storage_Mm3, na.rm = TRUE),
      sim_outflow_Mm3 = sum(sim_outflow_Mm3, na.rm = TRUE),
      .groups = "drop"
    )
  
  obs_m <- obs_data %>%
    mutate(
      obs_outflow_Mm3 = outflow * 86400 / 1e6,
      ym = floor_date(date, "month")
    ) %>%
    group_by(ym) %>%
    summarise(
      obs_storage_Mm3 = mean(storage, na.rm = TRUE),
      obs_outflow_Mm3 = sum(obs_outflow_Mm3, na.rm = TRUE),
      .groups = "drop"
    )
  
  left_join(sim_m, obs_m, by = "ym") %>%
    arrange(ym)
}

# Coefficient of determination
r2_cor <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred)
  if (sum(ok) < 2) return(NA_real_)
  cor(obs[ok], pred[ok])^2
}

# Aggregate daily series to monthly inputs for the solver
build_monthly_inputs <- function(daily_df) {
  daily_df %>%
    arrange(date) %>%
    mutate(
      inflow_Mm3_day = inflow * 86400 / 1e6,
      evap_mm = replace_na(ET_mm, 0),
      pcp_mm = replace_na(PCP_mm, 0),
      error_balance_Mm3 = replace_na(daily_err_Mm3, 0),
      year = year(date),
      month = month(date)
    ) %>%
    group_by(year, month) %>%
    summarise(
      inflow_month_Mm3 = sum(inflow_Mm3_day, na.rm = TRUE),
      evap_month_mm = sum(evap_mm, na.rm = TRUE),
      pcp_month_mm = sum(pcp_mm, na.rm = TRUE),
      error_month_Mm3 = sum(error_balance_Mm3, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(year, month)
}

# Build a seasonal monthly target vector from the demand table
# IMPORTANT: demand is monthly Mm3/month; it is not divided by days here.
build_monthly_demand_target <- function(grand_id,
                                        demand_table,
                                        demand_adj = 1) {
  x <- demand_table %>%
    filter(GRAND_ID == grand_id) %>%
    select(Month, Total_demand_res_Mm3_month) %>%
    mutate(Month = as.integer(Month)) %>%
    arrange(Month) %>%
    tidyr::complete(Month = 1:12)
  
  if (all(is.na(x$Total_demand_res_Mm3_month))) {
    stop(paste0("No demand data found for GRAND_ID ", grand_id))
  }
  
  if (any(is.na(x$Total_demand_res_Mm3_month))) {
    x$Total_demand_res_Mm3_month <- zoo::na.approx(
      x$Total_demand_res_Mm3_month,
      x = x$Month,
      na.rm = FALSE,
      rule = 2
    )
    if (any(is.na(x$Total_demand_res_Mm3_month))) {
      x$Total_demand_res_Mm3_month[is.na(x$Total_demand_res_Mm3_month)] <-
        mean(x$Total_demand_res_Mm3_month, na.rm = TRUE)
    }
  }
  
  as.numeric(x$Total_demand_res_Mm3_month) * demand_adj
}

# Pick initial storage from the year before the simulation period
get_initial_storage <- function(full_series, sim_period) {
  min_year <- min(sim_period)
  
  initial_stor <- full_series %>%
    filter(year(date) == min_year - 1) %>%
    arrange(date) %>%
    slice_tail(n = 1) %>%
    pull(storage)
  
  if (length(initial_stor) == 0 || is.na(initial_stor)) {
    initial_stor <- full_series %>%
      filter(year(date) == min_year, !is.na(storage)) %>%
      arrange(date) %>%
      slice_head(n = 1) %>%
      pull(storage)
  }
  
  initial_stor
}

# Create a compact plot + performance table
make_sdp_plot <- function(sim_df,
                          obs_df,
                          res_name,
                          cap_mcm,
                          title_suffix = "") {
  
  tab <- sim_df %>%
    select(date, sim_storage_Mm3, sim_outflow_Mm3) %>%
    left_join(
      obs_df %>%
        select(date, outflow, storage) %>%
        mutate(
          obs_outflow_Mm3 = outflow * 86400 / 1e6,
          obs_storage_Mm3 = storage
        ) %>%
        select(date, obs_outflow_Mm3, obs_storage_Mm3),
      by = "date"
    ) %>%
    summarise(
      KGE_out = KGE(sim_outflow_Mm3, obs_outflow_Mm3),
      KGE_sto = KGE(sim_storage_Mm3, obs_storage_Mm3),
      R2_out = r2_cor(sim_outflow_Mm3, obs_outflow_Mm3),
      R2_sto = r2_cor(sim_storage_Mm3, obs_storage_Mm3)
    ) %>%
    mutate(across(everything(), ~ round(.x, 2)))
  
  p <- sim_df %>%
    left_join(
      obs_df %>%
        transmute(
          date,
          obs_outflow_Mm3 = outflow * 86400 / 1e6,
          obs_storage_Mm3 = storage
        ),
      by = "date"
    ) %>%
    pivot_longer(cols = c(sim_outflow_Mm3, obs_outflow_Mm3, sim_storage_Mm3, obs_storage_Mm3)) %>%
    mutate(
      type = if_else(str_starts(name, "sim"), "Simulated", "Observed"),
      var = if_else(str_detect(name, "outflow"), "Release (Mm3/d)", "Storage (Mm3)")
    ) %>%
    ggplot(aes(x = date, y = value, color = type)) +
    geom_line(linewidth = 0.6) +
    facet_wrap(~var, scales = "free", ncol = 1) +
    scale_color_manual(values = c("Observed" = "black", "Simulated" = "red")) +
    labs(
      title = paste0(res_name, " (", round(cap_mcm), " Mm3) - ", title_suffix),
      x = "Date"
    ) +
    theme_bw() +
    theme(
      axis.title = element_blank(),
      legend.title = element_blank(),
      text = element_text(size = 12)
    )
  
  list(plot = p, perf = tab)
}

##### 2.2. Monthly SDP policy solver ######

# Adapted from reservoir::sdp_multi(). Turner and Galelli (2016)
sdp_policy_solver_monthly <- function(
    Q,
    capacity,
    target,
    hypso,
    pcp = NULL,
    evap = NULL,
    wb_err = NULL,
    R_max = NULL,
    spill_targ = 0.95, # par
    vol_targ = 0.75, # par
    weights = c(0.7, 0.2, 0.1), #pars
    S_disc = 500,
    R_disc = 100,
    Q_disc = c(0, 0.2375, 0.475, 0.7125, 0.95, 1),
    loss_exp = c(2, 2, 2),
    S_initial = 1,
    plot = FALSE,
    tol = 0.99,
    rep_rrv = FALSE,
    pcp_fix = 0.6
) {
  frq <- frequency(Q)
  if (!is.ts(Q)) stop("Q must be a seasonal time series object")
  if (frq != 12) stop("This script uses a monthly SDP solver (frequency = 12)")
  
  if (missing(target)) stop("target must be provided")
  
  if (length(target) == 1) {
    targ_seas <- rep(as.numeric(target), frq)
    target_ts <- ts(rep(as.numeric(target), length(Q)), start = start(Q), frequency = frq)
  } else if (length(target) == frq && !is.ts(target)) {
    targ_seas <- as.numeric(target)
    target_ts <- ts(rep(targ_seas, length(Q) / frq), start = start(Q), frequency = frq)
  } else if (is.ts(target)) {
    target_ts <- window(target, start = start(Q), end = end(Q), frequency = frq)
    targ_seas <- as.vector(tapply(target_ts, cycle(target_ts), mean, na.rm = TRUE))
  } else {
    stop("target must be a scalar, a vector of length 12, or a ts aligned with Q")
  }
  
  if (missing(hypso) || nrow(hypso) != 1) {
    stop("hypso must be a single hypsometric curve row for the reservoir")
  }
  
  if (missing(pcp) || is.null(pcp)) {
    pcp <- ts(rep(0, length(Q)), start = start(Q), frequency = frq)
  }
  if (missing(evap)) {
    evap <- ts(rep(0, length(Q)), start = start(Q), frequency = frq)
  }
  if (missing(wb_err) || is.null(wb_err)) {
    wb_err <- ts(rep(0, length(Q)), start = start(Q), frequency = frq)
  }
  
  if (length(pcp) == 1) pcp <- ts(rep(pcp, length(Q)), start = start(Q), frequency = frq)
  if (length(evap) == 1) evap <- ts(rep(evap, length(Q)), start = start(Q), frequency = frq)
  if (length(wb_err) == 1) wb_err <- ts(rep(wb_err, length(Q)), start = start(Q), frequency = frq)
  
  if (length(pcp) != length(Q) && length(pcp) != frq) stop("pcp must be a time series of length Q, a length-12 vector, or a scalar")
  if (length(evap) != length(Q) && length(evap) != frq) stop("evap must be a time series of length Q, a length-12 vector, or a scalar")
  if (length(wb_err) != length(Q) && length(wb_err) != frq) stop("wb_err must be a time series of length Q, a length-12 vector, or a scalar")
  
  # Trim incomplete first/last years
  if (start(Q)[2] != 1) {
    message("NOTE: First incomplete year of time series removed")
    Q <- window(Q, start = c(start(Q)[1] + 1, 1), frequency = frq)
    pcp <- window(pcp, start = c(start(pcp)[1] + 1, 1), frequency = frq)
    evap <- window(evap, start = c(start(evap)[1] + 1, 1), frequency = frq)
    wb_err <- window(wb_err, start = c(start(wb_err)[1] + 1, 1), frequency = frq)
    target_ts <- window(target_ts, start = c(start(target_ts)[1] + 1, 1), frequency = frq)
  }
  if (end(Q)[2] != frq) {
    message("NOTE: Final incomplete year of time series removed")
    Q <- window(Q, end = c(end(Q)[1] - 1, frq), frequency = frq)
    pcp <- window(pcp, end = c(end(pcp)[1] - 1, frq), frequency = frq)
    evap <- window(evap, end = c(end(evap)[1] - 1, frq), frequency = frq)
    wb_err <- window(wb_err, end = c(end(wb_err)[1] - 1, frq), frequency = frq)
    target_ts <- window(target_ts, end = c(end(target_ts)[1] - 1, frq), frequency = frq)
  }
  
  if (length(pcp) == frq) pcp <- ts(rep(as.numeric(pcp), length(Q) / frq), start = start(Q), frequency = frq)
  if (length(evap) == frq) evap <- ts(rep(as.numeric(evap), length(Q) / frq), start = start(Q), frequency = frq)
  if (length(wb_err) == frq) wb_err <- ts(rep(as.numeric(wb_err), length(Q) / frq), start = start(Q), frequency = frq)
  
  if (any(is.na(pcp))) pcp[is.na(pcp)] <- 0
  if (any(is.na(evap))) evap[is.na(evap)] <- 0
  if (any(is.na(wb_err))) wb_err[is.na(wb_err)] <- 0
  pcp[pcp < 0] <- 0
  evap[evap < 0] <- 0
  
  if (length(loss_exp) == 1) loss_exp <- rep(loss_exp, 3)
  if (length(loss_exp) != 3) stop("loss_exp must be scalar or length-3 vector")
  if (length(weights) == 1) weights <- rep(weights, 3)
  if (length(weights) != 3) stop("weights must be length 3")
  
  evap_seas <- as.vector(tapply(evap, cycle(evap), FUN = mean, na.rm = TRUE))
  pcp_seas  <- as.vector(tapply(pcp,  cycle(pcp),  FUN = mean, na.rm = TRUE))
  err_seas  <- as.vector(tapply(wb_err, cycle(wb_err), FUN = mean, na.rm = TRUE))
  
  # Use the same storage -> area relationship as the daily simulator
  S_states <- seq(from = 0, to = capacity, by = capacity / S_disc)
  S_area_rel <- vapply(S_states, function(s) get_area(hypso = hypso, storage = s), numeric(1))
  
  Q_month_mat <- matrix(Q, byrow = TRUE, ncol = frq)
  Q.probs <- diff(Q_disc)
  Q_class_med <- apply(Q_month_mat, 2, quantile, type = 8,
                       probs = Q_disc[-1] - (Q.probs / 2))
  
  if (is.null(R_max)) {
    R_max <- 2 * max(targ_seas, na.rm = TRUE)
  }
  R_disc_x <- seq(from = 0, to = R_max, by = R_max / R_disc)
  
  Shell.array <- array(0, dim = c(length(S_states), length(R_disc_x), length(Q.probs)))
  Cost_to_go <- vector("numeric", length = length(S_states))
  R_policy <- matrix(0, nrow = length(S_states), ncol = frq)
  Bellman <- R_policy
  R_policy_test <- R_policy
  
  message(paste0("policy converging... (>", tol, ")"))
  
  repeat {
    for (t in frq:1) {
      
      # Monthly net supply potential using the same flux equations as mb_fluxes()
      # but vectorized over the storage grid.
      monthly_area <- S_area_rel
      
      pcp_vol  <- pcp_fix * (monthly_area * pcp_seas[t]) / 1e3
      evap_vol <- (monthly_area * evap_seas[t]) / 1e3
      mb_term  <- pcp_vol - evap_vol + err_seas[t]
      
      R.cstr <- sweep(Shell.array, 3, Q_class_med[, t], "+") +
        sweep(Shell.array, 1, S_states, "+") +
        sweep(Shell.array, 1, mb_term, "+")
      
      R.star <- aperm(apply(Shell.array, c(1, 3), "+", R_disc_x), c(2, 1, 3))
      R.star[, 2:(R_disc + 1), ][which(R.star[, 2:(R_disc + 1), ] > R.cstr[, 2:(R_disc + 1), ])] <- NaN
      
      targ_t <- targ_seas[t]
      if (targ_t <= 0) {
        Deficit.arr <- 0 * R.star
      } else {
        Deficit.arr <- (targ_t - R.star) / targ_t
        Deficit.arr[Deficit.arr < 0] <- 0
      }
      
      Release_cost <- (abs(Deficit.arr))^loss_exp[1]
      
      S.t_plus_1 <- R.cstr - R.star
      S.t_plus_1[S.t_plus_1 < 0] <- 0
      
      Spill_costs <- S.t_plus_1 - capacity
      Spill_costs[Spill_costs < 0] <- 0
      Spill_costs <- (Spill_costs / quantile(Q, probs = spill_targ))^loss_exp[2]
      
      S.t_plus_1[S.t_plus_1 > capacity] <- capacity
      Vol_costs <- abs((S.t_plus_1 - vol_targ * capacity) / (vol_targ * capacity))^loss_exp[3]
      
      Stage_cost <- weights[1] * Release_cost +
        weights[2] * Spill_costs +
        weights[3] * Vol_costs
      
      Implied_S_state <- round(1 + (S.t_plus_1 / capacity) * (length(S_states) - 1))
      Implied_S_state[Implied_S_state < 1] <- 1
      Implied_S_state[Implied_S_state > length(S_states)] <- length(S_states)
      
      Cost_to_go.arr <- array(Cost_to_go[Implied_S_state],
                              dim = c(length(S_states), length(R_disc_x), length(Q.probs)))
      Min_cost_arr <- Stage_cost + Cost_to_go.arr
      Min_cost_arr_weighted <- sweep(Min_cost_arr, 3, Q.probs, "*")
      Min_cost_expected <- apply(Min_cost_arr_weighted, c(1, 2), sum)
      
      Bellman[, t] <- Cost_to_go
      Cost_to_go <- apply(Min_cost_expected, 1, min, na.rm = TRUE)
      R_policy[, t] <- apply(Min_cost_expected, 1, which.min)
    }
    
    sim_ratio <- sum(R_policy == R_policy_test) / (frq * length(S_states))
    message(sim_ratio)
    if (sim_ratio > tol) break
    R_policy_test <- R_policy
  }
  
  R_policy_frac <- (R_policy - 1) / R_disc
  total_penalty <- sum(Bellman[is.finite(Bellman)], na.rm = TRUE)
  
  if (plot) {
    plot(R_policy_frac[1, ], type = "l", ylim = c(0, 1),
         ylab = "Policy fraction", xlab = "Month")
    plot(Bellman[, 1], type = "l", ylab = "Bellman cost", xlab = "Storage state")
  }
  
  results <- list(
    release_policy = R_policy,
    release_policy_frac = R_policy_frac,
    Bellman = Bellman,
    storage_states = S_states,
    release_states = R_disc_x,
    target_seasonal = targ_seas,
    target_ts = target_ts,
    flow_disc = Q_disc,
    total_penalty = total_penalty
  )
  
  return(results)
}

##### 2.3. Daily simulation driven by the monthly SDP policy ######

# Daily forward simulation.
# The monthly policy is chosen once per month and then distributed evenly across
# the days of that month.
sdp_daily_sim <- function(grand_id,
                          res_tseries,
                          sim_period,
                          demand_table,
                          demand_adj = 1,
                          spill_targ = 0.95,
                          vol_targ = 0.75,
                          weights = c(0.7, 0.2, 0.1),
                          S_disc = 500,
                          R_disc = 100,
                          Q_disc = c(0, 0.2375, 0.475, 0.7125, 0.95, 1),
                          loss_exp = c(2, 2, 2),
                          S_initial = NA_real_,
                          plot = TRUE,
                          tol = 0.99,
                          rep_rrv = FALSE,
                          pcp_fix = 0.6) {
  
  full_series <- res_tseries %>%
    filter(GRAND_ID == grand_id, year(date) %in% sim_period) %>%
    arrange(date) %>%
    mutate(
      inflow_Mm3_day = inflow * 86400 / 1e6,
      evap_mm = replace_na(ET_mm, 0),
      pcp_mm = replace_na(PCP_mm, 0),
      error_balance_Mm3 = replace_na(daily_err_Mm3, 0),
      ym = format(date, "%Y-%m")
    )
  
  if (nrow(full_series) == 0) {
    stop(paste0("No time series found for GRAND_ID ", grand_id))
  }
  
  res_props <- reservoir_properties %>% filter(GRAND_ID == grand_id)
  if (nrow(res_props) != 1) stop(paste0("Invalid reservoir properties for GRAND_ID ", grand_id))
  
  cap_val <- res_props$Cap_mcm[1]
  hypso <- hypso_curves %>% filter(GRAND_ID == grand_id)
  if (nrow(hypso) != 1) stop(paste0("Invalid hypsometric curve for GRAND_ID ", grand_id))
  
  if (is.na(S_initial)) {
    initial_stor <- get_initial_storage(full_series, sim_period)
    if (length(initial_stor) == 0 || is.na(initial_stor)) {
      initial_stor <- full_series$storage[which(!is.na(full_series$storage))[1]]
    }
  } else {
    initial_stor <- S_initial * cap_val
  }
  
  # Monthly inputs for the solver
  monthly_inputs <- build_monthly_inputs(full_series)
  monthly_Q <- ts(
    monthly_inputs$inflow_month_Mm3,
    start = c(monthly_inputs$year[1], monthly_inputs$month[1]),
    frequency = 12
  )
  monthly_evap <- ts(
    monthly_inputs$evap_month_mm,
    start = c(monthly_inputs$year[1], monthly_inputs$month[1]),
    frequency = 12
  )
  monthly_pcp <- ts(
    monthly_inputs$pcp_month_mm,
    start = c(monthly_inputs$year[1], monthly_inputs$month[1]),
    frequency = 12
  )
  monthly_err <- ts(
    monthly_inputs$error_month_Mm3,
    start = c(monthly_inputs$year[1], monthly_inputs$month[1]),
    frequency = 12
  )
  
  # Monthly seasonal demand target for the reservoir
  target_seas <- build_monthly_demand_target(
    grand_id = grand_id,
    demand_table = demand_table,
    demand_adj = demand_adj
  )
  
  # Solve the monthly SDP policy
  policy_fit <- sdp_policy_solver_monthly(
    Q = monthly_Q,
    capacity = cap_val,
    target = target_seas,
    hypso = hypso,
    evap = monthly_evap,
    pcp = monthly_pcp,
    wb_err = monthly_err,
    R_max = 2 * max(target_seas, na.rm = TRUE),
    spill_targ = spill_targ,
    vol_targ = vol_targ,
    weights = weights,
    S_disc = S_disc,
    R_disc = R_disc,
    Q_disc = Q_disc,
    loss_exp = loss_exp,
    S_initial = initial_stor / cap_val,
    plot = FALSE,
    tol = tol,
    rep_rrv = rep_rrv,
    pcp_fix = pcp_fix
  )
  
  S_states <- policy_fit$storage_states
  R_disc_x <- policy_fit$release_states
  R_policy <- policy_fit$release_policy
  
  # Daily simulation arrays
  tsteps <- nrow(full_series)
  sim_storage <- numeric(tsteps)
  sim_outflow <- numeric(tsteps)
  sim_area <- numeric(tsteps)
  sim_pcp <- numeric(tsteps)
  sim_evap <- numeric(tsteps)
  sim_spill <- numeric(tsteps)
  sim_monthly_release <- numeric(tsteps)
  sim_daily_policy <- numeric(tsteps)
  
  sim_series <- full_series %>%
    select(date, inflow_Mm3_day, evap_mm, pcp_mm, error_balance_Mm3, ym) %>%
    mutate(
      sim_storage_Mm3 = NA_real_,
      sim_outflow_Mm3 = NA_real_,
      sim_area_km2 = NA_real_,
      sim_evap_Mm3 = NA_real_,
      sim_pcp_Mm3 = NA_real_,
      sim_spill_Mm3 = NA_real_,
      policy_release_Mm3_month = NA_real_,
      policy_release_Mm3_day = NA_real_,
      GRAND_ID = grand_id
    )
  
  monthly_ids <- unique(full_series$ym)
  idx <- 0
  current_storage <- initial_stor
  
  for (g in seq_along(monthly_ids)) {
    month_df <- full_series[full_series$ym == monthly_ids[g], ]
    month_i <- month(month_df$date[1])
    ndays <- nrow(month_df)
    
    S_state <- which.min(abs(S_states - current_storage))
    release_idx <- R_policy[S_state, month_i]
    monthly_release_target <- R_disc_x[release_idx]
    daily_release_target <- monthly_release_target / ndays
    
    for (j in seq_len(ndays)) {
      idx <- idx + 1
      
      storage_before_release <- if (idx == 1) {
        current_storage + month_df$inflow_Mm3_day[j]
      } else {
        sim_storage[idx - 1] + month_df$inflow_Mm3_day[j]
      }
      
      release_day <- min(storage_before_release, daily_release_target)
      storage_after_release <- max(0, storage_before_release - release_day)
      
      # Use the same mb_fluxes() function as the other policies
      mb_flux <- mb_fluxes(
        stor = storage_after_release,
        hypso = hypso,
        pcp = month_df$pcp_mm[j],
        evap = month_df$evap_mm[j],
        err = month_df$error_balance_Mm3[j],
        pcp_fix = pcp_fix,
        err_fix = TRUE
      )
      
      storage_pre_constr <- storage_after_release + mb_flux$mb_flux
      constrains <- physical_constrains(
        storage = storage_pre_constr,
        max_capacity = cap_val,
        prov_release = release_day
      )
      
      sim_area[idx] <- mb_flux$area_km
      sim_pcp[idx] <- mb_flux$sim_pcp_Mm3
      sim_evap[idx] <- mb_flux$sim_evap_Mm3
      sim_storage[idx] <- constrains$final_storage
      sim_outflow[idx] <- constrains$final_release
      sim_spill[idx] <- max(0, storage_pre_constr - cap_val)
      sim_monthly_release[idx] <- monthly_release_target
      sim_daily_policy[idx] <- daily_release_target
    }
    
    current_storage <- sim_storage[idx]
  }
  
  sim_series$sim_area_km2 <- sim_area
  sim_series$sim_pcp_Mm3 <- sim_pcp
  sim_series$sim_evap_Mm3 <- sim_evap
  sim_series$sim_storage_Mm3 <- sim_storage
  sim_series$sim_outflow_Mm3 <- sim_outflow
  sim_series$sim_spill_Mm3 <- sim_spill
  sim_series$policy_release_Mm3_month <- sim_monthly_release
  sim_series$policy_release_Mm3_day <- sim_daily_policy
  
  results <- list(
    release_policy = R_policy,
    Bellman = policy_fit$Bellman,
    storage = sim_storage,
    releases = sim_outflow,
    evap_loss = sim_evap,
    water_level = sim_area,
    spill = sim_spill,
    sim_series = sim_series,
    monthly_policy = policy_fit,
    target_monthly = target_seas
  )
  
  if (plot) {
    plot(sim_series$sim_outflow_Mm3, type = "l", ylab = "Controlled release (Mm3/d)", xlab = "Time")
    plot(sim_series$sim_storage_Mm3, type = "l", ylab = "Storage (Mm3)", xlab = "Time")
    abline(h = vol_targ * cap_val, lty = 2)
    plot(sim_series$sim_spill_Mm3, type = "l", ylab = "Spill (Mm3/d)", xlab = "Time")
  }
  
  return(results)
}

##### 2.4. Optimization wrapper ######

# Run DEoptim on the SDP calibration parameters.
# The SDP solver itself is the inner optimization.
# DEoptim controls and parameter bounds are passed in as inputs.
fitting_sdp <- function(grand_id,
                        sim_period,
                        demand_table,
                        parms,
                        low_bound,
                        upp_bound,
                        DEoptim_conf,
                        S_disc = 1000,
                        R_disc = 10,
                        Q_disc = c(0, 0.2375, 0.475, 0.7125, 0.95, 1),
                        loss_exp = c(2, 2, 2),
                        pcp_fix = 0.6) {
  
  obs_data <- reservoirs_time_series %>%
    filter(GRAND_ID == grand_id,
           year(date) %in% sim_period)
  
  if (nrow(obs_data) == 0) {
    stop(paste0("No observations found for GRAND_ID ", grand_id))
  }
  
  obs_month <- obs_data %>%
    mutate(
      obs_outflow_Mm3 = outflow * 86400 / 1e6,
      ym = floor_date(date, "month")
    ) %>%
    group_by(ym) %>%
    summarise(
      obs_storage_Mm3 = mean(storage, na.rm = TRUE),
      obs_outflow_Mm3 = sum(obs_outflow_Mm3, na.rm = TRUE),
      .groups = "drop"
    )
  
  obj_fn_sdp <- function(params) {
    demand_adj <- params[1]
    spill_targ <- params[2]
    vol_targ <- params[3]
    w_raw <- params[4:6]
    weights <- w_raw / sum(w_raw)
    
    sim <- tryCatch(
      sdp_daily_sim(
        grand_id = grand_id,
        res_tseries = reservoirs_time_series,
        sim_period = sim_period,
        demand_table = demand_table,
        demand_adj = demand_adj,
        spill_targ = spill_targ,
        vol_targ = vol_targ,
        weights = weights,
        S_disc = S_disc,
        R_disc = R_disc,
        Q_disc = Q_disc,
        loss_exp = loss_exp,
        plot = FALSE,
        tol = 0.99,
        rep_rrv = FALSE,
        pcp_fix = pcp_fix
      ),
      error = function(e) NULL
    )
    
    if (is.null(sim)) return(1e6)
    
    tab <- tryCatch(
      obs_sim_tab_monthly(sim_data = sim$sim_series, obs_data = obs_data),
      error = function(e) NULL
    )
    if (is.null(tab)) return(1e6)
    
    ok_out <- is.finite(tab$sim_outflow_Mm3) & is.finite(tab$obs_outflow_Mm3)
    ok_sto <- is.finite(tab$sim_storage_Mm3) & is.finite(tab$obs_storage_Mm3)
    
    if (sum(ok_out) < 2 || sum(ok_sto) < 2) return(1e6)
    
    nrmse_out <- hydroGOF::nrmse(
      sim = tab$sim_outflow_Mm3[ok_out],
      obs = tab$obs_outflow_Mm3[ok_out],
      na.rm = TRUE,
      norm = "sd"
    )
    
    nrmse_sto <- hydroGOF::nrmse(
      sim = tab$sim_storage_Mm3[ok_sto],
      obs = tab$obs_storage_Mm3[ok_sto],
      na.rm = TRUE,
      norm = "sd"
    )
    
    0.5 * nrmse_out + 0.5 * nrmse_sto
  }
  
  fit <- DEoptim(
    fn = obj_fn_sdp,
    lower = low_bound,
    upper = upp_bound,
    control = DEoptim_conf
  )
  
  best_pars <- tibble(
    !!parms[1] := fit$optim$bestmem[1],
    !!parms[2] := fit$optim$bestmem[2],
    !!parms[3] := fit$optim$bestmem[3],
    !!parms[4] := fit$optim$bestmem[4],
    !!parms[5] := fit$optim$bestmem[5],
    !!parms[6] := fit$optim$bestmem[6]
  ) %>%
    mutate(
      weight_sum = .data[[parms[4]]] + .data[[parms[5]]] + .data[[parms[6]]],
      !!parms[4] := .data[[parms[4]]] / weight_sum,
      !!parms[5] := .data[[parms[5]]] / weight_sum,
      !!parms[6] := .data[[parms[6]]] / weight_sum
    ) %>%
    select(-weight_sum)
  
  perf_evol <- tibble(
    iteration = seq_along(fit$member$bestvalit),
    obj_f = fit$member$bestvalit
  )
  
  best_ind <- tibble(
    iteration = seq_len(nrow(fit$member$bestmemit)),
    !!parms[1] := fit$member$bestmemit[, 1],
    !!parms[2] := fit$member$bestmemit[, 2],
    !!parms[3] := fit$member$bestmemit[, 3],
    !!parms[4] := fit$member$bestmemit[, 4],
    !!parms[5] := fit$member$bestmemit[, 5],
    !!parms[6] := fit$member$bestmemit[, 6]
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
    all_pop = all_pop,
    fit = fit
  )
}

#### 3. Running simulations for one reservoir ####

##### 3.1. Default simulation #####

gid <- 41
obs_data <- reservoirs_time_series %>%
  filter(GRAND_ID == gid)

# Default calibration parameters
parms <- c("demand_adj", "spill_targ", "vol_targ", "w1", "w2", "w3")
def_parms <- tibble(
  demand_adj = 1.5,
  spill_targ = 0.90,
  vol_targ = 0.90,
  w1 = 0.50, # demand
  w2 = 0.20, # spill
  w3 = 0.30  # storage
)

# Run default simulation
def_sim <- sdp_daily_sim(
  grand_id = gid,
  res_tseries = reservoirs_time_series,
  sim_period = 1990:2014,
  demand_table = wat_demand,
  demand_adj = def_parms$demand_adj,
  spill_targ = def_parms$spill_targ,
  vol_targ = def_parms$vol_targ,
  weights = c(def_parms$w1, def_parms$w2, def_parms$w3),
  S_disc = 500,
  R_disc = 100,
  Q_disc = c(0, 0.2375, 0.475, 0.7125, 0.95, 1),
  loss_exp = c(2, 2, 2),
  plot = FALSE,
  tol = 0.99,
  rep_rrv = FALSE,
  pcp_fix = 0.6
)

# Daily obs-sim table
obssim_tab <- obs_sim_tab(
  sim_data = def_sim$sim_series,
  obs_data = obs_data
)

# Performance
obssim_tab %>%
  summarise(
    KGE_stor = KGE(sim_storage_Mm3, obs_storage),
    KGE_out = KGE(sim_outflow_Mm3, obs_outflow_Mm3),
    PBIAS_stor = pbias(sim_storage_Mm3, obs_storage),
    PBIAS_out = pbias(sim_outflow_Mm3, obs_outflow_Mm3),
    R2_stor = r2_cor(sim_storage_Mm3, obs_storage),
    R2_out = r2_cor(sim_outflow_Mm3, obs_outflow_Mm3)
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
  scale_color_manual(values = c("Observed" = "black", "Simulated" = "steelblue")) +
  theme_bw() +
  theme(axis.title = element_blank(), legend.title = element_blank())

##### 3.2. Parameters optimization ######

# Optimization settings are defined outside fitting_sdp()
parms <- c("demand_adj", "spill_targ", "vol_targ", "w1", "w2", "w3")

low_bound <- c(0.50, 0.50, 0.50, 0.3, 0.01, 0.01)
upp_bound <- c(2.00, 1.00, 1.00, 1.00, 0.6, 0.60)

n_pars <- length(parms)
max_itrs <- 100

DEoptim_conf <- DEoptim.control(
  NP = 10 * n_pars,
  itermax = max_itrs,
  F = 0.8,
  CR = 0.9,
  trace = F,
  reltol = 1e-4,
  steptol = 10,
  storepopfrom = 1,
  storepopfreq = 1
)

# Run optimization
stp_fit <- fitting_sdp(
  grand_id = gid,
  sim_period = 1990:2004,
  demand_table = wat_demand,
  parms = parms,
  low_bound = low_bound,
  upp_bound = upp_bound,
  DEoptim_conf = DEoptim_conf,
  S_disc = 500,
  R_disc = 100,
  Q_disc = c(0, 0.2375, 0.475, 0.7125, 0.95, 1),
  loss_exp = c(2, 2, 2),
  pcp_fix = 0.6
)

##### 3.3. Optimized simulation ######



opt_pars <- stp_fit$best_pars

opt_sim <- sdp_daily_sim(
  grand_id = gid,
  res_tseries = reservoirs_time_series,
  sim_period = 1990:2014,
  demand_table = wat_demand,
  demand_adj = opt_pars$demand_adj,
  spill_targ = opt_pars$spill_targ,
  vol_targ = opt_pars$vol_targ,
  weights = c(opt_pars$w1, opt_pars$w2, opt_pars$w3),
  S_disc = 500,
  R_disc = 100,
  Q_disc = c(0, 0.2375, 0.475, 0.7125, 0.95, 1),
  loss_exp = c(2, 2, 2),
  plot = FALSE,
  tol = 0.99,
  rep_rrv = FALSE,
  pcp_fix = 0.6
)

# Daily obs-sim table
obssim_tab_opt <- obs_sim_tab(
  sim_data = opt_sim$sim_series,
  obs_data = obs_data
)

# Performance
obssim_tab_opt %>%
  summarise(
    KGE_stor = KGE(sim_storage_Mm3, obs_storage),
    KGE_out = KGE(sim_outflow_Mm3, obs_outflow_Mm3),
    R2_stor = r2_cor(sim_storage_Mm3, obs_storage),
    R2_out = r2_cor(sim_outflow_Mm3, obs_outflow_Mm3)
  )

# Plot
obssim_tab_opt %>%
  pivot_longer(-date) %>%
  mutate(
    var = if_else(str_detect(name, "out"), "Release (Mm3/d)", "Storage (Mm3)"),
    type = if_else(str_detect(name, "sim"), "Simulated", "Observed")
  ) %>%
  ggplot(aes(x = date, y = value, color = type)) +
  geom_line() +
  facet_wrap(~var, scales = "free", ncol = 1) +
  scale_color_manual(values = c("Observed" = "black", "Simulated" = "steelblue")) +
  theme_bw() +
  theme(axis.title = element_blank(), legend.title = element_blank())

# Compare fitted vs default
rbind(
  (opt_sim$sim_series %>%
     select(date, sim_storage_Mm3, sim_outflow_Mm3) %>%
     mutate(pars = "fitted")),
  (def_sim$sim_series %>%
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
  scale_color_manual(values = c("Observed" = "black", "Sim. default" = "red4", "Sim. fitted" = "steelblue")) +
  theme_bw() +
  theme(axis.title = element_blank(), legend.title = element_blank())

