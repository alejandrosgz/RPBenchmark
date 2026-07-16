# ISTARF-CONUS (Inferred Storage Targets and Release Functions for Conterminous United States)
#

# Core functions have been extracted from the Demo script available in ISTARF_CONUS dataset 
#(Turner et al., 2021)

#
# Simulation uses the common hydrologic functions shared with the other methods.


library(tidyverse)
library(lubridate)
library(starfit)
library(patchwork)
library(hydroGOF)


#### 1. Input datasets ####

# Reservoirs to simulate
studied_reservoirs <- read_table("data/studied_reservoirs.txt")$studied_reservoirs

# Reservoir properties
reservoir_properties <- read_csv("data/Processed_data/reserv_properties.csv", show_col_types = FALSE)

# Clean time series used to recalculate mean inflow and release bounds
obs_inflow_clean_tseries <- read_csv("data/Processed_data/time_series_for_simulations.csv", show_col_types = FALSE) %>%
  group_by(GRAND_ID) %>%
  summarise(
    mean_inflow_m3s = mean(inflow, na.rm = TRUE),
    min_outflow_m3s  = quantile(outflow, 0.05, na.rm = TRUE),
    max_outflow_m3s  = quantile(outflow, 0.95, na.rm = TRUE),
    .groups = "drop"
  )

# STARFIT parameters fitted for all reservoirs
starfit_pars_default <- read_csv("data/release_policies/starfit/ISTARF-CONUS.csv", show_col_types = FALSE) %>%
  filter(GRanD_ID %in% studied_reservoirs) %>%
  filter(!duplicated(GRanD_ID))

# Time series for simulations (includes water balance error derived from observations)
reservoirs_time_series <- read_csv("data/Processed_data/time_series_for_simulations.csv", show_col_types = FALSE) %>%
  fill(inflow, storage, .direction = "down") %>%
  fill(inflow, storage, .direction = "up")

# Hypsometric curves (used in the evaporation simulation)
hypso_curves <- read_csv("data/Processed_data/hypso_curves.csv", show_col_types = FALSE)


#### 2. Functions ####

##### 2.1. Common functions #####

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

# Daily version of observed vs simulated tabulation
obs_sim_tab <- function(sim_data,
                        obs_data){
  
  sim <- sim_data %>%
    select(date, inflow_Mm3_day, sim_storage_Mm3, sim_outflow_Mm3_day)
  
  obs <- obs_data %>%
    select(date, outflow, storage) %>%
    mutate(
      obs_outflow_m3s = outflow,
      obs_storage = storage
    ) %>%
    select(date, obs_outflow_m3s, obs_storage)
  
  left_join(sim, obs, by = "date")
}

r2_cor <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred)
  if (sum(ok) < 2) return(NA_real_)
  cor(obs[ok], pred[ok])^2
}

##### 2.2. STARFIT fit functions #####

# Weekly aggregation helper
starfit_weekly_series <- function(tsries = tibble(),
                                  sim_period,
                                  grand_id){
  
  tsries %>%
    filter(GRAND_ID == grand_id,
           year(date) %in% sim_period) %>%
    arrange(date) %>%
    mutate(
      cal_year = year(date),
      week_start = floor_date(date, "week", week_start = 1)
    ) %>%
    group_by(cal_year) %>%
    mutate(wk = dense_rank(week_start)) %>%
    ungroup() %>%
    group_by(cal_year, wk) %>%
    mutate(days_in_week = n()) %>%
    ungroup() %>%
    group_by(cal_year) %>%
    mutate(
      first_week_days = first(days_in_week[wk == 1]),
      week_in_year = if_else(first_week_days < 4,
                             pmax(1L, wk - 1L),
                             wk),
      max_week = max(week_in_year),
      week_in_year = if_else(max_week > 52 & week_in_year == max_week, 52L, week_in_year)
    ) %>%
    ungroup() %>%
    group_by(cal_year, week_in_year) %>%
    summarise(
      inflow_Mm3_week  = sum(inflow  * 86400, na.rm = TRUE) / 1e6,
      outflow_Mm3_week = sum(outflow * 86400, na.rm = TRUE) / 1e6,
      pcp_mm           = sum(PCP_mm, na.rm = TRUE),
      evap_mm          = sum(ET_mm, na.rm = TRUE),
      storage_Mm3      = first(storage),
      error_Mm3_week   = sum(daily_err_Mm3, na.rm = TRUE),
      .groups = "drop"
    )
}

# STARFIT target fitting logic, adapted to local time series inputs
fit_targets <- function(grand_id,
                        fit_period,
                        tsries = reservoirs_time_series,
                        reservoir_attributes){
  
  res <- grand_id
  
  message(paste0("Fitting targets for dam ", res, ": ", reservoir_attributes$Name))
  
  storage_capacity_MCM <- reservoir_attributes$Cap_mcm
  
  storage_daily <- tsries %>%
    filter(GRAND_ID == res) %>%
    filter(year(date) %in% fit_period) %>%
    select(date, s_MCM = storage) %>%
    filter(!is.na(s_MCM))
  
  if(nrow(storage_daily) == 0){
    stop(paste0("No storage data available for GRAND_ID ", res))
  }
  
  start_date <- storage_daily$date %>% first()
  end_date <- storage_daily$date %>% last()
  
  storage_daily_clipped <- left_join(
    tibble(date = seq.Date(start_date, end_date, by = 1)),
    storage_daily,
    by = "date"
  )
  
  storage_weekly <- storage_daily_clipped %>%
    mutate(
      year = year(date),
      epiweek = epiweek(date)
    ) %>%
    filter(year >= 1989) %>%
    group_by(year, epiweek) %>%
    summarise(
      s_pct = round(100 * median(s_MCM, na.rm = TRUE) / storage_capacity_MCM, 2),
      .groups = "drop"
    ) %>%
    ungroup() %>%
    filter(epiweek %in% 1:52)
  
  capacity_violations <- storage_weekly %>% filter(s_pct > 100)
  minimum_violations <- storage_weekly %>% filter(s_pct < 0)
  
  if(nrow(capacity_violations) > 0){
    message(paste0(nrow(capacity_violations), " capacity violations found for dam ", res))
  }
  if(nrow(minimum_violations) > 0){
    message(paste0(nrow(minimum_violations), " minimum violations found for dam ", res))
  }
  
  storage_weekly_for_fitting <- storage_weekly %>%
    mutate(
      s_pct = if_else(s_pct > 100, 100, s_pct),
      s_pct = if_else(s_pct < 0, 0, s_pct)
    )
  
  # Days per week used for the NOR storage bounds (larger smaller days)
  n_points <- 3
  
  data_for_flood_harmonic <- storage_weekly_for_fitting %>%
    group_by(epiweek) %>%
    mutate(rank = rank(-s_pct, ties.method = "first", na.last = "keep")) %>%
    filter(rank <= n_points) %>%
    ungroup() %>%
    select(epiweek, s_pct) %>%
    arrange(epiweek)
  
  data_for_conservation_harmonic <- storage_weekly_for_fitting %>%
    group_by(epiweek) %>%
    mutate(rank = rank(s_pct, ties.method = "first", na.last = "keep")) %>%
    filter(rank <= n_points) %>%
    ungroup() %>%
    select(epiweek, s_pct) %>%
    arrange(epiweek)
  
  fit_constrained_harmonic <- get("fit_constrained_harmonic", envir = asNamespace("starfit"))
  
  p_flood_harmonic <- fit_constrained_harmonic(data_for_flood_harmonic) %>%
    .[["solution"]] %>%
    round(3)
  
  p_conservation_harmonic <- fit_constrained_harmonic(data_for_conservation_harmonic) %>%
    .[["solution"]] %>%
    round(3)
  
  targets_flood <- convert_parameters_to_targets(p_flood_harmonic, target_name = "NORhi")
  targets_cons  <- convert_parameters_to_targets(p_conservation_harmonic, target_name = "NORlo")
  
  max_flood_target <- max(targets_flood$NORhi)
  min_flood_target <- min(targets_flood$NORhi)
  max_cons_target <- max(targets_cons$NORlo)
  min_cons_target <- min(targets_cons$NORlo)
  
  if(p_flood_harmonic[4] > max_flood_target) p_flood_harmonic[4] <- Inf
  if(p_flood_harmonic[5] < min_flood_target) p_flood_harmonic[5] <- -Inf
  if(p_conservation_harmonic[4] > max_cons_target) p_conservation_harmonic[4] <- Inf
  if(p_conservation_harmonic[5] < min_cons_target) p_conservation_harmonic[5] <- -Inf
  
  list(
    id = res,
    weekly_storage = storage_weekly,
    `NSR upper bound` = p_flood_harmonic,
    `NSR lower bound` = p_conservation_harmonic
  )
}

# STARFIT original release fitting logic, adapted to local time series inputs
fit_release_function <- function(grand_id,
                                 fit_period,
                                 tsries = reservoirs_time_series,
                                 targets_path = NULL,
                                 reservoir_attributes){
  
  res <- grand_id
  storage_capacity_MCM <- reservoir_attributes$Cap_mcm
  
  if(is.null(targets_path)){
    targets_path <- fit_targets(
      grand_id = res,
      fit_period = fit_period,
      tsries = tsries,
      reservoir_attributes = reservoir_attributes
    )
  }
  
  storage_target_parameters <- tibble(
    pf = targets_path[["NSR upper bound"]],
    pm = targets_path[["NSR lower bound"]]
  )
  
  if(all(is.na(storage_target_parameters))){
    message("Storage targets unavailable due to lack of data!")
    return(list())
  }
  
  daily_ops <- tsries %>%
    filter(GRAND_ID == res) %>%
    filter(year(date) %in% fit_period) %>%
    arrange(date) %>%
    mutate(
      i = inflow * 86400 / 1e6 + daily_err_Mm3, # To consider error during fitting
      r = outflow * 86400 / 1e6
    ) %>%
    select(date, s = storage, i, r) %>%
    mutate(
      year = year(date),
      epiweek = epiweek(date)
    ) %>%
    filter(year >= 1989)
  
  daily_ops_non_spill_periods <- daily_ops %>%
    filter(s + i < storage_capacity_MCM)
  
  aggregate_to_epiweeks <- get("aggregate_to_epiweeks", envir = asNamespace("starfit"))
  back_calc_missing_flows <- get("back_calc_missing_flows", envir = asNamespace("starfit"))
  
  weekly_ops_NA_removed <- daily_ops %>%
    aggregate_to_epiweeks() %>%
    back_calc_missing_flows() %>%
    filter(!is.na(i) & !is.na(r), i >= 0, r >= 0)
  
  min_r_i_datapoints <- 260
  min_r_maxmin_days <- 365
  r_st_max_quantile <- 0.95
  r_st_min_quantile <- 0.05
  r_sq_tol <- 0.2
  weeks_per_year <- 365.25 / 7
  
  if(nrow(weekly_ops_NA_removed) <= min_r_i_datapoints){
    message("Insufficient data to build release function")
    return(list(
      id = res,
      `mean inflow from GRAND. (MCM / wk)` = reservoir_attributes[["i_MAF_MCM"]] / weeks_per_year,
      `mean inflow from obs. (MCM / wk)` = NA_real_,
      `release harmonic parameters` = rep(NA_real_, 4),
      `release residual model coefficients` = rep(NA_real_, 3),
      `release constraints` = c(NA_real_, NA_real_)
    ))
  }
  
  i_daily <- daily_ops %>% filter(!is.na(i)) %>% pull(i)
  if(length(i_daily) > min_r_i_datapoints * 7){
    i_mean <- mean(i_daily) * 7
  } else {
    i_mean <- mean(weekly_ops_NA_removed[["i"]])
  }
  
  training_data_unfiltered <- weekly_ops_NA_removed %>%
    left_join(convert_parameters_to_targets(storage_target_parameters[["pf"]], target_name = "NORhi"), by = "epiweek") %>%
    left_join(convert_parameters_to_targets(storage_target_parameters[["pm"]], target_name = "NORlo"), by = "epiweek") %>%
    mutate(
      avail_pct = 100 * (s_start / storage_capacity_MCM),
      availability_status = (avail_pct - NORlo) / (NORhi - NORlo),
      i_st = (i / i_mean) - 1,
      r_st = (r / i_mean) - 1
    )
  
  r_daily <- daily_ops_non_spill_periods %>%
    filter(!is.na(r)) %>%
    pull(r)
  
  if(length(r_daily) > min_r_maxmin_days){
    r_st_max <- ((quantile(r_daily, r_st_max_quantile, na.rm = TRUE) %>% unname() %>% round(4) * 7) / i_mean) - 1
    r_st_min <- ((quantile(r_daily, r_st_min_quantile, na.rm = TRUE) %>% unname() %>% round(4) * 7) / i_mean) - 1
  } else {
    r_st_vector <- training_data_unfiltered %>%
      filter(s_start + i < storage_capacity_MCM) %>%
      pull(r_st)
    r_st_min <- quantile(r_st_vector, r_st_min_quantile, na.rm = TRUE) %>% unname() %>% round(4)
    r_st_max <- quantile(r_st_vector, r_st_max_quantile, na.rm = TRUE) %>% unname() %>% round(4)
  }
  
  # Keep only weeks where storage lies within the normal operating range (NOR).
  # availability_status is standardized between:
  #   0 -> storage at NORlo
  #   1 -> storage at NORhi
  # Weeks outside this interval are excluded from fitting because they are
  # considered abnormal operating conditions (e.g., spill or severe depletion).
  training_data <- training_data_unfiltered %>%
    filter(availability_status <= 1, availability_status > 0)
  
  # Fit the seasonal harmonic release component following the STARFIT formulation.
  # The standardized release (r_st) is represented as a combination of:
  #   - annual sine/cosine harmonics
  #   - semiannual sine/cosine harmonics
  #
  # This captures the mean seasonal release pattern of the reservoir.
  #
  # The model is fit without an intercept because the release anomaly is already
  # standardized relative to mean inflow.
  st_r_harmonic <- lm(
    data = training_data,
    r_st ~ 0 + sin(2 * pi * epiweek / 52) + cos(2 * pi * epiweek / 52) +
      sin(4 * pi * epiweek / 52) + cos(4 * pi * epiweek / 52)
  ) %>%
    .$coefficients %>%
    unname() %>%
    round(4)
  
  # Reconstruct the fitted harmonic release signal and compute residuals.
  # r_st_resid represents the part of the standardized release that is not
  # explained by the seasonal harmonic behavior alone.
  # These residuals are later modeled as a function of:
  #   - storage availability within the NOR range
  #   - standardized inflow anomaly
  data_for_linear_model_of_release_residuals <- training_data %>%
    mutate(
      st_r_harmonic = st_r_harmonic[1] * sin(2 * pi * epiweek / 52) +
        st_r_harmonic[2] * cos(2 * pi * epiweek / 52) +
        st_r_harmonic[3] * sin(4 * pi * epiweek / 52) +
        st_r_harmonic[4] * cos(4 * pi * epiweek / 52),
      r_st_resid = r_st - st_r_harmonic
    )
  
  # Fit the residual release model.
  # The residual component of release is explained using:
  #   availability_status -> storage position between NORlo and NORhi
  #   i_st                -> standardized inflow anomaly
  st_r_residual_model <- lm(
    data = data_for_linear_model_of_release_residuals,
    r_st_resid ~ availability_status + i_st
  )
  
  # Extract and round regression coefficients:
  #   [1] intercept
  #   [2] storage sensitivity coefficient
  #   [3] inflow sensitivity coefficient
  st_r_residual_model_coef <- st_r_residual_model$coefficients %>%
    unname() %>%
    round(3)
  
  # If the storage coefficient becomes negative while the inflow coefficient
  # remains positive, refit the residual model using inflow only
  # Negative storage sensitivity is considered physically inconsistent because
  # releases are generally expected to increase with increasing storage
  if(st_r_residual_model_coef[2] < 0 & st_r_residual_model_coef[3] >= 0){
    st_r_residual_model <- lm(
      data = data_for_linear_model_of_release_residuals,
      r_st_resid ~ i_st
    )
    
    st_r_residual_model_coef <- c(
      st_r_residual_model$coefficients[[1]],
      0,
      st_r_residual_model$coefficients[[2]]
    ) %>% round(3)
  }
  
  # If the inflow coefficient becomes negative while the storage coefficient
  # remains positive, refit the residual model using storage only.
  #
  # Negative inflow sensitivity is considered physically inconsistent because
  # releases are generally expected to increase during high inflow periods.
  if(st_r_residual_model_coef[3] < 0 & st_r_residual_model_coef[2] >= 0){
    st_r_residual_model <- lm(
      data = data_for_linear_model_of_release_residuals,
      r_st_resid ~ availability_status
    )
    
    st_r_residual_model_coef <- c(
      st_r_residual_model$coefficients[[1]],
      st_r_residual_model$coefficients[[2]],
      0
    ) %>% round(3)
  }
  
  # Discard the residual model entirely if:
  #   - adjusted R² is below the minimum threshold
  #   - storage sensitivity remains negative
  #   - inflow sensitivity remains negative
  #
  # In this case, release behavior is represented only by the harmonic seasonal
  # component, following the original STARFIT logic.
  if(summary(st_r_residual_model)$adj.r.squared < r_sq_tol ||
     st_r_residual_model_coef[2] < 0 ||
     st_r_residual_model_coef[3] < 0){
    
    message("Release residual model will be discarded; release will be based on harmonic function only")
    
    st_r_residual_model_coef <- c(0, 0, 0)
  }
  
  list(
    id = res,
    `mean inflow from GRAND. (MCM / wk)` = reservoir_attributes[["i_MAF_MCM"]] / weeks_per_year,
    `mean inflow from obs. (MCM / wk)` = i_mean,
    `release harmonic parameters` = st_r_harmonic,
    `release residual model coefficients` = st_r_residual_model_coef,
    `release constraints` = c(r_st_min, r_st_max)
  )
}

# Builds a STARFIT parameter row from the fitted targets and release parameters
build_starfit_parameters <- function(grand_id,
                                     reservoir_attributes,
                                     targets_stor,
                                     release_pars,
                                     starfit_pars_default,
                                     obs_inflow_clean_tseries){
  
  default_row <- starfit_pars_default %>%
    filter(GRanD_ID == grand_id) %>%
    slice(1)
  
  default_row %>%
    transmute(
      GRanD_ID = reservoir_attributes$GRAND_ID,
      GRanD_NAME = reservoir_attributes$Name,
      GRanD_CAP_MCM = reservoir_attributes$Cap_mcm,
      GRanD_MEANFLOW_CUMECS,
      Obs_MEANFLOW_CUMECS = obs_inflow_clean_tseries %>%
        filter(GRAND_ID == grand_id) %>%
        pull(mean_inflow_m3s),
      fit = "full",
      match = NA,
      NORhi_mu    = targets_stor[["NSR upper bound"]][1],
      NORhi_alpha = targets_stor[["NSR upper bound"]][2],
      NORhi_beta   = targets_stor[["NSR upper bound"]][3],
      NORhi_max    = targets_stor[["NSR upper bound"]][4],
      NORhi_min    = targets_stor[["NSR upper bound"]][5],
      NORlo_mu    = targets_stor[["NSR lower bound"]][1],
      NORlo_alpha = targets_stor[["NSR lower bound"]][2],
      NORlo_beta   = targets_stor[["NSR lower bound"]][3],
      NORlo_max    = targets_stor[["NSR lower bound"]][4],
      NORlo_min    = targets_stor[["NSR lower bound"]][5],
      Release_alpha1 = release_pars[["release harmonic parameters"]][1],
      Release_beta1  = release_pars[["release harmonic parameters"]][2],
      Release_alpha2 = release_pars[["release harmonic parameters"]][3],
      Release_beta2  = release_pars[["release harmonic parameters"]][4],
      Release_c  = release_pars[["release residual model coefficients"]][1],
      Release_p1 = release_pars[["release residual model coefficients"]][2],
      Release_p2 = release_pars[["release residual model coefficients"]][3],
      Release_min = release_pars[["release constraints"]][1],
      Release_max = release_pars[["release constraints"]][2]
    )
}

##### 2.3. Daily simulation functions #####

# Returns the daily preprocessed time series and STARFIT parameters for a reservoir
starfit_prep_day <- function(tsries = tibble(),
                             sim_period,
                             grand_id,
                             starfit_pars,
                             reservoir_properties,
                             obs_inflow_clean_tseries) {
  
  res <- grand_id
  
  daily_time_series <- tsries %>%
    filter(GRAND_ID == res, year(date) %in% sim_period) %>%
    arrange(date) %>%
    mutate(
      cal_year = year(date),
      epiweek  = epiweek(date),
      epiweek  = if_else(epiweek > 52, 52L, epiweek),
      inflow_Mm3_day  = inflow * 86400 / 1e6,
      outflow_Mm3_day = outflow * 86400 / 1e6,
      pcp_mm          = replace_na(PCP_mm, 0),
      evap_mm         = replace_na(ET_mm, 0),
      storage_Mm3     = storage,
      error_Mm3_day   = replace_na(daily_err_Mm3, 0)
    ) %>%
    select(GRAND_ID, date, cal_year, epiweek,
           inflow_Mm3_day, outflow_Mm3_day, pcp_mm, evap_mm, storage_Mm3, error_Mm3_day)
  
  hypso_curve <- hypso_curves %>% filter(GRAND_ID == res)
  
  pars <- starfit_pars %>%
    filter(GRanD_ID == res) %>%
    mutate(
      Obs_MEANFLOW_CUMECS = if_else(
        is.na(Obs_MEANFLOW_CUMECS),
        GRanD_MEANFLOW_CUMECS,
        Obs_MEANFLOW_CUMECS
      )
    )
  
  if(nrow(pars) != 1){
    stop(paste0("STARFIT parameter row not found for GRAND_ID ", res))
  }
  
  NOR_upper_bound <- convert_parameters_to_targets(
    parameters = c(pars[["NORhi_mu"]],
                   pars[["NORhi_alpha"]],
                   pars[["NORhi_beta"]],
                   pars[["NORhi_max"]],
                   pars[["NORhi_min"]]),
    target_name = "NORhi"
  )
  
  NOR_lower_bound <- convert_parameters_to_targets(
    parameters = c(pars[["NORlo_mu"]],
                   pars[["NORlo_alpha"]],
                   pars[["NORlo_beta"]],
                   pars[["NORlo_max"]],
                   pars[["NORlo_min"]]),
    target_name = "NORlo"
  )
  
  NOR_bounds <- NOR_upper_bound %>%
    left_join(NOR_lower_bound, by = "epiweek")
  
  harmonic_release <- convert_parameters_to_release_harmonic(c(
    pars[["Release_alpha1"]],
    pars[["Release_beta1"]],
    pars[["Release_alpha2"]],
    pars[["Release_beta2"]]
  ))
  
  obs_mean_flow <- pars$Obs_MEANFLOW_CUMECS[1]
  if(is.na(obs_mean_flow)){
    obs_mean_flow <- mean(daily_time_series$inflow_Mm3_day, na.rm = TRUE) * 1e6 / 86400
  }
  
  inflow_Mm3_week <- 7 * obs_mean_flow * 86400 / 1e6
  
  # Inflow predictor uses the next 7 days, consistent with the daily simulator design
  n_days <- nrow(daily_time_series)
  daily_time_series$week_ahead_inflow_Mm3 <- vapply(
    1:n_days,
    function(i) {
      idx_end <- min(i + 6, n_days)
      sum(daily_time_series$inflow_Mm3_day[i:idx_end], na.rm = TRUE)
    },
    FUN.VALUE = numeric(1)
  )
  
  stand_inflows <- daily_time_series %>%
    select(date, cal_year, epiweek, week_ahead_inflow_Mm3) %>%
    mutate(
      mean_weekly_inflow = inflow_Mm3_week,
      standrzd_inflow = (week_ahead_inflow_Mm3 - mean_weekly_inflow) / mean_weekly_inflow
    ) %>%
    select(date, cal_year, epiweek, standrzd_inflow)
  
  list(
    daily_time_series = daily_time_series,
    hypso_curve = hypso_curve,
    pars = pars,
    NOR_bounds = NOR_bounds,
    harm_release = harmonic_release,
    stand_inflows = stand_inflows,
    inflow_Mm3_week = inflow_Mm3_week,
    reservoir_attributes = reservoir_properties %>% filter(GRAND_ID == res)
  )
}

# Returns the storage state relative to NOR bounds
stor_State <- function(storage,
                       pars,
                       epiweek,
                       input){
  
  NOR_vols <- input$NOR_bounds[input$NOR_bounds$epiweek == epiweek, ] %>%
    mutate(
      NORhi_vol = 0.01 * NORhi * pars$GRanD_CAP_MCM,
      NORlo_vol = 0.01 * NORlo * pars$GRanD_CAP_MCM
    )
  
  (storage - NOR_vols$NORlo_vol) / (NOR_vols$NORhi_vol - NOR_vols$NORlo_vol)
}

# Returns the daily target release using STARFIT harmonic and residual terms
target_Release <- function(pars,
                           epi_week,
                           harm_release,
                           stor_state,
                           stand_inflow,
                           av_inflow){
  
  release_adjustment <- pars$Release_c +
    pars$Release_p1 * stor_state +
    pars$Release_p2 * stand_inflow
  
  harm_release_adjusted <- harm_release[harm_release$epiweek == epi_week, ]$release_harmonic + release_adjustment
  
  (harm_release_adjusted * av_inflow) + av_inflow
}

# Runs the daily STARFIT simulation
starfit_sim_day <- function(input,
                            pcp_fix = 0.6,
                            err_fix = TRUE){
  
  pars <- input$pars
  
  # Maximum and minimum release standardized with mean inflow
  max_min_release_stand <- tibble(
    max_release_Mm3_day = (1 + as.numeric(pars$Release_max)) * input$inflow_Mm3_week / 7,
    min_release_Mm3_day = (1 + as.numeric(pars$Release_min)) * input$inflow_Mm3_week / 7
  )
  
  # Simulated time steps
  tsteps <- nrow(input$daily_time_series)
  
  # Tibble to fill with simulations
  sim_series <- tibble(
    date = input$daily_time_series$date,
    epiweek = input$daily_time_series$epiweek,
    inflow_Mm3_day = input$daily_time_series$inflow_Mm3_day,
    stand_inflow = input$stand_inflows$standrzd_inflow,
    evap_mm = input$daily_time_series$evap_mm,
    pcp_mm = input$daily_time_series$pcp_mm,
    error_balance_Mm3 = input$daily_time_series$error_Mm3_day
  ) %>%
    mutate(
      sim_storage_Mm3 = NA_real_,
      sim_outflow_Mm3_day = NA_real_,
      sim_area_km2 = NA_real_,
      sim_evap_Mm3_day = NA_real_,
      sim_pcp_Mm3_day = NA_real_,
      sim_time_step = NA_integer_,
      GRAND_ID = pars$GRanD_ID
    )
  
  sim_series$error_balance_Mm3[is.na(sim_series$error_balance_Mm3)] <- 0
  
  # Simulation for each time step
  for(tstp in 1:tsteps){
    
    tstep <- tstp
    epi_week <- sim_series[tstep, ]$epiweek
    
    init_stor <- if(tstep == 1){
      input$daily_time_series[1, ]$storage_Mm3 + sim_series[tstep, ]$inflow_Mm3_day
    } else {
      sim_series[tstep - 1, ]$sim_storage_Mm3 + sim_series[tstep, ]$inflow_Mm3_day
    }
    
    init_stor_state <- stor_State(
      storage = init_stor,
      pars = pars,
      epiweek = epi_week,
      input = input
    )
    
    stand_inflow <- sim_series[tstep, ]$stand_inflow
    
    target_release <- target_Release(
      pars = pars,
      epi_week = epi_week,
      harm_release = input$harm_release,
      stor_state = init_stor_state,
      stand_inflow = stand_inflow,
      av_inflow = input$inflow_Mm3_week
    ) / 7
    
    # Constrain the daily release to the STARFIT limits
    release_Mm3_day <- if(all(is.na(max_min_release_stand))){
      target_release
    } else if(target_release < max_min_release_stand$min_release_Mm3_day){
      max_min_release_stand$min_release_Mm3_day
    } else if(target_release > max_min_release_stand$max_release_Mm3_day){
      max_min_release_stand$max_release_Mm3_day
    } else {
      target_release
    }
    
    # Apply release before evaporation and precipitation
    storage_after_release <- init_stor - release_Mm3_day
    
    # Calculating area, mass balance fluxes (P, E, and unknown fluxes extracted from observed data)
    mb_flux <- mb_fluxes(
      stor = storage_after_release,
      hypso = input$hypso_curve,
      pcp = sim_series[tstep, ]$pcp_mm,
      evap = sim_series[tstep, ]$evap_mm,
      err = sim_series[tstep, ]$error_balance_Mm3,
      pcp_fix = pcp_fix,
      err_fix = err_fix
    )
    
    storage_pre_constr <- storage_after_release + mb_flux$mb_flux
    
    # Apply physical constrains to prevent overflow and negative storage
    constr <- physical_constrains(
      storage = storage_pre_constr,
      max_capacity = pars$GRanD_CAP_MCM,
      prov_release = release_Mm3_day
    )
    
    sim_series[tstep, ]$sim_storage_Mm3 <- constr$final_storage
    sim_series[tstep, ]$sim_outflow_Mm3_day <- constr$final_release
    sim_series[tstep, ]$sim_area_km2 <- mb_flux$area_km
    sim_series[tstep, ]$sim_evap_Mm3_day <- mb_flux$sim_evap_Mm3
    sim_series[tstep, ]$sim_pcp_Mm3_day <- mb_flux$sim_pcp_Mm3
    sim_series[tstep, ]$sim_time_step <- tstep
  }
  
  sim_series
}

# Compares observations and simulations at daily scale
starfit_perf <- function(sim_series,
                         obs_series,
                         grand_id,
                         sim_period){
  
  obs_tab <- obs_series %>%
    filter(GRAND_ID == grand_id,
           year(date) %in% sim_period) %>%
    arrange(date) %>%
    select(date, obs_storage = storage, obs_outflow_m3s = outflow)
  
  sim_tab <- sim_series %>%
    select(date, inflow_Mm3_day, sim_storage_Mm3, sim_outflow_Mm3_day) %>%
    mutate(sim_outflow_m3s = 1e6 * sim_outflow_Mm3_day / 86400)
  
  obs_sim_daily <- left_join(sim_tab, obs_tab, by = "date") %>%
    mutate(
      sim_storage = sim_storage_Mm3,
      sim_outflow = sim_outflow_m3s
    )
  
  perf <- obs_sim_daily %>%
    summarise(
      NSE_stor = hydroGOF::NSE(sim_storage, obs_storage),
      KGE_stor = hydroGOF::KGE(sim_storage, obs_storage),
      PBIAS_stor = hydroGOF::pbias(sim_storage, obs_storage),
      RMSE_stor = hydroGOF::rmse(sim_storage, obs_storage),
      R2_stor = r2_cor(sim_storage, obs_storage),
      NSE_out = hydroGOF::NSE(sim_outflow, obs_outflow_m3s),
      KGE_out = hydroGOF::KGE(sim_outflow, obs_outflow_m3s),
      PBIAS_out = hydroGOF::pbias(sim_outflow, obs_outflow_m3s),
      RMSE_out = hydroGOF::rmse(sim_outflow, obs_outflow_m3s),
      R2_out = r2_cor(sim_outflow, obs_outflow_m3s),
      nRMSE_stor = hydroGOF::nrmse(sim = sim_storage, obs = obs_storage, na.rm = TRUE, method = "sd"),
      nRMSE_out = hydroGOF::nrmse(sim = sim_outflow, obs = obs_outflow_m3s, na.rm = TRUE, method = "sd")
    )
  
  list(obs_sim_daily = obs_sim_daily, perf = perf)
}

# Plot observed and simulated storage and outflow
plot_sim_perf <- function(obssimtab,
                          perf,
                          props,
                          res_p){
  
  daily_plot <- obssimtab %>%
    pivot_longer(-c(date, inflow_Mm3_day), names_to = "name", values_to = "value") %>%
    mutate(type = str_sub(name, 1, 3),
           var = if_else(str_detect(name, "flow"), "Outflow", "Storage")) %>%
    ggplot(aes(x = date, y = value, color = type)) +
    geom_line(linewidth = 0.6) +
    facet_wrap(~var, scales = "free", ncol = 1) +
    theme_bw() +
    scale_color_manual(values = c("grey30", "red")) +
    labs(x = "Date", title = paste0(props$Name, " (", props$Cap_mcm, ") ", " Uses: ", props$Uses)) +
    theme(
      legend.position = c(0.9, 0.9),
      legend.background = element_blank(),
      legend.title = element_blank(),
      axis.title = element_blank(),
      text = element_text(size = 15),
      title = element_text(size = 9)
    )
  
  tab <- perf %>%
    mutate(across(everything(), ~ifelse(is.numeric(.x), round(.x, 2), .x))) %>%
    pivot_longer(everything()) %>%
    mutate(Metric = str_sub(name, 1, str_locate(name, "_")[, 1] - 1),
           Var = str_sub(name, str_locate(name, "_")[, 1] + 1, 100)) %>%
    mutate(Var = if_else(str_detect(name, "out"), "Outflow", "Storage")) %>%
    rename(Value = value) %>%
    select(Var, Metric, Value) %>%
    pivot_wider(names_from = Metric, values_from = Value)
  
  tabplot <- ggplot() +
    annotate(geom = "table", x = 1, y = 1, label = list(tab), size = 4.8) +
    theme_bw() +
    theme(
      axis.line = element_blank(),
      line = element_blank(),
      text = element_blank(),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.background = element_rect(fill = "transparent", color = NA)
    )
  
  complot <- daily_plot +
    inset_element(
      tabplot,
      left = 0.02,
      bottom = 0.75,
      right = 0.30,
      top = 0.98,
      align_to = "full"
    )
  
  complot
}


#### 3. Try in one reservoir ####

fit_period <- 1990:2004
eval_period <- 2005:2014
studied_period <- 1990:2014

res_p <- 132
props <- reservoir_properties %>% filter(GRAND_ID == res_p)

# Prepare STARFIT inputs for one reservoir using the default parameter set
starf_input <- starfit_prep_day(
  tsries = reservoirs_time_series,
  sim_period = studied_period,
  grand_id = res_p,
  starfit_pars = starfit_pars_default,
  reservoir_properties = reservoir_properties,
  obs_inflow_clean_tseries = obs_inflow_clean_tseries
)

# Run the daily simulation with the original STARFIT parameters
sim_series <- starfit_sim_day(input = starf_input, pcp_fix = 0.6, err_fix = TRUE)

# Compare against observations
perf_default_one <- starfit_perf(
  sim_series = sim_series,
  obs_series = reservoirs_time_series,
  grand_id = res_p,
  sim_period = studied_period
)

obsims <- perf_default_one$obs_sim_daily %>%
  select(date, sim_storage_Mm3, obs_storage, sim_outflow_m3s, obs_outflow_m3s) %>% 
  
  pivot_longer(-date) %>%
  mutate(
    type = str_sub(name, 1, 3),
    var = if_else(str_detect(name, "out"), "Outflow", "Storage"))

perf <- perf_default_one$obs_sim_daily %>%
  select(date, sim_storage_Mm3, obs_storage, sim_outflow_m3s, obs_outflow_m3s) %>%
  rename(
    sim_storage = sim_storage_Mm3,
    obs_storage = obs_storage,
    sim_outflow = sim_outflow_m3s,
    obs_outflow = obs_outflow_m3s
  ) %>%
  summarise(
    NSE_stor = hydroGOF::NSE(sim_storage, obs_storage),
    KGE_stor = hydroGOF::KGE(sim_storage, obs_storage),
    PBIAS_stor = hydroGOF::pbias(sim_storage, obs_storage),
    RMSE_stor = hydroGOF::rmse(sim_storage, obs_storage),
    R2_stor = r2_cor(sim_storage, obs_storage),
    NSE_out = hydroGOF::NSE(sim_outflow, obs_outflow),
    KGE_out = hydroGOF::KGE(sim_outflow, obs_outflow),
    PBIAS_out = hydroGOF::pbias(sim_outflow, obs_outflow),
    RMSE_out = hydroGOF::rmse(sim_outflow, obs_outflow),
    R2_out = r2_cor(sim_outflow, obs_outflow),
    nRMSE_stor = hydroGOF::nrmse(sim = sim_storage, obs = obs_storage, na.rm = TRUE, method = "sd"),
    nRMSE_out = hydroGOF::nrmse(sim = sim_outflow, obs = obs_outflow, na.rm = TRUE, method = "sd"),
    .groups = "drop"
  ) %>%
  mutate(GRAND_ID = res_p)

daily_plot <- ggplot(obsims, aes(x = date, y = value, color = type)) +
  geom_line(linewidth = 0.6) +
  geom_vline(xintercept = as.Date(paste0(max(fit_period), "-12-31")),
             linetype = "dashed",
             linewidth = 1,
             color = "blue4") +
  facet_wrap(~var, scales = "free", ncol = 1) +
  theme_bw() +
  scale_color_manual(values = c("grey30", "red")) +
  labs(
    x = "Date",
    title = paste0(props$Name, " (", props$Cap_mcm, ")  Uses: ", props$Uses)
  ) +
  theme(
    legend.position = c(0.9, 0.9),
    legend.background = element_blank(),
    legend.title = element_blank(),
    axis.title = element_blank(),
    text = element_text(size = 15),
    title = element_text(size = 9)
  )
daily_plot

tab <- perf %>%
  select(starts_with(c("KGE", "PBI", "R2", "NSE", "RMSE", "nRMSE")), GRAND_ID) %>%
  mutate(across(where(is.numeric), ~round(.x, 2))) %>%
  pivot_longer(-c(GRAND_ID)) %>%
  mutate(
    Metric = str_sub(name, 1, str_locate(name, "_")[, 1] - 1),
    Var = str_sub(name, str_locate(name, "_")[, 1] + 1, 100)
  ) %>%
  mutate(Var = if_else(str_detect(name, "out"), "Outflow", "Storage")) %>%
  rename(Value = value) %>%
  select(Var, Metric, Value) %>%
  pivot_wider(names_from = Metric, values_from = Value) %>%
  arrange(Var) %>%
  rename(`R²` = R2)

tab
