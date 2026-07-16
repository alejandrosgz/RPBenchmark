# Hanasaki reservoir policy

# Based on the Original H06 algorithm (Hanasaki et al., 2006)

# Main differences in water demand estimation and in water balance simulator


library(tidyverse)
library(hydroGOF)
library(DEoptim)
library(lubridate)
library(parallel)
library(patchwork)


#### 1. Input datasets ####

# Reservoirs to simulate
studied_reservoirs <- read_table('data/studied_reservoirs.txt')$studied_reservoirs

# Properties (from GDW and adapted considering clean time series)
reservoir_properties <- read_csv('data/Processed_data/reserv_properties.csv')

# Time series for simulations (includes water balance error derived from observations)
reservoirs_time_series <- read_csv('data/Processed_data/time_series_for_simulations.csv')

# Hypsometric curves (used to compute area from simulated storage)
hypso_curves <- read_csv('data/Processed_data/hypso_curves.csv')

# Water demand by reservoir and month (required by Hanasaki policy)
# Expected columns: GRAND_ID, Month, Total_demand_res_Mm3_month
wat_demand <- read_csv('data/Processed_data/reservoirs_demands_month.csv')


#### 2. Functions ####

##### 2.1. Common functions #####
###### 2.1.1. Functions used in simulations ######

# Use of hypsometric curve to calculate area from simulated storage
get_area <- function(hypso, storage){
  # Piecewise polynomial (linear/quadratic) depending on curve definition
  area <- if (hypso$poly2 == 0) {
    hypso$intercept + hypso$poly1 * storage
  } else if (!is.na(hypso$max_stor_rel) && storage > hypso$max_stor_rel) {
    # Cap area when exceeding valid storage range
    hypso$intercept + hypso$poly1 * hypso$max_stor_rel +
      hypso$poly2 * hypso$max_stor_rel^2
  } else {
    hypso$intercept + hypso$poly1 * storage +
      hypso$poly2 * storage^2
  }
  max(area, 0.0005) # Avoid zero area
}

# Calculating area, mass balance fluxes (P, E, and unknown fluxes extracted from observed data)
mb_fluxes <- function(stor,
                      hypso,
                      pcp,
                      evap,
                      err,
                      pcp_fix,
                      err_fix){
  
  # Calculate area (km2) with updated storage
  sim_area <- get_area(hypso = hypso, storage = stor)
  
  # Evaporation and precipitation net volume
  sim_pcp_Mm3 <- pcp_fix * (sim_area * pcp) / 1e3 # mm*km2 to Mm3
  sim_evap_Mm3 <- (sim_area * evap) / 1e3         # mm*km2 to Mm3
  
  et_pcp_flux <- sim_pcp_Mm3 - sim_evap_Mm3
  
  wb_err <- if (isTRUE(err_fix)) {
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
  # Calculate spill if overflow
  spill <- max(0, storage - max_capacity)
  
  # Add spill to current release
  release <- prov_release + spill
  
  # Subtract spill to get final storage
  final_storage <- max(0, storage - spill)
  
  tibble::new_tibble(
    list(
      final_release = release,
      final_storage = final_storage
    ),
    nrow = 1L
  )
}

###### 2.1.2. Other functions (helpers) ######

# Observed/simulated table for releases in m3/s
obs_sim_hanasaki <- function(sim_data,
                             obs_data){
  
  sim <- sim_data %>%
    select(date, sim_storage_Mm3, sim_outflow_m3s)
  
  obs <- obs_data %>%
    select(date, outflow, storage) %>%
    mutate(obs_outflow_m3s = outflow,
           obs_storage = storage) %>%
    select(date, obs_outflow_m3s, obs_storage)
  
  left_join(sim, obs, by = 'date')
}

# Observed/simulated table for release in Mm3/day
obs_sim_tab <- function(sim_data,
                        obs_data){
  
  sim <- sim_data %>%
    select(date, sim_storage_Mm3, sim_outflow_Mm3_day)
  
  obs <- obs_data %>%
    select(date, outflow, storage) %>%
    mutate(obs_outflow_Mm3_day = outflow * 86400 / 1e6,
           obs_storage = storage) %>%
    select(date, obs_outflow_Mm3_day, obs_storage)
  
  left_join(sim, obs, by = 'date')
}

r2_cor <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred)
  if (sum(ok) < 2) return(NA_real_)
  cor(obs[ok], pred[ok])^2
}

##### 2.2. Specific functions #####
###### 2.2.1. Hanasaki static preparation ######

# Static object: everything that does not depend on calibration parameters
hanasaki_static_prep <- function(tsries = tibble(),
                                 sim_period = vector(),
                                 structure_period = sim_period,
                                 wat_demand = tibble(),
                                 res_properties = tibble(),
                                 hypso_curves = tibble(),
                                 grand_id = numeric()) {
  
  # Input series for the selected reservoir
  ts_struct <- tsries %>%
    filter(GRAND_ID == grand_id, year(date) %in% structure_period) %>%
    arrange(date)
  
  if (nrow(ts_struct) == 0) return(NULL)
  
  # Annual mean inflow over the structure period
  inflow_av_annual <- ts_struct %>%
    summarise(av_ann_inflow_m3s = mean(inflow, na.rm = TRUE)) %>%
    mutate(av_ann_disch_Mm3 = av_ann_inflow_m3s * 86400 * 365 / 1e6)
  
  i_mean <- inflow_av_annual$av_ann_inflow_m3s
  if (!is.finite(i_mean) || i_mean <= 0) return(NULL)
  
  # Monthly inflow climatology
  inflow_av_month <- ts_struct %>%
    mutate(mon = month(date)) %>%
    group_by(mon) %>%
    summarise(av_mon_inflow_m3s = mean(inflow, na.rm = TRUE), .groups = "drop") %>%
    arrange(mon)
  
  # Month classification: Recharge / Release
  # Months 3--10 that are isolated are smoothed to match their predecessor
  month_types <- inflow_av_month %>%
    mutate(type_month = if_else(av_mon_inflow_m3s >= i_mean, "Recharge", "Release")) %>%
    mutate(type_month = if_else(
      mon %in% 3:10 &
        lag(type_month)    != type_month &
        lag(type_month, 2) != type_month &
        lead(type_month)   != type_month &
        lead(type_month, 2) != type_month,
      lag(type_month),
      type_month
    )) %>%
    select(mon, type_month)
  
  # Monthly time series (year prior to sim included for init storage)
  ts_month <- tsries %>%
    filter(GRAND_ID == grand_id,
           year(date) %in% (min(sim_period) - 1L):max(sim_period)) %>%
    mutate(monyear = floor_date(date, "month"),
           month_key = as.character(floor_date(date, "month"))) %>%
    group_by(monyear, month_key) %>%
    summarise(
      inflow_mon_m3s  = mean(inflow, na.rm = TRUE),
      storage_mon_Mm3 = mean(storage, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(monyear)
  
  if (nrow(ts_month) == 0) return(NULL)
  
  # Operational years
  # Count how many distinct release periods exist in the annual cycle
  release_periods_count <- month_types %>%
    arrange(mon) %>%
    mutate(prev = lag(type_month, default = last(type_month))) %>%
    summarise(n = sum(type_month == "Release" & prev == "Recharge")) %>%
    pull(n)
  
  if (release_periods_count <= 1) {
    # Single release period: op_yr starts at first release month each year
    op_yr_start <- ts_month %>%
      mutate(mon = month(monyear), yr = year(monyear)) %>%
      left_join(month_types, by = "mon") %>%
      filter(type_month == "Release", yr %in% sim_period) %>%
      arrange(monyear) %>%
      group_by(yr) %>%
      slice_min(monyear, with_ties = FALSE) %>%
      ungroup() %>%
      transmute(monyear, op_yr = row_number())
  } else {
    # Multiple release periods: start op_yr after the longest recharge run
    op_start_mon <- month_types %>%
      arrange(mon) %>%
      { bind_rows(., mutate(., mon = mon + 12)) } %>%
      mutate(grp = cumsum(type_month != lag(type_month, default = first(type_month))) + 1L) %>%
      group_by(grp) %>%
      summarise(type = first(type_month), n = n(), end = max(mon), .groups = "drop") %>%
      filter(type == "Recharge", lead(type) == "Release") %>%
      mutate(rel_start_mon = (end %% 12) + 1) %>%
      arrange(desc(n), rel_start_mon) %>%
      slice(1) %>%
      pull(rel_start_mon)
    
    op_yr_start <- ts_month %>%
      mutate(mon = month(monyear), yr = year(monyear)) %>%
      filter(mon == op_start_mon, yr %in% sim_period) %>%
      arrange(monyear) %>%
      transmute(monyear, op_yr = row_number())
  }
  
  if (nrow(op_yr_start) == 0) return(NULL)
  
  # Fill the operational year down through the monthly series
  ts_month <- ts_month %>%
    left_join(op_yr_start, by = "monyear") %>%
    fill(op_yr, .direction = "down") %>%
    mutate(op_yr = replace_na(op_yr, 0L),
           mon   = month(monyear))
  
  # Initial storage
  first_sim_day <- as.Date(op_yr_start$monyear[1])
  
  stor_init <- tsries %>%
    filter(
      GRAND_ID == grand_id,
      !is.na(storage),
      date < first_sim_day
    ) %>%
    arrange(desc(date)) %>%
    slice_head(n = 1) %>%
    pull(storage)
  
  # If not available, take the closest date
  if (length(stor_init) == 0 || is.na(stor_init[1])) {
    stor_init <- tsries %>%
      filter(GRAND_ID == grand_id, !is.na(storage)) %>%
      mutate(dist = abs(as.numeric(date - first_sim_day))) %>%
      arrange(dist, date) %>%
      slice_head(n = 1) %>%
      pull(storage)
  }
  
  if (length(stor_init) == 0 || is.na(stor_init[1])) {
    stop("No valid initial storage found for GRAND_ID ", grand_id)
  }
  
  stor_init <- stor_init[1]
  
  # Reservoir properties
  res_props   <- res_properties %>% filter(GRAND_ID == grand_id)
  hypso_curve <- hypso_curves %>% filter(GRAND_ID == grand_id)
  
  if (nrow(res_props) != 1 || nrow(hypso_curve) != 1) return(NULL)
  
  # Capacity ratio used by the Hanasaki regime
  c_val <- res_props$Cap_mcm / inflow_av_annual$av_ann_disch_Mm3
  c_val <- if (!is.finite(c_val)) 0 else c_val
  
  # Reservoir purpose
  purpose <- if (isTRUE(str_detect(res_props$Uses, "rriga"))) "irrigation" else "other"
  
  # Water demand (Mm3/month -> m3/s)
  res_demand <- wat_demand %>%
    filter(GRAND_ID == grand_id) %>%
    select(Month, Total_demand_res_Mm3_month)
  
  if (nrow(res_demand) == 0) {
    demand_m3s <- tibble(mon = 1:12, demand_m3s = 0)
  } else {
    demand_m3s <- res_demand %>%
      transmute(
        mon = Month,
        demand_m3s = Total_demand_res_Mm3_month * 1e6 /
          (86400 * days_in_month(as.Date(paste(2001, Month, 1, sep = "-"))))
      ) %>%
      arrange(mon)
  }
  
  d_mean_m3s <- mean(demand_m3s$demand_m3s, na.rm = TRUE)
  if (!is.finite(d_mean_m3s)) d_mean_m3s <- 0
  
  aa_demand_Mm3 <- sum(res_demand$Total_demand_res_Mm3_month, na.rm = TRUE)
  if (!is.finite(aa_demand_Mm3)) aa_demand_Mm3 <- 0
  
  # Precompute the non-parameter-dependent daily data
  ts_day <- tsries %>%
    filter(GRAND_ID == grand_id,
           floor_date(date, "month") %in% ts_month$monyear) %>%
    arrange(date) %>%
    mutate(
      inflow        = replace_na(inflow, 0),
      ET_mm         = replace_na(ET_mm, 0),
      PCP_mm        = replace_na(PCP_mm, 0),
      daily_err_Mm3 = replace_na(daily_err_Mm3, 0),
      monyear       = floor_date(date, "month"),
      month_key     = as.character(floor_date(date, "month"))
    )
  
  if (nrow(ts_day) == 0) return(NULL)
  
  idx_list <- split(seq_len(nrow(ts_day)), ts_day$month_key)
  
  list(
    res_props         = res_props,
    hypso_curve       = hypso_curve,
    month_tseries_base = ts_month,
    init_stor         = stor_init,
    aa_inflow         = inflow_av_annual,
    i_mean            = i_mean,
    c_val             = c_val,
    purpose           = purpose,
    d_mean_m3s        = d_mean_m3s,
    aa_demand_Mm3     = aa_demand_Mm3,
    demand_m3s        = demand_m3s,
    month_types       = month_types,
    ts_day            = ts_day,
    idx_list          = idx_list
  )
}

###### 2.2.2. Hanasaki parameter preparation ######

# Parameter-dependent object: everything that changes with alpha, c_thr, dmd_thr and min_release
hanasaki_param_prep <- function(static,
                                alpha = 0.85,
                                c_thr = 0.5,
                                dmd_thr = 0.5,
                                min_release = 0.5,
                                dmd_factor = 1) {
  
  # Make sure the capacity threshold is valid
  c_thr <- if (is.finite(c_thr) && c_thr > 0) c_thr else NA_real_
  
  # Small-reservoir blending weight
  w <- if (is.finite(c_thr)) min((static$c_val / c_thr)^2, 1) else 0
  if (!is.finite(w)) w <- 0
  
  # Demand-to-inflow ratio used by the irrigation regime
  d_mean_m3s_adj <- static$d_mean_m3s * dmd_factor
  
  dmd_ratio <- if (is.finite(static$i_mean) && static$i_mean > 0) {
    d_mean_m3s_adj / static$i_mean
  } else {
    0
  }
  if (!is.finite(dmd_ratio)) dmd_ratio <- 0
  
  #dmd_ratio <- if (is.finite(static$i_mean) && static$i_mean > 0) static$d_mean_m3s / static$i_mean else 0
  #if (!is.finite(dmd_ratio)) dmd_ratio <- 0
  
  # Provisional release equation class
  prov_release_eq <- dplyr::case_when(
    static$purpose != "irrigation" | d_mean_m3s_adj <= 0 ~ "non_irr",
    dmd_ratio < dmd_thr                                      ~ "irr_low_demand",
    TRUE                                                     ~ "irr_high_demand"
  )
  
  # Large or small reservoir logic
  final_release_eq <- if (isTRUE(static$c_val >= c_thr)) "large_reservoir" else "small_reservoir"
  
  # Monthly series used by the simulator
  ts_month <- static$month_tseries_base %>%
    left_join(static$demand_m3s, by = "mon") %>%
    mutate(demand_m3s = replace_na(demand_m3s, 0)*dmd_factor)
  
  # Build the provisional monthly release series r0
  if (prov_release_eq == "irr_low_demand") {
    prov_release_vec <- static$i_mean + (ts_month$demand_m3s - d_mean_m3s_adj)
  } else if (prov_release_eq == "irr_high_demand") {
    if (is.finite(d_mean_m3s_adj) && d_mean_m3s_adj > 0) {
      prov_release_vec <- (static$i_mean * min_release) * (1 + ts_month$demand_m3s / d_mean_m3s_adj)
    } else {
      prov_release_vec <- rep(static$i_mean, nrow(ts_month))
    }
  } else {
    prov_release_vec <- rep(static$i_mean, nrow(ts_month))
  }
  
  # Avoid negative provisional releases
  prov_release_vec <- pmax(0, prov_release_vec)
  
  # Attach provisional release to the monthly table
  ts_month <- ts_month %>%
    mutate(prov_release_m3s = prov_release_vec)
  
  list(
    alpha            = alpha,
    c_thr            = c_thr,
    dmd_thr          = dmd_thr,
    min_release      = min_release,
    d_mean_m3s       = d_mean_m3s_adj,
    w                = w,
    dmd_ratio        = dmd_ratio,
    prov_release_eq  = prov_release_eq,
    final_release_eq = final_release_eq,
    month_tseries    = ts_month
  )
}

###### 2.2.3. Hanasaki simulation ######

# Simulation: daily water balance using the static object and the parameter object
hanasaki_sim_fast <- function(static,
                              param_prep,
                              pcp_fix = 0.6,
                              err_fix = TRUE) {
  
  # Preparing input data
  precip_fix <- pcp_fix
  mberr_fix  <- err_fix
  
  # Monthly table with demand and provisional release
  ts_month <- param_prep$month_tseries
  
  # Keep only the months that belong to the simulated operational years
  ts_mon <- ts_month %>% filter(op_yr > 0)
  if (nrow(ts_mon) == 0) {
    stop("No simulation months found (op_yr > 0) for GRAND_ID ", static$res_props$GRAND_ID[1])
  }
  
  # Daily data restricted to the simulated months only
  valid_months <- ts_mon$monyear
  ts_day <- static$ts_day %>%
    filter(monyear %in% valid_months) %>%
    arrange(date)
  
  if (nrow(ts_day) == 0) {
    stop("No daily data found for GRAND_ID ", static$res_props$GRAND_ID[1])
  }
  
  # Initial storage and release scaling factor
  stor     <- static$init_stor
  krls     <- stor / (param_prep$alpha * static$res_props$Cap_mcm)
  cur_opyr <- NA_integer_
  
  # Preallocate output list
  results <- vector("list", nrow(ts_day))
  k <- 0L
  
  # Run month by month, then day by day inside each month
  for (i in seq_len(nrow(ts_mon))) {
    
    mon_row <- ts_mon[i, ]
    
    # Recompute krls at the start of each operational year
    opy <- mon_row$op_yr
    if (!identical(opy, cur_opyr)) {
      cur_opyr <- opy
      krls     <- stor / (param_prep$alpha * static$res_props$Cap_mcm)
    }
    
    # Provisional release depends on reservoir class
    release_m3s <- if (param_prep$final_release_eq == "large_reservoir") {
      krls * mon_row$prov_release_m3s
    } else {
      param_prep$w * krls * mon_row$prov_release_m3s +
        (1 - param_prep$w) * mon_row$inflow_mon_m3s
    }
    
    # Convert monthly release from m3/s to Mm3/day
    release_Mm3_day <- release_m3s * 86400 / 1e6
    
    # Daily records inside the current month
    days <- ts_day[ts_day$monyear == mon_row$monyear, ]
    
    for (j in seq_len(nrow(days))) {
      
      d <- days[j, ]
      inflow_d <- d$inflow * 86400 / 1e6
      
      # Available water before release
      avail   <- stor + inflow_d
      rel_day <- min(release_Mm3_day, avail)
      stor    <- avail - rel_day
      
      # Compute evaporation, precipitation and mass-balance fluxes
      mb <- mb_fluxes(
        stor    = stor,
        hypso   = static$hypso_curve,
        pcp     = d$PCP_mm,
        evap    = d$ET_mm,
        err     = d$daily_err_Mm3,
        pcp_fix = precip_fix,
        err_fix = mberr_fix
      )
      
      # Update storage after mass-balance fluxes
      stor <- stor + mb$mb_flux
      
      # Apply physical constraints to avoid overflow and negative storage
      constrains <- physical_constrains(
        storage      = stor,
        max_capacity = static$res_props$Cap_mcm,
        prov_release = rel_day
      )
      
      stor    <- constrains$final_storage
      rel_day <- constrains$final_release
      
      k <- k + 1L
      results[[k]] <- tibble(
        date                = d$date,
        sim_storage_Mm3     = stor,
        sim_outflow_Mm3_day = rel_day,
        sim_outflow_m3s     = rel_day * 1e6 / 86400,
        sim_area_km2        = mb$area_km,
        sim_pcp_Mm3         = mb$sim_pcp_Mm3,
        sim_evap_Mm3        = mb$sim_evap_Mm3,
        GRAND_ID            = static$res_props$GRAND_ID[1]
      )
    }
  }
  
  list(
    sim = bind_rows(results[1:k]),
    metadata = tibble(
      GRAND_ID         = static$res_props$GRAND_ID[1],
      alpha            = param_prep$alpha,
      c_thr            = param_prep$c_thr,
      dmd_thr          = param_prep$dmd_thr,
      min_release      = param_prep$min_release,
      c_val            = static$c_val,
      w                = param_prep$w,
      purpose          = static$purpose,
      i_mean_m3s       = static$i_mean,
      d_mean_m3s       = param_prep$d_mean_m3s,
      aa_demand_Mm3    = static$aa_demand_Mm3,
      dmd_ratio        = param_prep$dmd_ratio,
      prov_release_eq  = param_prep$prov_release_eq,
      final_release_eq = param_prep$final_release_eq
    )
  )
}

###### 2.2.4. Functions used for optimization ######

# Function to run the optimization
fitting_hanasaki <- function(grand_id,
                             sim_period,
                             structure_period = sim_period,
                             reservoirs_time_series,
                             wat_demand,
                             reservoir_properties,
                             hypso_curves,
                             low_bound,
                             upp_bound,
                             DEoptim_conf,
                             parms){
  
  # Used in obj. function
  static <- hanasaki_static_prep(
    tsries           = reservoirs_time_series,
    sim_period       = sim_period,
    structure_period = structure_period,
    wat_demand       = wat_demand,
    res_properties   = reservoir_properties,
    hypso_curves     = hypso_curves,
    grand_id         = grand_id
  )
  
  if (is.null(static)) {
    return(list(
      best_pars = NULL,
      perf_evol = NULL,
      best_ind  = NULL,
      all_pop   = NULL,
      static    = NULL,
      fit       = NULL
    ))
  }
  
  static$obs_data <- reservoirs_time_series %>%
    filter(GRAND_ID == grand_id,
           year(date) %in% sim_period)
  
  # Function to use inside deoptim --> Run the simulation and calculate obj. function
  obj_fn_hanasaki <- function(params) {
    
    # Parameter names --> Order!
    alpha       <- params[1]
    c_thr       <- params[2]
    dmd_thr     <- params[3]
    min_release <- params[4]
    dmd_factor  <- params[5]
    
    # Run simulation (objects should be on environment)
    sim <- tryCatch(
      {
        param_prep <- hanasaki_param_prep(
          static       = static,
          alpha        = alpha,
          c_thr        = c_thr,
          dmd_thr      = dmd_thr,
          min_release  = min_release,
          dmd_factor   = dmd_factor
        )
        
        hanasaki_sim_fast(
          static      = static,
          param_prep  = param_prep
        )
      },
      error = function(e) NULL
    )
    
    # If error, return 1e6
    if (is.null(sim)) return(1e6)
    
    # Obs. sim. tab
    tab <- tryCatch(
      obs_sim_hanasaki(sim_data = sim$sim,
                       obs_data = static$obs_data),
      error = function(e) NULL
    )
    
    # If error, return 1e6
    if (is.null(tab)) return(1e6)
    
    # Valid rows
    ok_out <- is.finite(tab$sim_outflow_m3s) & is.finite(tab$obs_outflow_m3s)
    ok_sto <- is.finite(tab$sim_storage_Mm3) & is.finite(tab$obs_storage)
    
    if (sum(ok_out) < 2 || sum(ok_sto) < 2) return(1e6)
    
    # nRMSE for each variable
    nrmse_out <- hydroGOF::nrmse(
      sim = tab$sim_outflow_m3s[ok_out],
      obs = tab$obs_outflow_m3s[ok_out],
      na.rm = TRUE,
      norm = "sd"
    )
    nrmse_sto <- hydroGOF::nrmse(
      sim = tab$sim_storage_Mm3[ok_sto],
      obs = tab$obs_storage[ok_sto],
      na.rm = TRUE,
      norm = "sd"
    )
    
    # Objective funcion value
    0.5 * nrmse_out + 0.5 * nrmse_sto
    
  }
  
  fit <- DEoptim(fn = obj_fn_hanasaki,
                 lower = low_bound,
                 upper = upp_bound,
                 control = DEoptim_conf)
  
  best_pars <- tibble(!!parms[1] := fit$optim$bestmem[1],
                      !!parms[2] := fit$optim$bestmem[2],
                      !!parms[3] := fit$optim$bestmem[3],
                      !!parms[4] := fit$optim$bestmem[4],
                      !!parms[5] := fit$optim$bestmem[5])
  
  perf_evol <- tibble(iteration = seq_along(fit$member$bestvalit),
                      obj_f = fit$member$bestvalit)
  
  best_ind <- tibble(iteration = seq_len(nrow(fit$member$bestmemit)),
                     !!parms[1] := fit$member$bestmemit[,1],
                     !!parms[2] := fit$member$bestmemit[,2],
                     !!parms[3] := fit$member$bestmemit[,3],
                     !!parms[4] := fit$member$bestmemit[,4],
                     !!parms[5] := fit$member$bestmemit[,5])
  
  all_pop <- purrr::map_dfr(seq_along(fit$member$storepop), function(i) {
    tt <- tibble::as_tibble(fit$member$storepop[[i]])
    colnames(tt) <- parms
    dplyr::mutate(tt, iteration = i, .before = 1)
  })
  
  list(best_pars = best_pars,
       perf_evol = perf_evol,
       best_ind = best_ind,
       all_pop = all_pop,
       static = static,
       fit = fit)
}


#### 3. Running simulations for one reservoir #### 
##### 3.1. Default simulation #####

gid <- 132

# Observed data
obs_data <- reservoirs_time_series %>%
  filter(GRAND_ID == gid)

# Simulation parameters
parms <- c('alpha', 'c_thr', 'dmd_thr', 'min_release', 'dmd_factor')
def_parms <- tibble(alpha = 0.85,
                    c_thr = 0.5,
                    dmd_thr = 0.5,
                    min_release = 0.5,
                    dmd_factor = 1)

# Run simulation
def_sim <- tryCatch(
  hanasaki_sim_fast(
    static = {
      st <- hanasaki_static_prep(
        tsries           = reservoirs_time_series,
        sim_period       = 1990:2015,
        structure_period = 1990:2015,
        wat_demand       = wat_demand,
        res_properties   = reservoir_properties,
        hypso_curves     = hypso_curves,
        grand_id         = gid
      )
      st$obs_data <- obs_data
      st
    },
    param_prep = hanasaki_param_prep(
      static       = hanasaki_static_prep(
        tsries           = reservoirs_time_series,
        sim_period       = 1990:2015,
        structure_period = 1990:2015,
        wat_demand       = wat_demand,
        res_properties   = reservoir_properties,
        hypso_curves     = hypso_curves,
        grand_id         = gid
      ),
      alpha        = def_parms$alpha,
      c_thr        = def_parms$c_thr,
      dmd_thr      = def_parms$dmd_thr,
      min_release  = def_parms$min_release,
      dmd_factor   = def_parms$dmd_factor
    )
  ),
  error = function(e) NULL
)

# Obs.sim. table
obssim_tab <- obs_sim_hanasaki(sim_data = def_sim$sim,
                               obs_data = obs_data)

# Performance
obssim_tab %>%
  summarise(KGE_stor = KGE(sim_storage_Mm3, obs_storage),
            KGE_out = KGE(sim_outflow_m3s, obs_outflow_m3s))


# Plot
obssim_tab %>%
  pivot_longer(-date) %>%
  mutate(var = if_else(str_detect(name, 'out'), 'Release (m3/s)', 'Storage (Mm3)'),
         type = if_else(str_detect(name, 'sim'), 'Simulated', 'Observed')) %>%
  ggplot(aes(x = date, y = value, color = type)) +
  geom_line() +
  facet_wrap(~var, scales = 'free', ncol = 1) +
  scale_color_manual(values = c('black', 'steelblue')) +
  theme_bw() +
  theme(axis.title = element_blank(), legend.title = element_blank())


##### 3.2. Parameters optimization #####

gid <- 132

# Observed data
obs_data <- reservoirs_time_series %>%
  filter(GRAND_ID == gid)

# Simulation parameters
parms <- c('alpha', 'c_thr', 'dmd_thr', 'min_release', 'dmd_factor')
def_parms <- tibble(alpha = 0.85,
                    c_thr = 0.5,
                    dmd_thr = 0.5,
                    min_release = 0.5,
                    dmd_factor = 1)

n_pars <- length(parms)
max_itrs <- 200
low_bound <- c(0.40, 0.10, 0.10, 0.00, 0.2) # alpha, c_thr, dmd_thr, min_release, demand factor
upp_bound <- c(1.50, 2.00, 2.00, 0.60, 2)

# DE adjustments --> NP, F, CR default
DEoptim_conf <- DEoptim.control(
  NP = 10 * n_pars, # Number of individuals first pop
  itermax = max_itrs, # Max. number of sims. per generation
  F = 0.8, # Mutation extent
  CR = 0.9, # Mutation probability
  trace = FALSE,
  # To converge:
  reltol = 1e-4, # Minimum obj. funct. improvement
  steptol = 10, # In last n iterations
  # Keep obj. functions and population info
  storepopfrom = 1,
  storepopfreq = 1
)

# Run optimization
deoptim_fit <- fitting_hanasaki(
  grand_id = gid,
  sim_period = 1990:2004,
  structure_period = 1990:2004,
  reservoirs_time_series = reservoirs_time_series,
  wat_demand = wat_demand,
  reservoir_properties = reservoir_properties,
  hypso_curves = hypso_curves,
  low_bound = low_bound,
  upp_bound = upp_bound,
  DEoptim_conf = DEoptim_conf,
  parms = parms
)


##### 3.3. Optimized simulation #####

opt_pars <- deoptim_fit$best_pars

# Run simulation
opt_static <- hanasaki_static_prep(
  tsries           = reservoirs_time_series,
  sim_period       = 1990:2015,
  structure_period = 1990:2015,
  wat_demand       = wat_demand,
  res_properties   = reservoir_properties,
  hypso_curves     = hypso_curves,
  grand_id         = gid
)
opt_static$obs_data <- obs_data

opt_sim <- hanasaki_sim_fast(
  static = opt_static,
  param_prep = hanasaki_param_prep(
    static       = opt_static,
    alpha        = opt_pars$alpha,
    c_thr        = opt_pars$c_thr,
    dmd_thr      = opt_pars$dmd_thr,
    min_release  = opt_pars$min_release,
    dmd_factor =   opt_pars$dmd_factor
  )
)

# Obs.sim. table
obssim_tab <- obs_sim_hanasaki(sim_data = opt_sim$sim,
                               obs_data = obs_data)

# Performance
obssim_tab %>%
  summarise(KGE_stor = KGE(sim_storage_Mm3, obs_storage),
            KGE_out = KGE(sim_outflow_m3s, obs_outflow_m3s))


# Plot
obssim_tab %>%
  pivot_longer(-date) %>%
  mutate(var = if_else(str_detect(name, 'out'), 'Release (m3/s)', 'Storage (Mm3)'),
         type = if_else(str_detect(name, 'sim'), 'Simulated', 'Observed')) %>%
  ggplot(aes(x = date, y = value, color = type)) +
  geom_line() +
  facet_wrap(~var, scales = 'free', ncol = 1) +
  scale_color_manual(values = c('black', 'steelblue')) +
  theme_bw() +
  theme(axis.title = element_blank(), legend.title = element_blank())

rbind(
  (opt_sim$sim %>%
     select(date, sim_storage_Mm3, sim_outflow_m3s) %>%
     mutate(pars = 'fitted')),
  (def_sim$sim %>%
     select(date, sim_storage_Mm3, sim_outflow_m3s) %>%
     mutate(pars = 'default'))
) %>%
  left_join(., obs_data[, c('date', 'outflow', 'storage')]) %>%
  mutate(obs_storage = storage, obs_outflow_m3s = outflow) %>%
  select(date:pars, starts_with('obs')) %>%
  pivot_longer(., -c(date, pars)) %>%
  mutate(var = if_else(str_detect(name, 'out'), 'Release (m3/s)', 'Storage (Mm3)'),
         type = if_else(str_detect(name, 'sim'), 'Simulated', 'Observed')) %>%
  mutate(type = case_when(type == 'Observed' ~ 'Observed',
                          pars == 'default' ~ 'Sim. default',
                          .default = 'Sim. fitted')) %>%
  ggplot(aes(x = date, y = value, color = type)) +
  geom_line() +
  facet_wrap(~var, scales = 'free', ncol = 1) +
  scale_color_manual(values = c('black', 'steelblue', 'red4')) +
  theme_bw() +
  theme(axis.title = element_blank(), legend.title = element_blank())


