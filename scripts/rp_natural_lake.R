# Natural lake policy
library(tidyverse)
library(hydroGOF)
library(DEoptim)



# Based on the equation developed by Döll et al. (2003)
# Implemented in the common simulation framework developed for the benchmark

#### 1. Input datasets ####

# Reservoir to simulate
studied_reservoirs <- read_table('data/studied_reservoirs.txt')$studied_reservoirs

# Properties (from GDW and adapted considering clean time series)
reservoir_properties <- read_csv('data/Processed_data/reserv_properties.csv')

# Time series for simulations (includes water balance error derived from observations)
reservoirs_time_series <- read_csv('data/Processed_data/time_series_for_simulations.csv') #%>% 
 # fill(inflow, storage, .direction = "down") %>%
 # fill(inflow, storage, .direction = "up")

# Hypsometric curves (used in the evaporation simulation)
hypso_curves <- read_csv('data/Processed_data/hypso_curves.csv')


#### 2. Functions ####

##### 2.1. Common functions #####
###### 2.1.1. Functions used in simulations ######
# Use of hypsometric curve to calculate area from simulated storage
get_area <- function(hypso, storage){
  # Piecewise polynomial (linear/quadratic) depending on curve definition
  area <- if(hypso$poly2 == 0){
    hypso$intercept + hypso$poly1 * storage
  } else if(!is.na(hypso$max_stor_rel) && storage > hypso$max_stor_rel){
    # Cap area when exceeding valid storage range
    hypso$intercept + hypso$poly1 * hypso$max_stor_rel +
      hypso$poly2 * hypso$max_stor_rel^2
  } else{
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
  sim_pcp_Mm3 <- pcp_fix*(sim_area * pcp)/ 1e3 # mm*km2 to Mm3
  sim_evap_Mm3 <- (sim_area * evap) / 1e3 # mm*km2 to Mm3
  
  et_pcp_flux <- sim_pcp_Mm3 - sim_evap_Mm3
  
  wb_err <- if(isTRUE(err_fix)){
    err
  } else {
    0
  }
  
  mb_flux <- et_pcp_flux + wb_err
  
  tibble::new_tibble(
    list(area_km = sim_area,
         sim_pcp_Mm3 = sim_pcp_Mm3,
         sim_evap_Mm3 = sim_evap_Mm3,
         mb_flux = mb_flux),
    nrow = 1L
  )
}

# Apply physical constrains to prevent overflow and negative storage
physical_constrains <- function(storage,
                                max_capacity,
                                prov_release){
  # Calculate spill if overflow
  spill <- max(0, storage - max_capacity) 
  # Add to current release
  release <- prov_release + spill
  # Subtract spill to get final storage
  final_storage <- max(0, storage - spill)
  
  tibble::new_tibble(
    list(
      final_release = release,
      final_storage = max(0, storage - spill)
    ),
    nrow = 1L
  )
  
}

###### 2.1.2. Other functions (helpers) ######

obs_sim_tab <- function(sim_data,
                        obs_data){
  
  sim <- sim_data %>% 
    select(date, sim_storage_Mm3, sim_outflow_Mm3)
  
  obs <- obs_data %>% 
    select(date, outflow, storage) %>% 
    mutate(obs_outflow_Mm3 = outflow*86400/1e6, obs_storage = storage) %>% 
    select(date, obs_outflow_Mm3, obs_storage)
  
  left_join(sim, obs, by = 'date')
  
}

r2_cor <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred)
  if (sum(ok) < 2) return(NA_real_)
  cor(obs[ok], pred[ok])^2
}

##### 2.2. Specific functions #####
###### 2.2.1. Release policie functions
natural_lake_sim <- function(grand_id,
                             res_tseries,
                             sim_period,
                             r_factor_C = pars_natural_lake$r_factor_C, # Fraction of the storage to release
                             exp_C = pars_natural_lake$exp_C, # Exponent of non-linear function
                             pcp_fix = 0.6,
                             err_fix = TRUE) {
  
  # Preparing input data
  precip_fix <- pcp_fix
  mberr_fix <- err_fix
  
  # Input series
  full_series <- res_tseries[res_tseries$GRAND_ID == grand_id,]
  
  # Maximum capacity of the reservoir
  max_cap <- reservoir_properties[reservoir_properties$GRAND_ID == grand_id, 'Cap_mcm'] %>% pull(.)
  
  # Hypsometric curve
  hypso_curve <- hypso_curves[hypso_curves$GRAND_ID == grand_id, ]
  
  # Get storage for the beginning of simulation
  min_year <- min(sim_period)
  initial_stor <- full_series[year(full_series$date) == min_year - 1, ] %>%
    slice_tail(n = 1) %>%
    pull(storage)
  
  # In case no data for previous year, use current observed value
  initial_stor <- if (length(initial_stor) == 0 || is.na(initial_stor)) {
    full_series[year(full_series$date) == min_year & !is.na(full_series$storage), ] %>%
      slice_head(n = 1) %>%
      pull(storage)
  } else {
    initial_stor
  }
  
  # Time series for simulation
  tseries <- full_series[year(full_series$date) %in% sim_period, ]
  
  # Simulation
  # Simulated time steps
  tsteps <- length(tseries$date)
  
  sim_storage <- numeric(tsteps)
  sim_outflow <- numeric(tsteps)
  sim_area <- numeric(tsteps)
  sim_pcp <- numeric(tsteps)
  sim_evap <- numeric(tsteps)
  
  # Tibble to fill with simulations
  sim_series <-
    # Input vars
    tibble(
      date = tseries$date,
      inflow_Mm3_day = tseries$inflow * 86400 / 1e6,
      evap_mm = tseries$ET_mm,
      pcp_mm = tseries$PCP_mm,
      error_balance_Mm3 = tseries$daily_err_Mm3
    ) %>%
    # Variable to fill with simulations
    mutate(
      sim_storage_Mm3 = NA_real_,
      sim_outflow_Mm3 = NA_real_,
      sim_area_km2 = NA_real_,
      sim_evap_Mm3 = NA_real_,
      sim_pcp_Mm3 = NA_real_,
      GRAND_ID = grand_id
    ) %>%
    mutate(error_balance_Mm3 = if_else(is.na(error_balance_Mm3), 0, error_balance_Mm3))
  
  # Run simulation for each time step
  for (tstp in 1:tsteps) {
    tstep <- tstp
    
    # Storage to determine release --> Previous storage + Current inflow
    if (tstep == 1) {
      init_stor <- initial_stor + sim_series$inflow_Mm3_day[tstep]
    } else {
      # Others --> Previous simulated value
      init_stor <- sim_storage[tstep - 1] + sim_series$inflow_Mm3_day[tstep]
    }
    
    # In fraction
    stor_frac <- max(0, init_stor / max_cap)
    
    # Natural lake release function --> In-lined for faster run
    # natural_lake_release <- function(stor,
    #                                  stor_frac,
    #                                  r_factor_C,
    #                                  exp_C) {
    #   release <- r_factor_C * stor * stor_frac^exp_C
    #   min(stor, release) # If no available water to full release --> Release available
    # }
    # release_Mm3 <- natural_lake_release(stor = init_stor,
    #                                     stor_frac = stor_frac,
    #                                     r_factor_C = r_factor_C,
    #                                     exp_C = exp_C)
    
    # Calculate release (inlined function)
    release <- min(
      init_stor, # If no available water to full release --> Release available
      (r_factor_C * init_stor * stor_frac^exp_C) # Regular release
    )
    
    # Update storage with release
    init_stor <- max(0, init_stor - release)
    
    # Calculate P and E and mass balance flux
    mb_flux <- mb_fluxes(
      stor = init_stor,
      hypso = hypso_curve,
      pcp = sim_series$pcp_mm[tstep],
      evap = sim_series$evap_mm[tstep],
      err = sim_series$error_balance_Mm3[tstep],
      pcp_fix = precip_fix,
      err_fix = mberr_fix
    )
    
    # Update storage with fluxes
    init_stor <- init_stor + mb_flux$mb_flux
    
    # Apply physical constrains to get final release and storage
    constrains <- physical_constrains(
      storage = init_stor,
      max_capacity = max_cap,
      prov_release = release
    )
    
    sim_area[tstep] <- mb_flux$area_km
    sim_pcp[tstep] <- mb_flux$sim_pcp_Mm3
    sim_evap[tstep] <- mb_flux$sim_evap_Mm3
    sim_storage[tstep] <- constrains$final_storage
    sim_outflow[tstep] <- constrains$final_release
  }
  
  # Filling the tibble with the simulated time step
  sim_series$sim_area_km2 <- sim_area
  sim_series$sim_pcp_Mm3 <- sim_pcp
  sim_series$sim_evap_Mm3 <- sim_evap
  sim_series$sim_storage_Mm3 <- sim_storage
  sim_series$sim_outflow_Mm3 <- sim_outflow
  
  return(sim_series)
}

###### 2.2.2. Functions used for optimization ######

# Function to run the optimization
fitting_natural_lake <- function(grand_id, sim_period){
  
  # Used in obj. function
  obs_data <- reservoirs_time_series %>% 
    filter(GRAND_ID == grand_id,
           year(date) %in% sim_period)
  
  # Function to use inside deoptim --> Run the simulation and calculate obj. function
  obj_fn_natlake <- function(params) {
    
    # Parameter names --> Order!
    r_factor_C <- params[1]
    exp_C <- params[2]
    
    # Run simulation (onjects should be on environment)
    sim <- tryCatch(
      natural_lake_sim(grand_id = grand_id,
                       res_tseries = obs_data,
                       sim_period = sim_period, 
                       r_factor_C = r_factor_C, 
                       exp_C = exp_C),
      error = function(e) NULL
    )
    # If error, return 1e6
    if (is.null(sim)) return(1e6)
    
    # Obs. sim. tab
    tab <- tryCatch(obs_sim_tab(sim_data = sim, 
                                obs_data = obs_data),
                    error = function(e) NULL)
    # If error, return 1e6
    if (is.null(tab)) return(1e6)
    
    # Valid rows
    ok_out <- is.finite(tab$sim_outflow_Mm3) & is.finite(tab$obs_outflow_Mm3)
    ok_sto <- is.finite(tab$sim_storage_Mm3) & is.finite(tab$obs_storage)
    
    if (sum(ok_out) < 2 || sum(ok_sto) < 2) return(1e6)
    
    # nRMSE for each variable
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
    
    # Objective funcion value
    0.5 * nrmse_out + 0.5 * nrmse_sto
  }
  
  fit <- DEoptim(fn = obj_fn_natlake, 
                 lower = low_bound, 
                 upper = upp_bound,
                 control = DEoptim_conf)
  
  best_pars <- tibble(!!parms[1] := fit$optim$bestmem[1], 
                      !!parms[2] := fit$optim$bestmem[2])
  
  perf_evol <- tibble(iteration = seq_along(fit$member$bestvalit),
                      obj_f = fit$member$bestvalit)
  
  best_ind <- tibble(iteration = seq_len(nrow(fit$member$bestmemit)),
                     !!parms[1] := fit$member$bestmemit[,1],
                     !!parms[2] := fit$member$bestmemit[,2])
  
  all_pop <- purrr::map_dfr(seq_along(fit$member$storepop), function(i) {
    tt <- tibble::as_tibble(fit$member$storepop[[i]])
    colnames(tt) <- parms
    dplyr::mutate(tt, iteration = i, .before = 1)
  })
  
  list(best_pars = best_pars,
       perf_evol = perf_evol,
       best_ind = best_ind,
       all_pop = all_pop)
  
}


#### 3. Running simulations for one reservoir #### 
##### 3.1. Default simulation #####
gid <- 132

# Observed data
obs_data <- reservoirs_time_series %>%
  filter(GRAND_ID == gid)

# Simulation parameters
parms <- c('r_factor_C', 'exp_C')
def_parms <- tibble(0.01, 1.5)
colnames(def_parms) <- parms

# Run simulation
def_sim <- natural_lake_sim(grand_id = gid, 
                 res_tseries = reservoirs_time_series, 
                 sim_period = 1990:2015, 
                 r_factor_C = def_parms$r_factor_C, 
                 exp_C = def_parms$exp_C)

# Obs.sim. table
obssim_tab <- obs_sim_tab(sim_data = def_sim,
            obs_data = obs_data)

# Performance
obssim_tab %>% 
summarise(KGE_stor = KGE(sim_storage_Mm3, obs_storage),
          KGE_out = KGE(sim_outflow_Mm3, obs_outflow_Mm3))


# Plot
obssim_tab %>% 
  pivot_longer(-date) %>% 
  mutate(var = if_else(str_detect(name, 'out'), 'Release (Mm3/d)', 'Storage (Mm3)'),
         type = if_else(str_detect(name, 'sim'), 'Simulated', 'Observed')) %>% 
  ggplot(., aes(x = date, y = value, color = type))+
  geom_line()+facet_wrap(~var, scales  = 'free', ncol = 1)+
  scale_color_manual(values = c('black', 'steelblue'))+
  theme_bw()+ 
  theme(axis.title = element_blank(), legend.title = element_blank())


##### 3.2. Parameters optimization #####

gid <- 132

# Observed data
obs_data <- reservoirs_time_series %>%
  filter(GRAND_ID == gid)

# Simulation parameters
parms <- c('r_factor_C', 'exp_C')
def_parms <- tibble(0.01, 1.5)
colnames(def_parms) <- parms


n_pars <- length(parms) # r_factor_c and exp_C
max_itrs <- 200
low_bound <- c(0.001, 0.5)
upp_bound <- c(0.1, 2.5)

# DE adjustments --> NP, F, CR default
DEoptim_conf <- DEoptim.control(NP = 10*n_pars, # Number of individuals first pop
                           itermax = max_itrs, # Max. number of sims. per generation
                           F = 0.8, # Mutation extent
                           CR = 0.9, # Mutation probability
                           trace = FALSE, 
                           # To converge:
                           reltol = 1e-4, # Minimum obj. funct. improvement
                           steptol = 10, # In last n iterations 
                           # Keep obj. functions and population info
                           storepopfrom = 1,
                           storepopfreq = 1)

# Run optimization
deoptim_fit <- fitting_natural_lake(grand_id = gid, 
                                    sim_period = 1990:2004)


##### 3.3. Optimized simulation #####
opt_pars <- deoptim_fit$best_pars


# Run simulation
opt_sim <- natural_lake_sim(grand_id = gid, 
                            res_tseries = reservoirs_time_series, 
                            sim_period = 1990:2015, 
                            r_factor_C = opt_pars$r_factor_C, 
                            exp_C = opt_pars$exp_C)

# Obs.sim. table
obssim_tab <- obs_sim_tab(sim_data = opt_sim,
                          obs_data = obs_data)

# Performance
obssim_tab %>% 
  summarise(KGE_stor = KGE(sim_storage_Mm3, obs_storage),
            KGE_out = KGE(sim_outflow_Mm3, obs_outflow_Mm3))


# Plot
obssim_tab %>% 
  pivot_longer(-date) %>% 
  mutate(var = if_else(str_detect(name, 'out'), 'Release (Mm3/d)', 'Storage (Mm3)'),
         type = if_else(str_detect(name, 'sim'), 'Simulated', 'Observed')) %>% 
  ggplot(., aes(x = date, y = value, color = type))+
  geom_line()+facet_wrap(~var, scales  = 'free', ncol = 1)+
  scale_color_manual(values = c('black', 'steelblue'))+
  theme_bw()+ 
  theme(axis.title = element_blank(), legend.title = element_blank())


rbind(
(opt_sim %>% 
  select(date, sim_storage_Mm3, sim_outflow_Mm3) %>% 
  mutate(pars = 'default')),
(def_sim %>% 
  select(date, sim_storage_Mm3, sim_outflow_Mm3) %>% 
  mutate(pars = 'fitted'))) %>% 
  left_join(., obs_data[, c('date', 'outflow', 'storage')]) %>% 
  mutate(obs_storage = storage, obs_outflow_Mm3 = outflow*86400/1e6) %>% 
  select(date:pars, starts_with('obs')) %>% 
  pivot_longer(., -c(date, pars)) %>% 
  mutate(var = if_else(str_detect(name, 'out'), 'Release (Mm3/d)', 'Storage (Mm3)'),
         type = if_else(str_detect(name, 'sim'), 'Simulated', 'Observed')) %>%
  mutate(type = case_when(type == 'Observed' ~ 'Observed',
                          pars == 'default'~ 'Sim. default',
                          .default = 'Sim. fitted')) %>% 
  ggplot(., aes(x = date, y = value, color = type))+
  geom_line()+facet_wrap(~var, scales  = 'free', ncol = 1)+
  scale_color_manual(values = c('black', 'steelblue', 'red4'))+
  theme_bw()+ 
  theme(axis.title = element_blank(), legend.title = element_blank())

