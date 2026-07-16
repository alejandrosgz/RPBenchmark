library(tidyverse)
library(sf)
library(slider)
library(climateR)
library(readxl)
library(tidymodels)
library(patchwork)


r2_cor <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred)
  if (sum(ok) < 2) return(NA_real_)
  cor(obs[ok], pred[ok])^2
}

studied_period <- c(1989:2014)


gdw_data <- read_sf('C:/ASG/UCDavis/NSF_PROJECT_FOLDER/Datasets/GlobalDamWatcher/CONUS_reservoirs.shp') %>% 
  st_drop_geometry(.)

hucs <- read_sf('C:/ASG/UCDavis/NSF_PROJECT_FOLDER/Datasets/SpatialData/HUCS/HUCS_02_CONUS.shp')
reservoirs_polygons <- read_sf('C:/ASG/UCDavis/NSF_PROJECT_FOLDER/Datasets/GlobalDamWatcher/CONUS_reservoirs.shp')


# Table with reservoir properties
reservs_data <- reservoirs_polygons %>% 
  select(GDW_ID,
         GRAND_ID,
         RES_NAME,     # Name of the reservoir or impounded waterbody
         DAM_NAME,     # Name of the dam or barrier structure
         CAP_MCM,      # Representative max storage capacity in million cubic meters
         DIS_AVG_LS,   # Average long-term discharge at the barrier location (liters/sec)
         DOR_PC,       # Degree of regulation: ratio of capacity to annual flow (%)
         ELEV_MASL,    # Elevation of the reservoir surface (meters above sea level)
         RIVER,
         CATCH_SKM,    # Upstream catchment area draining into the reservoir (sq. km)
         starts_with("USE"), # Reservoir uses: irrigation, electricity, supply, etc. (multiple columns)
         MAIN_USE,      # Main stated purpose of the reservoir (e.g., irrigation, flood control)
         LONG_DAM, 
         LAT_DAM
  ) %>% 
  # Keeping the reservoir name
  mutate(Name = if_else(!is.na(RES_NAME),
                        RES_NAME,
                        DAM_NAME),
         .before = 2) %>% 
  # 2 had no names, I search in Gmaps and I use a name of close features
  mutate(Name = case_when(GDW_ID == 36198 ~ "FMI Henderson mining raft",
                          GDW_ID == 7978 ~ "Redrock Canyon",
                          .default = Name)) %>% 
  select(-c(RES_NAME, DAM_NAME)) %>% 
  # To simplify the use information
  mutate(USE_IRRI = if_else(is.na(USE_IRRI), "-", paste0("Irrigation", " (", USE_IRRI, ")")),
         USE_ELEC = if_else(is.na(USE_ELEC), "-", paste0("Hydropower", " (", USE_ELEC, ")")),
         USE_SUPP = if_else(is.na(USE_SUPP), "-", paste0("Water supply", " (", USE_SUPP, ")")),
         USE_FCON = if_else(is.na(USE_FCON), "-", paste0("Flood control", " (",USE_FCON, ")")),
         USE_RECR = if_else(is.na(USE_RECR), "-", paste0("Recreation", " (",USE_RECR, ")")),
         USE_NAVI = if_else(is.na(USE_NAVI), "-", paste0("Navigation", " (",USE_NAVI, ")")),
         USE_FISH = if_else(is.na(USE_FISH), "-", paste0("Fisheries", " (",USE_FISH, ")")),
         USE_PCON = if_else(is.na(USE_PCON), "-", paste0("Pollution control", " (",USE_PCON, ")")),
         USE_LIVE = if_else(is.na(USE_LIVE), "-", paste0("Livestock water supply", " (",USE_LIVE, ")")),
         USE_OTHR = if_else(is.na(USE_OTHR), "-", paste0("Other", " (",USE_OTHR, ")"))) %>% 
  mutate(Uses = paste(USE_IRRI, USE_ELEC, USE_SUPP, USE_FCON, USE_RECR, USE_NAVI, USE_FISH, USE_PCON, USE_LIVE, USE_OTHR, sep = ", "),
         .before = 3) %>% 
  mutate( Uses = str_remove_all(Uses, "-, "),
          Uses = str_remove_all(Uses, ", -")) %>% 
  select(-c(USE_IRRI:USE_OTHR)) %>% 
  # 2 columns kept: Once with all the uses, the main the first. The other a column with only the main
  mutate(
    Uses = map_chr(str_split(Uses, ",\\s*"), \(parts) {
      parts <- str_trim(parts)
      mains <- parts[str_detect(parts, "\\(Main\\)")]
      secs  <- parts[!str_detect(parts, "\\(Main\\)")]
      paste(c(mains, secs), collapse = ", ")
    })) %>% 
  st_intersection(., hucs[, c("HUC2", "NAME")]) %>% 
  mutate(area_sqkm = round(as.numeric(st_area(.)/1e6), 3),
         perim_km = round(as.numeric(st_perimeter(.)/1e3), 3)) %>% 
  select(-CATCH_SKM) %>% 
  select(GRAND_ID, GDW_ID, Name, MAIN_USE, Uses, 
         CAP_MCM, RIVER, DIS_AVG_LS, DOR_PC,  
         LONG_DAM, LAT_DAM, ELEV_MASL, area_sqkm, perim_km, HUC2, geometry)


colnames(reservs_data) <- c("GRAND_ID", "GDW_ID", "Name", 'Main_use', "Uses", "Cap_mcm", 'River', "Disch_ls", 
                            "Reg_degre", 'Long', 'Lat', "Elev", "area_km", "perm_km", 'HUC2', "geometry")



gridmet_et_data <- read_csv('data/Processed_data/gridmet_downloaded_series/gridmet_et_download.csv')
gridmet_pcp_data <- read_csv('data/Processed_data/gridmet_downloaded_series/gridmet_pcp_download.csv')

studied_reservoirs <- read_table('data/studied_reservoirs.txt')$studied_reservoirs

reservoirs_properties <- list()
reservoirs_available_data <- list()
reservoirs_clean_data <- list()
hypsometrics_curves <- list()
reservoirs_area_series <- list()
reservoirs_evaporation_series <- list()
reservoirs_precipitation_series <- list()


time_series_simulations <- list()

initial_storages <- list()
reservoir_lack_data_daily <- list()



for(i in 1:length(studied_reservoirs)){
  
  
  it_n <- i
  
  res_ID_grandID <- studied_reservoirs[i]
  
  res_properties <- reservs_data[reservs_data$GRAND_ID == res_ID_grandID,]
  
  resname <- res_properties$Name
  
  # Reservoir properties and shapefile
  
  r_gdw <- gdw_data[gdw_data$GRAND_ID == res_ID_grandID,]
  r_pol <- reservoirs_polygons[reservoirs_polygons$GRAND_ID == res_ID_grandID, c('GDW_ID', 'RES_NAME', 'GRAND_ID')]
  
  
  # ResOpsUS time series
  
  resops_ts_raw <- read_csv(paste0('C:/ASG/UCDavis/NSF_PROJECT_FOLDER/Datasets/ResOpsUS/time_series_all/ResOpsUS_', res_ID_grandID, '.csv')) %>% 
    filter(year(date) %in% studied_period) %>% 
    select(date:outflow)
  
  
  # Storage for these reservoirs?? --> Seems there is an unit error in ResOpsUS
  if(res_ID_grandID %in% c(767, 768)){
    
    # These reservoirs have wrong storage units in ResOpsUS
    resops_ts_raw$storage <- 1e3*resops_ts_raw$storage
  }
  
  # Inflow in this reservoir looks not ok
  if(res_ID_grandID == 1916){
    
    resops_ts_raw$inflow <- NA
    
  }
  
  # Removing negative values
  resops_ts_raw[resops_ts_raw$inflow < 0 &
                  !is.na(resops_ts_raw$inflow),]$inflow <- NA
  resops_ts_raw[resops_ts_raw$outflow < 0 &
                  !is.na(resops_ts_raw$outflow),]$outflow <- NA
  resops_ts_raw[resops_ts_raw$storage < 0 &
                  !is.na(resops_ts_raw$storage),]$storage <- NA
  
  # Detecting/Removing outliers
  
  # Storage outliers detection
  
  # Determining Reservoir total capacity
  obs_max <- max(resops_ts_raw$storage, na.rm = T)
  q95_obs <- quantile(resops_ts_raw$storage, 0.95, na.rm = T)
  q99_obs <- quantile(resops_ts_raw$storage, 0.99, na.rm = T)
  gdw_max <- r_gdw$CAP_MAX
  
  max_tib <- tibble(obs_max,  q95_obs, q99_obs, gdw_max)
  
  remove(obs_max, q95_obs, q99_obs, gdw_max)
  
  max_cap_res <- max_tib %>% 
    mutate(obs_max = if_else(obs_max > 5*q95_obs, NA, obs_max), # Maximum observed not reliable if larger than 5 times q95
           gdw_max = if_else(gdw_max > 5*q95_obs, NA, gdw_max), # GDW not reliable if larger than 5 times q95
           gdw_max = if_else(gdw_max < q95_obs, NA, gdw_max), # GDW not reliable if smaller than q95
           q99_obs = if_else(q99_obs > 5*q95_obs, NA, q99_obs)) %>% # q99 not reliable if larger than 5 times q95
    mutate(Max_cap = max(obs_max, q95_obs, q99_obs, gdw_max, na.rm = T)) %>% 
    pull(Max_cap)

  # In case all NAs --> use q95  
  if(!is.finite(max_cap_res)){
    max_cap_res <- max_tib$q95_obs
  }
  
  res_properties$Cap_mcm <- max_cap_res
  
  # List of reservoirs properties
  reservoirs_properties[[it_n]] <- res_properties
  
  # Thresholds
  pervar <- 3 # Storage variation (% relative to max capacity) considered large
  ntimes_stor_inout <- 2 # ntimes storage variation / inf out difference to suspect
  
  stor_out_det <- 
    resops_ts_raw %>% 
    arrange(., date) %>% 
    mutate(stor_var_perc = 100*(storage-lag(storage))/max_cap_res, # Relative variation 
           stor_var_hm = storage-lag(storage),
           inoutdiff_hm = 86400*(inflow-outflow)/1e6) %>% 
    mutate(largevar = abs(stor_var_perc) > pervar, # Large variation
           vlargevar = abs(stor_var_perc) > 2*pervar, # Very large variation
           pot_err = abs(inoutdiff_hm) > abs(stor_var_hm)*ntimes_stor_inout) %>% 
    mutate(error_detected = case_when(largevar == T & pot_err == T ~ 'out',
                                      vlargevar == T & is.na(pot_err) ~ 'out',
                                      .default = 'norm')) %>% 
    mutate(error_detected = if_else(lag(error_detected) == 'out' |
                                      lag(error_detected, 2) == 'out' |
                                      # lag(error_detected, 3) == 'out' |
                                      # lag(error_detected, 4) == 'out' |
                                      lead(error_detected) == 'out' |
                                      lead(error_detected, 2) == 'out',#|
                                    # lead(error_detected, 3) == 'out'|
                                    # lead(error_detected, 4) == 'out',
                                    'out', error_detected)) %>% 
    mutate(error_detected = if_else(is.na(error_detected), 'norm', error_detected)) %>% 
    ungroup(.)
  
  stor_data_avail <- stor_out_det %>% 
    mutate(storage = if_else(error_detected == 'out', NA, storage)) %>% 
    mutate(obs_S = if_else(!is.na(storage), 1, 0)) %>% 
    select(date, obs_S) %>% 
    mutate(GRAND_ID = res_ID_grandID)
  
  
  # Storage outliers removal (linear regression previous to next available data)
  resops_stor_fixed <- stor_out_det %>% 
    mutate(storage = if_else(error_detected == 'out', NA, storage)) %>% 
    mutate(storage = zoo::na.approx(storage, na.rm = F)) %>% 
    select(date, storage) %>% 
    mutate(GRAND_ID = res_ID_grandID)
  
  
  # Inflow and Outflow outliers detection
  
  outliers_threshold <- resops_ts_raw %>% 
    summarise(q99_i = quantile(inflow, 0.99, na.rm = TRUE),
              q99_o = quantile(outflow, 0.99, na.rm = TRUE)) %>% 
    mutate(inf_outlier_threshold = 5*q99_i,
           out_outlier_threshold = 5*q99_o) %>% 
    select(ends_with('hold')) %>% 
    mutate(GRAND_ID = res_ID_grandID)
  
  InOut_out_det <- resops_ts_raw %>% 
    mutate(GRAND_ID = res_ID_grandID) %>% 
    left_join(., outliers_threshold) %>% 
    mutate(inflow = if_else(inflow > inf_outlier_threshold, NA, inflow),
           outflow = if_else(outflow > out_outlier_threshold, NA, outflow),
           inflow = if_else(inflow > 15000, NA, inflow)) %>% # Extremely high, not realistic in the studied reservoirs
    select(GRAND_ID, date, inflow, outflow)
  
  InOut_data_avail <- InOut_out_det %>% 
    mutate(obs_I = if_else(!is.na(resops_ts_raw$inflow), 1, 0),
           obs_O = if_else(!is.na(resops_ts_raw$outflow), 1, 0)) %>% 
    select(date, GRAND_ID, obs_I, obs_O)
  
  # Locating lack data (1 observed, 0 infered)
  dat_avail_res <- left_join(stor_data_avail, InOut_data_avail) %>% 
    select(GRAND_ID, date, obs_I, obs_O, obs_S)
  
  
  reservoir_lack_data_daily[[it_n]] <- dat_avail_res
  
  
  # Fixing lack of data and outliers for Inflow and outflow
  
  # Short gaps (5 days or less --> Linear regression from last to next available data)
  inout_shortgaps_fixed <- InOut_out_det %>% 
    mutate(
      inflow  = if(all(is.na(inflow))) inflow else zoo::na.approx(inflow, maxgap = 5, na.rm = FALSE),
      outflow = if(all(is.na(outflow))) outflow else zoo::na.approx(outflow, maxgap = 5, na.rm = FALSE)
    ) %>% 
    mutate(
      inflow  = if(all(is.na(inflow))) inflow else zoo::na.locf(inflow, na.rm = FALSE),
      outflow = if(all(is.na(outflow))) outflow else zoo::na.locf(outflow, na.rm = FALSE)
    ) %>%
    mutate(
      inflow  = if(all(is.na(inflow))) inflow else zoo::na.locf(inflow, fromLast = TRUE),
      outflow = if(all(is.na(outflow))) outflow  else zoo::na.locf(outflow, fromLast = TRUE)
    )
  
  # Large gaps (> 5 days)
  
  if(res_ID_grandID %in% studied_reservoirs){
    
    # Use water balance when possible
    recalc_res_data <- inout_shortgaps_fixed %>%
      left_join(resops_stor_fixed) %>%
      mutate(avail_S = !is.na(storage),
             avail_I = !is.na(inflow),
             avail_O = !is.na(outflow)) %>%
      mutate(nvars_avail = avail_S + avail_I + avail_O) %>%
      mutate(storvar_day_m3s = ((storage-lag(storage))*1e6)/86400) %>%
      
      mutate(inflow = case_when(!is.na(inflow) ~ inflow,
                                is.na(inflow) & nvars_avail == 2 ~ storvar_day_m3s + outflow,
                                .default = NA),
             outflow = case_when(!is.na(outflow) ~ outflow,
                                 is.na(outflow) & nvars_avail == 2 ~ inflow - storvar_day_m3s,
                                 .default = NA))
    
    # When not possible, use na.approx
    recalc_res_data <- recalc_res_data %>%
      mutate(inflow = zoo::na.approx(inflow, na.rm = FALSE),
             outflow = zoo::na.approx(outflow, na.rm = FALSE)) %>%
      mutate(inflow = if_else(inflow < 0, 0, inflow),
             outflow = if_else(outflow < 0, 0, outflow)) 
    
    
    # This reservoir lacks storage from 1990
    if(res_ID_grandID == 1127){
      
      last_stor <- recalc_res_data$storage[which(is.na(recalc_res_data$storage))[1]-1]
      
      recalc_res_data <-  recalc_res_data %>% 
        mutate(cum_inout = cumsum(if_else(avail_S == T, 0, 0.0864*(inflow-outflow)))) %>% 
        mutate(storage = last_stor + cum_inout) %>% 
        select(GRAND_ID:storvar_day_m3s)
      
    }
    
    recalc_res_data <- recalc_res_data %>%
      # Final NAs filled use closest value
      fill(inflow, storage, .direction = "down") %>%
      fill(inflow, storage, .direction = "up") %>% 
      select(GRAND_ID, date, inflow, outflow, storage)
    
    # Sould never happen
  } else { # For reservoirs included in the project
    print('ERROR')
  }
  
  # List of reservoirs clean data
  reservoirs_clean_data[[it_n]] <- recalc_res_data
  
  # List of reservoirs available data
  dat_avail_sum <- dat_avail_res %>% 
    group_by(GRAND_ID) %>% 
    summarise(avail_inf = sum(obs_I)/length(date),
              avail_out = sum(obs_O)/length(date),
              avail_stor = sum(obs_S)/length(date),
    )
  
  reservoirs_available_data[[it_n]] <- dat_avail_sum
  
  
  # Time series fix plots
  dat_avail_long <- dat_avail_res %>% 
    pivot_longer(., -c(GRAND_ID, date),
                 names_to = 'var',
                 values_to = 'obs') %>% 
    mutate(var = case_when(str_detect(var, 'I') ~ 'inflow',
                           str_detect(var, 'O') ~ 'outflow',
                           str_detect(var, 'S') ~ 'storage'))
  
  series_fixed_long <- recalc_res_data %>%
    pivot_longer(., -c(GRAND_ID, date),
                 names_to = 'var',
                 values_to = 'value')
  
  ts_fix_plot <-  left_join(dat_avail_long, series_fixed_long) %>% 
    ggplot(., aes(x = date, y = value, color = factor(obs)))+
    geom_point()+
    facet_wrap(~var, ncol = 1, scales = 'free')+
    theme_bw()+
    scale_color_manual(
      values = c("1" = "steelblue", "0" = "brown")) +
    labs(x = 'Date',
         y = 'Value', 
         title = paste0('Time series fixed', ' GRAND_ID = ', res_ID_grandID, ', ', resname, ' (', round(max_cap_res, 0), ' Mm3)'))+
    theme(legend.position = 'nuop')
  
  ggsave(ts_fix_plot,
         filename = paste0('figures/time_series_figs/1_cleaning_series/', it_n, '_ts_clean_', res_ID_grandID, '.png'),
         device = 'png',
         width = 10, height = 8,
         dpi = 300)
  
  
  # Hypsometry data
  
  base_dir <- "C:/ASG/UCDavis/NSF_PROJECT_FOLDER/Datasets/Hypsometric_data/ZaoGao19/Reservoir_evaporation721/"
  
  file <- paste0(base_dir, res_ID_grandID, ".txt")
  exists_vec <- file.exists(file)
  
  
  # Only reservoirs in Zao Gao to be consistent 
  if(exists_vec){ # If available in Zao and Gao 2019, use it
    
    areas_gao1 <- read_table(paste0(base_dir, res_ID_grandID, ".txt")) %>%
      mutate(
        GRAND_ID = res_ID_grandID,
        Areas = as.numeric(Areas)) %>% 
      mutate(monyear = paste(year(ymd(Month)), month(ymd(Month)), sep = '/')) %>% 
      select(GRAND_ID, Areas, monyear) 
    
    month_storage_series <- recalc_res_data %>% 
      group_by(month(date), year(date), GRAND_ID) %>% 
      summarise(storage = mean(storage, na.rm = T)) %>% 
      ungroup(.) %>% 
      mutate(monyear = paste(`year(date)`, `month(date)`, sep = '/')) %>% 
      select(GRAND_ID, monyear, storage)
    
    area_stor_relationships <- areas_gao1 %>% 
      left_join(month_storage_series, .,c('GRAND_ID', 'monyear')) %>%
      mutate(Areas = round(Areas, 1)) %>%
      group_by(Areas, GRAND_ID) %>% 
      summarise(storage = mean(storage)) %>% 
      # mutate(Area_stor_rel = Areas/storage) %>% 
      mutate(
        Area_stor_rel = if_else(storage > 0, Areas / storage, NA_real_)
      ) %>% 
      group_by(GRAND_ID) %>% 
      mutate(z_area_stor = scale(Area_stor_rel)) %>% 
      # mutate(
      #   z_area_stor = if (sd(Area_stor_rel, na.rm = TRUE) > 0)
      #     as.numeric(scale(Area_stor_rel))
      #   else NA_real_
      # ) %>% 
      filter(., abs(z_area_stor) <= 2.5) %>% 
      select(-c(z_area_stor, Area_stor_rel))
    
    
    poly_lm_curves <- area_stor_relationships %>% 
      group_by(GRAND_ID) %>% 
      nest(data = c(storage, Areas)) %>% 
      mutate(
        model = map(
          data,
          ~ lm(Areas ~ poly(storage, 2, raw = TRUE), data = .x)
        )
      ) %>%
      mutate(
        coefs = map(model, broom::tidy),
        fit   = map(model, broom::glance)
      ) %>%
      unnest(coefs, fit) %>%
      select(-c(r.squared, std.error, statistic, p.value, sigma, statistic1, p.value1, df, logLik, AIC,
                BIC, deviance, df.residual,  nobs)) %>%
      pivot_wider(names_from = term, values_from = estimate) %>%
      rename(intercept = `(Intercept)`, 
             poly1 = `poly(storage, 2, raw = TRUE)1`,
             poly2 = `poly(storage, 2, raw = TRUE)2`) %>% 
      select(GRAND_ID, intercept, poly1, poly2, adj.r.squared, data)
    
    
    
    # Method added to improve hypsometric curves
    
    if(poly_lm_curves$poly2 > 0){ # Wrong curvature, use linear model
      
      poly_lm_curves <- area_stor_relationships %>% 
        group_by(GRAND_ID) %>% 
        nest(data = c(storage, Areas)) %>% 
        mutate(
          model = map(
            data,
            ~ (lm(Areas ~ storage, data = .x))
          )) %>%
        mutate(
          coefs = map(model, broom::tidy),
          fit   = map(model, broom::glance)
        ) %>%
        unnest(coefs, fit) %>%
        select(-c(r.squared, std.error, statistic, p.value, sigma, statistic1, p.value1, df, logLik, AIC,
                  BIC, deviance, df.residual,  nobs)) %>%
        pivot_wider(names_from = term, values_from = estimate) %>%
        rename(intercept = `(Intercept)`, 
               poly1 = `storage`) %>% 
        mutate(poly2 = 0) %>% 
        select(GRAND_ID, intercept, poly1, poly2, adj.r.squared, data)
      
      max_stor_reliable <- NA
      
      poly_lm_curves <- poly_lm_curves %>% 
        mutate(max_stor_rel = max_stor_reliable)
      
      
    } else { # Good curvature, check if fails at large values
      
      minvol <- quantile(recalc_res_data$storage, 0.01, na.rm = T)
      maxvol <- quantile(recalc_res_data$storage, 0.99, na.rm = T)
      
      stor_tries <- seq(minvol, maxvol, length.out = 3000)
      
      fitted <- poly_lm_curves$intercept + stor_tries*poly_lm_curves$poly1 + poly_lm_curves$poly2*stor_tries^2
      
      firstneg <- which((diff(fitted) / diff(stor_tries))<0)[1]
      
      if(is.na(firstneg)){ # No negatives, the curve is ok
        
        poly_lm_curves <- poly_lm_curves %>% 
          mutate(max_stor_rel = NA)
        
        
      } else { # See when it fails
        
        max_stor_reliable <- if(length(stor_tries[firstneg-1])==0){
          (stor_tries[firstneg])
        }else{
          (stor_tries[firstneg-1])}
        
        poly_lm_curves <- poly_lm_curves %>% 
          mutate(max_stor_rel = max_stor_reliable)
        
      }
      
    }
    

    # Hypso curve plot
    
    # Options
    # Wrong curve --> Use poly 2 only when storage is below max storage, Use max storage to calculate max area if is above  

    fitted_areas <- 
      area_stor_relationships %>% 
      left_join(poly_lm_curves, by = "GRAND_ID") %>% 
      mutate(area_km = 
               case_when(
                 
                 # Linear regression --> NA max storage (use poly 2 = to 0)
                 # Good curve --> NA max storage (use regular poly 2)
                 # Use regular equation
                 is.na(max_stor_rel) ~ intercept + poly1 * storage + poly2 * storage^2,
                 
                 # Wrong curve, but storage below maximum
                 # Use regular equation
                 storage <= max_stor_rel ~ intercept + poly1 * storage + poly2 * storage^2,
                 
                 # Wrong curve, and storage above maximum
                 # Use equation with maximum reliable storage
                 storage > max_stor_rel ~ intercept + poly1 * max_stor_rel + poly2 * max_stor_rel^2)
             
      )
    
    r2s <- fitted_areas %>%
      # unnest(data) %>% 
      group_by(GRAND_ID) %>% 
      summarise(xpos = quantile(storage, 0.35, na.rm = T),
                ypos = quantile(Areas, 0.95, na.rm = T),
                r2 = round(r2_cor(obs = Areas, area_km), 2))
    
    plylmmcurves <- fitted_areas %>%
      ggplot(., aes(y = Areas, x = storage))+
      geom_point(color = 'blue')+
      geom_line(aes(y = area_km),
                color = 'red')+
      geom_label(data = r2s, 
                 aes(label = paste0('R2 = ', r2), x = xpos, y = ypos),
                 fill = NA,        # no background
                 label.size = NA)+
      facet_wrap(facets = 'GRAND_ID',
                 scales = 'free', 
                 ncol = 5)+
      theme_bw()
    
    hypsometric_curve <- poly_lm_curves %>% 
      select(-data) %>% ungroup()

    # Should never happen
    
  } else {  # For reservoirs not included in Zao Gao --> InfeRes is used
    print('Error in hypsometry')
  }
  
  
  # List of hypsometrics curve
  hypsometrics_curves[[it_n]] <- hypsometric_curve
  
  
  ggsave(plylmmcurves,
         filename = paste0('figures/time_series_figs/2_hipso_curves/', it_n, '_hypso_curve_', res_ID_grandID, '.png'),
         device = 'png',
         width = 4, height = 4,
         dpi = 300)
  
  
  # Estimate reservoir areas from storage and hypsometric curve
  
  res_area <- recalc_res_data %>% 
    arrange(GRAND_ID, date) %>%
    left_join(hypsometric_curve, by = "GRAND_ID") %>% 
    group_by(GRAND_ID) %>% 
    mutate(
      area_km = 
        
        case_when(
          # Linear regression --> NA max storage (use poly 2 = to 0)
          # Good curve --> NA max storage (use regular poly 2)
          # Use regular equation
          is.na(max_stor_rel) ~ intercept + poly1 * storage + poly2 * storage^2,
          
          # Wrong curve, but storage below maximum
          # Use regular equation
          storage <= max_stor_rel ~ intercept + poly1 * storage + poly2 * storage^2,
          
          # Wrong curve, and storage above maximum
          # Use equation with maximum reliable storage
          storage > max_stor_rel ~ intercept + poly1 * max_stor_rel + poly2 * max_stor_rel^2)
      
    ) %>% 
    mutate(area_km = if_else(area_km < 0, NA_real_, area_km),
           
           d_area = abs(area_km - lag(area_km, default = first(area_km))),
           
           sum_2prev_2next = slide_dbl(
             d_area,
             ~ sum(.x, na.rm = TRUE),
             .before = 2, .after = 2,
             .complete = TRUE
           ),
           
           thr95 = quantile(sum_2prev_2next, 0.95, na.rm = TRUE),
           
           # flag spikes BUT never the first 2 rows in each GRAND_ID
           area_km = if_else(
             row_number() < 3,
             area_km,
             if_else(sum_2prev_2next > 2 * thr95 & !is.na(area_km), NA_real_, area_km)
           )
    ) %>% 
    ungroup() %>% 
    select(GRAND_ID, date, area_km) %>% 
    mutate(area_km = zoo::na.approx(area_km, na.rm = F))
  
  # List with estimated area series
  reservoirs_area_series[[it_n]] <- res_area
  
  
  resarea_plot <- res_area %>% 
    left_join(., recalc_res_data) %>% 
    select(date, area_km, storage) %>% 
    pivot_longer(-date) %>% 
    ggplot(., aes(x = date, y = value, color = name))+
    geom_line()+
    facet_wrap(~name, scales = 'free_y', ncol = 1)+
    labs(x = 'Date', y = 'Value', 
         title = paste0('Area and Storage series', 
                        ' GRAND_ID = ', res_ID_grandID, ', ', 
                        resname, ' (', round(max_cap_res, 0), ' Mm3)'))+
    
    theme_bw()+
    scale_color_manual(values = c('green4', 'orange3'))+
    theme(legend.position = 'nope')  
  
  ggsave(resarea_plot,
         filename = paste0('figures/time_series_figs/3_Area_series/', it_n, '_storarea_', res_ID_grandID, '.png'),
         device = 'png',
         width = 10, height = 6,
         dpi = 300)
  
  
  # Evaporation data
  
  # Better bulk download (saved in  gridmet_et_data)
  
  #start_date <- "1989-01-01"
  #end_date   <- "2020-12-31"
  #
  #  res_pol <- r_pol %>% 
  #    st_make_valid() %>%
  #    filter(!st_is_empty(.)) %>%
  #    st_collection_extract("POLYGON") %>%
  #    st_cast("MULTIPOLYGON", warn = FALSE)
  #  
  #  pet_res <- getGridMET(
  #    AOI = res_pol,
  #    varname = "pet",
  #    startDate = start_date,
  #    endDate  = end_date
  #  )
  #  
  #  # If getGridMET returns a tibble/data.frame (date, pet) ----
  #  if (inherits(pet_res, "data.frame") && all(c("date", "pet") %in% names(pet_res))) {
  #    
  #    daily_pet_w <- pet_res %>%
  #      transmute(
  #        date = as.Date(date),
  #        ET_mm = round(pet, 2),
  #        GRAND_ID = grand_id
  #      )
  #    
  #  } else {
  #    # If getGridMET returns a list with SpatRaster ----
  #    pet_r <- pet_res$daily_mean_reference_evapotranspiration_grass
  #    
  #    if (is.null(pet_r) || !inherits(pet_r, "SpatRaster")) {
  #      stop("Unexpected getGridMET output for GRAND_ID = ", grand_id)
  #    }
  #    
  #    # match CRS and convert to SpatVector for terra
  #    res_pol2 <- st_transform(res_pol, st_crs(pet_r))
  #    res_v    <- terra::vect(res_pol2)
  #    
  #    # crop to raster extent (helps when polygon barely overlaps / tiny rasters)
  #    res_v <- terra::crop(res_v, terra::ext(pet_r))
  #    
  #    daily_pet_w <- terra::extract(pet_r, res_v, fun = mean, na.rm = TRUE, ID = TRUE) %>%
  #      as_tibble() %>%
  #      pivot_longer(-ID, names_to = "name", values_to = "ET_mm") %>%
  #      mutate(
  #        # safer than separate(): works even if name formatting changes
  #        date = ymd(stringr::str_extract(name, "\\d{4}-\\d{2}-\\d{2}")),
  #        ET_mm = round(ET_mm, 2),
  #        GRAND_ID = res_ID_grandID
  #      ) %>%
  #      select(date, ET_mm, GRAND_ID)
  #  }
  
  
  # Precipitation data  bulk download (saved in  gridmet_pcp_data)
  
  
  #  start_date <- "1989-01-01"
  #  end_date   <- "2020-12-31"
  #  
  #  studied_reservoirs <- read_table('data/studied_reservoirs.txt')$studied_reservoirs
  #  reservoirs_polygons <- read_sf('C:/ASG/UCDavis/NSF_PROJECT_FOLDER/Datasets/GlobalDamWatcher/CONUS_reservoirs.shp')
  #  
  #  
  #  reservoirs_properties <- read_csv('data/Processed_data/reserv_properties.csv')
  #  
  #  
  #  
  #  reservoirs_pcp <- c()
  #  for(i in 1:length(studied_reservoirs)){
  #    
  #    
  #    it_n <- i
  #    
  #    res_ID_grandID <- studied_reservoirs[i]
  #    
  #    res_properties <- reservoirs_properties[reservoirs_properties$GRAND_ID == res_ID_grandID,]
  #    
  #    resname <- res_properties$Name
  #    
  #    # Reservoir properties and shapefile
  #    
  #    r_gdw <- reservoirs_properties[reservoirs_properties$GRAND_ID == res_ID_grandID,]
  #    r_pol <- reservoirs_polygons[reservoirs_polygons$GRAND_ID == res_ID_grandID, c('GDW_ID', 'RES_NAME', 'GRAND_ID')]
  #    
  #    
  #    
  #    res_pol <- r_pol %>% 
  #      st_make_valid() %>%
  #      filter(!st_is_empty(.)) %>%
  #      st_collection_extract("POLYGON") %>%
  #      st_cast("MULTIPOLYGON", warn = FALSE)
  #    
  #    pcp_res <- getGridMET(
  #      AOI = res_pol,
  #      varname = "pr",
  #      startDate = start_date,
  #      endDate  = end_date
  #    )
  #    
  #    # If getGridMET returns a tibble/data.frame (date, pet) ----
  #    if (inherits(pcp_res, "data.frame") && all(c("date", "precipitation_amount") %in% names(pcp_res))) {
  #      
  #      daily_pcp_w <- pcp_res %>%
  #        transmute(
  #          date = as.Date(date),
  #          PCP_mm = round(precipitation_amount, 2),
  #          GRAND_ID = grand_id
  #        )
  #      
  #    } else {
  #      # If getGridMET returns a list with SpatRaster ----
  #      pcp_r <- pcp_res$precipitation_amount
  #      
  #      if (is.null(pcp_r) || !inherits(pcp_r, "SpatRaster")) {
  #        stop("Unexpected getGridMET output for GRAND_ID = ", grand_id)
  #      }
  #      
  #      # match CRS and convert to SpatVector for terra
  #      res_pol2 <- st_transform(res_pol, st_crs(pcp_r))
  #      res_v    <- terra::vect(res_pol2)
  #      
  #      # crop to raster extent (helps when polygon barely overlaps / tiny rasters)
  #      res_v <- terra::crop(res_v, terra::ext(pcp_r))
  #      
  #      daily_pcp_w <- terra::extract(pcp_r, res_v, fun = mean, na.rm = TRUE, ID = TRUE) %>%
  #        as_tibble() %>%
  #        pivot_longer(-ID, names_to = "name", values_to = "PCP_mm") %>%
  #        mutate(
  #          # safer than separate(): works even if name formatting changes
  #          date = ymd(stringr::str_extract(name, "\\d{4}-\\d{2}-\\d{2}")),
  #          PCP_mm = round(PCP_mm, 2),
  #          GRAND_ID = res_ID_grandID
  #        ) %>%
  #        select(date, PCP_mm, GRAND_ID)
  #    }
  #    
  #    reservoirs_pcp <- rbind(reservoirs_pcp, daily_pcp_w)
  #    
  #  }
  #  
  #  reservoirs_pcp %>% 
  #    group_by(GRAND_ID) %>% 
  #    summarise(mean = mean(PCP_mm),
  #              totpcp = sum(PCP_mm)/32,
  #              min = min(PCP_mm),
  #              max = max(PCP_mm) ) %>%
  #    pivot_longer(mean:max) %>% 
  #    
  #    ggplot(., aes(x = name, y = value))+
  #    geom_boxplot()+
  #    facet_wrap(~name, scales ='free')
  #  
  #  
  #  write_csv(reservoirs_pcp, 
  #            'data/Processed_data/5_et_series/gridmet_pcp_download.csv')
  #  
  
  # GET Evaporation data in volume (m3d)

  daily_pet_w <- gridmet_et_data %>% 
    filter(GRAND_ID == res_ID_grandID) %>% 
    mutate(ET_mm = replace_na(ET_mm, 0))
  
  
  res_et <- res_area %>% 
    left_join(., daily_pet_w) %>% 
    mutate(ET_m3d = ET_mm*area_km*1e6/1e3) %>% 
    select(GRAND_ID, date, ET_m3d, ET_mm)
  
  
  # List with gridmet E (mm) and derived E volume (m3d)
  reservoirs_evaporation_series[[it_n]] <- res_et
  
  
  
  # GET Precipitation data in volume (m3d)
  
  daily_pcp_w <- gridmet_pcp_data %>% 
    filter(GRAND_ID == res_ID_grandID)%>% 
    mutate(PCP_mm  = replace_na(PCP_mm , 0))
  
  res_pcp <- res_area %>% 
    left_join(., daily_pcp_w) %>% 
    mutate(PCP_m3d = PCP_mm*area_km*1e6/1e3) %>% 
    select(GRAND_ID, date, PCP_m3d, PCP_mm)
  
  
  # List with gridmet E (mm) and derived E volume (m3d)
  reservoirs_precipitation_series[[it_n]] <- res_pcp
  
  
  # Final check of the water balance and calculation of flux to close mass balance
  
  stor_init <- recalc_res_data %>%
  arrange(date) %>%
  slice(1) %>%
  pull(storage)
  
  
  wbal_res <- recalc_res_data %>% 
    left_join(res_et, by = c("GRAND_ID", "date")) %>%
    left_join(res_pcp, by = c("GRAND_ID", "date")) %>%
    group_by(GRAND_ID) %>%
    arrange(date, .by_group = TRUE) %>%
    
    mutate(
      Inflow_Mm3  = inflow  * 86400 / 1e6,
      Outflow_Mm3 = outflow * 86400 / 1e6,
      
      PCP_corr_Mm3 = 0.6 * PCP_m3d / 1e6, # Raw PCP found generally overstimated for most reservoirs
      ET_Mm3       = ET_m3d / 1e6,
      
      # Daily flux
      DayVar_corr = Inflow_Mm3 - Outflow_Mm3 + PCP_corr_Mm3 - ET_Mm3,
      DayVar_io   = Inflow_Mm3 - Outflow_Mm3,
      
      # Cumulative
      Cum_corr = cumsum(replace_na(DayVar_corr, 0)),
      Cum_io   = cumsum(replace_na(DayVar_io, 0)),
      
      Storage_cum_corr = stor_init + Cum_corr,
      Storage_cum_io   = stor_init + Cum_io,
      
      # --- ERROR ---
      dS_obs = storage - lag(storage),
      extr   = dS_obs - DayVar_corr
    )
  
  storage_plot <- wbal_res %>% 
   #select(date, storage, stor_perc, Storage_cum, 
   #       Not_accounted, extr) %>% 
    ggplot(., aes(x = date))+
    geom_line(aes(y = Storage_cum_corr, color = 'Corrected (P and E)'), linewidth = .8) +
    geom_line(aes(y = Storage_cum_io, color = 'I-O only'), linewidth = .8) +
    geom_line(aes(y = storage, color = 'Observed'), linewidth = .8)+
    geom_hline(aes(yintercept = max_cap_res), color = 'black', linetype = 2)+
    scale_color_manual(values = c('blue4', 'grey', 'red4'))+
    scale_y_continuous(sec.axis = sec_axis(transform = ~100*./max_cap_res, name = 'Storage (%)'))+
    theme_bw()+
    labs(x = 'Date', y = 'Storage (Mm³)', color = '',
         title = paste0('WB', 
                        ' GRAND_ID = ', res_ID_grandID, ', ', 
                        resname, ' (', round(max_cap_res, 0), ' Mm³)  Uses: ', 
                        res_properties$Uses)) +
    theme(legend.position = c(0.15, 0.75),
          text = element_text(size = 18),
          axis.text.x = element_blank(),
          legend.background = element_blank(),
          title = element_text(size = 15))
  
  cumerror_plot <- wbal_res %>% 
    mutate(
      cum_extr = -cumsum(replace_na(extr, 0))
    ) %>% 
    ggplot(aes(x = date)) +
    geom_area(aes(y = cum_extr), fill = 'deepskyblue') +
    geom_hline(yintercept = 0, linetype = 2) +
    labs(
      x = 'Date', 
      y = 'Cumulative imbalance (Mm³)'
    ) +
    theme_bw() +
    theme(
      legend.position = c(0.15, 0.1),
      text = element_text(size = 15)
    )
  
  inout_plot <- wbal_res %>% 
    select(date, inflow, outflow) %>% 
    ggplot(., aes(x = date))+
    geom_line(aes(y = inflow, 
                  color = 'Inflow'))+
    geom_line(aes(y = outflow, 
                  color = 'Outflow'))+
    scale_color_manual(values = c('dodgerblue', 'brown4'))+
    labs(x = 'Date', y = 'Inflow and Outflow (m³/s)') +
    theme_bw()+
    theme(legend.position = c(0.15, 0.75),
          text = element_text(size = 15),
          axis.title.x = element_blank(),
          axis.text.x = element_blank())
  
  ylim_upper <- quantile(abs(wbal_res$extr), 0.975, na.rm = TRUE)
  
  extr_plot <- wbal_res %>% 
    ggplot(., aes(x = date))+
    geom_area(aes(y = extr),
              fill = 'red4')+
    geom_line(aes(y = ET_m3d/1e6), color = 'green4')+
    coord_cartesian(ylim = c(-ylim_upper, ylim_upper)) +
    labs(x = 'Date', y = 'Daily difference\n and Evaporation (Mm³)') +
    theme_bw()+
    theme(legend.position = c(0.15, 0.1),
          text = element_text(size = 15))
  
  
  complot <- storage_plot/ inout_plot/ cumerror_plot / extr_plot +
    plot_layout(heights = c(2, 3, 1, 3))
  
  
  ggsave(plot = complot, 
         filename = paste0('figures/time_series_figs/4_InOuStEv/', 
                           it_n, '_wbcheck_', res_ID_grandID, '.png'),
         device = 'png',
         width = 12, height = 14, dpi = 300)
  
  
  
  # Net inflow calculation and daily error calculation
  
  
  time_series_simulation <- wbal_res %>%
    # First extr always NA --> 0
    group_by(GRAND_ID) %>%
    mutate(
      extr = if_else(row_number() == 1, 0, extr)
    ) %>%
    ungroup() %>% 
    select(GRAND_ID, date, inflow, outflow, storage, ET_mm, PCP_mm, extr) %>%
    rename(daily_err_Mm3 = extr)
  
  time_series_simulations[[it_n]] <- time_series_simulation
  
  
  print(paste0('Reservoir ', i, ' of ', length(studied_reservoirs), ' done (', (100*i/length(studied_reservoirs)), ' %)'))
  
}




# Adding real coordinates, extracted from the GDW dams shapefile

reservoirs_coords <- read_sf('C:/ASG/UCDavis/NSF_PROJECT_FOLDER/Datasets/GlobalDamWatcher/CONUS_barriers.shp') %>% 
  select(GRAND_ID, LAT_DAM, LONG_DAM) %>% 
  st_transform(., crs = 4326) %>% 
  mutate(Lat = st_coordinates(.)[,2],
         Long = st_coordinates(.)[,1]) %>% 
  select(GRAND_ID, Lat, Long)


bind_rows(reservoirs_properties) %>% 
  # Duplicated and wrong area calculated (spatial error e.g., invalid geometry)
  filter(!(GRAND_ID == 1753 & area_km < 1)) %>% 
  mutate(Uses = if_else(is.na(Main_use), 'Irrigation', Uses),
         Main_use = if_else(is.na(Main_use), 'Irrigation', Main_use) # Vallecito Reservoir (575) purpose lacks 
  ) %>% 
  select(-c(Lat, Long)) %>% 
  st_drop_geometry(.) %>% 
  left_join(., st_drop_geometry(reservoirs_coords), by = 'GRAND_ID') %>% 
  write_csv(., 'data/Processed_data/reserv_properties.csv')


# Fixing initial storage for reservoir 530 (starts in 1990)
fix_530 <- time_series_simulations %>% bind_rows(.) %>% 
  filter(GRAND_ID == 530) %>% 
  select(GRAND_ID, date, inflow, outflow, storage, ET_mm, PCP_mm, daily_err_Mm3) %>% 
  slice(1) %>% 
  mutate(date = date-1,
         storage = 51.5, 
         daily_err_Mm3 = 0)


bind_rows(reservoirs_available_data) %>% 
  write_csv(., 'data/Processed_data/data_avail.csv')
bind_rows(hypsometrics_curves) %>% 
  write_csv(., 'data/Processed_data/hypso_curves.csv')
bind_rows(reservoirs_area_series) %>% 
  write_csv(., 'data/Processed_data/res_area_series.csv')
bind_rows(reservoirs_evaporation_series) %>% 
  write_csv(., 'data/Processed_data/res_et_series.csv')
bind_rows(reservoirs_precipitation_series) %>% 
  write_csv(., 'data/Processed_data/res_pcp_series.csv')
bind_rows(time_series_simulations) %>% 
  rbind(fix_530) %>% 
  arrange(., GRAND_ID, date) %>% 
  write_csv(., 'data/Processed_data/time_series_for_simulations.csv')
bind_rows(reservoir_lack_data_daily) %>% 
  write_csv(., 'data/Processed_data/reservoir_availability_timeseries.csv')

