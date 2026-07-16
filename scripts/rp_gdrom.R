# GDROM v2 release policy
#
# Script developed to implement GDROM model (v2, Zheng et al., 2025), originally
# developed by Chen et al., (2022)

# Fitting step is included in phyton scripts available in the repository
# (gdrom_01_prepare.py and gdrom_02_train.py)

#### 1. Input datasets ####

BASE <- normalizePath(getwd(), mustWork = FALSE)
DATA <- file.path(BASE, "data")

# input data 
studied_reservoirs <- read_table(file.path(DATA, "studied_reservoirs.txt"))$studied_reservoirs
reservoir_properties <- read_csv(file.path(DATA, "Processed_data/reserv_properties.csv"), show_col_types = FALSE)
reservoirs_time_series <- read_csv(file.path(DATA, "Processed_data/time_series_for_simulations.csv"), show_col_types = FALSE) %>%
  mutate(date = as.Date(date))
hypso_curves <- read_csv(file.path(DATA, "Processed_data/hypso_curves.csv"), show_col_types = FALSE)
pdsi_daily <- read_csv(file.path(DATA, "release_policies/gdrom/fit/PDSI/reservoir_pdsi_daily.csv"), show_col_types = FALSE) %>%
  mutate(date = as.Date(date))

FITTED_DIR <- "C:/ASG/UCDavis/Research/Papers/RPB/data/Release_policies/gdrom/data/release_policies/gdrom/fit/gdrom_fitted"
MOD_FITTED <- file.path(FITTED_DIR, "gdrom_module_params.csv")
CART_FITTED <- file.path(FITTED_DIR, "gdrom_cart_params.csv")
META_FITTED <- file.path(FITTED_DIR, "gdrom_metadata.csv")

# Default GDROM rules from the original paper release.
DEFAULT_RULES <- file.path(
  "C:/ASG/UCDavis/NSF_PROJECT_FOLDER/Datasets/GDROM_v2/data/contents/Operation Rules - GDROMs"
)

# Output folders.
OUT_BASE <- file.path(DATA, "release_policies/gdrom")
OUT_SIMS <- file.path(OUT_BASE, "simulations")
OUT_PERF <- file.path(OUT_BASE, "performance")
OUT_PARAMS <- file.path(OUT_BASE, "parameters")
OUT_CV <- file.path(OUT_BASE, "convergence")
FIG_BASE <- file.path(BASE, "figures/release_policies/gdrom")
FIG_DEFAULT <- file.path(FIG_BASE, "default/sims")
FIG_FITTED <- file.path(FIG_BASE, "fitted/sims")
FIG_PARAMS <- file.path(FIG_BASE, "fitted/params")

for (d in c(OUT_SIMS, OUT_PERF, OUT_PARAMS, OUT_CV, FIG_DEFAULT, FIG_FITTED, FIG_PARAMS)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Unit conversions.
M3S_TO_MCM <- 86400 / 1e6       # m3/s -> Mm3/day
MCM_TO_ACFT <- 810.714          # Mm3 -> acre-ft
ACFT_TO_MCM <- 1 / MCM_TO_ACFT  # acre-ft -> Mm3
M3S_TO_ACFT <- M3S_TO_MCM * MCM_TO_ACFT
ACFT_TO_M3S <- 1 / M3S_TO_ACFT

# Simulation periods used in the project.
CALIB_START <- as.Date("1990-01-01")
CALIB_END   <- as.Date("2004-12-31")
VAL_START   <- as.Date("2005-01-01")
VAL_END     <- as.Date("2015-12-31")
FULL_START  <- as.Date("1990-01-01")
FULL_END    <- as.Date("2015-12-31")

MODULE_COLS <- c(
  "#2166AC", "#D6604D", "#4DAC26", "#8073AC",
  "#F4A582", "#92C5DE", "#1A9850", "#762A83"
)

#### 2. Functions ####

##### 2.1. Common helpers ####

read_required <- function(path, ...) {
  if (!file.exists(path)) stop(sprintf("Missing required file: %s", path))
  readr::read_csv(path, show_col_types = FALSE, ...)
}

read_required_table <- function(path) {
  if (!file.exists(path)) stop(sprintf("Missing required file: %s", path))
  readr::read_table(path, show_col_types = FALSE)
}

safe_first <- function(x, default = NA_real_) {
  x <- x[!is.na(x)]
  if (length(x) == 0) default else x[[1]]
}

# Returns the first matching column among several options.
find_col <- function(df, candidates, required = TRUE, label = "column") {
  hit <- candidates[candidates %in% names(df)][1]
  if (is.na(hit) || length(hit) == 0) {
    if (required) stop(sprintf("Could not find %s. Available: %s", label, paste(names(df), collapse = ", ")))
    return(NA_character_)
  }
  hit
}

m3s_to_mcm_day <- function(x) x * M3S_TO_MCM

##### 2.2. Parameter loading ####

# Reads GDROM rules for one reservoir.
# fitted:
#   - direct CSV files written by Python
# default:
#   - original text rules from the GDROM paper release
load_params <- function(grand_id, source = c("fitted", "default")) {
  source <- match.arg(source)
  
  if (source == "fitted") {
    mod <- mod_fit %>% filter(GRAND_ID == grand_id)
    cart <- cart_fit %>% filter(GRAND_ID == grand_id)
    if (nrow(mod) == 0) return(NULL)
    return(list(mod = mod, cart = cart, smax_acft = mod$smax_acft[1]))
  }
  
  mod_dir <- file.path(DEFAULT_RULES, "modules")
  cond_dir <- file.path(DEFAULT_RULES, "module_conditions")
  mod_files <- list.files(mod_dir, pattern = paste0("^", grand_id, "_.*\\.txt$"), full.names = TRUE)
  if (length(mod_files) == 0) return(NULL)
  
  cap <- reservoir_properties %>% filter(GRAND_ID == grand_id) %>% pull(Cap_mcm)
  if (length(cap) == 0 || is.na(cap[1])) return(NULL)
  smax_acft <- cap[1] * MCM_TO_ACFT
  
  # Parse a module file into one row per leaf rule.
  parse_mod_file <- function(txt, mod_idx) {
    lines <- readLines(txt, warn = FALSE)
    rows <- list()
    
    # Single-leaf tree exported as one constant release line.
    if (length(lines) == 1 && grepl("^Release:", lines[1])) {
      return(tibble(
        GRAND_ID = grand_id, module_id = mod_idx, node_id = 0,
        smax_acft = smax_acft,
        Inflow_lo = -1e15, Inflow_hi = 1e15,
        Storage_lo = -1e15, Storage_hi = 1e15,
        release_acft = as.numeric(sub("Release:\\s*", "", lines[1]))
      ))
    }
    
    for (ln in lines) {
      if (!grepl("^if", ln)) next
      conds <- str_extract_all(ln, "\\([^)]+\\)")[[1]]
      
      # Some text files contain malformed leaves like "if  then Release: ...".
      # Keep them as broad rules instead of failing.
      if (length(conds) == 0) {
        m <- regmatches(ln, regexpr("Release: [0-9.\\-]+", ln))
        if (length(m) > 0) {
          rows[[length(rows) + 1]] <- tibble(
            GRAND_ID = grand_id, module_id = mod_idx, node_id = length(rows),
            smax_acft = smax_acft,
            Inflow_lo = -1e15, Inflow_hi = 1e15,
            Storage_lo = -1e15, Storage_hi = 1e15,
            release_acft = as.numeric(sub("Release: ", "", m))
          )
        }
        next
      }
      
      il <- -1e15; ih <- 1e15; sl <- -1e15; sh <- 1e15
      for (cond in conds) {
        p <- str_split(trimws(gsub("[()]", "", cond)), "\\s+")[[1]]
        if (length(p) < 3) next
        v <- suppressWarnings(as.numeric(p[3]))
        if (is.na(v)) next
        if      (p[1] == "Inflow"  & p[2] == "<=") ih <- min(ih, v)
        else if (p[1] == "Inflow"  & p[2] == ">" ) il <- max(il, v)
        else if (p[1] == "Storage" & p[2] == "<=") sh <- min(sh, v)
        else if (p[1] == "Storage" & p[2] == ">" ) sl <- max(sl, v)
      }
      m <- regmatches(ln, regexpr("Release: [0-9.\\-]+", ln))
      if (length(m) == 0) next
      rows[[length(rows) + 1]] <- tibble(
        GRAND_ID = grand_id, module_id = mod_idx, node_id = length(rows),
        smax_acft = smax_acft,
        Inflow_lo = il, Inflow_hi = ih, Storage_lo = sl, Storage_hi = sh,
        release_acft = as.numeric(sub("Release: ", "", m))
      )
    }
    
    bind_rows(rows)
  }
  
  # Parse the CART file. CART chooses the active module from inflow, storage,
  # PDSI and day-of-year.
  parse_cart_file <- function(txt) {
    lines <- readLines(txt, warn = FALSE)
    rows <- list()
    for (ln in lines) {
      if (!grepl("^if", ln)) next
      conds <- str_extract_all(ln, "\\([^)]+\\)")[[1]]
      il <- -1e15; ih <- 1e15; sl <- -1e15; sh <- 1e15
      pl <- -1e15; ph <- 1e15; dl <- -1e15; dh <- 1e15
      for (cond in conds) {
        p <- str_split(trimws(gsub("[()]", "", cond)), "\\s+")[[1]]
        if (length(p) < 3) next
        v <- suppressWarnings(as.numeric(p[3]))
        if (is.na(v)) next
        if      (p[1] == "Inflow"  & p[2] == "<=") ih <- min(ih, v)
        else if (p[1] == "Inflow"  & p[2] == ">" ) il <- max(il, v)
        else if (p[1] == "Storage" & p[2] == "<=") sh <- min(sh, v)
        else if (p[1] == "Storage" & p[2] == ">" ) sl <- max(sl, v)
        else if (p[1] == "PDSI"    & p[2] == "<=") ph <- min(ph, v)
        else if (p[1] == "PDSI"    & p[2] == ">" ) pl <- max(pl, v)
        else if (p[1] == "DOY"     & p[2] == "<=") dh <- min(dh, v)
        else if (p[1] == "DOY"     & p[2] == ">" ) dl <- max(dl, v)
      }
      m <- regmatches(ln, regexpr("module: [0-9]+", ln))
      if (length(m) == 0) next
      rows[[length(rows) + 1]] <- tibble(
        GRAND_ID = grand_id, node_id = length(rows), smax_acft = smax_acft,
        Inflow_lo = il, Inflow_hi = ih, Storage_lo = sl, Storage_hi = sh,
        PDSI_lo = pl, PDSI_hi = ph, DOY_lo = dl, DOY_hi = dh,
        module_id = as.integer(sub("module: ", "", m))
      )
    }
    bind_rows(rows)
  }
  
  mod <- map_dfr(seq_along(mod_files), function(i) parse_mod_file(mod_files[i], i - 1))
  if (is.null(mod) || nrow(mod) == 0) return(NULL)
  
  cart_file <- file.path(cond_dir, paste0(grand_id, ".txt"))
  cart <- if (file.exists(cart_file)) parse_cart_file(cart_file) else NULL
  if (is.null(cart) || nrow(cart) == 0) {
    cart <- tibble(
      GRAND_ID = grand_id, node_id = 0L, smax_acft = smax_acft,
      Inflow_lo = -1e15, Inflow_hi = 1e15,
      Storage_lo = -1e15, Storage_hi = 1e15,
      PDSI_lo = -1e15, PDSI_hi = 1e15,
      DOY_lo = -1e15, DOY_hi = 1e15,
      module_id = 0L
    )
  }
  
  list(mod = mod, cart = cart, smax_acft = smax_acft)
}

##### 2.3. Hypsometry and rule lookup ####

# Convert storage to surface area.
# Needed to compute PCP and ET volumes from depths.
get_area <- function(h, storage) {
  area <- if (h$poly2 == 0) {
    h$intercept + h$poly1 * storage
  } else if (!is.na(h$max_stor_rel) && storage > h$max_stor_rel) {
    h$intercept + h$poly1 * h$max_stor_rel + h$poly2 * h$max_stor_rel^2
  } else {
    h$intercept + h$poly1 * storage + h$poly2 * storage^2
  }
  max(area, 0.0005)
}

# Given the current state, find the GDROM release rule.
# Inputs are converted to acre-ft because the trained thresholds are in acre-ft.
get_release_acft <- function(inf_acft, sto_acft, pdsi, doy, mod, cart) {
  mod_id <- 0L
  if (!is.null(cart) && nrow(cart) > 0) {
    for (i in seq_len(nrow(cart))) {
      r <- cart[i, ]
      if (inf_acft >= r$Inflow_lo && inf_acft < r$Inflow_hi &&
          sto_acft >= r$Storage_lo && sto_acft < r$Storage_hi &&
          pdsi >= r$PDSI_lo && pdsi < r$PDSI_hi &&
          doy >= r$DOY_lo && doy < r$DOY_hi) {
        mod_id <- as.integer(r$module_id)
        break
      }
    }
  }
  
  rows <- mod[mod$module_id == mod_id, ]
  if (nrow(rows) == 0) rows <- mod
  
  for (i in seq_len(nrow(rows))) {
    r <- rows[i, ]
    if (inf_acft >= r$Inflow_lo && inf_acft < r$Inflow_hi &&
        sto_acft >= r$Storage_lo && sto_acft < r$Storage_hi) {
      return(list(release_acft = as.numeric(r$release_acft), module_id = mod_id))
    }
  }
  
  list(release_acft = as.numeric(rows$release_acft[1]), module_id = mod_id)
}

# Apply the spill and nonnegative-storage constraints.
physical_constraints <- function(storage, max_capacity, provisional_release) {
  spill <- max(0, storage - max_capacity)
  final_release <- provisional_release + spill
  final_storage <- max(0, storage - spill)
  tibble(final_release = final_release, final_storage = final_storage)
}

##### 2.4. Performance metrics ####

compute_perf <- function(df, grand_id, period_label) {
  r2 <- function(s, o) {
    ok <- is.finite(s) & is.finite(o)
    if (sum(ok) < 5) return(NA_real_)
    hydroGOF::rPearson(s[ok], o[ok])^2
  }
  
  df %>%
    summarise(
      NSE_out   = hydroGOF::NSE(sim_outflow_Mm3,  obs_outflow_Mm3),
      KGE_out   = hydroGOF::KGE(sim_outflow_Mm3,  obs_outflow_Mm3),
      PBIAS_out = hydroGOF::pbias(sim_outflow_Mm3, obs_outflow_Mm3),
      RMSE_out  = hydroGOF::rmse(sim_outflow_Mm3,  obs_outflow_Mm3),
      R2_out    = r2(sim_outflow_Mm3, obs_outflow_Mm3),
      NSE_sto   = hydroGOF::NSE(sim_storage_Mm3,  obs_storage_Mm3),
      KGE_sto   = hydroGOF::KGE(sim_storage_Mm3,  obs_storage_Mm3),
      PBIAS_sto = hydroGOF::pbias(sim_storage_Mm3, obs_storage_Mm3),
      RMSE_sto  = hydroGOF::rmse(sim_storage_Mm3,  obs_storage_Mm3),
      R2_sto    = r2(sim_storage_Mm3, obs_storage_Mm3)
    ) %>%
    mutate(GRAND_ID = grand_id, period = period_label, .before = 1)
}

##### 2.5. GDROM simulation ####

# The simulation follows the same structure as the other reservoir policies:
#   1) add observed inflow to previous storage
#   2) choose release using the GDROM rules
#   3) subtract release
#   4) compute area using hypsometry
#   5) add PCP, subtract ET, then add the balance correction
#   6) apply spill/capacity constraints
simulate_gdrom <- function(grand_id, params, ts_sub, pdsi_sub, h,
                           cap_mcm, init_stor,
                           date_start, date_end,
                           pcp_fix = 0.6,
                           err_fix = TRUE) {
  mod <- params$mod
  cart <- params$cart
  
  ts_p <- ts_sub %>%
    filter(date >= date_start, date <= date_end) %>%
    left_join(pdsi_sub %>% select(date, pdsi), by = "date") %>%
    mutate(
      inflow_Mm3 = replace_na(m3s_to_mcm_day(inflow), 0),
      outflow_Mm3 = replace_na(m3s_to_mcm_day(outflow), 0),
      et_mm = replace_na(ET_mm, 0),
      pcp_mm = replace_na(PCP_mm, 0),
      err_Mm3 = replace_na(daily_err_Mm3, 0),
      pdsi = replace_na(pdsi, 0),
      doy = yday(date)
    )
  
  if (nrow(ts_p) == 0) return(NULL)
  
  n <- nrow(ts_p)
  sim_out <- numeric(n)
  sim_sto <- numeric(n)
  sim_area <- numeric(n)
  sim_pcp <- numeric(n)
  sim_et <- numeric(n)
  sim_rule_rel <- numeric(n)
  stor <- init_stor
  
  for (t in seq_len(n)) {
    rw <- ts_p[t, ]
    
    # 1) GDROM rule lookup uses observed inflow and previous-day storage.
    lookup <- get_release_acft(
      inf_acft = rw$inflow_Mm3 * MCM_TO_ACFT,
      sto_acft = stor * MCM_TO_ACFT,
      pdsi = rw$pdsi,
      doy = rw$doy,
      mod = mod,
      cart = cart
    )
    rel_Mm3 <- max(lookup$release_acft * ACFT_TO_MCM, 0)
    sim_rule_rel[t] <- rel_Mm3
    
    # 2) Add inflow.
    init_s <- stor + rw$inflow_Mm3
    
    # 3) Apply release, limited by available water.
    rel_Mm3 <- min(rel_Mm3, init_s)
    upd_s <- init_s - rel_Mm3
    
    # 4) Calculate area from post-release storage.
    area <- get_area(h, upd_s)
    
    # 5) Convert PCP and ET depths into volumes using the current surface area.
    pcp_M <- (pcp_fix * area * rw$pcp_mm) / 1e3
    et_M  <- (area * rw$et_mm) / 1e3
    
    # Balance correction is added as in natural_lake.
    mb_flux <- pcp_M - et_M + if (err_fix) rw$err_Mm3 else 0
    final_s <- upd_s + mb_flux
    
    # 6) Spill and capacity constraint.
    spill <- max(0, final_s - cap_mcm)
    final_rel <- rel_Mm3 + spill
    final_s <- max(min(final_s, cap_mcm), 0)
    
    sim_out[t] <- final_rel
    sim_sto[t] <- final_s
    sim_area[t] <- area
    sim_pcp[t] <- pcp_M
    sim_et[t] <- et_M
    stor <- final_s
  }
  
  ts_p %>%
    mutate(
      GRAND_ID = grand_id,
      sim_outflow_Mm3 = sim_out,
      sim_storage_Mm3 = sim_sto,
      sim_area_km2 = sim_area,
      sim_pcp_Mm3 = sim_pcp,
      sim_et_Mm3 = sim_et,
      gdrom_release_rule_Mm3 = sim_rule_rel,
      obs_outflow_Mm3 = outflow_Mm3,
      obs_storage_Mm3 = storage
    ) %>%
    select(GRAND_ID, date, doy,
           obs_outflow_Mm3, obs_storage_Mm3,
           sim_outflow_Mm3, sim_storage_Mm3,
           gdrom_release_rule_Mm3,
           sim_area_km2, sim_pcp_Mm3, sim_et_Mm3,
           inflow_Mm3, outflow_Mm3, pdsi, err_Mm3)
}

##### 2.6. Convergence table ####

convergence_table <- function(sim_df, grand_id, period_label, cap_mcm) {
  if (is.null(sim_df) || nrow(sim_df) == 0) return(NULL)
  
  sim_df %>%
    summarise(
      GRAND_ID = grand_id,
      period = period_label,
      n_days = n(),
      cap_mcm = cap_mcm,
      start_storage_Mm3 = first(sim_storage_Mm3),
      end_storage_Mm3   = last(sim_storage_Mm3),
      delta_storage_Mm3 = end_storage_Mm3 - start_storage_Mm3,
      inflow_Mm3  = sum(inflow_Mm3, na.rm = TRUE),
      outflow_Mm3 = sum(sim_outflow_Mm3, na.rm = TRUE),
      pcp_Mm3     = sum(sim_pcp_Mm3, na.rm = TRUE),
      et_Mm3      = sum(sim_et_Mm3, na.rm = TRUE),
      err_Mm3     = sum(err_Mm3, na.rm = TRUE),
      balance_Mm3 = inflow_Mm3 - outflow_Mm3 + pcp_Mm3 - et_Mm3 + err_Mm3,
      closure_Mm3 = delta_storage_Mm3 - balance_Mm3,
      closure_pct_cap = 100 * closure_Mm3 / cap_mcm
    )
}

##### 2.7. Hydrograph plot ####

# Observed vs simulated outflow and storage, with a small performance table.
make_hydro_plot <- function(sim_list, grand_id, res_name, cap_mcm,
                            show_vline = TRUE, extra_title = "") {
  obsims <- bind_rows(sim_list) %>%
    pivot_longer(-c(date, period, GRAND_ID, gdrom_release_rule_Mm3, sim_area_km2, sim_pcp_Mm3, sim_et_Mm3, inflow_Mm3, outflow_Mm3, pdsi, err_Mm3)) %>%
    filter(name %in% c("obs_outflow_Mm3", "sim_outflow_Mm3", "obs_storage_Mm3", "sim_storage_Mm3")) %>%
    mutate(
      type = if_else(str_starts(name, "obs"), "obs", "sim"),
      var = if_else(str_detect(name, "outflow"), "Outflow (Mm3/day)", "Storage (Mm3)")
    )
  
  p <- ggplot(obsims, aes(x = date, y = value, color = type)) +
    geom_line(linewidth = 0.6) +
    { if (show_vline) geom_vline(xintercept = CALIB_END, linetype = "dashed", linewidth = 1, color = "blue4") else list() } +
    facet_wrap(~var, scales = "free", ncol = 1) +
    scale_color_manual(values = c("grey30", "red")) +
    labs(title = paste0(res_name, " (", round(cap_mcm), " Mm3)", extra_title), x = "Date") +
    theme_bw() +
    theme(
      legend.position = c(0.9, 0.9),
      legend.background = element_blank(),
      legend.title = element_blank(),
      axis.title = element_blank(),
      text = element_text(size = 15),
      title = element_text(size = 9)
    )
  
  perf_all <- bind_rows(lapply(sim_list, function(d) compute_perf(d, grand_id, unique(d$period))))
  
  tab <- perf_all %>%
    select(Period = period, GRAND_ID, KGE_out, PBIAS_out, R2_out, KGE_sto, PBIAS_sto, R2_sto) %>%
    mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
    pivot_longer(-c(Period, GRAND_ID)) %>%
    mutate(
      Metric = str_extract(name, "^[^_]+"),
      Var = if_else(str_detect(name, "out"), "Outflow", "Storage")
    ) %>%
    select(Var, Metric, Value = value, Period) %>%
    pivot_wider(names_from = Metric, values_from = Value) %>%
    arrange(Var) %>%
    rename(`R²` = R2)
  
  tab_plot <- ggplot() +
    annotate(geom = "table", x = 1, y = 1, label = list(tab), size = 3.2) +
    theme_bw() +
    theme(
      axis.line = element_blank(), line = element_blank(), text = element_blank(),
      panel.border = element_blank(), panel.background = element_rect(fill = "transparent", color = NA),
      plot.background = element_rect(fill = "transparent", color = NA)
    )
  
  p + patchwork::inset_element(tab_plot, left = 0.75, bottom = 0.9, right = 0.99, top = 0.99, align_to = "full")
}

#### 3. One reservoir test ####

# This section checks one reservoir first, comparing default and fitted GDROM.
cat("Running one-reservoir test...\n")

gid_test <- 132
obs_test <- reservoirs_time_series %>% filter(GRAND_ID == gid_test)
pr_test  <- reservoir_properties %>% filter(GRAND_ID == gid_test)
h_test   <- hypso_curves %>% filter(GRAND_ID == gid_test)
pd_test  <- pdsi_daily %>% filter(GRAND_ID == gid_test)

# Default GDROM rules.
def_params <- load_params(gid_test, "default")
if (!is.null(def_params)) {
  def_sim <- simulate_gdrom(
    grand_id = gid_test,
    params = def_params,
    ts_sub = reservoirs_time_series,
    pdsi_sub = pdsi_daily,
    h = h_test,
    cap_mcm = pr_test$Cap_mcm[1],
    init_stor = obs_test %>% filter(date >= FULL_START) %>% slice(1) %>% pull(storage),
    date_start = FULL_START,
    date_end = FULL_END,
    pcp_fix = 0.6,
    err_fix = TRUE
  )
  def_obs <- def_sim %>%
    select(date, sim_storage_Mm3, sim_outflow_Mm3) %>%
    left_join(obs_test %>% transmute(date, obs_outflow_Mm3 = outflow * M3S_TO_MCM, obs_storage_Mm3 = storage), by = "date")
  cat("Default test metrics:\n")
  print(def_obs %>% summarise(KGE_sto = KGE(sim_storage_Mm3, obs_storage_Mm3), KGE_out = KGE(sim_outflow_Mm3, obs_outflow_Mm3)))
}

# Fitted GDROM rules.
fit_params <- load_params(gid_test, "fitted")
if (!is.null(fit_params)) {
  fit_cal <- simulate_gdrom(
    grand_id = gid_test,
    params = fit_params,
    ts_sub = reservoirs_time_series,
    pdsi_sub = pdsi_daily,
    h = h_test,
    cap_mcm = pr_test$Cap_mcm[1],
    init_stor = obs_test %>% filter(date >= CALIB_START) %>% slice(1) %>% pull(storage),
    date_start = CALIB_START,
    date_end = CALIB_END,
    pcp_fix = 0.6,
    err_fix = TRUE
  )
  fit_val <- simulate_gdrom(
    grand_id = gid_test,
    params = fit_params,
    ts_sub = reservoirs_time_series,
    pdsi_sub = pdsi_daily,
    h = h_test,
    cap_mcm = pr_test$Cap_mcm[1],
    init_stor = obs_test %>% filter(date >= VAL_START) %>% slice(1) %>% pull(storage),
    date_start = VAL_START,
    date_end = VAL_END,
    pcp_fix = 0.6,
    err_fix = TRUE
  )
  fit_obs <- bind_rows(fit_cal, fit_val) %>%
    select(date, period, sim_storage_Mm3, sim_outflow_Mm3) %>%
    left_join(obs_test %>% transmute(date, obs_outflow_Mm3 = outflow * M3S_TO_MCM, obs_storage_Mm3 = storage), by = "date")
  cat("Fitted test metrics:\n")
  print(fit_obs %>% summarise(KGE_sto = KGE(sim_storage_Mm3, obs_storage_Mm3), KGE_out = KGE(sim_outflow_Mm3, obs_outflow_Mm3)))
}
