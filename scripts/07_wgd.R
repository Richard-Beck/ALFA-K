## Whole‑genome‑doubling (WGD) analysis – mixed‑model edition 2025‑04‑28 (rev‑E)
if(!basename(getwd())=="ALFA-K") stop("Ensure working directory set to ALFA-K root")

source("R/utils_karyo.R")
source("R/utils_theme.R")
source("R/utils_env.R")
ensure_packages(c("ggplot2","pbapply","cowplot","lme4","lmerTest"))
base_text_size <- 5
common_theme <- make_base_theme("classic",base_text_size)

# ───────────────────────────────────────── 1 • helpers
is_wgd   <- function(k, cut = 3) mean(k) > cut
get_mode <- function(x) { ux <- unique(x); ux[which.max(tabulate(match(x, ux)))] }

# ───────────────────────────────────────── 2 • load & QC fits
procDir <- "data/processed/salehi/alfak_outputs_proc/"
files   <- list.files(procDir, full.names = TRUE)

x_list <- lapply(files, readRDS)
x_list <- lapply(x_list, head, 1)

fits_raw <- do.call(rbind, x_list)
fits_raw <- subset(fits_raw, !is.na(xv) & xv > 0 & grepl("SA906", fi))
fits_raw$traj <- ifelse(grepl("906a", fits_raw$fi, TRUE), "trajA", "trajB")

ord <- with(fits_raw, order(fi, min_obs, -ntrain, -xv))
fits_sorted <- fits_raw[ord, ]
fits <- fits_sorted[!duplicated(fits_sorted[c("fi", "min_obs")]), ]
rownames(fits) <- NULL
cat("Loaded", nrow(fits), "unique fits after QC\n")

# 3A •  PASSAGE‑LEVEL dk WITHOUT FIT POOLING  (PlotA)
#       We take the two raw count matrices directly, ignoring any landscape QC.
## done because these two files contain all the raw data
raw_ids <- c("SA906a_X57_l_7_d1_0_d2_0", "SA906b_X55_l_8_d1_0_d2_0")

dfdiff <- do.call(rbind, lapply(raw_ids, function(fi){
  xin <- readRDS(file.path("data/processed/salehi/alfak_inputs", paste0(fi, ".Rds")))$x
  traj_lbl <- ifelse(grepl("906a", fi, TRUE), "trajA", "trajB")
  
  k  <- do.call(rbind, lapply(rownames(xin), s2v))
  dk <- apply(k, 1, function(ki) sum(ki != get_mode(ki)))
  pk <- sapply(rownames(xin), function(j) get_mode(s2v(j)))
  wgd_flag <- pk >= 3
  
  xin_freq <- sweep(xin, 2, colSums(xin), "/")
  pop_mat  <- xin
  
  do.call(rbind, lapply(c(FALSE, TRUE), function(flag){
    rows <- which(wgd_flag == flag); if(!length(rows)) return(NULL)
    xin_freq <- sweep(xin[rows,, drop=FALSE], 2, colSums(xin[rows,, drop=FALSE]), "/")
    data.frame(passage  = as.integer(colnames(xin)),
               avg_diff = colSums(xin_freq * dk[rows]),
               pop      = colSums(pop_mat[rows, , drop=FALSE]),
               wgd_str  = ifelse(flag, "WGD+", "WGD-"),
               traj     = traj_lbl,
               id       = fi)
  }))
}))

dfdiff$plot_id <- ifelse(dfdiff$traj == "trajA", "p53 k.o. A", "p53 k.o. B")


# ───────────────────────────────────────── 3 • build tables

dfl_list   <- list()
dffit_list <- list()

for (i in seq_len(nrow(fits))) {
  print(i)
  fi       <- fits$fi[i]
  min_obs  <- fits$min_obs[i]
  traj_lbl <- fits$traj[i]
  
  ## 3a – passage‑level dk & population
  xin_path <- file.path("data/processed/salehi/alfak_inputs", paste0(fi, ".Rds"))
  if (!file.exists(xin_path)) next
  xin <- readRDS(xin_path)$x
  
  k  <- do.call(rbind, lapply(rownames(xin), s2v))
  dk <- apply(k, 1, function(ki) sum(ki != get_mode(ki))) ## number of modifications
  pk <- sapply(rownames(xin), function(j) get_mode(s2v(j)))
  wgd_flag <- pk >= 3
  
  pop_mat  <- xin
  
  ## 3b – landscape predictions
  land_path <- file.path("data/processed/salehi/alfak_outputs", fi, paste0("minobs_", min_obs), "landscape.Rds")
  if (!file.exists(land_path)) next
  land <- readRDS(land_path)
  
  k_land  <- do.call(rbind, lapply(land$k, s2v))
  dk_land <- apply(k_land, 1, function(ki) sum(ki != get_mode(ki)))
  lut     <- land$mean; names(lut) <- land$k
  
  keep <- land$fq | land$nn
  dffit_tmp <- data.frame(
    fitness = lut[land$k[keep]],
    dk      = dk_land[keep],
    wgd     = sapply(land$k[keep], function(j) is_wgd(s2v(j))),
    id      = fi,
    k_id=land$k[keep],
    traj    = traj_lbl
  )
  dffit_list[[length(dffit_list) + 1]] <- dffit_tmp
  
  ## 3c – neighbour Δf
  fq <- land$k[land$fq]
  dfl_tmp <- do.call(rbind, lapply(fq, function(fqi) {
    nn <- apply(gen_all_neighbours(fqi, TRUE), 1, paste, collapse = ".")
    data.frame(df = lut[fqi] - lut[nn], wgd = is_wgd(s2v(fqi)), id = fi, traj = traj_lbl,k_id=fqi)
  }))
  dfl_list[[length(dfl_list) + 1]] <- dfl_tmp
}

dfl    <- do.call(rbind, dfl_list)
dffit  <- do.call(rbind, dffit_list)

# add WGD flag strings
dfl$wgd_str   <- ifelse(dfl$wgd,  "WGD+", "WGD-")
dffit$wgd_str <- ifelse(dffit$wgd, "WGD+", "WGD-")

# ───────────────────────────────────────── 4 • statistics

## 4a – permutation KS (fit‑level collapse)
set.seed(1)
cl_collapsed <- aggregate(df ~ id + wgd_str, dfl, mean)
obs_stat <- ks.test(split(cl_collapsed$df, cl_collapsed$wgd_str)[[1]],
                    split(cl_collapsed$df, cl_collapsed$wgd_str)[[2]])$statistic
perm_p <- mean(replicate(1e4, {
  lab <- sample(cl_collapsed$wgd_str)
  ks.test(split(cl_collapsed$df, lab)[[1]], split(cl_collapsed$df, lab)[[2]])$statistic
}) >= obs_stat)

library(permute)

if(FALSE){
  ##takes longer and p_hier still zero
  set.seed(1)
  
  ids  <- unique(dfl$id)                # one row per fitted landscape
  n_id <- length(ids)
  
  perm_D <- replicate(1e4, {
    ## 1 · assign a NEW random class to each landscape
    lab_map <- setNames(sample(c(TRUE, FALSE), n_id, replace = TRUE), ids)
    ## 2 · propagate to every Δf row
    perm_lab <- lab_map[dfl$id]
    
    ## 3 · make sure both classes are present; if not, resample once
    if (length(unique(perm_lab)) < 2) {
      lab_map <- setNames(sample(c(TRUE, FALSE), n_id, replace = TRUE), ids)
      perm_lab <- lab_map[dfl$id]
    }
    
    ks.test(dfl$df[perm_lab], dfl$df[!perm_lab])$statistic
  })
  
  obs_D  <- ks.test(dfl$df[dfl$wgd], dfl$df[!dfl$wgd])$statistic
  p_hier <- mean(perm_D >= obs_D)
}



## 4b – landscape shape (piecewise mixed model)
dffit$cut3 <- dffit$dk > 3
mod_piece <- lmer(fitness ~ dk * wgd_str + dk:cut3 * wgd_str + (1|id), data = dffit)

# --- Fit Additive Model with nls() & Generate Predictions ---

# 1. Prepare 0/1 indicator variables in your data frame 'dfdiff'
#    (Assumes 'wgd_str' and 'traj' are factors; uses 2nd level as contrast)
other_wgd_level <- levels(factor(dfdiff$wgd_str))[2]
other_traj_level <- levels(factor(dfdiff$traj))[2]
dfdiff$is_wgd_plus <- ifelse(dfdiff$wgd_str == other_wgd_level, 1, 0)
dfdiff$is_trajB <- ifelse(dfdiff$traj == other_traj_level, 1, 0)

# 2. Define starting values for the 6 parameters
#    (A0=Asym[ref], A_wgd=WGD+ effect, A_traj=TrajB effect, k0=k[ref], k_wgd=WGD+ effect, k_traj=TrajB effect)
start_params <- list(A0 = 1.2, A_wgd = 5.3, A_traj = 1.9, k0 = 0.03, k_wgd = 0.01, k_traj = -0.01)

# 3. Fit the nls model (basic error check)
fit_nls_additive <- tryCatch(nls(
  formula = avg_diff ~ (A0 + A_wgd * is_wgd_plus + A_traj * is_trajB) *
    (1 - exp(-(k0 + k_wgd * is_wgd_plus + k_traj * is_trajB) * passage)),
  data = dfdiff, start = start_params
), error = function(e) { message("NLS fitting failed: ", e$message); NULL })

# 4. Generate predictions (if fit worked)
predictions <- NULL
if (!is.null(fit_nls_additive)) {
  pred_grid <- expand.grid( # Create grid for all factor combinations over passage range
    passage = seq(min(dfdiff$passage), max(dfdiff$passage), length.out = 50), # Fewer points ok for lines
    wgd_str = levels(factor(dfdiff$wgd_str)), traj = levels(factor(dfdiff$traj)) )
  # Add indicators needed for prediction based on grid factor levels
  pred_grid$is_wgd_plus <- ifelse(pred_grid$wgd_str == other_wgd_level, 1, 0)
  pred_grid$is_trajB <- ifelse(pred_grid$traj == other_traj_level, 1, 0)
  # Predict and select relevant columns
  pred_grid$pred_avg_diff <- predict(fit_nls_additive, newdata = pred_grid)
  predictions <- pred_grid[, c("passage", "wgd_str", "traj", "pred_avg_diff")]
  # print(summary(fit_nls_additive)) # Optional: view model summary
  # print(head(predictions))          # Optional: view predictions head
}
predictions$plot_id <- ifelse(predictions$traj == "trajA", "p53 k.o. A", "p53 k.o. B")

pa <- ggplot(dfdiff, aes(passage, avg_diff, colour = wgd_str, group = wgd_str)) +
  facet_grid(rows = vars(plot_id)) +
  geom_point() + geom_line(data=predictions,aes(y=pred_avg_diff)) +
  scale_size_continuous(range = c(0.4, 2), guide = "none") +
  labs(x = "passage #", y = "chromosomes altered") +
  scale_colour_manual(values = c("WGD-" = "#8E44AD", "WGD+" = "#00724D"), "") +
  #scale_x_log10()+
  #scale_y_log10()+
  common_theme+
  theme(legend.position = "top")
pa

# ───────────────────────────────────────── 5 • plots

dffit$plot_id <- ifelse(dffit$traj == "trajA", "p53 k.o. A", "p53 k.o. B")
dffit_means <- aggregate(fitness ~ wgd_str + dk + plot_id, dffit, mean)
pb <- ggplot(dffit, aes(dk, fitness)) +
  facet_grid(cols = vars(wgd_str), rows = vars(plot_id),scales="free_x") +
  geom_jitter(height = 0, width = 0.15, colour = "grey70", size = 0.6) +
  geom_point(data = dffit_means, colour = "red", size = 5, shape = 95) +
  labs(x = "chromosomes altered", y = "fitness") +
  common_theme
pb
dfl$plot_id <- ifelse(dfl$traj == "trajA", "p53 k.o. A", "p53 k.o. B")
pc <- ggplot(dfl, aes(df, colour = wgd_str)) +
  facet_grid(rows=vars(plot_id))+
  stat_ecdf(size = 0.8) +
  labs(x = "fitness effect", y = "cum. dist.") +
  scale_colour_manual(values = c("WGD-" = "#8E44AD", "WGD+" = "#00724D"), "") +
  common_theme+
  scale_x_continuous(limits=c(-0.05,0.05))+
  theme(legend.position = "top")
pc

plt <- plot_grid(pa, pb, pc, nrow = 1, labels = letters[1:3],
                 rel_widths=c(2,3,2),label_size = base_text_size+2)

ggsave("figs/wgd.pdf",device=cairo_pdf, plt, width = 180, height = 60, units = "mm",dpi=600)

# ───────────────────────────────────────── 6 • console summary
cat("\n== Permutation KS test ==\n")
cat("Observed D =", round(obs_stat, 3), "; p_empirical =", perm_p, "\n")

cat("\n== Piecewise landscape model (lmer) ==\n")
print(summary(mod_piece)$coefficients)


##############################################################################
# 8. EXPORT SOURCE DATA (Nature Requirement) - FIGURE 5
##############################################################################

dir.create("data/source_data", showWarnings = FALSE, recursive = TRUE)

# --- Panel A: Divergence Kinetics (Points & Model Lines) ---
# Source 1: Raw Points (dfdiff)
sd_5a_points <- dfdiff %>%
  select(Experiment_ID = plot_id, Passage = passage, 
         WGD_Status = wgd_str, Chromosomes_Altered = avg_diff) %>%
  mutate(Chromosomes_Altered = round(Chromosomes_Altered, 3))

# Source 2: Model Lines (predictions)
# We only export this if the model converged
sd_5a_lines <- if (!is.null(predictions)) {
  predictions %>%
    select(Experiment_ID = plot_id, Passage = passage, 
           WGD_Status = wgd_str, Predicted_Alterations = pred_avg_diff) %>%
    mutate(Predicted_Alterations = round(Predicted_Alterations, 3))
} else {
  data.frame(Note = "Model fitting did not converge; no line data available.")
}

# Save as a list containing both components
saveRDS(list(Points = sd_5a_points, Model_Fit = sd_5a_lines), "data/source_data/Fig5a.Rds")


# --- Panel B: Landscape Shape (Fitness vs Karyotype Distance) ---
# Source 1: Raw Points (dffit)
sd_5b_points <- dffit %>%
  select(Experiment_ID = plot_id, WGD_Status = wgd_str, 
         Chromosomes_Altered = dk, Fitness = fitness) %>%
  mutate(Fitness = round(Fitness, 4))

# Source 2: Red Mean Bars (dffit_means)
sd_5b_means <- dffit_means %>%
  select(Experiment_ID = plot_id, WGD_Status = wgd_str, 
         Chromosomes_Altered = dk, Mean_Fitness = fitness) %>%
  mutate(Mean_Fitness = round(Mean_Fitness, 4))

saveRDS(list(Raw_Points = sd_5b_points, Mean_Values = sd_5b_means), "data/source_data/Fig5b.Rds")


# --- Panel C: ECDF of Fitness Effects ---
# Source: 'dfl'
sd_5c <- dfl %>%
  select(Experiment_ID = plot_id, WGD_Status = wgd_str, 
         Fitness_Effect = df, Karyotype_ID = k_id) %>%
  mutate(Fitness_Effect = round(Fitness_Effect, 5))

saveRDS(sd_5c, "data/source_data/Fig5c.Rds")

message("Figure 5 Source Data saved to data/source_data/")