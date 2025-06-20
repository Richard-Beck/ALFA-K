if(!basename(getwd())=="ALFA-K") stop("Ensure working directory set to ALFA-K root")
source("R/utils_env.R")
source("R/utils_karyo.R") ## for gen_randscape
ensure_packages("parallel","alfakR","stringr")

setup_and_run_abm <- function(rep_id,wavelength,outDir){
  library(alfakR)
  
  resample_sim <- function(sim, n_samples) {
    mat <- t(as.matrix(sim[,-1]))
    colnames(mat) <- sim$time
    mat <- apply(mat, 2, function(p) rmultinom(1, n_samples, prob = p / sum(p)))
    keep <- rowSums(mat) > 0
    rownames(mat) <- colnames(sim)[-1]
    mat <- mat[keep, , drop = FALSE]
    mat <- mat[order(rowSums(mat),decreasing = T),]
  }
  
  
  founder <- rep(2,22)
  Nwaves <- 10
  
  l <- gen_randscape(founder,Nwaves,wavelength = wavelength)
  times <- c(0,300)
  x0=c(1)
  names(x0) <- paste(founder,collapse=".")
  pmis <- 0.00005
  sim <- run_abm_simulation_grf(
    centroids = l, lambda = wavelength,
    p = pmis, times = times, x0 = x0,
    abm_pop_size = 5e4,abm_max_pop = 2e6, abm_delta_t = 0.1,
    abm_culling_survival = 0.01,
    abm_record_interval = -1, abm_seed = 42,normalize_freq = F
  )
  
  sim_rs <- resample_sim(sim,1000)
  yi <- list(x=data.frame(sim_rs,check.names = F),dt=1)
  out <- list(abm_output=yi,true_landscape=l)
  
  filePath <- paste0(outDir,
                     "w_",gsub("[.]","p",wavelength),
                     "_m_",gsub("[.]","p",pmis),
                     "_rep_",stringr::str_pad(rep_id,width = 2,pad = 0),".Rds")
  saveRDS(out,filePath)
  return(0)
}

outDir <- "data/raw/ABM/"
dir.create(outDir,showWarnings = F,recursive = TRUE)

reps <- 1:100
w <- c(0.2,0.4,0.8,1.6)


# Define parameter grid
param_grid <- expand.grid(rep_id = reps, wavelength = w, KEEP.OUT.ATTRS = FALSE)

# Convert to list of argument lists
param_list <- split(param_grid, seq_len(nrow(param_grid)))

# Create cluster
n_cores <- 50
cl <- makeCluster(n_cores)

# Export necessary variables and load packages
clusterExport(cl, varlist = c("setup_and_run_abm", "outDir"))
clusterEvalQ(cl, library(alfakR))

# Run in parallel
parLapplyLB(cl, param_list, function(params) {
  setup_and_run_abm(params$rep_id, params$wavelength, outDir)
})

# Stop cluster
stopCluster(cl)

