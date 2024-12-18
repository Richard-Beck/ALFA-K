## various functions for aggregating ABM sweep output.

## this function is used when we have data structure
##  inDir
##    subdir
##      target
## where target is a data frame.
aggregate_fit_summaries <- function(summaryName="fit_summaries.Rds",inDir="data/main/",outPath="data/proc/summaries/"){
  ## make sure directory exists
  dir.create(outPath,recursive = T,showWarnings = F)
  ff <- list.files(inDir)
  ff <- paste0(inDir,ff,"/",summaryName)
  ff <- ff[file.exists(ff)]
  x <- lapply(ff,readRDS)
  saveRDS(x,paste0(outPath,summaryName))
  return(x)
}

get_train_test_paths <- function(base_dir) {
  dirs_A <- list.files(base_dir, full.names = TRUE)
  
  do.call(rbind, lapply(dirs_A, function(dir_A) {
    train_path <- file.path(dir_A, "train", "00000")
    test_base <- file.path(dir_A, "test_v2")
    
    if (dir.exists(train_path) && dir.exists(test_base)) {
      test_dirs <- list.dirs(test_base, full.names = TRUE, recursive = TRUE)
      bottom_level_dirs <- test_dirs[sapply(test_dirs, function(d) {
        length(list.dirs(d, recursive = FALSE)) == 0
      })]
      
      base_path <- dir_A
      relative_test_paths <- sub(paste0(base_path, "/"), "", bottom_level_dirs)
      relative_train_path <- sub(paste0(base_path, "/"), "", train_path)
      
      data.frame(
        base_path = rep(base_path, length(relative_test_paths)),
        test_path = relative_test_paths,
        train_path = rep(relative_train_path, length(relative_test_paths)),
        stringsAsFactors = FALSE
      )
    } else {
      NULL
    }
  }))
}

get_train_train_paths <- function(base_dir) {
  dirs_A <- list.files(base_dir, full.names = TRUE)
  
  do.call(rbind, lapply(dirs_A, function(dir_A) {
    train_base <- file.path(dir_A, "train")
    
    if (dir.exists(train_base)) {
      train_dirs <- list.dirs(train_base, full.names = TRUE, recursive = TRUE)
      bottom_level_dirs <- train_dirs[sapply(train_dirs, function(d) {
        length(list.dirs(d, recursive = FALSE)) == 0
      })]
      
      base_path <- dir_A
      relative_train_paths <- sub(paste0(base_path, "/"), "", bottom_level_dirs)
      
      expand.grid(
        base_path = base_path,
        train_path1 = relative_train_paths,
        train_path2 = relative_train_paths,
        stringsAsFactors = FALSE
      )
    } else {
      NULL
    }
  }))
}


compute_population_metrics <- function(metrics=c("angle", "wass"), eval_times=seq(2000,2500,100), inDir="data/main/", outPath="data/proc/summaries/train_test_metrics.Rds", delta_t = 2000, cores = 70) {
  # Load required libraries
  library(parallel)
  
  # Validate input
  if (!all(metrics %in% c("angle", "wass"))) {
    stop("Metrics must be a subset of c('angle', 'wass')")
  }
  
  # Retrieve train-test paths
  df <- get_train_test_paths(inDir)
  if (nrow(df) == 0) {
    stop("No valid train-test paths found in the specified directory.")
  }
  
  # Set up parallel cluster
  cl <- makeCluster(cores)
  
  # Source necessary scripts on all workers
  clusterCall(cl, function() {
    source("utils/comparison_functions.R")
    source("utils/ALFA-K.R")
    library(transport)
  })
  
  # Export variables to the cluster
  clusterExport(cl, varlist = c("metrics", "eval_times", "delta_t", "df"), envir = environment())
  
  # Parallel processing
  res <- do.call(rbind, parLapplyLB(cl, 1:nrow(df), function(i) {
    tryCatch({
      # Extract paths
      train_path <- paste(df$base_path[i], df$train_path[i], sep = "/")
      test_path <- paste(df$base_path[i], df$test_path[i], sep = "/")
      
      # Process simulations
      x0 <- proc_sim(train_path, times = eval_times)
      x1 <- proc_sim(test_path, times = eval_times - delta_t)
      colnames(x1$x) <- as.numeric(colnames(x1$x))+delta_t
      
      # Compute metrics
      a <- w <- c()
      if ("angle" %in% metrics) {
        a <- sapply(eval_times, function(t) angle_metric(x0, x1, t = t))
        names(a) <- paste0("a", eval_times)
      }
      if ("wass" %in% metrics) {
        w <- sapply(eval_times, function(t) wasserstein_metric(x0, x1, t = t))
        names(w) <- paste0("w", eval_times)
      }
      return(c(a, w))
    }, error = function(e) {
      # Return a vector of NAs to ensure consistent output
      return(rep(NA, length(metrics) * length(eval_times)))
    })
  }))
  
  # Stop the cluster
  stopCluster(cl)
  
  # Combine results with the data frame
  df <- cbind(df, res)
  saveRDS(df, outPath)
  return(df)
}

