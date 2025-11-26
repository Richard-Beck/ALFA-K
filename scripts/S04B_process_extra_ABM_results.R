## Extra ABM results incorporating the extra missegregation rate simulations, 
## also the accuracy of the initial fraction estimation, which were requested 
## by a reviewer at the 11th hour

## the following code compares the intiial fractions estimated to those
## present in the sample data. The results are excellent. We do not show these plots
inDir <- "data/processed/ABM"
simDir <- "data/raw/ABM/"
sims <- list.files(inDir)

df <- pbapply::pblapply(sims,function(fi){
  data_path <- file.path(simDir,paste0(fi,".Rds"))
  
  x <- readRDS(data_path)$abm_output$x
  
  pass_times_all <- as.numeric(colnames(x))
  pass_times <- pass_times_all[pass_times_all<120]
  pass_times <- tail(pass_times,8)
  
  x <- x[,as.numeric(colnames(x))%in%pass_times]
  
  tt <- pass_times[1]
  
  xv_path <- file.path(inDir,fi,"minobs_10_ntp_8","xval.Rds")
  if(!file.exists(xv_path)) return(NULL)
  xv <- readRDS(xv_path)
  
  boot_path <- file.path(inDir,fi,"minobs_10_ntp_8","bootstrap_res.Rds")
  if(!file.exists(boot_path)) return(NULL)
  fit <- readRDS(boot_path)
  
  f_mat <- fit$final_fitness
  x0_mat <- fit$final_frequencies
  
  proj_mat <- x0_mat * exp(f_mat * tt)
  freq_mat <- proj_mat / rowSums(proj_mat)
  
  est_freq_at_tt <- colMeans(freq_mat, na.rm = TRUE)
  
  true_counts_at_tt <- x[, 1] 
  names(true_counts_at_tt) <- rownames(x)
  common_karyos <- names(est_freq_at_tt)
  
  true_counts_subset <- true_counts_at_tt[common_karyos]
  true_freq_at_tt <- true_counts_subset / sum(true_counts_subset)
  
  df1 <- data.frame(true=true_freq_at_tt,est=est_freq_at_tt)
  df2 <- data.frame(correlation = cor(est_freq_at_tt, true_freq_at_tt),
                    rmse = sqrt(mean((est_freq_at_tt - true_freq_at_tt)^2)),
                    nfq=length(common_karyos),
                    xv = xv,
                    id=fi)
  
  return(list(df1=df1,df2=df2))
})

df1 <- do.call(rbind,lapply(df,function(di) di$df1))
df2 <- do.call(rbind,lapply(df,function(di) di$df2))

p <- ggplot(df1,aes(x=true,y=est))+
  geom_point()
p

p <- ggplot(df1[df1$true==0,],aes(x=est))+
  geom_histogram()
p


median(df2$correlation)

## this one generates data for the plot of varying missegregation rate:

abmdir <- "data/raw/ABM_varyMisseg/"
fitDir <- "data/processed/ABM_vary_misseg/"
source("R/utils_karyo.R")
source("R/utils_theme.R")
source("R/utils_env.R")
ensure_packages(c("ggplot2","stringr"))

base_text_size <- 5
common_theme <- make_base_theme("classic",base_text_size)

ff <- list.files(fitDir)

df <- do.call(rbind,pbapply::pblapply(ff,function(fi){
  lscape_path <- file.path(fitDir,fi,"minobs_10_ntp_8/landscape.Rds")
  if(!file.exists(lscape_path)) return(NULL)
  lscape <- readRDS(lscape_path)
  
  xv_path <- file.path(fitDir,fi,"minobs_10_ntp_8/xval.Rds")
  if(!file.exists(xv_path)) return(NULL)
  xv <- readRDS(xv_path)
  
  
  ## extract metadata from filenames
  abm_ids <- strsplit(fi,split="_") |> unlist()
  wavelength <- gsub("p",".",abm_ids[2]) |> as.numeric()
  pmis <- as.numeric(abm_ids[4])
  
  xi <- readRDS(file.path(abmdir,paste0(fi,".Rds")))
  tru_lscape <- xi$true_landscape
  #Rxv <- readRDS(fx)
  
  f <- sapply(lscape$k,function(k){
    kvec <- strsplit(k,split="[.]") |> unlist() |> as.numeric()
    getf(kvec,wavelength,tru_lscape)
  })
  
  pmis <- as.numeric(unlist(strsplit(fi,split="_"))[4])
  
  data.frame(pearson=cor(f,lscape$mean),Rxv=xv, pmis=pmis)
}))

# 1. Create the factor with levels ordered by the numeric 'pmis' value
#    This ensures 5e-05 comes before 1e-04, etc.
df$pmis_factor <- factor(format(df$pmis, scientific = TRUE),
                         levels = format(sort(unique(df$pmis)), scientific = TRUE))

saveRDS(df,"data/processed/ABM_vary_misseg_summaries.Rds")