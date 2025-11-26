if(!basename(getwd())=="ALFA-K") stop("Ensure working directory set to ALFA-K root")
source("R/utils_env.R")
ensure_packages("alfakR")

inDir <- "data/raw/ABM_varyMisseg/"
outDir_base <- "data/processed/ABM_vary_misseg/"

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1)
  stop("Usage: Rscript name_of_this_script.R <input_file>")

fi_full <- args[1]
if (!file.exists(fi_full)) stop("Input file not found: ", fi_full)

fi <- basename(fi_full)# keep only filename

misrate <- as.numeric((strsplit(fi,split="_")|>unlist())[4])
minobs <- 10
ntp <- 8

xi <- readRDS(paste0(inDir,fi))
yi <- xi$abm_output

outDir <- paste0(outDir_base,gsub(".Rds","",fi),"/minobs_",minobs,"_ntp_",ntp)

pass_times_all <- as.numeric(colnames(yi$x))
pass_times <- pass_times_all[pass_times_all<120]
pass_times <- tail(pass_times,ntp)

yi$x <- yi$x[,pass_times_all%in%pass_times]

R2 <- alfak(yi,outdir = outDir,passage_times = pass_times,minobs=minobs,n0 = 2e5,nb=2e7,pm=misrate)


