#setwd("~/projects/ALFA-K")
setwd("~/projects/008_birthrateLandscape/ALFA-K")

Nreps <- 1
gen_abm_landscape <- function(fit){
  knots <- fit$knots
  cc <- fit$c
  d <- fit$d
  fscape <- rbind(cbind(knots,cc),c(d))
  return(fscape)
}
source("utils/sim_setup_functions.R")
source("utils/ALFA-K.R")
options(scipen=999)
dir <- "data/salehi/alfak_fits_minobs_adaptive_exclude_final/"
out_dir <- "data/salehi/forward_sims/"

cfig_template <- readLines("data/salehi/config_template.txt")
conditions <- list.files(dir)


lapply(conditions, function(fi){
  xi <- readRDS(paste0(dir,fi))
  
  xi <- xi[!sapply(xi,is.null)]
  if(length(xi)==0) return(NULL)
  cor<-sapply(xi, function(xij) xij$vx_cor)
  n <-sapply(xi, function(xij) nrow(xij$xv_res))
  
  goodness <- cor*n
  xi <- xi[[which.max(goodness)]]
  
  cpp_cmd <- "ABM/bin/ABM.exe"
  
  
  id <- unlist(strsplit(fi,split=".Rds"))[1]
  dir.create(paste0(out_dir,id))
  abm_out_dir <- paste0(out_dir,id,"/output/")
  dir.create(abm_out_dir)
  
  x <- readRDS(paste0("data/salehi/alfak_inputs/",fi))
  dt <- x$dt
  x0 <- x$x[,ncol(x$x)-1]
  x0 <- x0[x0>0]
  pop0 <- cbind(do.call(rbind,lapply(names(x0),s2v)),x0)
  rownames(pop0) <- NULL
  colnames(pop0) <- NULL
  pop0[,ncol(pop0)] <- round(100000*pop0[,ncol(pop0)]/sum(pop0[,ncol(pop0)]))
  
  lscape <- gen_abm_landscape(xi$fit)
  
  cfig_path <- paste0(out_dir,id,"/cfig.txt")
  lscape_path <- paste0(out_dir,id,"/lscape.txt")
  pop_path <- paste0(out_dir,id,"/pop.txt")
  
  cfig <- modify_config("fitness_landscape_file",lscape_path,cfig_template)
  cfig <- modify_config("output_dir",abm_out_dir,cfig)
  cfig <- c(cfig,paste0("population_file,",pop_path))
  
  writeLines(cfig,cfig_path)
  write.table(pop0,pop_path,sep=",",quote = F,col.names = F,row.names = F)
  write.table(lscape,lscape_path,sep=",",quote = F,col.names = F,row.names = F)
  
  cmd <- paste(cpp_cmd,cfig_path)
  
  for(i in 1:Nreps) (system(cmd))
})



