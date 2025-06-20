if(!basename(getwd())=="ALFA-K") stop("Ensure working directory set to ALFA-K root")
source("R/utils_karyo.R")
source("R/utils_theme.R")
source("R/utils_env.R")
ensure_packages(c("parallel","ggplot2","purrr","dplyr","viridisLite","deldir","reshape2","ggdendro","patchwork","magick","cowplot"))

base_text_size <- 5
classic_theme <- make_base_theme("classic",base_text_size)
void_theme <- make_base_theme("void",base_text_size)

## generate the example 2d landscapes
set.seed(1)
nchrom <- 2
Nwaves <- 10
wls <- c(0.2,0.4,0.8,1.6)
kspace <- expand.grid(x=seq(0,10,0.2),y=seq(0,10,0.2))

df <- do.call(rbind,lapply(wls,function(wavelength){
  l1 <- gen_randscape(sample(1:10,2,replace=T),10,wavelength=wavelength)
  df <- kspace
  df$f <- apply(kspace,1,getf,wavelength=wavelength,tru_lscape=l1)
  df$wavelength <- wavelength
  df
}))

p_grf_complexity <- ggplot(df,aes(x=x,y=y,fill=f))+
  facet_grid(rows=vars(paste0("lambda==",wavelength)),labeller="label_parsed")+
  geom_raster()+
  scale_fill_viridis_c("fitness")+
  scale_x_continuous("chromosome 1",breaks=0:10)+
  scale_y_continuous("chromosome 2",breaks=0:10)+
  guides(fill = guide_colorbar(title = "fitness",
                               frame.colour = "black",
                               barwidth = 3.5,
                               barheight = .5))+
  coord_fixed()+
  classic_theme+
  theme(legend.position = "top")
p_grf_complexity


## generate these "is evolving" plots:

identify_fit_times <- function(x,ntp=c(2,4,8),cutoff=120){
  all_times <- as.numeric(colnames(x))
  tt <- all_times[all_times<cutoff]
  do.call(rbind,lapply(ntp,function(n){
    data.frame(ntp=n,start=min(tail(tt,n)),end=tail(tt,1))
  }))
}

identify_pred_times <- function(x,cutoff=120){
  all_times <- as.numeric(colnames(x))
  tt <- all_times[all_times>cutoff]
  data.frame(end=tail(head(tt,10),1))
}

compute_fitness <- function(simulation_path){
  sim_data  <- readRDS(simulation_path)
  
  w <- (basename(simulation_path) |> strsplit(split="_") |> unlist())[2]
  w <- gsub("p",".",w) |> as.numeric()
  
  x <- sim_data$abm_output$x
  
  f <- sapply(rownames(x),function(kstr){
    strsplit(kstr,split="[.]") |> unlist() |> as.numeric() |> getf(wavelength=w,tru_lscape=sim_data$true_landscape)
  })
  
  f_t <- c(apply(x,2,function(xi) sum(xi*f/sum(xi))))
  
  data.frame(f=as.numeric(f_t),time=as.numeric(names(f_t)),w=w)
}

identify_evolution_intervals <- function(simulation_path) {
  
  sim_data  <- readRDS(simulation_path)
  
  w <- (basename(simulation_path) |> strsplit(split="_") |> unlist())[2]
  w <- gsub("p",".",w) |> as.numeric()
  
  x <- sim_data$abm_output$x
  
  f <- sapply(rownames(x),function(kstr){
    strsplit(kstr,split="[.]") |> unlist() |> as.numeric() |> getf(wavelength=w,tru_lscape=sim_data$true_landscape)
  })
  
  f_t <- c(apply(x,2,function(xi) sum(xi*f/sum(xi))))
  
  # Calculate the difference in fitness
  fitness_diff <- c(NA, diff(f_t))
  
  # Define threshold for considering upward trend
  threshold <- 0  # Can adjust based on noise tolerance
  
  # Find periods where fitness_diff > threshold
  evolving <- which(fitness_diff > threshold)
  
  # Find contiguous intervals of evolution
  intervals <- split(evolving, cumsum(c(TRUE, diff(evolving) != 1)))
  
  # Identify the longest contiguous interval
  longest_interval <- intervals[[which.max(sapply(intervals, length))]]
  
  # Return the start and end times of the interval
  start <- head(names(longest_interval),1)
  end <- tail(names(longest_interval),1)
  df_evo <- data.frame(start = start, end = end, id =basename(simulation_path))
  df_train <- identify_fit_times(x)
  df_train$id <- basename(simulation_path)
  df_pred <- identify_pred_times(x)
  df_pred$id <- basename(simulation_path)
  return(list(df_evo=df_evo,df_train=df_train,df_pred=df_pred))
}
get_w_from_id <- function(id){
  id <- unlist(strsplit(id,split="_"))
  id[1+which(id=="w")]
}

get_rep_from_id <- function(id){
  id <- unlist(strsplit(id,split="_"))
  as.numeric(gsub(".Rds","",id[1+which(id=="rep")]))
}
fpaths <- list.files("data/raw/ABM/",full.names = T)
all_fpaths <- fpaths
wavelengths <- sapply(basename(fpaths),get_w_from_id)

## may want to downsample for visibility...
#fpaths <- split(fpaths,f=wavelengths)
#fpaths <- as.character(do.call(c,lapply(fpaths,tail,20)))

cl <- makeCluster(32)                                         
clusterCall(cl,function(dummyVar) source("R/utils_karyo.R"))
clusterExport(cl,c("identify_fit_times","identify_pred_times"),
              envir = environment())
x <- parLapplyLB(cl,fpaths,identify_evolution_intervals)
df_fit <- parLapplyLB(cl,fpaths,compute_fitness)
for(i in 1:length(df_fit)) df_fit[[i]]$rep_id <- i
df_fit <- do.call(rbind,df_fit)
stopCluster(cl)

# Combine the results into one data frame

y <- do.call(rbind, lapply(x, `[[`, "df_train"))
z <- do.call(rbind, lapply(x, `[[`, "df_pred"))

df <- do.call(rbind, lapply(x, `[[`, "df_evo"))


df$w <- sapply(df$id,get_w_from_id)
df$rep <- sapply(df$id,get_rep_from_id)

y$w <- sapply(y$id,get_w_from_id)
y$rep <- sapply(y$id,get_rep_from_id)

z$w <- sapply(z$id,get_w_from_id)
z$rep <- sapply(z$id,get_rep_from_id)


# Plot intervals using ggplot
p_sample_overview <- ggplot(df, aes(y = rep, x = as.numeric(start), xend = as.numeric(end))) +
  facet_grid(rows = vars(paste0("lambda==", gsub("p", ".", w))), labeller = label_parsed)+
  geom_segment(aes(x = as.numeric(start), xend = as.numeric(end), y = rep, yend = rep)) +
  geom_point(data=y,aes(x = as.numeric(start), y = rep,color=as.character(ntp)),shape=95,size=2)+
  geom_point(data=y,aes(x = as.numeric(end), y = rep),color="grey",shape=95,size=2)+
  geom_point(data=z,aes(x = as.numeric(end), y = rep),color="red",shape=95,size=2)+
  coord_flip()+  # Flip axes for better readability+
  scale_color_viridis_d("num.\nsamples")+
  scale_x_continuous("simulation time (days)")+
  scale_y_continuous("ABM sim. replicate ID",breaks=NULL)+
  classic_theme+
  theme( legend.direction = "horizontal",
         legend.position = c(0.0, 1.0),  # Top-left inside the plot
         legend.justification = c("left", "top"))
p_sample_overview

p_popf <- ggplot(df_fit,aes(x=time,y=f,group=rep_id))+
  facet_grid(cols=vars(paste0("lambda==",w)), labeller = label_parsed)+
  geom_line(alpha=0.5)+
  classic_theme+
  scale_x_continuous("simulation time (days)")+
  scale_y_continuous("mean fitness")
p_popf


melt_for_plotting <- function(x,n){
  x <- head(x[order(rowSums(x),decreasing=T),],n)
  x <- reshape2::melt(
    data.frame(ID = rownames(x), x, row.names = NULL,check.names = F),  # add row-name column
    id.vars = "ID")
  x$variable <- as.numeric(as.character(x$variable))
  x
}

x1 <- readRDS(fpaths[1])$abm_output$x
xm1 <- melt_for_plotting(x1,5)
x2 <- readRDS(fpaths[301])$abm_output$x
xm2 <- melt_for_plotting(x2,5)

p_x1 <- ggplot(xm1,aes(x=variable,y=value,color=ID))+
  geom_line()+
  classic_theme+
  scale_color_discrete("")+
  theme(legend.position = c(0.95, 0.95),
        legend.justification = c(0.95, 0.95),
        plot.margin= margin(5, 5, 5, 5))+
  scale_x_continuous("simulation time (days)")+
  scale_y_continuous("num. in sample")
p_x1

p_x2 <- ggplot(xm2,aes(x=variable,y=value,color=ID))+
  geom_line()+
  classic_theme+
  scale_color_discrete("")+
  theme(legend.position = c(0.95, 0.95),
        legend.justification = c(0.95, 0.95),
        plot.margin= margin(5, 5, 5, 5))+
  scale_x_continuous("simulation time (days)")+
  scale_y_continuous("num. in sample")
p_x2

topr <- cowplot::plot_grid(p_x1,p_x2,labels=c("A","B"),label_size = base_text_size+2,nrow=1)
rt <- cowplot::plot_grid(p_popf,p_sample_overview,labels=c("D","E"),nrow=2,label_size = base_text_size+2,rel_heights = c(1,3))
botr <- cowplot::plot_grid(p_grf_complexity,rt,labels=c("C",""),nrow=1,label_size = base_text_size+2,rel_widths = c(1,3))
plt <- cowplot::plot_grid(topr,botr,rel_heights = c(1,2.2),nrow=2)
ggsave("figs/supp_overview.png",plt,bg="white",width=180,height=160,units="mm")

