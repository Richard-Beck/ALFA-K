if(!basename(getwd())=="ALFA-K") stop("Ensure working directory set to ALFA-K root")
source("R/utils_karyo.R")
source("R/utils_theme.R")
source("R/utils_env.R")
ensure_packages(c("ggplot2","purrr","dplyr","viridisLite","deldir","reshape2","ggdendro","patchwork","magick","cowplot"))

base_text_size <- 5
classic_theme <- make_base_theme("classic",base_text_size)
void_theme <- make_base_theme("void",base_text_size)



simulate_simple_evo <- function(N, K, p = 0.1) {
  # N: pop size; K: # timesteps; p: mut rate per timestep
  birth_time <- rep(1, N)        # initial birth times
  fitness    <- runif(N)         # initial fitnesses
  out <- matrix(NA, nrow = K, ncol = N)
  out[1, ] <- birth_time
  
  for (t in 2:K) {
    # mutation step
    is_mut <- runif(N) < p
    birth_time[is_mut] <- t
    fitness[is_mut]   <- fitness[is_mut]+runif(sum(is_mut))
    
    # selection step (Wright–Fisher)
    keep <- sample(N, N, replace = TRUE, prob = fitness)
    birth_time <- birth_time[keep]
    fitness    <- fitness[keep]
    
    out[t, ] <- birth_time
  }
  
  out
}


do_gen <- function(mat,tt,nx,ny,tmax=14){
  nx <- 10; ny <- 10
  pts <- expand.grid(i=1:nx, j=1:ny)
  pts$x <- pts$i/nx + runif(nrow(pts), -0.5/nx, 0.5/nx)
  pts$y <- pts$j/ny + runif(nrow(pts), -0.5/ny, 0.5/ny)
  
  # 2. compute Voronoi tesselation
  d <- deldir(pts$x, pts$y)
  tiles <- tile.list(d)
  
  # 3. assemble into one data.frame
  df <- do.call(rbind, lapply(seq_along(tiles), function(i) {
    with(tiles[[i]], data.frame(id=i, x=x, y=y))
  }))
  
  ix <- round(nrow(mat)*tt/tmax)
  df$em <- mat[ix,df$id]
  df$time <- tt
  pts$time <- tt
  return(list(df=df,pts=pts))
}
nx=10;ny=10
mat <- simulate_simple_evo(N = nx*ny, K = 50, p = 0.1)
x <- lapply(seq(2,14,2),function(tt) do_gen(mat,tt,nx,ny))

df <- do.call(rbind,lapply(x,function(xi) xi$df))
pts <- do.call(rbind,lapply(x,function(xi) xi$pts))

# 4. plot
p_evo <- ggplot(df, aes(x, y)) +
  facet_grid(cols=vars(paste0("time=",stringr::str_pad(time,2))))+
  geom_polygon(aes(group=id,fill=em), colour="white",show.legend=FALSE) +
  coord_equal() +
  scale_fill_viridis_c()+
  geom_point(data=pts)+
  coord_equal(xlim = c(0.3,.7), ylim = c(0.3,.7), expand = FALSE)+
  void_theme+
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )
p_evo


x0 <- rep(2, 22)
xi <- x0
set.seed(42)
n_pheno     <- 20
time_points <- seq(0, 100, by = 1)
phenotypes  <- paste0("P", seq_len(n_pheno))

for(i in 1:(n_pheno-1)) {
  j <- sample(1:22, 1)
  xi[j] <- xi[j] + (as.numeric(runif(1) > 0.5)*2 - 1)
  x0 <- rbind(x0, xi)
}
rownames(x0) <- phenotypes

# 2) Simulate freq_mat with variable peak times/heights
peak_times   <- seq(20, 80, length.out = n_pheno)
peak_heights <- runif(n_pheno, min = 0.5, max = 1.5)

freq_mat <- sapply(seq_len(n_pheno), function(i) {
  g <- dnorm(time_points, mean = peak_times[i], sd = 10)
  g <- g / max(g) * peak_heights[i]
})
colnames(freq_mat) <- phenotypes
rownames(freq_mat) <- time_points
freq_mat <- freq_mat / rowSums(freq_mat)  # row‑normalize so sums = 1

# 3) Cluster on x0
dist_mat   <- dist(x0)
hcl        <- hclust(dist_mat, method = "average")
phen_order <- phenotypes[hcl$order]

# 4) Prepare the dendrogram panel
dend   <- dendro_data(hcl)
p_dend <- ggplot(segment(dend)) +
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
  coord_flip() +
  scale_y_reverse() +
  labs(x = "Frequent karyotype") +
  classic_theme+
  theme(
    plot.margin = unit(c(0,0,0,0), "cm"),
    axis.title.y = element_text(
      size    = base_text_size,
      angle   = 90,
      vjust   = 0.5
    ),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )
p_dend

# 5) Prepare the karyotype (CNA) heatmap
df_cna <- melt(x0)
names(df_cna) <- c("phenotype","pos","copy_number")
df_cna$phenotype <- factor(df_cna$phenotype, levels = phen_order)

p_cna <- ggplot(df_cna, aes(x = pos, y = phenotype, fill = copy_number)) +
  geom_raster(show.legend = F) +
  scale_fill_viridis_c(option="plasma") +
  scale_x_discrete("chromosome")+
  classic_theme +
  theme(
    axis.title.y      = element_blank(),
    axis.text       = element_blank(),
    axis.ticks      = element_blank(),
    axis.line = element_blank(),
    panel.grid        = element_blank(),
    plot.margin       = unit(c(0,0,0,0), "cm")
  ) +
  # horizontal separators between rows:
  geom_hline(
    yintercept = seq(1.5, n_pheno - 0.5, by = 1),
    color = "white", size = 0.2
  )

# 6) Prepare the frequency heatmap
df_heat <- melt(freq_mat)
names(df_heat) <- c("time","phenotype","frequency")
df_heat$time      <- as.numeric(as.character(df_heat$time))
df_heat$phenotype <- factor(df_heat$phenotype, levels = phen_order)

p_heat <- ggplot(df_heat, aes(x = time/7, y = phenotype, fill = frequency)) +
  geom_raster() +
  scale_fill_viridis_c(option="magma") +
  scale_x_continuous(breaks = seq(0, 14, 2)) +
  labs(x = "Time (weeks)", fill = "Freq.") +
  classic_theme +
  theme(
    axis.title.y      = element_blank(),
    axis.text.y       = element_blank(),
    axis.ticks.y      = element_blank(),
    axis.line=element_blank(),
    panel.grid        = element_blank(),
    plot.margin       = unit(c(0,0,0,0), "cm"),
    legend.position   = c(0.95, 0.95),
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = "white", color = NA)
  )
p_heat
# 7) Combine and align vertically
p_seq <- p_dend + p_cna +
  plot_layout(widths = c(0.5, 1))

mr_unlab <- p_dend+p_cna+p_heat+plot_layout(widths=c(0.5,1,1.5))


crop_fractional <- function(img, xmin = 0, xmax = 1, ymin = 0, ymax = 1) {
  stopifnot(inherits(img, "magick-image"),
            all(c(xmin, xmax, ymin, ymax) >= 0),
            all(c(xmin, xmax, ymin, ymax) <= 1),
            xmin < xmax,
            ymin < ymax)
  
  info <- image_info(img)
  w <- info$width
  h <- info$height
  
  # Convert fractions to pixel-based geometry
  x_px      <- round(xmin * w)
  y_px      <- round(ymin * h)
  width_px  <- round((xmax - xmin) * w)
  height_px <- round((ymax - ymin) * h)
  
  geometry_string <- sprintf("%dx%d+%d+%d", width_px, height_px, x_px, y_px)
  
  image_crop(img, geometry_string)
}

# Load and crop the image (e.g., crop 50px from each edge)
img <- image_read("figs/python_panels/lscape_cells.png")
img <- crop_fractional(img,0.2,0.8,0,1)

landscape_margin_pts <- 7

# Convert to ggdraw image object
img_gg <- ggdraw() + draw_image(img)+
  theme(plot.margin = margin(landscape_margin_pts,landscape_margin_pts,landscape_margin_pts,landscape_margin_pts)) 


img2 <- image_read("figs/python_panels/lscape_explored.png")
img2 <- crop_fractional(img2,0.2,.8,0,1)

# Convert to ggdraw image object
img_gg2 <- ggdraw() + draw_image(img2)+
  theme(plot.margin = margin(landscape_margin_pts,landscape_margin_pts,landscape_margin_pts,landscape_margin_pts)) 

## generate the example 2d landscapes
set.seed(1)
nchrom <- 2
Nwaves <- 10
wls <- c(0.4,0.8)
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
                               barwidth = 2,
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

fpaths <- list.files("data/raw/ABM/",full.names = T)
fpaths <- fpaths[grepl("0p4",fpaths)|grepl("0p8",fpaths)]

tmp <- grepl("0p4",fpaths)
fpaths <- split(fpaths,f=tmp)
fpaths <- as.character(do.call(c,lapply(fpaths,tail,20)))

x <- lapply(fpaths,identify_evolution_intervals)

# Combine the results into one data frame

y <- do.call(rbind, lapply(x, `[[`, "df_train"))
z <- do.call(rbind, lapply(x, `[[`, "df_pred"))

df <- do.call(rbind, lapply(x, `[[`, "df_evo"))


get_w_from_id <- function(id){
  id <- unlist(strsplit(id,split="_"))
  id[1+which(id=="w")]
}

get_rep_from_id <- function(id){
  id <- unlist(strsplit(id,split="_"))
  as.numeric(gsub(".Rds","",id[1+which(id=="rep")]))
}

df$w <- sapply(df$id,get_w_from_id)
df$rep <- sapply(df$id,get_rep_from_id)

y$w <- sapply(y$id,get_w_from_id)
y$rep <- sapply(y$id,get_rep_from_id)

z$w <- sapply(z$id,get_w_from_id)
z$rep <- sapply(z$id,get_rep_from_id)


# Plot intervals using ggplot
p1b <- ggplot(df, aes(y = rep, x = as.numeric(start), xend = as.numeric(end))) +
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
p1b


## plot the R^2 and prediction performance of the ABM sims
results_summary <- readRDS("data/processed/ABM_summaries.Rds")

dfcv <- results_summary$dfl
dfcv <- split(dfcv,f=interaction(dfcv$abmrep,dfcv$ntp,w=dfcv$w))
dfcv <- do.call(rbind,lapply(dfcv,function(dfi) head(dfi[dfi$Rxv==max(dfi$Rxv),],1)))

plt_cv <- ggplot(dfcv,aes(x=as.character(ntp),y=pmax(-1,Rxv)))+
  geom_violin(aes(group=ntp))+
  geom_jitter()+
  coord_flip()+
  scale_y_continuous("CV score")+
  scale_x_discrete("training samples")+
  classic_theme
plt_cv

dfp <- results_summary$dfp
dfp <- dfp[dfp$w%in%c(0.4,0.8),]
dfp <- split(dfp,f=interaction(dfp$abmrep,dfp$ntp,w=dfp$w))
dfp <- do.call(rbind,lapply(dfp,function(dfi) dfi[dfi$Rxv==max(dfi$Rxv),]))
dfp <- dfp[dfp$pos_xv,]
dfpaf <- aggregate(list(fwin = dfp$win),
                   by = list(metric = dfp$metric, passage = dfp$passage, ntp = dfp$ntp,w=dfp$w),
                   mean, na.rm = TRUE)
dfpaf$lambda <- paste0("lambda==",dfpaf$w)
plt_bb <- ggplot(dfpaf[dfpaf$metric=="euclidean"&dfpaf$passage>1,], aes(x = passage-1, y = fwin, color = ntp)) +
  facet_grid(cols = vars(lambda), labeller = label_parsed, scales = "free") +
  geom_line(size = 1) +
  geom_hline(yintercept = 0.5, color = "red", linetype = 2) +
  scale_y_continuous("beat baseline", limits = c(0, 1)) +
  scale_x_continuous("passage number", breaks = seq(2, 10, 2)) +
  scale_color_viridis_d("num.\nsamples")+
  classic_theme+
  theme(legend.position="top")
plt_bb

mr <- cowplot::plot_grid(mr_unlab,nrow=1,labels = c("B"), label_size = base_text_size+2)

r3 <- cowplot::plot_grid(img_gg,img_gg2,labels=c("C","D"), label_size = base_text_size+2,nrow=1)

#r4 <- cowplot::plot_grid(p_grf_complexity,p1b,labels=c("E","F"),nrow=1, label_size = base_text_size+2)

#r5 <- cowplot::plot_grid(plt_cv,plt_bb,labels=c("G","H"),nrow=1, label_size = base_text_size+2)

sec1 <- cowplot::plot_grid(p_evo,mr,r3,nrow=3,labels=c("A","",""), label_size = base_text_size+2,rel_heights = c(1,2,2))
#ggsave("figs/overview.png",sec1,width=120,height=100,units="mm",bg = "white")

#sec2 <- cowplot::plot_grid(p_grf_complexity,p1b,labels=c("E","F"),nrow=1, rel_widths = c(2,3), label_size = base_text_size+2)
#ggsave("figs/overview.png",sec2,width=60,height=50,units="mm",bg = "white")
sec2 <- cowplot::plot_grid(p_grf_complexity,plt_bb,labels=c("E","G"),nrow=1, rel_widths = c(2,2), label_size = base_text_size+2)


#sec3 <-cowplot::plot_grid(plt_cv,plt_bb,labels=c("G","H"),nrow=1, label_size = base_text_size+2)
#ggsave("figs/overview.png",sec3,width=60,height=50,units="mm",bg = "white")

sec4 <- cowplot::plot_grid(sec2,p1b,nrow=2,labels=c("","F"), label_size = base_text_size+2)

plt <- cowplot::plot_grid(sec1,sec4,nrow=1,rel_widths = c(3,2))

#p1 <- cowplot::plot_grid(p_evo,mr,r3,r4,r5,nrow=5,labels=c("A","","","",""), label_size = base_text_size+2,rel_heights = c(1,2,2,2,2))

ggsave("figs/overview.png",plt,width=180,height=100,units="mm",bg = "white")
