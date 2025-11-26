if(!basename(getwd())=="ALFA-K") stop("Ensure working directory set to ALFA-K root")
# --- Utility for s2v ---
source("R/utils_karyo.R")
source("R/utils_env.R")
source("R/utils_theme.R")

ensure_packages(c("ggplot2","tidyverse","igraph","tidygraph","ggraph","ggrepel","reshape2","transport","pbapply","cowplot"))
base_text_size <- 5
base_theme <- make_base_theme()
void_theme <- make_base_theme("void",base_text_size)

# 0. How many unique karyotypes did we predict?
all_landscapes <- list.files("data/processed/salehi/alfak_outputs",pattern = "landscape.Rds",full.names = T,recursive = T)
all_xval <- gsub("landscape.Rds","xval.Rds",all_landscapes)
xv_scores <- sapply(all_xval,function(fi) {
  if(!file.exists(fi)) return(-Inf)
  return(readRDS(fi))
})

all_landscapes <- all_landscapes[xv_scores>0 & !is.na(xv_scores)]
nk <- count_unique_k(all_landscapes)

# 1. Read & filter prediction outputs
procDir <- "data/processed/salehi/alfak_outputs_proc/"
files <- list.files(procDir, full.names = TRUE)
x_list <- lapply(files, readRDS)
x <- do.call(rbind, x_list)
x0 <- x[!duplicated(interaction(x$fi,x$min_obs)),] 
x0 <- x0[!is.na(x0$xv),] ## keep one row per fit (NA vals indicate fit failed), for all CV scores
x0 <- do.call(rbind, lapply(split(x0, x0$fi), function(df) df[df$xv == max(df$xv), ]))

x <- x[!is.na(x$xv) & x$xv > 0 & x$test_treat != "NaN", ]
x <- do.call(rbind, lapply(split(x, x$fi), function(df) df[df$xv == max(df$xv), ]))
rownames(x) <- NULL
x$clean_id <- sub("_l_\\d+_d1_\\d+_d2_\\d+$", "", x$fi)
x <- do.call(rbind, lapply(split(x, x$clean_id), function(df) df[df$xv == max(df$xv), ]))

# 2. Merge metadata
m <- read.csv("data/raw/salehi/metadata.csv")
lut <- m$datasetname
names(lut) <- paste(m$datasetname, m$timepoint, sep = "_")
x$pdx <- lut[sapply(x$fi, function(i) {
  parts <- strsplit(i, "_")[[1]]
  paste(parts[1:(length(parts) - 6)], collapse = "_")
})]

x0$pdx <- lut[sapply(x0$fi, function(i) {
  parts <- strsplit(i, "_")[[1]]
  paste(parts[1:(length(parts) - 6)], collapse = "_")
})]

# 3. Map sample lineages
sample_lineage_map <- c(
  SA1035T                   = "SA1035",
  SA1035U                   = "SA1035",
  SA532                     = "SA532",
  SA609                     = "SA609",
  SA609UnBU                 = "SA609",
  SA906a                    = "p53 k.o",
  SA906b                    = "p53 k.o",
  SA039U                    = "p53 w.t",
  SA535_CISPLATIN_CombinedT = "SA535",
  SA535_CISPLATIN_CombinedU = "SA535",
  SA000 = "SA609",
  SA609R2 = "SA609"
)
x$type <- sample_lineage_map[x$pdx]
x0$type <- sample_lineage_map[x0$pdx]

# 4. Identify evaluation time point
x <- split(x, x$fi)
x <- x[sapply(x, function(df) 1 %in% df$tt)]
x <- lapply(x, function(df) {
  maxp <- if (df$type[1] %in% c("p53 k.o", "p53 w.t")) 25 else 15
  df <- df[df$tt %in% 1:maxp, ]
  df$eval_tp <- FALSE
  df$eval_tp[df$tt == maxp] <- TRUE
  df
})
x <- do.call(rbind, x)
rownames(x) <- NULL

# 5. Define win variable
x$win <- x$base > x$pred
sel <- x$metric %in% c("overlap", "cosine")
x$win[sel] <- (x$base < x$pred)[sel]

# 6. Aggregate for line & bar plots
z <- aggregate(win ~ type + metric + tt, data = x, FUN = mean)
names(z) <- c("lineage", "metric", "day", "fwin")

x_eval <- subset(x, eval_tp)
z_bar <- aggregate(win ~ type + metric, data = x_eval, FUN = mean)
names(z_bar)[3] <- "fwin"
n_count <- aggregate(win ~ type + metric, data = x_eval, FUN = length)
names(n_count)[3] <- "n"
z_bar$n <- n_count$n



# 7. Salehi Plots
p_salehi_lineage <- ggplot(z, aes(day, fwin)) +
  facet_grid(cols = vars(metric)) +
  geom_line(aes(color = lineage, group = lineage)) +
  scale_x_continuous("time (days)") +
  scale_y_continuous("frac. beat.\nbaseline") +
  base_theme+theme(
    legend.justification = c("left"),
    legend.key.width = unit(0.3, "lines"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.position = "top",
    legend.background = element_rect(fill = "transparent", color = NA)
  )
p_salehi_lineage
p_bar <- ggplot(z_bar, aes(type, fwin)) +
  facet_grid(cols = vars(metric)) +
  geom_col(color="black",fill="grey80")+
  coord_flip() +
  scale_y_continuous("fraction beating baseline", breaks = c(0, .5, 1)) +
  scale_x_discrete("") +
  base_theme
p_bar

x0s <- split(x0,f=x0$type)
x0s <- lapply(x0s,function(xi) xi$xv[xi$ntrain>3])
nsamples <- 10000

ci <- sapply(1:nsamples,function(i){
  median(unlist(x0s[sample(1:length(x0s),length(x0s),replace=T)]))
})

quantile(ci,probs=c(0.025,0.5,0.975))

p_violin <- ggplot(x0, aes(pmax(-1,xv), ntrain, group=ntrain)) +
  geom_violin(size=0.5) +
  geom_jitter(width=0, height=0.2,size=0.5) +
  scale_x_continuous("CV score") +
  scale_y_continuous("training\nsamples", breaks=2:9) +
  base_theme
p_violin
p_scatter <- ggplot(subset(x_eval, metric=="wasserstein"),
                    aes(base, pred0)) +
  geom_point() +
  scale_x_continuous("actual distance") +
  scale_y_continuous("predicted\ndistance") +
  base_theme

cor.tst.df <- subset(x_eval, metric=="wasserstein")[,c("base","pred0")]
cor(cor.tst.df$base,cor.tst.df$pred0,method="spear")
cor(cor.tst.df$base,cor.tst.df$pred0,method="pear")
##############################
# Section 1b: Reference ECDFs & Null Bar Plot
##############################

# A) Null sphere CDF
z_euc <- subset(x_eval, metric=="euclidean")
dSphereAngle <- function(theta, N) {
  coef <- integrate(function(t) sin(t)^(N-2), 0, pi)$value
  sin(theta)^(N-2)/coef
}
cSphereAngle <- function(theta, N) {
  if(!is.finite(theta)) return(1)
  integrate(function(t) dSphereAngle(t, N), 0, theta)$value
}
nulldf <- data.frame(angle=0:180)
nulldf$rads <- pi * nulldf$angle/180
nulldf$CDF  <- sapply(nulldf$rads, cSphereAngle, N=22)

# B) Compute reference-angle pairs (exact from Script 2)
input_dir <- "data/processed/salehi/alfak_inputs/"
ff <- list.files(input_dir)

pairs_list <- lapply(m$uid, function(id) m[m$parent==id & !is.na(m$parent),])
pairs_list <- pairs_list[sapply(pairs_list,nrow)>1]
pairs_df <- do.call(rbind, lapply(pairs_list, function(xi){
  combos <- combn(nrow(xi),2)
  do.call(rbind, lapply(seq_len(ncol(combos)), function(ci){
    i <- combos[1,ci]; j <- combos[2,ci]
    fi <- paste(c(xi$datasetname[i], xi$timepoint[i], "l_2"), collapse="_")
    fi <- ff[grepl(fi, ff)]
    fj <- paste(c(xi$datasetname[j], xi$timepoint[j], "l_2"), collapse="_")
    fj <- ff[grepl(fj, ff)]
    if(length(fi)==0||length(fj)==0) return(NULL)
    data.frame(fi=fi, fj=fj,
               on_treat=paste0(xi$on_treatment[c(i,j)], collapse=""),
               pdx=xi$PDX_id[1],
               stringsAsFactors=FALSE)
  }))
}))

compute_metrics <- function(x1,x2){
  clones <- union(rownames(x1), rownames(x2))
  get_counts <- function(mat,col){
    v <- numeric(length(clones)); names(v)<-clones
    present <- intersect(rownames(mat), clones)
    v[present] <- mat[present,col]
    v
  }
  v11<-get_counts(x1,1); v12<-get_counts(x1,2)
  v21<-get_counts(x2,1); v22<-get_counts(x2,2)
  k <- do.call(rbind, lapply(clones, s2v))
  vec <- function(y) colSums(k*y)
  mag <- function(v) sqrt(sum(v^2))
  getangle <- function(a,b) 180*acos(sum(a*b)/(mag(a)*mag(b)))/pi
  overlap<-function(a,b){ya<-a/sum(a); yb<-b/sum(b); sum(pmin(ya,yb))/min(sum(ya),sum(yb))}
  euclid <-function(a,b){ya<-a/sum(a); yb<-b/sum(b); sqrt(sum((vec(ya)-vec(yb))^2))}
  cosine <-function(a,b){ya<-a/sum(a); yb<-b/sum(b); sum(ya*yb)/(sqrt(sum(ya^2))*sqrt(sum(yb^2)))}
  wass   <-function(a,b){p<-wpp(k,mass=a/sum(a)); q<-wpp(k,mass=b/sum(b)); wasserstein(p,q)}
  funs <- list(overlap=overlap, euclidean=euclid, cosine=cosine, wasserstein=wass)
  df <- do.call(rbind, lapply(names(funs), function(mn){
    data.frame(metric=mn,
               value=funs[[mn]](v12,v22),
               baseline1=funs[[mn]](v11,v12),
               baseline2=funs[[mn]](v21,v22),
               stringsAsFactors=FALSE)
  }))
  delta1 <- vec(v12/sum(v12)) - vec(v11/sum(v11))
  delta2 <- vec(v22/sum(v22)) - vec(v21/sum(v21))
  df <- rbind(df, data.frame(metric="angle",
                             value=getangle(delta1,delta2),
                             baseline1=NA, baseline2=NA,
                             stringsAsFactors=FALSE))
  rownames(df)<-NULL; df
}

y <- do.call(rbind, lapply(seq_len(nrow(pairs_df)), function(i){
  x1 <- readRDS(file.path(input_dir, pairs_df$fi[i]))$x
  x2 <- readRDS(file.path(input_dir, pairs_df$fj[i]))$x
  cbind(pairs_df[i,,drop=FALSE], compute_metrics(x1,x2))
}))
# Map to lineage for these reference pairs
sample_lineage_map2 <- c(
  SA1035 = "SA1035", SA532 = "SA532",
  SA609  = "SA609",  SA906 = "p53 k.o",
  SA039  = "p53 w.t", SA535 = "SA535"
)
y$lineage   <- sample_lineage_map2[y$pdx]
y$treatdiff <- sapply(y$on_treat, function(i) length(unique(strsplit(i,"")[[1]])))
y$metricDir <- 1; y$metricDir[y$metric %in% c("cosine","overlap")] <- -1
y$win1      <- y$metricDir * (y$value - y$baseline1)
y$win2      <- y$metricDir * (y$value - y$baseline2)

# Melt including lineage
y2 <- reshape2::melt(y,
                     id.vars       = c("lineage","fi","fj","on_treat","pdx","metric","value","baseline1","baseline2","treatdiff","metricDir"),
                     measure.vars = c("win1","win2"),
                     variable.name = "which_win",
                     value.name    = "win"
)

# C) Also compute PDX‐pair angle ECDF (script2's p_angle_x data)
mm <- m
mm$fi <- paste0(mm$datasetname,"_",mm$timepoint,"_l_2")
mm$fi <- sapply(mm$fi, function(ii) ff[grepl(ii, ff)])
mm <- mm[sapply(mm$fi, length)>0, ]
combos_mm <- combn(seq_len(nrow(mm)), 2)
combos_mm <- combos_mm[, apply(combos_mm,2, function(ci) length(unique(mm$PDX_id[ci]))>1)]
z_pairs <- do.call(rbind, pbapply::pblapply(seq_len(ncol(combos_mm)), function(i){
  j1 <- combos_mm[1,i]; j2 <- combos_mm[2,i]
  x1 <- readRDS(file.path(input_dir, mm$fi[j1]))$x
  x2 <- readRDS(file.path(input_dir, mm$fi[j2]))$x
  angle <- (function(a,b){
    clones <- union(rownames(a),rownames(b))
    get_col <- function(mat,j){
      v<-numeric(length(clones)); names(v)<-clones
      idx<-match(rownames(mat),clones); v[idx]<-mat[,j]; v
    }
    c11<-get_col(a,1); c12<-get_col(a,2)
    c21<-get_col(b,1); c22<-get_col(b,2)
    y11<-c11/sum(c11); y12<-c12/sum(c12)
    y21<-c21/sum(c21); y22<-c22/sum(c22)
    kmat<-do.call(rbind,lapply(clones,s2v))
    d1<-colSums(kmat*(y12-y11)); d2<-colSums(kmat*(y22-y21))
    acos(sum(d1*d2)/(sqrt(sum(d1^2))*sqrt(sum(d2^2))))*180/pi
  })(x1,x2)
  data.frame(PDX1=mm$PDX_id[j1], PDX2=mm$PDX_id[j2], angle=angle, stringsAsFactors=FALSE)
}))
#last_char <- function(s) substr(s, nchar(s), nchar(s))
#tmp <- z_euc[]
# --- Plot: combined ECDF ---
p_ecdf <- ggplot() +
  stat_ecdf(data=z_euc,      aes(angle,color="predicted"),     geom="step", size=1.5)   +
  stat_ecdf(data=data.frame(angle=subset(y, metric=="angle")$value),
            aes(angle, color="forked"),     geom="step", size=1) +
  stat_ecdf(data=z_pairs,     aes(angle,color="unrelated"),     geom="step", size=1) +
  geom_line(data=nulldf,      aes(angle, CDF,color="theoretical"), size=1, linetype="dashed") +
  scale_color_manual("", values = c("#999999", "#E69F00", "#0072B2", "#009E73"))+
  labs(x="Angle (degrees)", y="cumulative\nfrequency") +
  base_theme+theme(legend.position=c(.99,0.01),
                   legend.justification = c("right", "bottom"),
                   legend.key.size = unit(0.5, "lines"),
                   legend.spacing.y = unit(0, "pt"),
                   legend.background = element_rect(fill = "transparent", color = NA))
p_ecdf
# Optional: KS test

source("R/utils_stats.R")

pred_res  <- sphere_null_test(z_euc$angle, z_euc$type)  # predictions
sis_res   <- sphere_null_test(y$value[y$metric=="angle"],
                               y$lineage[y$metric=="angle"])   
sis_res_treat_match   <- sphere_null_test(y$value[y$metric=="angle"&y$on_treat%in%c("yy","nn")],
                              y$lineage[y$metric=="angle"&y$on_treat%in%c("yy","nn")])  

sis_res_treat_mismatch   <- sphere_null_test(y$value[y$metric=="angle"&y$on_treat%in%c("yn","ny")],
                                          y$lineage[y$metric=="angle"&y$on_treat%in%c("yn","ny")])   
# sisters
unrel_res <- sphere_null_test(z_pairs$angle,
                              paste(pmin(z_pairs$PDX1, z_pairs$PDX2),
                                    pmax(z_pairs$PDX1, z_pairs$PDX2),
                                    sep = "--")) # unrelated

print(pred_res)
print(sis_res)
print(unrel_res)

# --- Null‐reference bar plot formatted like salehi_bar ---
yagg_ref <- aggregate(list(win=y2$win),
                      by=list(lineage=y2$lineage,metric=y2$metric),
                      FUN=function(i) mean(i<0))
yagg_ref$n <- aggregate(list(n=y2$win),
                        by=list(lineage=y2$lineage,metric=y2$metric),
                        FUN=length)$n

lineage_levels <- levels(factor(z_bar$type))
combos_ref <- expand.grid(lineage=lineage_levels,
                          metric=unique(yagg_ref$metric[yagg_ref$metric!="angle"]),
                          stringsAsFactors=FALSE)
plot_df_ref <- merge(combos_ref,
                     yagg_ref[yagg_ref$metric!="angle",c("lineage","metric","win","n")],
                     by=c("lineage","metric"), all.x=TRUE)
plot_df_ref$win[is.na(plot_df_ref$win)] <- 0
plot_df_ref$n  [is.na(plot_df_ref$n)]   <- 0
plot_df_ref$lineage <- factor(plot_df_ref$lineage, levels=lineage_levels)

p_ref_null <- ggplot(plot_df_ref, aes(lineage, win)) +
  facet_grid(cols=vars(metric)) +
  geom_col(color="black",fill="grey80")+
  coord_flip() +
  scale_x_discrete("") +
  scale_y_continuous("fraction beating baseline",breaks=c(0,0.5,1))+
  base_theme
p_ref_null
z_bar_tmp <- z_bar
colnames(z_bar_tmp) <- colnames(plot_df_ref)
z_bar_tmp$type <- "predicted"
plot_df_ref$type <- "forked"

z_bar_tmp <- rbind(z_bar_tmp,plot_df_ref)
z_bar_tmp <- z_bar_tmp[z_bar_tmp$metric=="cosine",]## removes duplicates, 
z_bar_tmp$n[z_bar_tmp$type=="forked"] <- z_bar_tmp$n[z_bar_tmp$type=="forked"]/2

p_pass_count <- ggplot(z_bar_tmp, aes(lineage, n)) +
  facet_grid(cols=vars(type)) +
  geom_col(color="black",fill="grey80")+
  coord_flip() +
  scale_x_discrete("")+
  scale_y_continuous("number of passages")+
  scale_fill_viridis_c("num.\nsub-lins.", option="magma") +
  base_theme
p_pass_count

print("Sister baseline beating fraction:")
sum(plot_df_ref$win*plot_df_ref$n)/sum(plot_df_ref$n)
print("prediction baseline beating fraction:")
sum(z_bar_tmp$win*z_bar_tmp$n)/sum(z_bar_tmp$n)
##############################
# Section 2: Network/Dendrogram Plot 
##############################

source("R/utils_lineage.R")
# --- Process metadata for network plot ---
sample_lineage_map2 <- c(
  SA1035 = "SA1035",
  SA532  = "SA532",
  SA609  = "SA609",
  SA906  = "p53 k.o",
  SA039  = "p53 w.t",
  SA535  = "SA535"
)
df <- read.csv("data/raw/salehi/metadata.csv")
df$shortName <- paste(df$datasetname,df$timepoint,sep="_")
x0$shortName <- sapply(x0$fi,function(i){
  strsplit(i,split="_l_") |> unlist() |> head(1)
})
x0_tmp <- do.call(rbind,lapply(split(x0,f=x0$shortName),function(xi) xi[xi$ntrain==max(xi$ntrain),]))
df <- merge(df,x0_tmp[,c("xv","shortName","min_obs","fi")],all=T)
colnames(df)[colnames(df)=="xv"] <- "xval"
df$xval <- pmax(-1,df$xval)
df <- prune_children(df)
df$linlab <- sample_lineage_map2[df$PDX_id]

# Create dummy root to fix missing parent values
dummy_uid <- "ROOT"
dummy_row <- df[1, ]
dummy_row$uid <- dummy_uid
dummy_row$datasetname <- "ROOT"
df_fixed <- rbind(df, dummy_row)
df_fixed$parent[is.na(df_fixed$parent) | df_fixed$parent == ""] <- dummy_uid

# Define lineage for network; mark dummy nodes
lm <- c("SA039U" = "p53 wt", "SA906a" = "p53 ko", "SA906b" = "p53 ko")
df_fixed$lineage <- ifelse(df_fixed$datasetname %in% names(lm), lm[df_fixed$datasetname], "TNBC PDX")
df_fixed$lineage[grepl("x", df_fixed$uid)] <- "dummy"
df_fixed$lineage[df_fixed$uid == dummy_uid] <- "dummy"
df_fixed <- cbind(df_fixed[, c("uid", "parent")],
                  df_fixed[, !colnames(df_fixed) %in% c("uid", "parent")])

# Build edge list and graph ----------------------------------------------
edges_df <- df_fixed %>% 
  filter(uid != dummy_uid) %>% 
  mutate(from_dummy = (parent == dummy_uid)) %>% 
  select(from = parent, to = uid, treatment = on_treatment, from_dummy)
g <- graph_from_data_frame(d = edges_df, vertices = df_fixed, directed = TRUE)

tg <- as_tbl_graph(g) %>% 
  mutate(dummy = (lineage == "dummy"),
         node_type = ifelse(dummy, "dummy", ifelse(parent == dummy_uid, "root", "non-root")))

# Define layout and compute label positions -------------------------------
set.seed(42)
layout <- create_layout(tg, layout = "dendrogram", circular = TRUE,
                        height = -node_distance_to(1, mode = "all"))
# Swap x/y and flip y (sometimes necessary)
#tmp <- layout$x; layout$x <- layout$y; layout$y <- tmp
#layout$y <- -layout$y

# Compute label positions for nodes whose parent is the dummy root
radius_multipliers <- c(
  "p53 w.t" = 5,
  "SA1035"  = 6,
  "SA532"   = 10,
  "SA535"   = 8,
  "SA609"   = 5,
  "p53 k.o" = 8
)
lineage_labels <- layout %>% 
  filter(parent == dummy_uid) %>% 
  mutate(
    r = sqrt(x^2 + y^2),
    theta = atan2(y, x),
    r_multiplier = radius_multipliers[linlab],
    r_label = r * r_multiplier,
    x_label = r_label * cos(theta),
    y_label = r_label * sin(theta)
  )
lineage_labels <- layout %>% 
  filter(parent == dummy_uid) %>% 
  mutate(x_label = x,y_label = y)

lineage_positions <- data.frame(x_label=c(0.08,0.4,0.5,-0.1,-0.2,-0.3),
                                y_label=c(0.6,0.6,0.2,-0.3,-0.1,0.8),
                                row.names = lineage_labels$linlab)

lineage_labels$x_label <- lineage_positions[lineage_labels$linlab,"x_label"]
lineage_labels$y_label <- lineage_positions[lineage_labels$linlab,"y_label"]

# --- Plot 6: Dendrogram/Network Plot ---
p_network <- ggraph(layout) + 
  # Connector lines from node to label
  geom_segment(data = lineage_labels,
               aes(x = x, y = y, xend = x_label, yend = y_label),
               color = "gray80", linetype = "dashed") +
  geom_edge_link(aes(color = ifelse(from_dummy, "dummy", treatment)),
                 edge_width = 1) +
  # Plot dummy nodes as white points
  geom_node_point(data = filter(layout, dummy),
                  aes(x = x, y = y), color = "white") +
  # Plot non-dummy nodes colored by xval
  geom_node_point(data = filter(layout, !dummy),
                  aes(x = x, y = y, color = xval), size = 2) +
  scale_edge_color_manual(name = "",
                          values = c("y" = "red", "n" = "black", "dummy" = "white"),
                          labels = c("y" = "Treatment ON", "n" = "Treatment OFF", "dummy" = "")) +
  scale_color_viridis_c("CV\nscore") +
  # Place repelled text labels at computed positions
  geom_text(data = lineage_labels,
            aes(x = x_label, y = y_label, label = linlab),
            fontface = "bold",size = base_text_size / .pt,family = "sans") +
  void_theme+
  coord_equal()+
  theme(legend.position = c(0.0, 0.0),
        legend.justification = c("left", "bottom"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box = "vertical",
        legend.spacing.y = unit(0, "cm"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA)
  ) 
p_network

##############################
# Final Output
##############################
plots <- list(
  salehi_lineage = p_salehi_lineage,
  salehi_bar     = p_bar,
  salehi_violin  = p_violin,
  salehi_scatter = p_scatter,
  salehi_ecdf    = p_ecdf,
  ref_null       = p_ref_null,
  network_dendro = p_network,
  pass_count=p_pass_count
)

sec1 <- plot_grid(plots$salehi_violin,plots$salehi_ecdf,plots$salehi_scatter, labels = c("b","c", "d"),label_size = base_text_size+2,nrow=3)
sec2 <- plot_grid(plots$network_dendro,sec1, labels = c("a","", ""),label_size = base_text_size+2,nrow=1,rel_widths = c(2,1))

sec3 <- plot_grid(plots$salehi_lineage,plots$salehi_bar, plots$ref_null,plots$pass_count, labels = c("e","f","g","h"),
                  label_size = base_text_size+2, nrow=4,rel_heights = c(3,2,2,2))

plt <- plot_grid(sec2,sec3,nrow=1,rel_widths=c(2,1))

ggsave("figs/salehi_validation.png",plot=plt,width=180,height=80,units="mm", 
       bg = "white")
ggsave("figs/salehi_validation.pdf",device = cairo_pdf,plt,width=180,height=80,units="mm",bg = "white",dpi=600)


# [ ... Previous code remains exactly the same ... ]

##############################################################################
# 8. EXPORT SOURCE DATA (Nature Requirement)
#    Saves cleaned, minimal dataframes to data/source_data/
##############################################################################

dir.create("data/source_data", showWarnings = FALSE, recursive = TRUE)

# --- Panel A: Network (Nodes & Edges) ---
# FIX: 'uid' is converted to 'name' by igraph/ggraph during layout creation
sd_2a_nodes <- layout %>% 
  select(uid = name, lineage, CV_Score = xval, x, y) %>%
  mutate(CV_Score = round(CV_Score, 3))

sd_2a_edges <- edges_df %>%
  select(Parent = from, Child = to, Treatment = treatment) %>%
  mutate(Treatment = ifelse(Treatment == "y", "On_Treatment", "Off_Treatment"))

saveRDS(list(nodes=sd_2a_nodes, edges=sd_2a_edges), "data/source_data/Fig2a.Rds")

# --- Panel B: Violin Plot (Model Fit) ---
sd_2b <- x0 %>%
  select(Num_Training_Samples = ntrain, CV_Score = xv) %>%
  mutate(CV_Score = round(CV_Score, 4))
saveRDS(sd_2b, "data/source_data/Fig2b.Rds")

# --- Panel C: ECDF (Model vs Sisters vs Unrelated) ---
# Combine the 4 layers of the plot into one readable table
sd_2c <- rbind(
  data.frame(Angle = z_euc$angle, Group = "Predicted"),
  data.frame(Angle = subset(y, metric=="angle")$value, Group = "Forked_Sister"),
  data.frame(Angle = z_pairs$angle, Group = "Unrelated"),
  # For theoretical, we provide the calculated curve points
  data.frame(Angle = nulldf$angle, Group = "Theoretical_Null") 
) %>%
  mutate(Angle = round(Angle, 2))
saveRDS(sd_2c, "data/source_data/Fig2c.Rds")

# --- Panel D: Scatter Plot (Wasserstein) ---
sd_2d <- subset(x_eval, metric=="wasserstein") %>%
  select(Actual_Distance = base, Predicted_Distance = pred0) %>%
  mutate(
    Actual_Distance = round(Actual_Distance, 4),
    Predicted_Distance = round(Predicted_Distance, 4)
  )
saveRDS(sd_2d, "data/source_data/Fig2d.Rds")

# --- Panel E: Line Plot (Time Series) ---
sd_2e <- z %>%
  select(Lineage = lineage, Metric = metric, Time_Days = day, Fraction_Beating_Baseline = fwin) %>%
  mutate(Fraction_Beating_Baseline = round(Fraction_Beating_Baseline, 3))
saveRDS(sd_2e, "data/source_data/Fig2e.Rds")

# --- Panel F: Bar Plot (Aggregate Performance) ---
sd_2f <- z_bar %>%
  select(Lineage = type, Metric = metric, Fraction_Beating_Baseline = fwin) %>%
  mutate(Fraction_Beating_Baseline = round(Fraction_Beating_Baseline, 3))
saveRDS(sd_2f, "data/source_data/Fig2f.Rds")

# --- Panel G: Reference Null Bar Plot ---
sd_2g <- plot_df_ref %>%
  select(Lineage = lineage, Metric = metric, Fraction_Beating_Baseline = win) %>%
  mutate(Fraction_Beating_Baseline = round(Fraction_Beating_Baseline, 3))
saveRDS(sd_2g, "data/source_data/Fig2g.Rds")

# --- Panel H: Passage Counts ---
sd_2h <- z_bar_tmp %>%
  select(Lineage = lineage, Type = type, Count = n)
saveRDS(sd_2h, "data/source_data/Fig2h.Rds")

message("Source data saved to data/source_data/")

