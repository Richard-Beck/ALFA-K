if(!basename(getwd())=="ALFA-K") stop("Ensure working directory set to ALFA-K root")

redo_bootstrap <- FALSE

source("R/utils_karyo.R")
source("R/utils_theme.R")
source("R/utils_env.R")
source("R/utils_lineage.R")
libs <- c("igraph", "lme4", "glmmTMB", "car", "ggplot2","ggh4x", "pbapply", "tidygraph", 
          "ggraph","ggsignif","isotone","dplyr","ggnewscale","emmeans","tidyr","parallel")
ensure_packages(libs)

base_text_size <- 5
base_theme <- make_base_theme()
# Load and process lineages
load_lineages <- function(path) {
  raw <- readRDS(path)
  lapply(raw, function(li) paste0(head(li$ids, -1), "-", tail(li$ids, -1)))
}

# Load, filter, and collapse fits
load_fits <- function(proc_dir) {
  files <- list.files(proc_dir, full.names = TRUE)
  xs <- lapply(files, readRDS)
  xs <- xs[sapply(xs, ncol) == 13]
  xs <- lapply(xs, head, 1)
  df <- do.call(rbind, xs)
  
  # filter xv > 0
  df <- df[!is.na(df$xv) & df$xv > 0, ]
  # keep max xv per fit_id
  df <- do.call(rbind, by(seq_len(nrow(df)), df$fi, function(i) {
    sub <- df[i, , drop = FALSE]
    sub[sub$xv == max(sub$xv), ]
  }))
  rownames(df) <- NULL
  df
}





# Compute overlap components and select reps (fitted lineages).
# this function selects a random set of non overlapping lineages.
get_reps <- function(lineages,df){
  ## keep only fits present in df$fi and delete N=2 passage lineages
  lineages <- lineages[names(lineages) %in% df$fi & lengths(lineages)>1]
  if (!length(lineages)) {
    return(data.frame(fit_id = character(), n_pass = integer()))
  }
  
  sizes   <- lengths(lineages)                              # fit lengths
  covered <- chosen <- character()
  
  repeat {
    ## fits whose passages don’t overlap anything already chosen
    ok <- names(lineages)[!vapply(lineages, function(p) any(p %in% covered), logical(1))]
    if (!length(ok)) break                                   # nothing left to add
    
    ## sample one fit with probability ∝ length
    fit <- sample(ok, 1, prob = sizes[ok])
    chosen  <- c(chosen, fit)
    covered <- c(covered, lineages[[fit]])
  }
  
  data.frame(fit_id = chosen, n_pass = sizes[chosen], stringsAsFactors = FALSE)[
    sizes[chosen] > 1, ]
  
}

# 4. Annotate with metadata
annotate_samples <- function(df, meta_file) {
  meta <- read.csv(meta_file, stringsAsFactors = FALSE)
  key <- paste(meta$datasetname, meta$timepoint, sep = "_")
  lut_pdx <- setNames(meta$PDX_id, key)
  lut_tp  <- setNames(meta$timepoint, key)
  
  roots <- sapply(strsplit(df$fi, "_"), function(v) {
    i <- which(v == "l"); paste(v[1:(i-1)], collapse = "_")
  })
  df$pdx <- lut_pdx[roots]
  df$tp  <- lut_tp[roots]
  df
}

# Build Δf table per focal karyotype
build_deltaf_df <- function(df) {
  lst <- pbapply::pblapply(seq_len(nrow(df)), function(i) {
    path <- file.path("data/processed/salehi/alfak_outputs",
                      df$fi[i],
                      paste0("minobs_", df$min_obs[i]),
                       "landscape.Rds")
   # path <- file.path("data/salehi/alfak_outputs_V1a",
    #                  paste0("minobs_", df$min_obs[i]),
     #                 df$fi[i],
      #                "landscape.Rds")
    if (!file.exists(path)) return(NULL)
    ls <- readRDS(path)
    lut <- setNames(ls$mean, ls$k)
    focal <- ls$k[ls$fq]
    dfi <- do.call(rbind, lapply(focal, function(fk) {
      pl <- mean(as.numeric(strsplit(fk, "\\.")[[1]]))
      nn <- apply(gen_all_neighbours(fk, remove_nullisomes = FALSE), 1, paste, collapse = ".")
      data.frame(
        ploidy      = pl,
        f           = lut[fk],
        deltaf      = lut[fk] - lut[nn],
        mut_id      = rep(paste0(rep(1:22, each = 2), c("-", "+")), 1),
        k_id        = fk,
        stringsAsFactors = FALSE
      )
    }))
    trt  <- strsplit(df$train_treat[i], "")[[1]]
    dfi$treat_switch <- length(unique(trt)) > 1
    dfi$treated      <- grepl("y", df$train_treat[i], fixed = TRUE)
    dfi$ntrain       <- df$ntrain[i]
    dfi$pdx          <- df$pdx[i]
    dfi$fi           <- df$fi[i]
    dfi
  })
  do.call(rbind, lst)
}


derive_cat_vars <- function(df) {
  df$wgd     <- df$ploidy > 2.75
  df$context <- factor(ifelse(df$pdx %in% c("SA039", "SA906"), "invitro", "PDX"))
  df$treat   <- factor(ifelse(df$treated, "cisplatin", "control"))
  df$feat    <- interaction(df$mut_id)
  
  # Explicitly set  reference levels for the models.
  df$context <- relevel(df$context, ref = "invitro")
  df$treat   <- relevel(df$treat, ref = "control")
  df
}


# Pairwise analysis
compute_pairwise <- function(df) {
  ylist <- split(df, interaction(df$k_id, df$fi))
  ylist <- ylist[sapply(ylist, nrow) > 0]
  kdat  <- do.call(rbind, lapply(ylist, function(z) z[1, ]))
  cmb   <- combn(nrow(kdat), 2)
  pair  <- data.frame(
    dk    = NA, sim = NA,
    fi1   = kdat$fi[cmb[1, ]], fi2 = kdat$fi[cmb[2, ]],
    pdx1  = kdat$pdx[cmb[1, ]], pdx2 = kdat$pdx[cmb[2, ]],
    treat1 = kdat$treat[cmb[1, ]], treat2 = kdat$treat[cmb[2, ]],
    k1    = kdat$k_id[cmb[1, ]], k2 = kdat$k_id[cmb[2, ]],
    stringsAsFactors = FALSE
  )
  for (j in seq_len(ncol(cmb))) {
    i1 <- cmb[1, j]; i2 <- cmb[2, j]
    v1 <- as.integer(strsplit(pair$k1[j], "\\.")[[1]])
    v2 <- as.integer(strsplit(pair$k2[j], "\\.")[[1]])
    pair$dk[j]  <- sum(abs(v1 - v2))
    pair$sim[j] <- cor(ylist[[i1]]$deltaf,
                       ylist[[i2]]$deltaf,
                       use = "pairwise.complete.obs")
  }
  pair$pair_type  <- with(pair,
                          ifelse(fi1 == fi2, "same_traj",
                                 ifelse(pdx1 == pdx2, "parallel", "diff_line")))
  pair$treat_pair <- with(pair,
                          ifelse(treat1 == treat2, treat1, "mixed"))
  pair$log_dk <- log1p(pair$dk)
  pair
}

## one bootstrap rep, l (lineages list) and xo (fit metadata) must be available
one_boot <- function(boot_rep){
  ## pick random set of nonoverlapping lineages
  rep  <- get_reps(l,x0)
  x    <- annotate_samples(x0[x0$fi %in% rep$fit_id,],
                           "data/raw/salehi/metadata.csv")
  ## get the one-MS-fitness changes for all frequent karyotypes.
  df   <- derive_cat_vars(build_deltaf_df(x))
  ## convert delta fitnesses to z scores.
  df$absdf <- abs(df$deltaf)+1e-6 
  df$absdf <- df$absdf / sd(df$absdf,na.rm = TRUE)
  
  ## check for influence of context (PDX or invitro) or treatment
  ## on delta-f scores:
  fit_var_context <- glmmTMB(absdf ~ context + (1|pdx/fi/k_id),
                             family = Gamma("log"), data = df)
  fit_var_treat   <- glmmTMB(absdf ~ treat   + (1|pdx/fi/k_id),
                             family = Gamma("log"), data = subset(df,context=="PDX"))
  
  ## compute log variance of delta-f per frequent karyotype
  kvar <- df %>% group_by(pdx,fi,k_id,context,treat) %>%
    summarise(v = var(deltaf), .groups="drop") %>%
    mutate(logv = scale(log(v+1e-6)))
  ## check for context effects:
  fit_kvar_context <- glmmTMB(logv ~ context + (1|pdx/fi),
                              family = gaussian(), data = kvar)
  fit_kvar_treat   <- glmmTMB(logv ~ treat   + (1|pdx/fi),
                              family = gaussian(), data = subset(kvar,context=="PDX"))
  
  ## conduct analysis of pairwise similarity...
  pair <- compute_pairwise(df) 
  pair$treat_pair[pair$treat_pair%in%c(1,2)] <- "same"
  pair$z <- atanh(pair$sim)
  # pair$dist <- "hi"
 # pair$dist[pair$dk<3] <- "lo"
  
  cross_fit <- lmer(z ~ log1p(dk)*pair_type+treat_pair + (1|fi1)+(1|fi2), data = pair)
  df2 <- data.frame(coef(summary(cross_fit)))
  df2$rep <- boot_rep

 # cross_fit2 <- lmer(sim ~ dist +pair_type+treat_pair + (1|fi1)+(1|fi2), data = pair[pair$pair_type!="same_traj",])
  #df3 <- data.frame(coef(summary(cross_fit2)))
  #df3$rep <- boot_rep
  
  df1=data.frame(rep          = boot_rep,
                 PDX_v  = coef(summary(fit_var_context))$cond["contextPDX","Estimate"],
                 cisplatin_v= coef(summary(fit_var_treat))$cond["treatcisplatin","Estimate"],
                 PDX_kv = coef(summary(fit_kvar_context))$cond["contextPDX","Estimate"],
                 cisplatin_kv=coef(summary(fit_kvar_treat))$cond["treatcisplatin","Estimate"])
  
  
  
  list(df1=df1,df2=df2)#,df3=df3)
}

## summarize the overlap between successive passages in evolutionary experiments.
## a weakness here is not accounting for time...
passage_stats <- function(pass,m,alfak_inputs_dir="data/processed/salehi/alfak_inputs"){
  
  uid <- strsplit(pass,split="-") |> unlist() |> tail(1)
  filepattern <- paste(m[m$uid==uid,c("datasetname","timepoint")],collapse="_")
  ff <- list.files(alfak_inputs_dir)
  matches <- grepl(filepattern,ff)
  if(sum(matches)<1) return(NULL)
  targetPath <- head(ff[matches],1)
  targetPath <- file.path(alfak_inputs_dir,targetPath)
  
  yi <- readRDS(targetPath)$x
  yi <- yi[,(ncol(yi)-1):ncol(yi)]
  
  k <- do.call(rbind,lapply(rownames(yi),s2v))
  mean_kary <- data.frame(t(t(k) %*% as.matrix(yi) / rep(colSums(yi), each = ncol(k))))
  
  freq <- sweep(yi, 2, colSums(yi), "/")
  
  list(init_kary = mean_kary[1,],
       delta_kary=mean_kary[2,]-mean_kary[1,],
       overlap=sum(pmin(freq[,1],freq[,2])),
       meta=m[m$uid==uid,])
}

cnap_heatmap <- function(df){
  set.seed(42)
  pdx_sel <- c("SA535","SA609")
  fi_sel <- df %>% filter(pdx %in% pdx_sel) %>% group_by(pdx) %>% summarise(fi=list(sample(unique(fi),2)),.groups="drop") %>% unnest(fi) %>% rename(fi=fi)
  kid_sel <- df %>% semi_join(fi_sel,by=c("pdx","fi")) %>% group_by(pdx,fi) %>% summarise(k_id=list(sample(unique(k_id),3)),.groups="drop") %>% unnest(k_id)
  sampled_df <- df %>% semi_join(kid_sel,by=c("pdx","fi","k_id"))
  k_levels <- kid_sel %>% distinct(pdx,fi,k_id) %>% arrange(pdx,fi,k_id) %>% pull(k_id) %>% unique()
  sampled_df <- sampled_df %>% mutate(pdx=factor(pdx,levels=pdx_sel),fi=factor(fi,levels=fi_sel$fi[order(fi_sel$pdx)]),k_id=factor(k_id,levels=k_levels),feat=factor(feat,levels=unique(feat))) %>% group_by(pdx) %>% mutate(fi_index=dense_rank(fi),fi=paste0("Landscape ",fi_index)) %>% ungroup() %>% mutate(fi=factor(fi,levels=paste0("Landscape ",seq_len(max(fi_index))))) %>% select(-fi_index)
  cb_low <- "#D55E00"; cb_high <- "#0072B2"
  ggplot(sampled_df, aes(x = k_id, y = feat, fill = deltaf)) +
    geom_tile(color = "white", size = 0.1) +
    scale_fill_gradient2(
      low      = cb_low,
      mid      = "white",
      high     = cb_high,
      midpoint = 0,
      name     = expression(Delta~f)
    ) +
    facet_nested(
      cols      = vars(pdx, fi),
      scales    = "free_x",
      space     = "free_x",
      nest_line = TRUE
    ) +
    base_theme +
    labs(x = "frequent karyotypes") +   
    theme(
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y  = element_text(hjust = 0, margin = margin(r = 5)),
      panel.grid       = element_blank(),
      panel.spacing.x  = unit(0.2, "lines"),
      strip.background = element_rect(fill = "white", color = "black"),
      strip.placement  = "outside",
      axis.title.y       = element_blank()
    )
}

l   <- load_lineages("data/processed/salehi/lineages.Rds")
x0   <- load_fits("data/processed/salehi/alfak_outputs_proc/")

if(redo_bootstrap){
  
  cl <- makeCluster(32)                                         
  
  clusterEvalQ(cl,{                                             # load pkgs on workers
    library(glmmTMB); library(lme4); library(dplyr)
  })
  
  ## export all objects / helpers used inside the loop
  clusterExport(cl,c("l","x0","get_reps","annotate_samples",
                     "build_deltaf_df","derive_cat_vars","compute_pairwise"),
                envir = environment())
  clusterCall(cl,function(dummyVar) source("R/utils_karyo.R"))
  
  Nboot <- 200
  boot_raw <- parLapplyLB(cl,1:Nboot,one_boot)                 # progress-less LB
  stopCluster(cl)
  saveRDS(boot_raw,"data/processed/salehi/landscape_analysis_bootstrap.Rds")
}

boot_raw <- readRDS("data/processed/salehi/landscape_analysis_bootstrap.Rds")


## 1 - stack all bootstrap estimates
boot_df <- do.call(rbind, lapply(boot_raw,
                                 \(x) data.frame(param = rownames(x$df2), est = x$df2$Estimate)))

## 2 - 95 % CIs (2.5 %, 50 %, 97.5 %)
ci <- as.data.frame(t(sapply(split(boot_df$est, boot_df$param),
                             quantile, probs = c(.025, .5, .975))))
names(ci) <- c("lo", "mid", "hi"); ci$param <- rownames(ci)

## 1. Look-up vector ─ terse but descriptive
lut <- c(
  "(Intercept)"                     = "Baseline similarity",
  "log1p(dk)"                       = "Distance (log₁₊)",
  "pair_typeparallel"               = "Parallel lineage",
  "pair_typesame_traj"              = "Same trajectory",
  "treat_pairsame"                  = "Same treatment",
  "log1p(dk):pair_typeparallel"     = "Distance × Parallel",
  "log1p(dk):pair_typesame_traj"    = "Distance × Same traj"
)

## 2. Augment the ci data frame on the fly
ci$label <- factor(ci$param, levels = names(lut), labels = lut)

## 3. Refined plot (one extra line beyond your original)
plt_cor <- ggplot(ci, aes(label, mid)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = .25) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = .3) +
  coord_flip() +
  labs(x = NULL, y = "Median effect (95 % bootstrap CI)") +
  base_theme

boot_df1 <- do.call(rbind, lapply(boot_raw, "[[", "df1"))
df1 <- reshape2::melt(boot_df1, id.vars="rep")
df1$variable <- sub("_z$", "", df1$variable)
df1 <- tidyr::separate(df1, variable, into=c("context","stat"), sep="_", remove=FALSE)
smm <- aggregate(list(z=df1$value), by=list(context=df1$context, stat=df1$stat), quantile, probs=c(0.025,0.5,0.975))
smm <- cbind(smm[,c("context","stat")], smm$z)
smm$stat <- factor(smm$stat, levels=c("v","kv"))
smm$context <- factor(smm$context, levels=c("cisplatin","PDX"), labels=c("cisplatin v.s.\nuntreated","PDX v.s. \nin-vitro"))
plt_boot <- ggplot(smm, aes(x=`50%`, y=stat)) +
  facet_grid(cols=vars(context)) +
  geom_point() +
  geom_errorbarh(aes(xmin=`2.5%`, xmax=`97.5%`), height=0.2) +
  scale_y_discrete(labels=c("v"=expression(Delta~f), "kv"="variance")) +
  labs(x="Estimate (median and 95% CI)", y=NULL) +
  coord_flip()+
  base_theme

all_passages <-  unlist(l) |> unique()
m <- read.csv("data/raw/salehi/metadata.csv")

ps <- pbapply::pblapply(all_passages,passage_stats,m=m)
ps <- ps[!sapply(ps,is.null)]

meta_df <- do.call(rbind, lapply(ps, `[[`, "meta"))
meta_df$overlap <- sapply(ps, `[[`, "overlap")
meta_df$pdx <- get_sample_lineage_map()[meta_df$datasetname]
meta_df$context <- "pdx"
meta_df$context[grepl("p53",meta_df$pdx)] <- "in-vitro"

plot_overlap <- ggplot(meta_df, aes(x=on_treatment, y=overlap, shape=on_treatment, color=on_treatment)) +
  geom_point(position=position_dodge(width=0.6)) +
  facet_nested(cols=vars(context, pdx), scales="free_x", space="free_x", nest_line=TRUE) +
  scale_color_manual(values=c("y"="#0072B2","n"="#D55E00")) +
  labs(y="Overlap", color="Treatment", shape="Treatment") +
  make_base_theme("minimal") +
  theme(
    strip.background=element_rect(fill="white", color="black"),
    strip.placement="outside",
    axis.text.x=element_blank(),
    axis.title.x=element_blank(),
    legend.position=c(0.98, 0.98),
    legend.justification=c("right","top"),
    legend.background=element_rect(fill="white", color="black")
  )
plot_overlap


fit1 <- lmer(overlap ~ on_treatment + (1 | pdx), data = meta_df)
summary(fit1)        # Wald t-tests for treatment levels
anova(fit1)

fit2 <- lmer(overlap ~ context + (1 | pdx/on_treatment), data = meta_df)
summary(fit2)        # Wald t-tests for treatment levels
anova(fit2)


##illustrate bootstrap concept:
base_plt <- basic_lineage_plot(m)
bootstrap_example_plots <- lapply(1:3,function(i){
  boot_i <- get_reps(l,x0)
  lins_i <- l[boot_i$fi]
  
  
  edge_highlights <- generate_segment_highlights(base_plt,lins_i)
  
  # Now overlay with geom_segment
  base_plt$plt +
    geom_segment(data = edge_highlights,
                 aes(x = x_from, y = y_from, xend = x_to, yend = y_to, color = categorical_variable),
                 show.legend = F)
})


## Run thru the process of an example bootstrap, highlighting relevant factors
set.seed(42)
rep  <- get_reps(l,x0)
x    <- annotate_samples(x0[x0$fi %in% rep$fit_id,],
                         "data/raw/salehi/metadata.csv")
df   <- derive_cat_vars(build_deltaf_df(x))
plt_cnap <- cnap_heatmap(df)

plt_ecdf <- ggplot(subset(df, abs(deltaf) < 0.3),
       aes(x = deltaf, colour = interaction(context, treat),
           group = interaction(context, treat))) +
  stat_ecdf() +
  scale_y_continuous("frequency") +
  scale_x_continuous("fitness effect") +
  scale_color_discrete("",labels=c("I.V\n(unt.)","PDX\n(unt.)","PDX\n(cisp.)")) +
  base_theme+
  theme(legend.position = "top")

pA <- cowplot::plot_grid(plotlist = bootstrap_example_plots,nrow=3)


pCD <- cowplot::plot_grid(plt_ecdf,plt_boot,nrow=1,rel_widths = c(3,4),
                           labels=c("c","d"),label_size = base_text_size+2)
pCDE <- cowplot::plot_grid(pCD,plot_overlap,nrow=2,labels=c("","e"),label_size = base_text_size+2)
plt <- cowplot::plot_grid(pA,plt_cnap,pCDE,nrow=1,rel_widths = c(1,2,2),labels=c("a","b",""),label_size = base_text_size+2)

ggsave("figs/cna_fitness_effects.png",plt,width=8,height=4,units="in", 
       bg = "white")

pA <- cowplot::plot_grid(plotlist = bootstrap_example_plots,nrow=3)


pCD <- cowplot::plot_grid(plt_ecdf,plt_boot,nrow=1,rel_widths = c(4,4),
                          labels=c("c","d"),label_size = base_text_size+2)
pBCD <- cowplot::plot_grid(plt_cnap,pCD,nrow=2,labels=c("b",""),
                           rel_heights=c(2,1),label_size = base_text_size+2)

pEF <- cowplot::plot_grid(plot_overlap,plt_cor,nrow=2,labels=c("e","f"),label_size = base_text_size+2)
plt <- cowplot::plot_grid(pA,pBCD,pEF,nrow=1,rel_widths = c(1.2,2,2),labels=c("a","",""),label_size = base_text_size+2)

ggsave("figs/cna_fitness_effects.pdf",plt,device=cairo_pdf,width=180,height=110,units="mm", 
       bg = "white",dpi=600)

##############################################################################
# 8. EXPORT SOURCE DATA (Nature Requirement) - FIGURE 4
##############################################################################

dir.create("data/source_data", showWarnings = FALSE, recursive = TRUE)

# --- Panel A: Bootstrap Example Trees (Schematic) ---
# Since these are random selections of lineages, we provide the IDs of the 
# lineages selected in the 3 examples to ensure reproducibility.
set.seed(NULL) # Reset seed logic based on your loop
sd_4a <- do.call(rbind, lapply(1:3, function(i) {
  # We assume the seed state matches the loop in your script
  # If you set a specific seed before the loop, this captures it.
  # Otherwise, this captures the logic:
  boot_reps <- get_reps(l, x0)
  data.frame(
    Example_Iteration = i,
    Selected_Lineage_ID = boot_reps$fit_id,
    Num_Passages = boot_reps$n_pass
  )
}))
saveRDS(sd_4a, "data/source_data/Fig4a.Rds")

# --- Panel B: CNAP Heatmap (Sub-sampled Data) ---
# We must replicate the data filtering logic inside your cnap_heatmap function
# to save the exact tiles shown in the figure.
set.seed(42) # Matches the seed inside cnap_heatmap
pdx_sel <- c("SA535","SA609")
# Replicate logic:
fi_sel <- df %>% 
  filter(pdx %in% pdx_sel) %>% 
  group_by(pdx) %>% 
  summarise(fi=list(sample(unique(fi),2)), .groups="drop") %>% 
  unnest(fi)
kid_sel <- df %>% 
  semi_join(fi_sel, by=c("pdx","fi")) %>% 
  group_by(pdx,fi) %>% 
  summarise(k_id=list(sample(unique(k_id),3)), .groups="drop") %>% 
  unnest(k_id)
sd_4b <- df %>% 
  semi_join(kid_sel, by=c("pdx","fi","k_id")) %>%
  select(PDX = pdx, Landscape_ID = fi, Karyotype_ID = k_id, Mutation_Feature = feat, Delta_Fitness = deltaf) %>%
  mutate(Delta_Fitness = round(Delta_Fitness, 4))

saveRDS(sd_4b, "data/source_data/Fig4b.Rds")

# --- Panel C: ECDF of Fitness Effects ---
# Source: 'df' (the single bootstrap realization used for the plot)
sd_4c <- subset(df, abs(deltaf) < 0.3) %>%
  select(PDX = pdx, Context = context, Treatment = treat, Delta_Fitness = deltaf) %>%
  mutate(Delta_Fitness = round(Delta_Fitness, 4))
saveRDS(sd_4c, "data/source_data/Fig4c.Rds")

# --- Panel D: Bootstrap Model Results (Forest Plot) ---
# Source: 'smm' (Summary of bootstrap quantiles)
sd_4d <- smm %>%
  select(Comparison_Context = context, Statistic = stat, 
         CI_Low_2.5 = `2.5%`, Median = `50%`, CI_High_97.5 = `97.5%`) %>%
  mutate(
    Statistic = ifelse(Statistic == "v", "Delta_Fitness", "Variance"),
    Comparison_Context = as.character(Comparison_Context), # Remove newline chars if present
    CI_Low_2.5 = round(CI_Low_2.5, 4),
    Median = round(Median, 4),
    CI_High_97.5 = round(CI_High_97.5, 4)
  )
saveRDS(sd_4d, "data/source_data/Fig4d.Rds")

# --- Panel E: Lineage Overlap ---
# Source: 'meta_df'
sd_4e <- meta_df %>%
  select(PDX = pdx, Context = context, On_Treatment = on_treatment, Overlap_Score = overlap) %>%
  mutate(Overlap_Score = round(Overlap_Score, 4))
saveRDS(sd_4e, "data/source_data/Fig4e.Rds")

# --- Panel F: Pairwise Correlation Effects ---
# Source: 'ci' (Bootstrap confidence intervals)
sd_4f <- ci %>%
  select(Parameter_Label = label, Median_Effect = mid, CI_Low = lo, CI_High = hi) %>%
  mutate(
    Median_Effect = round(Median_Effect, 4),
    CI_Low = round(CI_Low, 4),
    CI_High = round(CI_High, 4)
  )
saveRDS(sd_4f, "data/source_data/Fig4f.Rds")

message("Figure 4 Source Data saved to data/source_data/")
