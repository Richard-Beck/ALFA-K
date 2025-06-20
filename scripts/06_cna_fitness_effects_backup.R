if(!basename(getwd())=="ALFA-K") stop("Ensure working directory set to ALFA-K root")


source("R/utils_karyo.R")
source("R/utils_theme.R")
source("R/utils_env.R")
libs <- c("igraph", "lme4", "glmmTMB", "car", "ggplot2", "pbapply", "tidygraph", 
          "ggraph","ggsignif","isotone","dplyr","ggnewscale","emmeans")
ensure_packages(libs)

base_text_size <- 5
# 1. Load and process lineages
load_lineages <- function(path) {
  raw <- readRDS(path)
  lapply(raw, function(li) paste0(head(li$ids, -1), "-", tail(li$ids, -1)))
}

# 2. Load, filter, and collapse fits
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





# 3. Compute overlap components and select reps
get_reps <- function(lineages, df) {
  lineages <- lineages[names(lineages) %in% df$fi]
  fits <- names(lineages)
  all_p <- sort(unique(unlist(lineages)))
  mat <- matrix(FALSE, length(fits), length(all_p),
                dimnames = list(fits, all_p))
  for (i in seq_along(lineages)) mat[i, lineages[[i]]] <- TRUE
  adj <- (mat %*% t(mat)) > 0
  diag(adj) <- FALSE
  comps <- components(graph_from_adjacency_matrix(adj))$membership
  res <- data.frame(fit_id = fits, overlap_id = comps,
                    n_pass = lengths(lineages), stringsAsFactors = FALSE)
  reps <- by(seq_len(nrow(res)), res$overlap_id, function(ids) ids[which.max(res$n_pass[ids])])
  res[unlist(reps), ]
}

get_reps <- function(lineages, df) {
  # 1) restrict to fits present in df$fi
  lineages <- lineages[names(lineages) %in% df$fi]
  if (length(lineages) == 0) {
    return(data.frame(fit_id = character(0), n_pass = integer(0), stringsAsFactors = FALSE))
  }
  
  # 2) compute size (# of passages) for each fit
  sizes <- lengths(lineages)  # named integer vector
  
  # 3) order fits by descending size
  sorted_fits <- names(sort(sizes, decreasing = TRUE))
  
  # 4) greedily pick non-overlapping fits
  covered <- character(0)
  chosen  <- character(0)
  for (fit in sorted_fits) {
    passages <- lineages[[fit]]
    if (length(intersect(passages, covered)) == 0) {
      chosen  <- c(chosen, fit)
      covered <- c(covered, passages)
    }
  }
  
  # 5) return data.frame of chosen fits and their sizes
  res <- data.frame(
    fit_id = chosen,
    n_pass = sizes[chosen],
    stringsAsFactors = FALSE
  )
  res[res$n_pass>1,]
}

get_reps <- function(lineages,df){
  ## keep only fits present in df$fi
  lineages <- lineages[names(lineages) %in% df$fi]
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

# 5. Build Δf table per focal karyotype
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

# 6. Derive categorical vars
derive_cat_vars <- function(df) {
  df$wgd     <- df$ploidy > 2.75
  df$context <- factor(ifelse(df$pdx %in% c("SA039", "SA906"), "invitro", "PDX"))
  df$treat   <- factor(ifelse(df$treated, "cisplatin", "control"))
  df$feat    <- interaction(df$mut_id)
  df
}

# 9. Fit mixed model
fit_glmm <- function(df,fitvar) {
  df$deltaf <- as.numeric(scale(df$deltaf))
  res <- NULL
  if(fitvar=="context"){
    res <- glmmTMB(deltaf ~ context + (1 | pdx/fi/k_id),
                   #dispformula = ~ context + treat,
                   data = df, family = gaussian(),
                   control = glmmTMBControl(optArgs = list(
                     iter.max = 2e4, eval.max = 2e4
                   )))
  }
  if(fitvar=="treat"){
    res <- glmmTMB(deltaf ~ treat + (1 | pdx/fi/k_id),
                   #dispformula = ~ context + treat,
                   data = df, family = gaussian(),
                   control = glmmTMBControl(optArgs = list(
                     iter.max = 2e4, eval.max = 2e4
                   )))
  }
  return(res)
}

# 10. ECDF plot
plot_ecdf <- function(df, theme) {
  ggplot(subset(df, abs(deltaf) < 0.3),
         aes(x = deltaf, colour = interaction(context, treat),
             group = interaction(context, treat))) +
    stat_ecdf() +
    scale_y_continuous("frequency") +
    scale_x_continuous("fitness effect") +
    scale_color_discrete("condition") +
    theme_bw() + theme
}

# 12. Pairwise analysis
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

# 14. Cross-trajectory model
fit_cross <- function(pair) {
  lmer(sim ~ log_dk * pair_type + treat_pair + (1 | fi1) + (1 | fi2), data = pair)
}

# 15. Plot helpers for decay & cross
sig_code <- function(p){
  ifelse(p < .001, "***", ifelse(p < .01, "**",
                                 ifelse(p < .05,  "*",  "n.s.")))
}

# 2. Boxplot of Δf distributions by context & treatment
plot_deltaf_box <- function(df, theme) {
  ggplot(df, aes(x = context, y = deltaf, fill = treat)) +
    geom_boxplot(outlier.size = 1, position = position_dodge(0.8)) +
    ylab("Fitness effect (Δf)") + xlab("Context") +
    scale_fill_brewer("Treatment", palette = "Set2") +
    theme_bw() + theme
}

plot_deltaf_var <- function(kvar, fit_kvar_context, fit_kvar_treat){
  pvals <- c(coef(summary(fit_kvar_context))$cond["contextPDX","Pr(>|z|)"],
             coef(summary(fit_kvar_treat))$cond["treatcontrol","Pr(>|z|)"])
  
  pvals <- format.pval(pvals, digits = 2, eps = 0.001)
  
  # Your base plot
  p <- ggplot(kvar, aes(x = interaction(treat, context,sep = "\n"), y = logv)) +
    geom_jitter(width = 0.2, height = 0, alpha = 0.6) +
    theme_minimal()+ 
    geom_signif(
      annotations = c(paste0("p = ",pvals[1]), paste0("p = ",pvals[2])),  
      y_position = c(2.3, 2.0),                   
      xmin = c(1, 2),                             
      xmax = c(2.5, 3),                             
      tip_length = 0.02,
      textsize = 3.5
    )+
    scale_x_discrete("")+
    scale_y_continuous("log var.", expand = c(0, 1.15))
  p
}

plot_absdf_var <- function(df, fit_var_context, fit_var_treat) {
  pvals <- c(coef(summary(fit_var_context))$cond["contextPDX","Pr(>|z|)"],
             coef(summary(fit_var_treat))$cond["treatcontrol","Pr(>|z|)"])
  pvals <- format.pval(pvals, digits=2, eps=0.001)
  
  ggplot(df, aes(x=interaction(treat, context, sep="\n"), y=absdf)) +
    geom_jitter(width=0.2, height=0, alpha=0.6) +
    geom_signif(annotations=c(paste0("p = ", pvals[1]), paste0("p = ", pvals[2])),
                y_position = c(2.5, 2.0),
                xmin=c(1,2), xmax=c(2.5,3), tip_length=0.02, textsize=3.5) +
    scale_color_brewer(name="Context", palette="Set2") +
    scale_shape_discrete(name="Treatment") +
    scale_y_log10("log₁₀(abs Δf)", expand = c(0, 1.15)) +
    scale_x_discrete("")+
    theme_minimal()
}

## isotonic regression on corr-dist slope.
get_iso_slope <- function(dat) {
  if (nrow(dat) < 2) return(NA_real_)
  o   <- order(dat$dk)
  x   <- dat$dk[o]
  y   <- dat$sim[o]
  
  iso <- stats::isoreg(x, -y)       # isoreg enforces *increasing*; −y ⇒ decreasing
  yhat <- -iso$yf
  
  # first two distinct distance values
  uniq <- which(!duplicated(x))
  if (length(uniq) < 2) return(NA_real_)
  (yhat[uniq[2]] - yhat[uniq[1]]) / (x[uniq[2]] - x[uniq[1]])
}
# 16. Main orchestration
main <- function() {
  
  theme <- make_base_theme()
  
  l   <- load_lineages("data/processed/salehi/lineages.Rds")
  x   <- load_fits("data/processed/salehi/alfak_outputs_proc/")
  rep <- get_reps(l, x)
  x   <- x[x$fi %in% rep$fit_id, ]
  x   <- annotate_samples(x, "data/raw/salehi/metadata.csv")
  
  df  <- build_deltaf_df(x)
  df  <- derive_cat_vars(df)
  
  ##explicit variance model:
  
  df$absdf <- abs(df$deltaf) + 1e-6              # protect against exact zeros
  df$absdf <- df$absdf / sd(df$absdf,na.rm=T)           # scale for model stability
  
  fit_var_context <- glmmTMB(
    absdf ~ context  + (1 | pdx/fi/k_id),  # same hierarchy
    family  = Gamma(link = "log"),                # log-Gamma model for positive values
    data    = df
  )
  
  print(summary(fit_var_context))
  print(car::Anova(fit_var_context, type = "III"))
  print(confint(fit_var_context, parm = "beta_", level = 0.95))
  
  fit_var_treat <- glmmTMB(
    absdf ~ treat  + (1 | pdx/fi/k_id),  # same hierarchy
    family  = Gamma(link = "log"),                # log-Gamma model for positive values
    data    = subset(df,context=="PDX")
  )
  
  print(summary(fit_var_treat))
  print(car::Anova(fit_var_treat, type = "III"))
  print(confint(fit_var_treat, parm = "beta_", level = 0.95))
  
  kvar <- df %>% 
    group_by(pdx, fi, k_id, context, treat) %>% 
    summarise(v = var(deltaf), .groups = "drop")
  
  kvar$logv <- log(kvar$v + 1e-6)
  kvar$logv <- as.numeric(scale(kvar$logv))
  
  fit_kvar_context <- glmmTMB(
    logv ~ context + (1 | pdx/fi),
    family = gaussian(),
    data   = kvar
  )
  
  print(summary(fit_kvar_context))
  
  fit_kvar_treat <- glmmTMB(
    logv ~ treat + (1 | pdx/fi),
    family = gaussian(),
    data   = subset(kvar,context=="PDX")
  )
  print(summary(fit_kvar_treat))

  ## these give warnings:
  performance::check_residuals(fit_var_context)
  performance::check_residuals(fit_var_treat)
  
  performance::check_overdispersion(fit_var_context)
  performance::check_overdispersion(fit_var_treat)
  
  p_var <- plot_absdf_var(df,fit_var_context,fit_var_treat)
  p_kvar <- plot_deltaf_var(kvar,fit_kvar_context,fit_kvar_treat)
  p_ecdf <- plot_ecdf(df, theme)
  
  
  # ------------------------------------------------------------------
  # 1.  Pairwise table and Fisher-z transform          (same as before)
  # ------------------------------------------------------------------
  pair      <- compute_pairwise(df)
  pair      <- subset(pair, dk < 5)
  pair <- pair[pair$pair_type%in%c("paralell","same_traj"),]
  pair$z <- atanh(pmax(pmin(pair$sim,  0.999), -0.999))
  if(FALSE){
    m0 <- lmer(z ~ log1p(dk)                       + (1|fi1) + (1|fi2), data = pair)
    m1 <- lmer(z ~ log1p(dk) + pair_type + treat_pair              + (1|fi1) + (1|fi2), data = pair)
    m2 <- lmer(z ~ log1p(dk) + pair_type*treat_pair                + (1|fi1) + (1|fi2), data = pair)
    m3 <- lmer(z ~ log1p(dk)*pair_type*treat_pair                  + (1|fi1) + (1|fi2), data = pair)
    anova(m0, m1, m2, m3)   # step‑wise LRTs with valid df
    emtrends(m3, pairwise ~ pair_type*treat_pair, var = "log1p(dk)")
    emmeans(m3, ~ pair_type*treat_pair)
    
    same  <- subset(pair, pair_type == "same_traj")
    pair$cisplatin <- "no"
    pair$cisplatin[pair$treat_pair!="2"] <- "yes"
    
    same$z <- atanh(pmax(pmin(same$sim,  0.999), -0.999))
    
    
    cross_fit <- lmer(z~log1p(dk)*pair_type:treat_pair+ (1 | fi1) + (1 | fi2),
                      data=pair)
    car::Anova(cross_fit,type=3)
    
    emm <- emmeans(cross_fit,
                   ~ treat_pair | pair_type)
    plot(emm)
    emm <- emmeans(cross_fit, ~ z| treat_pair*pair_type)
    plot(emm)
    # Pairwise treatment comparisons within each pair_type
    pairs(emm)
    plot(emm)
  }
  
  
  # ------------------------------------------------------------------
  
  # C. Raw similarity vs distance by pair_type with LOESS
  renamr <- c(diff_line="diff. line",parallel="parallel",same_traj="same traj.")
  pC <- ggplot(pair, aes(x = dk, y = sim)) +
    facet_grid(rows=vars(renamr[pair_type]))+
    geom_point(alpha = 0.15, size = 0.5) +
    stat_summary(fun = mean, geom = "point", shape = 95, size = 5) +
    labs(x = "Manhattan distance", y = "Pearson similarity", color = "Pair type") +
    scale_color_brewer(palette = "Dark2") +
    theme_minimal() + theme
  pC
  

  list(plots=list(p_ecdf=p_ecdf,p_kvar=p_kvar,p_var=p_var,pC=pC),
       lineages=rep)
}

# Execute
main_output <- main()
plotList <- main_output$plots
##if theres a kvar context error run line by line (scope issue?)
#5.1 Get mean karyotype per timepoint/trajectory
get_centroids <- function(x){
  do.call(rbind,lapply(1:nrow(x),function(i){
    l <- readRDS(paste0("data/processed/salehi/alfak_inputs/",x$fi[i],".Rds"))$x
    k <- do.call(rbind,lapply(rownames(l),s2v))
    mean_kary <- data.frame(t(t(k) %*% as.matrix(l) / rep(colSums(l), each = ncol(k))))
    mean_kary$time <- colnames(l)
    mean_kary$fi <- x$fi[i]
    rownames(mean_kary) <- NULL
    mean_kary
  }))
}

# 1. Define a helper to compute “abundance overlap” between successive passages
get_overlap <- function(meta_df){
  do.call(rbind, lapply(seq_len(nrow(meta_df)), function(i){
    # load clone × time count matrix
    mat  <- readRDS(
      sprintf("data/processed/salehi/alfak_inputs/%s.Rds", meta_df$fi[i])
    )$x
    # convert to relative frequencies
    freq <- sweep(mat, 2, colSums(mat), "/")
    times <- colnames(mat)
    # compute sum(minimum) for each consecutive pair
    sim <- sapply(seq_len(ncol(freq)-1), function(j){
      sum(pmin(freq[, j], freq[, j+1]))
    })
    data.frame(
      fi    = meta_df$fi[i],
      time  = times[-1],
      sim   = sim,
      stringsAsFactors = FALSE
    )
  }))
}

getVectors <- function(cm){
  cm <- split(cm,f=cm$fi)
  do.call(rbind,lapply(cm,function(cmi){
    tmp <- cmi[2:nrow(cmi),1:22]-cmi[1:(nrow(cmi)-1),1:22]
    tmp$time <- cmi$time[-1]
    tmp$fi <- cmi$fi[-1]
    tmp
  }))
}
mag <- function(v) sqrt(sum(v^2))
getangle <- function(a,b) 180*acos(sum(a*b)/(mag(a)*mag(b)))/pi

dSphereAngle <- function(theta, N) {
  coef <- integrate(function(t) sin(t)^(N-2), 0, pi)$value
  sin(theta)^(N-2)/coef
}
cSphereAngle <- function(theta, N) {
  if(!is.finite(theta)) return(1)
  integrate(function(t) dSphereAngle(t, N), 0, theta)$value
}

xmeta <- function(x,linpath="data/processed/salehi/lineages.Rds"){
  l <- readRDS(linpath)
  
  do.call(rbind,lapply(1:nrow(x),function(i){
    uid <- l[[x$fi[i]]]$ids
    merge(data.frame(uid),x[i,])    
  }))
}

text_size_theme <- make_base_theme()

l   <- load_lineages("data/processed/salehi/lineages.Rds")
x   <- load_fits("data/processed/salehi/alfak_outputs_proc/")
rep <- main_output$lineages

x   <- x[x$fi %in% rep$fit_id, ]
x   <- annotate_samples(x, "data/raw/salehi/metadata.csv")
xm <- xmeta(x)

lut_pdx <- x$pdx
names(lut_pdx) <- x$fi
lut_pdx <- gsub("SA906","P53ko",lut_pdx)
lut_pdx <- gsub("SA039","P53wt",lut_pdx)
lut_traj <- split(lut_pdx,f=lut_pdx)
names(lut_traj) <- NULL
lut_traj <- unlist(lapply(lut_traj,function(li) {
  tmp <- paste0(li,"-",LETTERS[1:length(li)])
  names(tmp) <- names(li)
  tmp
  }))


cm <- get_centroids(x)
cmv <- getVectors(cm)
cmv <- merge(cmv,x[,c("fi","train_treat","pdx")])
cmv$d <- apply(cmv[,paste0("X",1:22)],1,function(ci) sqrt(sum(ci^2)))
cmv$treat <- grepl("y",cmv$train_treat)
cmv$traj <- lut_traj[cmv$fi]

ggplot(cmv, aes(x=d, y=traj, color=pdx,shape=treat))+
  geom_jitter(height=0)+
  scale_color_brewer(name="PDX", palette="Set2")+
  theme_minimal()+
  scale_y_discrete()+
  scale_x_continuous("passage_distance")


# 2. Run it on your filtered 'x' (contains fi, pdx, train_treat, etc)
sim_df <- get_overlap(x)
# 3. Annotate with PDX & treatment
sim_df <- merge(sim_df,
                x[, c("fi", "pdx", "train_treat")],
                by = "fi", all.x = TRUE)
sim_df$treat <- grepl("y", sim_df$train_treat)
sim_df$traj <- lut_traj[sim_df$fi]
p_sim <- ggplot(sim_df, aes(x=sim, y=traj, color=pdx,shape=treat))+
  geom_jitter(height=0)+
  scale_color_brewer(name="cell line", palette="Set2")+
  theme_minimal()+
  scale_shape_discrete("",labels=c("control","cisplatin"))+
  scale_y_discrete("")+
  scale_x_continuous("overlap coef.")+
  text_size_theme
p_sim

cmv <- getVectors(cm)

amat <- sapply(1:nrow(cmv),function(i) sapply(1:nrow(cmv),function(j){
  getangle(cmv[i,1:22],cmv[j,1:22])
}))
cmv <- merge(cmv,x[,c("fi","pdx","tp","train_treat")])
cmv$treat <- grepl("y",cmv$train_treat)

colnames(amat) <- paste(cmv$fi,cmv$time,sep="-")
rownames(amat) <- paste(cmv$fi,cmv$time,sep="-")

amat[is.na(amat)] <- 0

ycoords <- order(cmv$pdx,cmv$treat)
amat <- amat[ycoords,ycoords]
lut_coord <- 1:ncol(amat)
names(lut_coord) <- colnames(amat)

names(ycoords) <- paste0(cmv$fi,cmv$time,sep="-")
a <- reshape2::melt(amat)
colnames(a)[1:2] <- c("id1","i")
#a <- merge(a,x[,c("fi","pdx","tp")])
a$xcoord <- lut_coord[a$id1]
a$ycoord <- lut_coord[a$i]


nulldf <- data.frame(angle=0:180)
nulldf$rads <- pi * nulldf$angle/180
nulldf$CDF  <- sapply(nulldf$rads, cSphereAngle, N=22)

a$pval <- nulldf$CDF[floor(a$value)+1]
a$pval <- p.adjust(a$pval)

lut_treat <- grepl("y",x$train_treat)
names(lut_treat) <- x$fi



a$pdx1 <- lut_pdx[sapply(as.character(a$id1),function(ii) strsplit(ii,split="-") |> unlist() |> head(1))]
a$pdx2 <- lut_pdx[sapply(as.character(a$i),function(ii) strsplit(ii,split="-") |> unlist() |> head(1))]
a$treat1 <- lut_treat[sapply(as.character(a$id1),function(ii) strsplit(ii,split="-") |> unlist() |> head(1))]
a$treat2 <- lut_treat[sapply(as.character(a$i),function(ii) strsplit(ii,split="-") |> unlist() |> head(1))]


dtop <- a[a$xcoord==a$ycoord,]
dtop$ycoord <- max(a$ycoord+1)
dside <- a[a$xcoord==a$ycoord,]
dside$xcoord <- min(a$ycoord-1)

ui <- unique(a[, c("pdx1","xcoord")])

# for each pdx1, find the max xcoord
pdxs     <- unique(ui$pdx1)
maxs     <- sapply(pdxs, function(p) max(ui$xcoord[ui$pdx1==p]))
# sort them so they’re in plotting order
breaks   <- sort(maxs) + 0.5   # +0.5 puts the line *between* the integer coords

ui_pt <- a[a$id1==a$i,]
ui_pt <- ui_pt[c(FALSE,diff(ui_pt$treat1)>0),]
treat_breaks <- ui_pt$xcoord-0.5


perm_test_group_diff <- function(a, value_col = "pval", treat1 = "treat1", treat2 = "treat2", 
                                 pdx1 = "pdx1", pdx2 = "pdx2", id1 = "id1", id2 = "i",
                                 n_perm = 1000, seed = 1) {
  set.seed(seed)
  
  # Helper to assign treatments per PDX
  permute_treatments <- function(df) {
    df <- df
    pdxs <- unique(c(df[[pdx1]], df[[pdx2]]))
    for (p in pdxs) {
      trajs <- unique(c(df[[id1]][df[[pdx1]] == p], df[[id2]][df[[pdx2]] == p]))
      if (length(trajs) <= 1) next
      lab <- sample(c(TRUE, FALSE), length(trajs), replace = TRUE)
      names(lab) <- trajs
      df[[treat1]][df[[pdx1]] == p] <- lab[df[[id1]][df[[pdx1]] == p]]
      df[[treat2]][df[[pdx2]] == p] <- lab[df[[id2]][df[[pdx2]] == p]]
    }
    df
  }
  
  # Compute test statistic: mean between - mean within
  compute_stat <- function(df) {
    same <- df[[treat1]] == df[[treat2]]
    mean(df[[value_col]][!same], na.rm = TRUE) - mean(df[[value_col]][same], na.rm = TRUE)
  }
  
  # Observed statistic
  obs_stat <- compute_stat(a)
  
  # Permutations
  perm_stats <- replicate(n_perm, {
    a_perm <- permute_treatments(a)
    compute_stat(a_perm)
  })
  
  # Empirical p-value (greater or equal)
  perm_stats <- perm_stats[!is.na(perm_stats)]
  p_empirical <- (1 + sum(perm_stats >= obs_stat)) / (1 + length(perm_stats))
  
  
  list(
    observed_statistic = obs_stat,
    permuted_statistics = perm_stats,
    p_empirical = p_empirical
  )
}

tmp <- perm_test_group_diff(a[a$pdx1==a$pdx2&a$xcoord<a$ycoord,])
print(tmp$p_empirical)

perm_test_within_vs_across_pdx <- function(a, value_col = "value", n_perm = 1000, seed = 1) {
  set.seed(seed)
  a <- a[a$xcoord < a$ycoord, ]
  same_pdx <- a$pdx1 == a$pdx2
  
  obs_stat <- mean(a[[value_col]][!same_pdx], na.rm = TRUE) - 
    mean(a[[value_col]][same_pdx], na.rm = TRUE)
  
  perm_stats <- replicate(n_perm, {
    perm <- sample(same_pdx)
    mean(a[[value_col]][!perm], na.rm = TRUE) - 
      mean(a[[value_col]][perm], na.rm = TRUE)
  })
  
  p_empirical <- (1 + sum(perm_stats >= obs_stat)) / (1 + n_perm)
  
  list(
    observed_statistic = obs_stat,
    permuted_statistics = perm_stats,
    p_empirical = p_empirical
  )
}

tmp <- perm_test_within_vs_across_pdx(a)
print(tmp$p_empirical)


p_angles <- ggplot(a, aes(x=xcoord, y=ycoord,fill=pval))+
  geom_raster()+
  # first legend: Treatment
  geom_point(data=dtop,  aes(color=treat1))+
  geom_point(data=dside, aes(color=treat1))+
  scale_fill_viridis_c(expression(Pr(theta[null]<theta[obs])),direction=-1)+
  scale_color_manual("",values=c("black","red"),labels=c("cisplatin","control"))+
  # switch to a new colour scale
  new_scale_color()+
  # second legend: PDX model
  geom_point(data=dtop,  aes(y=ycoord+1,color=pdx1))+
  geom_point(data=dside, aes(x=xcoord-1,color=pdx1))+
  scale_color_brewer(name="cell line", palette="Set2")+
  geom_vline(xintercept = breaks, colour="white", size=0.4) +
  geom_hline(yintercept = treat_breaks, colour="white", size=0.2)+
  geom_vline(xintercept = treat_breaks, colour="white", size=0.2)+
  geom_hline(yintercept = breaks, colour="white", size=0.4) +
  coord_equal()+
  theme_void()+
  text_size_theme+
  theme(
    axis.line   = element_blank(),
    axis.ticks  = element_blank(),
    axis.text   = element_blank(),
    axis.title  = element_blank(),
    panel.grid  = element_blank(),
    legend.box.just      = "left",
    legend.box.margin     = margin(0, 5, 0, 0, "pt"),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA)
  )
p_angles


source("R/utils_lineage.R")

df <- read.csv("data/raw/salehi/metadata.csv") 
df <- prune_children(df)

# Create dummy root to fix missing parent values
dummy_uid <- "ROOT"
dummy_row <- df[1, ]
dummy_row$uid <- dummy_uid
dummy_row$datasetname <- "ROOT"
df_fixed <- rbind(df, dummy_row)
df_fixed$parent[is.na(df_fixed$parent) | df_fixed$parent == ""] <- dummy_uid

df_fixed <- cbind(df_fixed[, c("uid", "parent")],
                  df_fixed[, !colnames(df_fixed) %in% c("uid", "parent")])

df_fixed$highlight <- df_fixed$uid%in%xm$uid & df_fixed$parent%in%xm$uid 

traj_col <- do.call(rbind,lapply(split(x$fi,f=x$pdx),function(tmp){
  dfi <- data.frame(tname=tmp)
  dfi$col <- letters[1:nrow(dfi)]
  dfi
}))

tmp <- traj_col$col
names(tmp) <- traj_col$tname
traj_col <- tmp

traj_lut <- traj_col[xm$fi]
names(traj_lut) <- paste0("x",xm$uid)
df_fixed$traj <- traj_lut[paste0("x",df_fixed$uid)]
# Build edge list and graph ----------------------------------------------
edges_df <- df_fixed %>% 
  filter(uid != dummy_uid) %>% 
  mutate(from_dummy = (parent == dummy_uid)) %>% 
  select(from = parent, to = uid,  from_dummy,traj,highlight)

g <- graph_from_data_frame(d = edges_df, vertices = df_fixed, directed = TRUE)

tg <- as_tbl_graph(g) %>% 
  mutate(dummy = (name == dummy_uid),
         node_type = ifelse(dummy, "dummy", ifelse(parent == dummy_uid, "root", "non-root")))

# Define layout and compute label positions -------------------------------
set.seed(42)
layout <- create_layout(tg, layout = "dendrogram", circular = TRUE,
                        height = -node_distance_to(1, mode = "all"))
# the plot gets rotated for certain versions of the ggtree package.
# this flip reorients it to how it was originally setup:
do_switch <- FALSE
if(do_switch){
  tmp <- layout$x; layout$x <- layout$y; layout$y <- tmp
  layout$y <- -layout$y
}


lineage_labels <- layout %>% 
  filter(parent == dummy_uid) %>%
  filter(name != dummy_uid) %>% 
  mutate(x_label = x,y_label = y)

sample_lineage_map2 <- c(
  SA1035 = "SA1035",
  SA532  = "SA532",
  SA609  = "SA609",
  SA906  = "p53 k.o",
  SA039  = "p53 w.t",
  SA535  = "SA535"
)

lineage_labels$linlab <- sample_lineage_map2[lineage_labels$PDX_id]

lineage_positions <- data.frame(x_label=c(0.08,0.4,0.5,-0.1,-0.2,-0.3),
                                y_label=c(0.6,0.6,0.2,-0.3,-0.1,0.8),
                                row.names = lineage_labels$linlab)
lineage_labels$x_label <- lineage_positions[lineage_labels$linlab,"x_label"]
lineage_labels$y_label <- lineage_positions[lineage_labels$linlab,"y_label"]

# --- Plot 6: Dendrogram/Network Plot ---
cb_palette <- c("#E69F00", "#56B4E9", "#009E73", 
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7","pink")


p_network <- ggraph(layout) + 
  geom_segment(data = lineage_labels,
               aes(x = x, y = y, xend = x_label, yend = y_label),
               color = "gray80", linetype = "dashed") +
  geom_edge_link(aes(linetype=ifelse(from_dummy, "blank", "solid")),
                 edge_width = 0.5,show.legend = F,color="grey80") +
  geom_edge_link(aes(color=toupper(traj),
                     linetype=ifelse(!highlight, "blank", "solid")),
                 edge_width = 1.5) +
  # Plot dummy nodes as white points
  geom_node_point(data = filter(layout, dummy),
                  aes(x = x, y = y), color = "white") +
  # Plot non-dummy nodes colored by xval
  geom_node_point(data = filter(layout, !dummy),
                  aes(x = x, y = y), size = 2) +
  scale_edge_linetype_manual(values = c(solid = "solid", blank = "blank"),
                             guide = "none") +
  scale_edge_color_manual(name = "trajectory", values = cb_palette, na.translate = FALSE) +
  geom_text(data = lineage_labels,
            aes(x = x_label, y = y_label, label = linlab),
            fontface = "bold",size = base_text_size / .pt,family = "sans") +
  text_size_theme+
  theme(legend.position = c(0.01, 0.01),
        legend.justification = c("left", "bottom"),
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

plotList <- lapply(plotList,function(ppi) ppi+ 
                     theme(
                       panel.background = element_rect(fill = "white", colour = NA),
                       plot.background  = element_rect(fill = "white", colour = NA)
                     ))


p_sim <- p_sim+ 
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA)
  )

list(p_network,plotList[[1]],plotList[[2]],plotList[[3]],p_sim,plotList[[4]],p_angles)
px <- cowplot::plot_grid(plotList[[1]],plotList[[2]],labels=c("B","C"),ncol=1)
pxx <- cowplot::plot_grid(plotList[[3]],p_sim,labels=c("D","E"),ncol=1)

ptop <- cowplot::plot_grid(p_network,p_angles,ncol=2,labels=c("A","G"),rel_widths = c(1,1))

pbot <- cowplot::plot_grid(px,pxx,plotList[[4]],nrow=1,labels=c("","","F"),rel_widths=c(3,3,2))

plt <- cowplot::plot_grid(ptop,pbot,ncol=1,rel_heights = c(3,4))
ggsave("figs/cna_fitness_effects.png",width=7,height=7,units="in",bg="white")
