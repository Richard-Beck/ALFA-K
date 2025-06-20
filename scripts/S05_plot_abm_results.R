if(!basename(getwd())=="ALFA-K") stop("Ensure working directory set to ALFA-K root")
pkgs <- c("cli","ggplot2","reshape2","scales","stringr","cowplot","ggalluvial")

## set true if the ABM summaries need to be compiled, false if not to speed up a bit
REAGGREGATE_DATA <- FALSE

source("R/utils_env.R")
ensure_packages(pkgs)
source("R/utils_theme.R")
source("R/utils_karyo.R")

base_text_size <- 5
common_theme <- make_base_theme("classic",base_text_size)
minimal_theme <- make_base_theme("minimal",base_text_size)
project_forward_log <- function(x0, f, timepoints) {
  K <- length(x0)
  out <- matrix(NA, nrow = K, ncol = length(timepoints))
  log_x0 <- log(x0)
  for (i in seq_along(timepoints)) {
    lv <- log_x0 + f * timepoints[i]
    denom <- logSumExp(lv)
    out[, i] <- exp(lv - denom)
  }
  out
}

logSumExp <- function(v) {
  m <- max(v)
  m + log(sum(exp(v - m)))
}



## --- Data Loading ---------------------

if(REAGGREGATE_DATA){
  all_files <- list.files(
    path       = "data/processed/ABM/",
    pattern    = "abm_preds_summary.Rds",
    recursive  = TRUE,
    full.names = TRUE
  )
  
  all_results <- lapply(all_files,readRDS)
  
  dfl <- do.call(rbind,lapply(all_results,function(xij){
    cbind(xij$summary,xij$info)
  }))
  
  reformat_df <- function(df_,metric){
    colnames(df_) <- c("prediction","baseline")
    df_$passage <- 1:nrow(df_)
    df_$metric <- metric
    df_
  }
  
  dfp <- do.call(rbind,lapply(all_results,function(xij){
    tryCatch({
      tmp <- rbind(cbind(reformat_df(xij$df_so,"overlap"),xij$df_a),
                   cbind(reformat_df(xij$df_sc,"cosine"),xij$df_a),
                   cbind(reformat_df(xij$df_de,"euclidean"),xij$df_a),
                   cbind(reformat_df(xij$df_dw,"wasserstein"),xij$df_a))
      
      cbind(tmp,xij$info,xij$summary)
    },error=function(e) return(NULL))
  }))
  
  dfp$win <- FALSE
  dfp$win[dfp$metric%in%c("overlap","cosine")] <- (dfp$prediction>dfp$baseline)[dfp$metric%in%c("overlap","cosine")] 
  dfp$win[!dfp$metric%in%c("overlap","cosine")] <- (dfp$prediction<dfp$baseline)[!dfp$metric%in%c("overlap","cosine")]
  dfp$pos_xv <- dfp$Rxv>0
  
  res <- list(dfl=dfl,dfp=dfp)
  saveRDS(res,"data/processed/ABM_summaries.Rds")
}

results_summary <- readRDS("data/processed/ABM_summaries.Rds")
results_summary$dfp$passage <- results_summary$dfp$passage-1
dfl <- results_summary$dfl


## --- Helper Function for Aggregated Data (for plot_pa) ---------------------
prepare_aggregated_data <- function(dfl) {
  z1 <- melt(dfl, measure.vars = c("r", "R", "rho"))
  z1$value <- pmax(-1, z1$value)
  z2 <- aggregate(list(smm = z1$value), 
                  by = list(w = z1$w, minobs = z1$minobs, ntp = z1$ntp, metric = z1$variable),
                  quantile, probs = c(0.1, 0.5, 0.9))
  smm <- z2$smm
  z2 <- cbind(z2[,-5], smm)
  colnames(z2)[5:7] <- c("lo", "med", "hi")
  return(z2)
}

agg_data <- prepare_aggregated_data(dfl)
col_labels <- setNames(
  paste0("lambda==", as.numeric(gsub("p", ".", unique(agg_data$w)))),
  unique(agg_data$w)
)
custom_labels <- c(
  r   = "Pearson~r",
  R   = "Rescaled~R^2",
  rho = "Spearman~rho"
)


## --- Plotting Functions ---------------------

plot_pa <- function(z2, col_labels, custom_labels) {
  ggplot(z2, aes(x = stringr::str_pad(minobs, width = 2), color = as.character(ntp))) +
    facet_grid(rows = vars(metric), cols = vars(w), scales = "free",
               labeller = labeller(
                 w = as_labeller(col_labels, default = label_parsed),
                 metric = as_labeller(custom_labels, default = label_parsed))) +
    geom_errorbar(aes(ymin = lo, ymax = hi, group = interaction(ntp, minobs)),
                  position = position_dodge(width = 0.8)) +
    geom_point(aes(y = med), position = position_dodge(width = 0.8)) +
    scale_y_continuous("Metric value") +
    scale_x_discrete("Minimum Observations", expand = expansion(add = c(0.2, 0.8))) +
    scale_color_viridis_d("longitudinal\nsamples") +
    common_theme
}

plot_pb <- function(dfl) {
  df <- dfl
  df$nfq_bins <- ifelse(df$nfq > 8, ">8", "<=8")
  df$failVar <- ifelse(df$Rfq > 0, "R[~f]^2>0", "R[~f]^2<0")
  ggplot(df, aes(x = nfq_bins, y = pmax(-1, R))) +
    facet_grid(rows = vars(failVar), labeller = label_parsed) +
    geom_violin() +
    scale_x_discrete("num. freq. karyotypes", labels = c(">8" = ">8", "<=8" = bquote("≤" ~ 8))) +
    scale_y_continuous(expression(rescaled~R^2~(landscape))) +
    common_theme +
    coord_flip()
}

plot_pc <- function(base_dir, id, k1, k2, k3) {
 
  ypath <- file.path("data/raw/ABM",paste0(id,".Rds"))
  yi <- readRDS(ypath)
  
  wavelength <- (strsplit(id,split="_") |> unlist())[2]
  wavelength <- gsub("p",".",wavelength) |> as.numeric()
  
  dirs <- list.files(file.path(base_dir,id))
  
  
  z_list <- lapply(dirs, function(conds) {
    parts <- strsplit(conds, "_")[[1]]
    minobs <- parts[2]
    ntp <- parts[4]
    lscape_path <- file.path(base_dir, id, conds, "landscape.Rds")
    lscape <- readRDS(lscape_path)
    lscape <- lscape[lscape$fq, ]
    k <- do.call(rbind, lapply(lscape$k, s2v))
    f <- apply(k, 1, function(ki) getf(ki, wavelength, yi$true_landscape))
    lscape$f_tru <- f - mean(f)
    lscape$mean  <- lscape$mean - mean(lscape$mean)
    lscape$minobs <- minobs
    lscape$ntp <- ntp
    lscape
  })
  z <- do.call(rbind, z_list)
  z$hi <- z$k %in% c(k1, k2, k3)
  ggplot(z, aes(x = mean, y = f_tru, color = hi)) +
    facet_grid(rows = vars(paste0(stringr::str_pad(ntp, 2), " samples")),
               cols = vars(paste0("N=", stringr::str_pad(minobs, width = 2))),
               scales = "free") +
    geom_abline() +
    geom_point(show.legend = FALSE) +
    scale_color_manual("", values = c("FALSE" = "black", "TRUE" = "#FFA500")) +
    scale_x_continuous("true fitness", breaks = scales::pretty_breaks(n = 2)) +
    scale_y_continuous("estimated fitness") +
    common_theme
}

## Modified prepare_data now accepts a keys argument.
prepare_data <- function(cond, id, keys) {
  parts <- strsplit(cond, "_")[[1]]
  n_final_cols <- as.numeric(parts[4])
  bootstrap_path <- file.path("data/processed/ABM",id, cond, "bootstrap_res.Rds")
  x1 <- readRDS(bootstrap_path)
  
  ypath <- file.path("data/raw/ABM",paste0(id,".Rds"))
  yi <- readRDS(ypath)
  
  # Process yp
  yp <- yi$abm_output$x
  yp <- yp[, as.numeric(colnames(yp)) < 120, drop = FALSE]
  yp <- yp[, (ncol(yp) - n_final_cols + 1):ncol(yp), drop = FALSE]
  yp <- yp[order(rowSums(yp), decreasing = TRUE), ]
  yp <- head(yp, 8)
  for (i in seq_len(ncol(yp))) {
    yp[, i] <- yp[, i] / sum(yp[, i])
  }
  yp$Var1 <- rownames(yp)
  yp <- melt(yp,id.vars="Var1")
  yp$Var2 <- as.numeric(as.character(yp$variable))
  yp$hi <- yp$Var1 %in% keys
  
  # Process timepoints
  tt <- as.numeric(colnames(yi$abm_output$x)) 
  tt <- round(tail(tt[tt < 120], n_final_cols))
  tt <- seq(min(tt), max(tt), length.out = 100)
  
  # Process xp by projecting forward
  xp <- do.call(rbind, lapply(1:nrow(x1$final_frequencies), function(i) {
    xp_mat <- project_forward_log(x1$final_frequencies[i, ], x1$final_fitness[i, ], timepoints = tt)
    colnames(xp_mat) <- tt
    rownames(xp_mat) <- colnames(x1$final_frequencies)
    xp_mat <- xp_mat[order(rowSums(xp_mat), decreasing = TRUE), , drop = FALSE]
    xp_mat <- xp_mat[rownames(xp_mat) %in% unique(as.character(yp$Var1)), , drop = FALSE]
    for (j in seq_len(ncol(xp_mat))) {
      xp_mat[, j] <- xp_mat[, j] / sum(xp_mat[, j])
    }
    xp_mat <- melt(xp_mat)
    xp_mat$rep <- i
    xp_mat
  }))
  xp$hi <- xp$Var1 %in% keys
  
  list(xp = xp, yp = yp)
}

plot_pd_pe <- function(data8, data2, lut_label) {
  p_common <- function(xp, yp) {
    ggplot(xp, aes(x = as.numeric(Var2), y = value)) +
      facet_wrap(~Var1, nrow = 2, scales = "free",
                 labeller = as_labeller(lut_label)) +
      geom_line(aes(group = rep, color = hi), alpha = 0.5) +
      geom_point( data = yp, size = 2,color="#56B4E9") +
      scale_color_manual(values = c("FALSE" = "black", "TRUE" = "#FFA500")) +
      scale_x_continuous("sim time (days)", breaks = scales::pretty_breaks(n = 2)) +
      scale_y_continuous("kary. freq.", breaks = scales::pretty_breaks(n = 2)) +
      common_theme +
      theme(legend.position = "none")
  }
  list(pd = p_common(data8$xp, data8$yp),
       pe = p_common(data2$xp, data2$yp))
}

plot_p2a <- function(dfl, col_labels, custom_labels) {
  y <- dfl
  y$Rxv[!is.finite(y$Rxv)] <- -Inf
  y$Rxv <- pmax(-1, y$Rxv)
  z <- aggregate(list(value = y$Rxv),
                 by = list(ntp = y$ntp, minobs = y$minobs, w = y$w),
                 quantile, probs = c(0.1, 0.5, 0.9))
  dfz <- data.frame(z$value)
  colnames(dfz) <- c("lo", "med", "hi")
  z <- cbind(z[,-ncol(z)], dfz)
  
  ggplot(z, aes(x = stringr::str_pad(minobs, width = 2), color = as.character(ntp))) +
    facet_grid(cols = vars(w), scales = "free",
               labeller = labeller(
                 w = as_labeller(col_labels, default = label_parsed),
                 metric = as_labeller(custom_labels, default = label_parsed))) +
    geom_errorbar(aes(ymin = lo, ymax = hi, group = interaction(ntp, minobs)),
                  position = position_dodge(width = 0.8)) +
    geom_point(aes(y = med), position = position_dodge(width = 0.8)) +
    scale_y_continuous(expression(R[X]^2)) +
    scale_x_discrete("Observation threshold (N)", expand = expansion(add = c(0.2, 0.8))) +
    scale_color_viridis_d("longitudinal\nsamples") +
    common_theme
}

plot_p2c <- function(dfl, custom_labels) {
  y <- melt(dfl, measure.vars = c("rho", "r", "R"))
  y_split <- split(y, interaction(y$ntp, y$w, y$abmrep))
  y_filtered <- do.call(rbind, lapply(y_split, function(xi) xi[xi$Rxv == max(xi$Rxv), , drop = FALSE]))
  z <- aggregate(list(value = y_filtered$value),
                 by = list(filter = y_filtered$Rxv > 0, ntp = y_filtered$ntp, metric = y_filtered$variable),
                 quantile, probs = c(0.1, 0.5, 0.9))
  dfz <- data.frame(z$value)
  colnames(dfz) <- c("lo", "med", "hi")
  z <- cbind(z[,-ncol(z)], dfz)
  
  ggplot(z, aes(x = ntp, color = filter)) +
    facet_grid(rows = vars(metric), scales = "free",
               labeller = labeller(metric = as_labeller(custom_labels, default = label_parsed))) +
    geom_point(aes(y = med, group = filter), position = position_dodge(width = 0.8)) +
    geom_errorbar(aes(ymin = lo, ymax = hi, group = filter), position = position_dodge(width = 0.8)) +
    scale_y_continuous("metric value") +
    scale_x_discrete("num. longitudinal samples") +
    scale_color_viridis_d("", labels = c(expression(R[X]^2<0),
                                         expression(R[X]^2>0))) +
    common_theme +
    theme(legend.position = "top")
}

plot_p4a <- function(dfp) {
  dfpaf <- aggregate(list(fwin = dfp$win),
                     by = list(metric = dfp$metric, passage = dfp$passage, ntp = dfp$ntp, xv = dfp$Rxv > 0),
                     mean, na.rm = TRUE)
  dfpaf$label_filter <- ifelse(dfpaf$xv, "R[X]^2>0", "R[X]^2<0")
  ggplot(dfpaf, aes(x = passage, y = fwin, color = ntp)) +
    facet_grid(cols = vars(metric), rows = vars(label_filter), labeller = label_parsed, scales = "free") +
    geom_line(size = 1) +
    geom_hline(yintercept = 0.5, color = "red", linetype = 2) +
    scale_y_continuous("beat baseline", limits = c(0, 1)) +
    scale_x_continuous("passage number", breaks = seq(2, 10, 2)) +
    scale_color_viridis_d("longitudinal\nsamples") +
    common_theme
}

plot_p4b <- function(dfp, nulldf) {
  dfa <- dfp[dfp$metric == dfp$metric[1] & dfp$passage %in% c(1, 10), ]
  ggplot(dfa, aes(x = angle)) +
    facet_grid(rows = vars(paste("passage", passage)),
               cols = vars(paste(stringr::str_pad(ntp, width = 2), "samples")),
               scales = "free") +
    stat_ecdf(aes(color = as.character(Rxv > 0))) +
    geom_line(data = nulldf, aes(x = angle, y = CDF), linetype = 2) +
    scale_color_viridis_d("", labels = c(expression(R[X]^2<0), expression(R[X]^2>0))) +
    scale_y_continuous("cum. dist.") +
    scale_x_continuous("angle metric") +
    common_theme
}

## --- Alluvium Plot Functions ---------------------
# First alluvium plot (pf): using R, nfq, Rfq, etc.
plot_pf_alluvium <- function(dfl) {
  df <- dfl
  df$R[!is.finite(df$R)] <- -1
  z <- aggregate(list(n = df$r),
                 by = list(
                   ntp   = df$ntp,
                   nfq   = df$nfq > 8,
                   fq_r  = df$Rfq > 0,
                   w     = df$w,
                   adjR2 = df$R > 0),
                 length)
  z$fillvar <- "v1"
  z$fillvar[!z$fq_r] <- "v2"
  z$fillvar[z$fq_r & !z$nfq] <- "v3"
  z$w <- gsub("p", ".", z$w)
  
  ggplot(z, aes(axis1 = w, axis2 = ntp, axis3 = fq_r, axis4 = nfq, axis5 = adjR2, y = n)) +
    geom_alluvium(aes(fill = fillvar)) +
    geom_stratum(width = 1/3, color = "black") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)),
              size = common_theme$text$size * 0.3) +
    scale_fill_manual(
      name = "",
      values = c(
        "v1" = "grey70",  # greyscale fill
        "v2" = "#21908CFF",  # or use viridis::viridis(3)[2]
        "v3" = "#440154FF"   # or viridis::viridis(3)[1]
      ),
      labels = c(
        "v2" = expression(R[f]^2 < 0),
        "v3" = expression(N[f] <= 8)
      ),
      breaks = c("v2", "v3")  # omit v1 from legend
    ) +
    scale_x_discrete(
      labels = c("axis1" = expression(lambda),
                 "axis2" = "longitudinal\nsamples",
                 "axis3" = expression(R[f]^2 > 0),
                 "axis4" = expression(N[f] > 8),
                 "axis5" = expression(R^2 > 0)),
      limits = c("axis1", "axis2", "axis3", "axis4", "axis5"),
      expand = c(0.1, 0.1)) +
    minimal_theme +
    theme(
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      panel.grid = element_blank()
    )
  
}

# Second alluvium plot (p2b): using Rxv and R.
plot_p2b_alluvium <- function(dfl) {
  df <- dfl
  df$Rxv[!is.finite(df$Rxv)] <- -Inf
  df$Rxv <- pmax(-1, df$Rxv)
  z <- aggregate(list(n = df$r),
                 by = list(
                   ntp   = df$ntp,
                   xv    = df$Rxv > 0,
                   w     = df$w,
                   adjR2 = df$R > 0),
                 length)
  z$w <- gsub("p", ".", z$w)
  
  ggplot(z, aes(axis1 = w, axis2 = ntp, axis3 = xv, axis4 = adjR2, y = n)) +
    geom_alluvium(aes(fill = adjR2), show.legend = FALSE) +
    geom_stratum(width = 1/3, color = "black") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)),
              size = common_theme$text$size * 0.3) +
    scale_fill_viridis_d() +
    scale_x_discrete(
      labels = c("axis1" = expression(lambda),
                 "axis2" = "longitudinal\nsamples",
                 "axis3" = expression(R[X]^2 > 0),
                 "axis4" = expression(R[R]^2 > 0)),
      limits = c("axis1", "axis2", "axis3", "axis4"),
      expand = c(0.1, 0.1)) +
    theme_void() + common_theme +
    theme(
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      panel.grid = element_blank()
    )
}


## --- Additional Data for Some Plots ---------------------

## If using new data the manual highlights below won't work or may
## be meaningless. In general need to try a few different id values, run the next
## few lines of code (as far as plot_pd_pe(result8, result2, lut_label))
## then look for cases where the final two points of any panels in pd go against
## the trend for the highlighting.

base_dir <- "data/processed/ABM"

id <- "w_1p6_m_5e-05_rep_74"

# Define keys for highlighting in prepare_data and plot_pc.
k1 <- "2.2.2.2.2.2.2.2.2.2.2.2.2.2.2.2.2.2.2.1.2.2"
k2 <- "2.2.2.2.2.2.2.2.2.2.2.2.2.2.2.2.2.2.3.2.2.2"
k3 <- NULL#"2.2.2.2.2.2.2.2.1.2.2.2.2.2.2.2.2.2.2.2.2.2"
highlight_keys <- c(k1, k2, k3)

# For p4a and p4b, process dfp from the results_summary.
dfp <- results_summary$dfp
dfp <- split(dfp, interaction(dfp$abmrep, dfp$ntp, dfp$w))
dfp <- do.call(rbind, lapply(dfp, function(di) {
  if(sum(is.finite(di$Rxv)) < 1) return(NULL)
  di <- di[is.finite(di$Rxv), ]
  di[di$Rxv == max(di$Rxv), ]
}))

# For p4b, prepare a nulldf for spherical angle CDF.
dSphereAngle <- function(theta, N) {
  coef <- integrate(function(t) sin(t)^(N - 2), lower = 0, upper = pi)$value
  sin(theta)^(N - 2) / coef
}
cSphereAngle <- function(theta, N) {
  if (!is.finite(theta)) return(1)
  integrate(function(t) dSphereAngle(t, N), lower = 0, upper = theta)$value
}
nulldf <- data.frame(angle = 0:180)
nulldf$rads <- pi * nulldf$angle / 180
nulldf$CDF <- sapply(nulldf$rads, cSphereAngle, N = 22)


## --- Assemble and Return Plot List ---------------------
plots <- list(
  pa    = plot_pa(agg_data, col_labels, custom_labels),
  pb    = plot_pb(dfl),
  pc    = plot_pc(base_dir, id,  k1, k2, k3)
)

# Prepare data for pd and pe plots (with highlighting based on keys).
result8 <- prepare_data("minobs_20_ntp_8", id, highlight_keys)
result2 <- prepare_data("minobs_20_ntp_2", id, highlight_keys)
all_names <- unique(c(as.character(result8$yp$Var1), as.character(result2$yp$Var1)))
ordered_levels <- sort(all_names)
lut_label <- setNames(LETTERS[seq_along(ordered_levels)], ordered_levels)
lut_label[lut_label%in%c("C","G")]
# Reset factor levels for proper facet ordering.
result8$xp$Var1 <- factor(result8$xp$Var1, levels = ordered_levels)
result8$yp$Var1 <- factor(result8$yp$Var1, levels = ordered_levels)
result2$xp$Var1 <- factor(result2$xp$Var1, levels = ordered_levels)
result2$yp$Var1 <- factor(result2$yp$Var1, levels = ordered_levels)

pd_pe <- plot_pd_pe(result8, result2, lut_label)
plots$pd  <- pd_pe$pd
plots$pe  <- pd_pe$pe

plots$p2a <- plot_p2a(dfl, col_labels, custom_labels)
plots$p2c <- plot_p2c(dfl, custom_labels)
plots$p4a <- plot_p4a(dfp)
plots$p4b <- plot_p4b(dfp, nulldf)

# Add the two alluvium plots:
plots$pf   <- plot_pf_alluvium(dfl)
plots$p2b  <- plot_p2b_alluvium(dfl)

topr <- cowplot::plot_grid(plots$pa,plots$pb,labels=c("A","B"), label_size = base_text_size+2,rel_widths = c(2.5,1))
midrt <- cowplot::plot_grid(plots$pe,plots$pd,labels=c("D","E"), label_size = base_text_size+2,nrow=2)
midr <- cowplot::plot_grid(plots$pc,midrt,labels=c("C",""), label_size = base_text_size+2,nrow=1,rel_widths = c(2,3))

plt <- cowplot::plot_grid(topr,midr,plots$pf,labels=c("","","F"), label_size = base_text_size+2,nrow=3,rel_heights = c(4,4,3))
ggsave("figs/ABM_validation_p1.png",plt,width=150,height=175,units="mm",bg="white")
## Final output: a list of ggplot objects

left <- cowplot::plot_grid(plots$p2a,plots$p2b,labels=c("A","B"), label_size = base_text_size+2,nrow=2,rel_heights = c(3,4))
plt <- cowplot::plot_grid(left,plots$p2c,nrow=1,labels=c("","C"), label_size = base_text_size+2,rel_widths = c(3,2))
ggsave("figs/ABM_validation_p2.png",plt,width=150,height=100,units="mm",bg="white")

plt <- cowplot::plot_grid(plots$p4b,plots$p4a,labels=c("A","B"), label_size = base_text_size+2,nrow=2)
ggsave("figs/ABM_validation_p3.png",plt,width=150,height=100,units="mm",bg="white")
  