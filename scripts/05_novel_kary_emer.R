if(!basename(getwd())=="ALFA-K") stop("Ensure working directory set to ALFA-K root")

proc_data <- function(lscape, x) {
  ## 1) figure out which karyotypes “just appeared” vs. “old” 
  last_col <- ncol(x)
  prev_col <- last_col - 1
  
  present     <- x[, last_col] > 0
  just_appeared <- present & (x[, prev_col] == 0)
  
  newk <- rownames(x)[ just_appeared ]
  oldk <- rownames(x)[ present ]
  
  ## 2) build numeric matrices oldv (M × p) and vn (N × p) exactly as before
  #    (where M = length(oldk), N = length(lscape$k))
  oldv <- do.call(rbind, lapply(oldk, s2v))      # M rows
  vn   <- do.call(rbind, lapply(lscape$k, s2v))   # N rows
  
  ## 3) compute the full pairwise Manhattan‐distance matrix via dist(...)
  combined <- rbind(vn, oldv)
  d_all    <- as.matrix(dist(combined, method = "manhattan"))
  N        <- nrow(vn)
  M        <- nrow(oldv)
  # slice out the N × M block of “distances from each vn[i,] to each oldv[j,]”
  d_mat    <- d_all[ 1:N, (N + 1):(N + M) ]
  
  ## 4) grab “previous‐timepoint counts” exactly as in your old code
  #    Note: this is x[oldk, prev_col], but the normalization constant should be sum(x[,prev_col])
  prev_counts_oldk <- x[oldk, prev_col]   # length‐M
  total_prev       <- sum(x[, prev_col])   # THIS is the exact same denominator your original used
  
  ## 5) for d = 0..5, do one matrix‐multiply per “distance‐band”
  #    result[i, d+1] = sum(prev_counts_oldk[j] for all j where d_mat[i,j] == d)
  result <- matrix(0, nrow = N, ncol = 6)
  for (d in 0:5) {
    mask <- (d_mat == d)                # logical N×M
    # (mask %*% prev_counts_oldk) is an N‐vector: for each i, sum prev_counts_oldk[j] where d_mat[i,j] == d
    result[, d + 1] <- mask %*% prev_counts_oldk
  }
  
  ## 6) normalize exactly the same way your original did
  result <- result / total_prev
  
  ## 7) turn into a data.frame, add “f” and “y”, then filter on d0 == 0 and drop that column
  df <- as.data.frame(result, row.names = NULL)
  names(df) <- paste0("d", 0:5)
  
  df$f <- lscape$mean
  df$y <- lscape$k %in% newk
  
  ## keep only rows where d0 == 0, then drop “d0” (the first column)
  df <- df[ df$d0 == 0, -1, drop = FALSE ]
  
  ## (optional) reset rownames so they go 1,2,3,...
  rownames(df) <- NULL
  
  return(df)
}

get_descendents_ids <- function(id,lineages,meta){

  lin <- lineages[[id]]
  if (length(lin$dec1) == 0) return(NULL)
  
  daughter_names <- meta$datasetname[meta$uid %in% lin$dec1]
  daughter_passages <- meta$timepoint[meta$uid %in% lin$dec1]
  id_patterns <- paste(daughter_names,daughter_passages,"l",2,sep="_")
  
  as.character(sapply(id_patterns,function(idi) names(lineages)[grepl(idi,names(lineages))]))
  
  
}
source("R/utils_karyo.R")
source("R/utils_lineage.R")
source("R/utils_env.R")
source("R/utils_theme.R")
ensure_packages(c("ggplot2","reshape2","dplyr","viridis","gridExtra","ggforce"))
## --- Global Theme Definition ---------------------
base_text_size <- 5
common_theme <- make_base_theme("classic",base_text_size)
void_theme <-  make_base_theme("void",base_text_size)

REGENERATE_DATA <- FALSE ## this analysis takes about 5-10 mins
if(REGENERATE_DATA){
  df_meta <- do.call(rbind,lapply(list.files("data/processed/salehi/alfak_outputs_proc/",full.names = T),readRDS))
  df_meta <- df_meta[!is.na(df_meta$xv)&df_meta$tt==1&df_meta$metric=="cosine",]
  xv_thresh <- 0.0
  df_meta <- df_meta[df_meta$xv>xv_thresh,]
  df_meta$tail_id <- sapply(df_meta$fi,function(fij){
    strsplit(fij,split = "_l_") |> unlist() |> head(1)
  })
  df_meta <- split(df_meta,f=df_meta$tail_id)
  df_meta <- do.call(rbind,lapply(df_meta,function(ii) ii[ii$xv==max(ii$xv),]))
  
  ## we need to find the ids of the passages we will predict. It would have been better if these
  ## had originally been stored in the alfak_outputs. We do that here:
  meta <- read.csv("data/raw/salehi/metadata.csv")
  lineages <- readRDS("data/processed/salehi/lineages.Rds")
  dec_df <- do.call(rbind,lapply(unique(df_meta$fi),function(ii){
    decids <- get_descendents_ids(ii,lineages,meta)
    data.frame(fi=ii,dec_id=decids)
  }))
  
  ## Previously we did store summary data for predictions of these passages - resulting in
  ## duplicate entires differing only in their prediction summary scores (cosine in this case).
  ## we need to remove these before merging to avoid spurious duplicates:
  df_meta <- unique(df_meta[,c("fi","min_obs","ntrain","t_train","train_treat","xv","tail_id")])
  df_meta <- merge(df_meta,dec_df)
  
  ## Now ready to do the emergent karyotype checks
  df <- do.call(rbind,pbapply::pblapply(1:nrow(df_meta),function(i){
    res <- tryCatch({
      lscape <- readRDS(paste0("data/processed/salehi/alfak_outputs/",df_meta$fi[i],"/minobs_",df_meta$min_obs[i],"/landscape.Rds"))
      x <- readRDS(paste0("data/processed/salehi/alfak_inputs/",df_meta$dec_id[i],".Rds"))$x
      
      test <- proc_data(lscape,x)
      mod <- glm(y~.,data=test,family = binomial)
      res <- summary(mod)$coefficients
      ids <- rownames(res)
      res <- data.frame(res,row.names = NULL)
      res$ids <- ids
      res <- merge(res,df_meta[i,])
      res
    },error=function(e) return(NULL))
    return(res)
  }))
  
  i <- 10
  lscape <- readRDS(paste0("data/processed/salehi/alfak_outputs/",df_meta$fi[i],"/minobs_",df_meta$min_obs[i],"/landscape.Rds"))
  x <- readRDS(paste0("data/processed/salehi/alfak_inputs/",df_meta$dec_id[i],".Rds"))$x
  test <- proc_data(lscape,x)
  i1 <- unlist(test[which.max(test$d1),1:5])
  i2 <- unlist(test[which.min(test$d1),1:5])
  
  dfxmpl <- data.frame(frac=c(i1,i2),d=c(names(i1),names(i2)),
                       id=c(rep("1",length(i1)),rep("2",length(i2))))
  
  
 
  res <- list(df=df,dfxmpl=dfxmpl)
  saveRDS(res,"data/processed/salehi/novel_kary_predictions.Rds")
}

res <- readRDS("data/processed/salehi/novel_kary_predictions.Rds")
pthresh <- 0.05



df <- res$df
dfxmpl <- res$dfxmpl

meta <- read.csv("data/raw/salehi/metadata.csv")
df$pdx <- sapply(df$fi,assign_labels,meta=meta)
print("unique training lineages:")
print(table(df$pdx[!duplicated(df$fi)]))
print(length(df$pdx[!duplicated(df$fi)]))
print("unique test instances:")
print(table(df$pdx[df$ids==df$ids[1]]))
print(length(df$pdx[df$ids==df$ids[1]]))

pc <- ggplot(dfxmpl,aes(x=d,y=frac,fill=id))+
  geom_col(position="dodge")+
  scale_fill_discrete("novel\nkaryotype")+
  common_theme+
  scale_x_discrete("distance from novel karyotype")+
  scale_y_continuous("population fraction")+
  theme(legend.position = c(0.7,0.9),legend.key.size = unit(0.1,"in"))
pc


w <- split(df,f=df$fi)


w <- do.call(rbind,lapply(w,function(wi){
  wi <- wi[-1,]
  wi <- wi[wi$Estimate>0,]
  wi <- wi[which.min(wi$Pr...z..),]
  wi
}))
z <- split(df,f=df$ids)


z <- do.call(rbind,lapply(z,function(zi){
  data.frame(var=zi$ids[1],frac=mean(zi$Estimate>0 & zi$Pr...z..<pthresh))
}))

renamr <- c(d1="d[1]",d2="d[2]",d3="d[3]",d4="d[4]",d5="d[5]",f="f")
z <- z[!z$var=="(Intercept)",]
z$var <- renamr[z$var]
w$id <- renamr[w$id]



# Define the data manually
data <- data.frame(
  Time = rep(c(1, 2, 3, 4), times = 5),
  Clone = rep(c("Clone1", "Clone2", "Clone3", "Clone4", "Clone5"), each = 4),
  Frequency = c(10, 5, 3, 1,
                8, 4, 2, 1,
                4, 2, 0, 0,
                0, 2, 4, 8,
                0, 0, 1, 3)
)

# Normalize frequencies at each timepoint to sum to 1
data <- data %>%
  group_by(Time) %>%
  mutate(Frequency = Frequency / sum(Frequency))

# Add group information
data <- data %>%
  mutate(Group = case_when(
    Clone %in% c("Clone1", "Clone2") ~ "Group A",
    Clone == "Clone3" ~ "Group B",
    Clone %in% c("Clone4", "Clone5") ~ "Group C"
  ))

clone_names <- c("Clone1", "Clone2", "Clone3", "Clone4", "Clone5")
zeta <- clone_names[1:3]
psi <- clone_names[c(1:2,4:5)]
psiNzeta <- clone_names[c(4:5)]

groups <- list(zeta=zeta,
               psi=psi,
               psiNzeta=psiNzeta)
tmp <- do.call(rbind,lapply(unique(data$Group),function(di) {
  tmp <- data
  tmp$Group2 <- di
  tmp
}))

tmp$Clone2 <- FALSE
tmp$Clone2[tmp$Group2=="Group A" & tmp$Clone%in%groups[[1]]] <- TRUE
tmp$Clone2[tmp$Group2=="Group B" & tmp$Clone%in%groups[[2]]] <- TRUE
tmp$Clone2[tmp$Group2=="Group C" & tmp$Clone%in%groups[[3]]] <- TRUE

lablr <- c('Group B' = " \u03A8 ",
           'Group A' = "  \u03B6  ",
           'Group C' = "\u0060\u03B6\u2229\u03A8")

tmp$Group2 <- lablr[tmp$Group2]

# Create a Muller plot using ggplot2
p0 <- ggplot(tmp, aes(x = Time, y = Frequency, group=Clone,fill = Clone,alpha=Clone2)) +
  geom_area(color="black") +
  scale_fill_viridis_d("karyotype",labels=letters[1:5])+
  guides(alpha="none")+
  scale_alpha_manual(values=c(0,1))+
  facet_wrap(~ Group2, ncol = 3) +
  common_theme+
  scale_x_continuous("",breaks=c(1,4),labels=c(expression(S[0],S[t])))+
  scale_y_continuous("clone frequency")
p0


circles <- data.frame(
  x0 = c(1,4),
  y0 = c(0,0),
  r = c(2,2),
  id = c("\u03B6","\u03A8")
)

circles$id <- factor(circles$id,levels=c("\u03A8","\u03B6"))

anndf <- data.frame(x=c(1,4,-1.25),y=c(-2.5,-2.5,2.5),anno=c("\u03B6","\u03A8","\u0398"))
anndf2 <- data.frame(x=c(1,4,-1.1),y=c(-2.5,-2.5,2.5),
                     anno=c("\u03B6","\u0060\u03B6\u2229\u03A8","\u0060\u03B6\u2229\u0060\u03A8"))

noaxistheme <- theme(axis.text = element_blank(),axis.title = element_blank(),
                     axis.line = element_blank(),axis.ticks = element_blank())
margin_pts <- 5
margintheme <- theme(plot.margin = margin(margin_pts,margin_pts,margin_pts,margin_pts)) 
# Behold some circles
pa <- ggplot() +
  geom_rect(aes(xmin = -2, xmax = 7, ymin = -3, ymax = 3),alpha=0,color="black")+
  geom_circle(aes(x0 = x0, y0 = y0, r = r), data = circles)+
  geom_text(data=anndf,aes(x=x,y=y,label=anno))+
  void_theme+
  noaxistheme+
  margintheme
pa

pb <- ggplot() +
  geom_rect(aes(xmin = -2, xmax = 7, ymin = -3, ymax = 3,fill="\u0398"),
            color="black", show.legend = F)+
  geom_circle(aes(x0 = x0, y0 = y0, r = r,fill=id), 
              data = circles,show.legend = F)+
  geom_text(data=anndf2,aes(x=x,y=y,label=anno),parse=F)+
  void_theme+
  noaxistheme+
  margintheme
pb



pd <- ggplot(z,aes(x=var,y=frac))+
  geom_col()+
  common_theme+
  scale_x_discrete("variable predicting novel karyotype",labels = scales::parse_format())+
  scale_y_continuous(paste0("fraction significant\n(p=",pthresh,")"))
pd


pe <- ggplot(w,aes(x=ids))+ 
  stat_count(aes(y=..count../sum(..count..)))+
  common_theme+
  scale_x_discrete("variable predicting novel karyotype",labels = scales::parse_format())+
  scale_y_continuous("fraction most\nsignificant")
pe

frac_f_signif <- z["f","frac"]
paste("fitness significant in",round(nrow(w)*frac_f_signif),"/",nrow(w),"tests")
n_f_most_signif <- table(w$ids)["f"]
paste("fitness most significant in",n_f_most_signif,"/",nrow(w),"tests")

sec1 <- cowplot::plot_grid(pc,pd,pe,labels=c("D","E","F"),ncol=1,label_size = base_text_size+2)
sec2 <- cowplot::plot_grid(pa,pb,labels=c("B","C"),ncol=2,label_size = base_text_size+2)

sec3 <- cowplot::plot_grid(p0,sec2,labels=c("A",""),ncol=1,label_size = base_text_size+2)

plt <- cowplot::plot_grid(sec3,sec1,nrow=1,rel_widths=c(5,2))

ggsave("figs/novel_kary_emer.png",plot=plt,width=180,height=90,units="mm",bg="white")




