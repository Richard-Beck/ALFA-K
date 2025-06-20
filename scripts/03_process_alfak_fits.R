if(!basename(getwd())=="ALFA-K") stop("Ensure working directory set to ALFA-K root")
dataDir <- "data/processed/salehi/alfak_inputs/"
fitDir <- "data/processed/salehi/alfak_outputs/"
outputDir <- "data/processed/salehi/alfak_outputs_proc/"

source("R/utils_env.R")
ensure_packages(c("transport","parallel","pbapply"))


lineages <- readRDS("data/processed/salehi/lineages.Rds")
m <- read.csv("data/raw/salehi/metadata.csv")

treat_lut <- m$on_treatment
names(treat_lut) <- paste0("u",m$uid)

ff <- list.files(fitDir)

inputs <- list.files(dataDir)

## for each fit read the xval score
xval_res <- function(fi,mi){
  fpath <- paste(c(fitDir,fi,mi),collapse="/")
  avail_data <- list.files(fpath)
  ## if no xval score assume pipeline failed (e.g. there are no karyotypes above min_obs threshold)
  if(!"xval.Rds"%in%avail_data) return(NaN)
  xval <- as.numeric(readRDS(paste(c(fpath,"xval.Rds"),collapse="/")))
  return(xval)
}

prediction_res <- function(dj,fi,mi){
  
  fpath <- paste(c(fitDir,fi,mi),collapse="/")
  avail_data <- list.files(fpath)
  ## if no xval score assume pipeline failed (e.g. there are no karyotypes above min_obs threshold)
  if(!"abm_preds"%in%avail_data) return(NULL)
  
  abm_outs <- paste0(paste0(fpath,"/abm_preds/"),list.files(paste0(fpath,"/abm_preds")))
  preds <- lapply(abm_outs,readRDS)
  
  # Extract matrix parts (drop 'time'), stack, and average
  mat_avg <- Reduce(`+`, lapply(preds, \(x) as.matrix(x[,-1]))) / length(preds)
  
  # Recover time column
  time <- preds[[1]][,1]
  
  # Recombine
  preds <- as.data.frame(mat_avg)
  
  
  fj <- paste(m[m$uid==dj,c("datasetname","timepoint")],collapse="_")
  if(sum(grepl(fj,inputs))<1) return(NULL)
  fj <- inputs[grepl(fj,inputs)][1]
  print(fj)
  

  
  k <- do.call(rbind,lapply(colnames(preds),s2v))
  
  test_treat <- as.character(treat_lut[paste0("u",dj)])
  

  tru <- readRDS(paste(c(dataDir,fj),collapse="/"))$x
  ## the way this function is setup we're assuming we're only comparing against the next
  ## passage. So the last column of tru corresponds to the first row of preds
  ## or tail(as.numeric(colnames(tru$x))*15,1) == head(rownames(preds),1)
  
  preds <- preds[-1,]
  
  yt <- 0*preds[1,]
  yb <- 0*preds[1,]
  
  tru <- tru[rownames(tru)%in%colnames(preds),(ncol(tru)-1):ncol(tru)]
  
  if(ncol(tru)!=2) stop("expected 2 columns here")
  yb[rownames(tru)] <- tru[,1]/sum(tru[,1])
  yt[rownames(tru)] <- tru[,2]/sum(tru[,2])
  yt <- as.numeric(yt)
  yb <- as.numeric(yb)
  if(!is.finite(sum(yt))) return(NULL)
  if(!is.finite(sum(yb))) return(NULL)
  
  
  df <- do.call(rbind,pbapply::pblapply(1:nrow(preds),function(tt){
    yp <- as.numeric(preds[tt,])
    if(!is.finite(sum(yp))) return(NULL)
    
    
    if(sum(yt) <= 0) return(NULL)
    if(min(yp)<0) return(NULL)
    ## calculate the performance metrics
    
    ##euclidean:
    kt <- colSums(yt*k)
    kb <- colSums(yb*k)
    kp <- colSums(yp*k)
    
    mag <- function(v) sqrt(sum(v^2))
    getangle <- function(a,b) 180*acos(sum(a*b)/(mag(a)*mag(b)))/pi
    angle=getangle(kt-kb,kp-kb)

    
    df_de <- data.frame(pred0=sqrt(sum((kp-kb)^2)),
                        pred=sqrt(sum((kp-kt)^2)),
                        base=sqrt(sum((kb-kt)^2)),metric="euclidean")
    
    ##wasserstein:
    st <- wpp(k,mass=yt)
    sb <- wpp(k,mass=yb)
    sp <- wpp(k,mass=yp)
    df_dw <- data.frame(pred0=wasserstein(sb,sp),
                        pred=wasserstein(st,sp),
                        base=wasserstein(st,sb),metric="wasserstein")
    
    ## cosine:
    sc_0 <- sum(yb*yt)/prod(sqrt(sum(yb^2)) , sqrt(sum(yt^2)))
    sc_p <- sum(yp*yt)/prod(sqrt(sum(yp^2)) , sqrt(sum(yt^2)))
    sc_p0 <- sum(yp*yb)/prod(sqrt(sum(yp^2)) , sqrt(sum(yb^2)))
    df_sc <- data.frame(pred0=sc_p0,pred=sc_p,base=sc_0,metric="cosine")
    
    ##overlap:
    so_0 <- sum(pmin(yb, yt))/min(sum(yb), sum(yt))
    so_p <- sum(pmin(yp, yt))/min(sum(yp), sum(yt))
    so_p0 <- sum(pmin(yp, yb))/min(sum(yp), sum(yb))
    df_so <- data.frame(pred0=so_p0,pred=so_p,base=so_0,metric="overlap")
    
    df <- rbind(df_de,df_dw,df_sc,df_so)
    df$tt <- tt
    df$angle <- angle
    df
  }))
  
  
  test_treat <- as.character(treat_lut[paste0("u",dj)])
  
  df$test_treat <- test_treat
  df
}


cl <- makeCluster(getOption("cl.cores", 50))

clusterExport(cl,varlist=c("lineages","m","outputDir","fitDir","dataDir","xval_res","prediction_res","inputs","treat_lut"))
clusterEvalQ(cl,{
  library(transport)
  source("R/utils_karyo.R")
})
dir.create(outputDir,recursive = T,showWarnings = F)
df <- parLapplyLB(cl,ff,function(fi){
  print(fi)
  min_obs <- list.files(paste0(fitDir,fi))
  dfi <- do.call(rbind,lapply(min_obs,function(mi){
    li <- lineages[[fi]]
    xv <- data.frame(fi=fi,min_obs=gsub("minobs_","",mi),ntrain = length(li$ids),
                      t_train = m$timepoint[m$uid==tail(li$ids,1)],
                      train_treat = paste(treat_lut[paste0("u",li$ids)],collapse=""),
                     xv=xval_res(fi,mi))
    
    xp0 <- data.frame(pred0=NaN,pred=NaN,base=NaN,metric=NaN,tt=NaN,angle=NaN,test_treat=NaN)
    
    
    
    xp <- lapply(li$dec1,prediction_res,fi=fi,mi=mi)
    xp <- xp[!sapply(xp,is.null)]
    xp <- do.call(rbind,xp)
    
    res <- NULL
    if(is.null(xp)){
      res <- cbind(xv,xp0)
    }else{
        res <- cbind(xv,xp)
    }
    saveRDS(res,paste0(outputDir,fi,"_",mi,".Rds"))
    return(res)
  }))
})












