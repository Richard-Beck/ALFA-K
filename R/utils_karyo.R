## helper function convert string karyotype 2.2.2.2 -> numeric
s2v <- function(charvec) as.numeric(unlist(strsplit(charvec,split="[.]")))

## generates all one missegregation neighbours of inputs karyotypes.
gen_all_neighbours <- function(ids, as.strings = TRUE, remove_nullisomes = TRUE) {
  if (as.strings) 
    ids <- lapply(ids, s2v)
  nkern <- do.call(rbind, lapply(1:length(ids[[1]]), function(i) {
    x0 <- rep(0, length(ids[[1]]))
    x1 <- x0
    x0[i] <- -1
    x1[i] <- 1
    rbind(x0, x1)
  }))
  n <- do.call(rbind, lapply(ids, function(ii) t(apply(nkern, 1, function(i) i + ii))))
  n <- unique(n)
  nids <- length(ids)
  n <- rbind(do.call(rbind, ids), n)
  n <- unique(n)
  n <- n[-(1:nids), ]
  if (remove_nullisomes) 
    n <- n[apply(n, 1, function(ni) sum(ni < 1) == 0), ]
  n
}

## returns fitness value using one of the GRF random landscapes as input (tru_lscape)
getf <- function(k,wavelength,tru_lscape){
  Nwaves <- nrow(tru_lscape)
  scalef <- 1/(pi*sqrt(Nwaves))
  d <- apply(tru_lscape,1,function(ci){
    sqrt(sum((k-ci)^2))
  })
  sum(sin(d/wavelength)*scalef)
}

## this function takes a vector of paths to landscapes and counts the 
## number of unique karyotypes. 
## files pointed to by paths are data frames with a column 'k' which contains
## karyotype strings
count_unique_k <- function(paths) {
  unique_env <- new.env(hash = TRUE, parent = emptyenv())
  pb <- txtProgressBar(min = 0, max = length(paths), style = 3)
  
  for (i in seq_along(paths)) {
    df <- readRDS(paths[i])
    ks <- df$k
    for (s in ks) unique_env[[s]] <- TRUE
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  length(ls(unique_env))
}

## generate a random fitness landscape.
gen_randscape <- function(founder,Nwaves,scalef=NULL,wavelength=0.8){
  if(is.null(scalef)) scalef <- 1/(pi*sqrt(Nwaves))
  f0 <- 0
  ## we want to have simulations where the diploid founder is fit enough to survive.
  ## if scalef <- 1/(pi*sqrt(30))
  ## then this would make the diploid cell in the top 10% fittest
  ## clones:
  while(f0<=0.4){
    pk <- lapply(1:Nwaves, function(i){
      pk <- sample((-10):20,length(founder),replace=T)
    })
    
    d <- sapply(pk,function(ci) {
      sqrt(sum((founder-ci)^2))
    })
    f0=sum(sin(d/wavelength)*scalef)
    
  }
  do.call(rbind,pk)
}

##get the vector angle between two karyotype shifts:
mag <- function(v) sqrt(sum(v^2))
getangle <- function(a,b) 180*acos(sum(a*b)/(mag(a)*mag(b)))/pi