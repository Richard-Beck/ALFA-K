## generate a random field landscape
gen_rf_landscape <- function(founder,Nwaves,scalef=NULL,wavelength=1){
  
  if(is.null(scalef)) scalef <- 1/(pi*sqrt(Nwaves))
  pk <- lapply(1:Nwaves, function(i){
    #pk <- sample(0:10,length(founder),replace=T)
    pk <- sample((-10):20,length(founder),replace=T)
  })
  
  d <- sapply(pk,function(ci) {
    sqrt(sum((founder-ci)^2))
  })

  return(pk)
}

##get fitness from random field landscape
get_rf_fitness <- function(k,pk,scalef=NULL,wavelength=1){
  if(is.null(scalef)) scalef <- 1/(pi*sqrt(length(pk)))
  d <- sapply(pk,function(ci) {
    sqrt(sum((k-ci)^2))
  })
  f0=sum(sin(d/wavelength)*scalef)
  return(f0)
}

get_ruggedness <- function(p,pk,scalef=NULL,wavelength=1){
  Nchrom <- length(p)
  n <- expand.grid(lapply(1:Nchrom, function(novar) c(-1,1)))
  n <- n[!apply(n,1,function(ni) sum(abs(ni))==0),]
  n <- t(apply(n,1,function(ni) ni+p))
  fn <- apply(n,1,get_rf_fitness,pk=pk,scalef=scalef,wavelength=wavelength)
  fp <- get_rf_fitness(p,pk,scalef=scalef,wavelength=wavelength)
  mean(abs(fp-fn))
}

get_complexity <- function(v,pk,scalef=NULL,wavelength=1){
  f0 <- get_rf_fitness(v,pk,scalef=scalef,wavelength=wavelength)
  mean(sapply(1:length(v), function(i){
    vi <- v
    vi[i] <- vi[i]+1
    ff <- get_rf_fitness(vi,pk,scalef=scalef,wavelength=wavelength)
    
    vi <- v
    vi[i] <- vi[i]-1
    fr <- get_rf_fitness(vi,pk,scalef=scalef,wavelength=wavelength)
    
    sd(c(ff-f0,f0-fr))
  }))
}

assess_complexity <- function(wavelength,Nwaves,Nchrom=4,Nreps= 10){
  print(c(wavelength,Nwaves,Nchrom))
  x <- sapply(1:Nreps, function(novar){
    pk <- gen_rf_landscape(founder=rep(2,Nchrom),Nwaves = Nwaves,
                           wavelength=wavelength)
    p <- rep(2,Nchrom)
    get_complexity(p,pk,wavelength=wavelength)
  })
  
  data.frame(wavelength=wavelength,Nwaves=Nwaves,Nchrom=Nchrom,
             mean.complexity=mean(x),sd.complexity=sd(x))
}

assess_ruggedness <- function(wavelength,Nwaves,Nchrom=4,Nreps= 10){
  print(c(wavelength,Nwaves,Nchrom))
  x <- sapply(1:Nreps, function(novar){
    pk <- gen_rf_landscape(founder=rep(2,Nchrom),Nwaves = Nwaves,
                           wavelength=wavelength)
    p <- rep(2,Nchrom)
    get_ruggedness(p,pk,wavelength=wavelength)
  })
  
  data.frame(wavelength=wavelength,Nwaves=Nwaves,Nchrom=Nchrom,
             mean.ruggedness=mean(x),sd.ruggedness=sd(x))
}

abm_from_krig <- function(fit,dir,pars=NULL,cpp_cmd="ABM/bin/ABM.exe"){
  options(scipen = 999)
  dir.create(dir,recursive = T)
  dir.create(paste0(dir,"/train"))
  knots <- fit$knots
  cc <- fit$c
  d <- fit$d
  fscape <- rbind(cbind(knots,cc),c(d))
  fitted_landscape_path <- paste0(dir,"/landscape.txt")
  write.table(fscape, fitted_landscape_path,row.names = FALSE,col.names=FALSE,sep=",")
  
  parnames <- c("init_kary","fitness_landscape_type","fitness_landscape_file",
                "dt","p","pgd","Nsteps","output_dir","init_size")
  parvals <- c("2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2","krig",
               paste0(dir,"/landscape.txt"),"0.1","0.00005","0.0","3000",
               paste0(dir,"/train"),"100000")

  config <- cbind(parnames,parvals)
  for(i in 1:length(pars)){
    parname <- names(pars)[i]
    parval <- as.character(pars[i])
    config[config[,1]==parname,2] <- parval
  }
  config <- apply(config,1,paste,collapse=",")
  config_path <- paste0(dir,"/config.txt")
  writeLines(config,config_path)
  cmd <- paste(cpp_cmd,config_path)
  return(cmd)
}


roughness_meas <- function(landscape,only_fq=F){
  m <- landscape$fit
  x <- landscape$xo
  if(!is.null(ncol(x))) x <- list(x)
  if(only_fq) x <- lapply(x, function(xi) xi[xi$id=="fq",])
  k <- unlist(lapply(x, rownames))
  k <- unique(k)
  
  xmat <- do.call(rbind,lapply(k, function(i) as.numeric(unlist(strsplit(i,split="[.]")))))
  nn <- lapply(k,gen_all_neighbours)
  f0 <- c(predict(m,xmat))
  roughness <- sapply(1:length(nn), function(i){
    f0 <- f0[i]
    fn <- predict(m,nn[[i]])
    mean(abs(f0-fn))
  })
  
  ploidy <- apply(xmat,1,mean)
  
  data.frame(ploidy=ploidy,roughness=roughness,f_est=f0)
}

nviable_meas <- function(landscape,only_fq=F){
  m <- landscape$fit
  x <- landscape$xo
  if(!is.null(ncol(x))) x <- list(x)
  if(only_fq) x <- lapply(x, function(xi) xi[xi$id=="fq",])
  k <- unlist(lapply(x, rownames))
  k <- unique(k)
  
  xmat <- do.call(rbind,lapply(k, function(i) as.numeric(unlist(strsplit(i,split="[.]")))))
  nn <- lapply(k,gen_all_neighbours)
  f0 <- c(predict(m,xmat))
  fviable <- sapply(1:length(nn), function(i){
    f0 <- f0[i]
    fn <- predict(m,nn[[i]])
    mean(fn>f0)
  })
  
  ploidy <- apply(xmat,1,mean)
  
  data.frame(ploidy=ploidy,fviable=fviable,f_est=f0)
}

moran_freq <- function(karyotypes,fit=NULL,fitness=NULL,m=2){
  if(is.null(fit)&is.null(fitness)) stop("supply either landscape or fitness estimates")
  
  xmat <- do.call(rbind,lapply(karyotypes, function(i){
    as.numeric(unlist(strsplit(i,split="[.]")))
  }))
  if(is.null(fitness)) fitness <- predict(fit,xmat)
  ploidy <- apply(xmat,1,mean)
  xbar <- mean(fitness)
  N <- length(fitness)
  d <- as.matrix(dist(xmat))
  w <- (1/d^m)#d*(d==1)
  diag(w) <- 0
  w <- apply(w,1,function(wi) wi/sum(wi))
  zi <- fitness-xbar
  
  m2 <- sum(zi^2)/N
  
  Ii <- sapply(1:length(fitness), function(i){
    zi[i]/m2*sum(w[i,]*zi)
  })  
  data.frame(Ii,ploidy)
}