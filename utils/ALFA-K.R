library(lhs)
library(fields)
##prepare ABM data for input to pipeline
proc_sim <- function(dir,times){
  dt <- read.csv(paste0(dir,"/log.txt"),sep=",")$dt
  
  filenames <- list.files(dir)
  filenames <- filenames[!filenames%in%c("proc_data","log.txt")]
  tx <- as.numeric(substr(filenames,1,5))
  filenames <- paste0(dir,"/",sapply(times, function(ti) filenames[which.min(abs(ti-tx))]))
  tx <- sapply(times, function(ti) tx[which.min(abs(ti-tx))])
  
  x <- lapply(filenames,function(fi){
    xi <- read.csv(fi,header=F)
    nchrom <- ncol(xi)-2
    k <- apply(xi[,1:nchrom],1,paste,collapse=".")
    n <- xi[,nchrom+1]
    fitness <- xi[,nchrom+2]
    data.frame(karyotype=k,n=n,fitness=fitness)
  })
  
  
  
  pop.fitness <- sapply(x,function(xi) sum(xi$n*xi$fitness)/sum(xi$n))
  tmp <- do.call(rbind,x)  
  tmp <- tmp[!duplicated(tmp$karyotype),]
  unique.karyotypes <- tmp$karyotype
  clone.fitness <- tmp$fitness
  
  x <- do.call(rbind,lapply(as.character(unique.karyotypes), function(k){
    sapply(x, function(xi) sum(xi$n[xi$karyotype==k]))
  }))
  
  
  colnames(x) <- tx
  rownames(x) <- unique.karyotypes  
  list(x=x,pop.fitness=pop.fitness,clone.fitness = clone.fitness,dt=dt)
}

fitm <- function(par,xj,dt){
  N <- colSums(xj)
  x0i <- exp(head(par,length(par)/2))
  fi <- tail(par,length(par)/2)
  
  logx <- do.call(rbind,lapply(1:length(x0i), function(i){
    fi[i]*as.numeric(colnames(xj))*dt + log(x0i[i])
  }))
  pred <- apply(logx,2,function(lxi) exp(lxi)/sum(exp(lxi)))
  pred <- matrix(pred,ncol=ncol(xj))
  
  ll <- do.call(rbind,lapply(1:nrow(xj), function(i){
    -dbinom(xj[i,],size = N,prob=pred[i,],log=T)
  }))
  ll <- ll[,!N==0]
  ll[!is.finite(ll)] <- 1e+9
  sum(ll)
}

inifitm <- function(inipar, opar,xj,dt){
  N <- colSums(xj)
  x0i <- c(exp(head(opar,length(opar)/2)),exp(inipar[1]))
  fi <- c(tail(opar,length(opar)/2),inipar[2])
  
  logx <- do.call(rbind,lapply(1:length(x0i), function(i){
    fi[i]*as.numeric(colnames(xj))*dt + log(x0i[i])
  }))
  pred <- apply(logx,2,function(lxi) exp(lxi)/sum(exp(lxi)))
  
  ll <- do.call(rbind,lapply(1:nrow(xj), function(i){
    -dbinom(xj[i,],size = N,prob=pred[i,],log=T)
  }))
  ll <- ll[,!N==0]
  ll[!is.finite(ll)] <- 1e+9
  sum(ll)
  
}
opt_g_free <- function(x,min_obs=5,mintp=0){
  
  dt <- x$dt
  x <- x$x
  x <- x[rowSums(x)>min_obs,]
  ntp <- apply(x,1,function(xi) sum(xi>0))
  x <- x[ntp>mintp,]
  x <- x[order(rowSums(x),decreasing=T),]
  ##some initial parameter guesses
  x0 <- log(1/10^9)
  f0 <- 0.5
  
  i <- 1
  xj <- x[1:i,,drop=F]
  par <- c(x0,f0)
  
  opt <- optim(par,fitm,xj=xj,dt=dt)
  par <- opt$par
  err <- opt$value
  x0 <- head(par,length(par)/2)
  f0 <- tail(par,length(par)/2)
  
  for(i in 2:nrow(x)){
    print(i)
    guesses <- randomLHS(n=100,k=2)
    xj <- x[1:i,,drop=F]
    rx0 <- max(x0+25)-min(x0-25)
    guesses[,1] <- (guesses[,1]-0.5)*rx0+median(x0)
    
    inipar <- apply(guesses,1,function(gi){
      iniopt <- optim(gi,fn = inifitm,opar=par,xj=xj,dt=dt)
      c(iniopt$value,iniopt$par)
    })
    
    inipar <- inipar[2:3,which.min(inipar[1,])]
    par <- c( c(x0,inipar[1]),c(f0,inipar[2]))
    opt <- optim(par,fitm,xj=xj,dt=dt)
    par <- opt$par
    err <- c(err,opt$value)
    x0 <- head(par,length(par)/2)
    f0 <- tail(par,length(par)/2)
  }
  min_f0 <- min(f0)
  delta_f0 <- 0.2-min_f0
  f0 <- f0+delta_f0
  xopt <- data.frame(f_est=f0,u0=x0,ll=err,n=rowSums(x),ntp=apply(x,1,function(xi) sum(xi>0)))
  
  rownames(xopt) <- rownames(xj)[1:nrow(xopt)]
  xopt
}


gen_all_neighbours <- function(ids,as.strings=T){
  if(as.strings) ids <- lapply(ids, function(ii) as.numeric(unlist(strsplit(ii,split="[.]"))))
  nkern <- do.call(rbind,lapply(1:length(ids[[1]]), function(i){
    x0 <- rep(0,length(ids[[1]]))
    x1 <- x0
    
    x0[i] <- -1
    x1[i] <- 1
    rbind(x0,x1)
  }))
  n <- do.call(rbind,lapply(ids, function(ii) t(apply(nkern,1,function(i) i+ii))))
  n <- unique(n)
  nids <- length(ids)
  n <- rbind(do.call(rbind,ids),n)
  n <- unique(n)
  n <- n[-(1:nids),]  
  n <- n[apply(n,1,function(ni) sum(ni<1)==0),]
  n
}


find_grad_dist <- function(x_opt){
  xmat <- do.call(rbind,lapply(rownames(x_opt), function(xi){
    as.numeric(unlist(strsplit(xi,split="[.]")))
  }))
  
  dx <- as.matrix(dist(xmat))
  dy <- as.matrix(dist(x_opt$f_est))
  dy <- dy[which(dx==1&upper.tri(dx))]
  dy <- c(dy,-dy)
  
  sd(dy)
}
pij<-function(i, j, beta){
  qij <- 0
  if(abs(i-j)>i){ ## not enough copies for i->j
    return(qij)
  }
  # code fails for j = 0, but we can get this result by noting i->0 always
  ## accompanies i->2i
  if(j==0) j <- 2*i
  s <- seq(abs(i-j), i, by=2 )
  for(z in s){
    qij <- qij + choose(i,z) * beta^z*(1-beta)^(i-z) * 0.5^z*choose(z, (z+i-j)/2)
  }
  ## if there is a mis-segregation: ploidy conservation implies that both daughter cells will emerge simultaneously
  #if(i!=j) qij <- 2*qij
  
  return(qij)
}
optim_neighbor_fitness <- function(f,tp,iflux,ftp,sdy,u,tt,pop_size,f_est,dt){
  xtp <- rep(0,length(iflux))
  xtp[1] <- iflux[1]
  for(j in 2:length(xtp)){
    ## in the following statement we want c1 to evaluate to -Inf when 
    ## xtp[j-1] = 1. Due to precision we can have cases where xtp[j-1] > 1
    ## which evaluates to NaN. Fixed by bounding input to log at zero
    c1 <- log(max(1/xtp[j-1]-1,0)) 
    t <- diff(tp)[1]*dt
    xi <- exp(f*t)/(exp(f*t)+exp(c1+ftp[j]*t))
    xtp[j] <- iflux[j]+xi
  }
  ux <- approx(tp,xtp,xout=tt)$y
  ux <- pmin(ux,0.9999)
  deltaf <- abs(f-f_est)
  negll <- c(-log(dbinom(u,pop_size,prob=ux)), -log(dnorm(deltaf,mean=0,sd=sdy)))
  negll[!is.finite(negll)] <- 10^9## to prevent warnings when probability is zero
  sum(negll)
}
get_neighbor_fitness <- function(ni,x_opt,x,pm0,ntp=100){
  ni_s <- paste(ni,collapse=".")
  
  sdy <- find_grad_dist(x_opt)
  
  ##for each 1ms neighbour, identify likely parents
  sources <- rownames(x_opt)[rownames(x_opt)%in%apply(gen_all_neighbours(ni_s),1,paste,collapse=".")]
  x_opt_i <- x_opt[sources,,drop=FALSE]
  fu <- x$x[sources,,drop=FALSE]/colSums(x$x)
  
  ## get the longitudinal frequency data for neighbor of interest
  u <- rep(0,ncol(x$x))
  names(u) <- colnames(x$x)
  if(ni_s%in%rownames(x$x)) u <- x$x[ni_s,]
  
  ##misseg probability from parent to neighbor
  pm <- sapply(sources,function(si){
    si <- as.numeric(unlist(strsplit(si,split="[.]")))
    prod(sapply(1:length(si), function(k) pij(si[k],ni[k],pm0)))
  })
  
  ##interpolate fitness, parent frequency and time
  ## essentially we are approximating an ODE from here on
  mintp <- min(as.numeric(colnames(x$x)))
  maxtp <- max(as.numeric(colnames(x$x)))
  tp <- seq(mintp,maxtp,(maxtp-mintp)/ntp)
  
  
  ftp <- approx(as.numeric(colnames(x$x)),x$pop.fitness,xout = tp)$y
  
  iflux <- colSums(do.call(rbind,lapply(1:length(pm),function(j) {
    aj <- approx(as.numeric(colnames(x$x)),as.numeric(fu[j,]),xout = tp)$y
    aj*diff(tp)[1]*x$dt*pm[j]*x_opt_i$f_est[j]
  })))
  
  oni <- optimise(optim_neighbor_fitness,interval=c(0,1),tp=tp,iflux=iflux,ftp=ftp,sdy=sdy,
                  u=u,tt=colnames(x$x),pop_size=colSums(x$x),f_est=x_opt_i$f_est,dt=x$dt)
  data.frame(f_est=oni$minimum,u0=NaN,ll=oni$objective,n=sum(u),ntp=sum(u>0))
}

infer_population_fitness <- function(x,x_opt){
  m <- x$x[rownames(x_opt),] 
  for(i in 1:ncol(m)){
    m[,i]<- m[,i]*x_opt$f_est/sum(m[,i])
  }
  colSums(m)
}

wrap_neighbor_fitness <- function(x,x_opt,pm0=0.00005,ntp=100){
  x$pop.fitness <- infer_population_fitness(x,x_opt)
  n <- gen_all_neighbours(rownames(x_opt))
  x2 <- do.call(rbind,lapply(1:nrow(n), function(i){
    get_neighbor_fitness(n[i,],x_opt,x,pm0,ntp)
  }))
  rownames(x2) <- apply(n,1,paste,collapse=".")
  x2
}


alfak <- function(x,min_obs=5,min_tp=0,misseg_rate=0.00005){
  x0 <- opt_g_free(x,min_obs,min_tp)
  x1 <- wrap_neighbor_fitness(x,x0,pm0=misseg_rate)
  x0$id <- "fq"
  x1$id <- "nn"
  xo <- rbind(x0,x1)
  xmat <- do.call(rbind,lapply(rownames(xo), function(i){
    as.numeric(unlist(strsplit(i,split="[.]")))
  }))
  y <- xo$f_est
  fit <- Krig(xmat,y,m=1)
  return(list(fit=fit,xo=xo))
}