How does missegregation rate modify effective population fitness?

``` r
source("utils/ALFA-K.R")
```

    ## Warning: package 'lhs' was built under R version 4.1.3

    ## Warning: package 'fields' was built under R version 4.1.3

    ## Loading required package: spam

    ## Warning: package 'spam' was built under R version 4.1.3

    ## Spam version 2.9-1 (2022-08-07) is loaded.
    ## Type 'help( Spam)' or 'demo( spam)' for a short introduction 
    ## and overview of this package.
    ## Help for individual functions is also obtained by adding the
    ## suffix '.spam' to the function name, e.g. 'help( chol.spam)'.

    ## 
    ## Attaching package: 'spam'

    ## The following objects are masked from 'package:base':
    ## 
    ##     backsolve, forwardsolve

    ## Loading required package: viridis

    ## Warning: package 'viridis' was built under R version 4.1.2

    ## Loading required package: viridisLite

    ## 
    ## Try help(fields) to get started.

``` r
source("utils/ode_functions.R")
```

    ## 
    ## Attaching package: 'Matrix'

    ## The following object is masked from 'package:spam':
    ## 
    ##     det

    ## Warning: package 'igraph' was built under R version 4.1.3

    ## 
    ## Attaching package: 'igraph'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     decompose, spectrum

    ## The following object is masked from 'package:base':
    ## 
    ##     union

``` r
library(ggplot2)
```

    ## Warning: package 'ggplot2' was built under R version 4.1.3

``` r
library(deSolve)
```

    ## Warning: package 'deSolve' was built under R version 4.1.3

``` r
library(rootSolve)
```

    ## Warning: package 'rootSolve' was built under R version 4.1.1

Example of the principles of how mis-segregation rate can shape
landscape exploration:

``` r
## ODE Model equation
chrmod <- function(time,state,parms){
  with(as.list(parms),{
    ds <- state%*%A
    ds <- ds-sum(ds)*state
    return(list(ds))
  })
}
reverse_state_index <- function(index,maxchrom,minchrom=1,ndim){
  ## reset maxchrom to behave as if minchrom was 1,1,...
  mc <- maxchrom-minchrom+1
  if(length(mc)==1 & ndim>1) mc <- rep(mc,ndim)
  ## how many elements does each part of the state represent?
  ## works as prod(numeric(0)) evaluates to 1:
  Nsites <- cumprod(mc)
  cp <- c(1,cumprod(mc)[-length(mc)])
  state <- c()
  for(j in ndim:1){
    ni <- floor((index-1)/cp[j])
    state <- c(ni+1,state)
    index <- index-cp[j]*ni
  }
  state + minchrom-1
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

get_A <- function(maxchrom,beta){
 
  Ndim <- length(maxchrom)
  Nstates <- prod(maxchrom) 
  A<- matrix(0,nrow=Nstates,ncol=Nstates)
  for(i in 1:nrow(A)){
    state_i <- reverse_state_index(i,maxchrom,minchrom=1,ndim=Ndim)
    for(j in 1:nrow(A)){
      state_j <- reverse_state_index(j,maxchrom,minchrom=1,ndim=Ndim)
      qij <- sapply(1:Ndim, function(k) pij(state_i[k], state_j[k], beta))
      ## joint probability (i1,i2,...)->(j1,j2,...)
      qij <- prod(qij)
      A[i,j] <- 2*qij
      
      # ## case when there is no mis-segregation:
      if(i==j) A[i,j] <- (2*qij-1)
      
    }
  }
  A
}
get_fitness <- function(k,pk1,pk2){
  f1 <- as.numeric(sqrt(sum((k-pk1$centre)^2))<=pk1$rad)*pk1$fitness
  f2 <- as.numeric(sqrt(sum((k-pk2$centre)^2))<=pk2$rad)*pk2$fitness
  max(f1,f2)
}
mean_kary <- function(pm,pk1,pk2,return_df=FALSE){
  A <- get_A(maxchrom=c(8,8),beta=pm)
  state_ids <- do.call(rbind,lapply(1:nrow(A), reverse_state_index,maxchrom=c(8,8),ndim=2))
  f <- apply(state_ids,1,get_fitness,pk1=pk1,pk2=pk2)
  #f <- f/max(f)
  A <- A*f
  parms <- list(A=A)
  out <- runsteady(y=rep(1,length(f))/length(f),func=chrmod,parms=parms)
  df <- data.frame(state_ids)
  colnames(df) <- c("c1","c2")
  df$f <- f
  df$x <- out$y
  df$pm <- pm
  df$deltaf <- abs(pk1$fitness-pk2$fitness)
  if(!return_df) return(sum((df$c1)*df$x))
  if(!return_df) return(sum((df$c1+df$c2)*df$x/2))
  return(df)
  #data.frame(cn=1:8,f=f,pm=pm,deltaf=deltaf,x=out[nrow(out),-1])
}
```

``` r
pk1 <- list(centre=c(6,6),rad=1,fitness=1)
pk2 <- list(centre=c(3,3),rad=1,fitness=0.9)

df1 <- rbind(mean_kary(0.1,pk1,pk2,return_df = T),
            mean_kary(0.0001,pk1,pk2,return_df = T))
df1$id <- "scenario A"

pk3 <- list(centre=c(3,6),rad=0,fitness=1)
pk4 <- list(centre=c(6,3),rad=1.99,fitness=0.9)
df2 <- rbind(mean_kary(0.1,pk3,pk4,return_df = T),
            mean_kary(0.0001,pk3,pk4,return_df = T))
df2$id <-"scenario B"

df <- rbind(df1,df2)

p1 <- ggplot(df,aes(x=c1,y=c2))+
  facet_grid(cols=vars(paste("misrate:",pm)),rows=vars(id))+
  geom_raster(aes(fill=f))+
  geom_point(aes(size=x))+
  scale_color_manual("",values=c("black"))+
  scale_size("karyotype \nfrequency")+
  scale_fill_viridis_c("fitness")+
  scale_x_discrete("chromosome 1 copy number")+
  scale_y_discrete("chromosome 2 copy number")
p1
```

![](README_files/figure-markdown_github/unnamed-chunk-3-1.png)

``` r
ggsave("figures/misseg_landscape_exploration/figures/example.png",width = 4,height=3,units="in")
```

We will restrict our analysis to fitted lineages with no parents or
descendents, and that had a cross validation score r^2 above 0.3.

``` r
## first generate the nonzero elements of the sparse transition matrices:
source("figures/misseg_landscape_exploration/sparseMatrixCoords.R")
## then find dominant eigenvalues at various missegregation rates:
source("figures/misseg_landscape_exploration/eigenScreen.R")
```

![Example of karyotype distribution affected by mis-segregation
rate.](screen_ims/SA535_CISPLATIN_CombinedH_X9_l_5_d1_0_d2_0.png)
![Example of karyotype distribution unaffected by mis-segregation
rate.](screen_ims/SA609UnBU_X7_l_6_d1_0_d2_0.png)

Most fitted landscapes did not show any evidence of peak-selection via
missegregation rate changes. However, a small number did and we ran
these with the ABM to validate the eigenvector-based predictions:

``` r
source("figures/misseg_landscape_exploration/sweep_abm_with_fits2.R")
```

``` r
proc_res <- function(fi){
  print(fi)
  df <- screenR(paste0(fi,".Rds"))
  
  mr <- list.files(paste0(dir,fi))
  z <- do.call(rbind,lapply(mr,function(mri){
    reps <- list.files(paste0(dir,fi,"/",mri,"/output/"))
   # print(reps)
    do.call(rbind,lapply(reps,function(ri){
      x <- proc_sim(paste0(dir,fi,"/",mri,"/output/",ri,"/"),times=c(seq(0,10000,500)))
      x <- x$x
      x <- x[order(rowSums(x),decreasing=T),]
      x <- data.frame(x,check.names = F)
      
      id <- rep("none",nrow(x))
      id[rownames(x)%in%rownames(df[df$id,])] <- "winlo"
      id[rownames(x)%in%rownames(df[!df$id,])] <- "winhi"
      
      z <- split(x,f=id)
      z <- do.call(rbind,lapply(z,colSums))
      z <- reshape2::melt(z)
      colnames(z) <- c("id","time","n")
      z$rep <- ri
      z$misrate <- mri
      z
    }))
  }))
  z$cellLine <- fi
  return(z)
}

dir <- "data/salehi/misseg_landscape_exploration/minobs_5/"
ff <- list.files(dir)

z <- do.call(rbind,lapply(ff,proc_res))

saveRDS(z,file = "figures/misseg_landscape_exploration/screen_validation.Rds")
```

``` r
options(scipen=999)
z <- readRDS("figures/misseg_landscape_exploration/screen_validation.Rds")
cellLines <- c("SA535_CISPLATIN_CombinedH_X9_l_5_d1_0_d2_0",
               "SA906a_X57_l_7_d1_0_d2_0")
z2 <- z[!z$id=="none"&z$cellLine%in%cellLines,]

transp <- function(p) as.numeric(gsub("p",".",p))

m <- readRDS("figures/salehi_data_fitting/labelled_metadata.Rds")
xx <- readRDS("figures/salehi_data_fitting/fit_summaries.Rds")

uids <- sapply(cellLines, function(ci) xx$uid[xx$filenames==paste0(ci,".Rds")][1])
linlabs <- sapply(uids, function(ui){
  paste(m$PDX_id,m$timepoint,m$linlab)[m$uid==ui]
})
names(linlabs) <- cellLines

z3 <- aggregate(list(n=z2$n),by=list(id=z2$id,
                                     time=z2$time,
                                     misrate=z2$misrate,
                                     cellLine=z2$cellLine),
                mean)

p <- ggplot(z3,aes(x=time/10,y=n,color=id))+
  facet_grid(cols=vars(paste0("p=",transp(misrate))),rows=vars(linlabs[cellLine]))+
  geom_point()+
  scale_y_log10("num. cells")+
  scale_x_continuous("days")+
  scale_color_discrete("favored\nwhen p is:",labels=c("high","low"))
p
```

    ## Warning: Transformation introduced infinite values in continuous y-axis

![](README_files/figure-markdown_github/unnamed-chunk-7-1.png)