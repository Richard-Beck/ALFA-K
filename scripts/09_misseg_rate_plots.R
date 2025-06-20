if(!basename(getwd())=="ALFA-K") stop("Ensure working directory set to ALFA-K root")
source("R/utils_env.R")
source("R/utils_karyo.R")
source("R/utils_theme.R")
source("R/utils_lineage.R")
ensure_packages(c("alfakR","ggplot2","ggrepel","rootSolve","reshape2","umap","cowplot")) ## is rootSolve actually used??
base_text_size <- 5
base_theme <- make_base_theme("classic",base_text_size)

make_example_data <- function(){
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
  
  df$pm[df$pm==0.0001] <- "0p0001"
  df
}

alfak_fit_dir <- "data/processed/salehi/alfak_outputs/"
df <- readRDS("data/processed/salehi/steadyStatePredictions.Rds")
m <- read.csv("data/raw/salehi/metadata.csv")

# --- 1. helper ---------------------------------------------------------------
heatmap_stat <- function(dfi){
  kdom  <- apply(dfi, 2, which.max)
  kvecs <- do.call(rbind, lapply(rownames(dfi)[kdom], s2v))
  d_dom <- as.matrix(dist(kvecs, method = "manh"))[1, ]
  data.frame(misrate = colnames(dfi), d_dom,  row.names = NULL)
}

# --- 2. stats per dataset ----------------------------------------------------
delta_dom <- lapply(df, heatmap_stat)          # list the same length as df

# --- 3. lineage distance matrix ---------------------------------------------
li              <- readRDS("data/processed/salehi/lineages.Rds")
lineage_members <- lapply(names(df), function(id) li[[id]]$ids)

# symmetric pair-wise distance: mean of elements unique to either set
dlin <- outer(seq_along(lineage_members), seq_along(lineage_members),
              Vectorize(function(i, j){
                a <- lineage_members[[i]]; b <- lineage_members[[j]]
                sum(!(a %in% b))/length(unique(c(a,b)))
              }))
colnames(dlin) <- rownames(dlin) <- names(df)

hc <- hclust(as.dist(dlin))                    # no error now
row_order <- hc$labels[hc$order]              # the desired ordering

# --- 4. long data frame for ggplot ------------------------------------------
delta_long <- do.call(rbind, Map(function(x, id)
  data.frame(lineage = id, x, stringsAsFactors = FALSE),
  delta_dom, names(delta_dom)))
delta_long$lineage <- factor(delta_long$lineage, levels = row_order)

mis_levels               <- unique(as.character(delta_long$misrate))
delta_long$misrate       <- factor(delta_long$misrate,
                                   levels = mis_levels[order(as.numeric(mis_levels))])

delta_long$cellLine <- sapply(as.character(delta_long$lineage),assign_labels,meta=m)



nk <- lapply(df,rownames) |> unlist() |> unique() 
cat("num, unique karyotypes:", length(nk), "\n")

ndom <- sapply(df,function(di) apply(di,2,which.max) |> unique() |>length())
cat("num with multiple dominators:", sum(ndom>1), "\n")
cat("out of landscapes:", length(ndom), "\n")



if(FALSE){
  ## this is a manual screen looking to see if there is a missegregation rate dependent switch in clonal dominance.
  pscreen <- lapply(names(df),function(i){
    tryCatch({
      print(i)
      dfi <- df[[i]]
      dfi <- dfi[apply(dfi,1,max)>0.001,]
      k <- do.call(rbind,lapply(rownames(dfi), s2v))
      print(nrow(k))
      u <- umap::umap(k)$layout
      u <- data.frame(u)
      u <- cbind(u,dfi)
      u <- reshape2::melt(u,id.vars=c("X1","X2"))
      p <- ggplot(u,aes(x=X1,y=X2,color=value))+
        facet_wrap(~variable)+
        geom_point()+
        scale_color_viridis_c(trans="log")+
        ggtitle(i)
      p
    },error=function(e) return(NULL))
    
    
  })
  stop("check results of manual screen: inspect plots, check for two distinct clusters and population transition with misseg rate")
}


get_plot_data <- function(dfi){
  min_obs <- attr(dfi,"min_obs")
  fi <- attr(dfi,"fi")
  
  dfi <- head(dfi,400) ## kind of an arbitrary threshold...
  ploidy <- sapply(rownames(dfi),function(i) mean(as.numeric(s2v(i))))
  lscape <- readRDS(paste0(alfak_fit_dir,fi,"/minobs_",min_obs,"/landscape.Rds"))
  
  
  lscape <- lscape[lscape$k%in%rownames(dfi),]
  
  lut <- lscape$mean
  names(lut) <- lscape$k
  fitness <- lut[rownames(dfi)]  
  
  ump <- umap::umap(do.call(rbind,lapply(rownames(dfi),s2v)))
  ump <- data.frame(ump$layout)
  
  cbind(dfi,ump,ploidy,fitness,fi)
}

ids <- c("SA906a_X50_l_4_d1_1_d2_0","SA532_X4_l_4_d1_1_d2_1")
df_list <- df[ids]

x <- do.call(rbind,lapply(df_list,get_plot_data))
x$grp <- x[,5]>x[,7]

lut <- c("p53 k.o. A\n(scenario A)","SA532\n(scenario B)")
names(lut) <- ids
x$fi <- lut[x$fi]
xre <- reshape2::melt(x,id.vars=c("X1","X2","ploidy","fitness","fi","grp"))
xre <- xre[order(xre$fitness),]

df_ex <- make_example_data()
transpchar <- function(x) as.character(gsub("p",".",x))
p_xmpl <- ggplot(df_ex,aes(x=c1,y=c2))+
  facet_grid(rows=vars(paste("misrate:",transpchar(pm))),cols=vars(id))+
  geom_raster(aes(fill=f))+
  geom_point(aes(size=x))+
  scale_color_manual("",values=c("black"))+
  scale_size("karyotype \nfrequency")+
  scale_fill_viridis_c("fitness",begin=0.1)+
  scale_x_discrete("chromosome 1 copy number")+
  scale_y_discrete("chromosome 2\ncopy number")+
  base_theme
p_xmpl

##setup labels for UMAP plots:
tmp <- xre[xre$variable%in%c(0.00032,0.00256),]
tmp <- split(tmp,f=interaction(tmp$variable,tmp$fi))
tmp <- do.call(rbind,lapply(tmp,function(ti) ti[ti$value==max(ti$value),]))
tmp <- tmp[,!colnames(tmp)%in%c("variable","value")]
tmp$lab <- c("y","x","y","x")

p_umap <- ggplot(xre[xre$variable%in%c(0.00032,0.00256),],aes(x=X1,y=X2))+
  facet_grid(rows=vars(variable),cols=vars(fi),scales="free")+
  geom_point(aes(color=pmax(10^-7,value)))+
  geom_label_repel(data=tmp,aes(label=lab),box.padding = .5, size = base_text_size / 2.845)+
  scale_color_viridis_c("freq.",trans="log",breaks=c(1e-1,1e-3,1e-5,1e-7))+
  base_theme
p_umap

xf_norm <- xre[xre$variable%in%c(0.00032),]
xf_norm <- split(xf_norm,f=xf_norm$fi)
xf_norm <- do.call(rbind,lapply(xf_norm, function(i){
  i$fitness <- i$fitness-mean(i$fitness)
  i
} ))

p_fit_umap <- ggplot(xf_norm,aes(x=X1,y=X2))+
  facet_grid(cols=vars(fi),scales="free")+
  geom_point(aes(color=fitness))+
  geom_label_repel(data=tmp,aes(label=lab),box.padding = .5, size = base_text_size / 2.845)+
  scale_color_viridis_c("fitness",option="plasma")+
  base_theme
p_fit_umap

p_fit <- ggplot(x,aes(x=fitness,fill=grp))+
  facet_grid(cols=vars(fi),scales="free")+
  geom_histogram(position = "identity",binwidth=0.02,alpha=0.5,color="black")+
  scale_fill_manual("group",labels=c("x","y"),values=c("#CC79A7","#00724D"))+
  base_theme
p_fit

p_ploidy <- ggplot(x,aes(x=ploidy,fill=grp))+
  facet_grid(cols=vars(fi),scales="free")+
  geom_histogram(position = "identity",binwidth=0.1,alpha=0.5,color="black")+
  scale_fill_manual("group",labels=c("x","y"),values=c("#CC79A7","#00724D"))+
  base_theme
p_ploidy

xf <- x[,colnames(x)%in%c("0.00064","0.00256","X1","X2","ploidy","fitness","fi")]

xf <- reshape2::melt(x,measure.vars=colnames(x)[1:8])
xfa <- aggregate(list(frequency=xf$value),by=list(p=xf$variable,fi=xf$fi,grp=xf$grp),sum)
p_freq <- ggplot(xfa,aes(x=p,y=frequency,color=grp,group=grp))+
  facet_grid(cols=vars(fi))+
  scale_color_manual("group",labels=c("x","y"),values=c("#CC79A7","#00724D"))+
  geom_point()+
  geom_line()+
  scale_x_discrete("misseg. rate")+
  base_theme+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
p_freq

highlight <- ids             # <- rows to outline, optional



## plot ---------------------------------------------------
p_screen <- ggplot(delta_long, aes(misrate, lineage, fill = d_dom)) +
  facet_grid(                                   # keep facet strips
    rows  = vars(cellLine),
    scales = "free_y",                          # each panel gets its own y scale
    space  = "free_y"                           # …and shrinks/grows to keep tile height constant
  ) +
  geom_raster() +
  geom_tile(                                    # red boxes
    data = subset(delta_long, lineage %in% highlight),
    colour = "red", fill = NA, linewidth = 0.35
  ) +
  scale_y_discrete(breaks = highlight,labels=lut) +        # only label highlighted rows
  scale_fill_viridis_c(expression(d[k]),trans = "sqrt") +
  base_theme+
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line=element_blank())

p_screen

pleft <- cowplot::plot_grid(p_xmpl,p_umap,labels=LETTERS[c(1,3)],nrow=2, label_size = 10)
pright <- cowplot::plot_grid(p_freq,p_fit_umap,p_fit,p_ploidy,labels=LETTERS[c(2,4:6)],ncol=1, label_size = 10)

plt <- cowplot::plot_grid(pleft,pright,rel_widths = c(4,3))

ggsave("figs/misseg.png",plt,width=125,height=125,units="mm")

