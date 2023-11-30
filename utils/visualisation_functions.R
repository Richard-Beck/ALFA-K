melt_for_plotting <- function(data_obj, nclones=5, fit_obj=NULL){
  dt <- data_obj$dt
  x <- data_obj$x[order(rowSums(data_obj$x),decreasing=T),]
  x <- head(x,nclones)
  for(i in 1:ncol(x)) x[,i] <- x[,i]/sum(x[,i])
  
  x <- reshape2::melt(x)
  colnames(x) <- c("karyotype","time","frequency")
  x$time <- x$time*dt
  
  out <- list(data=x,fit=NULL)
  
  if(!is.null(fit_obj)){
    y <- fit_obj$xo
    y <- y[unique(x$karyotype),]
    
    t <- seq(min(x$time),max(x$time),length.out=100)
    
    z <- y$u0+y$f_est%*%t(t)
    z <- apply(z,2,function(zi) exp(zi)/sum(exp(zi)))
    colnames(z) <- t
    rownames(z) <- rownames(y)
    
    z <- reshape2::melt(z)
    colnames(z) <- c("karyotype","time","frequency")    
    out$fit <- z
  }
  return(out)
}