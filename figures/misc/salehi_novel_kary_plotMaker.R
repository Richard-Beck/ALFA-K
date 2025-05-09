setwd("~/projects/ALFA-K/")
## data in res should be generated with salehi_novel_kary.R
res <- readRDS("figures/misc/data/salehi_novel_kary_plotData.Rds")
pthresh <- 0.05

## --- Global Theme Definition ---------------------
base_text_size <- 8
common_theme <- theme_classic() + theme(
  text         = element_text(size = base_text_size, family = "sans"),
  axis.title   = element_text(size = base_text_size, family = "sans"),
  axis.text    = element_text(size = base_text_size, family = "sans"),
  legend.title = element_text(size = base_text_size, family = "sans"),
  legend.text  = element_text(size = base_text_size, family = "sans"),
  strip.text   = element_text(size = base_text_size, family = "sans")
)

df <- res$df
dfxmpl <- res$dfxmpl


pc <- ggplot(dfxmpl,aes(x=d,y=frac,fill=id))+
  geom_col(position="dodge")+
  scale_fill_discrete("novel\nkaryotype")+
  theme_classic(base_size=8)+
  labs(tag="D")+
  scale_x_discrete("distance from\nnovel karyotype")+
  scale_y_continuous("population fraction")+
  theme(legend.position = c(0.3,0.9),legend.key.size = unit(0.1,"in"))
pc


w <- split(df,f=df$dec_id)


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


# Load necessary libraries
library(ggplot2)
library(reshape2)
library(dplyr)
library(viridis)

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
  theme_classic(base_size=8)+
  labs(tag="A")+
  scale_x_continuous("",breaks=c(1,4),labels=c(expression(S[0],S[t])))+
  scale_y_continuous("clone frequency")
p0

library(ggforce)
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

# Behold some circles
pa <- ggplot() +
  geom_rect(aes(xmin = -2, xmax = 7, ymin = -3, ymax = 3),alpha=0,color="black")+
  geom_circle(aes(x0 = x0, y0 = y0, r = r), data = circles)+
  geom_text(data=anndf,aes(x=x,y=y,label=anno))+
  labs(tag="B")+
  theme_void(base_size=8)
pa

pb <- ggplot() +
  geom_rect(aes(xmin = -2, xmax = 7, ymin = -3, ymax = 3,fill="\u0398"),
            color="black", show.legend = F)+
  geom_circle(aes(x0 = x0, y0 = y0, r = r,fill=id), 
              data = circles,show.legend = F)+
  geom_text(data=anndf2,aes(x=x,y=y,label=anno),parse=F)+
  labs(tag="C")+
  theme_void(base_size=8)
pb



pd <- ggplot(z,aes(x=var,y=frac))+
  geom_col()+
  theme_classic(base_size=8)+
  labs(tag="E")+
  scale_x_discrete("variable predicting\nnovel karyotype",labels = scales::parse_format())+
  scale_y_continuous(paste0("fraction significant (p=",pthresh,")"))
pd


pe <- ggplot(w,aes(x=ids))+ 
  stat_count(aes(y=..count../sum(..count..)))+
  theme_classic(base_size=8)+
  labs(tag="F")+
  scale_x_discrete("variable predicting\nnovel karyotype",labels = scales::parse_format())+
  scale_y_continuous("fraction most significant")
pe

library(gridExtra)

plt <- grid.arrange(p0,grid.arrange(pa,pb,ncol=2),
                    grid.arrange(pc,pd,pe,ncol=3),
                    nrow=3,heights=c(4,3.5,4))
ggsave("figures/misc/figures/salehi_novel_kary.png",plot=plt,width=5,height=6,units="in")