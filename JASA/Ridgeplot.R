rm(list=ls())

# Libraries
library(ggridges)
library(ggplot2)
library(viridis)
library(hrbrthemes)

# Data

source("measure.R")

# Import R workspaces
n <- 36
K <- 6
cs <- 1
phi <- 0.20
phi.list <- c(0.2,0.25,0.3)
K.list <- c(4,6,7,9)
n.list <- c(30,36,60,96)
cs.list <- c(1,3)
method.list <- c("pop2","boin","kb")

metrics <- "pcs"

dt <- data.frame(K=integer(),n=integer(),phi=double(),cohortsize=integer(),method=character(),y=double())
for(k in 1:length(K.list)){
  K <- K.list[k]
  for(phi in phi.list){
    for(cs in cs.list){
      for(m in method.list){
        load(paste0(m,"-",phi*100,"-",K,"-",cs,".RData"))
        if(m==method.list[1]){
          mtd <- rep(0,n.scene)
          tol.mtd <- matrix(nrow = n.scene, ncol = K)
          for(i in 1:n.scene){
            mtd[i] <- unname(which.min(abs(skeleton_list[i,]-phi)))
            for(j in 1:K){
              tol.mtd[i,j] <- ifelse(skeleton_list[i,j]>=phi-0.05 & skeleton_list[i,j]<=phi+0.05,1,0)
            }
          }
          ref <- measure(design.mtd,design.allo,design)
          K <- K.list[k]
        }else{
          out <- measure(design.mtd,design.allo,design)
          K <- K.list[k]
          temp <- data.frame(K=K,n=n,phi=phi,cohortsize=cs,method=m,y=out$pcs-ref$pcs)
          dt <- rbind(dt,temp)
        }
      }
    }
  }
}
dt$sign <- as.factor(sign(dt$y))

dt2 <- dt#[dt$K %in% c(4,6,9),]

# New facet label names for cohort size
cs.labs <- c("cohort size = 1", "cohort size = 3")
names(cs.labs) <- c("1", "3")

# New facet label names for dose levels
dose.labs <- c("K=4, N=30", "K=6, N=36", "K=6, N=60", "K=8, N=48", 
               "K=8, N=96", "K=10, N=96")
names(dose.labs) <- c("4", "6", "7", "8", "9", "10")

# New facet label names for phi
phi.labs <- c("\u03a6=0.2","\u03a6=0.25","\u03a6=0.3")
names(phi.labs) <- c("0.2", "0.25","0.3")

# p1 <- ggplot(dt, aes(x = y, y = method, group=method, fill = method)) +
#   geom_density_ridges_gradient(scale = 1.5, rel_min_height = 0.001,
#                                quantile_lines=TRUE,
#                                quantile_fun=function(x,...) 0) +
#   xlim(-.12,.25)+scale_fill_brewer(palette = 4)+
#   #scale_fill_viridis_d(direction = -1, guide = "none")+
#   #scale_fill_viridis(name = "PCS", option = "A",begin = 0.2,end = 0.8) +
#   #labs(title = 'Temperatures in Lincoln NE in 2016') +
#   theme_ipsum() +
#   theme(
#     legend.position="bottom",
#     panel.spacing = unit(0.8, "lines"),
#     strip.text.x = element_text(hjust=0.5),
#     strip.text.y.left = element_text(angle = 0),
#     axis.text.y = element_blank(),
#     axis.title.y = element_blank()
#   )+
#   facet_grid(rows = vars(K),cols = vars(cohortsize), switch = "y",
#              labeller = labeller(K = dose.labs, cohortsize = cs.labs))
# p1

p2 <- ggplot(dt2, aes(x = -y, y = method, group=method, fill = method)) +
  geom_density_ridges(scale = 1.5, 
                      #stat = "density",
                      rel_min_height = 0.002,alpha=0.8,
                               quantile_lines=TRUE,
                               quantile_fun=function(x,...) 0) +
  xlim(-.12,.23)+
  #ylim(0,0.2)+
  #scale_fill_brewer(palette = "set1")+
  #geom_hline(yintercept = Inf, color = "white", size = 4) +
  scale_fill_manual(values=c("#505EA1",
                             "#86A8E7"))+
  theme_ridges() +
  theme(
    legend.position="bottom",
    panel.spacing = unit(0.5, "lines"),
    strip.text.x = element_text(hjust=0.5,vjust = 1.2),
    strip.text.y.left = element_text(angle = 0),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    strip.background.x = element_rect(color = "white", size = 8, fill = NA),
    strip.background.y = element_blank(),
    strip.placement = "outside"
  )+
  facet_grid(phi+cohortsize~K, switch = "y",
             labeller = labeller(K = dose.labs, cohortsize = cs.labs, phi=phi.labs))

dat_text <- data.frame(
  label = rep(NA,24),
  K = rep(unique(dt2$K),6),
  n = rep(unique(dt2$n),6),
  phi = rep(unique(dt2$phi),each=8),
  cohortsize = rep(c(1,3,1,3,1,3),each=4)
)
dat_text <- rbind(dat_text,dat_text)
dat_text$method <- rep(c("boin","kb"),each=24)

for(i in 1:dim(dat_text)[1]){
  temp <- -dt2$y[dt2$K==dat_text$K[i] & dt2$phi==dat_text$phi[i] & dt2$cohortsize==dat_text$cohortsize[i] & dt2$method==dat_text$method[i]]
  x <- mean(temp>0)
  dat_text$label[i] <- paste(round(100*x, 0), "%", sep="")
}

p3 <- p2+  annotate(
  geom = "curve", x = 0.05, y = 1.1, xend = 0.13, yend = 1.5, 
  curvature = .3, arrow = arrow(length = unit(2, "mm"))
 ) + annotate(
   geom = "curve", x = 0.05, y = 2.2, xend = 0.13, yend = 2.5, 
   curvature = .3, arrow = arrow(length = unit(2, "mm")),
   color = "#5D9693"
 )+
  geom_text(
  data    = dat_text,
  mapping = aes(x = 0.20, y = method, label = label, color=method),
  vjust = -2
)+scale_color_manual(values = c("black","#5D9693"))
p3


library(gtable)
library(grid)

# function to remove selected elements from gtables, keeping widths
gtable_grob_remove <- function (g, what = "guide-box") {
  require(gtable)
  matches <- c(grepl(pattern = what, g$layout$name))
  g$layout <- g$layout[!matches, , drop = FALSE]
  g$grobs <- g$grobs[!matches]
  return(g)
}

# Create gtable object
gt <- ggplot_gtable(ggplot_build(p3))
print(gt)

# Get the panels
gl.panel <- gtable_filter(gt, "panel")
panels <-c(subset(gt$layout, grepl("panel", gt$layout$name), se=t:r))

# Get the x.strips
gl.strip.x <- gtable_filter(gt,"strip-t")
strips.x <- c(subset(gt$layout, grepl("strip-t", gt$layout$name), se=t:r))

strips.y <- c(subset(gt$layout, grepl("strip-l", gt$layout$name), se=t:r))

# Get the legend


gt$widths[4] <- unit(3,"cm")
a <- textGrob("  \u03a6       cohort")
tmp <- gtable_grob_remove(gt,what = "strip-l")
tmp <- gtable_add_grob(tmp,a,t=6,l=4,name = "strip.y")

strip.y.list <- paste(c("0.2  ","   ","0.25","   ","0.3  ","   "),c("1","    3"),sep = "         ")
for(i in 1:length(strips.y[[1]])){
  a <- textGrob(strip.y.list[i])
  tmp <- gtable_add_grob(tmp,a,t=strips.y$t[i],l=strips.y$l[i],name=paste0("strip-l-",i))
}

tmp <- gtable_grob_remove(tmp)
tmp <- gtable_add_grob(tmp,gt$grobs[[which(gt$layout$name=="guide-box")]],t=23,l=9,name="guide-box")
plot(tmp)

grid.newpage()
#tmp <- gtable_grob_remove(tmp)
grid.draw(tmp)
