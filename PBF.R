f <- function(x,a,b){
  x^(a-1)*(1-x)^(b-1)
}


ibeta <- function(x,a,b) # Incomplete beta function
{
  integrate(f,lower = 0,upper=x,a=a,b=b)$value
}

cbeta <- function(x,a,b)
{
  integrate(f,lower = x,upper=1,a=a,b=b)$value
}


prbf10 <- function(p,y,n){
  p1 <- ibeta(p,y+2,n-y+1)/ibeta(p,y+1,n-y+1)
  dbinom(y,n,prob = p1)/(dbinom(y,n,prob = p)*exp(1))
}

prbf20 <- function(p,y,n){
  p1 <- cbeta(p,y+2,n-y+1)/cbeta(p,y+1,n-y+1)
  dbinom(y,n,prob = p1)/(dbinom(y,n,prob = p)*exp(1))
}

pp <- function(x,n,p){ # posterior probabilities of Hypothesis
  l <- c(1,prbf10(p,x,n),prbf20(p,x,n))
  l/sum(l)
}

p <- 0.3
lower <- upper <- 0
lower.ex <- upper.ex <- 0
k <- 1000
for(i in 1:k){
  x <- 0:i
  y <- lapply(x,pp,n=i,p=p)
  y.lower <- y.upper <- 0
  for(j in x){
    y.lower[j+1] <- y[[j+1]][2]
    y.upper[j+1] <- y[[j+1]][3]
  }
  lower[i] <- which.max(y.lower<0.5)-2
  upper[i] <- which.max(y.upper>0.5)-1
  lower.ex[i] <- which.max(y.lower<0.95)-2
  upper.ex[i] <- which.max(y.upper>0.95)-1
}

plot(log10(1:k),lower/1:k,main = bquote("Consistency of PBF design"~phi==.(p)),
     xlab = bquote("number of subject treated at current dose"~(log[10])),
     ylab = "boundaries (rate)",col="coral2",cex=0.3,ylim = c(0,1),
     xaxt="n")
axis(1,at=c(0,1,2,3),labels = c(1,10,100,1000))
points(log10(1:k),upper/1:k,cex=0.3,col="coral2")

abline(h=p)

# BOIN 
#points(log10(1:1000),rep(0.236,1000),cex=0.3,col="darkgoldenrod2")
#points(log10(1:1000),rep(0.358,1000),cex=0.3,col="darkgoldenrod2")

points(log10(1:1000),rep(0.197,1000),cex=0.3,col="darkgoldenrod2")
points(log10(1:1000),rep(0.298,1000),cex=0.3,col="darkgoldenrod2")

# Keyboard
library(readxl)
#keyb <- read_excel("D:/study/Dissertation/Coding/Keyboard_v1.2.1.0_Decision Table.xlsx")
keyb <- read_excel("D:/study/Dissertation/Coding/Keyboard_v1.2.1.0_Decision Table - 0.25.xlsx")

points(log10(2:100),as.numeric(keyb[3,2:100])/as.numeric(keyb[2,2:100]),cex=0.3,col="royalblue")
points(log10(2:100),as.numeric(keyb[4,2:100])/as.numeric(keyb[2,2:100]),cex=0.3,col="royalblue")


legend(2.1, 1, legend=c("PBF design", "BOIN design", "Keyboard design","target toxicity rate"),
       col=c("coral2", "darkgoldenrod","royalblue","black"), lty=c(3,3,3,1), cex=0.6)


for(i in 1:10){
  x <- 0:(i*3)
  y <- lapply(x,pp,n=(i*3),p=p)
  y.lower <- y.upper <- 0
  for(j in x){
    y.lower[j+1] <- y[[j+1]][2]
    y.upper[j+1] <- y[[j+1]][3]
  }
  lower[i] <- which.max(y.lower<0.5)-2
  upper[i] <- which.max(y.upper>0.5)-1
  lower.ex[i] <- which.max(y.lower<0.95)-2
  upper.ex[i] <- which.max(y.upper>0.95)-1
}
