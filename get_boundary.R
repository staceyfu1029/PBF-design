######### Calculate decision boundary of PBF design 


## Calculate the PrBF
prbf01 <- function(n,y,phi){
  p <- beta(y+2,n-y+1)/beta(1+y,n-y+1)
  dbinom(y,n,prob=phi)*exp(1)/dbinom(y,n,prob=p)
}

bound <- function(n.cohort,cohortsize,phi){
  lower <- upper <- 0
  for(i in 1:n.cohort){
    t <- i*cohortsize
    x <- 0:(i*cohortsize)
    y <- lapply(x,prbf01,n=(i*cohortsize),phi=phi)
    a <- 0
    for(j in 1:length(x)){
      if(y[[j]]<exp(1)){ #e,exp(1)*(1/t*log(1+t))^(1/(2*t)),exp(1)*(1/t*log(1+t))^(1/(t))
        if(x[j]/(i*cohortsize)<phi){
          a[j] <- 1
        }else{
          a[j] <- -1
        }
      }else{
        a[j] <- 0
      }
    }
    lower[i] <- max(which(a==1))-1
    upper[i] <- min(which(a==-1))-1
  }
  b <- rbind(lower,upper)
}

#### Boundaries
bound(n.cohort,cohortsize,phi)
