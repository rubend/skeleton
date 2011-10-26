getNextSet <- function(n,k,set){
  chInd <- k-(zeros <- sum((seq(n-k+1,n)-set) ==0))
  wasLast <- (chInd==0)
  if(!wasLast){
    
    set[chInd] <- set[chInd]+1
    if(chInd<k){
      set[(chInd+1):k] <- seq(set[chInd]+1,set[chInd]+zeros)
    }
  }
  list(nextSet=set,wasLast=wasLast)
}

pcorOrder <- function(i,j,k,C,cut.at=0.9999999){
  if(length(k) ==0){
    r <- C[i,j]
  }else if(length(k)==1){
    r <- (C[i,j]-C[i,k]*C[j,k])/sqrt((1-C[j,k]^2)*(1-C[i,k]^2))
  }else{
    stopifnot(require("corpcor",quietly=TRUE))
    PM <- pseudoinverse(C[c(i,j,k),c(i,j,k)])
    r <- -PM[1,2]/sqrt(PM[1,1]*PM[2,2])
  }
  if(is.na(r))
    r <- 0
  min(cut.at,max(-cut.at,r))
}
