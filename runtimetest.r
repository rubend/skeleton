## need to be in correct dir to be able to compile the cpp files!
## setwd("/Users/rubendezeure/Documents/studies/hs_2011/semester_paper/skeleton")
## in my case on laptop
source('prepareruntime.r') ## get everything loaded and compiled
## testing runtimes in a structured way
en <- c(3,5)
p <- c(50,100,500,1000,5000)
alpha <- 0.01
nrep <- 5

## for each combination of en and p do 5 replicates (if p=5000 takes too long, just drop it)
## record the runtimes of both methods and the ratio (R/C++)
## finally we look at log(ratio) = a + b*log(p) + c*en
timecpp <- timeR <- matrix(numeric(length(en)*length(p)),length(en),length(p))

for (i in 1:length(en)){
  for (j in 1:length(p)){
    for (k in 1:nrep){
      tmp <- makeGraph(p[j], en[i], seed = k) ## can take a while :/
      timecpp[i,j] <- estSkeltime(tmp$corMat, alpha)/nrep
      timeR[i,j] <- estSkel2time(tmp$corMat,alpha)/nrep
      ## divide by nrep so we get averages
    }
  }
}
ratio <- timeR/timecpp
##look at log(ratio) = a + b*log(p) + c*en

