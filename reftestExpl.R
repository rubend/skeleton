setwd("/Users/rubendezeure/Documents/studies/hs_2011/semester_paper/skeleton")
source('prepareruntime.r') ## get everything loaded and compiled

#debug test cases

## TEST 0: Simple test case
##nreps <- 100
##ok <- rep(NA, nreps)
##p <- 5
##en <- 2
##i <- 1
##cat("i = ",i,"\n")
##tmp <- makeGraph(p, en, seed = i) ## generate skeleton and true cor. matrix
##res <- estSkel2(tmp$corMat,verbose=TRUE) ## estimate skel perfectly (b/c true cor. mat. used)
##res == tmp$skel
##res2 <- estSkel(tmp$corMat) ## estimate skel perfectly (b/c true cor. mat. used)
##res2 == tmp$skel ## COMPARE YOUR SOLUTION WITH TRUTH (tmp$amat)
##res == res2


## timer performance test
p <- 20
en <- 3
tmp <- makeGraph(p, en, seed = 17)
cat("Time running R code")
res2 <- estSkel2(tmp$corMat) ## estimate skel perfectly (b/c true cor. mat. used)
cat("Time running cpp code")
res <- estSkel(tmp$corMat) ## estimate skel perfectly (b/c true cor. mat. used)

all(res2 == res)

##TEST 1-4 see email