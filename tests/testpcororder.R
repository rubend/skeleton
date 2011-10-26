setwd("/Users/rubendezeure/Documents/studies/hs_2011/semester_paper/skeleton")
library(inline)
library(pcalg)
pcorOrder2 <- cxxfunction(signature(i="integer",j="integer",k="integer",C="numeric"),body = paste(readLines("pcororderwrap.cpp"),collapse="\n"),includes=paste(readLines("graphfuncts.cpp"),collapse="\n"),plugin="RcppArmadillo")

myTrueCor <- function(g)
{
  ## Purpose: Compute true correlation matrix of graph
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - g: graph object
  ## ----------------------------------------------------------------------
  ## Value: True correlation matrix corresponding to graph
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 20 Sep 2011, 09:05
  a <- t(as(g, "matrix"))
  p <- ncol(a)
  m <- solve(diag(p) - a)
  cov2cor( m %*% t(m) )
}

makeGraph <- function(p, en, seed = NULL)
{
  ## Purpose: Generate a graph and produce corresponding skeleton and
  ## true covariance matrix of nodes
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - p: number of nodes (= gaussian random variable)
  ## - en: expected number of neighbors for each node
  ## - seed: seed for random number generator
  ## ----------------------------------------------------------------------
  ## Value: List with two elements
  ## - skel: adjacency matrix of skeleton as matrix (should be recovered)
  ## - corMat: True correlation matrix of all nodes (input for skeleton function)
  ## - g: generated graph object
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 20 Sep 2011, 08:36

s <- en/(p-1)
if (is.null(seed)) seed <- sample(10^4,1)
set.seed(seed)
rDAG <- randomDAG(p, prob = s, lB=0.1, uB=1)
amat <- as(rDAG, "matrix") + t(as(rDAG, "matrix"))
amat[amat != 0] <- 1
corMat <- myTrueCor(rDAG)

list(skel = amat, corMat = corMat, g = rDAG)
}


## TEST 1: Some random graphs of rather small size - true cor. mat
nreps <- 100
ok <- rep(NA, nreps)
p <- 7
en <- 3
i <- nreps
cat("i = ",i,"\n")
tmp <- makeGraph(p, en, seed = i) ## generate skeleton and true cor. matrix


## cases to test
## x y S
## 3 4 2
## 3 5 2
## 4 5 2

i <- 7
j <- 1
k <- c(2,3,4)
cat("i = ",i,"j = ",j,"\n")

t <- pcorOrder2(i,j,k,tmp$corMat)


cat("Correct pcor value = ",pcorOrder(i,j,k,tmp$corMat),"\n")

                                     
