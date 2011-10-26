library(pcalg)

library(inline)
## test!
funskeleton <- cxxfunction(signature(pt="integer",alphat="numeric",m_maxt="integer",C="numeric",nt="long"),body=paste(readLines("skeletonquick.cpp"),collapse="\n"),includes=paste(readLines("graphfuncts.cpp"),collapse="\n"),plugin="RcppArmadillo")

## NOTE: Usually, we would ESTIMATE the correlation matrix; but since we are
## only interested in the runtime performance, we can use the true correlation
## matrix

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

estSkel <- function(corMat, n = 10^15, alpha = 0.05, verbose = FALSE,m_max = 100000)
{
  ## Purpose: Estimate the skeleton using rcpp function
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - corMat: Correlation Matrix of true graph
  ## - n: Number of samples that were used to estimate corMat
  ## (since we use the true cor. matrix, n = infinity; for numerical reasons,
  ## we use n = 10^15; this is a bit quick and dirty, but should work for
  ## almost all cases...)
  ## - alpha: Sign. level to test the presence of an edge
  ## ----------------------------------------------------------------------
  ## Value: Estimated adjacency matrix
  ## ----------------------------------------------------------------------
  ## Authors: Markus Kalisch, Date: 20 Sep 2011, 08:40
  ## Modified by Ruben Dezeure

  p <- ncol(corMat)
  ## define independence test (partial correlations)
  indepTest <- gaussCItest
  ## define sufficient statistics
  suffStat <- list(C = corMat, n = n)
  ##skeleton.fit <- skeleton(suffStat, indepTest, p, alpha, verbose = verbose)
  ##skeleton.fit <- funskeleton(pt=p,alphat=alpha,m_maxt=10000,C=corMat)
  ##for testing purposes
  ## IMPORTANT change m_maxt according to en
  cat(system.time(skeleton.fit <- funskeleton(pt=p,alphat=alpha,m_maxt=m_max,C=corMat,nt=n)))
  ##as(skeleton.fit@graph, "matrix")
  as(skeleton.fit, "matrix")
}
estSkelTime <- function(corMat, n = 10^15, alpha = 0.05, verbose = FALSE, m_max = 100000)
{
  ## Purpose: Estimate the skeleton using rcpp function
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - corMat: Correlation Matrix of true graph
  ## - n: Number of samples that were used to estimate corMat
  ## (since we use the true cor. matrix, n = infinity; for numerical reasons,
  ## we use n = 10^15; this is a bit quick and dirty, but should work for
  ## almost all cases...)
  ## - alpha: Sign. level to test the presence of an edge
  ## ----------------------------------------------------------------------
  ## Value: Estimated adjacency matrix
  ## ----------------------------------------------------------------------
  ## Authors: Markus Kalisch, Date: 20 Sep 2011, 08:40
  ## Modified by Ruben Dezeure

  p <- ncol(corMat)
  ## define independence test (partial correlations)
  indepTest <- gaussCItest
  ## define sufficient statistics
  suffStat <- list(C = corMat, n = n)
  ##skeleton.fit <- skeleton(suffStat, indepTest, p, alpha, verbose = verbose)
  ##skeleton.fit <- funskeleton(pt=p,alphat=alpha,m_maxt=10000,C=corMat)
  ##for testing purposes
  ## IMPORTANT change m_maxt according to en
  time <- system.time(skeleton.fit <- funskeleton(pt=p,alphat=alpha,m_maxt=m_max,C=corMat,nt=n))
  time[3]
}

estSkel2 <- function(corMat, n = 10^15, alpha = 0.05, verbose = FALSE,m_max = Inf)
{
  ## Purpose: Estimate the skeleton using R algorithm
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - corMat: Correlation Matrix of true graph
  ## - n: Number of samples that were used to estimate corMat
  ## (since we use the true cor. matrix, n = infinity; for numerical reasons,
  ## we use n = 10^15; this is a bit quick and dirty, but should work for
  ## almost all cases...)
  ## - alpha: Sign. level to test the presence of an edge
  ## ----------------------------------------------------------------------
  ## Value: Estimated adjacency matrix
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 20 Sep 2011, 08:40

  p <- ncol(corMat)
  ## define independence test (partial correlations)
  indepTest <- gaussCItest
  ## define sufficient statistics
  suffStat <- list(C = corMat, n = n)
  cat(system.time(skeleton.fit <- skeleton(suffStat, indepTest, p, alpha, verbose = verbose,m.max=m_max)))
  as(skeleton.fit@graph, "matrix")
 }
estSkel2time <- function(corMat, n = 10^15, alpha = 0.05, verbose = FALSE,m_max = Inf)
{
  ## Purpose: Estimate the skeleton using R algorithm
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - corMat: Correlation Matrix of true graph
  ## - n: Number of samples that were used to estimate corMat
  ## (since we use the true cor. matrix, n = infinity; for numerical reasons,
  ## we use n = 10^15; this is a bit quick and dirty, but should work for
  ## almost all cases...)
  ## - alpha: Sign. level to test the presence of an edge
  ## ----------------------------------------------------------------------
  ## Value: Estimated adjacency matrix
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 20 Sep 2011, 08:40

  p <- ncol(corMat)
  ## define independence test (partial correlations)
  indepTest <- gaussCItest
  ## define sufficient statistics
  suffStat <- list(C = corMat, n = n)
  time <- system.time(skeleton.fit <- skeleton(suffStat, indepTest, p, alpha, verbose = verbose,m.max=m_max))
  time[3]
}
