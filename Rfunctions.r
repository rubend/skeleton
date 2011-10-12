##################################################
## function for computing the skeleton
## see help-file in pcalg package for more details on arguments
##################################################
skeleton <- function(suffStat, indepTest, p, alpha, verbose = FALSE,
                     fixedGaps = NULL, fixedEdges = NULL, NAdelete = TRUE,
                     m.max = Inf) {

  ## Purpose: Perform undirected part of PC-Algorithm, i.e.,
  ## estimate skeleton of DAG given data
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 09.12.2009

  ## x,y,S konstruieren
  ##-   tst <- try(indepTest(x,y,S, obj))
  ##-   if(inherits(tst, "try-error"))
  ##-     stop("the 'indepTest' function does not work correctly with 'obj'")

  stopifnot((p <- as.integer(p)) >= 2)
  cl <- match.call()
  ## start skeleton

  ## fixed gaps
  if (is.null(fixedGaps)) {
    ## G := complete graph :
    G <- matrix(TRUE, p,p)
    diag(G) <- FALSE
  } else {
    if (!identical(dim(fixedGaps),c(p,p))) {
      stop("Dimensions of the dataset and fixedGaps do not agree.")
    } else {
      if (!all(fixedGaps == t(fixedGaps)))
        stop("fixedGaps must be symmetric")
      G <- !fixedGaps
    }
  } ## if(is.null(G))

  ## fixed edges
  if (is.null(fixedEdges)) {
    fixedEdges <- matrix(FALSE, p,p)
  } else {
    if (!(identical(dim(fixedEdges),c(p,p))))
      stop("Dimensions of the dataset and fixedEdges do not agree.")
    if (fixedEdges != t(fixedEdges))
      stop("fixedEdges must be symmetric")
  }

  seq_p <- seq_len(p)
  sepset <- pl <- vector("list",p)
  for (i in seq_p) sepset[[i]] <- pl
  ## save maximal p value
  pMax <- matrix(-Inf, p,p)
  diag(pMax) <- 1

  done <- FALSE
  ord <- 0L
  n.edgetests <- numeric(1)# final length = max { ord}

  while (!done && any(G) && ord <= m.max) {
    n.edgetests[ord1 <- ord+1L] <- 0
    done <- TRUE
    ind <- which(G, arr.ind = TRUE)

    ind <- ind[order(ind[,1]) ,]
    remainingEdgeTests <- nrow(ind)
    if(verbose)
      cat("Order=",ord,"; remaining edges:",remainingEdgeTests,"\n", sep='')
    for (i in 1:remainingEdgeTests) {
      if(verbose) { if(i%%100==0) cat("|i=",i,"|iMax=",nrow(ind),"\n") }
      x <- ind[i,1]
      y <- ind[i,2]
      if (G[y,x] && !fixedEdges[y,x]) {
        nbrsBool <- G[,x]
        nbrsBool[y] <- FALSE
        nbrs <- seq_p[nbrsBool]
        length_nbrs <- length(nbrs)
        if (length_nbrs >= ord) {
          if (length_nbrs > ord) done <- FALSE
          S <- seq_len(ord)
          repeat { ## condition w.r.to all  nbrs[S] of size 'ord'
            n.edgetests[ord1] <- n.edgetests[ord1]+1
            pval <- indepTest(x,y, nbrs[S], suffStat)
            ## pval <- dsepTest(x,y,nbrs[S],gTrue,jp = jp)
            if (verbose) cat("x=",x," y=",y," S=",nbrs[S],": pval =",pval,"\n")
            if (is.na(pval)) pval <- if(NAdelete) 1 else 0
            if (pval > pMax[x,y]) pMax[x,y] <- pval
            
            if(pval >= alpha) { # independent
              G[x,y] <- G[y,x] <- FALSE
              sepset[[x]][[y]] <- nbrs[S]
              break
            } else {
              nextSet <- getNextSet(length_nbrs, ord, S)
              if(nextSet$wasLast)
                break
              S <- nextSet$nextSet
            }
          } ## repeat
        } ## if (length_nbrs >= ord)
      } ## if(!done)

    } ## for(i in 1:remainingEdgeTests)
    G
    ord <- ord1
  } ## while

  for (i in 1:(p-1)) {
    for (j in 2:p) {
      pMax[i,j] <- pMax[j,i] <- max(pMax[i,j],pMax[j,i])
    } ## for (j in 2:p)
  } ## for (i in 1:(p-1))

  ## transform matrix to graph object :
  nnms <- as.character(seq_p)
  Gobject <-
    if (sum(G) == 0) {
      new("graphNEL", nodes = nnms)
    } else {
      colnames(G) <- rownames(G) <- nnms
      as(G,"graphNEL")
    }

  ## final object
  new("pcAlgo",
      graph = Gobject,
      call = cl, n = integer(0), max.ord = as.integer(ord-1),
      n.edgetests = n.edgetests, sepset = sepset,
      pMax = pMax, zMin = matrix(NA,1,1))

}## end{ skeleton }

##################################################
## "indepTest" coubld be a function like "pcorOrder"
##################################################
pcorOrder <- function(i,j, k, C, cut.at = 0.9999999) {
  ## Purpose: Compute partial correlation
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - i,j,k: Partial correlation of i and j given k
  ## - C: Correlation matrix among nodes
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 26 Jan 2006; Martin Maechler
  if (length(k)==0) {
    r <- C[i,j]
  } else if (length(k)==1) {
    r <- (C[i, j] - C[i, k] * C[j, k])/sqrt((1 - C[j, k]^2) * (1 - C[i, k]^2))
  } else { ## length(k) >= 2
    ## If corpcor was decent and had a name space, we'd just use
    ## PM <- corpcor::pseudoinverse(C[c(i,j,k), c(i,j,k)])
    stopifnot(require("corpcor", quietly=TRUE))
    ## pseudoinverse to make sure no singularity problem comes up
    ## I guess, this is not optimal...
    PM <- pseudoinverse(C[c(i,j,k), c(i,j,k)]) 
    r <- -PM[1, 2]/sqrt(PM[1, 1] * PM[2, 2])
  }
  if(is.na(r)) r <- 0
  min(cut.at, max(-cut.at, r))
}

##################################################
## Compute tuples on the fly
##################################################
getNextSet <- function(n,k,set) {
  ## Purpose: Generate the next set in a list of all possible sets of size
  ##          k out of 1:n;
  ##  Also returns a boolean whether this set was the last in the list.
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - n,k: Choose a set of size k out of numbers 1:n
  ## - set: previous set in list
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 26 Jan 2006, 17:37

  ## chInd := changing Index
  chInd <- k - (zeros <- sum((seq(n-k+1,n)-set) == 0))
  wasLast <- (chInd == 0)
  if (!wasLast) {
    set[chInd] <- set[chInd] + 1
    if (chInd < k)
      set[(chInd+1):k] <- seq(set[chInd]+1, set[chInd]+zeros)
  }
  list(nextSet = set, wasLast = wasLast)
}
