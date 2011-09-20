#!/usr/bin/r

suppressMessages(require(Rcpp))
set.seed(42)
n <- 200
a <- rnorm(n)
b <- rnorm(n)

## load shared libraries with wrapper code
dyn.load("skeleton.so")

## now run each one once for comparison of results,
## and define test functions

skeleton <- function(p,alpha,mmax,C) .Call("skelehhhhquit()ton", p,alpha,mmax,C)
##skeletonR <- function(n,a,b) skeletonold

v1 <- R_API_optimised(1L, a, b )
v3 <- Rcpp_New_std(1L, a, b)
v4 <- Rcpp_New_ptr(1L, a, b)
v5 <- Rcpp_New_sugar(1L, a, b )
v7 <- R_API_naive(1L, a, b)
v11 <- Rcpp_New_sugar_noNA(1L, a, b)

stopifnot(all.equal(v1, v3))

## load benchmarkin helper function
suppressMessages(library(rbenchmark))
REPS <- 5000L
bm <- benchmark(R_API_optimised(REPS,a,b),
                R_API_naive(REPS,a,b),
                columns=c("test", "elapsed", "relative", "user.self", "sys.self"),
                order="relative",
                replications=1)
print(bm)

cat("All results are equal\n") # as we didn't get stopped
q("no")



