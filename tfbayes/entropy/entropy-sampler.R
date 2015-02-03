
library(logitnorm)
library(gtools)
library(KernSmooth)

################################################################################

source("entropy.R")

# simplex sampler for a general density f
################################################################################

simplex.sampler.dir <- function(n, f, k) {
  alpha <- rep(1, k)
  rdir  <- function( ) c(rdirichlet(1, alpha))
  ddir  <- function(x)   ddirichlet(x, alpha)
  # initial state
  result <- rbind(NULL, rdir())
  for (i in 1:n) {
    q <- result[i,]
    # draw a proposal
    p <- rdir()
    # accept or reject
    if (runif(1) <= min(1, f(p)/f(q)*ddir(q)/ddir(p))) {
      result <- rbind(result, p)
    }
    else {
      result <- rbind(result, q)
      row.names(result)[i+1] <- "x"
    }
  }
  result
}

simplex.sampler.norm <- function(n, f, k, sd=0.05, ...) {
  stopifnot(k > 1)
  rprop  <- function(qj) rnorm(1, mean=qj, sd=sd)
  # initial state
  result <- rbind(NULL, rdirichlet(1, rep(1, k)))
  for (i in 1:n) {
    q <- result[i,]
    r <- TRUE
    for (j1 in 1:k) {
      # copy the current state
      p  <- q
      # draw the second coordinate
      j2 <- sample((1:k)[-j1], 1)
      # compute the range
      ra <- p[j1] + p[j2]
      # draw a proposal
      p[j1] <- (q[j1] + ra*rprop(0)) %% ra
      p[j2] <- 1-sum(p[-j2])
      # accept or reject
      if (f(q) == 0 || runif(1) <= min(1, f(p)/f(q))) {
        q <- p
        r <- FALSE
      }
    }
    # propose a permutation
    p <- gtools::permute(q)
    # accept or reject
    ## if (f(q) == 0 || runif(1) <= min(1, f(p)/f(q))) {
    ##   q <- p
    ## }
    # save result
    result <- rbind(result, q)
    # mark if sample was rejected
    if (r) {
      row.names(result)[i+1] <- "x"
    }
    else {
      row.names(result)[i+1] <- ""
    }
  }
  result
}

simplex.sampler.ham.grad <- function(q) {
  n <- length(q)
  g <- rep(0.0, n)
  for (i in 1:(n-1)) {
    g[i] <- log(q[n]) - log(q[i])
  }
  g
}

simplex.sampler.ham <- function(n, f, k, epsilon=0.01, L=10, ...) {
  stopifnot(k > 1)
  # initial state
  result <- rbind(NULL, rdirichlet(1, rep(1, k)))
  for (i in 1:n) {
    current_q <- result[i,]
    current_p <- c(rnorm(k-1, 0, 1), 0)
    q <- current_q
    p <- current_p
    r <- TRUE
    # make a half step for the momentum
    p <- p - epsilon*simplex.sampler.ham.grad(q)
    for (l in 1:L) {
      # make a full step for the position
      q <- q + epsilon*p
      # renormalize
      q[k] <- 1-sum(q[-k])
      # check if still within the probability simplex
      if (q[k] < 0.0) {
        continue
      }
    }
    # accept or reject
    if (f(q) == 0 || runif(1) <= min(1, f(p)/f(q))) {
      q <- p
      r <- FALSE
    }
    # propose a permutation
    p <- gtools::permute(q)
    # save result
    result <- rbind(result, q)
    # mark if sample was rejected
    if (r) {
      row.names(result)[i+1] <- "x"
    }
    else {
      row.names(result)[i+1] <- ""
    }
  }
  result
}

simplex.sampler <- function(n, f, k, burnin=1000, method="normal", ...) {
  result <- switch(method,
                   normal    = simplex.sampler.norm(n+burnin, f, k, ...),
                   hamilton  = simplex.sampler.ham (n+burnin, f, k, ...),
                   dirichlet = simplex.sampler.dir (n+burnin, f, k, ...))
  result <- result[(burnin+1):(n+burnin),]
  print (sprintf("Acceptance rate: %f", sum(row.names(result) != "x")/n))
  result
}
