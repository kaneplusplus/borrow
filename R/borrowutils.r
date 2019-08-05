
betaESS <- function(meanV, varV)
{
  ess <- meanV * (1 - meanV) / varV - 1
  ess
}

boa.hpd <- function(x, alpha) {
  n <- length(x)
  m <- max(1, ceiling(alpha * n))
  y <- sort(x)
  a <- y[1:m]
  b <- y[(n - m + 1):n]
  i <- order(b - a)[1]
  structure(c(a[i], b[i]), names = c("Lower Bound", "Upper Bound"))
}

logMarg.DensSingle <- function(i, mvec, xvec, nvec, avec, bvec, betaV, prod.vec) {
  
  p.vec <- prod(prod.vec ^ (1 - mvec))
  marg.vec <-
    (beta(avec[i] + mvec %*% xvec, bvec[i] + mvec %*% (nvec - xvec)) /
       betaV[i]) * p.vec
  
  sum(log(marg.vec))
}

mem.PriorSingle <- function(mvec, pr.Inclus) {
  s.in <- pr.Inclus # source inclusion probability
  s.ex <- 1 - pr.Inclus # source exclusion probability
  prod(s.in^(mvec) * s.ex^(1 - mvec))
}


i.pep <- function(v, log.Marg, PRIOR) {
  ((exp(log.Marg) * PRIOR) %*% v) / sum(exp(log.Marg) * PRIOR)
}

compPEPSingle <- function(m, drug_index, log.Marg, PRIOR) {
  vecM <- matrix(0, 0, dim(m)[2] + 1)
  for(i in 1:dim(m)[1])
  {
    mvec <- insert(m[i,], drug_index, 1)
    vecM <- rbind(vecM, mvec)
  }
  apply(vecM, FUN = i.pep, MARGIN = 2, log.Marg, PRIOR)
}


post.Weights <- function(m, drug_index, log.Marg, PRIOR) {
  (exp(log.Marg) * PRIOR) / sum(exp(log.Marg) * PRIOR)
}

#' @importFrom stats pbeta
eval.Post <- function(p0, X, N, Omega, w, a, b, alternative = "greater") {
  alph <- a + Omega %*% X
  beta <- b + (Omega %*% N - Omega %*% X)
  if (alternative == "greater") {
    out <- sum((1 - pbeta(p0, alph, beta)) * w)
  } else {
    out <- sum((pbeta(p0, alph, beta)) * w)
  }
  return(out)
}

ESS <- function(X, N, Omega, a, b) {
  # a <- b <- 0.5
  alph <- a + Omega %*% X
  beta <- b + (Omega %*% N - Omega %*% X)
  alph + beta
}

#' @importFrom stats rbeta
gen.Post <- function(X, N, Omega, a, b) {
  alph <- a + Omega %*% X
  beta <- b + (Omega %*% N - Omega %*% X)
  return(rbeta(1, alph, beta))
}

#' @importFrom stats rmultinom
samp.Post <- function(X, N, Omega, w, a, b) {
  return(gen.Post(X, N, Omega[which(rmultinom(1, 1, w) == 1), ], a, b))
}

euc.dist <- function(x1, x2, w = c(1, 1)) {
  if (sum(is.na(x1)) > 1) {
    Inf
  } else {
    sqrt(sum(w * ((x1 - x2)^2)))
  }
}

#' @importFrom stats qbeta
dist.beta.HPD <- function(ess, fit, alpha, jj) {
  al <- fit$mean_est * ess
  al <- max(1e-2, al)
  be <- ess - al
  be <- max(1e-2, be)
  euc.dist(fit$HPD, qbeta(
    c(alpha / 2, 1 - alpha / 2),
    al, be
  ))
}

#' @importFrom GenSA GenSA
ESS.from.HPD.i <- function(jj, fit, alpha) {
  # library(GenSA)
  #browser()
  opt <-
    GenSA(
      par = 1,
      fn = dist.beta.HPD,
      lower = 0,
      upper = 10000000,
      # control=list(maxit=pars$DTW.maxit),
      fit = fit,
      alpha = alpha,
      jj = jj
    )
  opt$par
}

#' @importFrom foreach %dopar%
calc.ESS.from.HPD <- function(fit, alpha) {
  ## fit is list with median vec and HPD vec ##
  i <- fit$drug_index
  #res <-c()
  #foreach(i in 1:length(fit$mean_est))# = seq_along(fit$mean_est), .combine = c) %dopar% {
  ESS.from.HPD.i(i, fit, alpha)
  
}

## ESS3


#' @importFrom stats qbeta
dist.beta.HPDwid <- function(ess, fit, alpha, jj) {
  al <- fit$median_est * ess
  al <- max(1e-2, al)
  be <- ess - al
  be <- max(1e-2, be)
  return(abs((fit$HPD[2] - fit$HPD[1]) - (diff(
    qbeta(
      c(alpha / 2, 1 - alpha / 2),
      al, be
    )
  ))))
}





ESS.from.HPDwid.i <- function(jj, fit, alpha) {
  # library(GenSA)
  opt <-
    GenSA::GenSA(
      par = 1,
      fn = dist.beta.HPDwid,
      lower = 0,
      upper = 10000000,
      # control=list(maxit=pars$DTW.maxit),
      fit = fit,
      alpha = alpha,
      jj = jj
    )
  return(opt$par)
}



calc.ESS.from.HPDwid <- function(fit, alpha) {
  ## fit is list with median vec and HPD vec ##
  
  i <- fit$drug_index
  #res <-c()
  #foreach(i in 1:length(fit$mean_est))# = seq_along(fit$mean_est), .combine = c) %dopar% {
  ESS.from.HPDwid.i(i, fit, alpha)
}
