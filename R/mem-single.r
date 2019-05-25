#source("utils.r")

#library(R.utils)

#' @title Fit the Exact MEM Model
#'
#' @description Fit the MEM model using full Bayesian inference.
#' @param responses the number of responses in each basket.
#' @param size the size of each basket.
#' @param name the name of each basket.
#' @param p0 the null response rate for the poster probability calculation
#' (default 0.15).
#' @param shape1 the first shape parameter(s) for the prior of each basket
#' (default 0.5).
#' @param shape2 the second shape parameter(s) for the prior of each basket
#' (default 0.5).
#' @param prior the matrix giving the prior inclusion probability
#' for each pair of baskets. The default is on on the main diagonal and 0.5
#' elsewhere.
#' @param hpd_alpha the highest posterior density trial significance.
#' @param alternative the alternative case definition (default greater)
#' @param call the call of the function (default NULL).
#' @importFrom stats rbinom
#' @examples
#' # 3 baskets, each with enrollement size 5
#' trial_sizes <- rep(5, 3)
#' 
#' # The response rates for the baskets.
#' resp_rate <- 0.15
#' 
#' # The trials: a column of the number of responses and a column of the
#' # the size of each trial.
#' trials <- data.frame(
#'   responses = rbinom(trial_sizes, trial_sizes, resp_rate),
#'   size = trial_sizes,
#'   name = letters[1:3]
#' )
#' 
#' summary(mem_exact(trials$responses, trials$size, trials$name))
#' @importFrom foreach foreach %dopar% getDoParName getDoSeqName registerDoSEQ
#' %do%
#' @importFrom stats median
#' @importFrom R.utilsB insert

#' @importFrom crayon red
#' @importFrom itertools isplitRows
#' @export
mem_single <- function(responses,
                      size,
                      name,
                      drug_index,
                      p0 = 0.15,
                      shape1 = 0.5,
                      shape2 = 0.5,
                      prior = rep(1.0, length(responses)) / 2,
                      hpd_alpha = 0.05,
                      alternative = "greater",
                      call = NULL) {
  h <- mod_i <- NULL
  # if (is.null(getDoParName())) {
  #   registerDoSEQ()
  # }
  
  if (length(responses) != length(size)) {
    stop(red(
      "The length of the responses and size parameters",
      "must be equal."
    ))
  }
  
  if (length(shape1) == 1) {
    shape1 <- rep(shape1, length(responses))
  }
  if (length(shape2) == 1) {
    shape2 <- rep(shape2, length(responses))
  }


  prior[drug_index] <- 1
  prior_inclusion <- prior
  alp <- hpd_alpha
  
  xvec <- responses
  nvec <- size
  
  k <- length(xvec) - 1
  
  m <- as.matrix(expand.grid(rep(list(c(
    0, 1
  )), k)))

  betaV <- beta(shape1, shape2)
  prod.vec <- beta(xvec + shape1, nvec + shape2 - xvec) / beta(shape1, shape2)
  log.Marg <- c()
  for(i in 1:dim(m)[1])
  {
    mvec <- insert(m[i,], drug_index, 1)
    ld <- logMarg.DensSingle(drug_index, mvec, xvec, nvec, shape1,shape2, betaV, prod.vec)
    log.Marg <- c(log.Marg, ld)
  }
  #H <- length(mod.mat)

  max.i <- order(log.Marg, decreasing = TRUE)[1]
  
  MAX <- insert(m[max.i,], drug_index, 1)
  names(MAX) <- name
  
  
  ## compute prior + posterior MEM probabilities ##
  
  PRIOR <-c()
  for(i in 1:dim(m)[1])
  {
    mvec <- insert(m[i,], drug_index, 1)
    pr <- mem.PriorSingle(mvec, prior_inclusion)
    PRIOR <- c(PRIOR, pr)
  }

  POST <- (exp(log.Marg) * PRIOR) / sum(exp(log.Marg) * PRIOR)
  map.i <- order(POST, decreasing = TRUE)[1]
  
  MAP <- insert(m[map.i,], drug_index, 1)
  names(MAP) <- name
  
  
  
  ## Posterior Exchangeability Prob ##

 
  PEP <- compPEPSingle(m, drug_index, log.Marg, PRIOR) 
  names(prior_inclusion) <- name
  names(PEP) <- name
  
  
  ## Marginal CDF and ESS ##
  # Here's where all of the compute goes.
  pweights <- post.Weights(m, drug_index, log.Marg, PRIOR) 

  vecM <- matrix(0, 0, dim(m)[2] + 1)
  for(i in 1:dim(m)[1])
  {
    mvec <- insert(m[i,], drug_index, 1)
    vecM <- rbind(vecM, mvec)
  }
  
  post.prob <-
    eval.Post(
      p0,
      xvec,
      nvec,
      vecM,
      pweights,
      shape1[drug_index],
      shape2[drug_index],
      alternative
    )

  
  pESS <-pweights %*% ESS(xvec, nvec, vecM, shape1[drug_index], shape2[drug_index])
  samples <- replicate(
    10000,
    samp.Post(xvec, nvec, vecM, pweights, shape1[drug_index], shape2[drug_index])
  )
  HPD <- boa.hpd(
    samples,
    alp
  )  
#browser()
  #names(pESS) <- names(post.prob) <- name
  #rownames(HPD) <- c("lower", "upper")
  #colnames(HPD) <- name
  models <- vecM
  

  
  if (missing(name)) {
    name <- paste("basket", seq_along(size))
  } else {
    if (is.factor(name)) {
      name <- as.character(name)
    }
    if (!is.character(name) ||
        length(name) != length(size)) {
      stop(red(
        "The basket name argument must be a character vector with",
        "one name per basket."))
    }
  }
  
  if (is.null(call)) {
    call <- match.call()
  }
  ret <-
    list(
      MLE = MAX,
      PRIOR = prior_inclusion,
      MAP = MAP,
      PEP = PEP,
      drug_index = drug_index,
      post.prob = post.prob,
      ESS = round(pESS, 2),
      HPD = HPD,
      samples = samples,
      responses = responses,
      size = size,
      name = name,
      p0 = p0,
      alpha = hpd_alpha,
      alternative = alternative,
      pweights = pweights,
      shape1 = shape1,
      shape2 = shape2,
      models = models,
      call = call
    )
  #browser()
  ret$mean_est <- mean(ret$samples)
  ret$median_est <- median(ret$samples)
  #browser()
  ret$ESS <-round(calc.ESS.from.HPD(fit = ret, alpha = ret$alpha), 2)
  #ret$ESS3 <-round(calc.ESS.from.HPDwid(fit = ret, alpha = ret$alpha), 2)
  prec <- 1.0 / var(ret$samples)
  #ret$ESS4 <- comp.ESS(prec, responses[drug_index], responses[drug_index] / size[drug_index])
  #ret$ESS4 <- comp.ESS(prec, responses[drug_index], mean(ret$samples))
  #ret$ESS4 <- round(betaESS(mean(ret$samples), var(ret$samples)), 2)
  
  
    #ret$samples <- sample_posterior_model(ret)
  #ret$mean_est <- colMeans(ret$samples)
  #ret$median_est <- apply(ret$samples, 2, median)


  class(ret) <- c("mem_single", "exchangeability_model")
  return(ret)
}

#' @importFrom ggplot2 ggplot aes geom_density scale_fill_manual facet_grid
#' xlab ylab theme_minimal xlim geom_vline labeller label_wrap_gen
#' @export
plot_borrow_density <- function(x, ...) {
  dots <- list(...)
  Basket <- Density <- p0 <- NULL
  # if (length(x$p0) == 1 && length(x$size) > 1) {
  #   x$p0 <- rep(x$p0, length(x$size))
  # }
  d <- gather(as_tibble(x$samples), key = Basket, value = Density)
  d$p0 <- x$p0
  
  if ("basket_levels" %in% names(dots)) {
    basket_levels <- dots$basket_levels
    d$Basket <- factor(d$Basket, levels = basket_levels)
  } else {
    d$Basket <- as.factor(d$Basket)
  }
  
  if ("basket_colors" %in% names(dots)) {
    basket_colors <- dots$basket_colors
  } else {
    basket_colors <- rep("black", length(x$responses))
  }
  #d$p0 <- x$p0[match(d$Basket, x$name)]
  d$basket_name <- paste0(" (p0=", d$p0, ")")
  ggplot(d, aes(x = Density, fill = Basket)) +
    geom_density(alpha = 0.7) +
    facet_grid( basket_name~ ., labeller = label_wrap_gen(width = 10)) +
    geom_vline(aes(xintercept = p0)) +
    scale_fill_manual(values = basket_colors, guide = FALSE) +
    xlab("") +
    ylab("Density") +
    xlim(0, 1) +
    theme_minimal()
}

