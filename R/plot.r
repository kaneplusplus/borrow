#' @title Plot the posterior density of the fitting result
#'
#' @description Create the posterior density plot of the studied basket.
#' @param x the fitted model calculated from the borrow_single() and borrow_multiple() function.
#' @param ... other parameters
#' @importFrom ggplot2 ggplot 
#' @export
plot_borrow_density <- function(x, ...) {
  UseMethod("plot_borrow_density")
}

#' @importFrom crayon red
#' @export
plot_borrow_density.default <- function(x, ...) {
  stop(red(
    "Don't know how to make a density plot with an object of type",
    paste(class(x), collapse = ", "), "."))
}

#' @export
plot_borrow_density.borrow_single<- function(x, ...) {
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

#' @export
plot_borrow_density.borrow_multiple<- function(x, ...) {
  dots <- list(...)
  Basket <- Density <- p0 <- NULL
  d <- x$data
  allS <- d[[1]]$name[d[[1]]$drug_index]
  samp <- d[[1]]$samples
  if(length(d) > 1){
    for (i in 2:length(d)){
      samp <- cbind(samp, d[[i]]$samples)
      allS <- c(allS, d[[i]]$name[d[[i]]$drug_index])
    }
  }
  colnames(samp) <- allS
  # if (length(x$p0) == 1 && length(x$size) > 1) {
  #   x$p0 <- rep(x$p0, length(x$size))
  # }
  p0 <- d[[1]]$p0
  
  if (length(p0) == 1 && length(allS) > 1) {
    p0 <- rep(p0, length(allS))
  }
  d <- gather(as_tibble(samp), key = Basket, value = Density)
  d$p0 <- p0
  
  if ("basket_levels" %in% names(dots)) {
    basket_levels <- dots$basket_levels
    d$Basket <- factor(d$Basket, levels = basket_levels)
    d$Basket <- factor(d$Basket, levels = basket_levels)
  } else {
    d$Basket <- as.factor(d$Basket)
  }
  
  if ("basket_colors" %in% names(dots)) {
    basket_colors <- dots$basket_colors
  } else {
    basket_colors <- rep("black", length(x$data))
  }
  #d$p0 <- x$p0[match(d$Basket, x$name)]
  d$p0 <- d$p0[match(d$Basket, allS)]
  d$basket_name <- paste0(d$Basket, " (p0=", d$p0, ")")
  ggplot(d, aes(x = Density, fill = Basket)) +
    geom_density(alpha = 0.7) +
    facet_grid( basket_name ~ ., labeller = label_wrap_gen(width = 10)) +
    geom_vline(aes(xintercept = p0)) +
    scale_fill_manual(values = basket_colors, guide = FALSE) +
    xlab("") +
    ylab("Density") +
    xlim(0, 1) +
    theme_minimal()
}
