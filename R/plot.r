#' @title Plot the posterior density of the fitting result
#'
#' @description Create the posterior density plot of the studied basket.
#' @param x the fitted model calculated from the mem_single() function.

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

