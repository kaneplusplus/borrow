#' @export
summary.mem_single <- function(object, ...) {

  oct <- c(object$post.prob, object$HPD[1], object$HPD[2], object$ESS, object$mean_est, 
               object$median_est)
  
  names(oct) <- c("Post.Prob", "HPD LB", "HPD HB", "ESS", "Mean", "Median")
  return(oct)
}

