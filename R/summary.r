#' @export
summary.mem_single <- function(object, ...) {

  oct <- c(object$post.prob, object$HPD[1], object$HPD[2], object$ESS, object$mean_est, 
               object$median_est)
  
  names(oct) <- c("Post.Prob", "HPD LB", "HPD HB", "ESS", "Mean", "Median")
  return(oct)
}

#' @export
summary.borrow_multiple <- function(object, ...) {
  
  ind <- object$drug_index
  dat <- object$data
  numInd <- length(ind)
  MAP <- PEP <- matrix(0, 0, numInd)
  post.prob <- ESS <- mean_est <- median_est<-c()
  HPD <- matrix(0, 2, 0)
  allName <- dat[[1]]$name
  indexName <- allName[ind]
  # vemu_wide1$responders/vemu_wide1$evaluable
  for(i in 1:numInd)
  {
    # Full Bayes
    s <- dat[[i]]
    # print(exact_single)
    #print(vemu_wide1$responders / vemu_wide$evaluable)
    # print(exact_single$MAP)
    # print(exact_single$PEP)
    MAP <- rbind(MAP, s$MAP)
    PEP <- rbind(PEP, s$PEP)
    post.prob <- c(post.prob, s$post.prob)
    ESS <- c(ESS, s$ESS)
    HPD <- cbind(HPD, s$HPD)
    mean_est <- c(mean_est, s$mean_est)
    median_est <- c(median_est, s$median_est)
  }
  
  
  #browser()
  names(ESS) <- indexName
  names(post.prob) <- names(mean_est) <- names(median_est) <- indexName
  rownames(PEP) <- rownames(MAP) <- indexName
  oct <- list(MAP = MAP, PEP = PEP, post.prob = post.prob,HPD = HPD, mean_est = mean_est, median_est = median_est, ESS = ESS)
  
  #names(oct) <- c("Post.Prob", "HPD LB", "HPD HB", "ESS", "Mean", "Median")
  return(oct)
}

