#' @export
summary.borrow_single <- function(object, ...) {

  oct <- c(round(object$post.prob, 3), object$HPD[1], object$HPD[2], round(object$ESS,3),
           round(object$mean_est, 3), 
           round(object$median_est,3))
  oct <- round(oct, 3)
  
  names(oct) <- c("Post.Prob", "HPD LB", "HPD HB", "ESS", "Mean", "Median")
  return(oct)
}

#' @export
summary.borrow_multiple <- function(object, ...) {
  
  ind <- object$drug_index
  dat <- object$data
  numInd <- length(ind)

  post.prob <- ESS <- mean_est <- median_est<-c()
  HPD <- matrix(0, 2, 0)
  allName <- dat[[1]]$name
  indexName <- allName[ind]
  MAP <- PEP <- matrix(0, 0, length(allName))
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
  ESS <- round(ESS, 3)
  post.prob <- round(post.prob, 3)
  PEP <- round(PEP, 3)
  names(ESS) <- indexName
  names(post.prob) <- names(mean_est) <- names(median_est) <- indexName
  rownames(PEP) <- rownames(MAP) <- colnames(HPD) <- indexName
  oct <- list(MAP = MAP, PEP = PEP, post.prob = post.prob,HPD = round(HPD, 3), 
              mean_est = round(mean_est, 3), median_est = round(median_est, 3), ESS = ESS)
  
  #names(oct) <- c("Post.Prob", "HPD LB", "HPD HB", "ESS", "Mean", "Median")
  return(oct)
}

