MAP <- PEP <- matrix(0, 0, length(allName))

library(matrixStats)
calibrate <- function(simResult, prob = c(0.1, 0.8))
{
  data <- simResult$data
  numSim <- length(data)
  #r <- data[[1]]
  index <- simResult$drug_index
  numInd <- length(index)
  
  allPostProb <- matrix(0, 0, numInd)
  for (i in 1:numSim){
    allPostProb <- rbind(allPostProb, data[[i]]$post.prob)
  }
  thr <- colQuantiles(allPostProb, probs = prob)
  
  thr
  
}


library(ggplot2)
library(dplyr)
plot_sim <- function(simResult, threshold)
{
  data <- simResult$data
  numSim <- length(data)
  #r <- data[[1]]
  index <- simResult$drug_index
  numInd <- length(index)
  
  allPostProb <- matrix(0, 0, numInd)
  for (i in 1:numSim){
    allPostProb <- rbind(allPostProb, data[[i]]$post.prob)
  }
  
  allP <- allPostProb
  nDim <- dim(allP)
  dfAll <- data.frame(name=character(),
                      Post.prob=double())
  allName <- colnames(allP)
  for(i in 1:nDim[2])
  {
    nameS <- allName[i]
    df <- data.frame(name=rep(nameS, nDim[1]), Post.prob=allP[,i])
    dfAll <-rbind(dfAll, df)
  }
  
  #threshold <- c(0.1, 0.7)
  p <- ggplot(dfAll, aes(x = Post.prob)) +
    geom_density(fill = "lightblue") +  geom_density(alpha=0.4) +
    facet_wrap(. ~ name, ncol = 2) + 
    theme_minimal()
  for(i in 1:2)
  {
    p <- p + geom_vline(data = filter(dfAll, name == allName[i]), aes_string(xintercept=threshold[i]),
                        color="blue", size = 1)
  }
  p
}
