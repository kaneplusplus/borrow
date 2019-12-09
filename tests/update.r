
summary.borrow_simulate <- function(object, ...) {
  simResult <- object
  data <- simResult$data
  numSim <- length(data)
  #r <- data[[1]]
  index <- simResult$drug_index
  numInd <- length(index)
  
  allPostProb <- allESS  <- matrix(0, 0, numInd)
  for (i in 1:numSim){
    allPostProb <- rbind(allPostProb, data[[i]]$post.prob)
    allESS <- rbind(allESS, data[[i]]$ESS)
  }
  allMap <- data[[1]]$MAP
  for (i in 2:numSim){
    allMap <- allMap + data[[i]]$MAP
  }  
  allMap <- allMap / numSim
  post.bound <- colQuantiles(allPostProb, probs = c(0.25, 0.75))
  postmean <- colMeans(allPostProb)
  ess.bound <- colQuantiles(allESS, probs = c(0.25, 0.75))
  essmean <- colMeans(allESS)
  res <- cbind(postmean, post.bound, essmean, ess.bound)
  colnames(res) <- c("Post.prob Mean", "Post.prob 25%", "Post.prob 75%", "ESS Mean",  "ESS 25%", "ESS75%")
  list(num_sim = numSim, name = simResult$name,  drug_index = simResult$drug_index, 
       resp = simResult$resp, is.resp.rate = simResult$is.resp.rate, size = simResult$size, 
       allPostProb = allPostProb, allESS = allESS, Avg.MAP = allMap,  result = res)
}

# OC curve
ocCurve <- function(nullData, alterData)
{
  allData <- c(nullData, alterData)
  allData <- unique(allData)
  cutoff <- sort(allData)
  typeIError <- c()
  powerVal <- c()
  numNull <- length(nullData)
  numAlter <- length(alterData)
  for(i in 1:length(cutoff)){
    t1 <- sum(nullData >= cutoff[i]) / numNull
    po <- sum(alterData >= cutoff[i]) / numAlter
    typeIError <- c(typeIError, t1)
    powerVal <- c(powerVal, po) 
  }
  res <- data.frame(cutoff, typeIError, powerVal)
}

plot.occurve <- function(res)
{
  #plot(res$typeIError, res$powerVal, xlab = "Type I error Rate", ylab = "Power", type ="o")
  g <- ggplot(res, aes(x=typeIError, y=powerVal)) +
    geom_point(size=2, shape=23) +
    geom_path(size = 1)+
    xlab("Type I Error ") +
    ylab("Power")
  g
}

cali.onPower<- function(res, powerV = c(0.7, 0.8, 0.9))
{
  p <- powerV
  sm <- smooth.spline(res$powerVal, res$typeIError,  spar = 0.3)
  predTError <- predict(sm, x = p)$y
  x <- (1:1000)/1000
  predAll <- predict(sm, x = x)$y 
  smCurve <- data.frame(x, predAll)
  
  smCutoff <- smooth.spline(res$powerVal, res$cutoff,  spar = 0.3)
  
  Cutoff <- predict(smCutoff, x = p)$y
  # n <- 10
  # d <- data.frame(x = 1:n, y = rnorm(n))
  # ggplot(d,aes(x,y)) + geom_point() + 
  #   geom_line(data=data.frame(spline(d, n=n*10)))
  
  p <- plot.occurve(res) + geom_line(data=smCurve, aes(predAll, x), color="blue", size =1)
  print(p)
  data.frame(powerV, Cutoff, predTError)
}

cali.onTypeIError<- function(res, typeIError = c(0.1, 0.2, 0.3))
{
  p <- typeIError
  sm <- smooth.spline(res$typeIError, res$powerVal, spar = 0.3)
  predPower <- predict(sm, x = p)$y
  x <- (1:1000)/1000
  predAll <- predict(sm, x = x)$y 
  smCurve <- data.frame(x, predAll)
  
  smCutoff <- smooth.spline(res$typeIError, res$cutoff,  spar = 0.3)
  
  Cutoff <- predict(smCutoff, x = p)$y
  # n <- 10
  # d <- data.frame(x = 1:n, y = rnorm(n))
  # ggplot(d,aes(x,y)) + geom_point() + <- 
  #   geom_line(data=data.frame(spline(d, n=n*10)))
  
  p <- plot.occurve(res) + geom_line(data=smCurve, aes(x,predAll), color="blue", size =1)
  print(p)
  data.frame(typeIError, predCutoff, predPower)
}


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

# 
# library(pROC)
# data(aSAH)
# if(!require(DT)) install.packages("DT")
# DT::datatable(aSAH)
# 
# 
# rocobj <- plot.roc(aSAH$outcome, aSAH$s100b, percent = TRUE, main = "Smoothing")
# lines(smooth(rocobj), # smoothing (default: binormal)
#       col = "#1c61b6")
# lines(smooth(rocobj, 
#              method = "density"), # density smoothing
#       col = "#008600")
# lines(smooth(rocobj, 
#              method = "fitdistr", # fit a distribution
#              density = "lognormal"), # let the distribution be log-normal
#       col = "#840000")
