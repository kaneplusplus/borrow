
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
ocCurve <- function(nullData, alterData, digits = 3)
{
  allData <- c(nullData, alterData)
  if(!is.na(digits)){
    allData <-round(allData, digits)
  }
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
  data.frame(powerV, round(Cutoff, 3), round(predTError, 3))
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
  data.frame(typeIError, round(predCutoff , 3), round(predPower, 3))
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

plot_sim_violin <- function(simResult)
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
  
  #  = c(PP1, PP2, PP3)
  # Y=factor( c(rep("basket1",length(PP1)),rep("basket2",length(PP2)),rep("basket3",length(PP3))) )
  # #plot.XY(X, Y, x_name="Posterior Prob", y_name="Basket", pl.main=NULL, geom_violin.trim = TRUE, coord_flip=TRUE)
  # 
  # 
   
  plot.XY(dfAll$Post.prob, dfAll$name, y_name="Baskets", x_name="Posterior Prob", pl.main=NULL, geom_violin.trim = TRUE, coord_flip=FALSE)
}


# 
# library(pROC)
# data(aSAH)
# if(!require(DT)) install.packages("DT")


library(stringr)
library(ggplot2)
library(dplyr)

nLevels <- function(x){ return( sapply(X=levels(x), FUN=function(X,s) length(which(s==X)), s=x) ) }

append.n <- function(y){ n <- nLevels(y); y. <- as.character(y)
for(i in 1:length(levels(y))){ 
  u <- which(y==levels(y)[i]); y.[u] <- paste0(y[u]," (",n[i],")") }
return(factor(y.))
}

signDigit <- function(x){
  if(is.na(x)){ out <- NA 
  } else{ if(abs(x)>100){ out <- round(x)
  } else if(abs(x)>10){ out <- round(x,1) 
  } else if(abs(x)>1){ out <- round(x,2) 
  } else if(abs(x)>0.1){ out <- round(x,3) 
  } else{ s <- 0.01; k <- 1; while(abs(x)<s){ s <- s/10; k <- k+1 }; out <- round(x,k+2) }
  }
  return(out)
}

Kruskal <- function(U,cl.ORD=NA){
  g <- U$g; y <- U$y
  L <- levels(g)
  n <- nLevels(g); prop <- round(n/sum(!is.na(g)),3)
  N.perc <- Med.range <- rep(NA,length(L)) 
  for(i in 1:length(L)){ 
    N.perc[i] <- paste(n[i]," (",prop[i]*100,")",sep="") 
    Med.range[i] <- paste(signDigit(median(y[g==L[i]],na.rm=TRUE))," (",signDigit(min(y[g==L[i]],na.rm=TRUE)),"-",signDigit(max(y[g==L[i]],na.rm=TRUE)),")",sep="")
  }
  fit <- kruskal.test(x=y,g=g)
  if(fit$p.value>0.05){p. <- round(fit$p.value,3) } else{ p. <- signDigit(fit$p.value) }
  R <- cbind(N.perc,Med.range,rep(p.,length(L))); colnames(R) <- c("N (perc) g","Median-Range y","p-value"); rownames(R) <- L
  if(!is.na(cl.ORD)[1]){ R <- R[cl.ORD,]; rownames(R) <- paste0("cluster-",cl.ORD) }
  #write.csv(R,"R.csv")
  #print(paste("N = ",dim(na.omit(cbind(g,y)))[1],sep=""))
  return(R)
}

plot.XY <- function(X, Y, x_name="", y_name="", pl.main=NULL, Part.numeric.X=FALSE, 
                    X.colorPalette=colorRamps::matlab.like,
                    Y.colorPalette=colorRamps::matlab.like,
                    maxSplit=4, 
                    risk.table=function(x){ return(x<=10) },
                    risk.table.fontsize=function(x){
                      a <- solve(matrix(nrow=2,ncol=2,data = c(1,1,2,11)),c(5,3.3))
                      out <- max(min(a[1]+a[2]*x,5),1)
                      return(out) },
                    risk.table.height=0.4,
                    legend = "top",
                    surv.median.line="hv",
                    Y.order.level=NULL, 
                    CF=function(x,y){ 
                      x. <- c(30,2,30)
                      y. <- c(6,6,10)
                      z. <- c(0.22,1,0.34)
                      a <- solve(matrix(nrow=3,ncol=3,data = c(rep(1,3),x.,y.)),z.)
                      out <- max(min(a[1]+a[2]*x+a[3]*y,1),0.01)
                      return(out) }, 
                    truncLabs.g=7, 
                    truncLabs.y=7,
                    geom_tile.color="white",
                    scale_fill_gradient2.low="lightblue", 
                    scale_fill_gradient2.high="red",
                    scale_fill_gradient2.mid="white", 
                    midpoint=mean,
                    geom_text.color=1,
                    geom_text.size=function(x){
                      a <- solve(matrix(nrow=2,ncol=2,data = c(1,1,2,30)),c(4,3))
                      out <- max(min(a[1]+a[2]*x,4),1)
                      return(out) },
                    plot.title.color=1,
                    plot.title.face="bold",
                    plot.title.size=16,
                    plot.title.y.color=1,
                    plot.title.y.face="bold",
                    axis.title.y.size=function(x){
                      a <- solve(matrix(nrow=2,ncol=2,data = c(1,1,2,30)),c(12,10))
                      out <- max(min(a[1]+a[2]*x,10),4)
                      return(out) },
                    plot.title.x.color=1,
                    plot.title.x.face="bold",
                    axis.title.x.size=function(x){
                      a <- solve(matrix(nrow=2,ncol=2,data = c(1,1,2,30)),c(12,8))
                      out <- max(min(a[1]+a[2]*x,10),4)
                      return(out) },
                    axis.text.x.color=1,
                    axis.text.x.face="bold",
                    axis.text.x.size=function(x){
                      a <- solve(matrix(nrow=2,ncol=2,data = c(1,1,2,30)),c(10,6))
                      out <- max(min(a[1]+a[2]*x,10),4)
                      return(out) },
                    axis.text.y.color=1,
                    axis.text.y.face="bold",
                    axis.text.y.size=function(x){
                      a <- solve(matrix(nrow=2,ncol=2,data = c(1,1,2,10)),c(9,6))
                      out <- max(min(a[1]+a[2]*x,10),4)
                      return(out) },
                    geom_jitter.height=0,
                    geom_jitter.width=0.2,
                    legend.position="right",
                    legend.text.color=1,
                    legend.text.face="bold",
                    legend.text.size=10,
                    panel.background_blank=FALSE,
                    geom_violin.trim = FALSE, 
                    geom_violin.draw_quantiles = c(0.25, 0.5, 0.75),
                    coord_flip=FALSE ){
  out <- list()
  
  if(is.null(pl.main)){ 
    if(x_name!="" && y_name!=""){ pl.main <- paste0(y_name," by ",x_name)
    } else{ pl.main=c("Y by X") }
  }
  
  panel.back=element_rect(); if(panel.background_blank){ panel.back = element_blank() }
  
  Y.vec <- is.vector(Y) || length(ncol(Y))<1
  if( !Y.vec ){
    
    if( ncol(Y)>2 ){ stop("Improper Format Y: More than 2 cols in Y") 
    } else if( ncol(Y)==2 ){
      
      if(is.factor(X)){
        
        X. <- NULL; try( X. <- trunc.Factors(X=X, Start=truncLabs.g) )
        if( length(X.)>0 ){ X <- X. }
        TT <- as.data.frame( list(TTF=Surv(Y[,1], Y[,2]), x=X) )
        sdf <- survdiff( TTF ~ x, rho=0, data=TT)
        sft <- survfit( TTF ~ x, data=TT )
        
        thisPal <- X.colorPalette(length(levels(TT$x))) 
        if( "#FFFFFF" %in% thisPal ){
          temp <- X.colorPalette(length(levels(TT$x))+1); uu <- which( temp=="#FFFFFF" )
          if(length(uu)>0){ temp <- temp[-uu] } else{ temp <- temp[-length(temp)] }
          thisPal <- temp 
        }
        out. <- NULL
        
        try( out. <- survminer::ggsurvplot(sft, data=TT, palette=thisPal, ylab="proportion without event", title=pl.main, #subtitle=paste0("p-value = ",signDigit(1 - pchisq(sdf$chisq, length(sdf$n) - 1))),  
                                           risk.table=risk.table(length(levels(TT$x))), fontsize=risk.table.fontsize(length(levels(TT$x))), legend = legend, pval=TRUE, pval.size=4.5,
                                           risk.table.height=risk.table.height, surv.median.line=surv.median.line, legend.title="", legend.labs=levels(TT$x)) )
        
        try( out. <- out. + survminer::theme_survminer(
          font.main = c(plot.title.size, plot.title.face, plot.title.color),
          font.submain = c(plot.title.size-4, "plain", plot.title.color),
          font.x = c(axis.text.x.size(1), axis.text.x.face, axis.text.x.color),
          font.y = c(12, axis.text.y.face, axis.text.y.color) ) )
        
        if(length(out.)>0){ out[[length(out)+1]] <- out.; names(out)[length(out)] <- paste0(y_name,"..",x_name) }
        
      } else if(is.numeric(X)){
        if(Part.numeric.X){
          out. <- NULL
          try( out. <- cut.Cind(maxSplit=maxSplit,Dur=Y[,1],Event=Y[,2],X=X,tlab=y_name,xlab=x_name,ids=NULL,plot=TRUE,
                                pars=list(X.colorPalette=X.colorPalette, risk.table=risk.table, risk.table.fontsize=risk.table.fontsize, risk.table.height=risk.table.height, legend=legend, surv.median.line=surv.median.line),pl.main=pl.main ) )
          if(length(out.)>0){ for(u in 1:length(out.$pl)){ out[[length(out)+1]] <- out.$pl[[u]]; names(out)[length(out)] <- names(out.$pl)[u] }}
        }
      } else{ stop("X should be factor or numeric") }
      
    } else{ Y.vec <- TRUE }
  }
  
  if(Y.vec){
    
    if( is.numeric(X) && is.numeric(Y) ){
      out. <- NULL
      
      main=paste0(y_name," by ",x_name)
      df <- as.data.frame( list(x=X, y=Y) )
      PVAL <- 1.00; COR <- 0; cor3 <- NULL
      try( cor3 <- suppressWarnings( cor.test( df$y, df$x, method="spearman", alternative="two.sided" ) ) )
      if(length(cor3)>0){ COR <- cor3$estimate; PVAL <- signDigit(cor3$p.value) }
      
      out. <- ggplot2::ggplot(data=df, aes(x=X, y=Y)) + geom_point() + geom_smooth() +
        theme(plot.title = element_text(color=plot.title.color, face=plot.title.face, size=plot.title.size),
              axis.title.x = element_text(color=plot.title.y.color,face=plot.title.y.face,size=axis.title.y.size(1)),       
              axis.title.y = element_text(color=plot.title.x.color,face=plot.title.x.face,size=axis.title.x.size(1)),
              panel.grid.major = element_blank(),
              panel.background=panel.back,
              axis.text.x=element_text(color=axis.text.y.color,face=axis.text.y.face,size=axis.text.y.size(1)),
              axis.text.y=element_text(color=axis.text.x.color,face=axis.text.x.face,size=axis.text.x.size(1))) + 
        labs(y=paste0(y_name,"\n"),x=paste0("\n",x_name),title=main,subtitle=paste0("Spearman rank ",round(COR,3), "; p-value = ",PVAL) )
      
      if(length(out.)>0){ out[[length(out)+1]] <- out.; names(out)[length(out)] <- paste0(y_name,"..",x_name) }
      
    } else if( is.numeric(X) || is.numeric(Y) ){
      
      main=paste0(y_name," by ",x_name)
      panel.back=element_rect(); if(panel.background_blank){ panel.back = element_blank() }
      
      if( is.factor(X) ){ g <- X; y <- Y; g_lab=x_name; y_lab=y_name
      gPal <- X.colorPalette(length(levels(g)))
      if( "#FFFFFF" %in% gPal ){
        temp <- X.colorPalette(length(levels(g))+1); uu <- which( temp=="#FFFFFF" )
        if(length(uu)>0){ temp <- temp[-uu] } else{ temp <- temp[-length(temp)] }
        gPal <- temp 
      }
      } else{ g <- Y; y <- X; y_lab=x_name; g_lab=y_name
      gPal <- Y.colorPalette(length(levels(g)))
      if( "#FFFFFF" %in% gPal ){
        temp <- Y.colorPalette(length(levels(g))+1); uu <- which( temp=="#FFFFFF" )
        if(length(uu)>0){ temp <- temp[-uu] } else{ temp <- temp[-length(temp)] }
        gPal <- temp 
      }
      }
      PVAL <- 1.00
      if( sum(nLevels(g)>1)>1 ){
        K.test <- Kruskal( as.data.frame( list(g=g, y=y) ) )
        PVAL <- signDigit(as.numeric(K.test[1,ncol(K.test)]))
      }
      
      g. <- NULL; try( g. <- trunc.Factors(X=g, Start=truncLabs.g) )
      if( length(g.)>0 ){ g <- g. }
      #df <- as.data.frame( list(g=append.n(g), y=y) )
      df <- as.data.frame( list(g=g, y=y) )
      df2<-df%>%group_by(g) %>% summarise(n=n())
      out. <- NULL
      
      if(coord_flip){ 
        out. <- ggplot2::ggplot(aes(g,y), data=df) + 
          geom_violin(aes(fill=g), trim = geom_violin.trim, draw_quantiles = geom_violin.draw_quantiles) +
          geom_jitter(height = geom_jitter.height, width = geom_jitter.width) +
          scale_fill_manual(values = gPal) + scale_x_discrete(labels = paste0(names(nLevels(g)),"\n(",nLevels(g),")")) +
          theme(plot.title = element_text(color=plot.title.color, face=plot.title.face, size=plot.title.size),
                axis.title.x = element_text(color=plot.title.x.color,face=plot.title.x.face,size=axis.title.x.size(length(levels(g)))),       
                axis.title.y = element_text(color=plot.title.y.color,face=plot.title.y.face,size=axis.title.y.size(1)),
                panel.grid.major = element_blank(),
                panel.background=panel.back,
                legend.position=legend.position,
                legend.title=element_text(color=legend.text.color,face=legend.text.face,size=legend.text.size+1),
                legend.text=element_text(color=legend.text.color,face=legend.text.face,size=legend.text.size),
                axis.text.x=element_text(color=axis.text.x.color,face=axis.text.x.face,size=axis.text.x.size(length(levels(g)))),
                axis.text.y=element_text(color=axis.text.y.color,face=axis.text.y.face,size=axis.text.y.size(1))) + 
          labs(x=paste0(g_lab,"\n"),y=paste0("\n",y_lab),title=main,subtitle=paste0("p-value = ",PVAL) ) +
          guides(fill = guide_legend(title = g_lab)) + coord_flip() 
      } else{ 
        out. <- ggplot2::ggplot(aes(g,y), data=df) + 
          geom_violin(aes(fill=g), trim = geom_violin.trim, draw_quantiles = geom_violin.draw_quantiles) +
          geom_jitter(height = geom_jitter.height, width = geom_jitter.width) +
          scale_fill_manual(values = gPal) + scale_x_discrete(labels = paste0(names(nLevels(g)),"\n(",nLevels(g),")")) +
          theme(plot.title = element_text(color=plot.title.color, face=plot.title.face, size=plot.title.size),
                axis.title.x = element_text(color=plot.title.y.color,face=plot.title.y.face,size=axis.title.y.size(1)),       
                axis.title.y = element_text(color=plot.title.x.color,face=plot.title.x.face,size=axis.title.x.size(length(levels(g)))),
                panel.grid.major = element_blank(),
                panel.background=panel.back,
                legend.position=legend.position,
                legend.title=element_text(color=legend.text.color,face=legend.text.face,size=legend.text.size+1),
                legend.text=element_text(color=legend.text.color,face=legend.text.face,size=legend.text.size),
                axis.text.x=element_text(color=axis.text.y.color,face=axis.text.y.face,size=axis.text.y.size(1)),
                axis.text.y=element_text(color=axis.text.x.color,face=axis.text.x.face,size=axis.text.x.size(length(levels(g))))) + 
          labs(x=paste0("\n",g_lab),y=paste0(y_lab,"\n"),title=main,subtitle=paste0("p-value = ",PVAL) ) +
          guides(fill = guide_legend(title = g_lab))
      }
      
      if(length(out.)>0){ out[[length(out)+1]] <- out.; names(out)[length(out)] <- paste0(y_name,"..",x_name) }
      
    } else{
      
      out. <- NULL
      try( out. <- createContingTable(g=X, y=Y, g_lab=x_name, y_lab=y_name, main=pl.main, 
                                      Y.order.level=Y.order.level, CF=CF, 
                                      truncLabs.g=truncLabs.g, truncLabs.y=truncLabs.y,
                                      geom_tile.color=geom_tile.color,
                                      scale_fill_gradient2.low=scale_fill_gradient2.low, 
                                      scale_fill_gradient2.high=scale_fill_gradient2.high,
                                      scale_fill_gradient2.mid=scale_fill_gradient2.mid, 
                                      midpoint=midpoint,
                                      geom_text.color=geom_text.color,
                                      geom_text.size=geom_text.size,
                                      plot.title.color=plot.title.color,
                                      plot.title.face=plot.title.face,
                                      plot.title.size=plot.title.size,
                                      plot.title.y.color=plot.title.y.color,
                                      plot.title.y.face=plot.title.y.face,
                                      axis.title.y.size=axis.title.y.size,
                                      plot.title.x.color=plot.title.x.color,
                                      plot.title.x.face=plot.title.x.face,
                                      axis.title.x.size=axis.title.x.size,
                                      axis.text.x.color=axis.text.x.color,
                                      axis.text.x.face=axis.text.x.face,
                                      axis.text.x.size=axis.text.x.size,
                                      axis.text.y.color=axis.text.y.color,
                                      axis.text.y.face=axis.text.y.face,
                                      axis.text.y.size=axis.text.y.size ) )
      
      if(length(out.)>0){ out[[length(out)+1]] <- out.; names(out)[length(out)] <- paste0(y_name,"..",x_name) }
    }
  }
  
  return(out)
}




#### Nan: Try this function, 
## stack up the posterior probabilities into single vector,
# X = c(PP1, PP2, PP3)

# Mark these with factor level variable Y
# Y=factor( c(rep("basket1",length(PP1)),rep("basket1",length(PP12)),rep("basket1",length(PP3))) )

# PP1 <- rnorm(1000, 0, 1)
# 
# PP2 <- rnorm(1000, 8, 1)
# PP3 <- rnorm(1000, 38, 1)
# X = c(PP1, PP2, PP3)
# Y=factor( c(rep("basket1",length(PP1)),rep("basket2",length(PP2)),rep("basket3",length(PP3))) )
# #plot.XY(X, Y, x_name="Posterior Prob", y_name="Basket", pl.main=NULL, geom_violin.trim = TRUE, coord_flip=TRUE)
# 
# 
# 
# plot.XY(X, Y, y_name="Posterior Prob", x_name="Basket", pl.main=NULL, geom_violin.trim = TRUE, coord_flip=FALSE)
# 

