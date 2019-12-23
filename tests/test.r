#library(testthat)
library(basket)
library(tidyr)
library(tibble)
library(ggplot2)

library(borrow)
source("update.r")


data(vemu_wide)

baskets <- 1:6

vemu_wide1 <- vemu_wide[baskets, ]


allM <- allP <- matrix(0, 0, 6)
allPost <- ESS <-c()
aHPD <- matrix(0, 2, 0)
# vemu_wide1$responders/vemu_wide1$evaluable
for(i in 1:6)
{
  # Full Bayes
  exact_single <- mem_single(
    responses = vemu_wide1$responders,
    size = vemu_wide1$evaluable,
    name = vemu_wide1$baskets,
    drug_index = i, 
    p0 = 0.25
  )
  # print(exact_single)
  #print(vemu_wide1$responders / vemu_wide$evaluable)
  # print(exact_single$MAP)
  # print(exact_single$PEP)
  allM <- rbind(allM, exact_single$MAP)
  allP <- rbind(allP, exact_single$PEP)
  allPost <- c(allPost, exact_single$post.prob)
  ESS <- c(ESS, exact_single$ESS)
  
  aHPD <- cbind(aHPD, exact_single$HPD)
}




x <- exact_single
plot_borrow_density(x)
summary(x)




test2 <- mem_single(
  responses = c(3,4,10,9,2,11),
  size = rep(25, 6),
  name = c("1", "2", "3", "4", "5", "6"),
  drug_index = 2, 
  p0 = 0.22
)
print(test2$MAP)
summary(test2)

test3 <- mem_single(
  responses = c(3,4,10,9,2,11),
  size = rep(25, 6),
  name = c("1", "2", "3", "4", "5", "6"),
  drug_index = 3, 
  p0 = 0.22
)
print(test3$MAP)
summary(test3)



test4 <- borrow_multiple(
  responses = vemu_wide1$responders,
  size = vemu_wide1$evaluable,
  name = vemu_wide1$baskets,
  drug_index = 2:3, 
  p0 = 0.25
)
summary(test4)






testMult1 <- borrow_multiple(
  responses = vemu_wide1$responders,
  size = vemu_wide1$evaluable,
  name = vemu_wide1$baskets,
  drug_index = 1:6, 
  p0 = 0.25
)

res <- summary(testMult1)$PEP
summary(testMult1)$post.prob
dat <- data.frame(res)
library(knitr)
kable(dat, caption="Summary of Posterior Probabilities", align='lcc', row.names=F, format='html')

ess <- summary(testMult1)$ESS
df <- data.frame(baskets = vemu_wide1$baskets, ESS = ess)

p<-ggplot(df, aes(x=baskets, y=ESS, fill = baskets)) + labs(y = "ESS")+
  geom_bar(stat="identity")+theme(axis.title.x = element_blank())+
  geom_hline(yintercept=0.05, linetype="dashed", color = "red")
p




library(doParallel)
library(foreach)

numCores <- detectCores()
numCores
registerDoParallel(numCores)

resp <- vemu_wide1$responders
resp[1:2] <- rep(0.1, 2)
resp[4] <- 0.3
resp
test5 <- borrow_simulate(
  resp =  resp,
  is.resp.rate = c(TRUE, TRUE, FALSE, TRUE, FALSE, FALSE),
  size = vemu_wide1$evaluable,
  name = vemu_wide1$baskets,
  drug_index = c(2, 4), 
  p0 = 0.2,
  num_sim = 100
)
#saveRDS(test5, "simres.rds")

r <- summary(test5)
print(r)
# thr <- calibrate(test5, prob= (1:9)/ 10)
# thr
# plot_sim(test5, threshold = thr[,2])
plot_sim_violin(test5)


# OC curve
oc <- ocCurve(nullData = r$allPostProb[,1], alterData= r$allPostProb[,2])
oc
plot.occurve(oc) 
cali.onPower(oc, powerV = c(0.7, 0.8, 0.9))
cali.onTypeIError(oc, typeIError = c(0.1, 0.2, 0.3))


resp1 <- c(0.1, 0.1, 1.0, 0.3, 6.0, 2.0)
resp2 <- c(0.3, 0.1, 1.0, 0.3, 6.0, 2.0)
resp3 <- c(0.1, 0.2, 1.0, 0.3, 6.0, 2.0)
resp4 <- c(0.3, 0.3, 1.0, 0.3, 6.0, 2.0)

resps <- rbind(resp1, resp2, resp3, resp4)
test6 <- borrow_simulate_multiple(
  resp.scenarios = resps,
  is.resp.rate = c(TRUE, TRUE, FALSE, TRUE, FALSE, FALSE),
  size = vemu_wide1$evaluable,
  name = vemu_wide1$baskets,
  drug_index = c(2, 4), 
  p0 = 0.2,
  num_sim = 30,
  output.file = "testresult.RDS"
)


library(doParallel)
library(foreach)

numCores <- detectCores()
numCores
registerDoParallel(numCores)

test7 <- borrow_simulate(
  resp =  resp,
  is.resp.rate = c(TRUE, TRUE, FALSE, TRUE, FALSE, FALSE),
  size = vemu_wide1$evaluable,
  name = vemu_wide1$baskets,
  drug_index = c(2, 4), 
  interim_size = c(4,5),
  p0 = 0.2,
  num_sim = 24
)


r <- summary(test7)
print(r)
plot_sim_violin(test7)
plot_sim_interim_violin(test7, interim = 2)