#library(testthat)
library(basket)
library(tidyr)
library(tibble)
library(ggplot2)
library(borrow)



data(vemu_wide)


baskets <- 1:6

vemu_wide1 <- vemu_wide[baskets, ]


s <- mem_exact(
  responses = vemu_wide1$responders,
  size = vemu_wide1$evaluable,
  name = vemu_wide1$baskets, cluster_analysis = TRUE,
  p0 = 0.25
)
plot_density(s)

allM <- allP <- matrix(0, 0, 6)
allPost <- ESS <-c()
aHPD <- matrix(0, 2, 0)
# vemu_wide1$responders/vemu_wide1$evaluable
for(i in 1:6)
{
  # Full Bayes
  exact_single <- borrow_single(
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




test2 <- borrow_single(
  responses = c(3,4,10,9,2,11),
  size = rep(25, 6),
  name = c("1", "2", "3", "4", "5", "6"),
  drug_index = 2, 
  p0 = 0.22
)
print(test2$MAP)
summary(test2)

test3 <- borrow_single(
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
x <- test4
summary(test4)
plot_borrow_density(test4)


