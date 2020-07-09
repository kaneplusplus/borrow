  

#' @title borrow simulation
#' @description borrow trial simulations for multiple subgroups using borrow method.
#' @param resp the response outcome or the true response rates in each basket.
#' @param is.resp.rate the vector to indicate if the true_rate is response rate (FALSE means the resp value is response outcome)
#' @param size the size of each basket.
#' @param name the name of each basket.
#' @param drug_index the index vector of the basket to be studied.
#' @param interim_size basket sizes for interim analyses 
#' @param p0 the null response rate vector for the poster probability calculation
#' (default 0.15).
#' @param num_sim the number of simulationst
#' (default 100).
#' @param shape1 the first shape parameter(s) for the prior of each basket
#' (default 0.5).
#' @param shape2 the second shape parameter(s) for the prior of each basket
#' (default 0.5).
#' @param prior the matrix giving the prior inclusion probability
#' for each pair of baskets. The default is on on the main diagonal and 0.5
#' elsewhere.
#' @param hpd_alpha the highest posterior density trial significance.
#' @param alternative the alternative case definition (default greater)
#' @param call the call of the function (default NULL).
#' @importFrom stats rbinom
#' @examples
#' # 6 baskets, each with enrollement size 5
#' trial_sizes <- rep(25, 6)
#' 
#' # The response rates for the baskets.
#' resp_rate <- 0.15
#' 
#' # The trials: a column of the number of responses and a column of the
#' # the size of each trial.
#' trials <- data.frame(
#'   resprate = rep(0.15, 6),
#'   size = trial_sizes,
#'   name = letters[1:6]
#' )
#' 
#' borrow_simulate(trials$resprate, trials$size, trials$name, drug_index = 1:2)
#' 
#' @importFrom foreach foreach %dopar% getDoParName getDoSeqName registerDoSEQ
#' %do%
#' @importFrom stats median var
#' @importFrom R.utils insert
#' @importFrom crayon red
#' @importFrom foreach foreach
#' @importFrom itertools isplitRows
#' @export
borrow_simulate <- function(resp,
                      is.resp.rate,      
                      size,
                      name,
                      drug_index,
                      interim_size = NA,
                      p0 = 0.15,
                      num_sim =100,
                      shape1 = 0.5,
                      shape2 = 0.5,
                      prior = rep(1.0, length(size)) / 2,
                      hpd_alpha = 0.05,
                      alternative = "greater",
                      call = NULL) {
  
  if (length(resp) != length(size)) {
    stop(red(
      "The length of the response and size parameters",
      "must be equal."
    ))
  }

  numGroup <- length(size)
  allResp <- matrix(0, num_sim, numGroup)

#browser()
 
  allResu <- foreach(i = 1:num_sim, .combine = 'c') %dopar% {
    intRes <- list()
    if (!is.na(interim_size)){
      
      respC <- rep(0, numGroup)
      sizeC <- rep(0, numGroup)
      for(ii in 1:length(interim_size)){
        resS <- c()
        total <- c()
        for(j in 1:numGroup){
          if(is.resp.rate[j]){
            total <- c(total, interim_size[ii])
            resS <- c(resS, respC[j] + rbinom(1, interim_size[ii] - sizeC[j], resp[j]))
          } else {
            total <- c(total, size[j])
            resS <- c(resS, resp[j])           
          }
        }
        r <- borrow_multiple(responses =  resS,
                             size = total,
                             name = name,
                             drug_index = drug_index,
                             p0 = p0,
                             shape1 = shape1,
                             shape2 = shape2,
                             prior = prior,
                             hpd_alpha = hpd_alpha,
                             alternative = alternative)
        
        t <- summary(r)
        intRes[[ii]] <- list( res = t, size = total, resp = resS)
        respC <- resS
        sizeC <- total
      }
      allR <- c()
      for(j in 1:numGroup){
        if(is.resp.rate[j]){
          allR<- c(allR, respC[j] + rbinom(1, size[j] - sizeC[j], resp[j]))
        } else {
          allR <- c(allR, resp[j])
        }
      }
      
    }else{
      allR <- c()
      for(j in 1:numGroup){
        if(is.resp.rate[j]){
          allR<- c(allR, rbinom(1, size[j], resp[j]))
        } else {
          allR <- c(allR, resp[j])
        }
      }
    }
    
    r <- borrow_multiple(responses =  allR,
                         size = size,
                         name = name,
                         drug_index = drug_index,
                         p0 = p0,
                         shape1 = shape1,
                         shape2 = shape2,
                         prior = prior,
                         hpd_alpha = hpd_alpha,
                         alternative = alternative
    )
    t <- summary(r)
    res <- list(t = t, interim_res = intRes, finalResp = allR)
    list(res)
  }

  ret <- list(data = allResu, name = name,  drug_index = drug_index, interim_size = interim_size, resp = resp, is.resp.rate = is.resp.rate, size = size)
  class(ret) <- c("borrow_simulate", "exchangeability_model")
  return(ret)
}
