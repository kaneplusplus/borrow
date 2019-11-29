  

#' @title Fit the MEM Model for multiple subgroups using borrow
#'
#' @description Fit the MEM model for multiple subgroups using borrow method.
#' @param responses the number of responses in each basket.
#' @param size the size of each basket.
#' @param name the name of each basket.
#' @param drug_index the index vector of the basket to be studied.
#' @param p0 the null response rate vector for the poster probability calculation
#' (default 0.15).
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
#'   responses = rbinom(trial_sizes, trial_sizes, resp_rate),
#'   size = trial_sizes,
#'   name = letters[1:6]
#' )
#' 
#' borrow_multiple(trials$responses, trials$size, trials$name, drug_index = 1:2)
#' 
#' @importFrom foreach foreach %dopar% getDoParName getDoSeqName registerDoSEQ
#' %do%
#' @importFrom stats median var
#' @importFrom R.utils insert
#' @importFrom crayon red
#' @importFrom itertools isplitRows
#' @export
borrow_multiple <- function(responses,
                      size,
                      name,
                      drug_index,
                      p0 = 0.15,
                      shape1 = 0.5,
                      shape2 = 0.5,
                      prior = rep(1.0, length(responses)) / 2,
                      hpd_alpha = 0.05,
                      alternative = "greater",
                      call = NULL) {

  
  if (length(responses) != length(size)) {
    stop(red(
      "The length of the responses and size parameters",
      "must be equal."
    ))
  }
  
  if (length(shape1) == 1) {
    shape1 <- rep(shape1, length(responses))
  }
  if (length(shape2) == 1) {
    shape2 <- rep(shape2, length(responses))
  }
  
  if (length(p0) == 1) {
    p0 <- rep(p0, length(responses))
  }
  

  allResu <- list()
  for (i in 1:length(drug_index))
  {
    ind <- drug_index[i]
    p0_v = p0[ind]
    r <- mem_single(responses,
                    size,
                    name,
                    drug_index = ind,
                    p0 = p0_v,
                    shape1 = shape1,
                    shape2 = shape2,
                    prior = prior,
                    hpd_alpha = hpd_alpha,
                    alternative = "greater",
                    call = NULL)
    allResu[[i]] <- r
  }
  

  ret <- list(data = allResu, drug_index = drug_index)
  class(ret) <- c("borrow_multiple", "exchangeability_model")
  return(ret)
}
