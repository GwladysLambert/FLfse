#' Function to create FLIndices of the OM for MSE simulation
#'
#' @param idx object of class FLIndices
#' @param stk object of class FLStock
#' @param stk0 median of stk iterations created in \link[FLfse]{create_FLStock}
#' @param it iterations
#'
#' @details creates the FLIndices, script based on
#' http://www.flr-project.org/doc/An_introduction_to_MSE_using_FLR.html
#'
#' @return FLIndices
#'
#' @export

create_FLIndices <- function(idx, stk, stk0, it) {
  
  browser()
  
  # Estimate the index catchabilities from the a4a fit (without simulation)
  # Observation error is introduced through the index catchability-at-age
  # Set up the FLIndices object and populate it
  # (note, FLIndices potentially has more than one index, hence the for loop)
  
  #browser()
  
  idcs <- FLIndices()
  
  for (i in 1:length(idx)){
    
    #   Set up FLQuants and calculate mean and sd for catchability
    lst        <- mcf(list(index(idx[[i]]), stock.n(stk0))) # make FLQuants same dimensions
    
    # Temporary - if the above creates NAs du to difference in number of ages between stk and idx then copy last line of data over
    # This is only to run the cod example, would ahve to be carefully changed for real simulations!!!
    if (range(idx[[i]])[2] < range(stk0)[2]) {
      for (j in setdiff(seq(length=range(stk0)[2]), seq(length=range(idx[[i]])[2]))) {
        lst[[1]][j,] <- lst[[1]][c(j-1),]
      }
    }
    
    idx.lq     <- log(lst[[1]]/lst[[2]]) # log catchability of index
    
    #######################################################################################################
    ## IDEA - when we remove sites, we change catchability at age... should that go here somehow? ##
    #######################################################################################################
    
    idx.qmu    <- idx.qsig <- stock.n(iter(stk,1)) # create quants
    idx.qmu[]  <- yearMeans(idx.lq) # allocate same mean-at-age to every year
    idx.qsig[] <- sqrt(yearVars(idx.lq)) # allocate same sd-at-age to every year
    #   Build index catchability based on lognormal distribution with mean and sd calculated above
    idx.q      <- rlnorm(it, idx.qmu, idx.qsig)
    idx_temp   <- idx.q * stock.n(stk)
    idx_temp   <- FLIndex(index=idx_temp, index.q=idx.q) # generate initial index
    range(idx_temp)[c("startf", "endf")] <- c(0, 0) # timing of index (as proportion of year)
    idcs[[i]] <- idx_temp
  }
  names(idcs) <- names(idx)
  
  idx <- idcs #[1]
  
  return(idx)
}
