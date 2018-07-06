#' Function to create several iterations of FLIndices
#'
#' @param idx object of class FLIndices - can be NULL if one does not exist yet,
#' it will be created inside the function using \code{qmod} and \code{qmod_pars} 
#' @param stk object of class FLStock
#' @param qmod default model to create idx from scratch, i.e. logistic function
#' @param qmod_pars default model \code{qmod} parameters, i.e. logistic function with
#' inflection point of curve = 10% of max age
#' @param add.error whether or not the iterated indices should include some random noise, 
#' set at FALSE by default
#' @param assessment.yr latest assessment year for which an index is to be estimated,
#' default is NULL if indices for all years of \code{stk} are to be estimated
#' @param it iterations, set to NULL as not needed if \code{assessment.yr} 
#' is provided (see \code{Details})
#' @param seed.nb set seed number
#' 
#' @details creates the FLIndices, script based on
#' http://www.flr-project.org/doc/An_introduction_to_MSE_using_FLR.html
#' If asessment.yr is given, it will assume that indices were already produced up 
#' to taht year and that an index.q is available and only estimate FLIndices 
#' for that year
#'
#' @return FLIndices
#'
#' @export

create_FLIndices <- function(idx, stk, it=NULL, 
                             qmod='~max_q/(1+exp(-steepness*(age - age50)))', 
                             qmod_pars=FLPar(max_q = 1, steepness = 1, age50 = max(ages)/10),
                             add.error = FALSE,
                             assessment.yr =NULL, seed.nb = 321) {
  
  # Sets up the FLIndices object and populate it
  # (note, FLIndices potentially has more than one index, hence the for loop)
  # Estimate the index catchabilities from the stk fit 
  # Observation error is introduced through the index catchability-at-age (can be optional)
  
  #browser()
  
  idcs <- FLIndices()
  
  ### ~~~ if there was no idx to start with, create one for the first stk iteration~~~ ###
  
    if (is.null(idx)) { # assessment.yr should be NULL too here...
      
      ### get ages
      ages <- an(dimnames(stock.n(iter(stk,1)))[["age"]])
      
      ### define model for selectivity
      q_model <- FLModelSim(model = formula(qmod), params = qmod_pars)
      ### model selectivity
      q_modeled <- predict(q_model, age = ages)
      
      ### create index template
      idx <- iter(FLIndex(index = stock.n(iter(stk,1))),1)
      
      ### insert selectivity
      index.q(idx) <- c(q_modeled)
      
      ### calculate historical index values
      index(idx) <- index.q(idx) * iter(stock.n(iter(stk,1)),1) * iter(stock.wt(iter(stk,1)),1)
      
      ### save as FLIndices
      idx = FLIndices(idx = idx)
    }
  
  
  ### ~~~ create index for each iteration and add uncertainty if add.error = TRUE ~~~ ###
  
  for (i in 1:length(idx)){
    
    if (is.null(assessment.yr)) {
      
      #   Set up FLQuants and calculate mean and sd for catchability
      stk0       <- qapply(stk, iterMedians)
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
      set.seed(seed.nb)
      if (add.error ==T){
        idx.q      <- rlnorm(it, idx.qmu, idx.qsig)
      }  else { idx.q <- propagate(exp(idx.qmu), it) }    
      
    } else {
      idx.q    <- index.q(idx[[i]])
    }
    
    idx_temp   <- idx.q * stock.n(stk)
    idx_temp   <- FLIndex(index=idx_temp, index.q=idx.q) # generate initial index
    range(idx_temp)[c("startf", "endf")] <- c(0, 1) # timing of index (as proportion of year)
    idcs[[i]] <- idx_temp
  }
  names(idcs) <- names(idx)
  
  idx <- idcs #[1] - not sure why it was extrating [1] in example script
  
  # other outputs needed  
  stk.tmp <- NULL
  idx.tmp <- NULL
  
  if (!is.null(assessment.yr)) {
    stk.tmp <- window(stk, end=an(assessment.yr)) 
    idx.tmp <- window(idx, end=an(assessment.yr)) 
  }
  
  return(list(idx.om = idx, idx = idx.tmp, stk=stk.tmp))
}

