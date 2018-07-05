#' Function to include error in the projected index
#'
#' @param stk object of class FLStock
#' @param idx object of class FLindex
#' @param assessmentYear latest year of the assessment
#' @param dataYears position vector of years with data
#'
#' @return a list of stock, new index and index of the OM
#' @export

observation_error_proj <- function(stk, idx, assessmentYear, dataYears) {
  # dataYears is a position vector, not the years themselves
  stk.tmp           <- stk[, dataYears]
  # add small amount to avoid zeros
  catch.n(stk.tmp)  <- catch.n(stk.tmp) + 0.001
  # Generate the indices - just data years
  idx.tmp           <- lapply(idx, function(x) x[,dataYears])
  # Generate observed index
  for (i in 1:length(idx)) {
    index(idx[[i]])[, assessmentYear] <-
      stock.n(stk)[, assessmentYear]*index.q(idx[[i]])[, assessmentYear]
  }
  
  return(list(stk=stk.tmp, idx=idx.tmp, idx.om=idx))
}


#' Define control object for \link[FLash]{fwd}
#'
#' @param values 
#' @param quantity 
#' @param years 
#' @param it iterations
#'
#' @details The fwd method from FLash needs a control object, which is set by this function.
#'
#' @return ctrl for fwd
#' @export
#' 
getCtrl <- function(values, quantity, years, it){
  #browser()
  dnms <- list(iter=1:it, year=years, c("min", "value", "max"))
  arr0 <- array(NA, dimnames=dnms, dim=unlist(lapply(dnms, length)))
  arr0[,,"value"] <- unlist(values)
  arr0 <- aperm(arr0, c(2,3,1))
  #ctrl <- fwdControl(data.frame(year=years, quantity=quantity, val=NA))
  #ctrl@trgtArray <- arr0
  ctrl <- fwdControl(list(year=years, quant=quantity, value=values))
  ctrl
}

### ------------------------------------------------------------------------ ###
### additional functions (internal) ####
### ------------------------------------------------------------------------ ###


### ------------------------------------------------------------------------ ###
### adapt dimensions in 2 FLR objects ####
### ------------------------------------------------------------------------ ###
### returned as list
adapt_dims <- function(obj1, obj2, fill.iter = FALSE) {

  ### save objects in list
  res <- list(obj1, obj2)

  ### min year
  min_year <- min(as.numeric(sapply(lapply(res, dims), "[", "minyear")))
  ### max year
  max_year <- max(as.numeric(sapply(lapply(res, dims), "[", "maxyear")))

  ### adapt year range
  res <- lapply(res, window, start = min_year)
  res <- lapply(res, window, end = max_year)

  ### number of iterations
  n_iter <- as.numeric(sapply(lapply(res, dims), "[", "iter"))

  ### check comparability of iterations and stop if not
  if (!identical(n_iter[1], n_iter[2])) {

    if (!any(n_iter == 1)) stop("incompatible iter dimension")

    ### adapt
    res[[which(n_iter == 1)]] <- propagate(res[[which(n_iter == 1)]],
                                           n_iter[which(n_iter != 1)],
                                    fill.iter = fill.iter)

  }

  return(res)

}


### ------------------------------------------------------------------------ ###
### Functions attached to creating stocks ####
### ------------------------------------------------------------------------ ###

# oneWayTrip {{{
oneWayTrip <- function(stk, sr, brp, fmax=refpts(brp)['crash', 'harvest'] * 0.80,
                       years=seq(dims(stk)$minyear + 1, dims(stk)$maxyear),
                       residuals=FLQuant(1, dimnames=dimnames(rec(stk))),f0=NULL) {
  
  # limits
  if(is.null(f0)){ 
    f0 <- c(fbar(stk)[,1])
  } else{
    f0=as.vector(f0)
  }
  
  fmax <- c(fmax)
  rate <- exp((log(fmax) - log(f0)) / (length(years)))
  
  # linear trend: f <- rate ^ (seq(0, length(years))) * f0
  fs <- (matrix(rate, nrow=length(years) + 1, ncol=length(rate)) ^
           matrix(seq(0, length(years)), nrow=length(years) + 1, ncol=length(rate)) *
           matrix(f0, nrow=length(years) + 1, ncol=length(rate)))[-1,]
  
  # fwdControl
  ftar <- FLQuant(c(fs), dimnames=list(year=years, iter=seq(length(rate))), quant='age')
  
  
  # fwd
  res <- fwd(stk, control=as(FLQuants(f=ftar), "fwdControl"), sr=sr, residuals=residuals)
  
  return(res)
} # }}}
