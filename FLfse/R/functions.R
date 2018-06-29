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
  catch.n(stk.tmp)  <- catch.n(stk.tmp) + 0.1
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
  dnms <- list(iter=1:it, year=years, c("min", "val", "max"))
  arr0 <- array(NA, dimnames=dnms, dim=unlist(lapply(dnms, length)))
  arr0[,,"val"] <- unlist(values)
  arr0 <- aperm(arr0, c(2,3,1))
  ctrl <- fwdControl(data.frame(year=years, quantity=quantity, val=NA))
  ctrl@trgtArray <- arr0
  ctrl
}
