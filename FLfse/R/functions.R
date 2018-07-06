#' Define control object for \link[FLasher]{fwd}
#'
#' @param values target values
#' @param quantity 'f' or 'catch' or...
#' @param years years that ctrl needed for
#' @param it number of iterations
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
  values_org <- arr0[,"val",1:it] # for fwdControl with FLasher...
  if (is.matrix(values_org)) values_org <- unlist(c(as.data.frame(arr0[,"val",1:it]))) # for fwdControl with FLasher...
  ctrl <- fwdControl(list(year=years, quant=quantity, value=values_org))
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
### Functions attached to creating stocks #### NOT USED HERE YET !!!!!!!!!!!!!!!!!!
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
