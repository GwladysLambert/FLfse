#' Fit an SR model to FLStock and create residuals for MSE projection
#'
#' @param stk FLStock - if NULL the \code{sr_model} will be used to create 
#' projected residuals
#' @param srbh.res residuals of FLSR - default is NULL but if given, 
#' no need for \code{stk}
#' @param sr_model default is \code{"bevholt"}
#' @param method default is \code{"L-BFGS-B"}, not fit for other methods yet
#' @param it number of iterations
#' @param iy intermediate year (first year of projection)
#' @param fy final year of projection
#' @param seed.nb set seed number
#'
#' @return a list of "srbh" (FLSR), "srbh0", i.e. the median of all iterations, and
#' "srbh.res", i.e. recruitment residuals for projection period if \code{stk} given,
#' or only "srbh.res" if \code{srbh.res} given
#' 
#' @export

om_sr_model <- function(stk=NULL, srbh.res=NULL, sr_model="bevholt", method="L-BFGS-B", it, iy, fy, seed.nb = 321){
  
 # browser()
  
  srbh  <- NULL
  srbh0 <- NULL
  
  if (!is.null(stk)) {
  # Reduce to keep one iteration only for reference points - could be an input but easy enough to calculate again
  stk0 <- qapply(stk, iterMedians)
  
  # Fit the stock-recruit model
  if (method == "L-BFGS-B"){
    srbh  <- fmle(as.FLSR(stk, model=sr_model), method = method, lower=c(1e-6, 1e-6), upper=c(max(rec(stk)) * 3, Inf))
    srbh0 <- fmle(as.FLSR(stk0, model=sr_model), method= method, lower=c(1e-6, 1e-6), upper=c(max(rec(stk)) * 3, Inf))
    srbh.res <- residuals(srbh)
  }
  }
  
  # Generate stock-recruit residuals for the projection period
  set.seed(seed.nb)
  srbh.res  <- rnorm(it, FLQuant(0, dimnames=list(year=iy:fy)), mean(c(apply(srbh.res, 6, sd))))

  output_sr        <- list(srbh, srbh0, srbh.res)  
  names(output_sr) <- c("srbh", "srbh0", "srbh.res") 
  return(output_sr)
  
}


