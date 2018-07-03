#' Fit an SR model to FLStock and create residuals for MSE projection
#'
#' @param stk FLStock
#' @param sr_model default is "bevholt"
#' @param method default is "L-BFGS-B", not fit for other methods yet
#' @param it number of iterations
#' @param iy intermediate year (first year of projection)
#' @param fy final year of projection
#'
#' @return list of median FLSotck, FLSR for FLStock and median FLStock for reference points and
#' residuals for projection
#' 
#' @export

om_sr_model <- function(stk, sr_model="bevholt", method="L-BFGS-B", it, iy, fy){
  
  # Reduce to keep one iteration only for reference points - could be an input but easy enough to calculate again
  stk0 <- qapply(stk, iterMedians)
  
  # Fit the stock-recruit model
  if (method == "L-BFGS-B"){
    srbh  <- fmle(as.FLSR(stk, model=sr_model), method = method, lower=c(1e-6, 1e-6), upper=c(max(rec(stk_om)) * 3, Inf))
    srbh0 <- fmle(as.FLSR(stk0, model=sr_model), method= method, lower=c(1e-6, 1e-6), upper=c(max(rec(stk_om)) * 3, Inf))
  }
  
  # Generate stock-recruit residuals for the projection period
  srbh.res  <- rnorm(it, FLQuant(0, dimnames=list(year=iy:fy)), mean(c(apply(residuals(srbh), 6, sd))))
  
  output_sr        <- list(srbh, srbh0, srbh.res)  
  names(output_sr) <- c("srbh", "srbh0", "srbh.res") 
  return(output_sr)
  
}


