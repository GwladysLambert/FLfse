#' Function to create FLStock of the OM for MSE simulation
#'
#' @param stk object of class FLStock
#' @param idx object of class FLIndices
#' @param it iterations
#' @param qmod catchability submodel
#' @param fmod fishing mortality submodel
#' default ~te(replace(age, age>9,9), year, k=c(6,8))
#' @param mcsave mcmc parameter, default 100
#' @param seed.nb set seed number
#'
#' @details creates the FLStock, script based on
#' http://www.flr-project.org/doc/An_introduction_to_MSE_using_FLR.html
#' It uses a4aSCA to run the assessment on existing FLStock and FLIndices objects
#' @return a list of "stk", "dy", "stk0", "srbh","srbh0","srbh.res","brp","Fmsy","msy","Bmsy","Bpa","Blim"
#'
#'
#' @export
create_FLStock <- function (stk, idx, it, qmod = NULL, # = list(~s(age, k=6)) smooting spline
                            fmod = NULL, # = ~te(replace(age, age>9,9), year, k=c(6,8) tensor spline
                            mcsave=100, seed.nb = 321) {

  set.seed(seed.nb)
  
  require(FLa4a)

  mcmc <- it * mcsave

  # Fit the model - need to find better way to specify wether to use default fmodel and qmodel or not
  if (is.null(qmod) & !is.null(fmod))
  fit  <- FLa4a::a4aSCA(stk, idx, fmodel = fmod, fit = "MCMC",
                       mcmc = SCAMCMC(mcmc = mcmc, mcsave = mcsave, mcprobe = 0.4))
  if (!is.null(qmod) & is.null(fmod))
    fit  <- FLa4a::a4aSCA(stk, idx, qmodel = qmod, fit = "MCMC",
                          mcmc = SCAMCMC(mcmc = mcmc, mcsave = mcsave, mcprobe = 0.4))
  if (is.null(qmod) & is.null(fmod))
    fit  <- FLa4a::a4aSCA(stk, idx, fit = "MCMC",
                          mcmc = SCAMCMC(mcmc = mcmc, mcsave = mcsave, mcprobe = 0.4))
  if (!is.null(qmod) & !is.null(fmod))
    fit  <- FLa4a::a4aSCA(stk, idx, fit = "MCMC", qmodel = qmod, fmodel = fmod, 
                          mcmc = SCAMCMC(mcmc = mcmc, mcsave = mcsave, mcprobe = 0.4))
  
  # Update the FLStock object
  stk  <- stk + fit
  
  # Reduce to keep one iteration only for reference points
  stk0 <- qapply(stk, iterMedians)
  
  output_stk        <- list(stk,stk0)
  names(output_stk) <- c("stk","stk0")

  return(output_stk)
}


