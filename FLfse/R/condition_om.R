#' Function to create FLStock of the OM for MSE simulation
#'
#' @param stk object of class FLStock
#' @param idx object of class FLIndices
#' @param it iterations
#' @param ny number of years to project from initial year (iy = maxyear of stk + 1)
#' @param nsqy number of years to compute status quo metrics, default set at 3
#' @param qmod catchability submodel
#' @param fmod fishing mortality submodel
#' default ~te(replace(age, age>9,9), year, k=c(6,8))
#' @param mcsave mcmc parameter, default 100
#' @param sr stock-recruit model, default is bevholt
#'
#' @details creates the FLStock, script based on
#' http://www.flr-project.org/doc/An_introduction_to_MSE_using_FLR.html
#' It uses a4aSCA to run the assessment on existing FLStock and FLIndices objects
#' @return a list of "stk", "dy", "stk0", "srbh","srbh0","srbh.res","brp","Fmsy","msy","Bmsy","Bpa","Blim"
#'
#' @import FLa4a
#'
#' @export
create_FLStock <- function (stk, idx, it, qmod = NULL, # = list(~s(age, k=6)) smooting spline
                            fmod = NULL, # = ~te(replace(age, age>9,9), year, k=c(6,8) tensor spline
                            mcsave=100) {

  require(FLa4a)
  require(FLBRP)
  require(FLAssess)
  #require(FLash)

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

  return(list(stk,stk0))
}


