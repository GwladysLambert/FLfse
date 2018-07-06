#' Function to create several iterations of FLStock
#'
#' @param stk object of class FLStock
#' @param idx object of class FLIndices
#' @param it iterations
#' @param qmod catchability submodel
#' @param fmod fishing mortality submodel
#' @param mcsave mcmc parameter, default 100
#' @param seed.nb set seed number
#'
#' @details creates \code{it} iterations of the FLStock \code{stk} using MCMC,
#' calls a4aSCA to run the assessment on existing FLStock \code{stk} and 
#' FLIndices \code{idx} objects
#' script based on
#' http://www.flr-project.org/doc/An_introduction_to_MSE_using_FLR.html
#' 
#' @return a list of "stk" and "stk0", i.e. the median of all iterations
#'
#' @importFrom FLa4a a4aSCA SCAMCMC
#'
#' @export
create_FLStock <- function (stk, idx, it, qmod = NULL, fmod = NULL,
                            mcsave=100, seed.nb = 321) {

  set.seed(seed.nb)
  
  require(FLa4a)

  mcmc <- it * mcsave

  # Fit the model - !!! need to find better way to specify wether to use default fmodel and qmodel or not !!!
  if (is.null(qmod) & !is.null(fmod))
  fit  <- a4aSCA(stk, idx, fmodel = fmod, fit = "MCMC",
                       mcmc = SCAMCMC(mcmc = mcmc, mcsave = mcsave, mcprobe = 0.4))
  if (!is.null(qmod) & is.null(fmod))
    fit  <- a4aSCA(stk, idx, qmodel = qmod, fit = "MCMC",
                          mcmc = SCAMCMC(mcmc = mcmc, mcsave = mcsave, mcprobe = 0.4))
  if (is.null(qmod) & is.null(fmod))
    fit  <- a4aSCA(stk, idx, fit = "MCMC",
                          mcmc = SCAMCMC(mcmc = mcmc, mcsave = mcsave, mcprobe = 0.4))
  if (!is.null(qmod) & !is.null(fmod))
    fit  <- a4aSCA(stk, idx, fit = "MCMC", qmodel = qmod, fmodel = fmod, 
                          mcmc = SCAMCMC(mcmc = mcmc, mcsave = mcsave, mcprobe = 0.4))
  
  # Update the FLStock object
  stk  <- stk + fit
  
  # Reduce to keep one iteration only for reference points
  stk0 <- qapply(stk, iterMedians)
  
  output_stk        <- list(stk,stk0)
  names(output_stk) <- c("stk","stk0")

  return(output_stk)
}


