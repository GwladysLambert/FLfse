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
create_FLStock <- function (stk, idx, it, ny, nsqy = 3, qmod, # = list(~s(age, k=6)) smooting spline
                            fmod, # = ~te(replace(age, age>9,9), year, k=c(6,8) tensor spline
                            mcsave=100, sr="bevholt") {

  require(FLa4a)
  require(FLBRP)
  require(FLAssess)
  #require(FLash)

  y0 <- range(stk)["minyear"] # initial data year
  dy <- range(stk)["maxyear"] # final data year
  iy <- dy+1  # initial year of projection (also intermediate year)
  fy <- dy+ny # final year of projection

  mcmc <- it * mcsave

  # Fit the model - need to find better way to specify wether to use default fmodel and qmodel or not
  if (missing(qmod) & !missing(fmod))
  fit  <- FLa4a::a4aSCA(stk, idx, fmodel = fmod, fit = "MCMC",
                       mcmc = SCAMCMC(mcmc = mcmc, mcsave = mcsave, mcprobe = 0.4))
  if (!missing(qmod) & missing(fmod))
    fit  <- FLa4a::a4aSCA(stk, idx, qmodel = qmod, fit = "MCMC",
                          mcmc = SCAMCMC(mcmc = mcmc, mcsave = mcsave, mcprobe = 0.4))
  if (missing(qmod) & missing(fmod))
    fit  <- FLa4a::a4aSCA(stk, idx, fit = "MCMC",
                          mcmc = SCAMCMC(mcmc = mcmc, mcsave = mcsave, mcprobe = 0.4))
  if (!missing(qmod) & !missing(fmod))
    fit  <- FLa4a::a4aSCA(stk, idx, fit = "MCMC", qmodel = qmod, fmodel = fmod, 
                          mcmc = SCAMCMC(mcmc = mcmc, mcsave = mcsave, mcprobe = 0.4))
  
  # Update the FLStock object
  stk  <- stk + fit
  # Reduce to keep one iteration only for reference points
  stk0 <- qapply(stk, iterMedians)

  # Fit stock-recruit model
  # A Beverton-Holt stock-recruit model is fitted for each iteration, with residuals
  # generated for the projection window based on the residuals from the historic period.
  # A stock-recruit model is also fitted to the "median" stk for reference points.

  # Fit the stock-recruit model
  srbh  <- fmle(as.FLSR(stk, model=sr), method="L-BFGS-B",
               lower=c(1e-6, 1e-6), upper=c(max(rec(stk)) * 3, Inf))
  srbh0 <- fmle(as.FLSR(stk0, model=sr), method="L-BFGS-B",
                lower=c(1e-6, 1e-6), upper=c(max(rec(stk)) * 3, Inf))
  # Generate stock-recruit residuals for the projection period
  srbh.res <- rnorm(it, FLQuant(0, dimnames=list(year=iy:fy)), mean(c(apply(residuals(srbh), 6, sd))))

  # Calculate reference points and set up the operating model for the projection window
  # Reference points based on the "median" stk, assuming (for illustrative purposes only)
  # that Bpa=0.5Bmsy and Blim=Bpa/1.4. The stf method is applied to the operating model stk object
  # in order to have the necessary data (mean weights, etc.) for the projection window.

  # Calculate the reference points
   brp  <- FLBRP::brp(FLBRP(stk0, srbh0))
   Fmsy <- c(FLBRP::refpts(brp)["msy","harvest"])
   msy  <- c(FLBRP::refpts(brp)["msy","yield"])
   Bmsy <- c(FLBRP::refpts(brp)["msy","ssb"])
   Bpa  <- 0.5*Bmsy
   Blim <- Bpa/1.4
  ## Prepare the FLStock object for projections
   stk  <- stf(stk, fy-dy, nsqy, nsqy)

   output_list <- list(stk,dy, stk0, srbh,srbh0,srbh.res,brp,Fmsy,msy,Bmsy,Bpa,Blim)
   names(output_list) <- c("stk","dy","stk0", "srbh","srbh0","srbh.res","brp","Fmsy","msy","Bmsy","Bpa","Blim")

   return(output_list)
}



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

  # Estimate the index catchabilities from the a4a fit (without simulation)
  # Observation error is introduced through the index catchability-at-age
  # Set up the FLIndices object and populate it
  # (note, FLIndices potentially has more than one index, hence the for loop)

  idcs <- FLIndices()

  for (i in 1:length(idx)){

    #   Set up FLQuants and calculate mean and sd for catchability
    lst        <- mcf(list(index(idx[[i]]), stock.n(stk0))) # make FLQuants same dimensions
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

  idx<-idcs[1]

  return(idx)
}


