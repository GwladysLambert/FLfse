#' Run the full MSE
#'
#' @param stk object of class FLStock
#' @param idx object of class FLindex
#' @param dy final data year of original FLStock
#' @param ny number of years to project from initial year (iy = maxyear of stk + 1)
#' @param nsqy number of years to compute status quo metrics, default set at 3
#' @param srbh stock-recruit model
#' @param srbh.res stock-recruit model residuals
#' @param assessment assessment model to be used, default is sam (will also be set up for spict)
#' @param Bpa reference point
#' @param Fmsy reference points
#' @param seed.nb sets seed
#'
#' @return mse runs
#' @export

run_mse <- function(stk, idx, dy,  ny, nsqy = 3, srbh, srbh.res,
                    assessment= "sam", Bpa, Fmsy,
                    seed.nb=321) {

  # warnings and stops
  if (range(stk)[2] != range(idx)[2]) stop('different max age in FLIndices and FLStock will not work')
  
  # set yrs
  y0 <- range(stk)["minyear"] # initial data year
  iy <- dy+1      # 1st yr of projection
  fy <- dy+ny     # total number of yrs
  vy <- ac(iy:fy) # vector of position of projected years

  # Here TAC in the final year of data is assumed to be the realised
  #catch in stk for the same year, while the TAC in the intermediate year is set equal to the TAC in the final year of data
  TAC <- FLQuant(NA, dimnames=list(TAC="all", year=c(dy,vy), iter=1:dim(stk)[6]))
  TAC[,ac(dy)] <- catch(stk)[,ac(dy)]
  TAC[,ac(iy)] <- TAC[,ac(dy)] #assume same TAC in the first intermediate year
  ctrl   <- getCtrl(c(TAC[,ac(iy)]), "catch", iy, dim(stk)[6]) # dim(stk.om)[6] is "it" (iterations)
  # Set up the operating model FLStock object
  stk.om <- fwd(stk, control=ctrl, sr=srbh, sr.residuals = exp(srbh.res), sr.residuals.mult = TRUE)

  ## NOW LOOP
  set.seed(seed.nb) # set seed to ensure comparability between different runs
  for(i in vy[-length(vy)]){
    # set up simulations parameters
    ay  <- an(i)
    cat(i, " > ")
    flush.console()
    vy0 <- 1:(ay-y0) # data years (positions vector) - one less than current year
    sqy <- (ay-y0-nsqy+1):(ay-y0) # status quo years (positions vector) - one less than current year

    # apply observation error
    oem <- observation_error_proj(stk.om, idx, i, vy0)
    stk.mp <- oem$stk
    idx.mp <- oem$idx
    idx    <- oem$idx.om

    # perform assessment
    #out.assess <- xsa(stk.mp, idx.mp)
    #stk.mp <- out.assess$stk

    if (assessment=="sam") {
      require(stockassessment)
      out.assess <- FLR_SAM(stk.mp, idx.mp)
      stk.mp     <- SAM2FLStock(out.assess)
    }
    if (assessment=="spict") {
      require(stockassessment)
      out.assess <- FLR_SPiCT(stk.mp, idx.mp)
      stk.mp     <- out.assess
    }
    # apply ICES MSY-like Rule to obtain Ftrgt
    # (note this is not the ICES MSY rule, but is similar)
    flag <- ssb(stk.mp)[,ac(ay-1)]<Bpa
    Ftrgt <- ifelse(flag,ssb(stk.mp)[,ac(ay-1)]*Fmsy/Bpa,Fmsy)

    # project the perceived stock to get the TAC for ay+1
    fsq.mp <- yearMeans(fbar(stk.mp)[,sqy]) # Use status quo years defined above
    ctrl   <- getCtrl(c(fsq.mp, Ftrgt), "f", c(ay, ay+1), dim(stk.om)[6])
    stk.mp <- stf(stk.mp, 2)
    gmean_rec <- c(exp(yearMeans(log(rec(stk.mp)))))
    stk.mp    <- fwd(stk.mp, control=ctrl, sr=list(model="mean", params = FLPar(gmean_rec,iter=dim(stk.om)[6])))
    TAC[,ac(ay+1)] <- catch(stk.mp)[,ac(ay+1)]

    # apply the TAC to the operating model stock
    ctrl   <- getCtrl(c(TAC[,ac(ay+1)]), "catch", ay+1, dim(stk.om)[6])
    stk.om <- fwd(stk.om, control=ctrl,sr=srbh, sr.residuals = exp(srbh.res), sr.residuals.mult = TRUE)
  }

}

