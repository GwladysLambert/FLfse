#' Run the MSE loop
#'
#' @param stk object of class FLStock
#' @param idx object of class FLindex
#' @param dy final data year of original FLStock
#' @param nsqy number of years to compute status quo metrics, default set at 3
#' @param srbh stock-recruit model
#' @param srbh.res stock-recruit model residuals
#' @param assessment assessment model to be used, default is sam (will also be set up for spict)
#' @param Bpa reference point
#' @param Fmsy reference points
#' @param seed.nb sets seed
#' @param y0 intial data year
#' @param iy intermediate year
#' @param fy total number of years
#'
#' @importFrom FLasher fwdControl
#' 
#' @return mse runs
#' @export

mp_fn <- function(stk, idx, 
                  y0, iy, dy, fy, nsqy, 
                  srbh, srbh.res,
                  assessment= "sam", Bpa, Fmsy, 
                  seed.nb=321) {
  
 # browser()
  
  ############################
  ### ~~~ General Prep ~~~ ###
  ############################
  
  # set seed to ensure comparability between different runs
  set.seed(seed.nb) 
  # warnings and stops
  if (range(stk)[2] != range(idx)[2]) stop('different max age in FLIndices and FLStock will not work')
  # vector of position of projected years
  vy <- ac(iy:fy) 
  
  ###############################
  ### ~~~ Preparation HCR ~~~ ###
  ###############################
  
  # Empty FLQuant of TACs 
  TAC <- FLQuant(NA, dimnames=list(TAC="all", year=c(dy,vy), iter=1:dim(stk)[6])) #
  # First TAC is assumed to be the realised catch in stk for the same year
  TAC[,ac(dy)] <- catch(stk)[,ac(dy)]
  # Define stk.om
  stk.om <- stk
  
  ####################
  ### ~~~ Loop ~~~ ###
  ####################
  
  for(i in ac(an(vy)-1)){ #[-length(vy)] 
    
    # start running on last year of data to predict catch or F to use in intermediate year
    ## ~~~~~~~~~~~~~~~~~~~~ ##
    ## !!! NOT SURE HOW THE HCR WORKS REALLY SO MADE THIS UP... NEED TO BE CHECKED AND CORRECTED !!!
    ## +++ NNED TO RE ORGANISE + CREATE ARGUMENTS SO HCR CAN BE SCENARIOS 
    ## ~~~~~~~~~~~~~~~~~~~~ ##
    # set up simulations parameters
    ay  <- an(i)
    cat(i, " > ")
    flush.console()
    vy0 <- 1:(ay-y0) # data years (positions vector) - one less than current year
    sqy <- (ay-y0-nsqy+1):(ay-y0) # status quo years (positions vector) - one less than current year
    
    #################################
    ### ~~~ Observation model ~~~ ###
    #################################
    
    oem <- create_FLIndices(idx = idx, stk = stk.om, assessment.yr = i, seed.nb = seed.nb)
    stk.mp <- oem$stk
    idx.mp <- oem$idx
    idx    <- oem$idx.om
    
    ##########################
    ### ~~~ Assessment ~~~ ###
    ##########################
    
    # only "sam" is set up for now
    if (assessment=="sam") {
      out.assess <- FLR_SAM(stk.mp, idx.mp)
      stk.mp     <- SAM2FLStock(out.assess)
    }
    # if (assessment=="spict") {
    #   out.assess <- FLR_SPiCT(stk.mp, idx.mp)
    #   stk.mp     <- out.assess
    # } 
    # if (assessment=="xsa") {
    #   out.assess <- xsa(stk.mp, idx.mp)
    #   stk.mp <- out.assess$stk
    # }
    
    # keep only the iterations that did not fail
    if (dim(stk.mp)[6]!=1) # NEED TO AMKE SURE out.assess IS IN THE RIGHT LIST FORMAT EVEN IF JUST ONE ITERATION WORKING !!!! for now just set i == vy[1], might be good enough
    {
      flag <- unlist(lapply(out.assess, function(x) is(x)[1]=="sam")) 
      stk.mp <- iter(stk.mp, flag)
      stk.om <- iter(stk.om, flag)
      idx.mp <- iter(idx.mp, flag)
      idx    <- iter(idx, flag)
      TAC    <- iter(TAC, flag)
      srbh.res <- iter(srbh.res, flag)
    }
    
    ###################
    ### ~~~ HCR ~~~ ###
    ###################
    
    # apply ICES MSY-like Rule to obtain Ftrgt
    # (note this is not the ICES MSY rule, but is similar)
    flag  <- c(ssb(stk.mp)[,ac(ay-1)])<Bpa
    Ftrgt <- ifelse(flag,c(ssb(stk.mp)[,ac(ay-1)]*Fmsy/Bpa),Fmsy)
    
    # project the perceived stock to get the TAC for ay+1
    fsq.mp    <- yearMeans(fbar(stk.mp)[,sqy]) # Use status quo years defined above
    ctrl      <- getCtrl(c(fsq.mp, Ftrgt), "f", c(ay, ay+1),dim(stk.mp)[6])
    stk.mp    <- FLasher::stf(stk.mp, 2)
    gmean_rec <- c(exp(yearMeans(log(rec(stk.mp)))))
    stk.mp    <- fwd(stk.mp, control=ctrl, sr=list(model=mean, params = FLPar(gmean_rec,iter=dim(stk.mp)[6])))
    TAC[,ac(ay+1)] <- catch(stk.mp)[,ac(ay+1)]
    
    #################################
    ### ~~~ APPLY AND PROJECT ~~~ ###
    #################################
    
    # apply the TAC to the operating model stock
    ctrl   <- getCtrl(c(TAC[,ac(ay+1)]), "catch", ay+1, dim(stk.om)[6])
    stk.om <- fwd(stk.om, control=ctrl, sr=srbh, residuals = exp(srbh.res)) #, sr.residuals.mult = TRUE ## HAD TO REMOVE THE RESIDUALS MUL = TRUE HERE AS MOVING FROM FLash TO FLasher, IS THERE AN EQUIVALENT???
    
  }
  return(stk.om) # !!! need to set up a tracking FLQuant to extract all info needed along the way (i.e. TAC, Ftrgt etc) !!!
}


