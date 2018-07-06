#' Run the full MSE
#'
#' @param scens name of scenario
#' @param stock parameters controling the set up of the FLStock
#' \code{stock} = list(\code{stk}, \code{latin.name},
#' \code{common.name}, \code{region}, \code{biol.params}, \code{h})
#' \code{stk} if using an existing stock assessment - see \code{\link{create_FLStock}} 
#' \code{latin.name} if using biological parameters - see \code{\link{create_FLStock_biol}} 
#' \code{common.name} if using biological parameters - see \code{\link{create_FLStock_biol}} 
#' \code{region} if using biological parameters - see \code{\link{create_FLStock_biol}} 
#' \code{biol.params} if using biological parameters - see \code{\link{create_FLStock_biol}} 
#' \code{h}  if using biological parameters - see \code{\link{create_FLStock_biol}} 
#' @param ctrl.om parameters controling the set up of the Operating Model
#' \code{ctrl.om} = list(\code{sr_init}, \code{fmod_init}, \code{qmod_init}, \code{mcmc_init})
#' \code{sr_init} stock-recruit model for the OM set up, default is bevholt
#' \code{qmod_init} catchability submodel for the OM set up
#' \code{fmod_init} fishing mortality submodel for the OM set up
#' \code{mcmc_init} number of mcmc for the OM set up
#' @param ctrl.idx parameters controling the set up of the Observation Model
#' \code{ctrl.idx} = list(\code{idx}, \code{qmod_idx}, \code{qmod_pars_idx}, \code{add.error_idx})
#' \code{idx} if using an existing stock assessment - see \code{\link{create_FLStock}} 
#' \code{qmod_idx} defines idx selectivity function if idx does not exist as in 
#' \code{\link{create_FLStock_biol}} approach. Default is NULL, meaning to use default 
#' settings in \code{\link{create_FLIndices}}
#' \code{qmod_pars_idx} defines idx selectivity parameters for \code{qmod_idx}
#' \code{add.error_idx} whether or not the iterated indices should include some random noise, 
#' set at FALSE by default
#' @param ctrl.sims parameters controling the set up of the Simulations properties
#' \code{ctrl.sims} = list(\code{it}, \code{ny}, \code{nsqy}, \code{seed.nb})
#' \code{ny number} of years to project from initial year (iy = maxyear of stk + 1), default set at 3
#' \code{it number} of iterations, default set at 3
#' \code{nsqy number} of years to compute status quo metrics, default set at 3
#' \code{param} seed.nb set seed number
#' @param ctrl.mp parameters controling the set up of the Management Procedure
#' \code{assessment} assessment model to be used, default is \code{sam} 
#' (will also be set up for \code{spict})
#' 
#' @details can run an MSE based on existing stock, using arguments \code{stk} and \code{idx} 
#' (\code{\link{create_FLStock}}), or biological parameters using \code{latin.name}, \code{common.name}, 
#' \code{area} or \code{params} (see \code{\link{create_FLStock_biol}}). Although 
#' \code{\link{create_FLStock_biol}} can take several several stocks at once, this MSE function can 
#' only be run for one stock at a time for now (there will be an option to use scenarios later on to 
#' run several stocks in parallel)
#' 
#' @return mse runs
#' 
#' @importFrom FLasher stf
#' 
#' @export
#' @examples 
#' \dontrun{
#' 
#' ## Run ple4 example 
#' data(ple4)
#' data(ple4.index)
#' stk   = ple4[1:8,] # needed to adjust for different number of ages between stk and idx
#' idx   = FLIndices(idx=ple4.index)
#' test_ple  <-  mse_base(stock = list(stk = ple4), 
#'                        ctrl.idx=list(idx = FLIndices(idx=ple4.index)), 
#'                        ctrl.sims = list(it = 2, ny = 3))
#' 
#' ## Run cod example
#' data("cod4_stk")
#' cod4_stk <- window(cod4_stk, end=2016) # no catch in 2017, crashes the simulations as needed to set TAC
#' data("cod4_idx")
#' cod4_idx <- window(cod4_idx, end=2016)
#' # no uncertainty around FLIndices (perfect knowledge)
#' test_cod  <-  mse_base(stock = list(stk = cod4_stk), ctrl.idx=list(idx = cod4_idx),
#'                        ctrl.sims = list(it = 2, ny = 3))
#' # uncertainty around FLIndices (observation error)
#' test_cod_re  <-  mse_base(stock = list(stk = cod4_stk), ctrl.idx=list(idx = cod4_idx, add.error_idx=TRUE),
#'                           ctrl.sims = list(it = 2, ny = 3))
#' 
#' ## Run params example
#' # no uncertainty around FLIndices (perfect knowledge)
#' test_params <- mse_base(stock = list(biol.params=data.frame(linf=100)), 
#'                         ctrl.sims = list(it = 2, ny = 3, seed.nb=123))
#' # uncertainty around FLIndices (observation error)
#' test_params_re <- mse_base(stock = list(biol.params=data.frame(linf=100)), 
#'                            ctrl.idx=list(add.error_idx=TRUE), 
#'                            ctrl.sims = list(it = 2, ny = 3, seed.nb=123))
#' 
#' ## Run WKLife example (perfect knowledge)
#' test_wklife <- mse_base(stock = list(latin.name="Clupea harengus"), 
#'                         ctrl.sims = list(it=2, ny=3, seed.nb= 5))
#' 
#' }
#'

mse_base <- function(
  scens     = list(id = "default", stringsAsFactors = FALSE), # not used anywhere yet but will be attached to the outputs !!!
  stock     = list(stk = NULL, 
                   latin.name=NULL, common.name=NULL, region=NULL,
                   biol.params = NULL, h=0.75),
  ctrl.om   = list(sr_init="bevholt", fmod_init=NULL, qmod_init=NULL, mcmc_init =100),
  ctrl.idx  = list(idx = NULL, qmod_idx = NULL, qmod_pars_idx = NULL, add.error_idx = FALSE),
  ctrl.sims = list(it=3, ny=3, nsqy = 3, seed.nb= 321),
  ctrl.mp   = list(assessment="sam")) {
  
  ###########################################################################
  ###             ~~~ ADD IN VALUES BY GIVEN VALUES ~~~                   ###
  ###########################################################################
  
  ## This is needed to complete the list structure of the inputs with default values 
  ## (we can go back to "no list" but this looks tidier, i.e. easier  to visualise
  ## what paramters fits where in the scenario process and see where possible to add structure and options)
  
  names_args  <- unlist(lapply(formals(mse_base), function(x) names(formula(x))))
  values      <- unlist(lapply(formals(mse_base), formula), recursive = F)
  names(values) <- names_args
  
  new.args <- c(scens,stock,ctrl.om,ctrl.idx,ctrl.sims,ctrl.mp) 
  values[names(new.args)] <- new.args
  
  options(warn = -1)
  attach(values, warn.conflicts= FALSE)
  options(warn = 0)
  
  #############################################################
  ###             ~~~ OPERATING MODEL ~~~                   ###
  #############################################################
  
  ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
  ### ~~~ IF STARTING FROM AN EXISTING STOCK ASSESSMENT ~~~ ###
  ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
  
  if (!is.null(stk) & !is.null(idx)) {
    
    ### ~~ Set up the simulations ~~ ###
    
    y0 <- range(stk)["minyear"] # initial data year
    dy <- range(stk)["maxyear"] # final data year
    iy <- dy+1  # initial year of projection (also intermediate year)
    fy <- dy+ny # final year of projection
    
    
    ### ~~ Set up the OM ~~ ###
    
    ## This is to create an FLStock from a real stock
    # results in a list of 2 - one FLStock with iterations, one FLStock median for ref points calculations later
    stk_om       <-   create_FLStock(stk =  stk, idx = idx, it = it, 
                                     fmod = fmod_init, qmod = qmod_init, mcsave = mcmc_init, seed.nb = seed.nb) 
    
    ### ~~ Set up the Stock-Recruitment model ~~ ###
    
    # Fit stock-recruit model
    # A Beverton-Holt stock-recruit model is fitted for each iteration, with residuals
    # generated for the projection window based on the residuals from the historic period.
    # A stock-recruit model is also fitted to the "median" stk for reference points.
    sr_om        <-   om_sr_model(stk_om$stk, sr_model=sr_init, it =  it, iy = iy, fy = fy, seed.nb = seed.nb)
    
    ###   ~~ Calculate the reference points ~~   ###
    
    # Calculate reference points and set up the operating model for the projection window
    # Reference points based on the "median" stk, assuming (for illustrative purposes only)
    # that Bpa=0.5Bmsy and Blim=Bpa/1.4. The stf method is applied to the operating model stk object
    # in order to have the necessary data (mean weights, etc.) for the projection window.
    refpts_om    <-   ref_pts(stk_om$stk0, sr_om$srbh0)
    
  }
  
  ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
  ### ~~~ IF STARTING FROM PARAMETERS ~~~ ###
  ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
  
  if (is.null(stk) & is.null(idx)) {
    
    ### ~~ Set up the OM ~~ ###
    all_om <- create_FLStock_biol(latin.name=latin.name, common.name=common.name, region=region, 
                                  biol.params = biol.params, h=h, it = it, seed.nb=seed.nb) 
    
    ### ~~ Extract stk ~~ ###
    stk_om    <- lapply(all_om, '[[', "stk")
    stk0      <- lapply(stk_om, function(x) qapply(x, iterMedians))
    stk_om    <- list(stk=stk_om[[1]], stk0=stk0[[1]])
    
    ### ~~ Extract stock-recruit ~~ ###
    sr_om     <- lapply(all_om, '[[', "sr")
    srbh0     <- lapply(sr_om, function(x) qapply(x, iterMedians))
    sr_om     <- list(srbh=sr_om[[1]], srbh0=srbh0[[1]])
    
    ### ~~ Extract ref points ~~ ###
    refpts_om <- lapply(all_om, '[[', "brp")[[1]]
    refpts_om <- ref_pts(brps=refpts_om)
    
    ### ~~ Set up the simulations ~~ ###
    y0 <- range(stk_om$stk)["minyear"] # initial data year
    dy <- range(stk_om$stk)["maxyear"] # final data year
    iy <- dy+1  # initial year of projection (also intermediate year)
    fy <- dy+ny # final year of projection
    
    ### ~~ Projected stock-recruit residuals ~~ ###
    # Generate stock-recruit residuals for the projection period
    srbh.res.tmp   <- lapply(all_om, '[[', "sr_res")[[1]]
    sr_om$srbh.res <- om_sr_model(srbh.res = srbh.res.tmp, it =  it, iy = iy, fy = fy, seed.nb = seed.nb)$srbh.res
    
    ### ~~ Indices ~~ ###
    ## specify that the idx does not exist yet so will have to be created
    idx = NULL
    
  }
  
  
  #############################################################
  ###             ~~~ OBSERVATION MODEL ~~~                 ###
  #############################################################
  
  ### ~~ Projected stk ~~ ###
  ## First prepare the FLStock object for projections - this provides wts and fbars for fy-dy yrs, 
  #  averageing over last nsqy yrs of data
  stk_om$stk   <-   stf(stk_om$stk, fy-dy, nsqy, nsqy)
  
  ### ~~ Indices ~~ ###
  ## This is to create an FLindices
  ## Not giving any options for how to create those indices but this is where most development/innovation should occur in FLfse
  ## !!!!!!!!THERE HAS GOT TO BE A BETTER WAY TO SET UP THE MODELS WITH USING DEFAULT PARAMETERS OR NOT.... !!!!!!
  if (is.null(qmod_idx)){
    idx_oem     <-   create_FLIndices(idx =  idx, stk = stk_om$stk, it = it, add.error = add.error_idx, seed.nb = seed.nb)$idx.om
  } else {
    idx_oem     <-   create_FLIndices(idx =  idx, stk = stk_om$stk, it = it, 
                                      qmod = qmod_idx, qmod_pars = qmod_pars_idx,
                                      add.error = add.error_idx, 
                                      seed.nb = seed.nb)$idx.om
  }
  
  
  #############################################################
  ###                   ~~~ MP LOOP ~~~                    ###
  #############################################################
  
  ## SImple MSE, no option for hcr,or implementation error etc yet, but can/will be incorporated at a later stage
  stock_om_proj    <-   mp_fn(stk = stk_om$stk, idx = idx_oem,                     # stk and idx
                              y0 = y0, iy = iy, dy = dy, fy =fy, nsqy = nsqy,      # info on yrs
                              srbh= sr_om$srbh, srbh.res= sr_om$srbh.res,          # stock recruit
                              assessment= assessment,                              # model to use in the forecasted years
                              Bpa=refpts_om$Bpa, Fmsy=refpts_om$Fmsy,              # reference points
                              seed.nb = seed.nb)              
  
  return(stock_om_proj)
  
  detach(values)
}

