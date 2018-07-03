#' Run the full MSE
#'
#' @param stk object of class FLStock
#' @param idx object of class FLindex
#' @param ny number of years to project from initial year (iy = maxyear of stk + 1)
#' @param it number of iterations
#' @param assessment assessment model to be used, default is sam (will also be set up for spict)
#' @param nsqy number of years to compute status quo metrics, default set at 3
#' @param qmod_init catchability submodel for the OM set up
#' @param fmod_init fishing mortality submodel for the OM set up
#' default ~te(replace(age, age>9,9), year, k=c(6,8))
#' @param mcmc_init number of mcmc for the OM set up
#' @param sr_init stock-recruit model for the OM set up, default is bevholt
#'
#' @return mse runs
#' @export
#' @examples 
#' \dontrun{
#' 
#' ## Run ple4 example - need to adjust for different number of ages between stk and idx
#' data(ple4)
#' data(ple4.index)
#' stk   = ple4[1:8,]
#' idx   = FLIndices(idx=ple4.index)
#' test_ple  <-  mse_base(stk = ple4, idx = FLIndices(idx=ple4.index), it = 2, ny = 3)
#' 
#' ## Run cod example
#' data("cod4_stk")
#' data("cod4_idx")
#' test_cod  <-  mse_base(stk = cod4_stk, idx = cod4_idx, it = 2, ny = 3)
#' }


mse_base <- function(stk, idx, it, ny, nsqy = 3, 
                     sr_init="bevholt", fmod_init=NULL, qmod_init=NULL, mcmc_init =100, 
                     assessment = "sam") {
  
  
  ######################################
  ### ~~~ Set up the simulations ~~~ ###
  ######################################
  
  y0 <- range(stk)["minyear"] # initial data year
  dy <- range(stk)["maxyear"] # final data year
  iy <- dy+1  # initial year of projection (also intermediate year)
  fy <- dy+ny # final year of projection
  
  
  #############################
  ### ~~~ Set up the OM ~~~ ###
  #############################
  
  ## This is to create an FLStock from a real stock
  stk_om       <-   create_FLStock(stk =  stk, idx = idx, it = it, ny = ny, fmod = fmod_init, qmod = qmod_init, mcmc = mcmc_init) 
  
  
  ##################################################
  ### ~~~ Set up the Stock-Recruitment model ~~~ ###
  ##################################################
  
  # Fit stock-recruit model
  # A Beverton-Holt stock-recruit model is fitted for each iteration, with residuals
  # generated for the projection window based on the residuals from the historic period.
  # A stock-recruit model is also fitted to the "median" stk for reference points.
  sr_om        <-   om_sr_model(stk_om$stk, sr_model=sr_init, it =  it, iy = iy, fy = fy)
  
  
  ##################################################
  ###   ~~~ Calculate the reference points ~~~   ###
  ##################################################
  
  # Calculate reference points and set up the operating model for the projection window
  # Reference points based on the "median" stk, assuming (for illustrative purposes only)
  # that Bpa=0.5Bmsy and Blim=Bpa/1.4. The stf method is applied to the operating model stk object
  # in order to have the necessary data (mean weights, etc.) for the projection window.
  refpts_om    <- ref_pts(stk_om$stk0, sr_om$srbh0)
  
  
  ############################################
  ### ~~~ Set up the Observation model ~~~ ###
  ############################################
  
  ## Prepare the FLStock object for projections
  stk_om$stk   <- stf(stk_om$stk, fy-dy, nsqy, nsqy)
  

  ## This is to create an FLindices from a real stock
  idx_info     <-   create_FLIndices(idx =  idx, stk = stk_om$stk, stk0= stk_om$stk0, it = it)
  
  
###################################
### ~~~ Set up the MSE loop ~~~ ###
###################################
  
  final_mse    <-   mse_fn(stk = stk_om$stk, idx = idx_info,                    # stk and idx
                           dy = dy, ny = ny,                                    # info on yrs
                           srbh= sr_om$srbh, srbh.res= sr_om$srbh.res,          # stock recruit
                           assessment= assessment,                              # model to use in the forecasted years
                           Bpa=refpts_om$Bpa, Fmsy=refpts_om$Fmsy)              # reference points
  
  
}

