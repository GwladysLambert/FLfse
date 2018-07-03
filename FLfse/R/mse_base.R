#' Run the full MSE
#'
#' @param stk object of class FLStock
#' @param idx object of class FLindex
#' @param ny number of years to project from initial year (iy = maxyear of stk + 1)
#' @param it 
#' @param assessment assessment model to be used, default is sam (will also be set up for spict)
#'
#' @return mse runs
#' @export
#' @examples 
#' \dontrun{
#' stk = ple4
#' idx = FLIndices(idx=ple4.index)
#' }


mse_base <- function(stk, idx, it, ny, assessment = "sam") {
  
  
#############################
### ~~~ Set up the OM ~~~ ###
#############################
  
  ## This is to create an FLStock from a real stock and get some reference points
  stk_info     <-   create_FLStock(stk =  stk, idx = idx, it = it, ny = ny) 
  

############################################
### ~~~ Set up the Observation model ~~~ ###
############################################

  ## This is to create an FLindices from a real stock
  idx_info     <-   create_FLIndices(idx =  idx, stk = stk_info$stk, stk0= stk_info$stk0, it = it)
  
  
###################################
### ~~~ Set up the MSE loop ~~~ ###
###################################
  
  final_mse    <-   mse_fn(stk = stk_info$stk, idx = idx_info,                # stk and idx
                           dy=stk_info$dy, ny = ny,                           # info on yrs
                           srbh= stk_info$srbh, srbh.res= stk_info$srbh.res,  # stock recruit
                           assessment= "sam",                                 # model to use in the forecasted years
                           Bpa=stk_info$Bpa, Fmsy=stk_info$Fmsy)              # reference points
  
  
}

