#' Simple list of reference points
#' 
#' @param stk FLStock (median value if several iterations)
#' @param srbh FLSR (median value if several iterations)
#'
#' @description needs to be tidied up for options
#' 
#' @return list of reference points
#' @export


ref_pts <- function(stk, srbh) {
  
  require(FLRP)
  
  # Calculate the reference points
  brp  <- FLRP::brp(FLBRP(stk, srbh))
  Fmsy <- c(FLRP::refpts(brp)["msy","harvest"])
  msy  <- c(FLRP::refpts(brp)["msy","yield"])
  Bmsy <- c(FLRP::refpts(brp)["msy","ssb"])
  Bpa  <- 0.5*Bmsy
  Blim <- Bpa/1.4

  output_list        <- list(brp,Fmsy,msy,Bmsy,Bpa,Blim)
  names(output_list) <- c("brp","Fmsy","msy","Bmsy","Bpa","Blim")

  return(output_list)
  
}


