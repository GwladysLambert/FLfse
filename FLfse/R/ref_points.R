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
  
  require(FLBRP)
  
  # Calculate the reference points
  brp  <- FLBRP::brp(FLBRP(stk, srbh))
  Fmsy <- c(FLBRP::refpts(brp)["msy","harvest"])
  msy  <- c(FLBRP::refpts(brp)["msy","yield"])
  Bmsy <- c(FLBRP::refpts(brp)["msy","ssb"])
  Bpa  <- 0.5*Bmsy
  Blim <- Bpa/1.4

  output_list <- list(brp,Fmsy,msy,Bmsy,Bpa,Blim)
  names(output_list) <- c("brp","Fmsy","msy","Bmsy","Bpa","Blim")

  return(output_list)
  
}


