#' Simple list of reference points
#' 
#' @param stk FLStock (provide median value if several iterations)
#' @param srbh FLSR (provide median value if several iterations)
#' @param brps FLBRP if already existing (no need for other \code{stk}
#' or \code{srbh} in that case)
#'
#' @description needs to be tidied up for options
#' 
#' @return list of reference points
#' @export


ref_pts <- function(stk=NULL, srbh=NULL, brps=NULL) {
  
  require(FLRP)
  
  # Calculate the reference points
  if (is.null(brps))  {
    brp  <- FLRP::brp(FLBRP(stk, srbh)) } else { brp <- brps }
  Fmsy <- c(FLRP::refpts(brp)["msy","harvest"])
  msy  <- c(FLRP::refpts(brp)["msy","yield"])
  Bmsy <- c(FLRP::refpts(brp)["msy","ssb"])
  Bpa  <- 0.5*Bmsy
  Blim <- Bpa/1.4

  output_list        <- list(brp,Fmsy,msy,Bmsy,Bpa,Blim)
  names(output_list) <- c("brp","Fmsy","msy","Bmsy","Bpa","Blim")

  return(output_list)
  
}


