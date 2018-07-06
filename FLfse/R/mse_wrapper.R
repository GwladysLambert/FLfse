#' Master function to run MSE simulations
#'
#' This is the main high-level wrapper function for running \pkg{FLfse}
#' simulations. 
#'
#' @param scenarios list of scenarios to run
#' @param parallel if TRUE, will run scenarios in parallel. Default is F
#' 
#' @return output of mse_base
#' 
#' @export
#' 
#' @examples
#' \dontrun{
#' simple_scenario1 <- list(scens = list(id = "C.hareng"), 
#'                          stock = list(latin.name="Clupea harengus"),
#'                          ctrl.sims = list(it = 2, ny = 3))
#' simple_scenario2 <- list(scens = list(id = "Linf=100"),
#'                          stock = list(biol.params=data.frame(linf=100)),
#'                          ctrl.sims = list(it = 2, ny = 3))
#' 
#' test1 <- run_mse(scenarios=list(simple_scenario1, simple_scenario2), parallel = F)
#' test2 <- run_mse(scenarios=list(simple_scenario1, simple_scenario2), parallel = T)
#' }
#' 


run_mse <- function(scenarios, parallel = FALSE) {
  
 # browser()
  # Get arguments for each scenario
  arg_list <- lapply(scenarios, function(scenario) {
    list(
      scens             = scenario$scens,
      stock             = scenario$stock,
      ctrl.om           = scenario$ctrl.om,
      ctrl.idx          = scenario$ctrl.idx,
      ctrl.sims         = scenario$ctrl.sims,
      ctrl.mp           = scenario$ctrl.mp)
  })
  
  # to satisfy R CMD check in the foreach() call below
  x <- NULL
  
  if (parallel) {
    message("Running scenarios in parallel.")
    output <- foreach::foreach(x = arg_list, .packages = "FLfse", .combine = list,
                               .init = NULL,
                               .verbose = TRUE) %dopar% #
      {do.call("mse_base", x)
      message(paste("Finished running",x$scens$id))
      }
  } else {
    message("Running scenarios sequentially.")
    output <- lapply(arg_list, function(x) {
      do.call("mse_base", x)
      message(paste("Finished running",x$scens$id))
    })
  }
  return(output)
}

