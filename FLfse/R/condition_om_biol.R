#' Function to create FLStock of the OM for MSE simulation
#' from biological data using existing dataset of WKLife or
#' just providing own parameters
#'
#' @param latin.name vector of latin names of species to be selected from WKLife dataset
#' (if given, no need for common.name)
#' @param common.name vector of common names of species to be selected from WKLife dataset
#' (if given, no need for latin.name)
#' @param region stock area from WKLife dataset
#' @param params dataframe of biological params if not using from WKLife stocks 
#' (i.e. if not using the 3 above arguments). No need to provide all variables as FLife is used to 
#' complete the set
#' @param h steepness default as 0.75 (can ba a vecor of same length as number of species or stocks)
#' @param it number of iterations, default is 5
#' @param seed.nb set seed number
#'
#' @details creates the FLStock, script based on
#' https://github.com/shfischer/wklifeVII/blob/master/R/OM1.R
#' @return a list of "sr", "brp", "stk"
#'
#' @export
#' 
#' @examples 
#' \dontrun{
#' 
#' # Example using given life history parameter(s)
#' test1 <- create_FLStock_biol(params=data.frame(linf=100))
#' stock.n(test1[[1]]$stk)
#' 
#' # Example using WKLife dataset
#' test2 <- create_FLStock_biol(latin.name="Clupea harengus")
#' stock.n(test2[[1]]$stk)
#' }


create_FLStock_biol <- function (latin.name=NULL, common.name=NULL, region=NULL, 
                                 params = NULL, h=0.75, it = 5, seed.nb=321) {

  ############################
  ### ~~~ Prepare data ~~~ ###
  ############################

  ### extended list
  #setwd("C:\\Users\\GL04\\OneDrive - CEFAS\\MA016N\\Package\\FLfse\\R")
  stocks_lh <- read.csv("FLfse/R/input/stock_list_full2.csv")
  names(stocks_lh)[1] <- "latin.name"
  
  ### subset requested stock
  if(!is.null(latin.name))  stock_sub <- stocks_lh[stocks_lh$latin.name %in% latin.name,]
  if(!is.null(common.name)) stock_sub <- stocks_lh[stocks_lh$common %in% common.name,]
  if(!is.null(region))      stock_sub <- stocks_lh[stocks_lh$area %in% region,]
  if(!is.null(params))      {stock_sub <- params; stock_sub$stock <- "fake"[1:nrow(stock_sub)]}
  
  ### set steepness
  stock_sub$s <- h
  
  ###############################
  ### ~~~ Create FLStocks ~~~ ###
  ###############################
  
  OMs <- foreach(i = split(stock_sub, 1:nrow(stock_sub)), .errorhandling = "pass", 
                 .packages = c("FLife", "FLasher")) %dopar% {
                   
                   ## create brp
                   ### get lh params
                   lh_res             <- c(dimnames(lhPar(FLPar(linf = 1)))$params, "l50")
                   lh_avail           <- intersect(lh_res, names(i))
                   lh_pars.tmp        <- data.frame(i[, lh_avail])
                   names(lh_pars.tmp) <- lh_avail
                   lh_pars            <- data.frame(lh_pars.tmp[, !is.na(lh_pars.tmp)])
                   names(lh_pars)     <- names(lh_pars.tmp[, !is.na(lh_pars.tmp)])
                   lh_pars            <- as(lh_pars, "FLPar")
                   ### create missing pars
                   lh_pars            <- lhPar(lh_pars)
                   
                   # Max age: age at l = 0.95 * linf
                   max_age            <- ceiling(log(0.05)/(-c(lh_pars["k"]))+c(lh_pars["t0"]))
                   
                   # if minfbar and maxfbar are not given use age 1 to max age
                   if (is.null(i$minfbar)) i$minfbar <- 1
                   if (is.null(i$maxfbar)) i$maxfbar <- max_age
                   
                   ### create brp
                   brp <- lhEql(lh_pars, range = c(min = 1, max = max_age,
                                                   minfbar = i$minfbar, maxfbar = i$maxfbar,
                                                   plusgroup = max_age))
                   
                   ### save life-history parameters in FLBRP
                   attr(brp, "lhpar") <- lh_pars
                   
                   #### coerce FLBRP into FLStock
                   ### keep only second year (first year with non zero catch)
                   stk <- as(brp, "FLStock")[, 2]
                   ### name first year "1"
                   stk <- qapply(stk, function(x) {
                     dimnames(x)$year <- "1"; return(x)
                   })
                   
                   ### extend object to year 100
                   stk <- fwdWindow(stk, brp, end = 100)
                   
                   ### propagate with requested number of iterations
                   stk <- propagate(stk, it)
                   
                   ### create FLSR object
                   stk_sr <- FLSR(params = params(brp), model = model(brp))
                   
                   ### create residuals for (historical) projection
                   set.seed(seed.nb)
                   residuals <- rlnoise(it, rec(stk) %=% 0, sd = 0.2, b = 0.3)
                   
                   ### project forward, targeting 0.5 F_MSY
                   years_target <- ((dims(stk)$minyear + 1):100) # set to 75 for example if different fishing pattern for last 25 years? !!!!!!!!!!
                   stk <- fwd(stk, sr = stk_sr, 
                              control = fwdControl(year = years_target, 
                                                   value = refpts(brp)['msy', 'harvest']*0.5, ## here we can add the option to be able to set whatever in the scenarios!!!
                                                   quant = "f"),
                              residuals = residuals[, ac(years_target)])
                   
                   ### THIS IS WHERE WE COULD CHANGE THE FISHING PATTERN, e.g. roller coaster or one way !!!!!!!!!!!! - could be in the scenarios again
                   
                   ### return list
                   return(list(sr = stk_sr, brp = brp, stk = stk))
                   
                 }
  
  ### name the FLStock objects
  OMs <- lapply(seq_along(OMs), function(x){
    name(OMs[[x]]$stk) <- ac(stock_sub$stock[x])
    return(OMs[[x]])
  })
  
  ### set names for list elements
  names(OMs) <- stock_sub$stock

  ### subset to last 26 years
  OMs <- lapply(OMs, function(x){
    names(x)[length(names(x))] <- "stk"
     x$stk <- window(x$stk, start = 75)
     return(x)
   })
  
 # stopCluster(cl)

  return(OMs)
    
}

