# Skeleton
scenario_skeleton <- list(scns    = data.frame(id = "default", stringsAsFactors = FALSE),
                      stock     = list(stk = NULL,  
                                    latin.name=NULL, common.name=NULL, region=NULL, 
                                    params = NULL, h=0.75),
                      ctrl.om   = list(sr_init="bevholt", fmod_init=NULL, qmod_init=NULL, mcmc_init =100),
                      ctrl.idx  = list(idx = NULL, qmod_idx = NULL, qmod_pars_idx = NULL, add.error_idx = FALSE),
                      ctrl.sims = list(it=3, ny=3, nsqy = 3, seed.nb= 321), 
                      ctrl.mp   = list(assessment="sam")
)
