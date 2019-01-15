
create_model_basic <- function(case.file.name = case.data,
                               covar.file.name = tem.air,
                               foi.name = c("MT1","AH"),
                               contact.name = c("RH2"),
                               full = T,  #### two arguements for process
                               report.name = c("MT3","RH3lag")){

  #######################################################################################
  ################################# Gather pomp structure ###############################
  #######################################################################################
  
  case.data <- case.file.name;  covar.data <- covar.file.name
  
  foi.param     <- paste("p",foi.name,sep = "");
  contact.param <- paste("p",contact.name,sep = "");
  report.param  <- paste("p",report.name,sep = ""); 
  focal.param   <- c(foi.param,contact.param,report.param)
  
  rproc    <- creat_process(foi.name = foi.name, 
                            contact.name = contact.name, 
                            full = full)
  initlz   <- create_estimation(param.name = focal.param)[[1]];
  fromEst  <- create_estimation(param.name = focal.param)[[2]];
  toEst    <- create_estimation(param.name = focal.param)[[3]];
  dmeas    <- create_measure(measure.name = report.name)[[1]];
  rmeas    <- create_measure(measure.name = report.name)[[2]]
  
  model.param <- c("R0","mu","sigma","gamma","alpha","iota",
                      "rho","sigmaSE","psi","cohort","amplitude",
                      "S_0","E_0","I_0","R_0")
  
  total.parameters <- c(model.param,focal.param)
  
  case.data %>%
    pomp(t0 = case.data$week[1],
         time = "week",
         rprocess = euler.sim(step.fun=rproc,delta.t=2/365.25),
         initializer = initlz,
         dmeasure = dmeas,
         rmeasure = rmeas,
         covar = covar.data,
         toEstimationScale = toEst,
         fromEstimationScale = fromEst ,
         tcovar = "week",
         zeronames = c("C","W"),
         statenames = c("S","E","I","R","C","W"),
         paramnames = total.parameters
    ) ->  origin.model
  
  result <- list(origin.model,model.param,focal.param); 
  names(result) <- c("model","model.parameter.name","focal.parameter.name")
  return(result)
}

#model <- create_model_basic()
