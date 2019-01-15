########################################################################################
#################################### Global Searching ##################################
########################################################################################

create_global_simulation <- function(case.file.name = mumps.hz,
                                     covar.file.name = covar,
                                     foi.name = c("MT1","AH"),
                                     contact.name = c("MT2","RH2"),
                                     report.name = c("MT3","RH3"),
                                     full = T,
                                     origin.parameter.value,
                                     focal.value.range = c(0.5,2.0),
                                     rw.sd_rp = 0.02,
                                     rw.sd_ivp = 0.03,
                                     Np = 3000,
                                     Nmif = 50,
                                     mcors = 2,
                                     random.number = 5000){
  
  basic_model <- create_model_basic(case.file.name = case.file.name,
                                    covar.file.name = covar.file.name,
                                    foi.name = foi.name,
                                    contact.name = contact.name,
                                    report.name = report.name,
                                    full = full)
  
  ####################################### initial parameter values ###########################
  
  state.parameter <- c("S_0","E_0","I_0","R_0")
  model.parameter <- setdiff(sort(basic_model$model.parameter.name),state.parameter);
  focal.parameter <- basic_model$focal.parameter.name
  origin.model    <- basic_model$model
  
  match.v <- identical(sort(row.names(origin.parameter.value)),
                       sort(c(model.parameter,state.parameter[1:3])))
  
  if(!isTRUE(match.v)){
    print("Parameter names in pomp model doesn't match the parameter in previous settings")
  }
  
  focal.value.low <- rep(focal.value.range[1],length(focal.parameter));
  focal.value.hih <- rep(focal.value.range[2],length(focal.parameter));
  focal.value <- cbind(theta.lo = focal.value.low, theta.hi = focal.value.hih)
  row.names(focal.value) <- focal.parameter
  initia.value <- rbind(origin.parameter.value,focal.value)
  
  fixed.parameter  <- c("gamma","mu","R0","sigma");   #### Those parameters will not vary in mif2
  #### Those parameters will vary in mif2
  varied.parameter <- c(setdiff(model.parameter,fixed.parameter),focal.parameter) 
  parameter.rw.sd  <- c(rep(rw.sd_rp,length(varied.parameter)),rep(rw.sd_ivp,length(state.parameter)))
  ### (length(names(parameter.rw.sd)) - 3): length(names(parameter.rw.sd))

  names(parameter.rw.sd) <- c(varied.parameter,state.parameter);
  parameter.rw.sd  <- do.call(rw.sd, as.list(parameter.rw.sd))
  
  cl <- makeCluster(mcors,type = "SOCK"); registerDoSNOW(cl)
  set.seed(998468235L,kind = "L'Ecuyer")
  mcor.set <- list(preschedule = FALSE,set.seed = TRUE)
  
  global.result <- foreach(i = 1:random.number, #.combine=rbind,
                           .packages = c("pomp","magrittr"),.errorhandling = "remove",
                           .inorder = FALSE, .options.multicore = mcor.set) %dopar%  {
                             start.parameter <- apply(initia.value,1,function(x)runif(1,x[1],x[2]))
                             start.parameter[["R_0"]] <- 1 - start.parameter[["S_0"]] - start.parameter[["E_0"]] - start.parameter[["I_0"]]
                             mif2(origin.model, start = start.parameter, Np = Np,Nmif = Nmif,
                                  cooling.type = "geometric",cooling.fraction.50 = 0.1, 
                                  transform = TRUE,rw.sd = parameter.rw.sd) %>%  mif2() -> mf
                             pf <- replicate(10, pfilter(mf, Np = Np))
                             ll <- sapply(pf,logLik)
                             ll <- logmeanexp(ll, se = TRUE)
                             nfail <- sapply(pf,getElement,"nfail")
                             data.frame(as.list(coef(mf)),
                                        loglik = ll[1],
                                        loglik.se = ll[2],
                                        nfail.min = min(nfail),
                                        nfail.max = max(nfail))
                           }
  global.result <- do.call("rbind",global.result)
  return(global.result)
  stopCluster(cl)
  registerDoSEQ()
}

########################################################################################
#################################### Global Searching ##################################
########################################################################################

create_profile_simulation <- function(case.file.name,
                                      covar.file.name,
                                      profile.name,
                                      profile.value.range,
                                      origin.parameter.value,
                                      rw.sd_rp = 0.02,
                                      rw.sd_ivp = 0.03,
                                      Np = 3000,
                                      Nmif = 50,
                                      mcors,
                                      profile.length =15,
                                      left.length = 20){
  
  basic_model <- create_model_basic( case.file.name = case.file.name,
                                    covar.file.name = covar.file.name)
  
  ####################################### initial parameter values ###########################
  
  state.parameter <- c("S_0","E_0","I_0","R_0")
  model.parameter <- setdiff(sort(basic_model$model.parameter.name),state.parameter);
  focal.parameter <- basic_model$focal.parameter.name
  origin.model    <- basic_model$model
  
  match.v <- identical(sort(row.names(origin.parameter.value)),
                       sort(c(model.parameter,state.parameter)))
  
  if(!isTRUE(match.v)){
    print("Parameter names in pomp model doesn't match the parameter in initial settings")
  }

  column.num <- which(profile.name == names(origin.parameter.value))
  
  if(length(column.num) == 0){
    print(paste("The target parameter",dQuote(profile.name), 
                "was not in the", dQuote("origin.parameter.value"),sep = ""))
  }
  
  if(length(profile.value.range) == 2){
    target.lo <- profile.value.range[1];
    target.hi <- profile.value.range[2];
  } else{
    target.lo <- origin.parameter.value[1,column.num];
    target.hi <- origin.parameter.value[2,column.num];
  }
  
  left.parameter.value <- origin.parameter.value[,-c(column.num)]
  left.lo <- unlist(left.parameter.value[1,]); 
  left.hi <- unlist(left.parameter.value[2,]);
  
  fixed.parameter <- c("gamma","mu","R0","sigma"); 
  varied.parameter <- c(setdiff(model.parameter,fixed.parameter),focal.parameter)
  parameter.rw.sd <- c(rep(rw.sd_rp,length(varied.parameter)),rep(rw.sd_ivp,length(state.parameter)))
  names(parameter.rw.sd) <- c(varied.parameter,state.parameter);
  parameter.rw.sd[which(profile.name == names(parameter.rw.sd))] <- 0
  parameter.rw.sd <- do.call(rw.sd, as.list(parameter.rw.sd))
  
  cl <- makeCluster(mcors,type = "SOCK"); registerDoSNOW(cl)
  set.seed(998468235L,kind = "L'Ecuyer")
  mcor.set <- list(preschedule = FALSE,set.seed = TRUE)
  
  startparameter <- profileDesign(targetparameter = seq(from = target.lo,to = target.hi,length = profile.length),
                                  lower = left.lo, upper = left.hi, nprof = left.length)
  startparameter[["R_0"]] <- 1 - startparameter[["S_0"]] - startparameter[["E_0"]] - startparameter[["I_0"]]
  names(startparameter) <- c(profile.name, names(startparameter)[-1])
 
  global.result.local <- foreach(p = iter(startparameter,"row"), .combine = rbind,
                                 .packages = c("pomp","magrittr"),.errorhandling = "remove",
                                 .inorder = FALSE, .options.multicore = mcopts) %dopar% {
                                   mif2(origin.model, start = unlist(p), Np = Np, Nmif = Nmif,
                                        cooling.type = "geometric",cooling.fraction.50 = 0.1, transform = TRUE,
                                        rw.sd = parameter.rw.sd) %>% mif2() -> mf    
                                   pf <- replicate(10, pfilter(mf, Np = Np))
                                   ll <- sapply(pf,logLik)
                                   ll <- logmeanexp(ll,se = TRUE)
                                   nfail <- sapply(pf,getElement,"nfail")
                                   data.frame(as.list(coef(mf)),
                                              loglik = ll[1],
                                              loglik.se = ll[2],
                                              nfail.min = min(nfail),
                                              nfail.max = max(nfail))
                                 }
  return(global.result.local)
  stopCluster(cl)
  registerDoSEQ()
}
