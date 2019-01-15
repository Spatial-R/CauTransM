

toEst.basic <- "     Talpha = log(alpha);
                     Tiota = log(40-iota);
                     Trho = logit(rho);
                     Tcohort = logit(cohort);
                     Tamplitude = logit(amplitude);
                     TsigmaSE = log(sigmaSE);
                     Tpsi = log(psi);
                     to_log_barycentric (&TS_0, &S_0, 4);"


fromEst.basic <- "    Talpha = exp(alpha);
                      Tiota = 40-exp(iota);
                      Trho = expit(rho);
                      Tcohort = expit(cohort);
                      Tamplitude = expit(amplitude);
                      TsigmaSE = exp(sigmaSE);
                      Tpsi = exp(psi);
                      from_log_barycentric (&TS_0, &S_0, 4);"


create_estimation <- function(param.name,full = T){
 
  create.log.name <- unlist(lapply(1:length(param.name),function(id){
    log.tem <- paste("T",param.name[id],"=log(",param.name[id],");",sep = "")
    return(log.tem)
  }))
  create.log.name <- paste(create.log.name,collapse = "\n ");#cat(create.log.name)
  
  
  create.exp.name <- unlist(lapply(1:length(param.name),function(id){
    log.tem <- paste("T",param.name[id],"=exp(",param.name[id],");",sep = "")
    return(log.tem)
  }))
  create.exp.name <- paste(create.exp.name,collapse = "\n ");#cat(create.exp.name)

  toEst   <- Csnippet(paste(toEst.basic,"\n",create.log.name,sep = ""))
  fromEst <- Csnippet(paste(fromEst.basic,"\n",create.exp.name,sep = ""))
  
  initlz <- Csnippet("
                      double m = pop/(S_0+E_0+I_0+R_0);
                      S = nearbyint(m*S_0);
                      E = nearbyint(m*E_0);
                      I = nearbyint(m*I_0);
                      R = nearbyint(m*R_0);
                      W = 0;
                      C = 0;
                      ")
  return(c(initlz,fromEst,toEst))
}
