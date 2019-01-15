
dmeas.basic <- "  double m = rhonew * C;
                  double v = m*(1.0-rhonew+psi*psi*m);
                  double tol = 1.0e-18;
                  if (cases > 0.0) {
                  lik = pnorm(cases+0.5,m,sqrt(v)+tol,1,0)-pnorm(cases-0.5,m,sqrt(v)+tol,1,0)+tol;
                  } else {
                  lik = pnorm(cases+0.5,m,sqrt(v)+tol,1,0)+tol;
                  }"

rmeas.basic <- "double m = rhonew*C;
                double v = m*(1.0-rhonew+psi*psi*m);
                double tol = 1.0e-18;
                cases = rnorm(m,sqrt(v)+tol);
                if (cases > 0.0) {
                cases = nearbyint(cases);
                } else {
                cases = 0.0;
                }"

create_measure <- function(measure.name){
  
  measure.v <- paste("p",measure.name,sep = "")
  measure.part <- unlist(lapply(1:length(measure.v),function(id){
    parm.tem <- paste("log(",measure.v[id],")","*",measure.v[id],sep = "")
    return(parm.tem)
  }))
  measure.fin.1 <- paste(measure.part,collapse = "+")    ##### the formulation for the force of infection
  measure.fin.2 <- paste("double rhonew = exp(",measure.fin.1,")*rho/(exp(",measure.fin.1,")*rho+1);",sep = "")
  
  dmeas <-  Csnippet(paste(measure.fin.2,"\n",dmeas.basic,sep = ""))
  rmeas <-  Csnippet(paste(measure.fin.2,"\n",rmeas.basic,sep = ""))
 
   return(c(dmeas,rmeas))  
}