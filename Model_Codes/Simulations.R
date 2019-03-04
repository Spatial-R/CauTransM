#######################################################################################
#################################### Simulation in MIF ################################
#######################################################################################

library(pomp);require(magrittr);
library(doSNOW)

region <- "杭州"

rw.sd_rp <- 0.02;rw.sd_ivp <- 0.02;

######################################### load the functions #########################

source("Model-Codes/CreateCovars.R");
source("Model-Codes/CreateMeasures.R");
source("Model-Codes/CreateEstimation.R");
source("Model-Codes/CreateProcess.R");
source("Model-Codes/CreateModels.R");
source("Model-Codes/CreateDataset.R");

####################################### create the covariates  ######################## 

origin.data <- create_covars(covar.file.name = tem.air,
                             case.file.name = case.data,
                             foi.threshold.name = "MT1",
                             foi.threshold.value = 40,
                             foi.threshold.type = "V",
                             rr.threshold.name = "MT3",
                             rr.threshold.value = c(9,20),
                             rr.threshold.type = "U",
                             foi.continue = F)
case.num <-  which("cases" == names(origin.data)); week.num <-  which("week" == names(origin.data))
mumps.hz <-  origin.data[,c(week.num,case.num)]; covar <- origin.data[,-c(case.num)]


theta.lo <- c(mu = 5.92/1000,sigma=365/25,gamma=365/10,alpha = 3.115105e-01,iota = 5.952648e+00,rho = 0.2,
              R0 = 3,sigmaSE = 0,psi = 0.1,amplitude = 0, S_0 = 0.02,I_0 = 2e-06, E_0 = 2e-06)
theta.hi <- c(mu = 5.92/1000,sigma=365/12,gamma=365/4,alpha = 1,iota = 40,rho = 0.7,R0 = 10,
              sigmaSE = 1,psi = 0.4,amplitude = 1, S_0 = 0.03,I_0 = 0.0001,E_0 = 0.0001)
guess.parameter <- cbind(theta.lo,theta.hi)


create_global_simulation(case.file.name = mumps.hz,
                                     covar.file.name = covar,
                                     origin.parameter.value = guess.parameter,
                                     focal.value.range = c(0.5,5.0),
                                     rw.sd_rp = 0.02,
                                     rw.sd_ivp = 0.03,
                                     Np = 3000,
                                     Nmif = 50,
                                     mcors = 2,
                                     random.number = 5000)
  

  
  