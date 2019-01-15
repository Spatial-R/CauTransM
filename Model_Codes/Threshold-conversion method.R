####################################################################################################
####################################################################################################
####################################################################################################

library(pomp);require(magrittr)

source("Model-Codes/CreateCovars.R");
source("Model-Codes/CreateMeasures.R");
source("Model-Codes/CreatEstimation.R");
source("Model-Codes/CreatProcess.R");
source("Model-Codes/CreateModels.R");

region <- "杭州"

city.data <- read.csv("Data/city.csv",header = T,stringsAsFactors = F);
site.data <- read.csv("Data/Site_Data.csv",header = T,stringsAsFactors = F)[,c(2:5)];
names(site.data) <- c("site","name","lat","lon");
site.num <- site.data[site.data$name == region,1]   
city.name <- city.data[city.data$name == region,1]

case.data <- read.csv("Data/Cases/dattran2.csv",header = T,stringsAsFactors = F) %>%
  mutate(date = as.Date(date),week = floor_date(date, unit = "week")) %>%
  group_by(week,city) %>%
  summarize(cases = n()) %>%
  filter(city == city.name & (substr(week,1,4) %in% c(2012:2014))) %>%
  select(week,cases) %>%
  data.frame() %>%
  mutate(week = (julian(week,origin = as.Date("2012-01-01")))/365.25 + 2012) 


tem.air <- climatic.data(dir = "Data/Meterological",pattern.type = "txt",site.code = site.num)
names(tem.air) <- c("MT","RH","CR","date"); tem.air <- arrange(tem.air,date)  
tem.air[,1:3] <- apply(tem.air[,1:3],2,na.approx)     ####### interplation for the missing data
tem.air$AH <- AH(WVP(tem.air$RH,SVP(C2K(tem.air$MT/10))),C2K(tem.air$MT/10))*1000 ## Absolute humidity


create_model_threshold <- function(covar.file.name = tem.air,
                                   case.file.name = case.data,
                                   foi.name,
                                   foi.threshold.name = "RH1",
                                   foi.threshold.value = c(seq(70,140,10)),
                                   foi.threshold.type = "V",
                                   cr.name,
                                   rr.name,
                                   rr.threshold.name = "NO",
                                   rr.threshold.value = c(9,20),
                                   rr.threshold.type = "U",
                                   foi.continue = T,
                                   rr.continue = F){


origin.data <- create_covars(covar.file.name = covar.file.name,
                             case.file.name = case.file.name,
                             foi.threshold.name = foi.threshold.name,
                             foi.threshold.value = foi.threshold.value,
                             foi.threshold.type = foi.threshold.type,
                             rr.threshold.name = rr.threshold.name,
                             rr.threshold.value = rr.threshold.value,
                             rr.threshold.type = rr.threshold.type,
                             foi.continue = foi.continue,
                             rr.continue = rr.continue)

case.num <-  which("cases" == names(origin.data)); week.num <-  which("week" == names(origin.data))
mumps.hz <-  origin.data[,c(week.num,case.num)]; covar <- origin.data[,-c(case.num)]

target.parm <- paste(foi.threshold.name,"ths",sep = "")
#first.cond <- stri_detect_fixed(names(covar),target.parm);
#second.cond <- stri_detect_regex(names(covar),pattern = "[0-9]$") 
pattern <- paste("^(",target.parm,")[0-9a-z_A-Z]+[0-9]$",sep = "")
target.column <- which(stri_detect_regex(names(covar),pattern = pattern))
#target.column <- which(first.cond & second.cond)
result.list <- list()

for (i in c(1:length(target.column))){
  covar[,target.parm] <- covar[,target.column[i]] 
  model.sructure <- create_model_basic(case.file.name = mumps.hz,
                                        covar.file.name = covar,
                                        foi.name = foi.name,
                                        contact.name = cr.name,
                                        report.name = rr.name,
                                        full = T);
  threshold.name <- gsub("ths","",target.parm);
  threshold.value <- as.numeric(gsub(target.parm,"",names(covar)[target.column[i]]));
  
  result <- list(model = model.sructure[[1]],
                 threshold.name = threshold.name, 
                 threshold.value = threshold.value,
                 parameter.name = model.sructure[[2]])
  result.list[[i]] <- result
   print(paste("The model for ",dQuote(gsub("ths","",target.parm)), " with the threshold value at ",
           threshold.value," was successfully created.",sep = ""))
}
return(result.list)
}

model.threshold <- create_model_threshold(covar.file.name = tem.air,
                       case.file.name = case.data,
                       foi.name = c("MT1","AH","CR1"),
                       foi.threshold.name = "MT1",
                       foi.threshold.value = c(seq(70,140,10)),
                       foi.threshold.type = "V",
                       cr.name = c("MT2","RH2","CR2"),
                       rr.name = c("MT3","RH3","CR3"),
                       rr.threshold.name = "NO",
                       rr.threshold.value = c(9,20),
                       rr.threshold.type = "U",
                       foi.continue = T,
                       rr.continue = F)