
##################################################################################################
######################### creat different types of model to explore    ###########################
####################### the effect of the climatic factors on Mumps dynamic ######################
##################################################################################################


library(readr);library(dplyr);library(lubridate);
library(reshape2);require(magrittr);library(zoo);library(stringi)
library(humidity)


climatic.data <- function(dir,pattern.type,site.code  = "58633"){
  
  #### @ dir: which directory has the dataset for climatic dataset
  #### @ pattern.type: which kind of dataset is: .csv,.txt or other formulation 
  #### @ city.code: unique code for the unique city name. i.e., 32700 is the code for hangzhou city 
  
  if(!dir.exists(dir)){
    print("The directory doesn't exist, please select the true directory")
  } else{
    pattern.type <- paste0(".",pattern.type,sep = "") 
    all.files <- list.files(dir, pattern = pattern.type,full.names = T)
    if (pattern.type == "csv") {
      mylist <- lapply(all.files,function(i) read.csv(i,header = T,stringsAsFactors = F))
    } else {
      mylist <- lapply(all.files,function(i) read.table(i,header = T,stringsAsFactors = F))
    } 
    dat.tem <- do.call('rbind',mylist)
    dat.tem$V13201 <- ifelse(dat.tem$V13201 == "32700",0,dat.tem$V13201) #### replace the missing with 0
    dat.tem$date <- paste(dat.tem$V04001,dat.tem$V04002,dat.tem$V04003,sep = "-")
    dat.tem$date <- as.Date(dat.tem$date)
    dat.tem <- dat.tem[dat.tem$V01000 == site.code,c("V01000","V12001","V13003","V13201","date")]
    return(dat.tem[,-1])
  }  
}


threshold.set <- function(data,v.name,t.value,continue = F,tModel = "V"){
  ## @  data : the dataset to calculate the value
  ## @  variable: which variable to set the threshold
  ## @  threshold.value
  ## @  continue: construate a continue threshold settings
  ## @  tModel: V or U
  
  if(isTRUE(continue) & (length(t.value) == 1)){
    print("The parameter t.value must be a vector with length two")
  }
  if((tModel == "U") & (length(t.value) == 1)){
    print("The parameter t.value must be a vector with length two")
  }
  if(isTRUE(continue) & (tModel == "U")){
    print("The paramter value of model.type must be V if the continue is true")
  }
  dat.tem <- data
  
  if(isTRUE(continue) & (tModel == "V")){
    for (i in t.value){
      t.name <- paste(v.name,"ths",i,sep = "")
      dat.tem[,t.name] <- unlist(lapply(1:nrow(dat.tem),function(id){
        if (dat.tem[id,v.name] >= i) {
          threshold.value <- dat.tem[id,v.name] - i
        } else {
          threshold.value <- i - dat.tem[id,v.name]
        }
          return(threshold.value)
      }))
    }
  } else if(isTRUE(continue) & (tModel == "U")){
    stop("The paramter value of model.type must be V if the continue is true")
  } else if(!isTRUE(continue) & (tModel == "V")){
    t.name <- paste(v.name,"ths",sep = "")
    dat.tem[,t.name] <- unlist(lapply(1:nrow(dat.tem),function(id){
      if (dat.tem[id,v.name] >= t.value) {
        threshold.value <- dat.tem[id,v.name] - t.value
      } else {
        threshold.value <- t.value - dat.tem[id,v.name]
      }}))
  } else {
    t.name <- paste(v.name,"ths",sep = "")
    dat.tem[,t.name] <- unlist(lapply(1:nrow(dat.tem),function(id){
      if (dat.tem[id,v.name] >= t.value[2]) {
        threshold.value <- dat.tem[id,v.name] - t.value[2]
      } else if(dat.tem[id,v.name] <= t.value[1]){
        threshold.value <- t.value[1] - dat.tem[id,v.name]
      } else{
        threshold.value = 0
      }
      return(threshold.value)
    }))
  }
  return(dat.tem)  
}

create_covars <- function(covar.file.name,
                          case.file.name,
                          foi.threshold.name = "MT1",
                          foi.threshold.value = 4,
                          foi.threshold.type = "V",
                          foi.continue = FALSE,
                          rr.threshold.name = "MT3",
                          rr.threshold.value = c(12,25),
                          rr.threshold.type = "U",
                          rr.continue = FALSE){
  
  data.air <- covar.file.name; data.case.name <- case.file.name
  data.air[,"MT1"] <- data.air[,"MT2"] <- data.air[,"MT3"] <- data.air[,"MT"];
  data.air[,"RH1"] <- data.air[,"RH2"] <- data.air[,"RH3"] <- data.air[,"RH"];
  data.air[,"CR1"] <- data.air[,"CR2"] <- data.air[,"CR3"] <- data.air[,"CR"];
  
  ########################## set the threhold ##############################

  if(foi.threshold.name == "MT1"){
       data.air <- threshold.set(data = data.air, v.name = "MT1", 
                             t.value = foi.threshold.value,
                             continue = foi.continue,
                             tModel = foi.threshold.type)  
  } else if(foi.threshold.name == "AH"){
       data.air <- threshold.set(data = data.air, v.name = "AH", 
                             t.value = foi.threshold.value,
                             continue = foi.continue,
                             tModel = foi.threshold.type)  
  } else if(foi.threshold.name == "CR") {
    data.air <- threshold.set(data = data.air, v.name = "CR1", 
                              t.value = foi.threshold.value,
                              continue = foi.continue,
                              tModel = foi.threshold.type)  
  } else {
    data.air <- data.air
  }
    
  
if(rr.threshold.name == "MT3"){
    data.air <- threshold.set(data = data.air, v.name = "MT3", 
                             t.value = rr.threshold.value,
                             continue = rr.continue,
                             tModel = rr.threshold.type)  
} else if(rr.threshold.name == "RH3"){
  data.air <- threshold.set(data = data.air, v.name = "RH3", 
                           t.value = rr.threshold.value,
                           continue = rr.continue,
                           tModel = rr.threshold.type)  
  } else if(rr.threshold.name == "CR3"){
    data.air <- threshold.set(data = data.air, v.name = "CR3", 
                              t.value = rr.threshold.value,
                              continue = rr.continue,
                              tModel = rr.threshold.type)  
  }else{
    data.air <- data.air
  }
  
  
  ############################# lag effect #######################  
  
  lag.name <- setdiff(names(data.air),c("date","MT"))
  
  for (j in lag.name){
    lag.v <- paste(j,"lag",sep = "")
    data.air[,lag.v] <- c(data.air[-c(1:7),j],rep(NA,7)) 
  }
  
  #################################################################################
  
  data.air <- mutate(data.air,week = floor_date(date, unit = "week"))
  data.air <- aggregate(data.air,by = list(data.air$week),mean)[,-1]
  data.air <- filter(data.air,substr(week,1,4) %in% c("2012","2013","2014"))
  data.air <- mutate(data.air,week = (julian(week,origin = as.Date("2012-01-01")))/365.25 + 2012) 
  omit.num <- which("week" == names(data.air))    ### this column would not be sclaed
  data.air[,-c(4,omit.num)] <- apply(data.air[,-c(4,omit.num)],2,scale)
  
  demo <- read.csv("Data/demo.csv",header = T,stringsAsFactors = F,check.names = F)
  demo$pop <- apply(demo[,-c(1:3)],1,sum); demo <- demo[,c(1,2,4,21)]
  names(demo) <- c("year","city","births","pop")
  demo <- filter(demo,city == city.name)[,-2]
  covar <- plyr::summarize(demo,
                           week = seq(from = min(year) - 1,to = max(year) + 1,by = 1/52),
                           pop = predict(smooth.spline(x = year,y = pop),x = week)$y,
                           birthrate = predict(smooth.spline(x = year + 0.5,y = births),x = week - 4)$y)
 
  covar <- covar[c((dim(covar)[1] - dim(data.case.name)[1] + 1):dim(covar)[1]),]
  
  covar$week <- data.case.name$week; covar <- merge(covar,data.air,by = "week",all.x = T)
  date.num <- which("date" == names(covar));covar <- covar[,-c(date.num)]
  covar <- merge(covar,data.case.name,by = "week",all.x = T)
  vaccine.coverage <- city.data[city.data$name == region,3]/100
  covar$vaccine <- rep((1 - vaccine.coverage*0.8),nrow(covar))
  return(covar)
}

#create_covars(data.set = tem.air,
#              foi.threshold.name = "MT",
#              foi.threshold.value = 40,
#              foi.threshold.type = "V",
#              rr.threshold.name = "MT",
#              rr.threshold.value = c(9,20),
#              rr.threshold.type = "U",
#              foi.continue = F)