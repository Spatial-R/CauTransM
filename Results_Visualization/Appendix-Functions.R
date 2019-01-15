##########################################################################################
################################### Codes for the appendix  ##############################
##########################################################################################


Inter.ci <- function(dir,p.value,type){
  #### @ dir: the directory of the dataset   
  #### @ p.value: the probability of value to calculate the CI
  #### @ type: which to show? the value with the maximum likilihood
  
  all.files <- list.files(path = dir,full.names = T,pattern = ".csv")
  mylist <- lapply(all.files, function(i){
    dat.tem  <- read.csv(i); max.loglik <- max(dat.tem$loglik);
    lik.99 <- max.loglik - 0.5*qchisq(p = p.value,df = 1);
    dat.99 <- filter(dat.tem, loglik > lik.99)[,c(15:23,25)];
    n <- which.max(dat.99$loglik);
    if (type == "quantile") {
      result <- data.frame(t(apply(dat.99[,-10],2, function(data){
        median.result <- round(quantile(log(data),probs = 0.5),2);
        low.result    <- round(quantile(log(data),probs = 0.025),2);
        high.result   <- round(quantile(log(data),probs = 0.975),2);
        paste(median.result,"(",low.result,",",high.result,")",sep = "")
      })))} else {
        result <- data.frame(t(apply(dat.99[,-10],2, function(data){
          median.result <- round(log(data[n]),2);
          low.result    <- round(quantile(log(data),probs = 0.025),2);
          high.result   <- round(quantile(log(data),probs = 0.975),2);
          paste(median.result,"(",low.result,",",high.result,")",sep = "")
        })))}
    name <- gsub(".result.csv","",gsub(paste(dir,"/","global",sep = ""),"",i))
    result <- mutate(result,type = name,loglik = max.loglik,number = dim(dat.99)[1],
                     AIC=AICc(max.loglik,K = 22,n = 157))
    return(result)
  })
  result.final <- do.call("rbind",mylist);return(result.final)
}


Global.Paramter.Estimate <- function(dir,log=TRUE){
  #### dir: the directory of the dataset   
  all.files <- list.files(path = dir,full.names = T,pattern = ".csv")
  mylist <- lapply(all.files,function(i){
    dat.tem  <- read.csv(i); max.loglik <- max(dat.tem$loglik);
    name <- gsub(".result.csv","",gsub(paste(dir,"/","global",sep = ""),"",i))
    dat.99 <- filter(dat.tem, loglik == max.loglik); 
    dat.99 <- mutate(dat.99,type = name, loglik = round(loglik,3),
                     AIC = round(AICc(max.loglik,K = 22,n = 157),3))
    if(log == TRUE)  
      dat.99[,15:23] <- apply(dat.99[,15:23],2,function(data){
        round(log(data),3)
      }) 
    return(dat.99)
  })
  result.final <- do.call("rbind",mylist);return(result.final)
}



panel.hist <- function(x, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}



Global.Paramter.Paris <- function(dir,climatic=TRUE, maxlogik = 20){
  
  #### dir: the directory of the dataset   
  #### all: whether show all the parameter value or just the climatic paramter
  #### maxloglik: max(loglik) - maxlogik for the global searching 
  
  all.files <- list.files(path = dir,full.names = T,pattern = ".csv")
  lapply(all.files,function(i){
    dat.tem  <- read.csv(i); max.loglik <- max(dat.tem$loglik);
    name <- gsub(".result.csv","",gsub(dir,"",i)); 
    name <- gsub("-globallinefull","",name); name <- gsub("/","",name)
    dat.pars <- filter(dat.tem,loglik > (max.loglik - maxlogik))[,-c(1,23:25)]
    
    if (climatic == TRUE) {
      
      dat.pars[,15:20] <- apply(dat.pars[,15:20],2,function(data){
        round(log(data),3)
      })
      cairo_pdf(filename = paste(figure.app,"/",name,"climatic.pdf",sep=""),
                width = 8, height = 8, onefile = TRUE, family = "serif")
      pairs(~loglik+bmt+bah+kmt+krh+gmt+grh,data = dat.pars, 
            cex = 0.8, pch = 1, bg = "light blue",
            diag.panel = panel.hist, cex.labels = 1.5, font.labels = 1.5,
            labels = c("Loglik",
                       expression(beta[MT]),expression(beta[AH]),
                       expression(alpha[MT]),expression(alpha[RH]),
                       expression(psi[MT]),expression(psi[RH])))
      dev.off()
    } else {
      cairo_pdf(filename = paste(figure.app,"/",name,".pdf",sep=""),
                width = 8, height = 8, onefile = TRUE, family = "serif")
      pairs(~loglik+alpha+iota+rho+sigmaSE+psi+amplitude,data = dat.pars, 
            cex = 0.8, pch = 1, bg = "light blue",
            diag.panel = panel.hist, cex.labels = 1.5, font.labels = 1.5,
            labels = c("Loglik",expression(alpha[0]),expression(iota),expression(rho[0]),
                       expression(beta[S.D]),expression(psi),expression(epsilon)))
      dev.off()
    }
  })
}



rmse <- function(data){
  #data <- result.plot.sim
  dat.f <- filter(data,data == FALSE)[,-3]; dat.f$id  <- rep(1:100,each = (nrow(dat.f)/100))
  dat.f2 <- reshape(dat.f, v.names = "cases", idvar = "time",
                    timevar = "id", direction = "wide")
  dat.t <- filter(data,data == TRUE)[,-3]
  dat <- merge(dat.f2,dat.t,all.x = T,by = 'time')[-1,-1]; result <- c()
  for (i in c(1:100)){
    result[i] <- sqrt(sum((dat[,i]-dat[,ncol(dat)])^2)/nrow(dat))
  }
  mean(result)
}

Profile.calculate <- function(dir,city.name.omit = T){
  all.files <- list.files(path = dir,full.names = T,pattern = ".csv");
  dat.all <- lapply(all.files,function(i) {
    dat.tem <- read.csv(i, header = T,stringsAsFactors = F);
    city.n <- gsub(dir,"",i);city.n <- gsub(".csv","",city.n);
    if(isTRUE(city.name.omit)){
      city.n <- gsub("quzhou.","",city.n)
    }
    city.n <- gsub("fullprofile","",city.n); city.n <- gsub("/","",city.n);
    target <- substr(city.n,(stri_locate_first(city.n, fixed = '.')[1] + 1),nchar(city.n))
    type   <- substr(city.n,1,(stri_locate_first(city.n, fixed = '.')[1]) - 1)
    dat.tem$type <- rep(type,dim(dat.tem)[1]);dat.tem$parm <- rep(target,dim(dat.tem)[1])
    dat.tem %>% 
      subset(is.finite(loglik) & nfail.max == 1,-c(nfail.max,nfail.min)) -> dat.tem
    dat.tem$target <- round(dat.tem[,c(target)],3)
    dat.tem$dll    <- dat.tem$loglik - max(dat.tem$loglik)
    dat.tem %>%
      plyr::ddply(~target,subset,loglik == max(loglik)) -> dat.pred
    return(dat.pred)
  })
  dat.all <- do.call("rbind",dat.all)
  profile.tem <- dat.all %>%
    mutate(id = paste(type,parm,target,sep = "")) %>%
    plyr::ddply(~ id,subset,loglik == max(loglik))
  return(profile.tem)
}


rmse_plot <- function(model,theta,ylab.plot = TRUE){
  origin.model <- model
  origin.model %>% 
    simulate(params = theta,nsim = 100,as.data.frame = TRUE,include.data = TRUE) %>%
    subset(select = c(time,sim,cases)) %>%
    mutate(data=sim == "data") %>%
    select(time,cases,data) -> result.plot.sim 
  
  rmse.1 <- rmse(result.plot.sim)
  
  filter(result.plot.sim,data == F)[,-3] %>%
    group_by(time) %>%
    summarise(low = quantile(cases,0.05),
              med = quantile(cases,0.5),
              hi  = quantile(cases,0.95))  -> lh.result
  cases <- filter(result.plot.sim, data == T)[,-c(1,3)]
  result.plot <- cbind(cases,lh.result)
  #y.pos <- as.numeric(quantile(result.plot$hi,probs = 0.95))
  if(isTRUE(ylab.plot)){
    ylab.set <- ylab("Number of weekly cases")
    y.pos <- max(c(result.plot$hi,result.plot$cases)) - 4   
  }else{
    ylab.set <- ylab("")
    y.pos <- max(c(result.plot$hi,result.plot$cases)) - 1
  }
  
  fig.tem <- ggplot(result.plot[-1,],aes(x = time,y = med,ymin = low,ymax = hi))+
    geom_ribbon(alpha = 0.2) + 
    geom_point(aes(y = cases),size = 0.8,color = "red") +
    geom_line(aes(y = med),color = "black") +
    theme_bw(base_family = "serif",base_size = 15) +
    xlab("Date") + ylab.set + 
    geom_label(aes(x = 2014.5,y = y.pos,label = paste("RMSE = ",round(rmse.1,3),sep = "")),
                  family = "serif") +
   theme_classic(base_size = 15,base_family = "serif")
  return(fig.tem)
}


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
    dat.tem <- dat.tem[dat.tem$V01000 == site.code,]
    dat.tem <- dat.tem[,c("V12001","V13003","V13201","date")]
    names(dat.tem) <- c("MT","RH","CR","date")
    dat.tem <- arrange(dat.tem,date) 
    dat.tem[,1:3] <- apply(dat.tem[,1:3],2,na.approx)     ##### interplation for the missing data
    return(dat.tem)
  }  
}

case_manupilation <- function(dir.1, dir.2,city.name){
  
  case.1 <- read_csv(dir.1) %>%
    data.frame() %>%
    group_by(date,city) %>%
    summarize(cases=n()) %>%
    filter((city %in% city.name) & (substr(date,1,4) %in% c(2011:2015))) %>%
    data.frame()
  case.2 <- read.csv(dir.2,header = T,stringsAsFactors = F) #fileEncoding = "gb18030",
  
  names(case.2)[which(names(case.2) == c('发病日期'))] <- c("date")
  names(case.2)[which(names(case.2) == c("现住详细地址"))] <- c("address")
  
  case.2 <- mutate(case.2, date = as.Date(date))
  case.2$city <- substr(case.2$address,4,6)
  case.2$city <- factor(case.2$city,
                        levels=city.ch,
                        labels=c("Hangzhou","Jiaxing","Ningbo","Huzhou","Jinhua","Shaoxing",
                                 "Lishui","Quzhou","Zhoushan","Taizhou","Wenzhou"))
  
  case.2 <- case.2[case.2$date %in% seq(as.Date("2015-01-01"),as.Date("2015-04-20"),by = "day"),
                   c("date","city")]
  case.2 <- summarise(group_by(case.2,date,city),cases = n()) %>%
    filter(city %in% city.name) %>%
    data.frame()
  case.1$date <- as.Date(case.1$date)
  case.dat <- rbind(case.1,case.2)
  return(case.dat)
}


AICc <- function(loglihood, K, n, cc = T) {
  if(isTRUE(cc)){
    aicc <- -2 * loglihood + 2 * K +  2 * K*(K+1)/(n-K-1)
  } else {
    aicc <- -2 * loglihood + 2 * K
  }
  return(aicc)
}