library(ggplot2);library(reshape);library(magrittr);library(dplyr);
library(stringi);library(scales);library(locfit);library(cowplot)
library(ggsci)
source("Results-Visualization/mcap.R")

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




Profile.Full.Plot <- function(dataset,smooth.type = F,profile.figure.dir= NULL,
                              alpha, p.value = 0.95, seq.long,omit = TRUE,
                              macp = T,
                              parm.name = c("bmt","bah","kmt","krh","gmt","grh")){
  
  ### dataset <- dat.all.easy; smooth.type =T; alpha = 0.8; p.value = 0.95; seq.long = 0.01  
  ### parm.name <- c("bah","bmt","kmt","krh","gmt","grh");
  
  cols <- pal_lancet("lanonc",alpha = 0.5)(8); #show_col(cols)
  
  type.name <- unique(dataset$type)
  
  if(is.null(profile.figure.dir)){
    profile.figure.dir <- getwd()
  } 
  
  for (j in c(1:length(type.name))) {
    
    plot.total <- list()

    for (i in (1:length(parm.name))){
      dat.all.1 <- filter(dataset, (parm == parm.name[i]) & (type == type.name[j]))  
      
    if(isTRUE(omit)){
      if (parm.name[i] %in% c("bmt","bah")) {
        dat.all.1 <- filter(dat.all.1,target > log(0.5) & (target < log(1.11))) 
      }
    }
      if(isTRUE(macp)){
        max.x  <-  dat.all.1[which.max(dat.all.1$loglik),"target"];
        dat.ci <- mcap(dat.all.1$loglik,parameter = dat.all.1$target)$ci
        dat.plot <- mcap(dat.all.1$loglik,parameter = dat.all.1$target)$fit
        names(dat.plot) <- c("x","y","z")
        x.value.1 <- dat.ci[1];  x.value.2 <- dat.ci[2];
        dat.all.1 <- mutate(dat.all.1,LL = (loglik + 2*loglik.se),UL = (loglik - 2*loglik.se))
        lik.95 <-  max(dat.plot$y) - 0.5 * qchisq(p = p.value,df = 1)
      }else{
      fit <- locfit.robust(loglik~target,data = dat.all.1,family = "qrgauss",alpha = alpha);
      dat.all.1 <- mutate(dat.all.1,LL = (loglik + 2*loglik.se),UL = (loglik - 2*loglik.se))
      dat.plot <- data.frame(y = predict(fit,dat.all.1$target),x = dat.all.1$target)
      
      if(smooth.type){
        max.x <- dat.plot[which.max(dat.plot$y),"x"]
        lik.95 <-  max(dat.plot$y) - 0.5 * qchisq(p = p.value,df = 1)
        x.l.value <- seq(range(dat.all.1$target)[1],max.x,seq.long);
        x.h.value <- seq(max.x,range(dat.all.1$target)[2],seq.long);
        n.num.1 <- which.min(abs(predict(fit,newdata = x.l.value) - lik.95));
        n.num.2 <- which.min(abs(predict(fit,newdata = x.h.value) - lik.95));
        x.value.1 <- x.l.value[n.num.1]; x.value.2 <- x.h.value[n.num.2]; 
      }else{
        max.x  <-  dat.all.1[which.max(dat.all.1$loglik),"target"];
        lik.95 <-  max(dat.all.1$loglik) - 0.5*qchisq(p = p.value,df = 1)
        x.l.value <- seq(range(dat.all.1$target)[1],max.x,seq.long);
        x.h.value <- seq(max.x,range(dat.all.1$target)[2],seq.long);
        n.num.1 <- which.min(abs(predict(fit,newdata = x.l.value) - lik.95));
        n.num.2 <- which.min(abs(predict(fit,newdata = x.h.value) - lik.95));
        x.value.1 <- x.l.value[n.num.1]; x.value.2 <- x.h.value[n.num.2];
      }
      }
      
      if (all(dat.all.1$parm == "bah")) {
      #  xlab.label <- xlab(expression(paste("AH on TR: ", eta[AH],sep = "")))
        xlab.label <- xlab(expression(paste("Absolute Humidity on transmission rate (", eta[AH],")",sep = "")))
      } else if(all(dat.all.1$parm == "bmt")){
       # xlab.label <- xlab(expression(paste("MT on TR: ",eta[MT],sep = "")))
        xlab.label <- xlab(expression(paste("Mean temperature on transmission rate (", eta[MT],")",sep = "")))
      } else if(all(dat.all.1$parm == "gmt")){
       # xlab.label <- xlab(expression(paste(" MT on RR: ",psi[MT],sep = "")))
        xlab.label <- xlab(expression(paste("Mean temperature on reporting rate (", psi[MT],")",sep = "")))
      } else{
       # xlab.label <- xlab(expression(paste(" RH on RR: ", psi[RH],sep = "")))
        xlab.label <- xlab(expression(paste("Relative humidity on reporting rate (", psi[RH],")",sep = "")))
      }
      
      if(parm.name[i] %in% c("kmt","krh")) {
        vline.code <- NULL
        title.code <- NULL
      } else {
        vline.code <-  geom_vline(xintercept = c(x.value.1,x.value.2,max.x),colour=cols[2],linetype = c(4,4,1)) 
        title.code <-  labs(title = paste("MLE: ",round(max.x,2)," (",round(x.value.1,2),", ",round(x.value.2,2),")",sep = "")) 
      }
      
      plot.total[[i]] <- ggplot(data=dat.plot,aes(x = x,y = y)) + geom_line(linetype = 1,alpha = 0.5) + 
        theme_bw(base_family = "serif",base_size = 10) + 
        geom_pointrange(data = dat.all.1,aes(x = target,y = loglik,ymin = LL,ymax = UL),colour = cols[3],
                        fill = "yellow",size=0.2,linetype = 1) +
        geom_hline(yintercept=lik.95,colour="black",linetype = 2,alpha = 0.5) +
        vline.code + title.code +
        xlab.label + ylab("Profile log likehood")
    }
    plot_grid(plot.total[[1]],plot.total[[2]],plot.total[[3]],
              plot.total[[4]],plot.total[[5]],plot.total[[6]],
              labels = LETTERS[1:length(parm.name)],ncol = 2)
    
    ggsave(paste(profile.figure.dir,type.name[j],".pdf",sep = ""),
           dpi = 300,width = 18,height = 14,units = "cm")
  }
}






Profile.Easy.Plot <- function(dataset,smooth.type = F,profile.figure.dir = NULL,
                              alpha = 0.95, p.value, seq.long,mcap.t =T,
                              parm.name = c("bah","bmt","gmt","grh")){
  
  ### dataset <- dat.all.easy; smooth.type =T; alpha = 0.8; p.value = 0.95; seq.long = 0.01  
  ### parm.name <- c("bah","bmt","kmt","krh","gmt","grh");
  
  cols <- pal_lancet("lanonc",alpha = 0.5)(8); #show_col(cols)
  
  type.name <- unique(dataset$type)
  
  if(is.null(profile.figure.dir)){
    profile.figure.dir <- getwd()
  } 
  
  for (j in c(1:8)) {
    
    plot.total <- list()
    
    for (i in (1:length(parm.name))){
      dat.all.1 <- filter(dataset, (parm == parm.name[i]) & (type == type.name[j]))  
      if (parm.name[i] %in% c("bmt","bah")) {
        dat.all.1 <- filter(dat.all.1,target > log(0.5)) 
      }
      if(isTRUE(mcap.t)){
        max.x  <-  dat.all.1[which.max(dat.all.1$loglik),"target"];
        dat.ci <- mcap(dat.all.1$loglik,parameter = dat.all.1$target)$ci
        dat.plot <- mcap(dat.all.1$loglik,parameter = dat.all.1$target)$fit
        names(dat.plot) <- c("x","y","z")
        x.value.1 <- dat.ci[1];  x.value.2 <- dat.ci[2];
        dat.all.1 <- mutate(dat.all.1,LL = (loglik + 2*loglik.se),UL = (loglik - 2*loglik.se))
        lik.95 <-  max(dat.plot$y) - 0.5 * qchisq(p = p.value,df = 1)
        
      } else{
        
        fit <- locfit.robust(loglik~target,data = dat.all.1,family = "qrgauss",alpha = alpha);
        dat.all.1 <- mutate(dat.all.1,LL = (loglik + 2*loglik.se),UL = (loglik - 2*loglik.se))
        dat.plot <- data.frame(y = predict(fit,dat.all.1$target),x = dat.all.1$target)
        
        if(isTRUE(smooth.type)){
          max.x <- dat.plot[which.max(dat.plot$y),"x"]
          lik.95 <-  max(dat.plot$y) - 0.5 * qchisq(p = p.value,df = 1)
          x.l.value <- seq(range(dat.all.1$target)[1],max.x,seq.long);
          x.h.value <- seq(max.x,range(dat.all.1$target)[2],seq.long);
          n.num.1 <- which.min(abs(predict(fit,newdata = x.l.value) - lik.95));
          n.num.2 <- which.min(abs(predict(fit,newdata = x.h.value) - lik.95));
          x.value.1 <- x.l.value[n.num.1]; x.value.2 <- x.h.value[n.num.2]; 
        }else{
          max.x  <-  dat.all.1[which.max(dat.all.1$loglik),"target"];
          lik.95 <-  max(dat.all.1$loglik) - 0.5*qchisq(p = p.value,df = 1)
          x.l.value <- seq(range(dat.all.1$target)[1],max.x,seq.long);
          x.h.value <- seq(max.x,range(dat.all.1$target)[2],seq.long);
          n.num.1 <- which.min(abs(predict(fit,newdata = x.l.value) - lik.95));
          n.num.2 <- which.min(abs(predict(fit,newdata = x.h.value) - lik.95));
          x.value.1 <- x.l.value[n.num.1]; x.value.2 <- x.h.value[n.num.2];
        }}
      
      if (all(dat.all.1$parm == "bah")) {
        #  xlab.label <- xlab(expression(paste("AH on TR: ", eta[AH],sep = "")))
        xlab.label <- xlab(expression(paste("Absolute Humidity on transmission rate (", eta[AH],")",sep = "")))
      } else if(all(dat.all.1$parm == "bmt")){
        # xlab.label <- xlab(expression(paste("MT on TR: ",eta[MT],sep = "")))
        xlab.label <- xlab(expression(paste("Mean temperature on transmission rate (", eta[MT],")",sep = "")))
      } else if(all(dat.all.1$parm == "gmt")){
        # xlab.label <- xlab(expression(paste(" MT on RR: ",psi[MT],sep = "")))
        xlab.label <- xlab(expression(paste("Mean temperature on reporting rate (", psi[MT],")",sep = "")))
      } else{
        # xlab.label <- xlab(expression(paste(" RH on RR: ", psi[RH],sep = "")))
        xlab.label <- xlab(expression(paste("Relative humidity on reporting rate (", psi[RH],")",sep = "")))
      }
      
      plot.total[[i]] <- ggplot(data=dat.plot,aes(x = x,y = y)) + geom_line(linetype = 1,alpha = 0.5) + 
        theme_bw(base_family = "serif",base_size = 10) + 
        geom_pointrange(data=dat.all.1,aes(x = target,y = loglik,ymin = LL,ymax = UL),colour = cols[3],
                        fill="yellow",size=0.2,linetype = 1) +
        geom_hline(yintercept=lik.95,colour="black",linetype = 2,alpha = 0.5) +
        geom_vline(xintercept = c(x.value.1,x.value.2,max.x),colour=cols[2],linetype = c(4,4,1)) +
        labs(title = paste("MLE: ",fill.character(round(max.x,2),4)," (",
                           fill.character(round(x.value.1,2),4),", ",
                           fill.character(round(x.value.2,2),4),")",sep = "")) + 
        xlab.label + ylab("Profile log likehood")
    }
    plot_grid(plot.total[[1]],plot.total[[2]],plot.total[[3]],
              plot.total[[4]],
              labels = LETTERS[1:length(parm.name)],ncol = 2)
    
    ggsave(paste(profile.figure.dir,type.name[j],".pdf",sep = ""),
           dpi = 300,width = 18,height = 14,units = "cm")
  }
}




Omit.extreate.value <- function(dataset,range.value.1= 8, range.value.2 = 5){
  dataset <- mutate(dataset,id = paste(type,parm,sep = ""))
  type.name <- unique(dataset$id)
  result.list <- lapply(1:length(type.name), function(id){
    dat.tem <- dataset[dataset$id == type.name[id],]
    dat.tem <- arrange(dat.tem,target); length.n <- nrow(dat.tem)
    max.value <- max(dat.tem$loglik); max.n <- which.max(dat.tem$loglik)
    for (i in 1:(length.n - 3)){
      frist.loglik  <- dat.tem[i,"loglik"]; 
      second.loglik <- dat.tem[i + 1,"loglik"]; 
      third.loglik  <- dat.tem[i + 2,"loglik"];
      fouth.loglik  <- dat.tem[i + 3,"loglik"];
      mean.loglik.1 <- (frist.loglik + third.loglik)/2
      mean.loglik.2 <- (second.loglik + fouth.loglik)/2
      
      value.1 <- abs(second.loglik) > (abs(mean.loglik.1) + range.value.1)
      value.2 <- (abs(third.loglik) - (abs(fouth.loglik)))  > range.value.2
      if(value.1 == TRUE & value.2 == TRUE){
        dat.tem[i + 2,"loglik1"] <-  mean.loglik.2  
        dat.tem[i + 2, "normal"] <- 1
      } else if(value.1 == TRUE ){
        dat.tem[i + 1,"loglik1"] <-  mean.loglik.1
        dat.tem[i + 1, "normal"] <- 1
      } else {
        dat.tem[i + 1,"loglik1"] <- dat.tem[i + 1,"loglik"];
        dat.tem[i + 2,"loglik1"] <- dat.tem[i + 2,"loglik"] 
        dat.tem[i + 1, "normal"] <- 0
      }}
    return(dat.tem)
  })
  result <- do.call("rbind",result.list)
  return(result)
}



Point_CI_Plot <- function(dataset,figure.dir,city = F,city.color = c("#925E9FB2","#42B540B2")){
  ##  dataset <- result.easy; type.name = "easy"
  plot.total <- list(); #dataset$parm <- as.character(dataset$parm)
  parameter.name <- c("bmt","bah","gmt","grh")
  for (i in c(1:length(parameter.name))){
    name.parm <- as.character(parameter.name[i])
    dat.tem <- filter(dataset, parm == name.parm)
    
    if (name.parm == "bah"){
     # ylab.label <- ylab(expression(paste("AH on TR: ",eta[AH],sep = "")))
      ylab.label <- ylab(expression(paste("Absolute Humidity on transmission rate (", eta[AH],")",sep = "")))
    } else if(name.parm == "bmt"){
      # ylab.label <- ylab(expression(paste("MT on TR: ",eta[MT],sep = "")))
      ylab.label <- ylab(expression(paste("Mean temperature on transmission rate (", eta[MT],")",sep = "")))
      
    } else if(name.parm == "gmt"){
      # ylab.label <- ylab(expression(paste("MT on RR: ",psi[MT],sep = "")))
      ylab.label <- ylab(expression(paste("Mean temperature on reporting rate (", psi[MT],")",sep = "")))
    } else{
     # ylab.label <- ylab(expression(paste("RH on RR: ", psi[RH],sep = "")))
      ylab.label <- ylab(expression(paste("Relative humidity on reporting rate (", psi[RH],")",sep = "")))
    }
    
    if (name.parm %in% c("bah","bmt")){
      theme.tem <- theme(axis.text.x = element_blank(),
                         axis.ticks.x = element_blank(),
                         legend.spacing.x = unit(8,"cm"),
                         panel.grid.minor.x = element_blank())
    } else {
      theme.tem <- theme(axis.text.x = element_text(angle = -90,vjust = 0.6),
                         legend.spacing.y = unit(2,"cm"),
                         panel.grid.minor.x = element_blank())
    }
    
    if(isTRUE(city)){
      plot.total[[i]] <-  ggplot() +
        geom_pointrange(data = dataset[dataset$parm == name.parm,], 
                        aes(x = type, y = max.x, 
                            ymin = x.value.1, ymax = x.value.2, group = factor(city),
                            colour = factor(city)),fatten = 0.6,
                         position=position_dodge(width=0.40),size = 0.6) +
        geom_hline(yintercept = 0,linetype = 2) +
        scale_color_manual(values = city.color,
                           guide = guide_legend(title = "",keywidth = 1.5,keyheight = 1,
                                                override.aes = list(size = 0.5))) +
     #   geom_point(data = dataset[dataset$parm == name.parm,],aes( x= type, y =max.x),color = "red") +
        theme_bw(base_size = 10, base_family = "serif") +
        xlab("") + ylab.label + theme.tem
    } else{
      plot.total[[i]] <-  ggplot() +
        geom_pointrange(data = dataset[dataset$parm == name.parm,], 
                        aes(x = type, y = max.x, 
                            ymin = x.value.1, ymax = x.value.2, group = factor(type)),size = 0.2) +
        geom_hline(yintercept = 0,linetype = 2) + 
      #  geom_point(data = dataset[dataset$parm == name.parm,],aes( x = type, y = max.x),color = "red") +
        theme_bw(base_size = 10, base_family = "serif") +
        xlab("") + ylab.label + theme.tem
      
    }
  }
  legend_b  <- get_legend(plot.total[[1]] + theme(legend.position = "bottom",
                                                  legend.direction = "horizontal",
                                                  legend.justification = c(0.5,-2.5)))
  figure.1 <- plot_grid(plot.total[[1]] + theme(legend.position="none"),
              plot.total[[2]] + theme(legend.position="none"),
              plot.total[[3]] + theme(legend.position="none"),
              plot.total[[4]] + theme(legend.position="none"),
              labels = c("A", "B","C","D"),ncol = 2,rel_heights = c(1,1.2))#,vjust = -0.3) #,rel_heights = c(1,1.5))
  figure.2 <- plot_grid(figure.1, legend_b, ncol = 1, rel_heights = c(2, .1))
  ggsave(paste(figure.dir,"/pointci.pdf",sep = ""),figure.2,
         dpi = 300,width = 18,height = 14,units = "cm")
  print(paste("Figure was ",paste(figure.dir,".pdf",sep = ""),sep = ""))
}



ci.calculate <- function(dataset,p.value,alpha,seq.long,smooth.type = TRUE,mcap.t =T){
  
  dataset <- mutate(dataset,idtp = paste(type,parm,sep = ""))
  type.name <- unique(dataset$idtp);
  result.list <- lapply(type.name, function(data){
    dat.tem <- filter(dataset, idtp == data)

    if(isTRUE(mcap.t)){
      max.x  <-  dat.tem[which.max(dat.tem$loglik),"target"];
      dat.ci <- mcap(dat.tem$loglik,parameter = dat.tem$target)$ci
      dat.plot <- mcap(dat.tem$loglik,parameter = dat.tem$target)$fit
      names(dat.plot) <- c("x","y","z")
      max.x.2 <- dat.plot[which(dat.plot$y == max(dat.plot$y)),1]
      x.value.1 <- dat.ci[1];  x.value.2 <- dat.ci[2];
      lik.95 <-  max(dat.plot$y) - 0.5 * qchisq(p = p.value,df = 1)
    } else {
    fit <- locfit.robust(loglik~target,data = dat.tem,family = "qrgauss",alpha = alpha);
    predict.x <- seq(range(dat.tem$target)[1],range(dat.tem$target)[2],seq.long)
    dat.plot <- data.frame(y = predict(fit,predict.x),x = predict.x)
    
    if(smooth.type){
      max.x  <-  dat.tem[which.max(dat.tem$loglik),"target"];
      max.x.2 <- max.x;
      lik.95 <-  max(dat.tem$loglik) - 0.5*qchisq(p = p.value,df = 1)
      x.l.value <- seq(range(dat.tem$target)[1],max.x,seq.long);
      x.h.value <- seq(max.x,range(dat.tem$target)[2],seq.long);
      n.num.1 <- which.min(abs(predict(fit,newdata = x.l.value) - lik.95));
      n.num.2 <- which.min(abs(predict(fit,newdata = x.h.value) - lik.95));
      x.value.1 <- x.l.value[n.num.1]; x.value.2 <- x.h.value[n.num.2]; 
    } else {
      max.x  <-  dat.tem[which.max(dat.tem$loglik),"target"];
      max.x.2 <- dat.plot[which.max(dat.plot$y),"x"]
      lik.95 <-  max(dat.plot$y) - 0.5 * qchisq(p = p.value,df = 1)
      x.l.value <- seq(range(dat.tem$target)[1],max.x.2,seq.long);
      x.h.value <- seq(max.x.2,range(dat.tem$target)[2],seq.long);
      n.num.1 <- which.min(abs(predict(fit,newdata = x.l.value) - lik.95));
      n.num.2 <- which.min(abs(predict(fit,newdata = x.h.value) - lik.95));
      x.value.1 <- x.l.value[n.num.1]; x.value.2 <- x.h.value[n.num.2]; 
    }
    }
    result <- data.frame(max.x,max.x.2,x.value.1,x.value.2,dat.tem[1,"type"],dat.tem[1,"parm"])
    return(result)
  })
  result.final <- do.call("rbind",result.list)
  return(result.final)
}

fill.character <- function(data,number){
  origin.length <- nchar(data)
  
  if (data > 0){
  if(stri_detect_fixed(data,".")){
    zero.c <- ifelse(origin.length < (number),
                     stri_sub("1000000000", (-1 - (number - origin.length - 1)),-1),
                     "")
    result <- paste(data,zero.c,sep = "")    
  } else {
    zero.c <- stri_sub("1000000000", (-1 - (number - origin.length - 2)),-1)
    result <- paste(data,".",zero.c,sep = "")
  }}
  
  if (data < 0){
    if(stri_detect_fixed(data,".")){
      zero.c <- ifelse(origin.length < (number+1),
                       stri_sub("1000000000", (-1 - (number - origin.length)),-1),
                       "")
      result <- paste(data,zero.c,sep = "")    
    } else {
      zero.c <- stri_sub("1000000000", (-1 - (number - origin.length - 1)),-1)
      result <- paste(data,".",zero.c,sep = "")
    }
  }
  
  
  if (data == 0){
      zero.c <- stri_sub("1000000000", (-1 - (number - origin.length - 2)),-1)
      result <- paste(data,".",zero.c,sep = "")
  }
  return(result)
}




Point_CI_Plot.2 <- function(dataset,figure.dir,city = F,city.color = c("#925E9FB2","#42B540B2")){
  ##  dataset <- result.easy; type.name = "easy"
  plot.total <- list(); #dataset$parm <- as.character(dataset$parm)
  parameter.name <- c("bmt","bah","gmt","grh")
  for (i in c(1:length(parameter.name))){
    name.parm <- as.character(parameter.name[i])
    dat.tem <- filter(dataset, parm == name.parm)
    
    if (name.parm == "bah"){
      # ylab.label <- ylab(expression(paste("AH on TR: ",eta[AH],sep = "")))
      ylab.label <- ylab(expression(paste("Absolute Humidity on transmission rate (", eta[AH],")",sep = "")))
    } else if(name.parm == "bmt"){
      # ylab.label <- ylab(expression(paste("MT on TR: ",eta[MT],sep = "")))
      ylab.label <- ylab(expression(paste("Mean temperature on transmission rate (", eta[MT],")",sep = "")))
      
    } else if(name.parm == "gmt"){
      # ylab.label <- ylab(expression(paste("MT on RR: ",psi[MT],sep = "")))
      ylab.label <- ylab(expression(paste("Mean temperature on reporting rate (", psi[MT],")",sep = "")))
    } else{
      # ylab.label <- ylab(expression(paste("RH on RR: ", psi[RH],sep = "")))
      ylab.label <- ylab(expression(paste("Relative humidity on reporting rate (", psi[RH],")",sep = "")))
    }
    
    if (name.parm %in% c("bah","bmt")){
      theme.tem <- theme(axis.text.x = element_blank(),
                         axis.ticks.x = element_blank(),
                         legend.spacing.x = unit(8,"cm"),
                         panel.grid.minor.x = element_blank())
    } else {
      theme.tem <- theme(axis.text.x = element_text(angle = -90,vjust = 0.6),
                         legend.spacing.y = unit(2,"cm"),
                         panel.grid.minor.x = element_blank())
    }
    
    if(isTRUE(city)){
      plot.total[[i]] <-  ggplot() +
        geom_pointrange(data = dataset[dataset$parm == name.parm,], 
                        aes(x = type, y = max.x, 
                            ymin = x.value.1, ymax = x.value.2, group = factor(city),
                            colour = factor(city)),fatten = 0.6,
                        position=position_dodge(width=0.40),size = 0.6) +
        geom_hline(yintercept = 0,linetype = 2) +
        scale_color_manual(values = city.color,
                           guide = guide_legend(title = "",keywidth = 1.5,keyheight = 1,
                                                override.aes = list(size = 0.5))) +
        #   geom_point(data = dataset[dataset$parm == name.parm,],aes( x= type, y =max.x),color = "red") +
        theme_bw(base_size = 10, base_family = "serif") +
        xlab("") + ylab.label + theme.tem
    } else{
      plot.total[[i]] <-  ggplot() +
        geom_pointrange(data = dataset[dataset$parm == name.parm,], 
                        aes(x = type, y = max.x, 
                            ymin = x.value.1, ymax = x.value.2, group = factor(type)),size = 0.2) +
        geom_hline(yintercept = 0,linetype = 2) + 
        #  geom_point(data = dataset[dataset$parm == name.parm,],aes( x = type, y = max.x),color = "red") +
        theme_bw(base_size = 10, base_family = "serif") +
        xlab("") + ylab.label + theme.tem
      
    }
  }
  legend_b  <- get_legend(plot.total[[1]] + theme(legend.position = "bottom",
                                                  legend.direction = "horizontal",
                                                  legend.justification = c(0.5,-2.5)))
  figure.1 <- plot_grid(plot.total[[1]] + theme(legend.position="none"),
                        plot.total[[2]] + theme(legend.position="none"),
                        plot.total[[3]] + theme(legend.position="none"),
                        plot.total[[4]] + theme(legend.position="none"),
                        labels = c("A", "B","C","D"),ncol = 2,rel_heights = c(1,1.2))#,vjust = -0.3) #,rel_heights = c(1,1.5))
  figure.2 <- plot_grid(figure.1, legend_b, ncol = 1, rel_heights = c(2, .1))
  ggsave(paste(figure.dir,"/pointci.pdf",sep = ""),figure.2,
         dpi = 300,width = 16,height = 14,units = "cm")
  print(paste("Figure was ",paste(figure.dir,".pdf",sep = ""),sep = ""))
}




Profile.BS.Plot <- function(dataset,smooth.type = F,
                              alpha = 0.95, p.value, seq.long,mcap.t =T,
                              parm.name = c("gmt","grh")){
  
    cols <- pal_lancet("lanonc",alpha = 0.5)(8); #show_col(cols)
    dataset$target <- log(dataset$target)
    
    plot.total <- list()
    
    for (i in (1:length(parm.name))){
      dat.all.1 <- filter(dataset, (parm == parm.name[i]))  

      if(isTRUE(mcap.t)){
        max.x  <-  dat.all.1[which.max(dat.all.1$loglik),"target"];
        dat.ci <- mcap(dat.all.1$loglik,parameter = dat.all.1$target)$ci
        dat.plot <- mcap(dat.all.1$loglik,parameter = dat.all.1$target)$fit
        names(dat.plot) <- c("x","y","z")
        x.value.1 <- dat.ci[1];  x.value.2 <- dat.ci[2];
        dat.all.1 <- mutate(dat.all.1,LL = (loglik + 2*loglik.se),UL = (loglik - 2*loglik.se))
        lik.95 <-  max(dat.plot$y) - 0.5 * qchisq(p = p.value,df = 1)
        
      } else{
        
        fit <- locfit.robust(loglik~target,data = dat.all.1,family = "qrgauss",alpha = alpha);
        dat.all.1 <- mutate(dat.all.1,LL = (loglik + 2*loglik.se),UL = (loglik - 2*loglik.se))
        dat.plot <- data.frame(y = predict(fit,dat.all.1$target),x = dat.all.1$target)
        
        if(isTRUE(smooth.type)){
          max.x <- dat.plot[which.max(dat.plot$y),"x"]
          lik.95 <-  max(dat.plot$y) - 0.5 * qchisq(p = p.value,df = 1)
          x.l.value <- seq(range(dat.all.1$target)[1],max.x,seq.long);
          x.h.value <- seq(max.x,range(dat.all.1$target)[2],seq.long);
          n.num.1 <- which.min(abs(predict(fit,newdata = x.l.value) - lik.95));
          n.num.2 <- which.min(abs(predict(fit,newdata = x.h.value) - lik.95));
          x.value.1 <- x.l.value[n.num.1]; x.value.2 <- x.h.value[n.num.2]; 
        }else{
          max.x  <-  dat.all.1[which.max(dat.all.1$loglik),"target"];
          lik.95 <-  max(dat.all.1$loglik) - 0.5*qchisq(p = p.value,df = 1)
          x.l.value <- seq(range(dat.all.1$target)[1],max.x,seq.long);
          x.h.value <- seq(max.x,range(dat.all.1$target)[2],seq.long);
          n.num.1 <- which.min(abs(predict(fit,newdata = x.l.value) - lik.95));
          n.num.2 <- which.min(abs(predict(fit,newdata = x.h.value) - lik.95));
          x.value.1 <- x.l.value[n.num.1]; x.value.2 <- x.h.value[n.num.2];
        }}
      
      if (all(dat.all.1$parm == "bah")) {
        #  xlab.label <- xlab(expression(paste("AH on TR: ", eta[AH],sep = "")))
        xlab.label <- xlab(expression(paste("Absolute Humidity on transmission rate (", eta[AH],")",sep = "")))
      } else if(all(dat.all.1$parm == "bmt")){
        # xlab.label <- xlab(expression(paste("MT on TR: ",eta[MT],sep = "")))
        xlab.label <- xlab(expression(paste("Mean temperature on transmission rate (", eta[MT],")",sep = "")))
      } else if(all(dat.all.1$parm == "gmt")){
        # xlab.label <- xlab(expression(paste(" MT on RR: ",psi[MT],sep = "")))
        xlab.label <- xlab(expression(paste("Mean temperature on reporting rate (", psi[MT],")",sep = "")))
      } else{
        # xlab.label <- xlab(expression(paste(" RH on RR: ", psi[RH],sep = "")))
        xlab.label <- xlab(expression(paste("Relative humidity on reporting rate (", psi[RH],")",sep = "")))
      }
      
      plot.total[[i]] <- ggplot(data=dat.plot,aes(x = x,y = y)) + geom_line(linetype = 1,alpha = 0.5) + 
        theme_bw(base_family = "serif",base_size = 10) + 
        geom_pointrange(data=dat.all.1,aes(x = target,y = loglik,ymin = LL,ymax = UL),colour = cols[3],
                        fill="yellow",size=0.2,linetype = 1) +
        geom_hline(yintercept=lik.95,colour="black",linetype = 2,alpha = 0.5) +
        geom_vline(xintercept = c(x.value.1,x.value.2,max.x),colour=cols[2],linetype = c(4,4,1)) +
        labs(title = paste("MLE: ",fill.character(round(max.x,2),4)," (",
                           fill.character(round(x.value.1,2),4),", ",
                           fill.character(round(x.value.2,2),4),")",sep = "")) + 
        xlab.label + ylab("Profile log likehood")
    }
    plot_grid(plot.total[[1]],plot.total[[2]],
              labels = LETTERS[1:length(parm.name)],ncol = 2)
  }
