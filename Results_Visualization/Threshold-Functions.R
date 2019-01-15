source("Results-Visualization/mcap.R")


ci.threshold <- function(dataset,parameter,seq.long,p.value,smooth.type,alpha,mcap.t =T){
  result.list <- lapply(unique(dataset$type),function(id){
    dat.tem <- filter(dataset,type == id) 
    
    if(isTRUE(mcap.t)){
      max.x  <-  dat.tem[which.max(dat.tem$loglik),"target"];
      dat.ci <- mcap(dat.tem$loglik,parameter = dat.tem$target)$ci
      dat.plot <- mcap(dat.tem$loglik,parameter = dat.tem$target)$fit
      names(dat.plot) <- c("x","y","z")
      max.x.2 <- dat.plot[which(dat.plot$y == max(dat.plot$y)),1]
      x.value.1 <- dat.ci[1];  x.value.2 <- dat.ci[2];
      lik.95 <-  max(dat.plot$y) - 0.5 * qchisq(p = p.value,df = 1)
    } else {
    
    formulation <- as.formula(loglik~target)
    fit <- locfit.robust(formulation,data = dat.tem,alpha = alpha);
    
    if(smooth.type){
      predict.x <- seq(range(dat.tem$target)[1],range(dat.tem$target)[2],seq.long)
      dat.plot <- data.frame(y = predict(fit,predict.x),x = predict.x)
      max.x  <-  dat.tem[which.max(dat.tem$loglik),"target"];
      max.x.2 <- dat.plot[which.max(dat.plot$y),"x"]
      lik.95 <-  max(dat.plot$y) - 0.5 * qchisq(p = p.value,df = 1)
      x.l.value <- seq(range(dat.tem$target)[1],max.x.2,seq.long);
      x.h.value <- seq(max.x.2,range(dat.tem$target)[2],seq.long);
      n.num.1 <- which.min(abs(predict(fit,newdata = x.l.value) - lik.95));
      n.num.2 <- which.min(abs(predict(fit,newdata = x.h.value) - lik.95));
      x.value.1 <- x.l.value[n.num.1]; x.value.2 <- x.h.value[n.num.2]; 
    }else{
      max.x  <-  dat.tem[which.max(dat.tem$loglik),"target"];
      max.x.2 <- max.x;
      lik.95 <-  max(dat.tem$loglik) - 0.5*qchisq(p = p.value,df = 1)
      x.l.value <- seq(range(dat.tem[,"target"])[1],max.x,seq.long);
      x.h.value <- seq(max.x,range(dat.tem[,"target"])[2],seq.long);
      n.num.1 <- which.min(abs(predict(fit,newdata = x.l.value) - lik.95));
      n.num.2 <- which.min(abs(predict(fit,newdata = x.h.value) - lik.95));
      x.value.1 <- x.l.value[n.num.1]; x.value.2 <- x.h.value[n.num.2]; 
    }
    }
    result <- data.frame(
      type = id,loglik = max(dat.tem$loglik),mle = round(max.x,2),li = round(x.value.1,2),hi = round(x.value.2,2))
    return(result)
  })
  result <- bind_rows(result.list);return(result)
}


theshold.plot <- function(dataset,figure.file,mcap.t = T) {
  dataset <- mutate(dataset, id = substr(type,1,6), 
                    number = as.numeric(stri_extract_first(type,regex = "\\d+")))
  type.id <- unique(dataset$id)
  for (i in (1:length(type.id))) {
    dat.tem <- filter(dataset, id == type.id[i])
    number.name <- sort(as.numeric(unique(dat.tem$number)))

    plot.total <- list()
    for (j in (1:length(number.name))){
      dat.tem.1 <- filter(dat.tem, number == number.name[j])
      dat.tem.1 %>%
        plyr::ddply(~target,subset,loglik == max(loglik)) -> dat.tem.1
      dat.tem.1 <- mutate(dat.tem.1,LL = (loglik + 2*loglik.se),UL = (loglik - 2*loglik.se))
      fit <- locfit.robust(loglik~target,data = dat.tem.1,family = "qgauss",alpha = 0.9);
      max.x  <-  dat.tem.1[which.max(dat.tem.1$loglik),"target"];
      
      if(isTRUE(mcap.t)){
        dat.ci <- mcap(dat.tem.1$loglik,parameter = dat.tem.1$target)$ci
        dat.plot <- mcap(dat.tem.1$loglik,parameter = dat.tem.1$target)$fit
        names(dat.plot) <- c("x","y","z")
        x.value.1 <- dat.ci[1];  x.value.2 <- dat.ci[2];
        lik.95 <-  max(dat.plot$y) - 0.5 * qchisq(p = 0.95,df = 1)
      }else{
      lik.95 <-  max(dat.tem.1$loglik) - 0.5*qchisq(p = 0.95,df = 1)
      x.l.value <- seq(range(dat.tem.1$target)[1],max.x,0.01);
      x.h.value <- seq(max.x,range(dat.tem.1$target)[2],0.01);
      n.num.1 <- which.min(abs(predict(fit,newdata = x.l.value) - lik.95));
      n.num.2 <- which.min(abs(predict(fit,newdata = x.h.value) - lik.95));
      x.value.1 <- x.l.value[n.num.1]; x.value.2 <- x.h.value[n.num.2]; 
      dat.plot <- data.frame(y = predict(fit,dat.tem.1$target),x = dat.tem.1$target)
      }
      if (stri_detect_fixed(dat.tem.1[1,"type"],"AH")) {
      #  xlab.label <- xlab(expression(paste("AH on TR: ",eta[AH],sep = "")))
        xlab.label <- xlab(expression(paste("Absolute humidity on transmission rate (",eta[AH],")",sep = "")))
      } else if(stri_detect_fixed(dat.tem.1[1,"type"],"MT-FOI")){
      #  xlab.label <- xlab(expression(paste("MT on TR: ",eta[MT],sep = "")))
        xlab.label <- xlab(expression(paste("Mean temperature on transmission rate (",eta[MT],")",sep = "")))
      } else if(stri_detect_fixed(dat.tem.1[1,"type"],"MT-RR")){
      #  xlab.label <- xlab(expression(paste("MT on RR: ",psi[MT],sep = "")))
        xlab.label <- xlab(expression(paste("Mean temperature on reporting rate (",psi[MT],")",sep = "")))
      } else {
        xlab.label <- NULL
      }
      
      plot.total[[j]] <- ggplot(data = dat.plot,aes(x = x,y = y)) + geom_line(linetype = 1,alpha = 0.5) + 
        theme_bw(base_family = "serif",base_size = 7) + 
        geom_pointrange(data = dat.tem.1,aes(x = target,y = loglik,ymin = LL,ymax = UL),colour=cols[3],
                        fill = "yellow",size = 0.2,linetype = 1) +
        geom_hline(yintercept = lik.95,colour = "black",linetype = 2,alpha = 0.5) +
        geom_vline(xintercept = c(x.value.1,x.value.2,max.x),colour = cols[2],linetype = c(4,4,1)) +
        labs(title = paste("MLE:",round(max.x,2),"(",round(x.value.1,2),",",round(x.value.2,2),")",sep = "")) + 
        xlab.label + ylab("Profile log likehood")
    }
    plot.tem <- plot_grid(plotlist = plot.total,labels = number.name, ncol = 3)
    ggsave(paste(figure.file,"Detials/",unique(dataset$id)[i],".pdf",sep = ""),plot.tem,
           dpi = 300,width = 18,height = 18,units = "cm")
  }}



theshold.plot.2 <- function(dataset,figure.file,ncol = 3,mcap.t = T) {
    dat.tem <- dataset
    number.name <- sort(unique(dat.tem$type))
    
    plot.total <- list()
    for (j in (1:length(number.name))){
      dat.tem.1 <- filter(dat.tem, type == number.name[j])
      dat.tem.1 %>%
        plyr::ddply(~target,subset,loglik == max(loglik)) -> dat.tem.1
      dat.tem.1 <- mutate(dat.tem.1,LL = (loglik + 2*loglik.se),UL = (loglik - 2*loglik.se))
      fit <- locfit.robust(loglik~target,data = dat.tem.1,family = "qgauss",alpha = 0.9);
      max.x  <-  dat.tem.1[which.max(dat.tem.1$loglik),"target"];

      if(isTRUE(mcap.t)){
        dat.ci <- mcap(dat.tem.1$loglik,parameter = dat.tem.1$target)$ci
        dat.plot <- mcap(dat.tem.1$loglik,parameter = dat.tem.1$target)$fit
        names(dat.plot) <- c("x","y","z")
        x.value.1 <- dat.ci[1];  x.value.2 <- dat.ci[2];
        lik.95 <-  max(dat.plot$y) - 0.5 * qchisq(p = 0.95,df = 1)
      }else{
        lik.95 <-  max(dat.tem.1$loglik) - 0.5*qchisq(p = 0.95,df = 1)
        x.l.value <- seq(range(dat.tem.1$target)[1],max.x,0.01);
        x.h.value <- seq(max.x,range(dat.tem.1$target)[2],0.01);
        n.num.1 <- which.min(abs(predict(fit,newdata = x.l.value) - lik.95));
        n.num.2 <- which.min(abs(predict(fit,newdata = x.h.value) - lik.95));
        x.value.1 <- x.l.value[n.num.1]; x.value.2 <- x.h.value[n.num.2]; 
        dat.plot <- data.frame(y = predict(fit,dat.tem.1$target),x = dat.tem.1$target)
      }
       # xlab.label <- xlab(expression(paste("MT on RR: ",psi[MT],sep = "")))
      xlab.label <- xlab(expression(paste("Mean temperature on reporting rate (",psi[MT],")",sep = "")))
      
      plot.total[[j]] <- ggplot(data = dat.plot,aes(x = x,y = y)) + geom_line(linetype = 1,alpha = 0.5) + 
        theme_bw(base_family = "serif",base_size = 7) + 
        geom_pointrange(data = dat.tem.1,aes(x = target,y = loglik,ymin = LL,ymax = UL),colour=cols[3],
                        fill = "yellow",size = 0.2,linetype = 1) +
        geom_hline(yintercept = lik.95,colour = "black",linetype = 2,alpha = 0.5) +
        geom_vline(xintercept = c(x.value.1,x.value.2,max.x),colour = cols[2],linetype = c(4,4,1)) +
        labs(title = paste("MLE:",round(max.x,2),"(",round(x.value.1,2),",",round(x.value.2,2),")",sep = "")) + 
        xlab.label + ylab("Profile log likehood")
    }
    plot.tem <- plot_grid(plotlist = plot.total,labels = number.name, ncol = ncol)
    ggsave(paste(figure.file,"Detials/","MT-RR",".pdf",sep = ""),plot.tem,
           dpi = 300,width = 16,height = 18,units = "cm")
}


threshold.caculate <- function(dir,parameter){
  
  ##  dir: the list for all the csv loaction
  ##  parameter : the target parameter to estimate
  all.files <- list.files(path = dir,full.names = T,pattern = ".csv");
  if (parameter == "bah"){
    all.files <- all.files[stri_detect_fixed(all.files,"AH")]
  } else{
    all.files <- all.files[stri_detect_fixed(all.files,"MT")]
  }
  
  dat.all <- lapply(all.files,function(data) {
    dat.tem <- read.csv(data, header = T,stringsAsFactors = F);
    name <- gsub(dir,"",data);name <- gsub(".csv","",name);name <- gsub("/","",name);
    name <- gsub("-N","",name)
    dat.tem$type <- rep(name,dim(dat.tem)[1]);
    dat.tem %>% 
      subset(is.finite(loglik) & (nfail.max < 2),-c(nfail.max,nfail.min)) -> dat.tem
    dat.tem$target <- round(dat.tem[,parameter],3)
    dat.tem %>%
      plyr::ddply(~target,subset,loglik == max(loglik)) -> dat.pred
    return(dat.pred)
  })
  dat.temp <- bind_rows(dat.all); return(dat.temp)
} 



Paper.Figure.Threshold <- function(dataset,v.name,y.lim = c(0.4,1.6),y.by = 0.1,
                                   figure.name,plot = TRUE){
  
  dat.tem <- dataset;
  x.break <- sort(as.numeric(dat.tem$type))
  
  if(v.name == "bmt") {
  #  x.lab <- xlab(expression(paste("Threshold for MT on TR (",degree,"C)",sep = "")))
    x.lab <- xlab(expression(paste("Threshold for mean temperature on transmission rate (",degree,"C)",sep = "")))
    label_2 <- expression(paste("MLE for ",eta[MT]," and 95%CI",sep = ""));
  } else if(v.name == "bah"){
  #  x.lab <- xlab(expression(paste("Threshold for AH on TR (g/",m^3,")",sep = "")
    x.lab <- xlab(expression(paste("Threshold for absolute humidity on transmission rate (g/",m^3,")",sep = "")                              ))
    label_2 <- expression(paste("MLE for ",eta[AH]," and 95%CI",sep = ""));
  } else if(v.name == "gmt"){
   # x.lab <- xlab(expression(paste("Threshold for MT on RR (",degree,"C)",sep = "")))
    x.lab <- xlab(expression(paste("Threshold for mean temperature on reporting rate (",degree,"C)",sep = "")))
    label_2 <- expression(paste("MLE for ",psi[MT]," and 95%CI",sep = ""))
  } else {
    x.lab <- NULL; label_2 <- NULL
  }
  
  graph_2   <-  ggplot() +
    geom_pointrange(data = dat.tem, 
                    aes(x = type, y = mle, ymin = li, ymax = hi), color = "#925E9FB2") +
    geom_smooth(data = dat.tem,aes(x = type, y = mle)) +  ylab("") + 
    geom_hline(yintercept = 0,linetype = 2) + 
    scale_x_continuous(breaks = x.break) + 
    scale_y_continuous(breaks = round(seq(y.lim[1],y.lim[2],by = y.by),2),limits = y.lim)+
    x.lab +
    theme_bw(base_family = "serif",base_size = 15) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(margin = margin(0,2,0,0)))
  
  graph_1	 <- ggplot() +
    geom_line(data = dat.tem,aes(x = type, y = loglik ), colour = "#E8924299",size = 1) +
    scale_x_continuous(breaks = x.break) + ylab("") +
    theme_bw(base_family = "serif",base_size = 15) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.ticks.length = unit(0.1,"cm"),
          axis.text.y = element_text(margin = margin(0,0,0,4)))
  
  label_1 <- "Profile log likelihood"
  
  plot.tem <- dual_axis_graph(graph_2, graph_1, label_2,label_1)
  
  if(isTRUE(plot)){
    grid.draw(plot.tem)
  }
  ggsave(file = figure.name,plot.tem,dpi = 300, width = 18, height = 15, units = "cm")
  return(plot.tem)
} 


basic_plot <- function(dataset){
  ggplot(dataset, mapping = aes(x = target,y = loglik)) +
    geom_point(color = "red",size = 1.5) +
    facet_wrap(~type,scales = "free") + theme_bw(base_size = 10) +
    theme_bw() + xlab("Temperature effect on reporting probability") + ylab("")
}





theshold.plot.single <- function(dataset,figure.file,mcap.t = T) {

    dat.tem.1 <- dataset;
    dat.tem.1 %>%
        plyr::ddply(~target,subset,loglik == max(loglik)) -> dat.tem.1
    dat.tem.1 <- mutate(dat.tem.1,LL = (loglik + 2*loglik.se),UL = (loglik - 2*loglik.se))
    fit <- locfit.robust(loglik~target,data = dat.tem.1,family = "qgauss",alpha = 0.90);
    max.x  <-  dat.tem.1[which.max(dat.tem.1$loglik),"target"];
      
      if(isTRUE(mcap.t)){
        dat.ci <- mcap(dat.tem.1$loglik,parameter = dat.tem.1$target)$ci
        dat.plot <- mcap(dat.tem.1$loglik,parameter = dat.tem.1$target)$fit
        names(dat.plot) <- c("x","y","z")
        x.value.1 <- dat.ci[1];  x.value.2 <- dat.ci[2];
        lik.95 <-  max(dat.plot$y) - 0.5 * qchisq(p = 0.95,df = 1)
      }else{
      lik.95 <-  max(dat.tem.1$loglik) - 0.5*qchisq(p = 0.95,df = 1)
      x.l.value <- seq(range(dat.tem.1$target)[1],max.x,0.01);
      x.h.value <- seq(max.x,range(dat.tem.1$target)[2],0.01);
      n.num.1 <- which.min(abs(predict(fit,newdata = x.l.value) - lik.95));
      n.num.2 <- which.min(abs(predict(fit,newdata = x.h.value) - lik.95));
      x.value.1 <- x.l.value[n.num.1]; x.value.2 <- x.h.value[n.num.2]; 
      dat.plot <- data.frame(y = predict(fit,dat.tem.1$target),x = dat.tem.1$target)
      }
     # xlab.label <-  xlab(expression(paste("MT on RR: ",psi[MT],sep = "")))
    xlab.label <-  xlab(expression(paste("Mean temperature on reporting rate (",psi[MT],")",sep = ""))) 
      
  plot.tem <- ggplot(data = dat.plot,aes(x = x,y = y)) + geom_line(linetype = 1,alpha = 0.5) + 
        theme_bw(base_family = "serif",base_size = 9) + 
        geom_pointrange(data = dat.tem.1,aes(x = target,y = loglik,ymin = LL,ymax = UL),colour=cols[3],
                        fill = "yellow",size = 0.2,linetype = 1) +
        geom_hline(yintercept = lik.95,colour = "black",linetype = 2,alpha = 0.5) +
        geom_vline(xintercept = c(x.value.1,x.value.2,max.x),colour = cols[2],linetype = c(4,4,1)) +
        labs(title = paste("MLE: ",round(max.x,2)," (",round(x.value.1,2),", ",signif(x.value.2,2),")",sep = "")) + 
        xlab.label + ylab("Profile log likehood")
  
    ggsave(paste(figure.file,"Detials/","target.pdf",sep = ""),plot.tem,
           dpi = 300,width = 12,height = 12,units = "cm")
    return(plot.tem)
  }


stri_insert <- function(data,number,symbol){
  dat.1 <- substr(data,1,number)
  dat.2 <- substr(data,(number + 1), length(data))
  paste(dat.1,symbol,dat.2,sep = "")
}

num_pad_right <- function(data,num,symbol){
                    ifelse(stri_detect_fixed(data,"-"),
                           stri_pad_right(data,num + 1,symbol),
                           stri_pad_right(data,num,symbol))
}

p.value.test <- function(low,hi){
  dat.1 <- stri_detect_fixed(low,"-") & (stri_detect_fixed(hi,"-"))
  dat.2 <- (!stri_detect_fixed(low,"-")) & (!stri_detect_fixed(hi,"-"))
  ifelse(dat.1 | dat.2,1,0)
}




threshold.caculate.gmt <- function(dir,stri){
  
  ##  dir: the list for all the csv loaction
  ##  parameter : the target parameter to estimate
  all.files <- list.files(path = dir,full.names = T,pattern = ".csv");
  
  dat.all <- lapply(all.files,function(data) {
    dat.tem <- read.csv(data, header = T,stringsAsFactors = F);
    name <- gsub(dir,"",data);name <- gsub(".csv","",name);name <- gsub("/","",name);
    #name <- gsub("hangzhou_MT3ths","",name)
    name <- gsub(stri,"",name)
    dat.tem$type <- rep(name,dim(dat.tem)[1]);
    dat.tem %>% 
      subset(is.finite(loglik) & (nfail.max < 2),-c(nfail.max,nfail.min)) -> dat.tem
    dat.tem$target <- round(dat.tem[,"gmt"],3)
    dat.tem %>%
      plyr::ddply(~target,subset,loglik == max(loglik)) -> dat.pred
    return(dat.pred)
  })
  dat.temp <- bind_rows(dat.all); return(dat.temp)
} 
