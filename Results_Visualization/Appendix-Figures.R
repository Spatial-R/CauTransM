#######################################################################################################
############################################## Basic Description ######################################
#######################################################################################################

library(dplyr);library(zoo);library(humidity);library(readr)
library(reshape);library(ggplot2);library(lubridate)
library(cowplot)
source("Results-Visualization/Appendix-Functions.R")
city.ch <- c("杭州市","嘉兴市","宁波市","湖州市", "金华市","绍兴市","丽水市","衢州市","舟山市","台州市","温州市")
source("Results-Visualization/Basic.R")
figure.app <- "Figures/Appendix"
source("Results-Visualization/Profile-Functions.R");

#####################################################################################################

mump.case <- case_manupilation(dir.1 = "Data/Cases/dattran2.csv",dir.2 = "Data/Cases/2015.csv",
                               city.name = c("Hangzhou","Quzhou"))
mump.case$date <- as.Date(mump.case$date); mump.case <- arrange(mump.case,date);
dat.all <- data.frame(summarise(group_by(mump.case,city),count=sum(cases)))
dat.all[1,2]/9.02/100/4; dat.all[2,2]/2.54/100/4;

mump.case <- mutate(mump.case,week = floor_date(date, unit = "week"))
mump.case <- data.frame(summarise(group_by(mump.case,week,city),cases =sum(cases)))
mump.case <- filter(mump.case,!substr(week,1,4) == "2015")

dir.tem <- "Data/Meterological"
tem.air.hz <- climatic.data(dir.tem,pattern.type = "txt",site.code = "58457")
tem.air.hz$city <- rep("Hangzhou",nrow(tem.air.hz))
tem.air.qz <- climatic.data(dir.tem,pattern.type = "txt",site.code = "58633")
tem.air.qz$city <- rep("Quzhou",nrow(tem.air.qz))
tem.air <- rbind(tem.air.hz,tem.air.qz)
tem.air$AH <- AH(WVP2(tem.air$RH,SVP(C2K(tem.air$MT/10))),C2K(tem.air$MT/10))*1000
tem.air <- mutate(tem.air,week = floor_date(date, unit = "week"))
tem.air <- data.frame(summarise(group_by(tem.air,week,city),MT = mean(MT,na.rm = T),
                                                            RH = mean(RH,na.rm = T),
                                                            AH = mean(AH,na.rm = T)))
tem.air <- filter(tem.air,substr(week,1,4) %in% c(2011:2014))

mump.case$rate <- ifelse(mump.case$city == "Hangzhou",mump.case$case/9368416,mump.case$case/2564898)

first.size <- 10

case.plot <- ggplot(mump.case,aes(x = week,y = rate,colour = factor(city))) +
  geom_line(size = 0.9,alpha = 1) +
#  scale_color_manual(values = c("Hangzhou" = "#008B45FF","Quzhou" = "#631879FF"),
  scale_color_manual(values = c("Hangzhou" = "red","Quzhou" = "green"),
                     guide = guide_legend(title = "")) +
  theme_classic(base_family = "serif",base_size = first.size) +
  xlab("") + ylab("Weekly incidence of Mumps") + 
  theme(legend.position = c(0.8,0.9))

ggsave(paste(figure.app,"/case.pdf",sep = ""),case.plot,dpi = 400, width = 15,height = 9,units = "cm")


figure.ah <- ggplot(tem.air,aes(x = week,y = AH,colour = factor(city))) +
  geom_line(size = 0.9,alpha = 1) +
  scale_color_manual(values = c("Hangzhou" = "red","Quzhou" = "green"),
                      guide = "none") +
  theme_bw(base_family = "serif",base_size = first.size) +
  xlab("") + ylab(expression(paste("Absolute humidity (g/",m^3,")",sep = "")
  ))

figure.mt <- ggplot(tem.air,aes(x = week,y = MT/10,colour = factor(city))) +
  geom_line(size = 0.9,alpha = 1) +
  scale_color_manual(values = c("Hangzhou" = "red","Quzhou" = "green"),
                    # guide = guide_legend(title = "")
                    guide = "none") +
  theme_bw(base_family = "serif",base_size = first.size) +
  xlab("") + ylab(expression(paste("Mean temperature (",degree,"C )",sep = "")
  ))

figure.rh <- ggplot(tem.air,aes(x = week,y = RH,colour = factor(city))) +
  geom_line(size = 0.9,alpha = 1) +
  scale_color_manual(values = c("Hangzhou" = "red","Quzhou" = "green"),
                     #guide = guide_legend(title = "")) +
                     guide = "none") +
  theme_bw(base_family = "serif",base_size = first.size) +
  xlab("") + ylab("Relative humidity (%)")

fig.cli <- plot_grid(figure.mt,figure.ah,figure.rh,ncol = 1)

ggsave(paste(figure.app,"/climatic.pdf",sep = ""),fig.cli,dpi = 400, width = 16,height = 15,units = "cm")

library(imager)
library(magick)


p2 <- ggdraw() + draw_image("E:\\Project\\Mumps\\POMP\\Codes\\Figures\\Flowmap\\flowchart5.tiff", 
                            scale = 1)

fig.1 <- plot_grid(case.plot,p2,ncol = 1,labels = c("B","C"))
fig.2 <- plot_grid(fig.cli,fig.1,ncol = 2,labels = c("A",""))
ggsave(paste(figure.app,"/all.pdf",sep = ""),fig.2,dpi = 400, width = 20,height = 15,units = "cm")


#################################################################################################
############################################## Paris Matrix #####################################
#################################################################################################

Global.Paramter.Paris(dir = "Results/Global",maxlogik = 10)
Global.Paramter.Paris(dir = "Results/Global",maxlogik = 10,climatic = F)

#################################################################################################
############################################## Prediction #######################################
#################################################################################################


"Hangzhou-Easy"

dir.1 <- "Results/Hangzhou/Easy"; dat.all.1 <- Profile.calculate(dir.1)
dat.all.1 <- mutate(dat.all.1, type = gsub("N","", type))
dir.2 <- "Results/Hangzhou/Easy-3"; dat.all.2 <- Profile.calculate(dir.2)
dir.3 <- "Results/Hangzhou/Easy-complication-2";dat.all.3 <- Profile.calculate(dir.3)
dir.4 <- "Results/Hangzhou/Easy-2";dat.all.4 <- Profile.calculate(dir.4)
Profile.Easy.HZ <- bind_rows(dat.all.1,dat.all.2,dat.all.3,dat.all.4)
Profile.Easy.HZ$type <- ifelse(Profile.Easy.HZ$type == "thresold1","easythresold1",
                               Profile.Easy.HZ$type)
loglik.max.hz <- plyr::ddply(Profile.Easy.HZ,~type, subset,loglik == max(loglik))
loglik.max.hz <- plyr::ddply(Profile.Easy.HZ[Profile.Easy.HZ$parm == "gmt",],
                             ~type, subset,loglik == max(loglik))

loglik.max.hz$type <- gsub("thresold","threshold",loglik.max.hz$type)
loglik.max.hz$type <- gsub("easy","",loglik.max.hz$type)

dir.profile.1 <- "Results/Quzhou/Profile-Easy/20170502"
dir.profile.2 <- "Results/Quzhou/Profile-Easy/20170505"
dir.profile.3 <- "Results/Quzhou/Profile-Easy/20170505-1"
dir.profile.4 <- "Results/Quzhou/Profile-Easy/20170508"
dir.profile.5 <- "Results/Quzhou/Profile-Easy/20170504"
dir.profile.6 <- "Results/Quzhou/Profile-Easy/20170601"

profile.5 <- Profile.calculate(dir.profile.5);
profile.1 <- Profile.calculate(dir.profile.1);
profile.2 <- Profile.calculate(dir.profile.2);
profile.3 <- Profile.calculate(dir.profile.3);
profile.4 <- Profile.calculate(dir.profile.4);
profile.6 <- Profile.calculate(dir.profile.6);

Profile.Easy.QZ <- rbind(profile.5, profile.1,profile.2,profile.3,profile.4,profile.6)
loglik.max.qz <- plyr::ddply(Profile.Easy.QZ,~type, subset,loglik == max(loglik))
loglik.max.qz$type <- gsub("thresold","threshold",loglik.max.qz$type)
loglik.max.qz$type <- gsub("easy","",loglik.max.qz$type)


identical(loglik.max.hz$type,loglik.max.qz$type)
load("Data/dataset-hangzhou.RData");load("Data/dataset-quzhou.RData")


for (i in c(1:8)){

theta.hz <- unlist(loglik.max.hz[i,1:21])
mumps.hz <- case.hz; covar <- covar.hz
source(paste("Models/Prediction/",loglik.max.qz[i,"type"],"full.R",sep = ""))  
fig.1 <- rmse_plot(theta = theta.hz, model = origin.model)

theta.qz <- unlist(loglik.max.qz[i,1:21])
mumps.hz <- case.qz; covar <- covar.qz
source(paste("Models/Prediction/",loglik.max.qz[i,"type"],"full.R",sep = ""))  
fig.2 <- rmse_plot(theta = theta.qz, model = origin.model,ylab.plot = F)
plot.final <- plot_grid(fig.1,fig.2,labels = c("Hangzhou","Quzhou"))

ggsave(paste(figure.app,"/Predictions/",loglik.max.hz[i,"type"],".pdf",sep = ""), plot.final,
       width = 18, height = 10,units = "cm", dpi = 300)
print(i)
}




