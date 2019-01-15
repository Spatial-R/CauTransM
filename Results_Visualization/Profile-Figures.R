################################################################################################
#########################################  local profile #######################################
################################################################################################

library(ggplot2);library(reshape);library(magrittr);library(dplyr);
library(stringi);library(scales);library(locfit);library(cowplot);library(ggsci)


source("Results-Visualization/Profile-Functions.R");figure.dir <- "Figures/Profiles/"


################################################################################################
#########################################  full model  #########################################
################################################################################################

"Hangzhou"

Profile.Full.HZ <- read.csv("Results/Hangzhou/FullProfile/Changed/profile-1-20170303.csv",
                    header = T, stringsAsFactors = F)

#dat.all.1 <-  Omit.extreate.value(dat.all,range.value.1 = 3, range.value.2 = 2)
# dat.all <- dat.all.1[!dat.all.1$normal %in% c(1),]

for (i in c(1:nrow(Profile.Full.HZ))){
  parameter.name <- Profile.Full.HZ[i,"parm"]
  Profile.Full.HZ[i,"target"] <- Profile.Full.HZ[i,parameter.name]
}
type.name <- unique(Profile.Full.HZ$type)
Profile.Full.HZ$target <- round(log(Profile.Full.HZ$target),3)
Profile.Full.HZ.Line <- filter(Profile.Full.HZ,type == "line" & loglik > (-575))
Profile.Full.HZ.Line.2 <- Profile.calculate(dir = "Results/Hangzhou/FullProfile/Changed/Need")
Profile.Full.HZ.Line.2$type <- "line"
Profile.Full.HZ.Line.2$target <- log(Profile.Full.HZ.Line.2$target)
Profile.Full.HZ.Line <- rbind(Profile.Full.HZ.Line[!Profile.Full.HZ.Line$parm %in% c("kmt","krh"),],
                              Profile.Full.HZ.Line.2)

Profile.Full.Plot(dataset = Profile.Full.HZ,
                  profile.figure.dir = paste(figure.dir,"Full/Hangzhou/",sep = ""),
                  smooth.type = F,alpha = 0.8, p.value = 0.95, seq.long = 0.01)

Profile.Full.Plot(dataset = Profile.Full.HZ.Line,
                  profile.figure.dir = paste(figure.dir,"Full/Hangzhou/",sep = ""),
                  smooth.type = F,alpha = 0.8, p.value = 0.95, seq.long = 0.01)


"Quzhou"


Profile.Full.QZ <- Profile.calculate(dir = "Results/Quzhou/Profile-Line")

type.name <- unique(Profile.Full.QZ$type)
Profile.Full.QZ$target <- round(log(Profile.Full.QZ$target),3)
Profile.Full.Plot(dataset = Profile.Full.QZ,profile.figure.dir = paste(figure.dir,"Full/Quzhou/",sep = ""),
                  smooth.type = F,alpha = 0.8,omit = F,
                  p.value = 0.95, seq.long = 0.01)


################################################################################################
#########################################  easy model  #########################################
################################################################################################

"Hangzhou-Easy"

dir.1 <- "Results/Hangzhou/Easy"; dat.all.1 <- Profile.calculate(dir.1)
dat.all.1 <- mutate(dat.all.1, type = gsub("N","", type))

dir.2 <- "Results/Hangzhou/Easy-3"; dat.all.2 <- Profile.calculate(dir.2)

dir.3 <- "Results/Hangzhou/Easy-complication-2";dat.all.3 <- Profile.calculate(dir.3)

dir.4 <- "Results/Hangzhou/Easy-2";dat.all.4 <- Profile.calculate(dir.4)



Profile.Easy.HZ <- bind_rows(dat.all.1,dat.all.2,dat.all.3,dat.all.4)

Profile.Easy.HZ$type <- ifelse(Profile.Easy.HZ$type == "thresold1","easythresold1",
                               Profile.Easy.HZ$type)

for (i in c(1:nrow(Profile.Easy.HZ))){
  parameter.name <- Profile.Easy.HZ[i,"parm"]
  Profile.Easy.HZ[i,"target"] <- Profile.Easy.HZ[i,parameter.name]
}

Profile.Easy.HZ$target <- round(log(Profile.Easy.HZ$target),3)
omit.true.1 <- !((Profile.Easy.HZ$parm == "bmt") & (Profile.Easy.HZ$target < -0.2))
omit.true.2 <- !((Profile.Easy.HZ$parm == "grh") & (Profile.Easy.HZ$loglik < -580))
Profile.Easy.HZ <- Profile.Easy.HZ[omit.true.1 & omit.true.2,]

Profile.Easy.Plot(dataset = Profile.Easy.HZ,mcap.t = T,
                  smooth.type = F,profile.figure.dir = paste(figure.dir,"Easy/Hangzhou/",sep = ""), 
                  alpha = 0.9, p.value = 0.95, seq.long = 0.01)

"Quzhou-Easy"

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
Profile.Easy.QZ$target <- round(log(Profile.Easy.QZ$target),3)

loglik.hz <- unlist(lapply(unique(Profile.Easy.HZ$type),function(data) {
  dat.tem <- filter(Profile.Easy.HZ,type == data)
  max(dat.tem$loglik)
}))
data.frame(type = unique(Profile.Easy.HZ$type),loglik = loglik.hz)

loglik.qz <- unlist(lapply(unique(Profile.Easy.QZ$type),function(data) {
  dat.tem <- filter(Profile.Easy.QZ,type == data)
  max(dat.tem$loglik)
}))
data.frame(type = unique(Profile.Easy.QZ$type),loglik = loglik.qz)


#Profile.Easy.QZ %>%
#  plyr::ddply(~target + type ,subset,loglik == max(loglik)) -> Profile.Easy.QZ


Profile.Easy.Plot(dataset = Profile.Easy.QZ,mcap.t = T,
                  smooth.type = F,profile.figure.dir = paste(figure.dir,"Easy/Quzhou/",sep = ""), 
                  alpha = 0.9, p.value = 0.95, seq.long = 0.01)

##################################  Figures in Paper  #################################

Profile.Easy.QZ$city <- rep("Quzhou",nrow(Profile.Easy.QZ))
Profile.Easy.HZ$city <- rep("Hangzhou",nrow(Profile.Easy.HZ))
Profile.Easy <- rbind(Profile.Easy.HZ,Profile.Easy.QZ)

result.list <- lapply(c("Hangzhou","Quzhou"),function(id){
  result <- ci.calculate(dataset = Profile.Easy[Profile.Easy$city == id,],mcap.t = T,
                         p.value = 0.95, alpha = 0.9,seq.long = 0.01,smooth.type = F)
  
  names(result) <- c("max.x" ,"max.x.2","x.value.1","x.value.2","type","parm")
  result$type <- gsub("easy","",result$type)
  result$type <- factor(result$type,levels = c("line","thresold1","thresold2","thresold3","lag","lagthresold1","lagthresold2","lagthresold3"),
                        labels = c("Linear","T-MT-FOI","T-AH-FOI","T-MT-RR",
                                   "Lag","L-T-MT-FOI","L-T-AH-FOI","L-T-MT-RR"))
  result$city <- rep(id,nrow(result))
  return(result)
})
result.profile.ci <- bind_rows(result.list);        

result.profile.ci$parm <-  factor(result.profile.ci$parm,levels = c("bah","bmt","gmt","grh"))

figure.dir <- "Figures/Profiles"
Point_CI_Plot(result.profile.ci,figure.dir = "Figures/Profiles",city = T)


######################################### Changed  #########################################


result.profile.ci$lg <- stri_detect_fixed(result.profile.ci$type,"L-")
