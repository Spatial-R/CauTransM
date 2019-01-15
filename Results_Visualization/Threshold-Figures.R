
library(ggplot2);library(reshape);library(magrittr);library(dplyr)
library(stringi);library("scales");library(locfit);library(cowplot)
library(ggsci)

figue.ah.dir  <- "Figures/Threshold-AH-FOI/"
figue.mt.dir  <- "Figures/Threshold-MT-FOI/"
figue.mt2.dir <- "Figures/Threshold-MT-RR/"

source("Results-Visualization/Threshold-Functions.R")
source("Results-Visualization/duaggplot2.R")

cols <- pal_lancet("lanonc", alpha = 0.7)(9)

mypal = pal_lancet("lanonc", alpha = 0.7)(9)
mypal = pal_rickandmorty("schwifty", alpha = 0.6)(9)
show_col(mypal)


####################################################################################################
#########################################      Hangzhou   ##########################################
####################################################################################################

dir <- "Results/Hangzhou/Easy-Theshold/O"

dat.all.bah <- threshold.caculate(dir,"bah")
dat.all.bah <- filter(dat.all.bah, target > 0.7)

ggplot(dat.all.bah, mapping=aes(x = log(target),y = loglik)) +
  geom_point(color = "red",size=1.5) +
  facet_wrap(~type,scales="free") + theme_bw(base_size = 10) +
  theme_bw() + xlab("Temperature effect on reporting probability") + ylab("")

##################################### log the parameter value #############################

dat.all.bah$target <- round(log(dat.all.bah$target),3)

dat.all.1 <- filter(dat.all.bah,substr(type,1,4) == "T-AH" & target > log(0.7))
dat.all.2 <- filter(dat.all.bah,substr(type,1,6) == "T-L-AH" & target > log(0.7))
dat.ah.threshold <- rbind(dat.all.1,dat.all.2)
dat.ah.threshold <- filter(dat.ah.threshold, !type %in% c("T-AH-FOI-2","T-AH-FOI-3","T-AH-FOI-16",
                                                          "T-L-AH-FOI-3","T-L-AH-FOI-2","T-L-AH-FOI-16"))
dat.theshold.ci <- ci.threshold(dataset = dat.ah.threshold,parameter = "bah",alpha = 0.9,
                                smooth.type = F,p.value = 0.95,seq.long = 0.01)
dat.theshold.ci$id <- substr(dat.theshold.ci$type,1,3)
dat.theshold.ci$type <- as.numeric(unlist(stri_match_all(dat.theshold.ci$type, regex = "[0-9]+")))

###################################  T-AH #########################################

figure.ah.t.hz  <- Paper.Figure.Threshold(dataset = dat.theshold.ci[dat.theshold.ci$id == "T-A",],
                       v.name = "bah",y.lim = c(-0.2,0.2),y.by = 0.1,
                       figure.name = paste(figue.ah.dir,"Hangzhou/hangzhou-T-AH.pdf",sep = ""),
                       plot = TRUE)

###################################  L-T-AH #########################################

figure.ah.lt.hz  <- Paper.Figure.Threshold(dataset = dat.theshold.ci[dat.theshold.ci$id == "T-L",],
                       v.name = "bah",y.lim = c(-0.2,0.2),y.by = 0.1,
                       figure.name = paste(figue.ah.dir,"Hangzhou/hangzhou-T-L-AH.pdf",sep = ""),
                       plot = TRUE)

plot.ah.final.hz <- plot_grid(figure.ah.t.hz,figure.ah.lt.hz,labels = c("A1","A2"),ncol = 2 )

ggsave(file = paste(figue.ah.dir,"Hangzhou/hangzhou-AH.pdf",sep = ""), plot.ah.final.hz,
       dpi = 300, width = 28, height = 17, units = "cm")

##########################################  T-MT ####################################

dat.all.bmt <- threshold.caculate(dir,"bmt")
dat.all.1 <- filter(dat.all.bmt,substr(type,1,4) == "T-MT" & target > 0.72)
dat.all.2 <- filter(dat.all.bmt,substr(type,1,6) == "T-L-MT" & target > 0.72)
dat.mt.threshold <- rbind(dat.all.1,dat.all.2)
dat.mt.threshold <- filter(dat.mt.threshold, !type %in% c("T-MT-FOI-1","T-MT-FOI-2",
                                                          "T-L-MT-FOI-1","T-L-MT-FOI-2"))
dat.mt.threshold$target <- round(log(dat.mt.threshold$target),3)

bmt.threshold.ci <- ci.threshold(dat.mt.threshold,parameter = "bmt",alpha = 0.9,
                                smooth.type = F,p.value = 0.95,seq.long = 0.01)
bmt.threshold.ci$id <- substr(bmt.threshold.ci$type,1,3)
bmt.threshold.ci$type <- as.numeric(unlist(stri_match_all(bmt.threshold.ci$type, regex = "[0-9]+")))


###################################  T- MT #########################################

figure.mt.t.hz <- Paper.Figure.Threshold(dataset = bmt.threshold.ci[bmt.threshold.ci$id == "T-M",],
                                       v.name = "bmt",y.lim = c(-0.2,0.2),y.by = 0.1,
                                       figure.name = paste(figue.mt.dir,"Hangzhou/hangzhou-T-MT.pdf",sep = ""),
                                       plot = TRUE)

###################################  L-T-MT  #########################################


figure.mt.lt.hz <- Paper.Figure.Threshold(dataset = bmt.threshold.ci[bmt.threshold.ci$id == "T-L",],
                                        v.name = "bmt",y.lim = c(-0.2,0.2),y.by = 0.1,
                                        figure.name = paste(figue.mt.dir,"Hangzhou/hangzhou-T-L-MT.pdf",sep = ""),
                                        plot = TRUE)

plot.mt.final.hz <- plot_grid(figure.mt.t.hz,figure.mt.lt.hz,labels = c("A1","A2"),ncol = 2 )

ggsave(file = paste(figue.mt.dir,"Hangzhou/hangzhou-MT.pdf",sep = ""), plot.mt.final.hz,
       dpi = 300, width = 28, height = 17, units = "cm")


threshold.hz.tr <- rbind(bmt.threshold.ci[bmt.threshold.ci$id == "T-M",],
                         dat.theshold.ci[dat.theshold.ci$id == "T-A",])

#################################################################################################
################################## Profile Deatials ############################################
#################################################################################################

theshold.plot(dat.mt.threshold,figure.file = paste(figue.mt.dir,"Hangzhou/",sep = ""))

dat.ah.threshold <- filter(dat.ah.threshold,!type == "T-L-AH-FOI-16" )
theshold.plot(dat.ah.threshold,figure.file = paste(figue.ah.dir,"Hangzhou/",sep = ""))

####################################################################################################
#########################################      Quzhou   ############################################
####################################################################################################

###########################################  AH on FOI  ###########################################

dir.ah.qz.1 <- "Results/Quzhou/Threshold-AH-FOI/20170508"
dir.ah.qz.2 <- "Results/Quzhou/Threshold-AH-FOI/Complicated"
dir.ah.qz.3 <- "Results/Quzhou/Threshold-AH-FOI/Complicated-1"


dat.bah.qz.1 <- threshold.caculate(dir.ah.qz.1,"bah")
dat.bah.qz.1$type <- gsub("-C","",dat.bah.qz.1$type)

dat.bah.qz.2 <- threshold.caculate(dir.ah.qz.2,"bah")
dat.bah.qz.2$type <- gsub("-C","",dat.bah.qz.2$type)

dat.bah.qz.3 <- threshold.caculate(dir.ah.qz.3,"bah")
dat.bah.qz.3$type <- gsub("-C","",dat.bah.qz.3$type)


dat.bah.qz <- rbind(dat.bah.qz.1,dat.bah.qz.2,dat.bah.qz.3)
basic_plot(dat.bah.qz)

dat.bah.qz <- filter(dat.bah.qz,
                     !type %in% c("T-AH-FOI-2","T-AH-FOI-3","T-AH-FOI-16",
                                  "L-T-AH-FOI-3","L-T-AH-FOI-2","L-T-AH-FOI-16"))
dat.bah.qz$target <- round(log(dat.bah.qz$target),3)
dat.bah.ci <- ci.threshold(dat.bah.qz,parameter = "bah",alpha = 0.9,
                                smooth.type = F,p.value = 0.95,seq.long = 0.01)

#dat.theshold.ci <- lapply(unique(dat.bah.qz$type),function(id){
#  dat.tem <- dat.bah.qz[dat.bah.qz$type == id,]
#  result <- mcap(dat.tem$loglik,dat.tem$target,Ngrid=100,lambda=0.8)
#  data.frame(mle = result$mle,li = result$ci[1],hi = result$ci[2],id = id)
#})
#dat.theshold.ci <- bind_rows(dat.theshold.ci)

dat.bah.ci$id <- substr(dat.bah.ci$type,1,3)
dat.bah.ci$type <- as.numeric(unlist(stri_match_all(dat.bah.ci$type, regex = "[0-9]+")))

###################################  T-AH     ##########################################

figure.ah.t.qz  <- Paper.Figure.Threshold(dataset = dat.bah.ci[dat.bah.ci$id == "T-A",],
                                       v.name = "bah",y.lim = c(-0.6,0.6),y.by = 0.2, 
                                       figure.name = paste(figue.ah.dir,"Quzhou/quzhou-T-AH.pdf",sep = ""),
                                       plot = TRUE)

###################################  L-T-AH   #########################################

figure.ah.lt.qz  <- Paper.Figure.Threshold(dataset = dat.bah.ci[dat.bah.ci$id == "L-T",],
                                        v.name = "bah",y.lim = c(-0.6,0.6),y.by = 0.2, 
                                        figure.name = paste(figue.ah.dir,"Quzhou/quzhou-T-L-AH.pdf",sep = ""),
                                        plot = TRUE)

plot.ah.final.qz <- plot_grid(figure.ah.t.qz,figure.ah.lt.qz,labels = c("B1","B2"),ncol = 2 )

ggsave(file = paste(figue.ah.dir,"Quzhou/quzhou-AH.pdf",sep = ""), plot.ah.final.qz,
       dpi = 300, width = 28, height = 17, units = "cm")

theshold.plot(dat.bah.qz,figure.file = paste(figue.ah.dir,"Quzhou/",sep = ""))


################################  MT   #############################################

dir.mt.qz <- paste("Results/Quzhou/Threshold-MT-FOI/Complicated-",1:4,sep = "")
  
dat.bmt.list <- lapply(dir.mt.qz,function(id){
  dat.tem <- threshold.caculate(id,"bmt")
  dat.tem$type <- gsub("-C","",dat.tem$type)
  return(dat.tem)
})
dat.bmt.qz <- bind_rows(dat.bmt.list)

omit.value <- !(dat.bmt.qz$type %in% c("L-T-MT-FOI-6","L-T-MT-FOI-4",
                                       "L-T-MT-FOI-7","L-T-MT-FOI-8") & (dat.bmt.qz$target > 1.16))
dat.bmt.qz <- dat.bmt.qz[omit.value,]

dat.bmt.qz$target <- round(log(dat.bmt.qz$target),3)

dat.bmt.qz.1 <- data.frame(summarise(group_by(dat.bmt.qz,type,target,bmt),loglik = mean(loglik)))

dat.bmt.qz.1 <- filter(dat.bmt.qz.1, !type %in% c("T-MT-FOI-2"))

basic_plot(dat.bmt.qz.1)

dat.bmt.ci <- ci.threshold(dat.bmt.qz.1,parameter = "bmt",alpha = 0.9,
                           smooth.type = F,p.value = 0.95,seq.long = 0.01)
dat.bmt.ci$id <- substr(dat.bmt.ci$type,1,3)
dat.bmt.ci$type <- as.numeric(unlist(stri_match_all(dat.bmt.ci$type, regex = "[0-9]+")))


#dat.theshold.ci <- lapply(unique(dat.bmt.qz.1$type),function(id){
#  dat.tem <- dat.bmt.qz.1[dat.bmt.qz.1$type == id,]
#  result <- mcap(dat.tem$loglik,dat.tem$target,Ngrid=100,lambda=0.8)
#  data.frame(mle = result$mle,li = result$ci[1],hi = result$ci[2],id = id)
#})
#dat.theshold.ci <- bind_rows(dat.theshold.ci)


###################################  T-MT    ##########################################

figure.mt.t.qz  <- Paper.Figure.Threshold(dataset = dat.bmt.ci[dat.bmt.ci$id == "T-M",],
                                       v.name = "bmt", y.lim = c(-0.6,0.6),y.by = 0.2,
                                       figure.name = paste(figue.mt.dir,"Quzhou/quzhou-T-MT.pdf",sep = ""),
                                       plot = TRUE)

###################################  L-T-MT   #########################################

figure.mt.lt.qz  <- Paper.Figure.Threshold(dataset = dat.bmt.ci[dat.bmt.ci$id == "L-T",],
                                        v.name = "bmt",y.lim = c(-0.6,0.6),y.by = 0.2,
                                        figure.name = paste(figue.mt.dir,"Quzhou/quzhou-T-L-MT.pdf",sep = ""),
                                        plot = TRUE)

plot.mt.final.qz <- plot_grid(figure.mt.t.qz,figure.mt.lt.qz,labels = c("B1","B2"),ncol = 2)

ggsave(file = paste(figue.mt.dir,"Quzhou/quzhou-MT.pdf",sep = ""), plot.mt.final.qz,
       dpi = 300, width = 28, height = 17, units = "cm")

theshold.plot(dat.bmt.qz,figure.file = paste(figue.mt.dir,"Quzhou/",sep = ""))

######################################  Merge the Hangzhou and Quzhou   ##########################

plot.mt.final <- plot_grid(plot.mt.final.hz,plot.mt.final.qz,labels = c("",""), ncol = 1)

ggsave(file = paste("Figures","/MT.pdf",sep = ""), plot.mt.final,
       dpi = 300, width = 28, height = 28, units = "cm")


plot.ah.final <- plot_grid(plot.ah.final.hz,plot.ah.final.qz,labels = c("",""), ncol = 1)

ggsave(file = paste("Figures","/AH.pdf",sep = ""), plot.ah.final,
       dpi = 300, width = 28, height = 28, units = "cm")

threshold.qz.tr <- rbind(dat.bmt.ci[dat.bmt.ci$id == "T-M",],
                         dat.bah.ci[dat.bah.ci$id == "T-A",])
threshold.qz.tr$city <- "QZ";threshold.hz.tr$city <- "HZ"
threshold.tr <- rbind(threshold.hz.tr,threshold.qz.tr)

###############################################################################################


figure_threshold <-  plot_grid(figure.ah.t.hz,figure.ah.t.qz,figure.mt.t.hz,figure.mt.t.qz,ncol = 2,
          labels = c("A(1)","A(2)","B(1)","B(2)"))
ggsave(file = paste("Figures","/TR-THRESHOLD.pdf",sep = ""), figure_threshold,
       dpi = 300, width = 32, height = 22, units = "cm")

###############################################################################################
################################### Threshold for temperature on RR  ##########################
###############################################################################################

dir.mt.qz.rr  <- "Results/Quzhou/Threshold-MT-RR/Basic"
dat.all.gmt   <- threshold.caculate(dir.mt.qz.rr,"gmt"); 

dir.mt.qz.rr.1  <- "Results/Quzhou/Threshold-MT-RR/Complicated"
dat.all.gmt.1   <- threshold.caculate(dir.mt.qz.rr.1,"gmt"); 

dat.all.gmt <- rbind(dat.all.gmt,dat.all.gmt.1)

dat.all.gmt$type <- as.numeric(unlist(stri_match_all(dat.all.gmt$type, regex = "[0-9]+"))) + 9

################## add thress threshold 7-9 C for MT on 

dir.mt.qz.rr.add  <- "Results/Quzhou/Threshold-MT-RR/Complicated-1"
dat.all.gmt.add   <- threshold.caculate(dir.mt.qz.rr.add,"gmt"); 
dat.all.gmt.add$type <- as.numeric(unlist(stri_match_all(dat.all.gmt.add$type, regex = "[0-9]+"))) + 6

dat.all.gmt <- rbind(dat.all.gmt,dat.all.gmt.add)

omit.true.19 <- !(dat.all.gmt$type %in% c(19:24) & (dat.all.gmt$target < 0.2 | dat.all.gmt$target > 1.0))
omit.true.18 <- !(dat.all.gmt$type %in% c(15:18) & (dat.all.gmt$target < 0.2 | dat.all.gmt$target > 2.0))
omit.true.10 <- !(dat.all.gmt$type %in% c(13:14) & (dat.all.gmt$target < 0.2))
omit.true.79 <- !(dat.all.gmt$type %in% c(7:9) & (dat.all.gmt$target < 1))

dat.all.gmt <- dat.all.gmt[omit.true.19 & omit.true.10 & omit.true.18 & omit.true.79,]
dat.all.gmt$target <- round(log(dat.all.gmt$target),3)
basic_plot(dat.all.gmt)


######################################################################################################
############################################# plot the dataset #######################################
######################################################################################################

dat.gmt.ci <- ci.threshold(dat.all.gmt,parameter = "gmt",alpha = 0.9,
                           smooth.type = F,p.value = 0.95,seq.long = 0.01)

figure.mt.lt  <- Paper.Figure.Threshold(dataset = dat.gmt.ci,
                                        v.name = "gmt",y.lim = c(-1.2,1.6),y.by = 0.4,
                                        figure.name = paste(figue.mt2.dir,"Quzhou/quzhou-MT-RR.pdf",sep = ""),
                                        plot = TRUE)

theshold.plot.2(dat.all.gmt,figure.file = paste(figue.mt2.dir,"Quzhou/",sep = ""))


#######################################################################################################
############################################### MT-RR: step two #######################################
#######################################################################################################


dir.gmt.ushape  <- "Results/Quzhou/Threshold-MT-RR/U-shape-2"  
dat.gmt.ushape   <- threshold.caculate(dir.gmt.ushape,"gmt");basic_plot(dat.gmt.ushape)
dat.gmt.ushape$type <- as.numeric(unlist(stri_match_all(dat.gmt.ushape$type, regex = "[0-9]+"))) 
dat.gmt.ushape$type <- stri_pad_left(dat.gmt.ushape$type, 4, pad = "0")
dat.gmt.ushape$type <- stri_insert(dat.gmt.ushape$type,2,"-")

dat.gmt.ushape$target <- round(log(dat.gmt.ushape$target),3)
dat.gmt.ushape.ci <- ci.threshold(dat.gmt.ushape,parameter = "gmt",alpha = 0.9,
                           smooth.type = F,p.value = 0.95,seq.long = 0.01)

gmt.qz.u <- filter(dat.gmt.ushape.ci,type == "07-15");gmt.qz.u$city <- "QZ" 

dat.gmt.target <- filter(dat.gmt.ushape,type == "07-15")
plot.target <- theshold.plot.single(dat.gmt.target,figure.file = paste(figue.mt2.dir,"Quzhou/",sep = ""))


basic.plot <- ggplot() +
  geom_pointrange(data = dat.gmt.ci, 
                  aes(x = type, y = mle, ymin = li, ymax = hi), color = "#925E9FB2") +
  geom_smooth(data = dat.gmt.ci,aes(x = type, y = mle)) +  ylab("") + 
  geom_hline(yintercept = 0,linetype = 2) + 
  scale_x_continuous(breaks = dat.gmt.ci$type) + 
  scale_y_continuous(breaks = round(seq(-1.2,1.6,by = 0.4),2))+
  xlab(expression(paste("Threshold for mean temperature on reporting rate (",degree,"C )",sep = ""))) +
  ylab(expression(paste(psi[MT],sep = ""))) +
  theme_bw(base_family = "serif",base_size = 15) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(margin = margin(0,2,0,0)))

basic.plot + annotation_custom(ggplotGrob(plot.target),
                               xmin = 15, xmax = 23, ymin = 0.1, ymax = 1.5)

ggsave(filename = "Figures/Threshold-MT-RR/Quzhou/MT-RR-Total.pdf",dpi = 300,
       width = 16 ,height = 14,units = "cm")


#################################################################################################
##################################### U-shape to test the result ################################
#################################################################################################


theshold.plot.2(dat.gmt.ushape,figure.file = "Figures/Threshold-MT-RR/Quzhou/U-shape/",ncol = 3)

dat.gmt.ushape.ci[,c(3:5)] <- apply(dat.gmt.ushape.ci[,c(3:5)],2,function(data){
  num_pad_right(data,num = 4, symbol = "0")
})
dat.gmt.ushape.ci <- mutate(dat.gmt.ushape.ci, type1 = substr(type,1,2),
                                               type2 = substr(type,4,5),
                                               p.value = p.value.test(li,hi),
                                               or.text = paste("(",li,",",hi,")",sep = ""))
dat.gmt.ushape.ci[,c(6,7)] <- apply(dat.gmt.ushape.ci[,c(6,7)],2,as.numeric)

ggplot(dat.gmt.ushape.ci, aes(x = type1, y = type2)) + 
  geom_point(aes(size = as.numeric(mle), color = factor(p.value))) + 
  scale_size_continuous(range = c(0.2,1.0)*10,
                        breaks = c(seq(-0.3,0.5,0.2)),
                        labels = seq(-0.3,0.5,0.2),
                        guide = guide_legend(title = expression(psi[MT]),keywidth = 2,
                                             override.aes = list(colour = "grey"))) +
  geom_text(aes(label = mle,parse = T,nudge_y = -0.1),size = 4) +
  geom_text(data = dat.gmt.ushape.ci,aes(x = type1, y = type2 - 0.1,label = or.text),size = 4) +
  theme_bw(base_family = "serif",base_size = 12) +
  scale_colour_manual(values = c("grey","red"),labels = c("Insiginificant","Siginificant"),
                      guide = guide_legend(title = "Statistical Level",keyheight = 1.3,override.aes = list(size=5)))+
  theme(legend.key.size = unit(1.2,"cm")) + 
  scale_x_continuous(limits = c(6.5,10.5),breaks = c(7:10),labels = c(7:10)) +
  scale_y_continuous(limits = c(14.8,17.2),breaks = c(15:17),labels =c(15:17)) +
  xlab(expression(paste("Cold threshold (",degree,"C)",sep = ""))) + 
  ylab(expression(paste("Hot threshold (",degree,"C)",sep = "")))

ggsave(filename = "Figures/Threshold-MT-RR/Quzhou/combination.pdf",
       width = 16,height = 14,dpi = 400,units = "cm")


#########################################################################################################
############################################## U-shaped for Hangzhou  ###################################
#########################################################################################################


########### V-shaped

gmt.result <- threshold.caculate.gmt(
  dir = "E:/Project/Mumps/POMP/Codes/Analysis-Codes/Threshold-MT-RR/Results",
  stri = "hangzhou_MT3ths")
gmt.result$target <- round(log(gmt.result$target),3)
gmt.result$type <- as.numeric(gmt.result$type)
theshold.plot.2(gmt.result,figure.file = "Figures/Threshold-MT-RR/Hangzhou/",ncol = 3)

dat.gmt.ci <- ci.threshold(gmt.result,parameter = "gmt",alpha = 0.9,
                                  smooth.type = F,p.value = 0.95,seq.long = 0.01)
basic.plot <- ggplot() +
  geom_pointrange(data = dat.gmt.ci, 
                  aes(x = type, y = mle, ymin = li, ymax = hi), color = "#925E9FB2") +
  geom_smooth(data = dat.gmt.ci,aes(x = type, y = mle)) +  ylab("") + 
  geom_hline(yintercept = 0,linetype = 2) + 
  scale_x_continuous(breaks = dat.gmt.ci$type) + 
  scale_y_continuous(breaks = round(seq(-1.2,1.6,by = 0.4),2))+
  xlab(expression(paste("Threshold for mean temperature on reporting rate (",degree,"C )",sep = ""))) +
  ylab(expression(paste(psi[MT],sep = ""))) +
  theme_bw(base_family = "serif",base_size = 15) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(margin = margin(0,2,0,0)))

############ U-shaped

gmt.result.u <- threshold.caculate.gmt(dir = "E:/Project/Mumps/POMP/Codes/Analysis-Codes/Threshold-MT-RR/U-shape",
                                       stri = "hangzhou_U-MT3ths")
gmt.result.u$target <- round(log(gmt.result.u$target),3)
gmt.result.u <- filter(gmt.result.u,target > (-0.25))
theshold.plot.2(gmt.result.u,figure.file = "Figures/Threshold-MT-RR/Hangzhou/U-shape/",ncol = 3)

dat.gmt.ushape.ci <- ci.threshold(gmt.result.u,parameter = "gmt",alpha = 0.9,
                             smooth.type = F,p.value = 0.95,seq.long = 0.01)
gmt.hz.u <- filter(dat.gmt.ushape.ci,type == "20&7");gmt.hz.u$city <- "HZ" 


dat.gmt.ushape.ci[,c(3:5)] <- apply(dat.gmt.ushape.ci[,c(3:5)],2,function(data){
ifelse(abs(data) < 0.0001,"0.00",num_pad_right(data,num = 4, symbol = "0"))
})

dat.gmt.ushape.ci <- mutate(dat.gmt.ushape.ci, 
                            type2 = lapply(stri_split_fixed(dat.gmt.ushape.ci$type,"&"),function(data)data[1]),
                            type1 = lapply(stri_split_fixed(dat.gmt.ushape.ci$type,"&"),function(data)data[2]),
                            p.value = p.value.test(li,hi),
                            or.text = paste("(",li,",",hi,")",sep = ""))
dat.gmt.ushape.ci[,c(6,7)] <- apply(dat.gmt.ushape.ci[,c(6,7)],2,as.numeric)

dat.gmt.ushape.ci <- filter(dat.gmt.ushape.ci,type2 > 17 & type1 > 6)
dat.gmt.ushape.ci <- filter(dat.gmt.ushape.ci,!type == "19&7")

ggplot(dat.gmt.ushape.ci, aes(x = type1, y = type2)) + 
  geom_point(aes(size = as.numeric(mle), color = factor(p.value))) + 
  scale_size_continuous(range = c(0.2,1.0)*10,
                        breaks = c(seq(-0.1,0.2,0.05)),
                        labels = seq(-0.1,0.2,0.05),
                        guide = guide_legend(title = expression(psi[MT]),keywidth = 2,
                                             override.aes = list(colour = "grey"))) +
  geom_text(aes(label = mle,parse = T,nudge_y = 0.2),size = 4) +
  geom_text(data = dat.gmt.ushape.ci,aes(x = type1, y = type2 - 0.2,label = or.text),size = 4) +
  theme_bw(base_family = "serif",base_size = 12) +
  scale_colour_manual(values = c("grey","red"),labels = c("Insiginificant","Siginificant"),
                      guide = guide_legend(title = "Statistical Level",keyheight = 1.3,override.aes = list(size=5)))+
  theme(legend.key.size = unit(1.2,"cm")) + 
  scale_x_continuous(limits = c(6.5,13.5),breaks = seq(7,13,2),labels = seq(7,13,2)) +
  scale_y_continuous(limits = c(17,24),breaks = seq(18,24,2),labels = seq(18,24,2)) +
  xlab(expression(paste("Cold threshold (",degree,"C)",sep = ""))) + 
  ylab(expression(paste("Hot threshold (",degree,"C)",sep = "")))

ggsave(filename = "Figures/Threshold-MT-RR/Hangzhou/combination.pdf",
       width = 24,height = 14,dpi = 400,units = "cm")


filter(dat.gmt.ushape.ci,type1 == 8)


dat.gmt.target <- filter(gmt.result.u,type == "22&7")
plot.target <- theshold.plot.single(dat.gmt.target,figure.file = paste(figue.mt2.dir,"Hangzhou/",sep = ""))

basic.plot + annotation_custom(ggplotGrob(plot.target),
                               xmin = 17, xmax = 23, ymin = 0.1, ymax = 0.5)

ggsave(filename = "Figures/Threshold-MT-RR/Hangzhou/MT-RR-Total.pdf",dpi = 300,
       width = 16 ,height = 14,units = "cm")


gmt.u <- rbind(gmt.hz.u,gmt.qz.u)

write.csv(threshold.tr,file = "Results/threshold-transmission.csv",row.names = F)
write.csv(gmt.u,file = "Results/threshold-reporting.csv",row.names = F)
