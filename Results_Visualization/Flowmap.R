##########################################################################################
################################### Codes for model mechanism ############################
##########################################################################################


library(diagram)
library(Cairo);library(pBrackets)
library(grid)
library(ggsci)

#setwd('/home/spatial-r/E盘/Project/Mumps/写作模板/figures')

size.state <- 1.5
size.true <- 1.2

#cairo_pdf("Figures/Flowmap/flowchart2.pdf",width = 16, height = 8, family = "serif",bg="transparent",pointsize = 15)
#cairo_pdf("flowchart2.pdf",width = 16, height = 8, family = "serif",bg="transparent",pointsize = 15)
tiff(file = paste("Figures/MS-Figures/","Figure 9.tiff",sep = ""),width = 1800, height = 1600, res = 300,units = "px")


openplotmat()

par(mai =c(0,0,0,0),mar =c(0,0,0,0))
elpos <- coordinates(c(4))
elpos[,2] <- elpos[,2] + 0.2

for (i in 1:(nrow(elpos) - 1)){
  straightarrow(from = c(elpos[i, 1],elpos[i, 2]),
                to  = c(elpos[i + 1, 1],elpos[i + 1, 2]),
                lwd = 1, arr.pos = 0.6, arr.length = 0.2,arr.width = 0.15)
}

curvedarrow(from = c(elpos[1, 1],elpos[1, 2]),
            to  = c(elpos[4, 1],elpos[4, 2]), lwd = 1, 
            arr.pos = 0.5, curve = -0.3, dr = 0.01,arr.length = 0.2,arr.width = 0.15)

curvedarrow(from = c(elpos[4, 1],elpos[4, 2]),
            to  = c(elpos[1, 1],elpos[1, 2]), lwd = 1, 
            arr.pos = 0.5, curve = 0.2, dr = 0.01,arr.length = 0.2,arr.width = 0.15)


labels <- c("S","E","I","R")

for (i in 1:4){
  textellipse(elpos[i,], 0.08, lab = labels[i], box.col = "white",
              shadow.col = "white", shadow.size = 0.005, cex = size.state)
}


v.text.pos <- as.matrix(data.frame(x = c(0.5,0.6,0.6),y = c(0.98,0.88,0.82)))

textplain(v.text.pos[1,], 0.06, lab = "Vaccine effect",cex = 1.2)
#textplain(c(v.text.pos[1,1],v.text.pos[1,2] - 0.12), 0.06, lab = expression(paste("Vaccine Coverage (",
#                                                       omega,") * Vaccine Effectivity (",xi,")",sep="")),cex = 1.2)


f.text.pos <- as.matrix(data.frame(x = c((elpos[1,1] + elpos[2,1])/2), y = c(0.3)))


textplain(c(f.text.pos[1,1],f.text.pos[1,2] + 0.44), 0.08, lab = expression(lambda(t)), cex = size.true)
textplain(c((elpos[2,1] + elpos[3,1])/2,f.text.pos[1,2] + 0.44), 0.08, lab = expression(sigma), cex = size.true)
textplain(c((elpos[3,1] + elpos[4,1])/2,f.text.pos[1,2] + 0.44), 0.08, lab = expression(gamma), cex = size.true)

textplain(c((elpos[1,1] + elpos[2,1])/2, elpos[2,2] - 0.15),
          lab = expression(paste("Holiday effect ",theta(t),sep = "")),
          cex = size.true)
textplain(c((elpos[1,1] + elpos[2,1])/2, elpos[2,2] - 0.20),
          lab = expression(paste("Transmission rate ",beta(t),sep = "")),
          cex = size.true,col = "#FFA07A")

textplain(c((elpos[1,1] + elpos[2,1])/2, elpos[2,2] - 0.25),
          lab = expression(paste("Effective contact rate ",kappa(t),sep = "")),
          cex = size.true)


###  lamda up

straightarrow(to = c(f.text.pos[1,1],elpos[2,2] - 0.01),
              from = c(f.text.pos[1,1], elpos[2,2] - 0.12),lty = 2,
              lwd = 1, arr.pos = 0.6, arr.length = 0.2,arr.width = 0.15)


straightarrow(from = c(elpos[1,1],elpos[1,2] + 0.20),
              to  = c(elpos[1,1], elpos[1,2] + 0.08),lty = 1,
              lwd = 1, arr.pos = 0.6, arr.length = 0.2,arr.width = 0.15)

textplain(c(elpos[1,1], elpos[1,2] + 0.23),lab = "Birth",cex = size.state)
textplain(c(elpos[1,1] + 0.015, elpos[1,2] + 0.14),lab = expression(v),cex = size.state)


for (i in c(1:4)){
  straightarrow(from = c(elpos[i,1],elpos[1,2] - 0.09),
                to  = c(elpos[i,1], elpos[1,2] - 0.21),lty = 1,
                lwd = 1, arr.pos = 0.6, arr.length = 0.2,arr.width = 0.15)
  textplain(c(elpos[i,1], elpos[1,2] - 0.22),lab = "Death",cex = size.state)
  textplain(c(elpos[i,1]+ 0.015,elpos[1,2] - 0.16),lab = expression(mu),cex = size.state)
}

####    case up

straightarrow(from = c((elpos[2,1] + elpos[3,1])/2 + 0.25, elpos[1,2]),
              to  = c((elpos[2,1] + elpos[3,1])/2 + 0.25, elpos[1,2] - 0.20),lty = 2,
              lwd = 1, arr.pos = 0.6, arr.length = 0.2,arr.width = 0.15)

#textplain(c((elpos[2,1] + elpos[3,1])/2, elpos[1,2] - 0.22),lab = "Cases",cex = size.state)
#textplain(c((elpos[2,1] + elpos[3,1])/2 + 0.03, elpos[1,2] - 0.11),lab = "Negbin",cex = size.true)
#textplain(c((elpos[2,1] + elpos[3,1])/2 + 0.043, elpos[1,2] - 0.145),lab = "distribution",cex = size.true)

textplain(c((elpos[2,1] + elpos[3,1])/2 + 0.25, elpos[1,2] - 0.22),lab = "Cases",cex = size.state)
textplain(c((elpos[2,1] + elpos[3,1])/2 + 0.284, elpos[1,2] - 0.11),lab = "Negbin",cex = size.true)
textplain(c((elpos[2,1] + elpos[3,1])/2 + 0.297, elpos[1,2] - 0.145),lab = "distribution",cex = size.true)


##### case blue

#straightarrow(from = c(elpos[3,1] - 0.05, elpos[1,2] - 0.26),
#              to  = c((elpos[2,1] + elpos[3,1])/2 + 0.05, elpos[1,2] - 0.19),lty = 3,lcol = "blue",
#              lwd = 1, arr.pos = 0.6, arr.length = 0.2,arr.width = 0.15)

#textplain(c(elpos[3,1] - 0.05, elpos[1,2] - 0.30),col ="#483D8B",
#          lab = expression(paste(C[n]," | ",P[EI]," ~ ","Negbin(",rho(t),P[EI],", ",rho(t)^2,psi^2,P[EI]^2,")",sep = "")),
#          cex = size.true)

straightarrow(from = c(elpos[4,1] - 0.2, elpos[1,2] - 0.26),
              to  = c((elpos[3,1] + elpos[4,1])/2 - 0.02, elpos[1,2] - 0.12),lty = 3,lcol = "blue",
              lwd = 1, arr.pos = 0.6, arr.length = 0.2,arr.width = 0.15)


textplain(c(elpos[3,1] + 0.05, elpos[1,2] - 0.30),col ="#483D8B",
          lab = expression(paste(C[n]," | ",Y[n]," ~ ","Negbin(",rho(t),Y[n],", ",rho(t)^2,psi^2,Y[n]^2,")",sep = "")),
          cex = size.true)

textplain(c(elpos[1,1]+ 0.15, elpos[1,2] - 0.40),lab = expression(paste("1 + 2(1-p)",epsilon,", during the holiday",sep = "")),cex = size.true)
textplain(c(elpos[1,1]+ 0.15, elpos[1,2] - 0.46),lab = expression(paste("1 - 2p",epsilon,", during the school term",sep = "")),cex = size.true)
#grid.locator(unit = "native")
grid.brackets(250, 390, 250, 350, lwd = 1, col = "black",curvature = 0.9,h = unit(0.01, 'native'))
textplain(c(elpos[1,1] + 0.015, elpos[1,2] - 0.43),lab = expression(paste(theta,"(t) =",sep = "")),cex = size.true)


textplain(c(elpos[1,1] + 0.07, elpos[1,2] - 0.54),lab = expression(paste(kappa,"(t) = ",(I[(t)] + iota)^{alpha},"/N(t)",sep = "")),cex = size.true)

textplain(c(elpos[1,1] + 0.074, elpos[1,2] - 0.63),lab = expression(paste(beta,"(t) = ",eta,"(t)","(",gamma,"+",mu,")",
                                                                          R[0],Delta,"",tau[t],sep = "")),cex = size.true)

grid.brackets(200, 460, 200, 360, lwd = 1, col = "black",curvature = 0.85,h =unit(0.02, 'native'))
textplain(c(elpos[1,1] - 0.029, elpos[1,2] - 0.52),lab = expression(paste(lambda,"(t)",sep = "")),cex = size.true)



textround(c(elpos[3,1]+0.05, elpos[1,2] - 0.5),radx = 0.18,rady = 0.05,box.col ="#8470FF",
          lab = expression(paste("Hypothesis 2 on reporting rate: ",rho,(t),
                                 "= 1/(exp(-(",psi[x[i]]*x[i],"+",rho[0],"))+1)",sep = "")),cex = size.true)

textround(c(elpos[3,1] + 0.052, elpos[1,2] - 0.63),radx = 0.17,rady = 0.05,box.col = "#FFA07A",
          lab = expression(paste("Hypothesis 1 on trasnsmission rate: ",eta,(t),
                                 "= exp(",eta[x[i]]*x[i],")",sep = "")),cex = size.true)

dev.off()

