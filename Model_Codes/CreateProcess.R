library(pomp);library(readr)

creat_formula <- function(v.name,effect = "foi"){
  ## @ effect: foi or cr
  
  parm.v.name <- paste("p",v.name,sep = "")
  parm.part <- unlist(lapply(1:length(v.name),function(id){
    parm.tem <- paste("log(",parm.v.name[id],")","*",v.name[id],sep = "")
    return(parm.tem)
  }))
  parm.fin.1 <- paste(parm.part,collapse = "+")    ##### the formulation for the force of infection
  
  if(effect == "foi"){
      parm.fin.2 <- paste("beta=R0*(gamma+mu)*seas*exp(",parm.fin.1,");",sep = "")
    }else{
      parm.fin.2 <- paste("ctr=exp(",parm.fin.1,"+alpha",")/(","exp(",parm.fin.1,"+alpha)","+1);",sep = "") 
      parm.fin.2 <- paste(parm.fin.2,"\n","foi = beta*pow(I+iota,ctr)/pop;",sep = "")
          }
  return(parm.fin.2)
} 

Part.1 <- read_file("Model-Codes/Process-1.txt") 
Part.2 <- read_file("Model-Codes/Process-2.txt") 

creat_process <- function(foi.name,contact.name,full = T){
  
  #### foi.name: which variables to constute the force of infection
  #### full: whether include the effect of climatic factors on contacting rate
  
  foi.v <- creat_formula(v.name = foi.name,effect = "foi")
  if(isTRUE(full)){
    contact.v <- creat_formula(v.name = contact.name,effect = "cr")
  } else{
    contact.v <- "foi = beta*pow(I+iota,alpha)/pop;"
  }
  proc.text <- paste(Part.1,"\n",foi.v,"\n",Part.2,sep = "");
  rproc <- Csnippet(proc.text)
  return(rproc)
}
#creat_process(foi.name = c("MT","AHths"),contact.name = c("RH","MT"),full = T)
