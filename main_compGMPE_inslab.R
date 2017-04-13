## to compare GMPEs
rm(list=ls())
#####--------------  main functions --------------#####
library(nga)

# wd = "C:\\Users\\FELi\\Documents\\others\\PSHA\\GMPEs"
wd = "C:\\2017\\Sources\\GMPE\\GMPE R Codes"
setwd(wd)

# source("siteamp.r")
# source("createtable.r")
# source("_Silva2002DC_Sat.r")
# source("_Toro1997.r")
# source("_Pezeshk2011.R")
# source("_AkB10.r")

source("AB03SZ.r")
source("Youngs97SZ.r")
source("Zhao2006SZ.r")
source("AbrahamsonEtAl16SZ.r")
source("MontalvaEtAl17SZ.r")
source("CB15SZ.r")
source("IdiniEtAl16SZ.r")
source("utility.r")

### calculation parameters
# list.p = c(0.01)
list.p = c(0.01,0.02,0.03, 0.04,0.05,0.075,##0.04,0.03,
           0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.75,
           1,1.5,2,3,4,5,7.5,10) ##

flag.scenario = "magnitude" #"distance" #"period" # "period" ##

flag.period = "PGA" #0.05 #"all" ## "0.2""PGA" ##"all" #
if(flag.period == "PGA") {
  list.p = c(0.01)
}
list.dist = 100 ##exp(seq(log(60),log(500),length.out=50)) ##c(50) #,50,100,200) ##100 #c(75) #
list.mag = seq(5,8.1,0.1) #c(9) # c(8)  ##c(7.5,6.5)
list.zh =  100 #20 #c(25,50,100,200) 

## parameters
flag.source = "inslab" #

vs30 = 560 #761 # avoid use 760 which is the boundary limit. Vs30 is site class C in AB03
F_FABA = 0  ## 0=forearc; 1=backarc
fltstyl = "strike slip" ##"normal" ##"normal" ##
##"reverse"
##"strike slip"  

#### read Chile database for comparison
fname = "C:\\2017\\Sources\\GM Data\\Data_from__Chile_Strong_Ground_Motion_flatfile__table.csv"
dataraw = read.csv(fname)
data.temp = dataraw[dataraw$TectEnvi =="inslab" 
                    & dataraw$Comp!="Z"
                     & dataraw$Vs30..m.s.>360,]# & dataraw$Vs30..m.s.<1500,]

# data.temp = data.temp[data.temp$Mag..Mw.>7.7,]
# temp = as.numeric(as.character(data.temp$HypDist..km.))

data.temp = data.temp[data.temp$HypDist..km.>=50 & data.temp$HypDist..km.<=150,]
data.temp = data.temp[data.temp$EvDepth..km.>=75 & data.temp$EvDepth..km.<=125,]
temp = as.numeric(as.character(data.temp$Mag..Mw.))

temp2 = ifelse(data.temp$Vs30..m.s.<180,"E",ifelse(data.temp$Vs30..m.s.<360,"D",ifelse(data.temp$Vs30..m.s.<760,"C",ifelse(data.temp$Vs30..m.s.<1500,"B","A"))))
data.recorded = data.frame(variable= temp, value=data.temp$PGA..cm.sec2./981,group = temp2)

rm("dataraw")
#### end

flag.AB03SZ = 1
flag.Y97SZ = 1
flag.Z06SZ = 1
flag.AEA16SZ = 1
flag.MBR17SZ = 1
flag.CB15SZ = 1
flag.IdiniEtAl16SZ = 1

flag.plot = 1
flag.testGMPE = 0

out.Sa = NULL
out.Sigma = NULL
list.GMPEs = NULL


if(flag.Y97SZ == 1) {
  ###CB15################################################################
  ###   
  title = paste0("YEA97") #, vs30.site,"_amp",flag.ampf)
  # datafile.name = paste0("data_",title,".csv")
  # fout.name = paste0(title,".csv")
  descrpt = "Youngs et al. (1997)"
  list.GMPEs = c(list.GMPEs, descrpt)
  region = "world"
  dist.metrics = "Rrup" ## Rrup 
  
  
  ### calcualte Sa and sigma that are amplified to site condition
  output = YEA97SZ.itr (list.mag, list.dist, list.p, list.zh, 
                              flag.source, vs30, 
                              AddMedian=0)
  
  ## for plotting purpose
  output.sa = output[[1]]  
  out.temp = GetSubdata(output.sa, list.p, flag.scenario, flag.period)
  out.Sa = cbind(out.Sa, out.temp)
  ## for test GMPE purpose
  if (flag.testGMPE ==1) {
    if(length(list.dist)==1) {
      write.table(output.sa[1,,,], sprintf("test_%s.csv",title), sep=",",row.names=F)  
    } else {
      write.table(output.sa[,,,], sprintf("test_%s.csv",title), sep=",",row.names=F)  
    }
    
  }
}


if(flag.AB03SZ == 1) {
  ###AB03#########################################################################
  ### Atkinson and Boore 2003
  
  title = paste0("AB03") #, vs30.site,"_amp",flag.ampf)
  # datafile.name = paste0("data_",title,".csv")
  # fout.name = paste0(title,".csv")
  descrpt = "Atkinson and Boore (2003)"
  
  list.GMPEs = c(list.GMPEs, descrpt)
  region = "world"
  dist.metrics = "Rrup"
  
  if(region == "world"){
    Zl=0
  }
  if(region == "Japan"){
    Zl=2
  }
  if(region == "Cascadia"){
    Zl=1
  }
  
  ### calcualte Sa and sigma that are amplified to site condition
  output = AB03SZ.itr (list.mag, list.dist, list.p, list.zh,
                       flag.source, vs30, Zl)
  output.sa = output[[1]]
  
  out.temp = GetSubdata(output.sa, list.p, flag.scenario, flag.period)
  out.Sa = cbind(out.Sa, out.temp)
  
  ## for test GMPE purpose
  if (flag.testGMPE ==1) {
    if(length(list.dist)==1) {
      write.table(output.sa[1,,,], sprintf("test_%s.csv",title), sep=",",row.names=F)  
    }else {
      write.table(output.sa[,,,], sprintf("test_%s.csv",title), sep=",",row.names=F)  
    }
    
  }
}


if(flag.Z06SZ == 1) {
  ###Zhao et al. 2006################################################################
  ### 
  
  title = paste0("Z06") #, vs30.site,"_amp",flag.ampf)
  # datafile.name = paste0("data_",title,".csv")
  # fout.name = paste0(title,".csv")
  descrpt = "Zhao et al. (2006)"
  list.GMPEs = c(list.GMPEs, descrpt)
  region = ""
  dist.metrics = "Rx" ## a mixture
  
  ### calcualte Sa and sigma that are amplified to site condition
  output = Z06SZ.itr (list.mag, list.dist, list.p, list.zh, 
                         flag.source, vs30, FR=0, MS=0)
  output.sa = output[[1]]
  
  out.temp = GetSubdata(output.sa, list.p, flag.scenario, flag.period)
  out.Sa = cbind(out.Sa, out.temp)
  ## for test GMPE purpose
  if (flag.testGMPE ==1) {
    if(length(list.dist)==1) {
      write.table(output.sa[1,,,], sprintf("test_%s.csv",title), sep=",",row.names=F)  
    }else {
      write.table(output.sa[,,,], sprintf("test_%s.csv",title), sep=",",row.names=F)  
    }
    
  }
}

if(flag.AEA16SZ == 1) {
  ###Abrahamson et al. 2016################################################################
  ### 
  
  title = paste0("AEA16") #, vs30.site,"_amp",flag.ampf)
  # datafile.name = paste0("data_",title,".csv")
  # fout.name = paste0(title,".csv")
  descrpt = "Abrahamson et al. (2016)"
  list.GMPEs = c(list.GMPEs, descrpt)
  region = "global"
  dist.metrics = "Rrup" ## Rrup for interface and Rhypo for inslab
  
  ### calcualte Sa and sigma that are amplified to site condition
  output = AEA16SZ.itr (list.mag, list.dist, list.p, list.zh, 
                      flag.source, vs30, F_DeltaC1=0, 
                      Delta_C1_input=-0.3, F_FABA,
                      F_sigma=0,
                      AddMedian=0)
    
  ## for plotting purpose
  output.sa = output[[1]]  
  out.temp = GetSubdata(output.sa, list.p, flag.scenario, flag.period)
  out.Sa = cbind(out.Sa, out.temp)
  ## for test GMPE purpose
  if (flag.testGMPE ==1) {
    if(length(list.dist)==1) {
      write.table(output.sa[1,,,], sprintf("test_%s.csv",title), sep=",",row.names=F)  
    }else {
      write.table(output.sa[,,,], sprintf("test_%s.csv",title), sep=",",row.names=F)  
    }
    
  }
  
}


if(flag.MBR17SZ == 1) {
  ###Abrahamson et al. 2016################################################################
  ###   
  title = paste0("MBR17") #, vs30.site,"_amp",flag.ampf)
  # datafile.name = paste0("data_",title,".csv")
  # fout.name = paste0(title,".csv")
  descrpt = "Montalva et al. (2017)"
  list.GMPEs = c(list.GMPEs, descrpt)
  region = "Chile"
  dist.metrics = "Rrup" ## Rrup for interface and Rhypo for inslab
  
  ### calcualte Sa and sigma that are amplified to site condition
  output = MBR17SZ.itr (list.mag, list.dist, list.p, list.zh, 
                        flag.source, vs30, F_DeltaC1=0, 
                        Delta_C1_input=-0.3, F_FABA,
                        F_sigma=0,
                        AddMedian=0)
  ## for plotting purpose
  output.sa = output[[1]]  
  out.temp = GetSubdata(output.sa, list.p, flag.scenario, flag.period)
  out.Sa = cbind(out.Sa, out.temp)
  ## for test GMPE purpose
  if (flag.testGMPE ==1) {
    if(length(list.dist)==1) {
      write.table(output.sa[1,,,], sprintf("test_%s.csv",title), sep=",",row.names=F)  
    }else {
      write.table(output.sa[,,,], sprintf("test_%s.csv",title), sep=",",row.names=F)  
    }
    
  }
}


if(flag.CB15SZ == 1) {
  ###CB15################################################################
  ### 
  
  title = paste0("CB15") #, vs30.site,"_amp",flag.ampf)
  # datafile.name = paste0("data_",title,".csv")
  # fout.name = paste0(title,".csv")
  descrpt = "Contreras and Boroschek (2015)"
  list.GMPEs = c(list.GMPEs, descrpt)
  region = "Chile"
  dist.metrics = "Rrup" ## Rrup for interface and Rhypo for inslab
  
  ### calcualte Sa and sigma that are amplified to site condition
  output = CB15SZ.itr (list.mag, list.dist, list.p, list.zh, 
                        flag.source, vs30, 
                        AddMedian=0)
  ## for plotting purpose
  output.sa = output[[1]]  
  out.temp = GetSubdata(output.sa, list.p, flag.scenario, flag.period)
  out.Sa = cbind(out.Sa, out.temp)
  ## for test GMPE purpose
  if (flag.testGMPE ==1) {
    if(length(list.dist)==1) {
      write.table(output.sa[1,,,], sprintf("test_%s.csv",title), sep=",",row.names=F)  
    }else {
      write.table(output.sa[,,,], sprintf("test_%s.csv",title), sep=",",row.names=F)  
    }
    
  }
}



if(flag.IdiniEtAl16SZ == 1) {
  ###CB15################################################################
  ### 
  
  title = paste0("IdiniEtAl16") #, vs30.site,"_amp",flag.ampf)
  # datafile.name = paste0("data_",title,".csv")
  # fout.name = paste0(title,".csv")
  descrpt = "Idini et al. (2016)"
  list.GMPEs = c(list.GMPEs, descrpt)
  region = "Chile"
  dist.metrics = "Rrup" ## Rrup for interface and Rhypo for inslab
  
  ### calcualte Sa and sigma that are amplified to site condition
  output = IdiniEtAl16SZ.itr (list.mag, list.dist, list.p, list.zh, 
                       flag.source, vs30, T_star=.1,
                       AddMedian=0)
  ## for plotting purpose
  output.sa = output[[1]]  
  out.temp = GetSubdata(output.sa, list.p, flag.scenario, flag.period)
  out.Sa = cbind(out.Sa, out.temp)
  ## for test GMPE purpose
  if (flag.testGMPE ==1) {
    if(length(list.dist)==1) {
      write.table(output.sa[1,,,], sprintf("test_%s.csv",title), sep=",",row.names=F)  
    } else {
      write.table(output.sa[,,,], sprintf("test_%s.csv",title), sep=",",row.names=F)  
    }
    
  }
}

ymax = round(max(out.Sa,na.rm = T),1)+0.1
if(flag.scenario=="period") { 
  out.Sa = cbind(list.p, out.Sa)
}

if(flag.scenario=="distance") { 
  out.Sa = cbind(list.dist, out.Sa)
}
if(flag.scenario=="magnitude") {  
  out.Sa = cbind(list.mag, out.Sa)
}
if(flag.scenario=="depth") {  
  out.Sa = cbind(list.zh, out.Sa)
}

figname = sprintf("comp_%s_%s_%s_vs%s", flag.source, flag.scenario, flag.period,vs30)
colnames(out.Sa) = c(flag.scenario, list.GMPEs)
write.table( out.Sa, sprintf("%s.csv",figname),
             row.names = F, sep=",")
### plot result figure 
if(flag.plot==1){
  ## plot figures
  
  datain = out.Sa
  plot.comp(figname, datain, data.recorded, list.GMPEs, flag.scenario, flag.period, 
            flag.source, vs30,
            list.mag, list.dist, list.zh, ymax)
  
}


