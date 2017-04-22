## to compare GMPEs
rm(list=ls())
#####--------------  main functions --------------#####
library(nga)

# wd = "C:\\Users\\FELi\\Documents\\others\\PSHA\\GMPEs"
# wd = "C:\\2017\\Sources\\GMPE\\GMPE R Codes"
# setwd(wd)

# source("siteamp.r")
# source("createtable.r")
# source("_Silva2002DC_Sat.r")
# source("_Toro1997.r")
# source("_Pezeshk2011.R")
# source("_AkB10.r")

source("AB03SZ.r")
source("Youngs97SZ.r")
source("Zhao2006SZ.r")
source("ZhaoEA16SZ.r")
source("AbrahamsonEtAl16SZ.r")
source("MontalvaEtAl17SZ.r")
source("CB15SZ.r")
source("IdiniEtAl16SZ.r")
source("utility.r")
source("GMPEcal.r")

# #### read Chile database for comparison
# fname = "C:\\2017\\Sources\\GM Data\\Data_from__Chile_Strong_Ground_Motion_flatfile__table.csv"
# dataraw = read.csv(fname)
# data.temp = dataraw[dataraw$TectEnvi =="inslab" 
#                     & dataraw$Comp!="Z"
#                      & dataraw$Vs30..m.s.>360,]# & dataraw$Vs30..m.s.<1500,]
# 
# # data.temp = data.temp[data.temp$Mag..Mw.>7.7,]
# # temp = as.numeric(as.character(data.temp$HypDist..km.))
# 
# data.temp = data.temp[data.temp$HypDist..km.>=50 & data.temp$HypDist..km.<=150,]
# data.temp = data.temp[data.temp$EvDepth..km.>=75 & data.temp$EvDepth..km.<=125,]
# temp = as.numeric(as.character(data.temp$Mag..Mw.))
# 
# temp2 = ifelse(data.temp$Vs30..m.s.<180,"E",ifelse(data.temp$Vs30..m.s.<360,"D",ifelse(data.temp$Vs30..m.s.<760,"C",ifelse(data.temp$Vs30..m.s.<1500,"B","A"))))
# data.recorded = data.frame(variable= temp, value=data.temp$PGA..cm.sec2./981,group = temp2)
# 
# rm("dataraw")
data.recorded = data.frame(variable=NULL)
#### ========== Main Inputs =================
list.gmpe = c(#"AB03SZ",
              # "Y97SZ",
              # "Z06SZ",
              "Z16SZ"
              # "AEA16SZ",
              # "MBR17SZ",
              # "CB15SZ",
              # "IdiniEtAl16SZ"
              )

### calculation parameters
# list.p = c(0.01)
list.p = c(0.01,0.02,0.03, 0.04,0.05,0.075,##0.04,0.03,
           0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.75,
           1,1.5,2,3,4,5,7.5,10) ##

flag.scenario = "magnitude" #"distance" # "period" ##"period" #

flag.period = "PGA" #0.05 #"all" ## "0.2""PGA" ##"all" #
# if(flag.period == "PGA") {
#   list.p = c(0.01,0.5,1,3)
# }
list.dist = 30 #exp(seq(log(20),log(500),length.out=50)) ##c(50) #,50,100,200) ##100 #c(75) #
list.mag = c(9) #seq(5,8.1,0.1) #c(9) # c(8)  ##c(7.5,6.5)
list.zh =  26 #25 #20 #c(25,50,100,200) 

vs30 = 761 #560 #761 # avoid use 760 which is the boundary limit. Vs30 is site class C in AB03

# F_FABA = 0  ## 0=forearc; 1=backarc
flag.source = "interface" #"inslab" #
fltstyl = "strike slip" ##"normal" ##"normal" ##
##"reverse"
##"strike slip"  

flag.plot = 1
# flag.testGMPE = 0

temp = GMPE.cal(list.gmpe,list.mag, list.dist, list.p, list.zh, 
                vs30,
                flag.source,  
                F_FABA=0, F_DeltaC1=0, Delta_C1_input=-0.3, F_sigma=0, ##AEA16 & MBR17SZ
                region="world",  ##AB03
                FR =0, MS =0, ##Zhao 06
                Xv = 0, ## Zhao 16
                T_star = 0.1, ## Idini16
                AddMedian=0,
                flag.testGMPE=1,
                flag.scenario, flag.period)
out.Sa = temp[[1]]
list.GMPEs = temp[[2]]
ymax = temp[[3]]

figname = sprintf("comp_%s_%s_%s_vs%s", flag.source, flag.scenario, flag.period,vs30)
colnames(out.Sa) = c(flag.scenario, list.GMPEs)
write.table( out.Sa, sprintf("%s.csv",figname),
             row.names = F, sep=",")
### plot result figure 
if(flag.plot==1){
  ## plot figures
  
  datain = out.Sa
  plot.comp(figname, datain, data.recorded, 
            list.GMPEs, flag.scenario, flag.period, 
            flag.source, vs30,
            list.mag, list.dist, list.zh, ymax)
  
}


