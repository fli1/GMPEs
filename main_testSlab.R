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
source("rank.r")

## note the sequence in the list of GMPEs has to be the same as the calculation sequence in GMPEcal.r
list.gmpe = c("Y97SZ",
              "AB03SZ",
  "Z06SZ",
  "Z16SZ",
  "AEA16SZ",
  "MBR17SZ",
  "CB15SZ",
  "IdiniEtAl16SZ"
)
list.p = c(0.01, 0.04, 0.1, 0.2, 1,2)
# list.p = c(0.01,0.02,0.03, 0.04,0.05,0.075,##0.04,0.03,
#            0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.75,
#            1,1.5,2,3,4,5,7.5,10) ##

readObserved = 0
computePred = 0
testGMPE = 1

flag.source = "inslab" 
flag.scenario = "period" #"magnitude" #"distance" # #"period" #
flag.period = "all" #"PGA" #0.05 ## "0.2""PGA" ##"all" #

fname.predicted = sprintf("predicted_%s_%s_%s", flag.source, flag.scenario, flag.period) #,vs30)
fname.observed = sprintf("observed_%s_%s_%s", flag.source, flag.scenario, flag.period) #,vs30)
fname.rank = sprintf("rank_%s_%s_%s", flag.source, flag.scenario, flag.period) #,vs30)
figname.base = sprintf("%s_%s_%s", flag.source, flag.scenario, flag.period) #,vs30)
if(readObserved == TRUE){
  #### read Chile database for comparison
  fname = "C:\\2017\\Sources\\GM Data\\Data_from__Chile_Strong_Ground_Motion_flatfile__table.csv"
  dataraw = read.csv(fname)
  
  ## target inslab events; horizontal components 
  data.temp = dataraw[dataraw$TectEnvi ==flag.source
                      & dataraw$Comp!="Z",]
  # & dataraw$Vs30..m.s.>360,]# & dataraw$Vs30..m.s.<1500,]
  
  # # data.temp = data.temp[data.temp$Mag..Mw.>7.7,]
  # # temp = as.numeric(as.character(data.temp$HypDist..km.))
  # 
  # data.temp = data.temp[data.temp$HypDist..km.>=50 & data.temp$HypDist..km.<=150,]
  # data.temp = data.temp[data.temp$EvDepth..km.>=75 & data.temp$EvDepth..km.<=125,]
  # temp = as.numeric(as.character(data.temp$Mag..Mw.))
  # 
  # temp2 = ifelse(data.temp$Vs30..m.s.<180,"E",ifelse(data.temp$Vs30..m.s.<360,"D",ifelse(data.temp$Vs30..m.s.<760,"C",ifelse(data.temp$Vs30..m.s.<1500,"B","A"))))
  # data.recorded = data.frame(variable= temp, value=data.temp$PGA..cm.sec2./981,group = temp2)
  
  temp = data.temp[,c(29,seq(32,82))]
  temp[temp<0] = -981
  temp=temp/981
  
  observed = cbind(data.temp[,c(1,2,5,15,17,seq(19,24,1))],
                   temp[c(1,8,18,20,35,37,39)])
  colnames(observed)[3:5] = c("depth","mag","Vs30")
  colnames(observed)[6:11] = c("EpiDist","HypDist" ,"RupDist" ,
                               "RupDist_PubliPlane",	"RupDist_CMT" ,"RupDist_CSN")
  # colnames(observed)[12:18] = c("PGA","S004","S")
  # list.period= c(0.01,0.015,0.02,0.05,0.025,0.03,0.035,
  #                0.04,0.045,0.055,0.06,0.065,0.07,0.075,0.08,
  #                0.085,0.09,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,
  #                0.5,0.55,0.6,0.65,0.7,
  #                0.75,0.8,0.85,0.9,1,1.5,2,
  #                2.5,3,3.5,4,4.5,5,5.5,6,
  #                6.5,7,7.5,8,8.5,9,10)
  rm("dataraw","temp")
  write.table(observed, sprintf("%s.csv",fname.observed),
              sep=",",row.names=FALSE)
  # observed = data.frame(variable=NULL)
}

if(computePred == TRUE) {
  observed = read.csv(sprintf("%s.csv",fname.observed))
  no.observed = nrow(observed)
  
  ### calculation parameters
  
  # if(flag.period == "PGA") {
  #   list.p = c(0.01,0.5,1,3)
  # }
  # F_FABA = 0  ## 0=forearc; 1=backarc
  
  predicted = NULL
  for (i in 1:no.observed){
    list.dist = observed$HypDist[i] #exp(seq(log(20),log(500),length.out=50)) ##c(50) #,50,100,200) ##100 #c(75) #
    list.mag =  observed$mag[i] #seq(5,8.1,0.1) #c(9) # c(8)  ##c(7.5,6.5)
    list.zh =  observed$depth[i] #25 #20 #c(25,50,100,200) 
    list.vs30 = observed$Vs30[i] #560 #761 # avoid use 760 which is the boundary limit. Vs30 is site class C in AB03
    
    temp = GMPE.cal(list.gmpe,list.mag, list.dist, list.p, list.zh, 
                    list.vs30,
                    flag.source,  
                    F_FABA=1, F_DeltaC1=0, Delta_C1_input=-0.3, F_sigma=0, ##AEA16 & MBR17SZ
                    region="world",  ##AB03
                    FR =0, MS =0, ##Zhao 06
                    Xv = 0, ## Zhao 16
                    T_star = 0.1, ## Idini16
                    AddMedian=0,
                    flag.testGMPE=0,
                    flag.scenario, flag.period)
    out.Sa = temp[[1]]
    list.GMPEs = temp[[2]]
    # ymax = temp[[3]]
    
    colnames(out.Sa) = c(flag.scenario, paste(rep(list.GMPEs,each=4),
                                              c("sa","sigmatotal","sigma","tau"),sep="_"))
    out.Sa = cbind(out.Sa, ReID = observed$ReID[i], EQID = observed$EQID[i])
    predicted = rbind(predicted, out.Sa)
  }
  
  write.table( predicted, sprintf("%s.csv",fname.predicted),
               row.names = F, sep=",")
}

#####=========== test the quality of GMPEs
if(testGMPE==TRUE) {
  
  out.rank = NULL
  # figname = sprintf("predicted_inslab_period_all.csv")
  predicted = read.csv(sprintf("%s.csv",fname.predicted))
  observed = read.csv(sprintf("%s.csv",fname.observed))
  
  predicted.nrow=nrow(predicted)
  no.gmpe = length(list.gmpe)
  no.periods = length(list.p)
  sum.diff = NULL
  sum.inter = NULL
  sum.stats = NULL
  for(i in 1:no.gmpe){
    temp = c(list.gmpe[i],"medLH", "medianDiff", "meanDiff", "stdDiff")
    sum.temp.diff=NULL
    sum.temp.inter = NULL
    
    for (j in 1:no.periods) {
      ## get subset of observed
      observed.sub = observed[,c(seq(1,11,1),j+11)]
      colnames(observed.sub)[12] = "sa"
      ## get subset of predicted
      temp.colno = c(1,seq((i-1)*4+2,(i-1)*4+5,1),no.gmpe*4+2,no.gmpe*4+3)
      predicted.sub = predicted[predicted$period==list.p[j],temp.colno]
      
      colnames(predicted.sub)[c(2,3,4,5)] = c("sa","totalsigma","sigma","tau")
      
      totalresidual = GetResidual(observed.sub, predicted.sub)
      diff =totalresidual[[1]]
      LH = totalresidual[[2]]
      Inter_residl=totalresidual[[3]]
      LH_inter = totalresidual[[4]]
      Intra_residl = totalresidual[[5]] 
      LH_intra = totalresidual[[6]]
      iflag = totalresidual[[7]]
      
      #### plot histogram here
      figname = sprintf("hist_%s_%s_%s.png",figname.base,list.gmpe[i],list.p[j])
      # plothist(figname, totalresidual)
      #### save the mean for future plot
      temp.diff = data.frame(cbind(diff, Intra_residl, LH, LH_intra))
      temp.inter = data.frame(cbind(Inter_residl, LH_inter))
      temp.diff$period = list.p[j]
      temp.inter$period = list.p[j]
      temp.diff$gmpe = list.gmpe[i]
      temp.inter$gmpe = list.gmpe[i]
      
      sum.temp.diff= rbind(sum.temp.diff, temp.diff)
      sum.temp.inter = rbind(sum.temp.inter,temp.inter)
      
      sum.stats = rbind(sum.stats, c(mean(diff),mean(LH),
                                     mean(Inter_residl), mean(LH_inter),
                                     mean(Intra_residl), mean(LH_intra),
                                     list.p[j], list.gmpe[i]))
      ####
      out = Rank(diff, LH)
      tmp.rank = out[1]
      if (iflag == 1) {
        tmp.rank = c(tmp.rank, Rank(Inter_residl, LH_inter)[1])
        tmp.rank = c(tmp.rank, Rank(Intra_residl, LH_intra)[1])
      }
      tmp.rank = paste(tmp.rank, collapse="/")
      ## write output into a single file
      out = c(tmp.rank,out[seq(2,5,1)])
      temp = cbind(temp, out)
      
    }
    out.rank = rbind(out.rank, temp)
    
    sum.diff = rbind(sum.diff, sum.temp.diff)
    sum.inter = rbind(sum.inter,sum.temp.inter)
    
  }
  
  colnames(out.rank) = c("gmpe", paste0("sa",list.p))
  write.table(out.rank, sprintf("%s.csv",fname.rank),sep=",",row.names=F)
  
  ## plot distribution of differences 
  ##1. plot vs periods
  
  iflag = c(0,1,1,1,1,1,0,1)
  plot.xperiod(iflag, figname.base, sum.diff, sum.inter, sum.stats,
               list.gmpe, list.p)
  
}



