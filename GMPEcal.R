GMPE.cal <- function (list.gmpe,list.mag, list.rupdist, list.hypodist,
                      list.p, list.zh, vs30,
                      flag.source,  
                      F_FABA, F_DeltaC1=0, Delta_C1_input=-0.3, F_sigma=0, ##AEA16 & MBR17SZ
                      region="world",  ##AB03
                      FR =0, MS =0, ##Zhao 06
                      Xv = 0, ## Zhao 16
                      T_star = 0.1, ## Idini16
                      AddMedian=0,
                      flag.testGMPE=0,
                      flag.scenario, flag.period){
  
  out.Sa = NULL
  out.Sigma = NULL
  list.GMPEs = NULL
  flag.AB03SZ = 0
  flag.Y97SZ = 0
  flag.Z06SZ = 0
  flag.Z16SZ = 0
  flag.AEA16SZ = 0
  flag.MBR17SZ = 0
  flag.CB15SZ = 0
  flag.IdiniEtAl16SZ = 0
  
  if("AB03SZ" %in% list.gmpe) {
    flag.AB03SZ = 1
  }
  if("Y97SZ" %in% list.gmpe) {
    flag.Y97SZ = 1
  }
  if("Z06SZ" %in% list.gmpe) {
    flag.Z06SZ = 1
  }
  if("Z16SZ" %in% list.gmpe) {
    flag.Z16SZ = 1
  }
  if("AEA16SZ" %in% list.gmpe) {
    flag.AEA16SZ = 1
  }
  if("MBR17SZ" %in% list.gmpe) {
    flag.MBR17SZ = 1
  }
  if("CB15SZ" %in% list.gmpe) {
    flag.CB15SZ = 1
  }
  if("IdiniEtAl16SZ" %in% list.gmpe) {
    flag.IdiniEtAl16SZ = 1
  }
  
  if(flag.Y97SZ == 1) {
    ###CB15################################################################
    ###   
    title = paste0("YEA97") #, vs30.site,"_amp",flag.ampf)
    # datafile.name = paste0("data_",title,".csv")
    # fout.name = paste0(title,".csv")
    descrpt = "Youngs et al. (1997)"
    list.GMPEs = c(list.GMPEs, descrpt)
    # region = "world"
    dist.metrics = "Rrup" ## Rrup 
    
    ### calcualte Sa and sigma that are amplified to site condition
    output = YEA97SZ.itr (list.mag, list.rupdist, list.p, list.zh, 
                          flag.source, vs30, 
                          AddMedian)
    
    ## for plotting purpose
    output.sa = output[[1]]  
    output.sigma = output[[2]]
    out.temp = GetSubdata(output.sa, output.sigma, list.p, flag.scenario, flag.period)
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
    # region = "world"
    dist.metrics = "Rrup"
    
    if(region == "World"){
      Zl=0
    }
    if(region == "Japan"){
      Zl=2
    }
    if(region == "Cascadia"){
      Zl=1
    }
    
    ### calcualte Sa and sigma that are amplified to site condition
    # print(region)
    output = AB03SZ.itr (list.mag, list.rupdist, list.p, list.zh,
                         flag.source, vs30, Zl)
    output.sa = output[[1]]
    output.sigma = output[[2]]
    out.temp = GetSubdata(output.sa, output.sigma, list.p, flag.scenario, flag.period)
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
    # region = ""
    dist.metrics = "Rx" ## a mixture
    
    
    ### calcualte Sa and sigma that are amplified to site condition
    output = Z06SZ.itr (list.mag, list.rupdist, list.hypodist, list.p, list.zh, 
                        flag.source, vs30, FR, MS)
    output.sa = output[[1]]
    output.sigma = output[[2]]
    out.temp = GetSubdata(output.sa, output.sigma, list.p, flag.scenario, flag.period)
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
  
  if(flag.Z16SZ == 1) {
    ###Zhao et al. 2006################################################################
    ### 
    
    title = paste0("Zhao16") #, vs30.site,"_amp",flag.ampf)
    # datafile.name = paste0("data_",title,".csv")
    # fout.name = paste0(title,".csv")
    descrpt = "Zhao et al. (2016)"
    list.GMPEs = c(list.GMPEs, descrpt)
    # region = ""
    dist.metrics = "Rx" ## a mixture
    
    ### calcualte Sa and sigma that are amplified to site condition
    output = Zhao16SZ.itr (list.mag, list.rupdist, list.hypodist, list.p, list.zh-30, 
                           flag.source, vs30, Xv, AddMedian)
    output.sa = output[[1]]
    output.sigma = output[[2]]
    out.temp = GetSubdata(output.sa, output.sigma, list.p, flag.scenario, flag.period)
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
    # region = "global"
    dist.metrics = "R mixed" ## Rrup for interface and Rhypo for inslab
    
    ### calcualte Sa and sigma that are amplified to site condition
    output = AEA16SZ.itr (list.mag, list.rupdist, list.hypodist, list.p, list.zh, 
                          flag.source, vs30, F_DeltaC1, 
                          Delta_C1_input, F_FABA,
                          F_sigma,
                          AddMedian)
    
    ## for plotting purpose
    output.sa = output[[1]]  
    output.sigma = output[[2]] 
    out.temp = GetSubdata(output.sa, output.sigma,
                          list.p, flag.scenario, flag.period)
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
    # region = "Chile"
    dist.metrics = "R mixed" ## Rrup for interface and Rhypo for inslab
    
    ### calcualte Sa and sigma that are amplified to site condition
    output = MBR17SZ.itr (list.mag, list.rupdist, list.hypodist, list.p, list.zh, 
                          flag.source, vs30, F_DeltaC1, 
                          Delta_C1_input, F_FABA,
                          F_sigma,
                          AddMedian)
    ## for plotting purpose
    output.sa = output[[1]]  
    output.sigma = output[[2]] 
    
    out.temp = GetSubdata(output.sa, output.sigma, list.p, flag.scenario, flag.period)
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
    # region = "Chile"
    dist.metrics = "Rrup" ## Rrup for interface and Rhypo for inslab
    
    ### calcualte Sa and sigma that are amplified to site condition
    output = CB15SZ.itr (list.mag, list.hypodist, list.p, list.zh, 
                         flag.source, vs30, 
                         AddMedian)
    ## for plotting purpose
    output.sa = output[[1]]  
    output.sigma = output[[2]] 
    
    out.temp = GetSubdata(output.sa, output.sigma, list.p, flag.scenario, flag.period)
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
    # region = "Chile"
    dist.metrics = "R mixed" ## Rrup for interface and Rhypo for inslab
    
    ### these were defined by FL for comparison purpose; NOT from the original paper
    if(vs30>760) {T_star = 0.01}   ## rock
    else if(vs30>360) {T_star = 0.1}
    else if(vs30>180) {T_star = 0.3}
    else if(vs30>100) {T_star = 0.5}
    else {T_star = 1}
    
    ### calcualte Sa and sigma that are amplified to site condition
    output = IdiniEtAl16SZ.itr (list.mag, list.rupdist, list.hypodist, list.p, list.zh, 
                                flag.source, vs30, T_star,
                                AddMedian)
    ## for plotting purpose
    output.sa = output[[1]]  
    output.sigma = output[[2]] 
    out.temp = GetSubdata(output.sa,output.sigma, list.p, flag.scenario, flag.period)
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
  
  return(list(out.Sa,list.GMPEs,ymax))
}