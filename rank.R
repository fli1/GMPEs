

erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1


GetResidual <- function (observed.sub, predicted.sub) {
  # 'read from spreadsheet "Interim Results" and get spectral ordinates of interest
  
  RecordNo = nrow(observed.sub) 
  
  sumYij = CalSumYij(RecordNo, observed.sub, predicted.sub)
  
  preID = predicted.sub$EQID[1]
  # eid = 1
  out.diff = NULL
  out.LH = NULL
  out.Inter_residl = NULL
  out.LH_inter = NULL
  out.Intra_residl = NULL
  out.LH_intra = NULL
  i=2
  for ( i in 1:RecordNo){
    sa_observed = observed.sub$sa[i]
    sa_pred = predicted.sub$sa[i]
    totalsigma = predicted.sub$totalsigma[i]
    sigma = predicted.sub$sigma[i]
    tau = predicted.sub$tau[i]
    eventID = predicted.sub$EQID[i]
    
    diff = (sa_observed - sa_pred) / totalsigma
    LH = 1 - erf(abs(diff) / sqrt(2))
    out.diff = c(out.diff, diff)
    out.LH = c(out.LH, LH)
    
    if (sigma > 0 & tau > 0) {
      temp = ((sigma ^ 2) * sumYij[i]) / (sumYij[i + RecordNo] * (sigma ^ 2) + (tau ^ 2))
      if (eventID != preID) {
        Inter_residl =  temp / sigma  ## one per each event
        LH_inter = 1 - erf(abs(Inter_residl) / sqrt(2))
        out.Inter_residl = c(out.Inter_residl, Inter_residl)
        out.LH_inter = c(out.LH_inter, LH_inter)
        # eid = eid + 1
        preID = eventID
      }
      Intra_residl = sa_observed - sa_pred - temp  ## one per each record
      Intra_residl = Intra_residl / tau
      LH_intra = 1 - erf(abs(Intra_residl) / sqrt(2))
      iflag = 1
      out.LH_intra = c(out.LH_intra, LH_intra) 
      out.Intra_residl = c(out.Intra_residl, Intra_residl)
    } else {
      Inter_residl = 0
      Intra_residl = 0
      LH_inter = 0
      LH_intra = 0
      iflag = 0
      if (eventID != preID) {
        out.Inter_residl = c(out.Inter_residl, Inter_residl)
        out.LH_inter = c(out.LH_inter, LH_inter)
        # eid = eid + 1
        preID = eventID
      }
      out.LH_intra = c(out.LH_intra, LH_intra) 
      out.Intra_residl = c(out.Intra_residl, Intra_residl)
    }
  }
  
  return(list(out.diff, out.LH,out.Inter_residl,out.LH_inter,
              out.Intra_residl,out.LH_intra, iflag))
}


CalSumYij <- function(RecordNo, observed.sub, predicted.sub) {
  temp = 0
  prevID = predicted.sub$EQID[1]
  stID = 1
  edID = 1
  out.sum = NULL
  out.no = NULL
  for(i in 1:RecordNo) {
    sa_observed = observed.sub$sa[i]
    sa_predicted = predicted.sub$sa[i]
    # 'sa_observed = Exp(sa_observed)
    eventID = predicted.sub$EQID[i]
    if (eventID == prevID) {
      temp = temp + sa_observed - sa_predicted
    } else {
      edID = i
      
      for (j in stID: edID){
        out.sum = c(out.sum, temp)
        out.no = c(out.no,edID - stID)
      }
      stID = edID
      temp = sa_observed - sa_predicted
    }
    prevID = eventID
  } 
  
  edID = i
  for (j in stID:edID){
    out.sum = c(out.sum, temp)
    out.no = c(out.no, edID - stID + 1)
  }
  return(c(out.sum, out.no))
}


Rank <- function(diff,LH) {
  
  # limits_diff = Array(-4, 4)
  # limits_LH = Array(0, 1)
  
  # ' for total residual
  meanDiff = mean(diff)
  medianDiff = median(diff)
  stdDiff = sd(diff)
  medLH = median(LH)
  
  # out = Hist(17, diff, limits_diff, row_no_st, 68, iperiod, shtname)
  # out = Hist(11, LH, limits_LH, row_no_st, 76, iperiod, shtname)
  
  # 'ranking criteria
  # print(medianDiff)
  # print(meanDiff)
  if(abs(medianDiff) == 0 & abs(meanDiff) == 0) {
    rank = "-"
  } else if (abs(medLH) >= 0.4 & abs(medianDiff) < 0.25 & abs(meanDiff) < 0.25 & abs(stdDiff) < 1.125) {
    rank = "A"
  } else if (abs(medLH) >= 0.3 & abs(medianDiff) < 0.5 & abs(meanDiff) < 0.5 & abs(stdDiff) < 1.25) {
    rank = "B"
  } else if (abs(medLH) >= 0.2 & abs(medianDiff) < 0.75 & abs(meanDiff) < 0.75 & abs(stdDiff) < 1.5) {
    rank = "C"
  } else {
    rank = "D"
  }
  return(c(rank, medLH, medianDiff, meanDiff, stdDiff))
  
}

# 
plothist <- function(figname, totalresidual, 
         diff.limits=c(-2,2),
         LH.limits=c(0,1)) {
  
  diff =totalresidual[[1]]
  LH = totalresidual[[2]]
  Inter_residl=totalresidual[[3]]
  LH_inter = totalresidual[[4]]
  Intra_residl = totalresidual[[5]] 
  LH_intra = totalresidual[[6]]
  iflag = totalresidual[[7]]
  
  
  diff.step = 0.2
  LH.step = 0.1
  diff.list = seq(-10,10,diff.step)
  diff.list = diff.list[diff.list>=diff.limits[1] & diff.list<=diff.limits[2]]
  
  # diff.list = c(-10,diff.list,10)
  LH.list = seq(0,1,LH.step)
  
  diff[diff>diff.limits[2]] = diff.limits[2]
  diff[diff<diff.limits[1]] = diff.limits[1]
  
  Inter_residl[Inter_residl>diff.limits[2]] = diff.limits[2]
  Inter_residl[Inter_residl<diff.limits[1]] = diff.limits[1]
  
  Intra_residl[Intra_residl>diff.limits[2]] = diff.limits[2]
  Intra_residl[Intra_residl<diff.limits[1]] = diff.limits[1]
  # # 'Counting the number of occurrences for each of the bins
  # no.total = length(diff)
  # freq = array(0, dim = c(length(diff.list)-1,3))
  # 
  # for (i in 1:no.total) {
  #   diffs = c(diff[i], Inter_residl[i], Intra_residl[i])
  #   for(j in 1:3){
  #     idiff = diffs[j]
  #     ibin = which(idiff<=diff.list)[1]-1
  #     freq[ibin,j] = freq[ibin,j]+1
  #   }
  # }
  # freq = freq/no.total
  # 
  # ## frequency for LH
  # no.total = length(LH)
  # freqLH = array(0, dim = c(length(LH.list)-1,3))
  # for (i in 1:no.total) {
  #   LHs = c(LH[i], LH_inter[i], LH_intra[i])
  #   for(j in 1:3){
  #     idiff = LHs[j]
  #     ibin = which(idiff<=LH.list)[1]-1
  #     freqLH[ibin,j] = freqLH[ibin,j]+1
  #   }
  # }
  # freqLH = freqLH/no.total
  
  ### plot histograms
  png(figname, width=8,height=10, res=300,unit="in")
  par(mfrow=c(3,2))
  plot(0,0,type="n",xlim=diff.limits,ylim=c(0,1/diff.step),
       xlab="Total Residual",ylab="Freq")
  p1 <- hist(diff,breaks = diff.list,  
             plot=FALSE
             )
  plot(p1, freq=FALSE,add=TRUE)
  box()
  
  plot(0,0,type="n",xlim=LH.limits,ylim=c(0,1/LH.step),
       xlab="Total LH",ylab="Freq")
  p1 <- hist(LH, breaks = LH.list, #freq=FALSE, 
             plot=FALSE)
  plot(p1, freq=FALSE,add=TRUE)
  box()
  
  if(iflag>0) {
    plot(0,0,type="n",xlim=diff.limits,ylim=c(0,1/diff.step),
         xlab="Inter-Residual",ylab="Freq")
    p1 <- hist(Inter_residl,breaks = diff.list, plot=FALSE)
    plot(p1, freq=FALSE, add=TRUE)
    box()
    
    plot(0,0,type="n",xlim=LH.limits,ylim=c(0,1/LH.step),
         xlab="LH for Inter-residual",ylab="Freq")
    p1 <- hist(LH_inter, breaks = LH.list, 
               plot=FALSE)
    plot(p1,freq=FALSE,  add=TRUE)
    box()
    
    plot(0,0,type="n",xlim=diff.limits,ylim=c(0,1/diff.step),
         xlab="Intra-Residual",ylab="Freq")
    
    p1 <- hist(Intra_residl, breaks = diff.list, plot=FALSE)
    plot(p1, freq=FALSE, add=TRUE)
    box()
    
    plot(0,0,type="n",xlim=LH.limits,ylim=c(0,1/LH.step),
         xlab="LH for Intra-residual",ylab="Freq")
    p1 <- hist(LH_intra, breaks = LH.list,  
               plot=FALSE)
    plot(p1, freq=FALSE, add=TRUE)
    box()
  }
  dev.off()
}


plot.xperiod <- function(iflag, figname.base, sum.diff, sum.inter, sum.stats,
                         list.gmpe, list.p,
                         diff.limits = c(-1,1)){
  
  
  no.gmpe = length(list.gmpe)
  xlim0 = c(min(list.p),max(list.p))
  i=1
  for (i in 1:no.gmpe) {
    # for (i in 1:1) {
    sum.diff.sub = sum.diff[sum.diff$gmpe == list.gmpe[i], ]
    sum.inter.sub = sum.inter[sum.inter$gmpe == list.gmpe[i], ]
    sum.stats.sub = sum.stats[sum.stats[,8]== list.gmpe[i],]
    
    figname = sprintf("scat_%s_%s.png",figname.base,list.gmpe[i])

    png(figname, width=8,height=10, res=300,unit="in")
    par(mfrow=c(3,2))
    
    plot(0.01,0,type="n",xlim=xlim0,ylim=diff.limits,
         log="x", #xaxs="i", yaxs="i",
         xlab="Period",ylab="Total Residual")
    temp = sum.diff.sub$diff
    points(sum.diff.sub$period, temp)
    lines(c(0.0001,10),c(0,0))
    lines(sum.stats.sub[,7],sum.stats.sub[,1],col="red")
    txt = sprintf("max=%s\nmin=%s",round(max(temp),1),
                  round(min(temp),1))
    text(0.02,.9,labels=txt)
    box()
    
    plot(0.01,0,type="n",xlim=xlim0,ylim=c(0,1),
         log="x", #xaxs="i", yaxs="i",
         xlab="Period",ylab="Total LH")
    temp = sum.diff.sub$LH
    points(sum.diff.sub$period,temp )
    # line(c(0.0001,10),c(0,0))
    lines(sum.stats.sub[,7],sum.stats.sub[,2],col="red")
    txt = sprintf("max=%s\nmin=%s",round(max(temp),1),
                  round(min(temp),1))
    text(0.02,.9,labels=txt)
    box()
    
    
    if(iflag[i]>0) {
      plot(0.01,0,type="n",xlim=xlim0,ylim=diff.limits,
           log="x", #xaxs="i", yaxs="i",
           xlab="Period",ylab="Inter-Residual")
      temp = sum.inter.sub$Inter_residl
      points(sum.inter.sub$period, temp)
      lines(c(0.0001,10),c(0,0))
      lines(sum.stats.sub[,7],sum.stats.sub[,3],col="red")
      txt = sprintf("max=%s\nmin=%s",round(max(temp),1),
                    round(min(temp),1))
      text(0.02,.9,labels=txt)
      box()
      
      plot(0.01,0,type="n",xlim=xlim0,ylim=c(0,1),
           log="x", #xaxs="i", yaxs="i",
           xlab="Period",ylab="LH for Inter-residual")
      temp = sum.inter.sub$LH_inter
      points(sum.inter.sub$period,temp )
      # line(c(0.0001,10),c(0,0))
      lines(sum.stats.sub[,7],sum.stats.sub[,4],col="red")
      txt = sprintf("max=%s\nmin=%s",round(max(temp),1),
                    round(min(temp),1))
      text(0.02,.9,labels=txt)
      box()
      
      plot(0.01,0,type="n",xlim=xlim0,ylim=diff.limits,
           log="x", #xaxs="i", yaxs="i",
           xlab="Period",ylab="Intra-Residual")
      temp = sum.diff.sub$Intra_residl
      points(sum.diff.sub$period, temp)
      lines(c(0.0001,10),c(0,0))
      lines(sum.stats.sub[,7],sum.stats.sub[,5],col="red")
      txt = sprintf("max=%s\nmin=%s",round(max(temp),1),
                    round(min(temp),1))
      text(0.02,.9,labels=txt)
      box()
      
      plot(0.01,0,type="n",xlim=xlim0,ylim=c(0,1),
           log="x", #xaxs="i", yaxs="i",
           xlab="Period",ylab="LH for Intra-residual")
      temp = sum.diff.sub$LH_intra
      points(sum.diff.sub$period,temp )
      # line(c(0.0001,10),c(0,0))
      lines(sum.stats.sub[,7],sum.stats.sub[,6],col="red")
      txt = sprintf("max=%s\nmin=%s",round(max(temp),1),
                    round(min(temp),1))
      text(0.02,.9,labels=txt)
      box()
    }
    
    
    dev.off()
    
    
  }
  
}
