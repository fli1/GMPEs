library(scales)

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
  i=1
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
  # for(i in 1:10) {
    sa_observed = observed.sub$sa[i]
    sa_predicted = predicted.sub$sa[i]
    # 'sa_observed = Exp(sa_observed)
    eventID = predicted.sub$EQID[i]
    if (eventID == prevID) {
      temp = temp + sa_observed - sa_predicted
    } else {
      edID = i-1
      # print(stID)
      # print(edID)
      for (j in stID: edID){
        out.sum = c(out.sum, temp)
        out.no = c(out.no,edID - stID+1)
      }
      stID = edID+1
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
         diff.limits=c(-4,4),
         LH.limits=c(0,1)) {
  
  diff =totalresidual[[1]]
  LH = totalresidual[[2]]
  Inter_residl=totalresidual[[3]]
  LH_inter = totalresidual[[4]]
  Intra_residl = totalresidual[[5]] 
  LH_intra = totalresidual[[6]]
  iflag = totalresidual[[7]]
  
  
  diff.step = 0.5
  LH.step = 0.1
  
  diff.ylim = c(0,.8)
  LH.ylim = c(0,4)
  
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
  plot(0,0,type="n",xlim=diff.limits,ylim=diff.ylim,
       xlab="Total Residual",ylab="Freq")
  p1 <- hist(diff,breaks = diff.list,  
             plot=FALSE
             )
  plot(p1, freq=FALSE,add=TRUE)
  box()
  
  plot(0,0,type="n",xlim=LH.limits,ylim=LH.ylim,
       xlab="Total LH",ylab="Freq")
  p1 <- hist(LH, breaks = LH.list, #freq=FALSE, 
             plot=FALSE)
  plot(p1, freq=FALSE,add=TRUE)
  box()
  
  if(iflag>0) {
    plot(0,0,type="n",xlim=diff.limits,ylim=diff.ylim,
         xlab="Inter-Residual",ylab="Freq")
    p1 <- hist(Inter_residl,breaks = diff.list, plot=FALSE)
    plot(p1, freq=FALSE, add=TRUE)
    box()
    
    plot(0,0,type="n",xlim=LH.limits,ylim=LH.ylim,
         xlab="LH for Inter-residual",ylab="Freq")
    p1 <- hist(LH_inter, breaks = LH.list, 
               plot=FALSE)
    plot(p1,freq=FALSE,  add=TRUE)
    box()
    
    plot(0,0,type="n",xlim=diff.limits,ylim=diff.ylim,
         xlab="Intra-Residual",ylab="Freq")
    
    p1 <- hist(Intra_residl, breaks = diff.list, plot=FALSE)
    plot(p1, freq=FALSE, add=TRUE)
    box()
    
    plot(0,0,type="n",xlim=LH.limits,ylim=LH.ylim,
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

plot.statsB <- function(iflag, figname.base, sum.stats,
                         list.gmpe, list.p,
                         diff.limits = c(-1.5,1.5)){
  
  no.gmpe = length(list.gmpe)
  xlim0 = c(min(list.p),max(list.p))
  i=1
 
  figname = sprintf("stats_total_%s.png",figname.base)
  
  pal <-rainbow_hcl(no.gmpe, c=100, l=50, start = 0, end = 360*(no.gmpe-1)/no.gmpe)
  
  list.col = c(pal)
  list.lty = c(rep(1:2,10))
  list.lwd = c(rep(c(1,1,2,2),3)) ##c(rep(c(1,3),each=3))
  
  png(figname, width=10,height=8, res=300,unit="in")
  par(mfrow=c(2,2))
  
  plot(0.01,0,type="n",xlim=xlim0,ylim=diff.limits,
       log="x", #xaxs="i", yaxs="i",
       xlab="Period",ylab="Mean Total Residual")
  
  alpha0= seq(0.1,0.9,0.1) #c(0.2,0.4,0.6,0.8)
  cols = rep("black",4) #c("green", "yellow","orange","red")
  rect(0.0001,-0.25,10,0.25,col=alpha(cols[1],alpha0[1]),border=NA)
  rect(0.0001,-0.5,10,-0.25,col=alpha(cols[2],alpha0[2]),border=NA)
  rect(0.0001,0.25,10,0.5,col=alpha(cols[2],alpha0[2]),border=NA)
  rect(0.0001,-0.75,10,-0.5,col=alpha(cols[3],alpha0[3]),border=NA)
  rect(0.0001,0.5,10,0.75,col=alpha(cols[3],alpha0[3]),border=NA)
  rect(0.0001,-3,10,-0.75,col=alpha(cols[4],alpha0[4]),border=NA)
  rect(0.0001,0.75,10,3,col=alpha(cols[4],alpha0[4]),border=NA)
  abline(h=c(0))
  
  for (i in 1:no.gmpe) {
    sum.stats.sub = sum.stats[sum.stats[,8]== list.gmpe[i],]
    lines(sum.stats.sub[,7],sum.stats.sub[,1],col=list.col[i],
          lty = list.lty[i],
          lwd = list.lwd[i])  
  }
  # txt = sprintf("max=%s\nmin=%s",round(max(temp),1),
  #               round(min(temp),1))
  # text(0.02,.9,labels=txt)
  # abline(h=seq(-0.75,0.75,0.25))
  box()
  
  plot(0.01,0,type="n",xlim=xlim0,ylim=c(0.1,.9),
       log="x", #xaxs="i", yaxs="i",
       xlab="Period",ylab="Median Total LH")
  # lines(c(0.0001,10),c(0.5,0.5))
  # abline(h=seq(0.2,0.4,0.1))
  rect(0.0001,0.4,10,1,col=alpha(cols[1],alpha0[1]),border=NA)
  rect(0.0001,0.3,10,0.4,col=alpha(cols[2],alpha0[2]),border=NA)
  rect(0.0001,0.2,10,0.3,col=alpha(cols[3],alpha0[3]),border=NA)
  rect(0.0001,0,10,.2,col=alpha(cols[4],alpha0[4]),border=NA)
  for (i in 1:no.gmpe) {
    sum.stats.sub = sum.stats[sum.stats[,8]== list.gmpe[i],]
    lines(sum.stats.sub[,7],sum.stats.sub[,2],col=list.col[i],
          lty = list.lty[i],
          lwd = list.lwd[i])  
  }
  legend('topright',list.gmpe,
         ncol=3,
         col=list.col, lty=list.lty,
         lwd=list.lwd)
  
  box()
  
  plot(0.01,0,type="n",xlim=xlim0,ylim=diff.limits,
       log="x", #xaxs="i", yaxs="i",
       xlab="Period",ylab="Median Total Residual")

  rect(0.0001,-0.25,10,0.25,col=alpha(cols[1],alpha0[1]),border=NA)
  rect(0.0001,-0.5,10,-0.25,col=alpha(cols[2],alpha0[2]),border=NA)
  rect(0.0001,0.25,10,0.5,col=alpha(cols[2],alpha0[2]),border=NA)
  rect(0.0001,-0.75,10,-0.5,col=alpha(cols[3],alpha0[3]),border=NA)
  rect(0.0001,0.5,10,0.75,col=alpha(cols[3],alpha0[3]),border=NA)
  rect(0.0001,-3,10,-0.75,col=alpha(cols[4],alpha0[4]),border=NA)
  rect(0.0001,0.75,10,3,col=alpha(cols[4],alpha0[4]),border=NA)
  abline(h=c(0))
  
  for (i in 1:no.gmpe) {
    sum.stats.sub = sum.stats[sum.stats[,8]== list.gmpe[i],]
    lines(sum.stats.sub[,7],sum.stats.sub[,9],col=list.col[i],
          lty = list.lty[i],
          lwd = list.lwd[i])  
  }
  # txt = sprintf("max=%s\nmin=%s",round(max(temp),1),
  #               round(min(temp),1))
  # text(0.02,.9,labels=txt)
  # abline(h=seq(-0.75,0.75,0.25))
  box()
  
  plot(0.01,0,type="n",xlim=xlim0,ylim=c(0.5,2),
       log="x", #xaxs="i", yaxs="i",
       xlab="Period",ylab="StDev of Total Residual")

  rect(0.0001,0,10,1.125,col=alpha(cols[1],alpha0[1]),border=NA)
  rect(0.0001,1.125,10,1.25,col=alpha(cols[2],alpha0[2]),border=NA)
  rect(0.0001,1.25,10,1.5,col=alpha(cols[3],alpha0[3]),border=NA)
  rect(0.0001,1.5,10,3,col=alpha(cols[4],alpha0[4]),border=NA)
  abline(h=c(1))
  
  for (i in 1:no.gmpe) {
    sum.stats.sub = sum.stats[sum.stats[,8]== list.gmpe[i],]
    lines(sum.stats.sub[,7],sum.stats.sub[,12],col=list.col[i],
          lty = list.lty[i],
          lwd = list.lwd[i])  
  }
  # txt = sprintf("max=%s\nmin=%s",round(max(temp),1),
  #               round(min(temp),1))
  # text(0.02,.9,labels=txt)
  # abline(h=seq(-0.75,0.75,0.25))
  box()
  
  dev.off()
  
  ###-------- for inter-residual
  figname = sprintf("stats_inter_%s.png",figname.base)

  png(figname, width=10,height=8, res=300,unit="in")
  par(mfrow=c(2,2))
  plot(0.01,0,type="n",xlim=xlim0,ylim=diff.limits,
       log="x", #xaxs="i", yaxs="i",
       xlab="Period",ylab="Mean Inter-Residual")
  # abline(h=seq(-0.75,0.75,0.25))
  rect(0.0001,-0.25,10,0.25,col=alpha(cols[1],alpha0[1]),border=NA)
  rect(0.0001,-0.5,10,-0.25,col=alpha(cols[2],alpha0[2]),border=NA)
  rect(0.0001,0.25,10,0.5,col=alpha(cols[2],alpha0[2]),border=NA)
  rect(0.0001,-0.75,10,-0.5,col=alpha(cols[3],alpha0[3]),border=NA)
  rect(0.0001,0.5,10,0.75,col=alpha(cols[3],alpha0[3]),border=NA)
  rect(0.0001,-3,10,-0.75,col=alpha(cols[4],alpha0[4]),border=NA)
  rect(0.0001,0.75,10,3,col=alpha(cols[4],alpha0[4]),border=NA)
  abline(h=c(0))
  for (i in 1:no.gmpe) {
    if(iflag[i]>0) {
      sum.stats.sub = sum.stats[sum.stats[,8]== list.gmpe[i],]
      lines(sum.stats.sub[,7],sum.stats.sub[,3],col=list.col[i],
            lty = list.lty[i],
            lwd = list.lwd[i])  
    }
  }
  box()
  
  plot(0.01,0,type="n",xlim=xlim0,ylim=c(0.1,.9),
       log="x", #xaxs="i", yaxs="i",
       xlab="Period",ylab="Median LH for Inter-residual")
  rect(0.0001,0.4,10,1,col=alpha(cols[1],alpha0[1]),border=NA)
  rect(0.0001,0.3,10,0.4,col=alpha(cols[2],alpha0[2]),border=NA)
  rect(0.0001,0.2,10,0.3,col=alpha(cols[3],alpha0[3]),border=NA)
  rect(0.0001,0,10,.2,col=alpha(cols[4],alpha0[4]),border=NA)
  
  for (i in 1:no.gmpe) {
    if(iflag[i]>0) {
      sum.stats.sub = sum.stats[sum.stats[,8]== list.gmpe[i],]
      lines(sum.stats.sub[,7],sum.stats.sub[,4],col=list.col[i],
            lty = list.lty[i],
            lwd = list.lwd[i])  
    }
  }
  legend('topright',list.gmpe,
         ncol=3,
         col=list.col, lty=list.lty,
         lwd=list.lwd)
  box()
  
  plot(0.01,0,type="n",xlim=xlim0,ylim=diff.limits,
       log="x", #xaxs="i", yaxs="i",
       xlab="Period",ylab="Median Inter-Residual")
  # abline(h=seq(-0.75,0.75,0.25))
  rect(0.0001,-0.25,10,0.25,col=alpha(cols[1],alpha0[1]),border=NA)
  rect(0.0001,-0.5,10,-0.25,col=alpha(cols[2],alpha0[2]),border=NA)
  rect(0.0001,0.25,10,0.5,col=alpha(cols[2],alpha0[2]),border=NA)
  rect(0.0001,-0.75,10,-0.5,col=alpha(cols[3],alpha0[3]),border=NA)
  rect(0.0001,0.5,10,0.75,col=alpha(cols[3],alpha0[3]),border=NA)
  rect(0.0001,-3,10,-0.75,col=alpha(cols[4],alpha0[4]),border=NA)
  rect(0.0001,0.75,10,3,col=alpha(cols[4],alpha0[4]),border=NA)
  abline(h=c(0))
  for (i in 1:no.gmpe) {
    if(iflag[i]>0) {
      sum.stats.sub = sum.stats[sum.stats[,8]== list.gmpe[i],]
      lines(sum.stats.sub[,7],sum.stats.sub[,10],col=list.col[i],
            lty = list.lty[i],
            lwd = list.lwd[i])  
    }
  }
  box()
  
  plot(0.01,0,type="n",xlim=xlim0,ylim=c(0.5,2),
       log="x", #xaxs="i", yaxs="i",
       xlab="Period",ylab="Stdev of Inter-Residual")
  # abline(h=seq(-0.75,0.75,0.25))
  rect(0.0001,0,10,1.125,col=alpha(cols[1],alpha0[1]),border=NA)
  rect(0.0001,1.125,10,1.25,col=alpha(cols[2],alpha0[2]),border=NA)
  rect(0.0001,1.25,10,1.5,col=alpha(cols[3],alpha0[3]),border=NA)
  rect(0.0001,1.5,10,3,col=alpha(cols[4],alpha0[4]),border=NA)
  abline(h=c(1))
  for (i in 1:no.gmpe) {
    if(iflag[i]>0) {
      sum.stats.sub = sum.stats[sum.stats[,8]== list.gmpe[i],]
      lines(sum.stats.sub[,7],sum.stats.sub[,13],col=list.col[i],
            lty = list.lty[i],
            lwd = list.lwd[i])  
    }
  }
  box()
  dev.off()
  
  ##---------- for intra 
  figname = sprintf("stats_intra_%s.png",figname.base)
  
  png(figname, width=10,height=8, res=300,unit="in")
  par(mfrow=c(2,2))
  
  plot(0.01,0,type="n",xlim=xlim0,ylim=diff.limits,
       log="x", #xaxs="i", yaxs="i",
       xlab="Period",ylab="Mean Intra-Residual")
  rect(0.0001,-0.25,10,0.25,col=alpha(cols[1],alpha0[1]),border=NA)
  rect(0.0001,-0.5,10,-0.25,col=alpha(cols[2],alpha0[2]),border=NA)
  rect(0.0001,0.25,10,0.5,col=alpha(cols[2],alpha0[2]),border=NA)
  rect(0.0001,-0.75,10,-0.5,col=alpha(cols[3],alpha0[3]),border=NA)
  rect(0.0001,0.5,10,0.75,col=alpha(cols[3],alpha0[3]),border=NA)
  rect(0.0001,-3,10,-0.75,col=alpha(cols[4],alpha0[4]),border=NA)
  rect(0.0001,0.75,10,3,col=alpha(cols[4],alpha0[4]),border=NA)
  abline(h=c(0))
  
  for (i in 1:no.gmpe) {
    if(iflag[i]>0) {
      sum.stats.sub = sum.stats[sum.stats[,8]== list.gmpe[i],]
      lines(sum.stats.sub[,7],sum.stats.sub[,5],col=list.col[i],
            lty = list.lty[i],
            lwd = list.lwd[i])  
    
    }
  }
  
  # abline(h=seq(-0.75,0.75,0.25))
  box()
  
  plot(0.01,0,type="n",xlim=xlim0,ylim=c(0.1,.9),
       log="x", #xaxs="i", yaxs="i",
       xlab="Period",ylab="Median LH for Intra-residual")
  rect(0.0001,0.4,10,1,col=alpha(cols[1],alpha0[1]),border=NA)
  rect(0.0001,0.3,10,0.4,col=alpha(cols[2],alpha0[2]),border=NA)
  rect(0.0001,0.2,10,0.3,col=alpha(cols[3],alpha0[3]),border=NA)
  rect(0.0001,0,10,.2,col=alpha(cols[4],alpha0[4]),border=NA)
  
  for (i in 1:no.gmpe) {
    if(iflag[i]>0) {
      sum.stats.sub = sum.stats[sum.stats[,8]== list.gmpe[i],]
      lines(sum.stats.sub[,7],sum.stats.sub[,6],col=list.col[i],
            lty = list.lty[i],
            lwd = list.lwd[i])  
    }
  }
  legend('topright',list.gmpe,
         ncol=3,
         col=list.col, lty=list.lty,
         lwd=list.lwd)
  
  box()
  
  plot(0.01,0,type="n",xlim=xlim0,ylim=diff.limits,
       log="x", #xaxs="i", yaxs="i",
       xlab="Period",ylab="Median Intra-Residual")
  rect(0.0001,-0.25,10,0.25,col=alpha(cols[1],alpha0[1]),border=NA)
  rect(0.0001,-0.5,10,-0.25,col=alpha(cols[2],alpha0[2]),border=NA)
  rect(0.0001,0.25,10,0.5,col=alpha(cols[2],alpha0[2]),border=NA)
  rect(0.0001,-0.75,10,-0.5,col=alpha(cols[3],alpha0[3]),border=NA)
  rect(0.0001,0.5,10,0.75,col=alpha(cols[3],alpha0[3]),border=NA)
  rect(0.0001,-3,10,-0.75,col=alpha(cols[4],alpha0[4]),border=NA)
  rect(0.0001,0.75,10,3,col=alpha(cols[4],alpha0[4]),border=NA)
  abline(h=c(0))
  
  for (i in 1:no.gmpe) {
    if(iflag[i]>0) {
      sum.stats.sub = sum.stats[sum.stats[,8]== list.gmpe[i],]
      lines(sum.stats.sub[,7],sum.stats.sub[,11],col=list.col[i],
            lty = list.lty[i],
            lwd = list.lwd[i])  
      
    }
  }
  
  # abline(h=seq(-0.75,0.75,0.25))
  box()
  
  plot(0.01,0,type="n",xlim=xlim0,ylim=c(0.5,2),
       log="x", #xaxs="i", yaxs="i",
       xlab="Period",ylab="StDev of Intra-Residual")
  rect(0.0001,0,10,1.125,col=alpha(cols[1],alpha0[1]),border=NA)
  rect(0.0001,1.125,10,1.25,col=alpha(cols[2],alpha0[2]),border=NA)
  rect(0.0001,1.25,10,1.5,col=alpha(cols[3],alpha0[3]),border=NA)
  rect(0.0001,1.5,10,3,col=alpha(cols[4],alpha0[4]),border=NA)
  abline(h=c(1))
  
  for (i in 1:no.gmpe) {
    if(iflag[i]>0) {
      sum.stats.sub = sum.stats[sum.stats[,8]== list.gmpe[i],]
      lines(sum.stats.sub[,7],sum.stats.sub[,14],col=list.col[i],
            lty = list.lty[i],
            lwd = list.lwd[i])  
      
    }
  }
  
  # abline(h=seq(-0.75,0.75,0.25))
  box()
  dev.off()
  
}

plot.stats <- function(iflag, figname.base, sum.stats,
                       list.gmpe, list.p,
                       diff.limits = c(-1.5,1.5)){
  
  no.gmpe = length(list.gmpe)
  xlim0 = c(min(list.p),max(list.p))
  i=1
  
  figname = sprintf("stats_%s.png",figname.base)
  
  pal <-rainbow_hcl(no.gmpe, c=100, l=50, start = 0, end = 360*(no.gmpe-1)/no.gmpe)
  
  list.col = c(pal)
  list.lty = c(rep(1:2,10))
  list.lwd = c(rep(c(1,1,2,2),3)) ##c(rep(c(1,3),each=3))
  
  png(figname, width=8,height=9, res=300,unit="in")
  par(mfrow=c(3,2))
  
  plot(0.01,0,type="n",xlim=xlim0,ylim=diff.limits,
       log="x", #xaxs="i", yaxs="i",
       xlab="Period",ylab="Mean Total Residual")
  
  alpha0= seq(0.1,0.9,0.1) #c(0.2,0.4,0.6,0.8)
  cols = rep("black",4) #c("green", "yellow","orange","red")
  rect(0.0001,-0.25,10,0.25,col=alpha(cols[1],alpha0[1]),border=NA)
  rect(0.0001,-0.5,10,-0.25,col=alpha(cols[2],alpha0[2]),border=NA)
  rect(0.0001,0.25,10,0.5,col=alpha(cols[2],alpha0[2]),border=NA)
  rect(0.0001,-0.75,10,-0.5,col=alpha(cols[3],alpha0[3]),border=NA)
  rect(0.0001,0.5,10,0.75,col=alpha(cols[3],alpha0[3]),border=NA)
  rect(0.0001,-3,10,-0.75,col=alpha(cols[4],alpha0[4]),border=NA)
  rect(0.0001,0.75,10,3,col=alpha(cols[4],alpha0[4]),border=NA)
  abline(h=c(0))
  
  for (i in 1:no.gmpe) {
    sum.stats.sub = sum.stats[sum.stats[,8]== list.gmpe[i],]
    lines(sum.stats.sub[,7],sum.stats.sub[,1],col=list.col[i],
          lty = list.lty[i],
          lwd = list.lwd[i])  
  }
  # txt = sprintf("max=%s\nmin=%s",round(max(temp),1),
  #               round(min(temp),1))
  # text(0.02,.9,labels=txt)
  
  # abline(h=seq(-0.75,0.75,0.25))
  
  box()
  
  plot(0.01,0,type="n",xlim=xlim0,ylim=c(0.1,.9),
       log="x", #xaxs="i", yaxs="i",
       xlab="Period",ylab="Median Total LH")
  # lines(c(0.0001,10),c(0.5,0.5))
  # abline(h=seq(0.2,0.4,0.1))
  rect(0.0001,0.4,10,1,col=alpha(cols[1],alpha0[1]),border=NA)
  rect(0.0001,0.3,10,0.4,col=alpha(cols[2],alpha0[2]),border=NA)
  rect(0.0001,0.2,10,0.3,col=alpha(cols[3],alpha0[3]),border=NA)
  rect(0.0001,0,10,.2,col=alpha(cols[4],alpha0[4]),border=NA)
  
  
  for (i in 1:no.gmpe) {
    sum.stats.sub = sum.stats[sum.stats[,8]== list.gmpe[i],]
    lines(sum.stats.sub[,7],sum.stats.sub[,2],col=list.col[i],
          lty = list.lty[i],
          lwd = list.lwd[i])  
  }
  legend('topright',list.gmpe,
         ncol=3,
         col=list.col, lty=list.lty,
         lwd=list.lwd)
  
  box()
  
  plot(0.01,0,type="n",xlim=xlim0,ylim=diff.limits,
       log="x", #xaxs="i", yaxs="i",
       xlab="Period",ylab="Mean Inter-Residual")
  # abline(h=seq(-0.75,0.75,0.25))
  rect(0.0001,-0.25,10,0.25,col=alpha(cols[1],alpha0[1]),border=NA)
  rect(0.0001,-0.5,10,-0.25,col=alpha(cols[2],alpha0[2]),border=NA)
  rect(0.0001,0.25,10,0.5,col=alpha(cols[2],alpha0[2]),border=NA)
  rect(0.0001,-0.75,10,-0.5,col=alpha(cols[3],alpha0[3]),border=NA)
  rect(0.0001,0.5,10,0.75,col=alpha(cols[3],alpha0[3]),border=NA)
  rect(0.0001,-3,10,-0.75,col=alpha(cols[4],alpha0[4]),border=NA)
  rect(0.0001,0.75,10,3,col=alpha(cols[4],alpha0[4]),border=NA)
  abline(h=c(0))
  for (i in 1:no.gmpe) {
    if(iflag[i]>0) {
      sum.stats.sub = sum.stats[sum.stats[,8]== list.gmpe[i],]
      lines(sum.stats.sub[,7],sum.stats.sub[,3],col=list.col[i],
            lty = list.lty[i],
            lwd = list.lwd[i])  
    }
  }
  
  box()
  
  plot(0.01,0,type="n",xlim=xlim0,ylim=c(0.1,.9),
       log="x", #xaxs="i", yaxs="i",
       xlab="Period",ylab="Median LH for Inter-residual")
  rect(0.0001,0.4,10,1,col=alpha(cols[1],alpha0[1]),border=NA)
  rect(0.0001,0.3,10,0.4,col=alpha(cols[2],alpha0[2]),border=NA)
  rect(0.0001,0.2,10,0.3,col=alpha(cols[3],alpha0[3]),border=NA)
  rect(0.0001,0,10,.2,col=alpha(cols[4],alpha0[4]),border=NA)
  
  for (i in 1:no.gmpe) {
    if(iflag[i]>0) {
      sum.stats.sub = sum.stats[sum.stats[,8]== list.gmpe[i],]
      lines(sum.stats.sub[,7],sum.stats.sub[,4],col=list.col[i],
            lty = list.lty[i],
            lwd = list.lwd[i])  
    }
  }
  box()
  
  plot(0.01,0,type="n",xlim=xlim0,ylim=diff.limits,
       log="x", #xaxs="i", yaxs="i",
       xlab="Period",ylab="Mean Intra-Residual")
  rect(0.0001,-0.25,10,0.25,col=alpha(cols[1],alpha0[1]),border=NA)
  rect(0.0001,-0.5,10,-0.25,col=alpha(cols[2],alpha0[2]),border=NA)
  rect(0.0001,0.25,10,0.5,col=alpha(cols[2],alpha0[2]),border=NA)
  rect(0.0001,-0.75,10,-0.5,col=alpha(cols[3],alpha0[3]),border=NA)
  rect(0.0001,0.5,10,0.75,col=alpha(cols[3],alpha0[3]),border=NA)
  rect(0.0001,-3,10,-0.75,col=alpha(cols[4],alpha0[4]),border=NA)
  rect(0.0001,0.75,10,3,col=alpha(cols[4],alpha0[4]),border=NA)
  abline(h=c(0))
  
  for (i in 1:no.gmpe) {
    if(iflag[i]>0) {
      sum.stats.sub = sum.stats[sum.stats[,8]== list.gmpe[i],]
      lines(sum.stats.sub[,7],sum.stats.sub[,5],col=list.col[i],
            lty = list.lty[i],
            lwd = list.lwd[i])  
      
    }
  }
  
  # abline(h=seq(-0.75,0.75,0.25))
  
  box()
  
  plot(0.01,0,type="n",xlim=xlim0,ylim=c(0.1,.9),
       log="x", #xaxs="i", yaxs="i",
       xlab="Period",ylab="Median LH for Intra-residual")
  rect(0.0001,0.4,10,1,col=alpha(cols[1],alpha0[1]),border=NA)
  rect(0.0001,0.3,10,0.4,col=alpha(cols[2],alpha0[2]),border=NA)
  rect(0.0001,0.2,10,0.3,col=alpha(cols[3],alpha0[3]),border=NA)
  rect(0.0001,0,10,.2,col=alpha(cols[4],alpha0[4]),border=NA)
  
  for (i in 1:no.gmpe) {
    if(iflag[i]>0) {
      sum.stats.sub = sum.stats[sum.stats[,8]== list.gmpe[i],]
      lines(sum.stats.sub[,7],sum.stats.sub[,6],col=list.col[i],
            lty = list.lty[i],
            lwd = list.lwd[i])  
    }
  }
  
  box()
  dev.off()
  
}


####
plothist4gmpes <- function(figname.base, sum.diff, 
                           list.gmpe, period = 0.01,
                     diff.limits=c(-4,4),
                     LH.limits=c(0,1)) {
  
  sum.diff = sum.diff[sum.diff$period==period,]
  no.gmpe = length(list.gmpe)
  diff.step = 0.5
  LH.step = 0.1
  
  diff.ylim = c(0,.8)
  LH.ylim = c(0,4)
  
  diff.list = seq(-10,10,diff.step)
  diff.list = diff.list[diff.list>=diff.limits[1] & diff.list<=diff.limits[2]]
  
  # diff.list = c(-10,diff.list,10)
  LH.list = seq(0,1,LH.step)
  
  sum.diff$diff[sum.diff$diff>diff.limits[2]] = diff.limits[2]
  sum.diff$diff[sum.diff$diff<diff.limits[1]] = diff.limits[1]
  
  ### plot histograms
  figname = sprintf("comp_%s_%s.png",figname.base,period)
  png(figname, width=10,height=10, res=300,unit="in")
  par(mfrow=c(4,4))
  
  ##plot diffs
  for (i in 1:no.gmpe) {
    diff = sum.diff[sum.diff$gmpe==list.gmpe[i],]$diff
    plot(0,0,type="n",xlim=diff.limits,ylim=diff.ylim,
         xlab="Total Residual",ylab="Freq", main=paste0(list.gmpe[i]," (",period,")"))
    p1 <- hist(diff,breaks = diff.list,  
               plot=FALSE)
    plot(p1, freq=FALSE,add=TRUE)
    box()
  }
  for (i in 1:no.gmpe) {
    LH = sum.diff[sum.diff$gmpe==list.gmpe[i],]$LH
    plot(0,0,type="n",xlim=LH.limits,ylim=LH.ylim,
         xlab="Total LH",ylab="Freq", main=paste0(list.gmpe[i]," (",period,")"))
    p1 <- hist(LH, breaks = LH.list, #freq=FALSE, 
               plot=FALSE)
    plot(p1, freq=FALSE,add=TRUE)
    box()
  }

  dev.off()
}
