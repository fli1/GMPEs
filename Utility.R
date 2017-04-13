# library(grid)
library(colorspace)
# library(ggplot2)
# library(gridExtra)
# library(reshape2)
# library(scales)


### get the data for plotting purpose


GetSubdata <- function(output.sa, list.p,
                       flag.scenario="distance", flag.period = "PGA") {
  
  # output.Sa <- array(NA, dim = c(n, m, le, he)) # distance, magnitude, period, depth
  
  # list.dist = exp(seq(log(1),log(500),length.out=50)) ##
  # list.mag = 8  ##c(7.5,6.5)
  # list.zh = 50
  # 
  if(flag.period == "all"){
    ind.p = 0
  }else if(flag.period == "PGA"){
    ind.p = 1
  } else{
    # flag.period=0.2
    ind.p = which(list.p ==flag.period)
  }
  
  if(flag.scenario=="period") {  ## by default, read the results for the 1st magnitude and depth
    out = output.sa[1,1,,1]
  }
  
  if(flag.scenario=="distance") {  ## by default, read the results for the 1st magnitude and depth
    out = output.sa[,1,ind.p,1]
  }
  if(flag.scenario=="magnitude") {  ## by default, read the results for the 1st magnitude and depth
    out = output.sa[1,,ind.p,1]
  }
  if(flag.scenario=="depth") {  ## by default, read the results for the 1st magnitude and depth
    out = output.sa[1,1,ind.p,]
  }
  return(out)
}

fmt <- function(){
  f <- function(x) as.character(round(x,2))
  f
}

gg_color_hue <- function(n) {
  hues <- seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}

### plot figures

plot.comp <- function (figname, datain, data.recorded, 
                       list.GMPEs, flag.scenario, flag.period, 
                       flag.source,vs30,
                       list.mag, list.dist, list.zh, ymax){
  
  n=length(list.GMPEs)
  pal <-rainbow_hcl(n, c=100, l=50, start = 0, end = 360*(n-1)/n)
  
  lwd <- c(rep(c(2,3),n))
  lty <- c(1,1,1,1,2,2,2) #c(rep(c(1,2),n))
  col <- pal
  
  flag.log = ""
  figtitle = ""
  if(flag.scenario == "period") {
    xticks <- c(0.01, 0.02,0.05,0.1, 0.2, 0.5, 1, 2,5,10)
    ygridlines <- c(seq(0.01,0.1,0.01), 
                    seq(0.2, 1, 0.1), seq(2,10,1)) #c(seq(0.02,0.1, 0.01), seq(0.2,1,0.1), 2)
    flag.log="x"
    xlim0 = c(0.01, 10)
    xlabel = "Spectral Period (second)"
    figtitle = sprintf("%sMw %s, Distance %skm, Depth %skm, Vs30 %sm/s",figtitle, list.mag[1], list.dist[1], list.zh[1],vs30)
  }
  if(flag.scenario == "distance") {
    dist.max = max(list.dist)
    dist.min = 10 #min(list.dist)
    xticks <- c(seq(10,100,10),seq(200,1000,100)) #seq(dist.min, dist.max, 50)
    ygridlines <- c(seq(10,100,10),seq(200,1000,100)) #seq(dist.min, dist.max, 50)
    flag.log="x"
    xlim0 = c(dist.min, dist.max)
    
    xlabel = "Rupture Distance (km)"
    
    if(flag.source=="inslab") {
      xlabel = "Hypocentral Distance (km)"
    }
    
    figtitle = sprintf("%sPeriod %s, Mw %s, Depth %skm, Vs30 %sm/s",figtitle, flag.period, list.mag[1], list.zh[1],vs30)
  }
  if(flag.scenario == "magnitude") {
    temp.max = max(list.mag)
    temp.min = min(list.mag)
    xticks <- seq(temp.min, temp.max, 0.5)
    ygridlines <- seq(temp.min, temp.max, .5)
    xlim0 = c(temp.min, temp.max)
    xlabel = "Magnitude (mw)"
    figtitle = sprintf("%sPeriod %s, Distance %skm, Depth %skm, Vs30 %sm/s",figtitle, flag.period, list.dist[1], list.zh[1],vs30)
  }
  if(flag.scenario == "depth") {
    temp.max = max(list.zh)
    temp.min = min(list.zh)
    xticks <- seq(temp.min, temp.max, 10)
    ygridlines <- seq(temp.min, temp.max, 10)
    xlim0 = c(dist.min, dist.max)
    xlabel = "Hypocentral Depth (km)"
    figtitle = sprintf("%sPeriod %s, Mw %s, Distance %skm, Vs30 %sm/s",figtitle, flag.period, list.mag[1], list.dist[1],vs30)
  }
  xlabels <- as.character(xticks)
  
  y1ticks <- seq(0,3,0.2) #c(1e-04, 1e-03, 1e-02, 1e-01)
  y1labels <- as.character(y1ticks)
  xgridlines <- seq(0,3, 0.2) 
  
  png(file = sprintf("%s.png", figname), width=7, height=6,
      units = 'in', res = 400)
  
  par(oma = c(0.5,0.5,0,0.3))
  # par(mfrow =c (1,2))
  
  #   par(new=TRUE)
  plot(c(0), c(0), type = 'n', xaxs = 'i', yaxs = 'i', 
       log = flag.log, 
       xlim = xlim0, 
       ylim = c(0, ymax), axes = FALSE, ann = FALSE)
  axis(1, at = xticks, label = xlabels, cex = 0.55)
  axis(2, at = y1ticks, label = y1labels, cex = 0.55)
  
  abline(v=ygridlines, col = 'grey')
  abline(h=xgridlines, col = 'grey')
  
  mtext(xlabel, side = 1, line = 2.5, cex=0.95)
  mtext("Spectral Acceleration (g)", side = 2, line = 2.5, cex=0.95)
  mtext(sprintf("%s", figtitle),
        side = 3, line = 1, cex=0.95)
  
  legend("topright", list.GMPEs,
         col = col,
         ncol = 1, bty = 'o', border = 'white', bg = 'white',
         box.col = 'white', cex = 0.75, xjust = 1, yjust = 1,
         lwd = lwd, lty=lty)
  for(j in 1:length(list.GMPEs))
  {
    lines(datain[,1], datain[,j+1], col = col[j], lwd = lwd[j], lty=lty[j])
  }
  if(nrow(data.recorded)>1) {
    data.recorded = data.recorded[order(data.recorded[,3]),]
    points(data.recorded$variable, data.recorded$value,col=data.recorded$group)   #heat.colors(2),pch=16)#
  }
  
  box()
  dev.off()
  graphics.off()
  
}
