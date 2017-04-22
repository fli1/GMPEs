# rm(list = ls())
# library(nga)
# library(ggplot2)
# library(scales)
AB03SZ.subCal <- function (ip, M, h, Df, Zt, Vs30, Zl)
{
#   
  
#   % Purpose: Computes the mean and standard deviation of the PGA
#   % or psuedoacceleration, PSA, 5% damping. Additional modifications are also
#   % included for the regions of Cascadia and Japan.
#   
#   % Citation: "Empirical Ground-Motion Relations for Subduction-Zone
#   % Earthquakes and Their Application to Cascadia and Other Regions" by
#   % Atkinson, Gail M., David M. Boore. Bulletin of the Seismological Society
#   % of America, Vol. 93, No. 4, pp. 1703-1729, August 2003
#   
#   % General Limitations: This equation is obtained from a global database of
#   % ground motions (~1200 horizontal records). It is highly recommended that
#   % a comparison be drawn using records that compare well with the region of
#   % interest.
#   %   The authors suggest applying the equation for these given bounds:
#     
#     %   For Interface Events:
#     %   5.5 <= M < 6.5 and Df <= 80 km
#   %   6.5 <= M < 7.5 and Df <= 150 km
#   %   M >= 7.5 and Df <= 300 km
#   
#   %   For Intraslab (In-slab) Events:
#     %   6.0 <= M < 6.5 and Df <= 100 km
#   %   M >= 6.5 and Df <= 200 km
#   
#   %   The established bounds optimize the equations for seismic-hazard
#   %   analysis.
#   
#   %------------------------------INPUTS-------------------------------------%
#   
#   % T  = Period (sec)
#   % M  = Moment Magnitude
#   % h  = Focal (Hypocentral) Depth (km)
#   % Df = closest distance to the fault surface (km)
#   % Zt = Subduction Type Indicator: Zt = 0 for interface events
#   %                                 Zt = 1 for intraslab (in-slab) events
#   % Vs30 = Shear Wave Velocity averaged over the top 30 meters of soil of the
#   %      soil profile (m/sec)
#   % Zl = Cascadia or Japan indicator: Zl = 0 for General Cases
#   %                                   Zl = 1 for Cascadia
#   %                                   Zl = 2 for Japan
#   
#   %------------------------------OUTPUTS------------------------------------%
#     
#     % Sa    = Median spectral acceleration prediction (g)
#   % sigma = Logarithmic standard deviation of spectral acceleration
#   %         prediction
  
  #-------------------Period----------------------------
  period <- array(c(0.01, 0.04, 0.1, 0.2, 0.4, 1, 2, 1/0.33))
#   %---------------------Interslab Events Coefficients-----------------------%
    
  c1_it <- array(c(2.991, 2.8753, 2.7789, 2.6638, 2.5249, 2.1442, 2.1907, 2.301))
  c1_it_jp <- array(c(3.14, 3.05, 2.95, 2.84, 2.58, 2.18, 2.14, 2.27))
  c1_it_cas <- array(c(2.79, 2.60, 2.50, 2.54, 2.50, 2.18, 2.33, 2.36))
  c2_it <- array(c(0.03525, 0.07052, 0.09841, 0.12386, 0.1477, 0.1345, 0.07148, 0.02237))
  c3_it <- array(c(0.00759, 0.01004, 0.00974, 0.00884, 0.00728, 0.00521, 0.00224, 0.00012))
  c4_it <- array(c(-0.00206, -0.00278, -0.00287, -0.0028, -0.00235, -0.0011, 0, 0))
  c5_it <- array(c(0.19, 0.15, 0.15, 0.15, 0.13, 0.1, 0.1, 0.1))
  c6_it <- array(c(0.24, 0.2, 0.23, 0.27, 0.37, 0.3, 0.25, 0.25))
  c7_it <- array(c(0.29, 0.2, 0.2, 0.25, 0.38, 0.55, 0.4, 0.36))
  sigma_it <- array(c(0.23, 0.26, 0.27, 0.28, 0.29, 0.34, 0.34, 0.36))
  sigma1_it <- array(c(0.2, 0.22, 0.25, 0.25, 0.25, 0.28, 0.29, 0.31))
  sigma2_it <- array(c(0.11, 0.14, 0.1, 0.13, 0.15, 0.19, 0.18, 0.18))
  
#   %---------------------In-slab Events Coefficients-------------------------%
    
  c1_in <- array(c(-0.04713, 0.50697, 0.43928, 0.51589, 0.005445, -1.02133, -2.39234, -3.70012))
  c1_in_jp <- array(c(0.10, 0.68, 0.61, 0.70, 0.07, -0.98, -2.44, -3.73))
  c1_in_cas <- array(c(-0.25, 0.23, 0.16, 0.40, -0.01, -0.98, -2.25, -3.64))
  c2_in <- array(c(0.6909, 0.63273, 0.66675, 0.69186, 0.7727, 0.8789, 0.9964, 1.1169))
  c3_in <- array(c(0.0113, 0.01275, 0.0108, 0.00572, 0.00173, 0.0013, 0.00364, 0.00615))
  c4_in <- array(c(-0.00202, -0.00234, -0.00219, -0.00192, -0.00178, -0.00173, -0.00118, -0.00045))
  c5_in <- array(c(0.19, 0.15, 0.15, 0.15, 0.13, 0.1, 0.1, 0.1))
  c6_in <- array(c(0.24, 0.2, 0.23, 0.27, 0.37, 0.3, 0.25, 0.25))
  c7_in <- array(c(0.29, 0.2, 0.2, 0.25, 0.38, 0.55, 0.4 ,0.36))
  sigma_in <- array(c(0.27, 0.25, 0.28, 0.28, 0.28, 0.29, 0.3, 0.3))
  sigma1_in <- array(c(0.23, 0.24, 0.27, 0.26, 0.26, 0.27, 0.28, 0.29))
  sigma2_in <- array(c(0.14, 0.07, 0.07, 0.1, 0.1, 0.11, 0.11, 0.08))
  
  #Preliminary Inital Conditions and Variables Computation
  if (h >= 100){
    h <- 100
  }
    
  if (Zt == 0 & M >= 8.5){
    M <- 8.5  }
#   if (Zt == 0 & M < 8.5) {
#     M  }
  if (Zt == 1 & M >= 8.0) {
    M <- 8.0   }
#   if (Zt == 1 & M < 8.0) {
#     M   }
 
  
  delta <- 0.00724*(10^(0.507*M))
  R <- sqrt(Df^2 + delta^2)
  
  if (Zt == 0) {
    g <- 10^(1.2 - 0.18*M) }  
  if (Zt == 1) {
    g <- 10^(0.301 - 0.01*M)}
  
  if (Vs30 > 760) {
    Sc <- 0
    Sd <- 0
    Se <- 0 }
  if (Vs30 > 360 & Vs30 <= 760) {
    Sc <- 1
    Sd <- 0
    Se <- 0}
  if (Vs30 >= 180 & Vs30 <= 360 ) {
    Sc <- 0
    Sd <- 1
    Se <- 0}
  if (Vs30 < 180) {
    Sc <- 0
    Sd <- 0
    Se <- 1}
  
 ## % Begin Computation of Ground Motion Parameter with the modifications
  if (Zt == 0) {
    if (Zl == 0) {
      c1 <- c1_it[1]}
    if (Zl == 1) {
      c1 <- c1_it_cas[1]}
    if (Zl == 2) {
      c1 <- c1_it_jp[1]}
  
    log_PGArx <- c1 + c2_it[1]*M + c3_it[1]*h + c4_it[1]*R - g*log10(R);
    PGArx <- 10^(log_PGArx);
    
    if (PGArx <= 100 | (1/ip) <= 1) {
      sl <- 1}
    
    if ((1/ip) > 1 & (1/ip) < 2) {
      if (PGArx < 500 & PGArx > 100) {
        sl <- 1 - ((1/ip)-1)*(PGArx - 100)/400}
      if (PGArx >= 500) {
        sl <- 1 - ((1/ip)-1)}
    }
    
    if ((1/ip) >= 2) {
      if (PGArx < 500 & PGArx > 100) {
        sl <- 1 - (PGArx - 100)/400}
      if (PGArx >= 500) {
        sl <- 0}
    }
    
  }
  if (Zt == 1) {
    if (Zl == 0) {
      c1 <- c1_in[1]}
    if (Zl == 1) {
      c1 <- c1_in_cas[1]}
    if (Zl == 2) {
      c1 <- c1_in_jp[1]}
  
    log_PGArx <- c1 + c2_in[1]*M + c3_in[1]*h + c4_in[1]*R - g*log10(R);
    PGArx <- 10^(log_PGArx);
  
    if (PGArx <= 100 | (1/ip) <= 1) {
      sl <- 1}
  
    if ((1/ip) >= 2) {
      if (PGArx < 500 & PGArx > 100) {
        sl <- 1 - (PGArx - 100)/400}
      if (PGArx >= 500) {
        sl <- 0}
    }
    if ((1/ip) < 2 & (1/ip) > 1) {
      if (PGArx < 500 & PGArx > 100) {
        sl <- 1 - ((1/ip)-1)*(PGArx - 100)/400}
      if (PGArx >= 500) {
        sl <- 1 - ((1/ip)-1)}
    }
  }
  
  if (Zt == 0) {
#     if (length(which(period == T)) == 0) {
#       i_lo <- sum(period<T)
#       T_lo <- period[i_lo]
#       T_hi <- period[i_lo + 1]
#       Sa_sigma_hi <- AB_2003_SZ(T_hi,M,h,Df,Zt,Vs30,Zl)
#       Sa_sigma_lo <- AB_2003_SZ(T_lo,M,h,Df,Zt,Vs30,Zl)
#       
#       x <- array(c(T_lo, T_hi))    
#       Y_Sa <- array(c(Sa_sigma_lo[1], Sa_sigma_hi[1]))
#       Y_sigma <- array(c(Sa_sigma_lo[2], Sa_sigma_hi[2]))
#       Sa <- interp1(x,Y_Sa,T)   #################
#       sigma <- interp1(x,Y_sigma,T)   #################
#     }
#     else {
      i <- which(period == ip)  
  
      if (Zl == 0) {
        c1 <- c1_it[i]}
      if (Zl == 1) {
        c1 <- c1_it_cas[i]}
      if (Zl == 2) {
        c1 <- c1_it_jp[i]}
      
      log_10_Y <- c1 + c2_it[i]*M + c3_it[i]*h + c4_it[i]*R - 
        g*log10(R) + c5_it[i]*sl*Sc + c6_it[i]*sl*Sd + 
        c7_it[i]*sl*Se   #### % Log10 Sa in cm/s^2
      sigma_10 <- sqrt((sigma1_it[i])^2 + (sigma2_it[i])^2)
      
      Sa <- 10.^(log_10_Y)/981 ## % Median Sa in g
      sigmatotal <- log(10.^sigma_10)
      sigma = log(10.^sigma1_it[i])
      tau = log(10.^sigma2_it[i])
#     }
  }
  if (Zt == 1) {
#     if (length(which(period == T)) == 0) {
#       i_lo <- sum(period<T)   ########
#       T_lo <- period[i_lo]
#       T_hi <- period[i_lo + 1]
#   
#       Sa_sigma_hi <- AB_2003_SZ(T_hi,M,h,Df,Zt,Vs30,Zl)
#       Sa_sigma_lo <- AB_2003_SZ(T_lo,M,h,Df,Zt,Vs30,Zl)
#   
#       x <- array(c(T_lo, T_hi))
#       Y_Sa <- array(c(Sa_sigma_lo[1], Sa_sigma_hi[1]))
#       Y_sigma <- array(c(Sa_sigma_lo[2], Sa_sigma_hi[2]))
#       Sa <- interp1(x,Y_Sa,T)    #########
#       sigma <- interp1(x,Y_sigma,T)
#     }
#     else {
      i <- which(period == ip)   
  
      if (Zl == 0) {
        c1 <- c1_in[i]}
      if (Zl == 1) {
        c1 <- c1_in_cas[i]}
      if (Zl == 2) {
        c1 <- c1_in_jp[i]}
      
      log_10_Y <- c1 + c2_in[i]*M + c3_in[i]*h + c4_in[i]*R - 
        g*log10(R) + c5_in[i]*sl*Sc + c6_in[i]*sl*Sd + 
        c7_in[i]*sl*Se  ### % Log10 Sa in cm/s^2
      sigma_10 <- sqrt((sigma1_in[i])^2 + (sigma2_in[i])^2)
      
      Sa <- 10.^(log_10_Y)/981   #### % Median Sa in g
      sigmatotal <- log(10.^sigma_10)
      sigma = log(10.^sigma1_it[i])
      tau = log(10.^sigma2_it[i])
#     }
  }
  temp <-array(c(Sa, sigmatotal, sigma, tau)) 
  return (temp)
}


AB03SZ.Cal <- function (ip, M, h, Df, Zt, Vs30, Zl)
{
  period <- array(c(0.01, 0.04, 0.1, 0.2, 0.4, 1, 2, 1/0.33))
  if (length(which(period == ip)) == 0) {
    i_lo <- sum(period<ip)   ########
    if (i_lo>=length(period))
    {
      i_lo <- length(period)
      T_lo <- period[i_lo]
      Sa_sigma <- AB03SZ.subCal(T_lo,M,h,Df,Zt,Vs30,Zl)
      
    }
    else
    {
      T_lo <- period[i_lo]
      T_hi <- period[i_lo + 1]      
      Sa_sigma_hi <- AB03SZ.subCal(T_hi,M,h,Df,Zt,Vs30,Zl)
      Sa_sigma_lo <- AB03SZ.subCal(T_lo,M,h,Df,Zt,Vs30,Zl)
      
  #     x <- array(c(T_lo, T_hi))
  #     Y_Sa <- array(c(Sa_sigma_lo[1], Sa_sigma_hi[1]))
  #     Y_sigma <- array(c(Sa_sigma_lo[2], Sa_sigma_hi[2]))
      Sa <- interpolate(ip, T_lo, T_hi, Sa_sigma_lo[1], Sa_sigma_hi[1])  ###(x,Y_Sa,T)    #########
      sigmatotal <- interpolate(ip, T_lo, T_hi, Sa_sigma_lo[2], Sa_sigma_hi[2])
      sigma <- interpolate(ip, T_lo, T_hi, Sa_sigma_lo[3], Sa_sigma_hi[3])
      tau <- interpolate(ip, T_lo, T_hi, Sa_sigma_lo[4], Sa_sigma_hi[4])
      Sa_sigma <- array(c(Sa, sigmatotal, sigma, tau))
    }
  }
  if (length(which(period == ip)) > 0) {
    Sa_sigma <- AB03SZ.subCal(ip,M,h,Df,Zt,Vs30,Zl)
  }
  return (Sa_sigma)
}

#####
AB03SZ.itr <- function (list.mag, list.dist, list.p, list.zh, 
                        flag.source, Vs30=761, Zl=0) {
  m <- length(list.mag)
  n <- length(list.dist)
  le <- length(list.p)
  he = length(list.zh)
  
  if (flag.source=="interface") {
    Zt = 0
  } 
  if (flag.source=="inslab") {
    Zt = 1
  } 
  output.Sa <- array(NA, dim = c(n, m, le, he))
  output.Sd <- array(NA, dim = c(le))
  output.sigma <- array(NA, dim = c(le))
  output.tau <- array(NA, dim = c(le))
  for (hh in 1:he)
  {
    for (j in 1:m)
    {
      for (k in 1:n)
      {
        for (t in 1:le) {
          #         results <- AB03SZ.Cal(T.list[t], ListMag[j], h.list[he], ListDist[k], Zt, Vs30, Zl)
#           dist.temp <- sqrt(list.dist[k]^2 + list.zh[hh]^2)
          dist.temp <- list.dist[k]
          if(list.p[t]>1/.33) {
            results=c(NA,NA,0,0)
          } else {
            results <- AB03SZ.Cal(list.p[t], list.mag[j], list.zh[hh], 
                                  dist.temp, Zt, Vs30, Zl)
          }
          output.Sa[k,j,t,hh] <- results[1] #*980
          output.Sd[t] <- results[2]
          output.sigma[t] = results[3]
          output.tau[t] = results[4]
        }
      }
    }
    # print(h.list[hh])
  }
  output.Sd = cbind(output.Sd, output.sigma, output.tau)
  return(list(output.Sa, output.Sd))
}


# 
# #--------------------------INPUTS-------------------------------#
# #   % T  = Period (sec)
# #   % M  = Moment Magnitude
# #   % h  = Focal (Hypocentral) Depth (km)
# #   % Df = closest distance to the fault surface (km)
# #   % Zt = Subduction Type Indicator: Zt = 0 for interface events
# #   %                                 Zt = 1 for intraslab (in-slab) events
# #   % Vs30 = Shear Wave Velocity averaged over the top 30 meters of soil of the
# #   %      soil profile (m/sec)
# #   % Zl = Cascadia or Japan indicator: Zl = 0 for General Cases
# #   %                                   Zl = 1 for Cascadia
# #   %                                   Zl = 2 for Japan
# # T <- 0.3
# # M <- 7
# h <- 80
# h.list <- as.array(seq(50, 400,10)) #as.array(c(50, 150)) #
# he <- length(h.list)
# # Df <- 50
# Zt <- 1   # Zt = 0 for interface events
# #   %       Zt = 1 for intraslab (in-slab) events
# Vs30 <- 760
# Zl <- 0  #general
# 
# #------------------------ Create the ATN file ------------------#
# # step 1: create the list of magnitudes
# MagMinMax <- array(c(5,8,13))
# ListMag <- as.array(seq(MagMinMax[1],MagMinMax[2], length.out=MagMinMax[3]))
# ListMag <- as.array(7.8)
# m <- dim(ListMag)
# # step 2: create the list of distances
# DistMinMax <- array(c(5, 200,21))
# temp <- seq(log(DistMinMax[1]), log(DistMinMax[2]), length.out=DistMinMax[3])  ####VERIFY THIS
# ListDist <- as.array(exp(temp))
# ListDist <- as.array(c(0,100,500))
# n <- dim(ListDist)
# # step 3: create the list of periods
# # T.list <- as.array(c(1, 3))
# # T.list <- as.array(c(0.01, 0.02,0.04,0.05, 0.075,0.1,0.15, 0.2,
# #                      0.25, 0.3, 0.4, 0.5, 0.6,0.75,1, 1.5,2, 2.5,
# #                      3))
# T.list <- as.array(c(0.01, 0.2,1))
# le <- dim(T.list)
# 
# output.Sa <- array(NA, dim = c(n, m, le, he))
# output.Sd <- array(NA, dim = c(le))
# # step 4: get matrix of ground motions 
# for (hh in 1:he)
# {
#   for (j in 1:m)
#   {
#     for (k in 1:n)
#     {
#       for (t in 1:le) {
# #         results <- AB03SZ.Cal(T.list[t], ListMag[j], h.list[he], ListDist[k], Zt, Vs30, Zl)
#         dist.temp <- sqrt(ListDist[k]^2 + h.list[hh]^2)
#         results <- AB03SZ.Cal(T.list[t], ListMag[j], h.list[hh], dist.temp, Zt, Vs30, Zl)
#         
#         output.Sa[k,j,t,hh] <- results[1]*980
#         output.Sd[t] <- results[2]
#       }
#     }
#   }
# print(h.list[hh])
# }
# 
# # # }
# # # step 4: export to an output file
# # if (Zt == 0)
# #   {listFiles <- array(c("AB03Itf.atn"))}
# # if (Zt == 1)
# # {listFiles <- array(c("AB03Intra.atn"))}
# # for (f in 1:1)
# # {
# #   filename = listFiles[f]
# #   write(MagMinMax, file = filename, ncolumns = 3,
# #         append = FALSE, sep = "   ")
# #   
# #   write(DistMinMax, file = filename, ncolumns = 3,
# #         append = TRUE, sep = "   ")
# #   for (i in 1:le)
# #   {
# #     temp <- array(c(T.list[i], output.Sd[i], 0.00))
# #     write(temp, file = filename, ncolumns = 3,
# #           append = TRUE, sep = "   ")
# #     
# # #     options(scipen=10)  
# #     write(output.Sa[,,i,1], file = filename, ncolumns = DistMinMax[3],
# #           append = TRUE, sep = "   ")
# # #     options(scipen=0)
# #     
# #   }
# # }
# 
# #-----------------------------------------------------------------#
# #         See the change of spectra VS depth                      #
# #-----------------------------------------------------------------#
# # #-------------------INPUTS------------------------#
# # # T <- 0.3
# # M <- 7
# # h <- 80
# # Df <- 50
# # Zt <- 1   # interface
# # Vs30 <- 760
# # Zl <- 0  #general
# # # 
# # 
# # # # step 2: create the list of depths
# # HMinMax <- array(c(40, 100,7))
# # ListH <- as.array(seq(HMinMax[1], HMinMax[2], length.out=HMinMax[3]))  ####VERIFY THIS)
# # n <- length(ListH)
# # # step 3: create the list of periods
# # T.list <- as.array(c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1, 2, 3, 4))
# # le <- dim(T.list)
# # # 
# # output.Sa <- array(NA, dim = c(n, le))
# # output.Sd <- array(NA, dim = c(le))
# # # step 4: get matrix of ground motions 
# # for (j in 1:n)
# # {
# #   for (t in 1:le) {
# #     results <- AB03SZ.Cal(T.list[t], M, ListH[j], Df, Zt, Vs30, Zl)
# #     output.Sa[j,t] <- results[1]
# #     output.Sd[t] <- results[2]
# #   }
# # }
# #                   
# #---------plot-------------
# cols <- cm.colors(le)
# df <- data.frame(x=rep(h.list, le*n), val=as.vector(t((output.Sa[1,1,,]))), 
#                  variable=rep(c("PGA","0.2s", "1s"), each=he),
#                  HDist=rep(c("0km"), each=le*n))
# 
# df <- rbind(df, data.frame(x=rep(h.list, le*n), val=as.vector(t((output.Sa[2,1,,]))), 
#                  variable=rep(c("PGA","0.2s", "1s"), each=he),
#                  HDist=rep(c("100km"), each=le*n)))
# df <- rbind(df, data.frame(x=rep(h.list, le*n), val=as.vector(t((output.Sa[3,1,,]))), 
#                            variable=rep(c("PGA","0.2s", "1s"), each=he),
#                            HDist=rep(c("500km"), each=le*n)))
# # plot
# if (Zt == 0)
# {
#   temp <- sprintf("Interface (M=%1.1f, AB03)", ListMag[1]) #, Df) , D=%dkm
# }
# if (Zt == 1)
# {
# #   temp <- sprintf("Intraslab (M=%1.1f , D=%dkm)", ListMag[1], Df)
#   temp <- sprintf("Intraslab (M=%1.1f, AB03)", ListMag[1]) #, Df) , D=%dkm
# }
# g <- ggplot(data = df, aes(x=x, y=val)) + 
#   geom_line(aes(colour=variable, lty=factor(HDist))) +
#   scale_color_discrete(name="Period") +
#   scale_linetype_manual(name="Horiz. Dist.", values=c(1,2,3)) +
#     xlab("Depth (km)") + ylab("SA (cm/s2)")  + ggtitle(temp)+ xlim(c(50, 400)) +
#   scale_y_log10(limits=c(1,1000),
#                 breaks = c(1,10,100,1000), #c(seq(1,10,1),seq(20,100,10),seq(200,1000,100)), #('log10', function(x) 10^x),
#                 labels = print(c(1,10,100,1000))) 
# #   theme(panel.grid.minor = element_line(colour="blue", size=0.5)) 
# g      
# ggsave("AB03.png")
# #                   
# # 
