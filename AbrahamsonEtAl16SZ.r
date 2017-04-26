# rm(list = ls())
# library(nga)
# library(ggplot2)
# library(scales)
AEA16SZ.subCal <- function (ip, M, Zh, R, F_event, Vs30, 
                            F_DeltaC1, Delta_C1_input, F_FABA,
                           F_sigma=0, AddMedian = 0) {
#   
  ### the following needs to be updated 
#   % Purpose: Computes the mean and standard deviation of the PGA
#   % or psuedoacceleration, PSA, 5% damping. 
#   
#   % Citation: 
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
#   % Zh  = Focal (Hypocentral) Depth (km)
#   % R = closest distance to the fault surface (km)
#   % Fevent = Subduction Type Indicator: Zt = 0 for interface events
#   %                                 Zt = 1 for intraslab (in-slab) events
#   % Vs30 = Shear Wave Velocity averaged over the top 30 meters of soil of the
#   %      soil profile (m/sec)
#   
#   
#   %------------------------------OUTPUTS------------------------------------%
#     
#     % Sa    = Median spectral acceleration prediction (g)
#   % sigma = Logarithmic standard deviation of spectral acceleration
#   %         prediction
  
  #-------------------Period----------------------------
#   period <- array(c(0.01, 0.04, 0.1, 0.2, 0.4, 1, 2, 1/0.33))
  # read coefficients
  file = sprintf("%s\\coeffs\\AEA16Coeffs.csv", getwd())
  coeffs <- read.csv(file, header = TRUE, sep=",")
  period <- as.array(coeffs$Period)
  
  #Preliminary Inital Conditions and Variables Computation
  C1 = 7.8
  # Zh = Zcoe;
  if (Zh>120) { #//in the 2015 paper
    Zh=120
  }
  
  i <- which(period == ip)   
  
  if (F_event == 0) {    #//interface
    if (F_DeltaC1 == -1)  {
      Delta_C1 = coeffs$Delta_C1_Lower[i]
      Delta_C10 = 0
    } else if (F_DeltaC1 == 0)  {
      Delta_C1 = coeffs$Delta_C1_Central[i]
      Delta_C10 = 0.2
    } else if (F_DeltaC1 == 1)  {
      Delta_C1 = coeffs$Delta_C1_Upper[i]
      Delta_C10 = 0.4
    } else {
      Delta_C1 = coeffs$Delta_C1_Central[i]
      Delta_C10 = 0.2
    }
  } else {
    Delta_C1 = Delta_C1_input
    Delta_C10 = Delta_C1_input
  }
  
  C4 = coeffs$C4[i]
  theta1 = coeffs$theta1[i]
  theta2 = coeffs$theta2[i]
  theta3 = coeffs$theta3[i]
  theta4 = coeffs$theta4[i]
  theta5 = coeffs$theta5[i]
  theta6 = coeffs$theta6[i]
  theta7 = coeffs$theta7[i]
  theta8 = coeffs$theta8[i]
  theta9 = coeffs$theta9[i]
  theta10= coeffs$theta10[i]
  theta11= coeffs$theta11[i]
  theta12= coeffs$theta12[i]
  theta13= coeffs$theta13[i]
  theta14= coeffs$theta14[i]
  theta15= coeffs$theta15[i]
  theta16= coeffs$theta16[i]
  Vlin= coeffs$Vlin[i]

  b= coeffs$b[i]
  c= coeffs$c[i]
  n= coeffs$n[i]
  
  fi= coeffs$fi[i]
  tau= coeffs$tau[i]
  sigma= coeffs$sigma[i]
  sigmaSS= coeffs$sigmaSS[i]
  
  
  if (M <= C1+Delta_C1) {
    f_mag_M = theta4*(M-(C1+Delta_C1)) + theta13*(10-M)^2
    f_mag_M0 = 0.9*(M-(C1+Delta_C10)) -0.0135*(10-M)^2
    
  }  else {
    f_mag_M = theta5*(M-(C1+Delta_C1)) + theta13*(10-M)^2
    f_mag_M0 = 0.0*(M-(C1+Delta_C10)) -0.0135*(10-M)^2
  }
  
  f_depth_Zh = theta11 * (Zh-60)*F_event;
  f_depth_Zh0 = 0.013 * (Zh-60)*F_event;
  
  if (F_event == 1) {
    f_FABA_R = (theta7+theta8*log(max(R,85)/40))*F_FABA
    f_FABA_R0 = (1.0988-1.42*log(max(R,85)/40))*F_FABA
  }   else {
    f_FABA_R = (theta15+theta16*log(max(R,100)/40))*F_FABA
    f_FABA_R0 = (0.9969 - 1 * log(max(R,100)/40))*F_FABA
  }
  
  if (Vs30>1000)  {  
    Vs_star  = 1000
  } else {
    Vs_star = Vs30
  }
  
  lnPGA1000  = 0.98*log(1000/865.1) -1.186*1.18* log(1000/865.1) 
  lnPGA1000  = lnPGA1000 + 4.2203 + 0.9 * Delta_C10 
  lnPGA1000  = lnPGA1000 + (-1.35 -0.4*F_event + 0.1*(M-7.8))*log(R+C4*exp((M-6)*0.4)) 
  lnPGA1000  = lnPGA1000 - 0.0012*R + 3.12*F_event + f_mag_M0 + f_depth_Zh0 + f_FABA_R0
  
  PGA1000  = exp(lnPGA1000)

  
  if (Vs30 < Vlin)  {
    f_site_PGA = theta12 * log(Vs_star/Vlin) - b* log(PGA1000 + c)+ b * log(PGA1000 + c * (Vs_star/Vlin)^n)
    # // f_site_PGA =0.0083;  //
  }  else {
    f_site_PGA = theta12*log(Vs_star/Vlin) + b*n* log(Vs_star/Vlin)
  }
  
  lnY =  theta1 + theta4 * Delta_C1 
  lnY = lnY + (theta2+ theta14*F_event + theta3*(M-7.8))*log(R+C4*exp((M-6)*theta9)) + theta6*R 
  lnY = lnY + theta10*F_event + f_mag_M + f_depth_Zh + f_FABA_R + f_site_PGA
  lnY = lnY + AddMedian
  
  if (F_sigma == 0) {
    sigmaLnY = sigma
  }  else {
    sigmaLnY = sigmaSS
  }
  
  Sa <- exp(lnY)   #### % Median Sa in g
  sigmatotal <- sigmaLnY
  
  temp <-array(c(Sa, sigmatotal, fi, tau)) 
  return (temp)
}


AEA16SZ.Cal <- function (ip, M, Zh, R, Fevent, Vs30,  
                         F_DeltaC1, Delta_C1_input, F_FABA,
                         F_sigma, AddMedian = 0)
{
  period <- array(c(0.01,0.02, 0.05,0.075, 0.1,0.15,
                    0.2,0.25,0.3, 0.4,0.5,0.6,0.75, 1,1.5,
                    2, 2.5,3,4,5,6,7.5,10)) 
  
  if (length(which(period == ip)) == 0) {
    i_lo <- sum(period<ip)   ########
    if (i_lo>=length(period))
    {
      i_lo <- length(period)
      T_lo <- period[i_lo]
      Sa_sigma <- AEA16SZ.subCal(T_lo,M,Zh, R, Fevent, Vs30, 
                                 F_DeltaC1, Delta_C1_input, F_FABA,
                                 F_sigma, AddMedian)
      
    }
    else
    {
      T_lo <- period[i_lo]
      T_hi <- period[i_lo + 1]      
      Sa_sigma_hi <- AEA16SZ.subCal(T_hi,M,Zh, R, Fevent, Vs30,  
                                    F_DeltaC1, Delta_C1_input, F_FABA,F_sigma,AddMedian)
      Sa_sigma_lo <- AEA16SZ.subCal(T_lo,M,Zh, R, Fevent, Vs30, 
                                    F_DeltaC1, Delta_C1_input, F_FABA,F_sigma,AddMedian)
      
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
    Sa_sigma <- AEA16SZ.subCal(ip, M, Zh, R, Fevent, Vs30, 
                               F_DeltaC1, Delta_C1_input, F_FABA,F_sigma,AddMedian)
  }
  return (Sa_sigma)
}

#####
AEA16SZ.itr <- function (list.mag, list.rupdist, list.hypodist, list.p, list.zh, 
                         flag.source, Vs30=760, F_DeltaC1=0, 
                         Delta_C1_input=-0.3, F_FABA=0,
                         F_sigma=0,
                         AddMedian=0) {
  m <- length(list.mag)
  n <- length(list.rupdist)
  le <- length(list.p)
  he = length(list.zh)
  
  if(flag.source=="interface") {
    Fevent=0
    list.dist = list.rupdist
  } else{
    Fevent=1
    list.dist = list.hypodist
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
          
          results <- AEA16SZ.Cal(list.p[t], list.mag[j], list.zh[hh], 
                                 list.dist[k], Fevent, Vs30,  
                                 F_DeltaC1, Delta_C1_input, F_FABA,
                                 F_sigma, AddMedian)
          
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

