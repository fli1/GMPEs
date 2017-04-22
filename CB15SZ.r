# rm(list = ls())

# library(nga)
# library(ggplot2)
# library(scales)

CB15SZ.subCal <- function (ip, M, Zh, R, F_event, Vs30, 
                            epsilon=0) {
#   // GMPE from V. Contreras and R. Boroschek (2015)
  
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

#   % Input:
#     % R  = Closest Distance to the Rupture (Rrup) or Hypocentral Distance
#   % M  = Moment Magnitude
#   % I  = 1:interface, 0:in slab
#   % Zh = hypocentral depth in km 
#   % F_FABA  = 0:forearc or unknown sites, 1: backarc sites 
#   % Vs30    = Average shear-wave velocity in the upper 30 m
#   % epsilon =  Number of standard deviations from median
#   %
#   % With the purpose of extending the utility of the ground motion prediction model 
#   % to other subduction regions, the provided Matlab script maintains the coefficients 
#   % for sites located in backarc zones fitted by the Abrahamson et al. (2016) model (BCHydro).
#   %
#   %  
#   % Montalva, G., Bastias, N., Rodriguez-Marek, A. (2016) Ground Motion
#   % Prediction Equation for the Chilean Subduction Zone. Submitted to BSSA
#   % 
#   %------------------------------OUTPUTS------------------------------------%

#     % Sa    = Median spectral acceleration prediction (g)
#   % sigma = Logarithmic standard deviation of spectral acceleration
#   %         prediction
  
  #-------------------Period----------------------------
#   period <- array(c(0.01, 0.04, 0.1, 0.2, 0.4, 1, 2, 1/0.33))
  # read coefficients

  if(F_event==0){
    file = sprintf("%s\\coeffs\\CB15InterfCoeffs.csv", getwd())  
  } else {
    file = sprintf("%s\\coeffs\\CB15InslabCoeffs.csv", getwd())
  }
  
  coeffs <- read.csv(file, header = TRUE, sep=",")
  period <- as.array(coeffs$Period)
  
  #Preliminary Inital Conditions and Variables Computation
  i <- which(period == ip)   

  c1 = coeffs$c1[i]
  c2 = coeffs$c2[i]
  c3 = coeffs$c3[i]
  c4 = coeffs$c4[i]
  c5 = coeffs$c5[i]
  c6 = coeffs$c6[i]
  c7 = coeffs$c7[i]
  c8 = coeffs$c8[i]
  c9 = coeffs$c9[i]

  sigma= coeffs$sigma[i]

  c1b = coeffs$c1b[i]
  c2b = coeffs$c2b[i]
  c3b = coeffs$c3b[i]
  c4b = coeffs$c4b[i]
  c5b = coeffs$c5b[i]
  c6b = coeffs$c6b[i]
  c7b = coeffs$c7b[i]
  c8b = coeffs$c8b[i]
  c9b = coeffs$c9b[i]
  
  sigmab= coeffs$sigmab[i]
  if (Vs30 >760) {
    z_siteclass = 0  
  }  else {
    z_siteclass = 1
  }
  
  H=Zh
  
  if (M>=6.5) {
    delta = c6 * 10^(c7 * M)
    r_cmb = sqrt( R^2 + delta^2)
    g_term = c8 + c9 * M
    log10Y = c1 + c2*M + c3*H + c4*r_cmb - g_term*log(r_cmb)/log(10) + c5 * z_siteclass
    lnY = log(10^log10Y)
    
    sigmaLnY = log(10^sigma)
  }  else {
    delta = c6b * 10^( c7b * M)
    r_cmb = sqrt(R^2 + delta^2)
    g_term = c8b + c9b*M
    log10Y = c1b + c2b*M + c3b*H + c4b*r_cmb - g_term*log(r_cmb)/log(10) + c5b * z_siteclass
    lnY = log(10^log10Y)
    
    sigmaLnY = log(10^sigmab)
  }
  
  Sa <- exp(lnY)   #### % Median Sa in g
  sigma <- sigmaLnY

  temp <-array(c(Sa, sigma, 0,0)) 
  return (temp)
}


CB15SZ.Cal <- function (ip, M, Zh, R, Fevent, Vs30,  
                          AddMedian = 0)
{
  period <- array(c(0.01,0.04, 0.1,
                    0.2,0.4, 1,
                    2)) 
  
  if (length(which(period == ip)) == 0) {
    i_lo <- sum(period<ip)   ########
    if (i_lo>=length(period))
    {
      i_lo <- length(period)
      T_lo <- period[i_lo]
      Sa_sigma <- CB15SZ.subCal(T_lo,M,Zh, R, Fevent, Vs30, 
                                  AddMedian)
      
    }
    else
    {
      T_lo <- period[i_lo]
      T_hi <- period[i_lo + 1] 
      
      Sa_sigma_hi <- CB15SZ.subCal(T_hi,M,Zh, R, Fevent, Vs30,  
                                    AddMedian)
      Sa_sigma_lo <- CB15SZ.subCal(T_lo,M,Zh, R, Fevent, Vs30, 
                                    AddMedian)
#       print(paste0("Sa_sigma_lo:",Sa_sigma_lo))
#       print(paste0("Sa_sigma_hi:",Sa_sigma_hi))
  #     x <- array(c(T_lo, T_hi))
  #     Y_Sa <- array(c(Sa_sigma_lo[1], Sa_sigma_hi[1]))
  #     Y_sigma <- array(c(Sa_sigma_lo[2], Sa_sigma_hi[2]))
      Sa <- interpolate(ip, T_lo, T_hi, Sa_sigma_lo[1], Sa_sigma_hi[1])  ###(x,Y_Sa,T)    #########
      sigma <- interpolate(ip, T_lo, T_hi, Sa_sigma_lo[2], Sa_sigma_hi[2])
      Sa_sigma <- array(c(Sa, sigma,0,0))
    }
  }
  if (length(which(period == ip)) > 0) {
    Sa_sigma <- CB15SZ.subCal(ip, M, Zh, R, Fevent, Vs30, 
                               AddMedian)
  }
  return (Sa_sigma)
}

#####
CB15SZ.itr <- function (list.mag, list.dist, list.p, list.zh, 
                        flag.source, Vs30=761, 
                        AddMedian=0) {
  m <- length(list.mag)
  n <- length(list.dist)
  le <- length(list.p)
  he = length(list.zh)
  
  if(flag.source=="interface") {
    Fevent=0
  } 
  if(flag.source=="inslab"){
    Fevent=1
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
#           print(list.p[t])
          if(list.p[t]>2) { ## only available up to 2 seconds
            results = c(NA,NA,0,0)
          } else {
            results <- CB15SZ.Cal(list.p[t], list.mag[j], list.zh[hh], 
                                  list.dist[k], Fevent, Vs30,  
                                  AddMedian)
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

