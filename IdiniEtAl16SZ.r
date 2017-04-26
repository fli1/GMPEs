

IdiniEA16SZ.subCal <- function (ip, M, Zh, R, F_event, Vs30, 
                                T_star,
                            epsilon=0) {
#   print(paste0("F_event:",F_event))
  
  
#   // GMPE from Idini et al. (2016)
  
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
#   % 
#   % Fevent = Subduction Type Indicator: = 0 for interface events
#   %                                  = 1 for intraslab (in-slab) events
#   % Vs30 = Shear Wave Velocity averaged over the top 30 meters of soil of the
#   %      soil profile (m/sec)
#   
#   % Input:
#     % R  = Closest Distance to the Rupture (Rrup) for Mw>=7.7 interface earthquakes 
#    or Hypocentral Distance for intraslab and Mw<7.7 interface earthquakes
#   % M  = Moment Magnitude
#   % 
#   % Zh = hypocentral depth in km 
#   % epsilon =  Number of standard deviations from median
#   %
#   % With the purpose of extending the utility of the ground motion prediction model 
#   % to other subduction regions, the provided Matlab script maintains the coefficients 
#   % for sites located in backarc zones fitted by the Abrahamson et al. (2016) model (BCHydro).
#   % 
#   %------------------------------OUTPUTS------------------------------------%

#     % Sa    = Median spectral acceleration prediction (g)
#   % sigma = Logarithmic standard deviation of spectral acceleration
#   %         prediction
  
#   
#   ip = 0.01
#   M=8
#   Zh=30
#   R=10
#   F_event=0
#   Vs30=500 
#   T_star = 0.01
#   
  
  #-------------------Period----------------------------
#   period <- array(c(0.01, 0.04, 0.1, 0.2, 0.4, 1, 2, 1/0.33))
  # read coefficients

  file = sprintf("%s\\coeffs\\IdiniEtAl2016Coeffs.csv", getwd())
  
  coeffs <- read.csv(file, header = TRUE, sep=",")
  period <- as.array(coeffs$Period)
  
  #Preliminary Inital Conditions and Variables Computation
  i <- which(period == ip)   

  c1 = coeffs$c1[i]
  c2 = coeffs$c2[i]
  c3 = coeffs$c3[i]
  c4 = 0.1 # #coeffs$c4[i]
  c5 = coeffs$c5[i]
  c6 = 5 #coeffs$c6[i]
  c7 = 0.35 #coeffs$c7[i]
  c8 = coeffs$c8[i]
  c9 = coeffs$c9[i]

  delta_c1 = coeffs$delta_c1[i]
  delta_c2 = coeffs$delta_c2[i]
  delta_c3 = coeffs$delta_c3[i]

  sigma_t= coeffs$sigma_t[i]
  sigma= coeffs$sigma_r[i]
  tau= coeffs$sigma_e[i]

  if(T_star<=0.01) {
    sT = 0
  }else if(T_star<=0.2) {
    sT = coeffs$sII[i]
  } else if(T_star<=0.4) {
    sT = coeffs$sIII[i]
  } else if(T_star<=0.8) {
    sT = coeffs$sIV[i]
  } else  {
    sT = coeffs$sV[i]
  }
  ## the SI and SVI cases were not coded here
  
  h0 = 50
  Mr = 5
  Vref = 1530
# b = -2.1396
  if(F_event == 0) {
    delta_fM = c9*M^2
  }
  if(F_event == 1){
    delta_fM = delta_c1 + delta_c2*M
  }
 
  F_F = c1 + c2*M + c8*(Zh-h0)*F_event + delta_fM

  R0 = (1-F_event) * c6 * 10^(c7*(M-Mr))
  g = c3 + c4*(M-Mr) + delta_c3*F_event
  F_D = g*log(R+R0)/log(10) + c5*R
  
  F_S = sT * log(Vs30/Vref)/log(10)

  log10Y = F_F + F_D + F_S
  lnY = log(10^log10Y)
    
  sigmaLnY = log(10^sigma_t)
  sigma = log(10^sigma)
  tau = log(10^tau)
  
  Sa <- exp(lnY)   #### % Median Sa in g
  sigmatotal <- sigmaLnY

  temp <-array(c(Sa, sigmatotal, sigma, tau)) 
  return (temp)
}


IdiniEtAl16SZ.Cal <- function (ip, M, Zh, R, Fevent, Vs30,T_star, 
                          AddMedian = 0)
{
  period <- array(c(0.01,0.02, 0.03,0.05, 0.07,0.1,
                    0.15,0.2,0.25,0.3,0.4,0.5,0.75, 1,
                    1.5,2,3,4,5,7.5,10)) 
  
  if (length(which(period == ip)) == 0) {
    i_lo <- sum(period<ip)   ########
    if (i_lo>=length(period))
    {
      i_lo <- length(period)
      T_lo <- period[i_lo]
      Sa_sigma <- IdiniEA16SZ.subCal(T_lo,M,Zh, R, Fevent, Vs30, T_star,
                                  AddMedian)
      
    }
    else
    {
      T_lo <- period[i_lo]
      T_hi <- period[i_lo + 1] 
      
      Sa_sigma_hi <- IdiniEA16SZ.subCal(T_hi,M,Zh, R, Fevent, Vs30,  
                                        T_star,AddMedian)
      Sa_sigma_lo <- IdiniEA16SZ.subCal(T_lo,M,Zh, R, Fevent, Vs30, 
                                        T_star,AddMedian)
#       print(paste0("Sa_sigma_lo:",Sa_sigma_lo))
#       print(paste0("Sa_sigma_hi:",Sa_sigma_hi))
  #     x <- array(c(T_lo, T_hi))
  #     Y_Sa <- array(c(Sa_sigma_lo[1], Sa_sigma_hi[1]))
  #     Y_sigma <- array(c(Sa_sigma_lo[2], Sa_sigma_hi[2]))
      Sa <- interpolate(ip, T_lo, T_hi, Sa_sigma_lo[1], Sa_sigma_hi[1])  ###(x,Y_Sa,T)    #########
      sigmatotal <- interpolate(ip, T_lo, T_hi, Sa_sigma_lo[2], Sa_sigma_hi[2])
      sigma <- interpolate(ip, T_lo, T_hi, Sa_sigma_lo[2], Sa_sigma_hi[2])
      tau <- interpolate(ip, T_lo, T_hi, Sa_sigma_lo[2], Sa_sigma_hi[2])
      Sa_sigma <- array(c(Sa, sigmatotal, sigma, tau))
    }
  }
  if (length(which(period == ip)) > 0) {
    Sa_sigma <- IdiniEA16SZ.subCal(ip, M, Zh, R, Fevent, Vs30, 
                                   T_star,AddMedian)
  }
  return (Sa_sigma)
}

#####
IdiniEtAl16SZ.itr <- function (list.mag, list.rupdist, list.hypodist, list.p, list.zh, 
                        flag.source, Vs30=760, T_star,
                        AddMedian=0) {
  m <- length(list.mag)
  n <- length(list.rupdist)
  le <- length(list.p)
  he = length(list.zh)
  
  if(flag.source=="interface") {
    Fevent = 0
    list.dist = list.rupdist
    list.dist[list.mag<7.7] = list.hypodist[list.mag<7.7]
    
  } 
  if(flag.source=="inslab"){
    Fevent = 1
    list.dist = list.hypodist
  }
  
  output.Sa <- array(NA, dim = c(n, m, le, he))
  output.Sd <- array(NA, dim = c(le))
  output.sigma <- array(NA, dim = c(le))
  output.tau <- array(NA, dim = c(le))
  icount=0
  for (hh in 1:he)
  {
    for (j in 1:m)
    {
      for (k in 1:n)
      {
        for (t in 1:le) {
#           print(list.p[t])
          if(list.p[t]>10) { ## only available up to 10 seconds
            results = c(NA,NA,0,0)
          } else {
            results <- IdiniEtAl16SZ.Cal(list.p[t], list.mag[j], list.zh[hh], 
                                  list.dist[k], Fevent, Vs30,  T_star,
                                  AddMedian)
            icount = icount+1
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
# print(icount)
# print(paste0("depth, mag,dist,period",c(he, m,n,le)))
  output.Sd = cbind(output.Sd, output.sigma, output.tau)
  return(list(output.Sa, output.Sd))
}

