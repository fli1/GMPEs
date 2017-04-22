
# library(nga)

Youngs_1997_SZ <- function(M, ip, r, H, Zt, Zr,AddMedian=0)
{

  # % Purpose: Computes the median and logaritmic standard deviation of a
  # %          subduction zone earthquake with 5% damping.
  # 
  # % Citation: "Strong Ground Motion Attentuation Relationships for Subduction
  # % Zone Earthquakes" by Youngs, R.R., S.J. Chiou, W.J. Silva, J.R. Humphrey.
  # % Seismological Research Letters, vol. 68, no. 1, p. 58-73
  # 
  # % General Limitations: According to the authors, "The attenuation
  # % relationships developed in this study are considered appropriate for
  # % earthquakes of magnitude M 5 and greater and for distances to the rupture
  # % surface of 10 to 500km."
  # 
  # %----------------------------------INPUTS---------------------------------%
  #   
  #   % M   = Moment Magnitude
  # % T   = Period of Interest (sec.). For PGA Computation, T = 0
  # %       For Rock Sites, 0<=T<=3
  # %       For Soil Sites, 0<=T<=4
  # % r   = Source to site distance to rupture surface (km)
  # % H   = Hypocentral (Focal) depth from surface to focus (km)
  # % Zt  = Subduction Zone Type: For Intraslab, Zt = 1
  # %                             For Interface, Zt = 0
  # % Zr  = Rock Site Indicator: For Rock Sites, Zr = 1
  # %                                 Otherwise, Zr = 0
  # 

  # %---------------------------------OUTPUTS---------------------------------%
  #   
  #   % Sa    = Median spectral acceleration prediction (g)
  # % sigma = Logarithmic standard deviation of spectral acceleration
  # %         prediction
  # 
  # %--------------------------Generic Rock Period----------------------------%
    period_rock <- array(c(0, 0.075, 0.10, 0.20, 0.30, 0.40, 0.50, 0.75, 1.00, 1.50, 2.00, 3.00))
  
  # %-------------------------Generic Rock Coefficients-----------------------%
    grc1  <- array(c(0.2418, 1.5168, 1.4298 ,0.9638, 0.4878, 0.1268, -0.1582, -0.9072, -1.4942, -2.3922, -3.0862, -4.2692))
  grc2 <- array(c(1.414, 1.414, 1.414, 1.414, 1.414, 1.414, 1.414, 1.414, 1.414, 1.414, 1.414, 1.414))
  grc3 <- array(c(0.0000, 0.0000, -0.0011, -0.0027, -0.0036, -0.0043, -0.0048, -0.0057, -0.0064, -0.0073, -0.0080, -0.0089))
  grc4 <- array(c(-2.552, -2.707, -2.655, -2.528, -2.454, -2.401, -2.360, -2.286, -2.234, -2.160 ,-2.107, -2.033))
  grc5 <- array(c(0.00607, 0.00607, 0.00607, 0.00607, 0.00607, 0.00607, 0.00607, 0.00607, 0.00607, 0.00607, 0.00607, 0.00607))
  grc6 <- array(c(0.3846, 0.3846, 0.3846, 0.3846, 0.3846, 0.3846, 0.3846, 0.3846, 0.3846, 0.3846, 0.3846, 0.3846))
  grc7 <- array(c(1.7818, 1.7818, 1.7818, 1.7818, 1.7818, 1.7818, 1.7818, 1.7818, 1.7818, 1.7818, 1.7818, 1.7818))
  grc8 <- array(c(0.554, 0.554, 0.554, 0.554, 0.554, 0.554, 0.554, 0.554, 0.554, 0.554, 0.554, 0.554))
  grc9 <- array(c(1.45, 1.45, 1.45, 1.45, 1.45, 1.45, 1.45, 1.45, 1.45, 1.50, 1.55, 1.65))
  grc10 <- array(c(-0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1))
  grc11 <- array(c(0.650, 0.650, 0.650, 0.650, 0.650, 0.650, 0.650, 0.650, 0.650, 0.700, 0.750, 0.850))
    
  
  # %--------------------------Generic Soil Period----------------------------%
    period_soil <- array(c(0.0, 0.075, 0.10, 0.20, 0.30, 0.40, 0.50, 0.75, 1.00, 1.50, 2.00, 3.00, 4.00))
  
  # %-------------------------Generic Soil Coefficients-----------------------%
  gsc1 <- array(c(-0.6687, 1.7313, 1.8473, 0.8803, 0.1243, -0.5247, -1.1067, -2.3727, -3.5387, -5.7697, -7.1017, -7.3407, -8.2867))
  gsc2 <- array(c(1.438, 1.438, 1.438, 1.438, 1.438, 1.438, 1.438, 1.438, 1.438, 1.438, 1.438, 1.438, 1.438))
  gsc3 <- array(c(0.0000, -0.0019, -0.0019, -0.0019, -0.0020, -0.0020, -0.0035, -0.0048, -0.0066, -0.0114, -0.0164, -0.0221, -0.0235))
  gsc4 <- array(c(-2.329, -2.697, -2.697, -2.464, -2.327, -2.230, -2.140, -1.952, -1.785, -1.470, -1.290, -1.347, -1.272))
  gsc5 <- array(c(0.00648, 0.00648, 0.00648, 0.00648, 0.00648, 0.00648, 0.00648, 0.00648, 0.00648, 0.00648, 0.00648, 0.00648, 0.00648))
  gsc6 <- array(c(0.3648, 0.3648, 0.3648, 0.3648, 0.3648, 0.3648, 0.3648, 0.3648, 0.3648, 0.3648, 0.3648, 0.3648, 0.3648))
  gsc7 <- array(c(1.0970, 1.0970, 1.0970, 1.0970, 1.0970, 1.0970, 1.0970, 1.0970, 1.0970, 1.0970, 1.0970, 1.0970, 1.0970))
  gsc8 <- array(c(0.617, 0.617, 0.617, 0.617, 0.617, 0.617, 0.617, 0.617, 0.617, 0.617, 0.617, 0.617, 0.617))
  gsc9 <- array(c(1.45, 1.45, 1.45, 1.45, 1.45, 1.45, 1.45, 1.45, 1.45, 1.50, 1.55, 1.65, 1.65))
  gsc10 <- array(c(-0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1))
  gsc11 <- array(c(0.650, 0.650, 0.650, 0.650, 0.650, 0.650, 0.650, 0.650, 0.650, 0.700, 0.750, 0.850, 0.850))
  
  
  # % Computation of Parameters
  # % For Sa Computation in Generic Soil
  if (Zr == 0) {
    if (length(which(period_soil == ip)) == 0) {
      i_low <- sum(period_soil < ip)
      if (i_low >= length(period_soil)){
        i_low <- length(period_soil)
        T_low <- period_soil[i_low]
        Sa_sigma <- Youngs_1997_SZ(M, T_low, r, H, Zt, Zr)
      }
      else {
        T_low <- period_soil[i_low]
        T_high <- period_soil[i_low + 1]
        Sa_sigma_high <- Youngs_1997_SZ(M, T_high, r, H, Zt, Zr)
        Sa_sigma_low <- Youngs_1997_SZ(M, T_low, r, H, Zt, Zr)
    
    #     x = [T_low T_high]
    #     Y_Sa = [Sa_low Sa_high]
    #     Y_sigma = [sigma_low sigma_high]
        Sa <- interpolate(ip, T_low, T_high, Sa_sigma_low[1], Sa_sigma_high[1])
        sigma <- interpolate(ip, T_low, T_high, Sa_sigma_low[2], Sa_sigma_high[2])
        Sa_sigma <- array(c(Sa, sigma, 0,0))
      }
    }
    else {
      i <- which(period_soil == ip)
      
      ln_Y <- gsc1[i] + gsc2[i]*M + gsc3[i]*(10-M)^3 + 
        gsc4[i]*(log((r)+gsc7[i]*exp(gsc8[i]*M))) + gsc5[i]*H + gsc6[i]*Zt
    
      sigma <- max((gsc9[i] + gsc10[i]*M), gsc11[i])
      Sa <- exp(ln_Y)
      Sa_sigma <- array(c(Sa, sigma, 0,0))
    }
  }
  
  # % For Sa Computation in Generic Rock
  if (Zr == 1) {
    if (length(which(period_rock == ip)) == 0) {
    
      i_low <- sum(period_rock < ip)
      if (i_low >= length(period_rock)){
        i_low <- length(period_rock)
        T_low <- period_rock[i_low]
        Sa_sigma <- Youngs_1997_SZ(M, T_low, r, H, Zt, Zr)
      }
      else {
        T_low <- period_rock[i_low]
        T_high <- period_rock[i_low + 1]
        Sa_sigma_high <- Youngs_1997_SZ(M, T_high, r, H, Zt, Zr)
        Sa_sigma_low <- Youngs_1997_SZ(M, T_low, r, H, Zt, Zr)
        
        #     x = [T_low T_high]
        #     Y_Sa = [Sa_low Sa_high]
        #     Y_sigma = [sigma_low sigma_high]
        Sa <- interpolate(ip, T_low, T_high, Sa_sigma_low[1], Sa_sigma_high[1])
        sigma <- interpolate(ip, T_low, T_high, Sa_sigma_low[2], Sa_sigma_high[2])
        Sa_sigma <- array(c(Sa, sigma, 0,0))
      }
    }
    else {
      i <- which(period_rock == ip)
      ln_Y <- grc1[i] + grc2[i]*M + grc3[i]*(10-M)^3 + 
        grc4[i]*(log((r)+grc7[i]*exp(grc8[i]*M))) + grc5[i]*H + grc6[i]*Zt
      sigma <- max((grc9[i] + grc10[i]*M), grc11[i])
  
      Sa <- exp(ln_Y)
      Sa_sigma <- array(c(Sa, sigma, 0,0))
    }
  }
  return (Sa_sigma)
}


#####
YEA97SZ.itr <- function (list.mag, list.dist, list.p, list.zh, 
                               flag.source, Vs30=760, 
                               AddMedian=0) {
  m <- length(list.mag)
  n <- length(list.dist)
  le <- length(list.p)
  he = length(list.zh)
  
  if(flag.source=="interface") {
    Zt = 0
  } 
  if(flag.source=="inslab"){
    Zt = 1
  }
  
  ### this is a simple rule set by Feng Li to differentiate rock or soil site 
  ### According to Youngs et al.(1997), rock is close to site class A defined by Boore et al. (1993), which is 
  ### close to NEHRP site Class B; soil is close to site class C in Boore et al. (1993), which is close to 
  ### NEHRP site class D.
  if(Vs30 > 760) {
    Zr = 1
  } else {
    Zr = 0
  }  
  
  output.Sa <- array(NA, dim = c(n, m, le, he))
  output.Sd <- array(NA, dim = c(le))
  output.sigma <- array(NA, dim = c(le))
  output.tau <- array(NA, dim = c(le))
#   icount=0
  for (hh in 1:he)
  {
    for (j in 1:m)
    {
      for (k in 1:n)
      {
        for (t in 1:le) {
          #           print(list.p[t])
          if(list.p[t]>3) { ## only available up to 3 seconds
            results = c(NA,NA,0,0)
          } else {
            results <- Youngs_1997_SZ(list.mag[j], list.p[t], list.dist[k],
                                      list.zh[hh], 
                                      Zt, Zr, AddMedian)

#             icount = icount+1
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
# # T <- 0.3
# # M <- 7
# # % H   = Hypocentral (Focal) depth from surface to focus (km)
# # % Zt  = Subduction Zone Type: For Intraslab, Zt = 1
# # %                             For Interface, Zt = 0
# # % Zr  = Rock Site Indicator: For Rock Sites, Zr = 1
# # %                                 Otherwise, Zr = 0
# # 
# r <- 50
# H <- 30
# Zt <- 0  #interface
# Zr <- 1   #rock site
# #------------------------ Create the ATN file ------------------#
# # step 1: create the list of magnitudes
# MagMinMax <- array(c(5,9,17))
# ListMag <- as.array(seq(MagMinMax[1],MagMinMax[2], length.out=MagMinMax[3]))
# # ListMag <-as.array(c(8.5,9))
# m <- dim(ListMag)
# # step 2: create the list of distances
# DistMinMax <- array(c(5, 500,50))
# temp <- seq(log(DistMinMax[1]), log(DistMinMax[2]), length.out=DistMinMax[3])  ####VERIFY THIS
# ListDist <- as.array(exp(temp))
# # ListDist<-as.array(c(52, 100))
# n <- dim(ListDist)
# # step 3: create the list of periods
# # T.list <- as.array(c(1, 3))
# T.list <- as.array(c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1, 2, 3, 4))
# le <- dim(T.list)
# 
# output.Sa <- array(NA, dim = c(n, m, le))
# output.Sd <- array(NA, dim = c(le))
# # step 4: get matrix of ground motions 
# # for (i in T.list)
# # {
# for (j in 1:m)
# {
#   for (k in 1:n)
#   {
#     for (t in 1:le) {
#       results <- Youngs_1997_SZ(ListMag[j], T.list[t], ListDist[k], H, Zt, Zr)
#       output.Sa[k,j,t] <- exp(log(results[1])+0*results[2]) *981# 84th prct
#       output.Sd[t] <- results[2]
#     }
#   }
# }
# 
# # }
# # step 4: export to an output file
# listFiles <- array(c("Youngs97SZ_itf.atn"))
# for (f in 1:1)
# {
#   filename = listFiles[f]
#   write(MagMinMax, file = filename, ncolumns = 3,
#         append = FALSE, sep = "   ")
#   
#   write(DistMinMax, file = filename, ncolumns = 3,
#         append = TRUE, sep = "   ")
#   for (i in 1:le)
#   {
#     temp <- array(c(T.list[i], output.Sd[i], 0.00))
#     write(temp, file = filename, ncolumns = 3,
#           append = TRUE, sep = "   ")
#     
#     #     options(scipen=10)  
#     write(output.Sa[,,i], file = filename, ncolumns = DistMinMax[3],
#           append = TRUE, sep = "   ")
#     #     options(scipen=0)
#     
#   }
# }
# 
