

Zhao_2006 <- function(ip,M,x,h,Vs30,FR,SI,SS,MS)
{
# % by 08/27/2013
# % Purpose: Computes the median and dispersion of the Peak Ground
# % Acceleration, PGA, or the Spectral Acceleration, Sa, at 5% Damping.
# 
# % Citation: "Attentuation Relations of Strong Ground Motion in Japan Using
# % Site Classification Based on Predominant Period" by Zhao, John X., Jian
# % Zhang, Akihiro Asano, Yuki Ohno, Taishi Oouchi, Toshimasa Takahashi,
# % Hiroshi Ogawa, Kojiro Irikura, Hong K. Thio, Paul G. Somerville, Yasuhiro
# % Fukushima, and Yoshimitsu Fukushima. Bulletin of the Seismological
# % Society of America, Vol. 96, No. 3, pp. 898-913, June 2006
# 
# % General Limitations: This particular model utilizes data primarily from
# % earthquakes from Japan, with some data from the western United States.
# % The model accounts for shallow crustal events (including reverse fault
#                                                  % type mechanisms), and subduction zone events, both interface and in-slab
# % (intraslab) types. The model yields results for focal depths no greater
# % than 120 km. The focal depth, h, is the greatest influence in this model.
# 
# %-------------------------------INPUTS------------------------------------%
#   
#   % T  = Period (sec)
# % M  = Moment Magnitude
# % x  = Source to Site distance (km); Defines as the shortest distance from
# %      the site to the rupture zone. NOTE: 'x' must be a positive,
# %      non-negative value.
# % h  = Focal (Hypocentral) Depth (km)
# % Vs30  = Average Shear Velocity in the first 30 meters of the soil profile
# %     (m/sec)
# % FR = Reverse-Fault Parameter: FR = 1, For Crustal Earthquakes ONLY IF a
# %                                       Reverse-Fault Exists
# %                               FR = 0, Otherwise
# % SI = Source-Type Indicator: SI = 1, For Interface Events
# %                             SI = 0, Otherwise
# % SS = Source-Type Indicator: SS = 1, For Subduction Slab Events
# %                             SS = 0, Otherwise
# % MS = Magnitude-Squared Term: MS = 1, Includes the Magnitude-squared term
# %                              MS = 0, Does not include magnitude-squared
# %                                      term

# %------------------------------OUTPUTS------------------------------------%
#   
#   % Sa    = Median spectral acceleration prediction (g)
# % sigma = Logarithmic standard deviation of spectral acceleration
# %         prediction
# 
# %-------------------------------Period------------------------------------%
#   
#   period <- array(c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.25, 1.5, 2, 2.5, 3, 4, 5))
# 
# # %-------------------General Equation Coefficients-------------------------%
#   a <- array(c(1.101, 1.076, 1.118, 1.134, 1.147, 1.149, 1.163, 1.2, 1.25, 1.293, 1.336, 1.386, 1.433, 1.479, 1.551, 1.621, 1.694, 1.748, 1.759, 1.826, 1.825))
#   b <- array(c(-0.00564, -0.00671 -0.00787, -0.00722, -0.00659, -0.0059, -0.0052, -0.00422, -0.00338, -0.00282, -0.00258, -0.00242, -0.00232, -0.0022, -0.00207, -0.00224, -0.00201, -0.00187, -0.00147, -0.00195, -0.00237))
#   c <- array(c(0.0055, 0.0075, 0.009, 0.01, 0.012, 0.014, 0.015, 0.01, 0.006, 0.003, 0.0025, 0.0022, 0.002, 0.002, 0.002, 0.002, 0.0025, 0.0028, 0.0032, 0.004, 0.005))
#   d <- array(c(1.08, 1.06, 1.083, 1.053, 1.014, 0.966, 0.934, 0.959, 1.008, 1.088, 1.084, 1.088, 1.109, 1.115, 1.083, 1.091, 1.055, 1.052, 1.025, 1.044, 1.065))
#   e <- array(c(0.01412, 0.01463, 0.01423, 0.01509, 0.01462, 0.01459, 0.01458, 0.01257, 0.01114, 0.01019, 0.00979, 0.00944, 0.00972, 0.01005, 0.01003, 0.00928, 0.00833, 0.00776, 0.00644, 0.0059, 0.0051))
#   Fr <- array(c(0.251, 0.251, 0.24, 0.251, 0.26, 0.269, 0.259, 0.248, 0.247, 0.233, 0.22, 0.232, 0.22, 0.211, 0.251, 0.248, 0.263, 0.262, 0.307, 0.353, 0.248))
#   Si <- array(c(0, 0, 0, 0, 0, 0, 0, -0.041, -0.053, -0.103, -0.146, -0.164, -0.206, -0.239, -0.256, -0.306, -0.321, -0.337, -0.331, -0.39, -0.498))
#   Ss <- array(c(2.607, 2.764, 2.156, 2.161, 1.901, 1.814, 2.181, 2.432, 2.629, 2.702, 2.654, 2.48, 2.332, 2.233, 2.029, 1.589, 0.966, 0.789, 1.037, 0.561, 0.225))
#   Ssl <- array(c(-0.528, -0.551, -0.42, -0.431, -0.372, -0.36, -0.45, -0.506, -0.554, -0.575, -0.572, -0.54, -0.522, -0.509, -0.469, -0.379, -0.248, -0.221, -0.263, -0.169 ,-0.12))
# 
# # %------------Site Class Coefficients & Prediction Uncertainty-------------%
#   CH <- array(c(0.293, 0.939 ,1.499, 1.462, 1.28, 1.121, 0.852, 0.365, -0.207, -0.705, -1.144, -1.609, -2.023, -2.451, -3.243, -3.888, -4.783, -5.444, -5.839, -6.598, -6.752))
#   C1 <- array(c(1.111, 1.684, 2.061, 1.916, 1.669, 1.468, 1.172, 0.655, 0.071, -0.429, -0.866, -1.325, -1.732, -2.152, -2.923, -3.548, -4.41, -5.049, -5.431, -6.181, -6.347))
#   C2 <- array(c(1.344, 1.793, 2.135, 2.168, 2.085, 1.942, 1.683, 1.127, 0.515, -0.003, -0.449, -0.928, -1.349, -1.776, -2.542, -3.169, -4.039, -4.698, -5.089, -5.882, -6.051))
#   C3 <- array(c(1.355, 1.747, 2.031, 2.052, 2.001, 1.941, 1.808, 1.482, 0.934, 0.394, -0.111, -0.62, -1.066, -1.523, -2.327, -2.979, -3.871, -4.496, -4.893, -5.698, -5.873))
#   C4 <- array(c(1.42, 1.814, 2.082, 2.113, 2.03, 1.937, 1.77, 1.397, 0.955, 0.559, 0.188, -0.246, -0.643, -1.084, -1.936, -2.661, -3.64, -4.341, -4.758, -5.588, -5.798))
#   s_t <- array(c(0.604, 0.64, 0.694, 0.702, 0.692, 0.682, 0.67, 0.659, 0.653, 0.653, 0.652, 0.647, 0.653, 0.657, 0.66, 0.664, 0.669, 0.671, 0.667, 0.647, 0.643))
#   t_t <- array(c(0.398, 0.444, 0.49, 0.46, 0.423, 0.391, 0.379, 0.39, 0.389, 0.401, 0.408, 0.418, 0.411, 0.41, 0.402, 0.408, 0.414, 0.411, 0.396, 0.382, 0.377))
#   S_T <- array(c(0.723, 0.779, 0.849, 0.839, 0.811, 0.786, 0.77, 0.766, 0.76, 0.766, 0.769, 0.77, 0.771, 0.775, 0.773, 0.779, 0.787, 0.786, 0.776, 0.751, 0.745))
# 
# # %----------Magnitude Squared Interevent Uncertainty Coefficents-----------%
#   Qc <- array(c(0, 0, 0, 0, 0, 0, 0, 0, -0.0126, -0.0329, -0.0501, -0.065, -0.0781, -0.0899, -0.1148 ,-0.1351, -0.1672, -0.1921, -0.2124, -0.2445, -0.2694))
#   Wc <- array(c(0, 0, 0, 0, 0, 0, 0, 0, 0.0116, 0.0202, 0.0274, 0.0336, 0.0391, 0.044, 0.0545, 0.063, 0.0764, 0.0869, 0.0954, 0.1088, 0.1193))
#   t_c <- array(c(0.303, 0.326, 0.342, 0.331, 0.312, 0.298, 0.3, 0.346, 0.338, 0.349, 0.351, 0.356, 0.348, 0.338, 0.313, 0.306, 0.283, 0.287, 0.278, 0.273, 0.275))
#   QI <- array(c(0, 0, 0, -0.0138, -0.0256, -0.0348, -0.0423, -0.0541, -0.0632, -0.0707, -0.0771, -0.0825, -0.0874, -0.0917, -0.1009, -0.1083, -0.1202, -0.1293, -0.1368, -0.1486, -0.1578))
#   WI <- array(c(0, 0, 0, 0.0286, 0.0352, 0.0403, 0.0445, 0.0511, 0.0562, 0.0604, 0.0639, 0.067, 0.0697, 0.0721, 0.0772, 0.0814, 0.088, 0.0931, 0.0972, 0.1038, 0.109))
#   t_I <- array(c(0.308, 0.343, 0.403, 0.367, 0.328, 0.289, 0.28, 0.271, 0.277, 0.296, 0.313, 0.329, 0.324, 0.328, 0.339, 0.352, 0.36, 0.356, 0.338, 0.307, 0.272))
#   Ps <- array(c(0.1392, 0.1636, 0.169, 0.1669, 0.1631, 0.1588, 0.1544, 0.146, 0.1381, 0.1307, 0.1239, 0.1176, 0.1116, 0.106, 0.0933, 0.0821, 0.0628, 0.0465, 0.0322, 0.0083, -0.0117))
#   Qs <- array(c(0.1584, 0.1932, 0.2057, 0.1984, 0.1856, 0.1714, 0.1573, 0.1309, 0.1078, 0.0878, 0.0705, 0.0556, 0.0426, 0.0314, 0.0093, -0.0062, -0.0235, -0.0287, -0.0261, -0.0065, 0.0246))
#   Ws <- array(c(-0.0529, -0.0841, -0.0877, -0.0773, -0.0644, -0.0515, -0.0395, -0.0183, -0.0008, 0.0136, 0.0254, 0.0352, 0.0432, 0.0498, 0.0612, 0.0674, 0.0692, 0.0622, 0.0496, 0.015, -0.0268))
#   t_s <- array(c(0.321, 0.378, 0.42, 0.372, 0.324, 0.294, 0.284, 0.278, 0.272, 0.285, 0.29, 0.299, 0.289, 0.286, 0.277, 0.282, 0.3, 0.292, 0.274, 0.281, 0.296))
#   
#   
  # read coefficients
  # file <- "C:\\Users\\FELi\\Documents\\others\\crisis2007\\CustomizedAttFunctions\\Zhao06Coeffs.csv"
  file = sprintf("%s\\coeffs\\Zhao06Coeffs.csv", getwd())
  coeffs <- read.csv(file, header = TRUE, sep=",")
  period <- as.array(coeffs$Period)
  
  # % Preliminary Computations & Constants
  if (h >= 125) {  h <- 125}
  
  hc <- 15
  if (h >= hc) {
    delh = 1}
  if (h < hc) {
    delh = 0}
  
  # % Ground Motion Prediction Computation
  if (length(which(period == ip)) == 0) {
  
    i_lo = sum(period<ip)
    if (i_lo >= length(period))
    {
      i_lo = length(period)
      T_lo = period[i_lo]
      Sa_sigma = Zhao_2006(T_lo,M,x,h,Vs30,FR,SI,SS,MS)
    }
    else {
      T_lo = period[i_lo]
      T_hi = period[i_lo + 1]
      
      Sa_sigma_lo = Zhao_2006(T_lo,M,x,h,Vs30,FR,SI,SS,MS)
      Sa_sigma_hi = Zhao_2006(T_hi,M,x,h,Vs30,FR,SI,SS,MS)
      
  #     X = [T_lo T_hi]
  #     Y_Sa = [Sa_lo Sa_hi]
  #     Y_sigma = [sigma_lo sigma_hi]
      Sa = interpolate(ip, T_lo, T_hi, Sa_sigma_lo[1], Sa_sigma_hi[1])
      sigmatotal = interpolate(ip, T_lo, T_hi, Sa_sigma_lo[2], Sa_sigma_hi[2])
      sigma = interpolate(ip, T_lo, T_hi, Sa_sigma_lo[3], Sa_sigma_hi[3])
      tau = interpolate(ip, T_lo, T_hi, Sa_sigma_lo[4], Sa_sigma_hi[4])
      Sa_sigma <- array(c(Sa, sigmatotal, sigma, tau))
    }
  }  else {
    i = which(period == ip)
  
    if (MS == 1) {
      if (SS == 1) {
        Mc = 6.5
        Pst = coeffs$Ps[i]
        Qst = coeffs$Qs[i]
        Wst = coeffs$Ws[i]
        t_t = coeffs$t_s[i]
      } else if (SI == 1) {
        Mc = 6.3
        Pst = 0.0
        Qst = coeffs$QI[i]
        Wst = coeffs$WI[i]
        t_t = coeffs$t_I[i]
      } else {
        Mc = 6.3
        Pst = 0.0
        Qst = coeffs$Qc[i]
        Wst = coeffs$Wc[i]
        t_t = coeffs$t_c[i]
      }
      lnS_MS = Pst*(M - Mc) + Qst*((M - Mc)^2) + Wst
    }  else if (MS == 0) {
      t_t = coeffs$t_t[i]
      lnS_MS = 0
    }
    
    if (Vs30 <= 200) {              Ck <- coeffs$C4[i] }
    if (Vs30 > 200 & Vs30 <= 300) { Ck <- coeffs$C3[i]}
    if (Vs30 > 300 & Vs30 <= 600) { Ck <- coeffs$C2[i] }
    if (Vs30 > 600 & Vs30 <= 1100) { Ck <- coeffs$C1[i] }
    if (Vs30 > 1100) {  Ck <- coeffs$CH[i]}
    
    r <- x + coeffs$c[i]*exp(coeffs$d[i]*M)  
    lnY <- coeffs$a[i]*M + coeffs$b[i]*x - log(r) + 
      coeffs$e[i]*(h-hc)*delh + FR*coeffs$Sr[i] + 
      SI*coeffs$Si[i] + SS*coeffs$Ss[i] + 
      SS*coeffs$SsL[i]*log(x) + Ck + lnS_MS  #### Log Sa in cm/s^2  
    sigmatotal <- sqrt(coeffs$s_t[i]^2 + t_t^2)
    sigma = coeffs$s_t[i]
    tau = t_t
    
    Sa <- exp (lnY)/981    ##### Median Sa in g
    Sa_sigma <- array(c(Sa, sigmatotal, sigma, tau))
  }
  
  return (Sa_sigma)
}


Z06SZ.itr <- function (list.mag, list.rupdist, list.hypodist, list.p, list.zh, 
                       flag.source="interface", Vs30=761, FR=0, MS=0) {
  
  # FR <- 0  # not reverse fault
  SI <- 0  # interface
  SS <- 1  # not intraslab
  
  if(flag.source=="interface"){
    SI <- 1  # interface
    SS <- 0  # not intraslab
    list.dist = list.rupdist
  }
  if(flag.source=="inslab"){
    SI <- 0  # not interface
    SS <- 1  # intraslab
    list.dist = list.hypodist
  }
  
  # MS <- 0  # include the magnitude item
  
  m <- length(list.mag)
  n <- length(list.dist)
  le <- length(list.p)
  he = length(list.zh)
  
  output.Sa <- array(NA, dim = c(n, m, le, he))
  output.Sd <- array(NA, dim = c(le))
  output.sigma <- array(NA, dim = c(le))
  output.tau <- array(NA, dim = c(le))
  hh=1
  j=1
  k=1
  t=2
  
  for (hh in 1:he)
  {
    for (j in 1:m)
    {
      for (k in 1:n)
      {
        for (t in 1:le) {
          if(list.p[t]>5) {
            results=c(NA,NA,0,0)
          } else {
            results <- Zhao_2006(list.p[t],list.mag[j],list.dist[k],
                                 list.zh[hh],Vs30,FR,SI,SS,MS)
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
#   
#   #--------------------------INPUTS-------------------------------#
#   
# #   T <- 0.01
# #   M <- 7
# #   x <- 50
#   h <- 80
#   Vs30 <- 760
#   FR <- 0  # not reverse fault
#   SI <- 0  # interface 
#   SS <- 1  # not intraslab
#   MS <- 0  # include the magnitude item
#   #------------------------ Create the ATN file ------------------#
#   # step 1: create the list of magnitudes
#   MagMinMax <- array(c(5,8,13))
#   ListMag <- as.array(seq(MagMinMax[1],MagMinMax[2], length.out=MagMinMax[3]))
#   m <- dim(ListMag)
#   # step 2: create the list of distances
#   DistMinMax <- array(c(5, 200,21))
#   temp <- seq(log(DistMinMax[1]), log(DistMinMax[2]), length.out=DistMinMax[3])  ####VERIFY THIS
#   ListDist <- as.array(exp(temp))
#   n <- dim(ListDist)
#   # step 3: create the list of periods
#   # T.list <- as.array(c(1, 3))
#   T.list <- as.array(c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1, 2, 3, 4))
#   le <- dim(T.list)
#   
#   output.Sa <- array(NA, dim = c(n, m, le))
#   output.Sd <- array(NA, dim = c(le))
#   # step 4: get matrix of ground motions 
#   # for (i in T.list)
#   # {
#   for (j in 1:m)
#   {
#     for (k in 1:n)
#     {
#       for (t in 1:le) {
#         results <- Zhao_2006(T.list[t],ListMag[j],ListDist[k],h,Vs30,FR,SI,SS,MS)
#         output.Sa[k,j,t] <- results[1]*980
#         output.Sd[t] <- results[2]
#       }
#     }
#   }
#   
#   # }
#   # step 4: export to an output file
# if (SI == 1){
#   listFiles <- array(c("Zhao06_Itf.atn"))}
# if (SS == 1){
#   listFiles <- array(c("Zhao06_Intra.atn"))}
# 
#   for (f in 1:1)
#   {
#     filename = listFiles[f]
#     write(MagMinMax, file = filename, ncolumns = 3,
#           append = FALSE, sep = "   ")
#     
#     write(DistMinMax, file = filename, ncolumns = 3,
#           append = TRUE, sep = "   ")
#     for (i in 1:le)
#     {
#       temp <- array(c(T.list[i], output.Sd[i], 0.00))
#       write(temp, file = filename, ncolumns = 3,
#             append = TRUE, sep = "   ")
#       
#       #     options(scipen=10)  
#       write(output.Sa[,,i], file = filename, ncolumns = DistMinMax[3],
#             append = TRUE, sep = "   ")
#       #     options(scipen=0)
#       
#     }
#   }
#   