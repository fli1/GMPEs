

Zhao16SZ_Interf.subCal <- function (ip, M, R, Zh, Vs30,  
                                  Xv=0)
{
  ##
  ## R :
  ## Zh: 
  #-------------------Period----------------------------
  file = sprintf("%s\\coeffs\\Zhao16InterfCoeffs.csv", getwd())
  coeffs <- read.csv(file, header = TRUE, sep=",")
  period <- as.array(coeffs$Period)
  # period <- array(c(0.01, 0.04, 0.1, 0.2, 0.4, 1, 2, 1/0.33))
  #   %---------------------Interslab Events Coefficients-----------------------%
  
  #Preliminary Inital Conditions and Variables Computation
  mc = 7.1
  hc = 25
  i <- which(period == ip)
  
  SC1         = coeffs$ SC1        [i]
  SC2         = coeffs$ SC2        [i]
  SC3         = coeffs$ SC3        [i]
  SC4         = coeffs$ SC4        [i]
  c1          = coeffs$ c1         [i]
  c2          = coeffs$ c2         [i]
  cintD       = coeffs$ cintD      [i]
  cintS       = coeffs$ cintS      [i]
  dint        = coeffs$ dint       [i]
  gammaintS   = coeffs$ gammaintS  [i]
  bint        = coeffs$ bint       [i]
  gint        = coeffs$ gint       [i]
  gintLD      = coeffs$ gintLD     [i]
  gintLS      = coeffs$ gintLS     [i]
  eintV       = coeffs$ eintV      [i]
  eintS       = coeffs$ eintS      [i]
  gammaint    = coeffs$ gammaint   [i]
  RSF         = coeffs$ RSF        [i]
  S2          = coeffs$ S2         [i]
  S3          = coeffs$ S3         [i]
  S4          = coeffs$ S4         [i]
  S5          = coeffs$ S5         [i]
  S6          = coeffs$ S6         [i]
  S7          = coeffs$ S7         [i]
  sigma0       = coeffs$ sigma      [i]
  tau0         = coeffs$ tau        [i]
  sT          = coeffs$ sT         [i]
  SC1sigmaS   = coeffs$ SC1sigmaS  [i]
  SC1tauS     = coeffs$ SC1tauS    [i]
  SC1ST       = coeffs$ SC1ST      [i]
  SC2sigmaS   = coeffs$ SC2sigmaS  [i]
  SC2tauS     = coeffs$ SC2tauS    [i]
  SC2ST       = coeffs$ SC2ST      [i]
  SC3sigmaS   = coeffs$ SC3sigmaS  [i]
  SC3tauS     = coeffs$ SC3tauS    [i]
  SC3ST       = coeffs$ SC3ST      [i]
  SC4sigmaS   = coeffs$ SC4sigmaS  [i]
  SC4tauS     = coeffs$ SC4tauS    [i]
  SC4ST       = coeffs$ SC4ST      [i]
  Amax_s1     = coeffs$ Amax_s1    [i]
  Src_s1      = coeffs$ Src_s1     [i]
  Ca_s1       = coeffs$ Ca_s1      [i]
  Cb_s1       = coeffs$ Cb_s1      [i]
  Ca_s2       = coeffs$ Ca_s2      [i]
  Cb_s2       = coeffs$ Cb_s2      [i]
  Ca_s3       = coeffs$ Ca_s3      [i]
  Cb_s3       = coeffs$ Cb_s3      [i]
  Amax_s2     = coeffs$ Amax_s2    [i]
  Src_s2      = coeffs$ Src_s2     [i]
  Amax_s3     = coeffs$ Amax_s3    [i]
  Src_s3      = coeffs$ Src_s3     [i]
  Amax_s4     = coeffs$ Amax_s4    [i]
  Src_s4   =coeffs$Src_s4  [i]
  Ca_s4    =coeffs$Ca_s4   [i]
  Cb_s4    =coeffs$Cb_s4   [i]
  
  
  if (Vs30> 1100) {
    siteclass = 0
    Imf = 0
  }  else if (Vs30 > 600) {
    siteclass = 1
    Imf = 0.91
  }  else if (Vs30 > 300) {
    siteclass = 2
    Imf = 1.023
  }  else if (Vs30 > 200) {
    siteclass = 3
    Imf = 1.034
  }  else{
    siteclass = 4
    Imf = 0.737
  }
  
  if (Zh > hc) {
    deltah = 1  
  }   else {
    deltah = 0
  }
  
  if (M > mc) {
    Cm = mc  
  }  else {
    Cm = M
  }
  
  # //calculate the hard rock motion first
  rij = 10 + R + exp(c1 + 1.151 * Cm)
  # print(deltah)
  if (deltah == 0) { ## shallow
    if (M > mc) {
      fmintS = bint*Zh + gammaintS + cintS*mc + dint*(M-mc)
    }    else {
      fmintS = bint*Zh + gammaintS + cintS*M
    }
    lnSr = fmintS + gint * log(rij) + 
      gintLS*log(R+200) + eintS*R + eintV*Xv + gammaint - log(RSF)
  }  else { ## deep
    if (M > mc) {
      fmintD = bint*Zh + cintD*mc + dint*(M-mc)  
    } else {
      fmintD = bint*Zh + cintD*M
    }
    lnSr = fmintD + gint*log(rij) + gintLD*log(R+200) + 
      eintV*Xv + gammaint - log(RSF) #+ log(60)
    # print(lnSr)
  } 
  
  Sr = exp(lnSr)
  Sreff = Sr*Imf
  
  # //calculate the site term
  alpha = 2
  beta = 0.6
  theta=10
  if (siteclass==1) {
    Anmax = RSF
    SreffC = Imf*Src_s1
    Sf = Anmax/Amax_s1
    if (Anmax >= 1.25) {
      Snc = (exp((log(Anmax)*log(SreffC^alpha + beta) - log(Sf)*log(beta))/log(Amax_s1))-beta)^(1/alpha)
      
    } else{
      if (Ca_s1 != 0) {
        temp_up = Ca_s1 * (alpha - 1) * log(beta) * log(theta * beta) - log(theta) * (Cb_s1 + log(Sf))
        temp_under = Ca_s1 * (alpha * log(theta * beta) - log(beta))
        Snc = exp(temp_up/temp_under)
        
      } else {
        Snc = 0
      }
    }      
    Smr = Sreff * (Snc / SreffC) *SC1
    lnA = -log(Amax_s1) * ((log(Smr ^ alpha + beta) - log(beta)) / (log(SreffC ^ alpha + beta) - log(beta)))
  } else if (siteclass == 2) {
    if (Zh <= hc) {
      Anmax=RSF*exp(S2)
      ln_linearterm = S2
    } else {
      Anmax = RSF*exp(S5)
      ln_linearterm = S5
    }
    SreffC = Imf*Src_s2
    Sf = Anmax/Amax_s2
    if (Anmax >= 1.25) {
      Snc = (exp((log(Anmax)*log(SreffC^alpha + beta) - log(Sf)*log(beta))/log(Amax_s2))-beta)^(1/alpha)
    } else{
      if (Ca_s2 != 0) {
        Snc = exp((Ca_s2*(alpha-1)*log(beta)*log(theta*beta)-log(theta)*(Cb_s2+log(Sf)))/(Ca_s2*(alpha*log(theta*beta)-log(beta))))
      } else{
        Snc = 0
      }
    }
    Smr = Sreff*(Snc/SreffC)*SC2
    lnA = ln_linearterm-log(Amax_s2)*((log(Smr^alpha + beta) - log(beta))/(log(SreffC^alpha + beta)-log(beta)))
  }  else if (siteclass==3) {
    if (Zh <= hc) {
      Anmax = RSF*exp(S3)
      ln_linearterm = S3
    } else {
      Anmax = RSF*exp(S6)
      ln_linearterm = S6
    }
    SreffC = Imf*Src_s3
    Sf = Anmax/Amax_s3
    if (Anmax >= 1.25) {
      Snc = (exp((log(Anmax)*log(SreffC^alpha+beta)-log(Sf)*log(beta))/log(Amax_s3))-beta)^(1/alpha)
    } else{
      if (Ca_s3 != 0) {
        Snc = exp((Ca_s3*(alpha-1)*log(beta)*log(theta*beta)-log(theta)*(Cb_s3+log(Sf)))/(Ca_s3*(alpha*log(beta*theta)-log(beta))))
      }  else{
      Snc = 0
      } 
    } 
    Smr = Sreff*(Snc/SreffC)*SC3
    lnA = ln_linearterm-log(Amax_s3)*((log(Smr^alpha+beta)-log(beta))/(log(SreffC^alpha+beta)-log(beta)))
    
  }  else if (siteclass == 4) {
    if (Zh <= hc) {
      Anmax = RSF*exp(S4)
    ln_linearterm = S4
    } else{
      Anmax = RSF*exp(S7)
      ln_linearterm = S7
    }
    SreffC = Imf*Src_s4
    Sf = Anmax/Amax_s4
    if (Anmax >= 1.25) {
      Snc = (exp((log(Anmax)*log(SreffC^alpha+beta)-log(Sf)*log(beta))/log(Amax_s4))-beta)^(1/alpha)
    } else {
      if (Ca_s4 != 0) {
        Snc = exp((Ca_s4*(alpha-1)*log(beta)*log(theta*beta)-log(theta)*(Cb_s4+log(Sf)))/(Ca_s4*(alpha*log(theta*beta)-log(beta))))
      }   else {
        Snc = 0
      }
    }
    Smr = Sreff*(Snc/SreffC)*SC4
    lnA = ln_linearterm - log(Amax_s4)*((log(Smr^alpha+beta)-log(beta))/(log(SreffC^alpha+beta)-log(beta)))
  }  else if (siteclass == 0) {
    lnA = -log(RSF)
  }
  
  # //calculate ground motions
  if (deltah == 0) {
    lnY = fmintS + gint*log(rij) + gintLS*log(R+200) + 
      eintS*R + eintV*Xv + gammaint + lnA
  }  else {
    lnY = fmintD + gint*log(rij) + gintLD*log(R+200) + 
      eintV*Xv  + gammaint + lnA #+log(60)+log(0.85) ##!!!! temporary solution 
  }
  
  # //standard deviation
  # if (siteclass  == 1) {
  #   sigmaLnY = SC1ST
  #   sigma = SC1sigmaS
  #   tau = SC1tauS
  # }  else if (siteclass  == 2) {
  #   sigmaLnY = SC2ST
  #   sigma = SC2sigmaS
  #   tau = SC2tauS
  # } else if (siteclass  == 3) {
  #   sigmaLnY = SC3ST
  #   sigma = SC3sigmaS
  #   tau = SC3tauS
  # } else if (siteclass  == 4) {
  #   sigmaLnY = SC4ST
  #   sigma = SC4sigmaS
  #   tau = SC4tauS
  # } else {
    sigmaLnY = sT
    sigma = sigma0
    tau = tau0
  # }
  Sa = exp(lnY)
  sigmatotal = sigmaLnY
  
  temp <-array(c(Sa, sigmatotal, sigma, tau)) 
  return (temp)
}

####
Zhao16SZ_Interf.Cal <- function (ip, M, R, Zh, Vs30, Xv,AddMedian)
{
  period <- array(c(0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08, 0.09,0.1,
                    0.12,0.14,0.15,0.16,0.18,0.2,0.25,0.3,0.35,0.4,0.45,
                    0.5, 0.6,0.7,0.8,0.9, 1,1.25,1.5,2,2.5,3, 3.5,4,4.5,5))
  if (length(which(period == ip)) == 0) {
    i_lo <- sum(period<ip)   ########
    if (i_lo>=length(period))
    {
      i_lo <- length(period)
      T_lo <- period[i_lo]
      Sa_sigma <- Zhao16SZ_Interf.subCal(T_lo,M,R, Zh,Vs30,Xv)
    }
    else
    {
      T_lo <- period[i_lo]
      T_hi <- period[i_lo + 1]      
      Sa_sigma_hi <- Zhao16SZ_Interf.subCal(T_hi,M,R, Zh,Vs30,Xv)
      Sa_sigma_lo <- Zhao16SZ_Interf.subCal(T_lo,M,R, Zh,Vs30,Xv)
      #     x <- array(c(T_lo, T_hi))
      #     Y_Sa <- array(c(Sa_sigma_lo[1], Sa_sigma_hi[1]))
      #     Y_sigma <- array(c(Sa_sigma_lo[2], Sa_sigma_hi[2]))
      Sa <- interpolate(ip, T_lo, T_hi, Sa_sigma_lo[1], Sa_sigma_hi[1])  ##(x,Y_Sa,T)    #########
      sigmatotal <- interpolate(ip, T_lo, T_hi, Sa_sigma_lo[2], Sa_sigma_hi[2])
      sigma = interpolate(ip, T_lo, T_hi, Sa_sigma_lo[3], Sa_sigma_hi[3])
      tau = interpolate(ip, T_lo, T_hi, Sa_sigma_lo[4], Sa_sigma_hi[4])
      Sa_sigma <- array(c(Sa, sigmatotal, sigma, tau))
    }
  }
  if (length(which(period == ip)) > 0) {
    Sa_sigma <- Zhao16SZ_Interf.subCal(ip,M,R, Zh,Vs30,Xv)
  }
  Sa_sigma[1] = Sa_sigma[1]+AddMedian
  return (Sa_sigma)
}



Zhao16SZ_Inslab.subCal <- function (ip, M, R, Zh, Vs30,  
                                  Xv=0)
{
  ##
  ## R :
  ## Zh: 
  #-------------------Period----------------------------
  file = sprintf("%s\\coeffs\\Zhao16InslabCoeffs.csv", getwd())
  coeffs <- read.csv(file, header = TRUE, sep=",")
  period <- as.array(coeffs$Period)
  # period <- array(c(0.01, 0.04, 0.1, 0.2, 0.4, 1, 2, 1/0.33))
  #   %---------------------Interslab Events Coefficients-----------------------%
  
  #Preliminary Inital Conditions and Variables Computation
  mc = 7.1
  msc = 6.3
  i <- which(period == ip)
  
  AmSC1 = coeffs$ AmSC1        [i]
  SC1         = coeffs$ SC1        [i]
  SC2         = coeffs$ SC2        [i]
  SC3         = coeffs$ SC3        [i]
  SC4         = coeffs$ SC4        [i]
  c1          = coeffs$ c1         [i]
  # c2          = coeffs$ c2         [i]
  cSL1  = coeffs$ cSL1      [i]
  cSL2 = coeffs$ cSL2      [i]
  dSL      =coeffs$dSL     [i]
  bSL      =coeffs$bSL     [i]
  gSL      =coeffs$gSL     [i]
  gSLL     =coeffs$gSLL    [i]
  eSLv     =coeffs$eSLv    [i]
  eSL      =coeffs$eSL     [i]
  eSLH     =coeffs$eSLH    [i]
  gamma    =coeffs$gamma   [i]
  S2       =coeffs$S2      [i]
  S3       =coeffs$S3      [i]
  S4       =coeffs$S4      [i]
  sigma0    =coeffs$sigma   [i]
  tau0      =coeffs$tau     [i]
  sT       =coeffs$sT      [i]
  Amax_s1  =coeffs$Amax_s1 [i]
  Src_s1   =coeffs$Src_s1  [i]
  Imf_s1   =coeffs$Imf_s1  [i]
  Ca_s1    =coeffs$Ca_s1   [i]
  Cb_s1    =coeffs$Cb_s1   [i]
  Amax_s2  =coeffs$Amax_s2 [i]
  Src_s2   =coeffs$Src_s2  [i]
  Ca_s2    =coeffs$Ca_s2   [i]
  Cb_s2    =coeffs$Cb_s2   [i]
  Amax_s3  =coeffs$Amax_s3 [i]
  Src_s3   =coeffs$Src_s3  [i]
  Ca_s3    =coeffs$Ca_s3   [i]
  Cb_s3    =coeffs$Cb_s3   [i]
  Amax_s4  =coeffs$Amax_s4 [i]
  Src_s4   =coeffs$Src_s4  [i]
  Ca_s4    =coeffs$Ca_s4   [i]
  Cb_s4    =coeffs$Cb_s4   [i]
  sigmaS1  =coeffs$sigmaS1 [i]
  tauS1    =coeffs$tauS1   [i]
  sigmaST1 =coeffs$sigmaST1[i]
  sigmaS2  =coeffs$sigmaS2 [i]
  tauS2    =coeffs$tauS2   [i]
  sigmaST2 =coeffs$sigmaST2[i]
  sigmaS3  =coeffs$sigmaS3 [i]
  tauS3    =coeffs$tauS3   [i]
  sigmaST3 =coeffs$sigmaST3[i]
  sigmaS4  =coeffs$sigmaS4 [i]
  tauS4    =coeffs$tauS4   [i]
  sigmaST4 =coeffs$sigmaST4[i]
  
  if (Vs30> 1100) {
    siteclass = 0
    Imf = 0
  }  else if (Vs30 > 600) {
    siteclass = 1
    Imf = 0.91
  }  else if (Vs30 > 300) {
    siteclass = 2
    Imf = 1.023
  }  else if (Vs30 > 200) {
    siteclass = 3
    Imf = 1.034
  }  else{
    siteclass = 4
    Imf = 0.737
  }
  
  if (M > mc) {
    Cm = mc  
  }  else {
    Cm = M
  }
  
  # //calculate the hard rock motion first
  rij = R + exp(c1 + 1.151 * Cm)
  
  if (M > mc) {
    fmSL = bSL*Zh + cSL1*mc + cSL2*(mc - msc)^2 + dSL*(M - mc)
  }    else {
    fmSL = bSL*Zh + cSL1*M + cSL2*(M - msc)^2
  }
  
  if (Zh<50) {
    qslh = eSLH * 0
  }  else {
    qslh = eSLH * (0.02 * Zh - 1)
  } 
  lnSr = fmSL + gSL * log(rij) + 
    gSLL * log(R + 200) + eSL * R + qslh * R + eSLv * Xv + gamma
  Sr = exp(lnSr - log(AmSC1))
  Sreff = Sr * Imf
  
  # //calculate the site term
  alpha = 2
  beta = 0.6
  theta=10
  if (siteclass==1) {
    Anmax = AmSC1
    SreffC = Imf*Src_s1
    Sf = Anmax/Amax_s1
    if (Anmax >= 1.25) {
      Snc = (exp((log(Anmax) * log(SreffC ^ alpha + beta) - log(Sf) * log(beta)) / log(Amax_s1)) - beta) ^ (1/alpha)
      
    } else{
      if (Ca_s1 != 0) {
        temp_up = Ca_s1 * (alpha - 1) * log(beta) * log(theta * beta) - log(theta) * (Cb_s1 + log(Sf))
        temp_under = Ca_s1 * (alpha * log(theta * beta) - log(beta))
        Snc = exp(temp_up/temp_under)
        
      } else {
        Snc = 0
      }
    }      
    Smr = Sreff * (Snc / SreffC) *SC1
    lnA = -log(Amax_s1) * 
      ((log(Smr^alpha + beta) - log(beta))/(log(SreffC^alpha + beta) - log(beta)))
  } else if (siteclass == 2) {

    Anmax=AmSC1*exp(S2)
    SreffC = Imf*Src_s2
    Sf = Anmax/Amax_s2
    if (Anmax >= 1.25) {
      Snc = (exp((log(Anmax)*log(SreffC^alpha + beta) - log(Sf)*log(beta))/log(Amax_s2))-beta)^(1/alpha)
    } else{
      if (Ca_s2 != 0) {
        Snc = exp((Ca_s2*(alpha-1)*log(beta)*log(theta*beta)-log(theta)*(Cb_s2+log(Sf)))/(Ca_s2*(alpha*log(theta*beta)-log(beta))))
      } else{
        Snc = 0
      }
    }
    Smr = Sreff*(Snc/SreffC)*SC2
    lnA = S2-log(Amax_s2)*((log(Smr^alpha + beta) - log(beta))/(log(SreffC^alpha + beta)-log(beta)))
  }  else if (siteclass==3) {
    Anmax = AmSC1*exp(S3)
    SreffC = Imf*Src_s3
    Sf = Anmax/Amax_s3
    if (Anmax >= 1.25) {
      Snc = (exp((log(Anmax)*log(SreffC^alpha+beta)-log(Sf)*log(beta))/log(Amax_s3))-beta)^(1/alpha)
    } else{
      if (Ca_s3 != 0) {
        Snc = exp((Ca_s3*(alpha-1)*log(beta)*log(theta*beta)-log(theta)*(Cb_s3+log(Sf)))/(Ca_s3*(alpha*log(beta*theta)-log(beta))))
      }  else{
        Snc = 0
      } 
    } 
    Smr = Sreff*(Snc/SreffC)*SC3
    lnA = S3-log(Amax_s3)*((log(Smr^alpha+beta)-log(beta))/(log(SreffC^alpha+beta)-log(beta)))
    
  }  else if (siteclass == 4) {
    Anmax = AmSC1*exp(S4)
    SreffC = Imf*Src_s4
    Sf = Anmax/Amax_s4
    if (Anmax >= 1.25) {
      Snc = (exp((log(Anmax)*log(SreffC^alpha+beta)-log(Sf)*log(beta))/log(Amax_s4))-beta)^(1/alpha)
    } else {
      if (Ca_s4 != 0) {
        Snc = exp((Ca_s4*(alpha-1)*log(beta)*log(theta*beta)-log(theta)*(Cb_s4+log(Sf)))/(Ca_s4*(alpha*log(theta*beta)-log(beta))))
      }   else {
        Snc = 0
      }
    }
    Smr = Sreff*(Snc/SreffC)*SC4
    lnA = S4 - log(Amax_s4)*((log(Smr^alpha+beta)-log(beta))/(log(SreffC^alpha+beta)-log(beta)))
  }  else if (siteclass == 0) {
    lnA = -log(AmSC1)
  }
  
  # //calculate ground motions
  lnY = lnSr + lnA
  
  # //standard deviation
  # if (siteclass  == 1) {
  #   sigmaLnY = sigmaST1
  #   sigma = sigmaS1
  #   tau = tauS1
  # }  else if (siteclass  == 2) {
  #   sigmaLnY = sigmaST2
  #   sigma = sigmaS2
  #   tau = tauS2
  # } else if (siteclass  == 3) {
  #   sigmaLnY = sigmaST3
  #   sigma = sigmaS3
  #   tau = tauS3
  # } else if (siteclass  == 4) {
  #   sigmaLnY = sigmaST4
  #   sigma = sigmaS4
  #   tau = tauS4
  # } else {
  #   sigmaLnY = sT
  #   sigma = sigma0
  #   tau = tau0
  # }
  ## the total sigma does not vary with site classes 
  sigmaLnY = sT
  sigma = sigma0
  tau = tau0
  
  Sa = exp(lnY)
  sigmatotal = sigmaLnY
  
  temp <-array(c(Sa, sigmatotal, sigma, tau))
  
  return (temp)
}
####
Zhao16SZ_Inslab.Cal <- function (ip, M, R, Zh, Vs30, Xv,AddMedian)
{
  period <- array(c(0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08, 0.09,0.1,
                    0.12,0.14,0.15,0.16,0.18,0.2,0.25,0.3,0.35,0.4,0.45,
                    0.5, 0.6,0.7,0.8,0.9, 1,1.25,1.5,2,2.5,3, 3.5,4,4.5,5))
  if (length(which(period == ip)) == 0) {
    i_lo <- sum(period<ip)   ########
    if (i_lo>=length(period))
    {
      i_lo <- length(period)
      T_lo <- period[i_lo]
      Sa_sigma <- Zhao16SZ_Inslab.subCal(T_lo,M,R, Zh,Vs30,Xv)
    }
    else
    {
      T_lo <- period[i_lo]
      T_hi <- period[i_lo + 1]      
      Sa_sigma_hi <- Zhao16SZ_Inslab.subCal(T_hi,M,R, Zh,Vs30,Xv)
      Sa_sigma_lo <- Zhao16SZ_Inslab.subCal(T_lo,M,R, Zh,Vs30,Xv)
      #     x <- array(c(T_lo, T_hi))
      #     Y_Sa <- array(c(Sa_sigma_lo[1], Sa_sigma_hi[1]))
      #     Y_sigma <- array(c(Sa_sigma_lo[2], Sa_sigma_hi[2]))
      Sa <- interpolate(ip, T_lo, T_hi, Sa_sigma_lo[1], Sa_sigma_hi[1])  ##(x,Y_Sa,T)    #########
      sigmatotal <- interpolate(ip, T_lo, T_hi, Sa_sigma_lo[2], Sa_sigma_hi[2])
      sigma <- interpolate(ip, T_lo, T_hi, Sa_sigma_lo[2], Sa_sigma_hi[2])
      tau <- interpolate(ip, T_lo, T_hi, Sa_sigma_lo[2], Sa_sigma_hi[2])
      
      Sa_sigma <- array(c(Sa, sigmatotal, sigma, tau))
    }
  }
  if (length(which(period == ip)) > 0) {
    Sa_sigma <- Zhao16SZ_Inslab.subCal(ip,M,R, Zh,Vs30,Xv)
  }
  Sa_sigma[1] = Sa_sigma[1]+AddMedian
  return (Sa_sigma)
}



#####
Zhao16SZ.itr <- function (list.mag, list.rupdist, list.hypodist, list.p, list.zh, 
                        flag.source, Vs30=761, Xv=0, AddMedian=0) {
  m <- length(list.mag)
  n <- length(list.rupdist)
  le <- length(list.p)
  he = length(list.zh)
  output.Sa <- array(NA, dim = c(n, m, le, he))
  output.Sd <- array(NA, dim = c(le))
  output.sigma <- array(NA, dim = c(le))
  output.tau <- array(NA, dim = c(le))
  
  
  for (hh in 1:he) {
    for (j in 1:m)  {
      for (k in 1:n)  {
        for (t in 1:le) {
          
          if(list.p[t]>5) {
            results=c(NA,NA,0,0)
          } else {
            if (flag.source=="interface") {
              list.dist = list.rupdist
              dist.temp <- list.dist[k]
              results <- Zhao16SZ_Interf.Cal(list.p[t], list.mag[j], 
                                           dist.temp, list.zh[hh], 
                                           Vs30, Xv, AddMedian)
            }
            if (flag.source=="inslab") {
              list.dist = list.hypodist
              dist.temp <- list.dist[k]
              results <- Zhao16SZ_Inslab.Cal(list.p[t], list.mag[j], 
                                             dist.temp, list.zh[hh], 
                                             Vs30, Xv,AddMedian)
            }
          }
          output.Sa[k,j,t,hh] <- results[1] #*980
          output.Sd[t] <- results[2]
          output.sigma[t] = results[3]
          output.tau[t] = results[4]
        }
      }
    # print(h.list[hh])
    }
  } 
  output.Sd = cbind(output.Sd, output.sigma, output.tau)
  return(list(output.Sa, output.Sd))
}

