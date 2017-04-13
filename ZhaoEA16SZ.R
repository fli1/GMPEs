

mc = 7.1;
hc = 25;


if (Vs30> 1100) then
siteclass = 0;
Imf = 0;
else if (Vs30 > 600) then
siteclass = 1;
Imf = 0.91;
else if (Vs30 > 300) then
siteclass = 2;
Imf = 1.023;
else if (Vs30 > 200) then
siteclass = 3;
Imf = 1.034;
else
  siteclass = 4;
Imf = 0.737;
endif

//if (Fx == 1) then
//   R = Rcl;
//else
  //   R = Rcoe;
//endif
R = Rcl;
Zh = Zcoe; 

if (Zh > hc) then
deltah = 1;
else
  deltah = 0;
endif

if (Mw > mc) then
Cm = mc;
else
  Cm = Mw;
endif

//calculate the hard rock motion first

rij = 10 + R + exp(c1 + 1.151 * Cm);

if (deltah == 0) then
if (Mw > mc) then
fmintS = bint*Zh + gammaintS + cintS*mc + dint*(Mw-mc);
else
  fmintS = bint*Zh + gammaintS + cintS*Mw;
endif
lnSr = fmintS + gint * ln(rij) + gintLS*ln(R+200) + eintS*R + eintV*Xv + gammaint - ln(RSF);
else
  ;
if (Mw > mc) then
fmintD = bint*Zh + gammaintS + cintD*mc + dint*(Mw-mc);
else
  fmintD = bint*Zh + gammaintS + cintD*Mw;
end if
lnSr = fmintD + gint*ln(rij) + gintLD*ln(R+200) + eintV*Xv + gammaint - ln(RSF);
endif

Sr = exp(lnSr);
Sreff = Sr*Imf;

//calculate the site term
alpha = 2;
beta = 0.6;
theta=10;
if (siteclass==1) then
Anmax = RSF;
SreffC = Imf*Src_s1;
Sf = Anmax/Amax_s1;
if (Anmax >= 1.25) then
Snc = (exp((ln(Anmax)*ln(SreffC^alpha + beta) - ln(Sf)*ln(beta))/ln(Amax_s1))-beta)^(1/alpha);
else
  ;
if (Ca_s1 != 0) then
temp_up = Ca_s1 * (alpha - 1) * ln(beta) * ln(theta * beta) - ln(theta) * (Cb_s1 + ln(Sf));
temp_under = Ca_s1 * (alpha * ln(theta * beta) - ln(beta));
Snc = exp(temp_up/temp_under);
else
  Snc = 0;
end if
end if         
Smr = Sreff * (Snc / SreffC) *SC1;
lnA = -ln(Amax_s1) * ((ln(Smr ^ alpha + beta) - ln(beta)) / (ln(SreffC ^ alpha + beta) - ln(beta)));   
else if (siteclass == 2) then
if (Zh <= hc) then
Anmax=RSF*exp(S2);
ln_linearterm = S2;
else
  ;
Anmax = RSF*exp(S5);
ln_linearterm = S5;
endif
SreffC = Imf*Src_s2;
Sf = Anmax/Amax_s2;
if (Anmax >= 1.25) then
Snc = (exp((ln(Anmax)*ln(SreffC^alpha + beta) - ln(Sf)*ln(beta))/ln(Amax_s2))-beta)^(1/alpha);
else
  ;
if (Ca_s2 != 0) then
Snc = exp((Ca_s2*(alpha-1)*ln(beta)*ln(theta*beta)-ln(theta)*(Cb_s2+ln(Sf)))/(Ca_s2*(alpha*ln(theta*beta)-ln(beta))));
else
  ;
Snc = 0;
endif
endif
Smr = Sreff*(Snc/SreffC)*SC2;
lnA = ln_linearterm-ln(Amax_s2)*((ln(Smr^alpha + beta) - ln(beta))/(ln(SreffC^alpha + beta)-ln(beta)));
else if (siteclass==3) then
if (Zh <= hc) then
Anmax = RSF*exp(S3);
ln_linearterm = S3;
else
  ;
Anmax = RSF*exp(S6);
ln_linearterm = S6;
endif
SreffC = Imf*Src_s3;
Sf = Anmax/Amax_s3;
if (Anmax >= 1.25) then
Snc = (exp((ln(Anmax)*ln(SreffC^alpha+beta)-ln(Sf)*ln(beta))/ln(Amax_s3))-beta)^(1/alpha);
else
  ;
if (Ca_s3 != 0) then
Snc = exp((Ca_s3*(alpha-1)*ln(beta)*ln(theta*beta)-ln(theta)*(Cb_s3+ln(Sf)))/(Ca_s3*(alpha*ln(beta*theta)-ln(beta))));
else
  ;
Snc = 0;
endif
endif
Smr = Sreff*(Snc/SreffC)*SC3;
lnA = ln_linearterm-ln(Amax_s3)*((ln(Smr^alpha+beta)-ln(beta))/(ln(SreffC^alpha+beta)-ln(beta)));
else if (siteclass == 4) then
if (Zh <= hc) then
Anmax = RSF*exp(S4);
ln_linearterm = S4;
else
  ;
Anmax = RSF*exp(S7);
ln_linearterm = S7;
endif
SreffC = Imf*Src_s4;
Sf = Anmax/Amax_s4;
if (Anmax >= 1.25) then
Snc = (exp((ln(Anmax)*ln(SreffC^alpha+beta)-ln(Sf)*ln(beta))/ln(Amax_s4))-beta)^(1/alpha);
else
  ;
if (Ca_s4 != 0) then
Snc = exp((Ca_s4*(alpha-1)*ln(beta)*ln(theta*beta)-ln(theta)*(Cb_s4+ln(Sf)))/(Ca_s4*(alpha*ln(theta*beta)-ln(beta))));
else
  ;
Snc = 0;
endif
endif
Smr = Sreff*(Snc/SreffC)*SC4;
lnA = ln_linearterm - ln(Amax_s4)*((ln(Smr^alpha+beta)-ln(beta))/(ln(SreffC^alpha+beta)-ln(beta)));
else if (siteclass == 0) then
lnA = -ln(RSF);
endif

//calculate ground motions
if (deltah == 0) then
lnY = fmintS + gint*ln(rij) + gintLS*ln(R+200) + eintS*R + eintV*Xv + gammaint + lnA;
else
  ;
lnY = fmintD + gint*ln(rij) + gintLD*ln(R+200) + eintV*Xv  + gammaint + lnA;
endif

//standard deviation

if (siteclass  == 1) then
sigmaLnY = SC1ST;
// sig1 = sigma;
//  sig2 = tau;
else if (siteclass  == 2) then
sigmaLnY = SC2ST;
// sig1 = sigma;
// sig2 = tau;
else if (siteclass  == 3) then
sigmaLnY = SC3ST;
else if (siteclass  == 4) then
sigmaLnY = SC4ST;
else
  sigmaLnY = sT;
