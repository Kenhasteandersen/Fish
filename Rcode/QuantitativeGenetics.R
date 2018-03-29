#
#  Code for quantitative genetics
#
library(pracma)  # for cumtrapz function

baseparamQG <- function(wm=1000, kr=NULL, p=baseparameters()) {
  p$wm <- wm
  
  # Set a default value of kr either from a specification of wm or W
  if (!is.null(p$W))
    p$kr <- p$A * p$W^(p$n-1)
  else
    if (is.null(kr))
      p$kr <- p$A * (p$wm/p$etaM)^(p$n-1)
  else
    p$kr <- kr
  p$W <- calcW(p)
  
  p$dFac <- 0.99  # factor used when evaluating numerical derivatives
  p$h2 <- 0.2     # Heritability
  p$cv <- 0.2     # Coefficient of variation
  p$fracSpawnerFishery <- 0 # Fraction of the F that's a spawner fishery
  p
}

calcW <- function(param)
  ( param$kr/param$A )^(1/(param$n-1))

dwdt <- function(p, w) 
  p$A * w^p$n - psi(w/p$wm)*p$kr*w

#
# Fitness (R0):
#
R <- function(param, F=0, nGrid=1000) {
  # Calc Winf
  W <- calcW(param)

  w = 10^seq(from=log10(param$w0), to=log10(W), length.out = nGrid+1)
  w = w[1:nGrid]
  # Calc. growth rate
  g <- dwdt(param,w)
  # Calc. mortality rates
  mup <- param$a*param$A*w^(param$n-1)
  param$F <- F
  muF <- (1-param$fracSpawnerFishery)*param$funcFishing(w,param) +
    param$fracSpawnerFishery*param$F*(w>param$wm)
  mu <- mup + muF
  # Calc. survival:
  S <- exp(-cumtrapz(w, t(mu/g)))
  S[length(S)] <- 0
  # Calc R
  param$epsEgg*param$kr*trapz(w, S*w*psi(w/param$wm)/g)/param$w0
}
#
# Gradients in R :
#
dRdwm <- function(p,F=0) {
  R0 <- R(p,F)
  p$wm <- p$dFac*p$wm
  R1 <- R(p,F)
  
  (R0-R1)/((1-p$dFac)*p$wm) / (0.5*(R0+R1))
}

dRdA <- function(p,F=0) {
  R0 <- R(p,F)
  p$A <- p$dFac*p$A
  R1 <- R(p,F)
  
  (R0-R1)/((1-p$dFac)*p$A) / (0.5*(R0+R1))
}

dRdkr <- function(p,F=0) {
  R0 <- R(p,F)
  p$kr <- p$dFac*p$kr
  R1 <- R(p,F)
  
  (R0-R1)/((1-p$dFac)*p$kr) / (0.5*(R0+R1))
}
# 
# Relative selection response
#
dSdt <- function(p,F=0) {
  # Generation time:
  T <- ageMaturation(p$W, p)
  
  # Traits:
  p$wm <- p$etaM*p$W
  p$kr <- p$A*p$W^(p$n-1)

  # Maturation:
  # no fishing:
  dRdwm0 <- dRdwm(p,0)
  # with fishing:
  dRdwmF <- dRdwm(p,F)
  dwmdt <- p$h2 * p$cv^2 * p$wm * (dRdwmF - dRdwm0)/T
  
  # A:
  dRdA0 <- dRdA(p,0)
  dRdAF <- dRdA(p,F)
  dAdt <- p$h2 * p$cv^2 * p$A * (dRdAF - dRdA0)/T

  # kr:
  dRdkr0 <- dRdkr(p,0)
  dRdkrF <- dRdkr(p,F)
  dkrdt <- p$h2 * p$cv^2 * p$kr* (dRdkrF - dRdkr0)/T
  
  # dWdt:
  p$wm <- p$wm*(1+dwmdt)
  p$A <- p$A*(1+dAdt)
  p$kr <- p$kr*(1+dkrdt)
  dWdt <- (calcW(p)-p$W)/p$W

  data.frame(dwmdt=dwmdt, dAdt=dAdt, dkrdt=dkrdt, dWdt=dWdt)
}
#
# Sel. response over a range of asymptotic sizes:
#
calcSelectionResponse <- function(p=baseparamQG(), 
                                 W=10^seq(log10(4),5,length.out=9), 
                                 F=0) {
  S <- data.frame(W = W)
  for (i in 1:length(W)) {
    p$W <- W[i]

    SS <- dSdt(p,F)
    
    S$dwmdt[i] <- SS$dwmdt
    S$dAdt[i] <- SS$dAdt
    S$dkrdt[i] <- SS$dkrdt
    S$dWdt[i] <- SS$dWdt
  }
  S
}
#
# Iterate selection over "years" years
#
calcIteratedSelectionResponse <- function(p=baseparamQG(), F=0.3, years) {
  for (year in 1:years){
    S <- calcSelectionResponse(p=p, F=F, W=calcW(p))
    p$wm <- p$wm + S$dwmdt*p$wm
    p$A <- p$A + S$dAdt*p$A
    p$kr <- p$kr + S$dkrdt*p$kr
  }
  p$W <- calcW(p)
  p$etaM <- p$wm/p$W
  return(p)
}
