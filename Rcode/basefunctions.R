#
# Basic functions for single-population calculations, used throughout all chapters
#

#
# Global flag:
#
bRecalcExpensiveFunctions <- TRUE  # Use pre-calculated values for heavy functions and previously downloaded databases from Fishbase

library(deSolve)

# Length to weight calculation
weight <- function(length)
{
  p <- baseparameters()
  length^3*p$c
}

# Weight to length calculation
weight2length <- function(weight)
{
  p<-baseparameters()
  (weight/p$c)^(1/3)
}

# Calculate the A-parameter from observations of K and Linf:
calcA <- function(K,Linf, p=baseparameters())
{
  3*p$c^0.25 * p$etaM^(-1/12) * K * Linf^(3/4)
  #-p$c^(1-p$n) * p$etaM^(1-p$n)/((1-p$n)*log(1-p$etaM^(1/3))) *
  #  K * Linf^(3*(1-p$n))
}
# Calculate K from A and Linf:
calcK <- function(A, Linf, p=baseparameters())
{
  A * Linf^(-3/4) / (3*p$c^0.25 * p$etaM^(-1/12) )
  #p$A*(-p$c^(1-p$n) * p$etaM^(1-p$n)/((1-p$n)*log(1-p$etaM^(1/3))))^(-1) *
  #  Linf^(-3*(1-p$n))
}

# The switching function
psi <- function(z, u=5)
{
  (1+z^(-u))^(-1)
}

# Somatic growth rate:
growth <- function(p, w)
{
  p$A*w^p$n*( 1 - psi(w/(p$etaM*p$W)) * (w/p$W)^(1-p$n)) 
}

# Calculate weight at age
weightatage <- function(p, ages=seq(0,5*ageMaturation(p$W, p),length.out=500))
{
  out <- ode(y=p$w0, times=ages, func=function(ages,w,p) list(growth(p,w)), parms=p)
  data.frame(age=ages, w=out[,2])
}

# Mortality
mortality <- function(p, w)
  p$a*p$A*w^(p$n-1) + p$funcFishing(w,p)

# Investment in reproduction
calcKr <- function(p)
{
  p$epsEgg * p$A * p$W^(p$n-1)
}

# Reproductive output:
calcRegg <- function(p)
{
  calcKr(p)/p$w0
}

# Age at maturation
ageMaturation <- function(W, p=baseparameters()) 
{
  p$etaM^(1-p$n) / (p$A*(1-p$n)) * W^(1-p$n)
  #log(1-p$etaM^0.25*p$etaA) / (-0.25*p$etaA*p$A ) * W^(1-p$n)
}

# Calc spectrum per recruit numerically using the full growth equation:
spectrum <- function(p, nGrid=1000)
{
  w <- 10^(seq( log10(p$w0), log10(p$W), length.out = nGrid+1))
  w <- w[1:nGrid]
  g <- growth(p,w)
  mu <- mortality(p,w)
  
  # solve
  a <- mu*w/g
  D <- w[2]/w[1]
  NprR <- cumprod( D^(-mu/g*w) ) / g
  
  N <- calcAll(p, data.frame(w=w, g=g, mu=mu, NprR=NprR))
  N
}

# Analytical solution to the spectrum per recruit using von B growth:
spectrumana <- function(p, nGrid=1000)
{
  w <- 10^(seq( log10(p$w0), log10(p$W), length.out = nGrid+1))
  w <- w[1:nGrid]
  g <- p$A*w^p$n * (1 - (w/p$W)^(1-p$n))
  mu <- p$a*p$A*w^(p$n-1)
  
  NprR <- 1/g[1]* (w/p$w0)^(-p$n-p$a) * (1 - (w/p$W)^(1-p$n))^(p$a/(1-p$n)-1)
  
  N <- calcAll(p, data.frame(w=w, g=g, mu=mu, NprR=NprR))
  N
}

# Analytical solution using only the juvenile growth equation:
spectrumsimple <- function(p, nGrid=1000 )
{
  w <- 10^(seq( log10(p$w0), log10(p$W), length.out = nGrid+1))
  w <- w[1:nGrid]
  g <- p$A*w^p$n
  mu <- p$a*p$A*w^(p$n-1)
  
  NprR <- 1/g[1]* (w/p$w0)^(-p$n-p$a)
  
  N <- calcAll(p, data.frame(w=w, g=g, mu=mu, NprR=NprR))
  N
}

# Add calculations of reproduction and recruitment to a spectrum frame:
calcAll <- function(p, N)
{
  N$SSBperR <- calcSSB(p,N)
  N$R0 <- calcR0(p,N)
  N$R <- 1-1/N$R0 # R/Rmax
  N$Rp <- N$R0-1  # Rp/Rmax
  N$YprR <- trapz(N$w, N$NprR * N$w * p$funcFishing(N$w, p))
  N
}

# SSB per R:
calcSSB <- function(p, N) 
{
  trapz(N$w, psi(N$w/(p$etaM*p$W))*N$w*N$NprR )
}

# R0 = Rp/R:
calcR0 <- function(p, N) {
  SSB <- calcSSB(p,N)
  p$epsR * calcKr(p) *SSB / p$w0
}

# Yield from fishing
calcYield <- function(p,spec)
{
  trapz(spec$w,  spec$w * spec$NprR * spec$R * p$funcFishing(spec$w, p))
}

# Reference points:
calcRefpoints <- function(p) 
{
  maxF <- 50 # The maximum F that the function will deal with
  # Fmax:
  opt2 <- function(F) 
  {
    p$F <- F
    spec <- spectrum(p)
    spec$YprR[1]
  }
  Fmax <- optimize(opt2, interval=c(0,100), maximum=TRUE)$maximum
  # Bmax:
  p$F <- Fmax
  spec <- spectrum(p)
  Bmax <- max(0, spec$SSBperR[1] * spec$R[1])
  
  # Check whether we're crashed at F=0:
  p$F <- 0
  spec0 <- spectrum(p)
  B0 <- max(0, spec0$SSBperR[1] * spec0$R[1])
  if (spec0$R[1] <=0)
  {
    Fmsy <- 0
    Ymsy <- 0
    Fcrash <- 0
    Flim <- 0
    Bmsy <- B0
    Blim <- B0
  }  else
  {
    # Fmsy:
    opt <- function(F) 
    {
      p$F <- F
      spec <- spectrum(p)
      yield <- calcYield(p,spec)
      yield
    }
    Fmsycalc <- optimize(opt, interval=c(0, maxF), maximum=TRUE)
    Fmsy <- Fmsycalc$maximum
    Ymsy <- Fmsycalc$objective
    
    # Bmsy:
    p$F <- Fmsy
    spec <- spectrum(p)
    Bmsy <- spec$SSBperR[1] * spec$R[1]
    
    # Flim:
    opt3 <- function(F) 
    {
      if (F>maxF)
        return(0)
      else
      {
        p$F <- F
        spec <- spectrum(p)
        #cat(F,",",spec$R[1],"\n")
        return(spec$R[1]-0.5)
      }
    }
    
    # Check whether the population's crashed at Fmax; if not set Flim and Fcrash to Fmax
    pFmax <- p
    pFmax$F <- maxF
    specFmax <- spectrum(pFmax)
    if (specFmax$R[1]>0) {
      Flim <- maxF
      Blim <- specFmax$SSBperR[1] * specFmax$R[1]
      Fcrash <- maxF
    }
    else
    {
      
      if (spec0$R[1]<=0.5)
      {
        Flim <- 0
        Blim <- 0
      }
      else
      {
        Flim <- uniroot(opt3, interval=c(0,maxF))$root
        
        # Blim
        p$F <- Flim
        spec <- spectrum(p)
        Blim <- spec$SSBperR[1] * spec$R[1]
      }
      # Fcrash
      p$F<-0
      spec<-spectrum(p)
      if (spec$R[1]<=0)
        Fcrash <- 0
      else
      {
        target <- function(F) 
        {
          if (F>maxF)
            0
          else {
            p$F <- F
            spec <- spectrum(p)
            spec$R[1]
          }
        }
        Fcrash <- uniroot(target, interval=c(0,maxF))$root
        
      }
    }
  }
  list(Fmsy=Fmsy, Ymsy=Ymsy, Fcrash=Fcrash, Fmax=Fmax,
       Flim=Flim, Bmsy=Bmsy, Blim=Blim, W=p$W, Bmax=Bmax, B0=B0)
}
#
# Population growth rate from analytical approximation #1
#
calc_rana1 <- function(W,p=baseparameters()) {
#  (p$A*(-1 + p$n)*p$w0^p$n*W^p$n*
#     log((p$a*p$w0^(1 - p$a)*W^(-1 + p$a))/(p$epsEgg*p$epsR)))/
#  (p$w0^p$n*W - p$w0*W^p$n)
  p$A*(1-p$n) * (W^(1-p$n) - p$w0^(1-p$n))^(-1) * 
    ( (1-p$a)*log(W/p$w0)+log(p$epsEgg*p$epsR))
}
#
# Population growth rate from analytical approximation #2
#
require(gsl) # For Lambert W function
calc_rana2 <- function(W, p=baseparameters()) { 
  n <- p$n
  a <- p$a
  w0 <- p$w0
  epsilon <- p$epsEgg*p$epsR
  (-(a*p$A*W^n*w0) +  p$A*W*w0^n*(a + (-1 + n)* 
                                    lambert_W0(((W^n*w0 - W*w0^n)* 
                                                  (exp(a - a*W^(-1 + n)*w0^(1 - n))*p$A^(1 - n)*
                                                     W^((a - n)*(-1 + n))*w0^(-1 + a + n - a*n)*epsilon^(1 - n))^ (1/(1 - n)))/((-1 + n)*p$A*W^n*w0^n))))/(W*w0 - W^(2 - n)*w0^n)
}
#
# Functions for various types of fishing gear:
#
fishingTrawl <- function(w,pp) {
  muF <- matrix(0, ncol=length(w), nrow=length(pp$W))
  for (i in 1:length(pp$W)) {
    muF[i,] <- pp$F[i]*psi(w/(pp$etaF*pp$W[i]), u=3)   # trawl
    if (!is.null(pp$etaFF))
      muF[i,] <- muF[i,] * (1-psi(w/(pp$etaFF*pp$W[i]), u=3))  # Possible upper cut-off
  }   
  muF
} 

fishingKnifeedge <- function(w,pp) 
  pp$F*(w>pp$W*pp$etaF)

fishingGillnet <- function(w,pp)
  pp$F*exp(-(log(w/(pp$W*pp$etaF)))^2/(1.5^2))

fishingNarrow <- function(w,pp)
  pp$F*exp(-(log(w/(pp$W*pp$etaF)))^2/.05)

