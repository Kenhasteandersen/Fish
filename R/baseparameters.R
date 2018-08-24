#
# Default values of all parameters
#

baseparameters <- function(W=NULL)
{
  p <- list()
  
  p$c <- 0.01  # weight-length relationship
  
  # Physiological parameters:
  p$W <- W
  p$n <- 3/4
  p$A <- 5.35 #4.39 # g^0.25 / yr
  p$epsA <- 0.6 # Assimilation efficiency
  p$epsEgg <- 0.22
  p$etaM <- 0.28
  p$f0 <- 0.6
  p$fc <- 0.2
  p$h <- p$A / (p$epsA*(p$f0-p$fc))
  
  # Pred-prey interaction parameters:
  p$beta <- 408
  p$sigma <- 1
  p$q <- 0.8
  
  # Reproduction:
  p$u <- 5
  p$w0 <- 0.001
  p$epsR <- 0.03 #0.015
  p$epsR_Elasmobranch = 0.3
  p$ElasmobranchEggMassRatio <- 0.0044
  
  # Mortality:
  
  p$a <- p$beta^(2*p$n-p$q-1) * exp((2*p$n-p$q-1)*(p$q-1)*p$sigma^2/2) * p$f0/(p$epsA*(p$f0-p$fc)) #0.425 # 0.34
  p$aMKconstant <- 1/4.83   # Relation between a and M/K
  
  # Fishing
  p$etaF = 0.05
  p$funcFishing = fishingTrawl
  p$F = 0
  p
}

