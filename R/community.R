#
# Code for the community model calculations
#

source("R/basetools.R")
source("R/basefunctions.R")
source("R/baseparameters.R")

require(limSolve)
require(reshape)
require(tictoc)
require(deSolve)
require(pracma)

#
# Base parameters for the trait-based model:
#
paramTraitbasedmodel <- function(
  nSpecies = 9,
  W=10^seq(log10(4), 5, length.out = nSpecies),
  facRmax = 0.25) 
{
  param <- baseparameters()
  
  param$bVerbose <- FALSE
  
  param$W <- W
  param$nSpecies <- length(param$W)
  param$GridExpand <- 0.1
  param$wrec <- param$w0
  param$wMax <- max(param$W)
  param$Tend <- 200
  param$dt <- 0.1
  param$stepSave <- 20
  #
  # Mortality
  #
  param$mu0Coefficient <- 1.2
  param$mu0Exponent <- -0.25
  param$muS0 <- 1
  param$F <- rep(0, param$nSpecies)
  param$SizeBasedMortCoef <- 0
  #
  # Width of feeding preference
  #
  param$sigma <- 1.3 # Increased from the default value (1) p$Ato represent a spread in different species' prey size preference
  #
  # Resource parameters
  #
  param$KR <- 0.005 # Carrying capacity
  param$lambdaR <- -2-param$q+param$n # Exponent
  param$rR <- 4 # Productivity
  param$wRstart <- param$wrec / param$beta^2
  param$wRend <- min(param$W)/2
  param$thetaFish <- 1 # Preference for eating fish
  #
  # Calculate Rmax:
  #
  param$aTheo <- param$f0*param$h / param$A * param$beta^(2*param$n-param$q-1) *
    exp( (2*param$n*(param$q-1)-param$q^2+1)*param$sigma^2/2)
  fac = param$W[2]/param$W[1]# Assume that the Winf's are exponentially distributed
  deltaW <- param$W*(sqrt(fac)-1/sqrt(fac)) # Assume that the Winf's are exponentially distributed
  param$facRmax <- facRmax
  param$Rmax <- facRmax * param$KR * param$A *
    param$W^(2*param$n-param$q-3+param$aTheo) * param$wrec^(-param$aTheo) * deltaW
  param$wMax <- 2*max(param$W)
  #
  # Calc. gamma:
  #
  alphae = sqrt(2*pi)*param$sigma*param$beta^(-param$lambdaR-2)* 
    exp((-param$lambdaR-2)^2*param$sigma^2/2);
  # Note: we add DeltaF0Expected to the feeding level under the expectation that the average feeding level
  # is lower than the one for which gamme is calculated:
  DeltaF0Expected <- 0.05
  param$gamma <- (DeltaF0Expected+param$f0)*param$h / (alphae*param$KR*(1-param$f0-DeltaF0Expected)) #2.9190e+03
  #
  # Do simulation with infinite fixed resource spectrum:
  #
  # param$rR <- 1000
  # param$wRend <- param$wMax
  # param$Rmax <- 0.00001*param$Rmax
  # param$SizeBasedMortCoef <- param$a*param$A
  # param$mu0Coefficient <- 0
  # param$F <-  0*param$F
  
  param
}
#
# Base parameters for the single-species model embedded in a static community:
#
paramConsumerResourcemodel <- function(W, facRmax) {
  param <- paramTraitbasedmodel(W=W)
  param$sigma <- 1 # Assume a single species with narrow feeding preference
  
  # Fix maximum recruitment. Set to high value to avoid the SR-relation
  param$Rmax <- facRmax * param$A * param$KR * param$w0^-param$a *
    W^(-2-param$q+2*param$n+param$a)
  
  param$wRend <- W  # Make an infinite resource spectrum
  
  param$SizeBasedMortCoef <- param$A*param$a # Metabolic mortality
  param$mu0Coefficient <- 0 # No constant background mortality
  
  # Fix gamma to use non-corrected f0
  alphae = sqrt(2*pi)*param$sigma*param$beta^(-param$lambdaR-2)* 
    exp((-param$lambdaR-2)^2*param$sigma^2/2);
  param$gamma <- param$f0*param$h / (alphae*param$KR*(1-param$f0)) 
  
  param$thetaFish <- 1 # Cannibalism
  
  return(param)
}
#
# Core function to run the community model
#
runspectrum <- function(param, prev_run=NULL)  
{
  # --------------------------------------
  # Set up grids
  # --------------------------------------
  tic()
  # Fish grid:
  grid = makegrid(param$wrec, param$wMax, param$GridExpand)
  w = grid$w
  dw = grid$dw
  nGrid = grid$nGrid
  
  # Resource grid:
  grid = makegrid(param$wRstart, param$wRend, param$GridExpand)
  wR = grid$w
  dwR = grid$dw
  nGridR = grid$nGrid
  
  # --------------------------------------
  # Init parameters
  # --------------------------------------
  nSpecies = param$nSpecies  # No. of species
  nIte = param$Tend/param$dt # No. of iterations
  nSave = nIte/param$stepSave
  W = param$W # Asymptotic size
  Imax = param$h * w^param$n # Maximum consumption
  Clearance = param$gamma * w^param$q # Clearance rate
  wMature = param$etaM * W # Size at maturation
  # Allocation to reproduction:
  psiRepro <- matrix(data=0, nrow=nSpecies, ncol=nGrid)
  for (i in 1:nSpecies) {
    wRel <- w/W[i]
    psiRepro[i,] = wRel^(1-param$n) * psi(wRel/param$etaM)
    psiRepro[i, wRel>=1] <- 1
  }
  mu0 = param$mu0Coefficient * W^param$mu0Exponent # Background mortality
  muF = param$funcFishing(w, param) # Fishing mortality
  # --------------------------------------
  # Initial conditions
  # --------------------------------------
  NR = rep(0, nGridR)
  N = matrix(data=0, nrow=nSpecies, ncol=nGrid) 
  NinfR <- param$KR * wR^param$lambdaR # Resource spectrum at carrying capacity
  if (is.null(prev_run)) {
    # Set up from scratch:
    NR <- NinfR
    for (iSpecies in 1:nSpecies) {
      N[iSpecies, ] = 0.01*param$Rmax[iSpecies] * w^(-param$a - param$n) * (1 - (w/W[iSpecies])^(1-param$n))^((param$a/(1-param$n))-1)
      N[iSpecies, w>=W[iSpecies]] = 0
    }
  } else {
    # Use previous run:
    NR <- prev_run$resource$NR
    N <- as.matrix(prev_run$N[prev_run$nSave,,])
    if (nSpecies==1)
      N <- t(N)
  }
  #
  # Set up predation kernel
  #
  predkernelR = matrix(data=0, ncol=nGridR, nrow=nGrid)
  predkernel = matrix(data=0, nrow=nGrid, ncol=nGrid)
  
  for (j in 1:nGrid) {  # Loop over predator sizes
    predkernelR[j,] = exp(-0.5*(log(w[j]/(param$beta*wR))/param$sigma)^2)
    predkernel[j,] = param$thetaFish * exp(-0.5*(log(w[j]/(param$beta*w))/param$sigma)^2)
  }
  
  zeros <- rep(0,nSpecies*nGrid*nIte)
  #Nsave = data.frame(iteration=zeros, species=zeros, w=zeros, N=zeros) #array(data=0, dim=c(nIte, nSpecies, nGrid)) # Dimensions: iteration, species, size class
  Nsave = array(0, c(nSave, nSpecies, nGrid))
  SSB <- matrix(0, nSave, nSpecies)
  Rsave <- 0*SSB
  Rpsave <- 0*SSB
  R0save <- 0*SSB
  Ysave <- 0*SSB
  #
  # Other init stuff
  #
  A = matrix(data=0, nrow=nSpecies, ncol=nGrid-1)
  B = matrix(data=0, nrow=nSpecies, ncol=nGrid)
  C = rep(0,nGrid-1)
  S = matrix(data=0, nrow=nSpecies, ncol=nGrid)
  
  x <- 0 # Only used for noisy recruitment
  # --------------------------------------
  # Main loop
  # --------------------------------------
  for (ite in 1:nIte) {
    #
    # Calculate feeding level
    #
    Nc <- colSums(N)  # Community spectrum
    phiprey <- predkernelR %*% (NR*wR*dwR) + predkernel %*% (Nc*w*dw)   # Available food
    f <- Clearance*phiprey/(Imax + Clearance*phiprey)
    # 
    # Predation mortality:
    #
    tmp <- (1-f)*Clearance*Nc*dw
    muPR <- t(crossprod(predkernelR, tmp))  # ...on resource
    muP <- t( crossprod(predkernel, tmp) )  # ...on fish
    # 
    # Calculate growth and reproduction:
    #
    
    # Available energy
    Eavailable <- param$epsA*(f - param$fc) * Imax
    
    # Starvation mortality:
    muS <- matrix(0,nrow=1,ncol=nGrid)
    ix <- Eavailable<0
    muS[ix] <- -param$muS0*Eavailable[ix]/w[ix]
    Eavailable[ix] <- 0
    
    # Energy used for reproduction:
    Erepro <- psiRepro * (matrix(1,nSpecies,1)%*%t(Eavailable))
    # Reproductive output (note that I use the parameter "a" and not "aTheo" to calculate survival)
    Rp <- (param$wrec/param$w0)^(-param$a) *
      param$epsR*param$epsEgg*colSums(t(N*Erepro* (rep(1,nSpecies) %*% t(dw))))/param$w0
    # ...and the rest for somatic growth:
    g <- (1-psiRepro) * (matrix(1,nSpecies,1)%*%t(Eavailable))
    
    # Recruitment:
    R <- param$Rmax*Rp/(Rp + param$Rmax)
    
    # Noise on recruitment?
    if (!is.null(param$sigmaRecruitment)) {
      x <- x - param$thetaRecruitment*x*param$dt + param$sigmaRecruitment*sqrt(param$dt)*rnorm(1)
      R <- R*exp(x)
    }
    
    # Total mortality:
    mu <- matrix(1,nSpecies,1)%*%(muP + muS + param$SizeBasedMortCoef*w^(param$n-1)) + mu0 + muF 
    
    #
    # Set up matrix:
    #
    for (iSpecies in 1:nSpecies) {
      A[iSpecies, ] <- -g[iSpecies, 1:(nGrid-1)]*param$dt/dw[2:nGrid]
      B[iSpecies, ] <- 1 + param$dt*(g[iSpecies,]/dw + mu[iSpecies,])
      S[iSpecies, ] <- N[iSpecies, ]
      S[iSpecies, 1] <- N[iSpecies, 1] + R[iSpecies]*param$dt/dw[1]
    }  
    #
    # Invert matrix:
    #
    for (iSpecies in 1:nSpecies) 
      N[iSpecies,] <- Solve.tridiag(A[iSpecies,],B[iSpecies,],C,S[iSpecies,])

    #
    # Old way:
    #
    #N[iSpecies,1] = (N[iSpecies,1]+S[iSpecies,1])/B[iSpecies,1];
    #for (j in 2:nGrid) 
    #  N[iSpecies,j] = (S[iSpecies,j]-A[iSpecies,j-1]*N[iSpecies,j-1])/B[iSpecies,j]

    #
    # Update resource
    #
    tmp <- param$rR*NinfR / (param$rR + muPR)
    NR <- t( tmp - (tmp - t(NR)) * exp( -(param$rR + muPR)*param$dt) )
    #NR[nR<realmin] <- realmin # Fix points where the background has vanished
    # 
    # Save results
    #
    if ((ite %% param$stepSave) == 0) {
      iSave <- ite/param$stepSave
      for (iSpecies in 1:nSpecies) {
        #tmp = 1+(iSave-1)*nSpecies*nGrid + (iSpecies-1)*nGrid
        #Nsave[tmp:(tmp+nGrid-1),] = 
        #  data.frame(iteration=ite, species=iSpecies, w=w, N=N[iSpecies,] )
        Nsave[iSave,iSpecies,] <- N[iSpecies,]
        SSB[iSave, iSpecies] <- sum(w*N[iSpecies,]*dw*psi(psi(w/W[iSpecies]/param$etaM)))
        Ysave[iSave, iSpecies] <- sum(w*N[iSpecies,]*dw*muF[iSpecies,])
      }
      Rpsave[iSave, ] <- Rp
      Rsave[iSave, ] <- R / param$Rmax
      R0save[iSave, ] <- Rp/R
    }
  }
  #
  # Assemble the output:
  #
  result <- list(
    nSpecies = nSpecies,
    nIte = nIte,
    nSave = nSave,
    nGrid = nGrid,
    muF = muF, 
    param = param,
    R = as.data.frame(melt(SSB)),
    N = Nsave,
    g = g,
    w = w,
    Eavailable = Eavailable, 
    psi = psiRepro,
    fish = data.frame(
      w = w,
      dw = dw,
      f = f,
      g = g[nSpecies,],
      muP = t(muP),
      Nc = Nc
    ),
    resource = data.frame(
      wR = wR,
      muPR = t(muPR),
      NR = NR
    )
  )
  # Recruitment
  names(result$R) <- c("time","species","SSB")
  result$R$time <- result$R$time*param$dt*param$stepSave
  result$R$R <- melt(Rsave)$value
  result$R$Rp <- melt(Rpsave)$value
  result$R$Y <- melt(Ysave)$value
  result$R$R0 <- melt(R0save)$value
  
  if (param$bVerbose)
    toc()
  return(result)
}
#
# Helper function to make the computational grid
#
makegrid <- function(wStart, wEnd, expansion) 
{
  i = 1
  w = wStart
  while (w[i] < wEnd) {
    w[i+1] = w[i]*(1+expansion)
    i = i + 1
  }
  nGrid = i
  
  tmp = log(1+expansion);
  dw = tmp*wStart*exp((0:(nGrid-1))*tmp);
  
  list(w=w, dw=dw, nGrid=nGrid)
}
#
# calc. yield vs. F:
#
calcYield_vs_F <- function(param, res, iSpecies, nPoints=10) {
  p <- param
  F <- seq(0,1.5, length.out = nPoints)
  Yield <- 0*F
  wgtmat = 0*F
  SSB = 0*F
  ress = list()
  p$Tend <- 50
  resF <- res
  for (i in 1:length(F)) {
    p$F[iSpecies] <- F[i]
    resF <- runspectrum(p, resF)
    N <- resF$N[resF$nSave,iSpecies,]
    muF <- p$funcFishing(resF$fish$w,p)
    Yield[i] <- trapz(resF$w+0.5*resF$fish$dw, N*resF$fish$w*muF[iSpecies,]) #sum(N*resF$fish$w*resF$fish$dw*muF[iSpecies,])
    
    wgt = calcWeightAtAge(p, resF)
    wgtmat[i] = wgt$w[which(wgt$ages > ageMaturation(p$W,p))[1]]
    
    SSB[i] = sum(N*resF$fish$w*resF$fish$dw*psi(resF$fish$w/(p$etaM*p$W)))
    
    ress[[i]] = resF
  }
  data.frame(F=F, Yield=Yield, wmat=wgtmat, SSB=SSB)
}
#
# Weight at age
#
calcWeightAtAge <-function(p,res) {
  ages <- seq(0, 30, length.out = 100)
  out <- ode(y=p$w0, times=ages, func=function(t,w,dummy) list(interp1(res$fish$w, res$fish$g, w)))
  return( list(ages=ages, w=out[,2]) )
} 

plotTraitbasedmodel<- function(param, res) 
{
  #
  # Plot spectra:
  #
  
  # Assemble data frame:
  iPlot <- res$nSave # Iteration to plot
  N <- data.frame(N=rep(0,res$nGrid*param$nSpecies))
  for (i in 1:param$nSpecies) {
    ix = ((i-1)*res$nGrid+1):(i*res$nGrid)
    N$w[ix] <- res$fish$w
    N$N[ix] <- res$N[iPlot,i,]
    N$species[ix] <- i
  }
  # Plot:
  figSpec <- ggplot() +
    geom_line(data=res$resource, aes(x=wR, y=NR*wR), size=thick, linetype="dashed") +
    geom_line(data=res$fish, aes(x=w, y=Nc*w), size=thick) +
    geom_line(data=res$fish, aes(x=w, y=param$KR*w^(-1-param$q+param$n)), 
              color="red", size=thin, linetype="dotted") +
    geom_line(data=res$resource, aes(x=wR, y=param$KR*wR^(-1-param$q+param$n)), 
              color="red", size=thin, linetype="dotted")
  
  for (i in 1:res$nSpecies) {
    #figSpec <- figSpec + geom_line(data=subset(res$N, species==i & iteration==res$nIte), aes(x=w, y=N*w),color="blue")
    figSpec <- figSpec + geom_line(data=subset(N, species==i), aes(x=w, y=N*w),color="blue")
  }
  figSpec <- loglog(figSpec, ylim=c(1e-10,10), xlim=c(0.00001, 1e5)) +
    xlab("")
  #
  # Plot feeding level:
  #
  figF <- ggplot() +
    geom_line(data=res$fish, aes(x=w, y=f)) +
    geom_hline(yintercept = param$f0, size=thin, linetype="dotted") +
    geom_hline(yintercept = param$fc, size=thin, linetype="dotted")
  
  w = res$fish$w
  res$fish$loss = (res$fish$muP + param$SizeBasedMortCoef*w^(param$n-1)) /
    (res$fish$f*param$h*w^(param$n-1))
  figF <- figF +
    geom_line(data=res$fish, aes(x=w, y=loss), color="red")
  
  figF <- semilogx(figF)
  #
  # Plot mortality:
  #
  figMu <- ggplot() +
    geom_line(data=res$fish, aes(x=w, y=muP)) +
    geom_line(data=res$resource, aes(x=wR, y=muPR), linetype="dashed") +
    geom_line(dat=res$fish, aes(x=w, y=param$a*param$A*w^(param$n-1)), linetype="dotted") +
    geom_line(aes(x=res$fish$w, y=param$SizeBasedMortCoef*res$fish$w^(param$n-1)), linetype="dashed", color="blue")
  # Fishing mortality:
  muF <- as.data.frame(melt(param$funcFishing(res$fish$w, param)))
  names(muF) <- c("species","ix","muF")
  for (i in 1:res$nSpecies) 
    figMu <- figMu + geom_line(dat=subset(muF, species==i), aes(x=res$fish$w, y=muF))
  # Background mortality:
  figMu <- figMu + geom_point(aes(x=param$W, y=param$mu0Coefficient*param$W^param$mu0Exponent)) 
  figMu <- figMu + geom_line(data=res$fish, aes(x=w, y=param$SizeBasedMortCoef*w^(-0.25)))
  figMu <- loglog(figMu, xlim=c(1e-3, 1e5), ylim=c(0.05,7)) 
  #
  # Plot Biomasses vs time:
  #
  figB <- ggplot()
  for (i in 1:param$nSpecies)
    figB <- figB + geom_line(data=subset(res$R, species==i), aes(x=time, y=SSB))
  figB <- semilogy(figB)
  #
  # Plot R0:
  #
  figR0 <- ggplot() +
    geom_line(data=subset(res$R, time==param$Tend), aes(x=param$W, y=R0), size=thick) +
    geom_hline(yintercept = 1, size=thin, linetype="dotted")
  figR0 <- loglog(figR0)
  #
  # Save plot
  #
  pdf("baserun.pdf", width=doublewidth, height=3*height)
  multiplot(figSpec,figF,figMu,figB,figR0, cols=1, rows=2)
  dev.off()
  
  save_plot("baserun2.pdf",
    plot_grid(figSpec,figF,figMu,figB,figR0, labels=c("a","b","c","d","e"), align="v",ncol=1))
  
}
#
# Make a yield vs F curves for all species
#
plotYields <- function(param, res) {
  #calcYield_vs_F(param,res,4)
  Y <- data.frame(iSpecies=NULL, F=NULL, Y=NULL)
  for (i in 1:param$nSpecies) {
    Ysp <- calcYield_vs_F(param, res, i)
    Y <- rbind(Y, cbind(iSpecies=rep(i ,dim(Ysp)[1]), Ysp))
  }
  fig <- ggplot()
  for (i in 1:param$nSpecies)
    fig <- fig + 
      geom_line(data=subset(Y, iSpecies==i), aes(x=F, y=Yield/max(Yield)), size=i/3)
  ggsave("plotYields.pdf", fig)
  fig
}
#
# Compare several runs:
#
plotComparison <- function(param, res0, listRes) {
  xlim <- c(1e-3, 1e5)
  nSim <- length(listRes)
  #
  # Spectra
  #
  figS <- ggplot()
  for (i in 1:nSim) 
    figS <- figS +
    geom_line(data=listRes[[i]]$fish, aes(x=w, y=Nc/res0$fish$Nc), size=i/nSim*2)
  
  figS <- loglog(figS, ylim=c(0.01,100), xlim=xlim) +
    geom_hline(yintercept = 1, linetype="dashed") +
    ylab("Rel. change") +
    xlab("")
  #
  # Feeding level:
  #
  figF <- ggplot()
  for (i in 1:nSim) 
    figF <- figF +
    geom_line(data=listRes[[i]]$fish, aes(x=w, y=f), size=i/2)
  figF <- semilogx(figF, xlim=xlim) +
    geom_line(data=res0$fish, aes(x=w, y=f), linetype="dashed") +
    geom_hline(yintercept = param$f0, linetype="dotted", size=thin) +
    ylab("feeding level") +
    xlab("") +
    ylim(c(0,1))
  #
  # Predation mortality:
  #
  figM <- ggplot()
  for (i in 1:nSim) 
    figM <- figM +
    geom_line(data=listRes[[i]]$fish, aes(x=w, y=muP), size=i/2)
  figM <- semilogx(figM, xlim=xlim, ylim=c(0.05,7)) +
    geom_line(data=res0$fish, aes(x=w, y=muP), linetype="dashed") +
    geom_line(data=res0$fish, aes(x=w, y=param$A*param$aTheo*w^(param$n-1)), linetype="dotted", size=thin) +
    ylab(TeX("Mortality ($yr^{-1}$)")) +
    xlab("Weight (g)")    
  #
  # Recruitment:
  #
  figR0 <- ggplot()
  for (i in 1:nSim)
    figR0 <- figR0 + 
      geom_line(data=subset(listRes[[i]]$R, time==param$Tend), aes(x=param$W, y=Rp/R), size=i/2)
  figR0 <- figR0 +
    geom_line(data=subset(res0$R, time==param$Tend), aes(x=param$W, y=Rp/R), linetype="dashed") +
    geom_hline(yintercept = 1, size=thin, linetype="dotted")
  
  figR0 <- loglog(figR0)
  #
  # Fishing mortality
  #
  figFM <- ggplot()
  muF <- data.frame(sim=integer(), species=integer(), muF=double(), w=double())
  for (i in 1:nSim)
    for (j in 1:res0$nSpecies) 
      muF <- rbind(muF, data.frame(sim=i, species=j, w=listRes[[i]]$fish$w, muF = listRes[[i]]$muF[j,]))
  for (i in 1:nSim)
    for (j in 1:res0$nSpecies) 
      figFM <- figFM + geom_line(dat=subset(muF, sim==i & species==j), aes(x=w, y=muF), size=i/2)

  figFM <- semilogx(figFM)
  #
  # Yield
  #
  figY <- ggplot()
  for (i in 1:nSim)
    figY <- figY + 
    geom_line(data=subset(listRes[[i]]$R, time==param$Tend), 
              aes(x=param$W, y=Y), size=i/2)
  figY <- figY +
    geom_line(data=subset(res0$R, time==param$Tend), 
              aes(x=param$W, y=Y), linetype="dashed")

  figY <- semilogx(figY)
  
  fig <- plot_grid(figS, figR0, figF, figFM, figM, figY, nrow=3, align="v")
  #ggsave("Comparison.pdf", fig)
  fig
}
#
# Calculate reference points for all species:
#
calcRefpointsTraitbasedmodel <- function(param, res0, iSpec=1:param$nSpecies) {
  p <- param
  maxF <- 10
  refpoints <- data.frame(species=NULL, W=NULL, Fmsy=NULL, Ymsy=NULL)
  
  # Fmsy:
  opt <- function(F) 
  {
    p$F[iSpecies] <- F
    res <- runspectrum(p, res0)
    return( subset(res$R, time==p$Tend & species==iSpecies)$Y )
  }
  
  for (iSpecies in iSpec) {
    Fmsycalc <- optimize(opt, interval=c(0, maxF), maximum=TRUE, tol=0.03)
    refpoints <- rbind(refpoints, 
                       data.frame(species=iSpecies, W=param$W[iSpecies], Fmsy=Fmsycalc$maximum, Ymsy=Fmsycalc$objective))
    #print(refpoints)
  }  
  return(refpoints)
}

testCommunitySingleSpeciesFishing <- function(W=20000, cannibalism=1) {
  p <- paramCommunitymodel(W,100) # essentially without SR relation
  p$thetaFish=cannibalism # no cannibalism
  Y <- calcYield_vs_F(p, runspectrum(p), 1)
  
  p <- paramCommunitymodel(W,0.001) # dominated by SR relation
  p$thetaFish=cannibalism
  Ystd <- calcYield_vs_F(p, runspectrum(p), 1)
  
  defaultplot(mfcol=c(3,1))
  #
  # Yield vs F
  #
  defaultpanel(xlim=range(Y$F), ylim=c(0,1),
               xlab="F", ylab="Relative yield")
  lines(Y$F, Y$Yield/max(Y$Yield))
  lines(Ystd$F, Ystd$Yield/max(Ystd$Yield), lty=2)
  #
  # Growth reduction
  #
  defaultpanel(xlim=range(Y$F), ylim=c(0, 1.2*p$etaM*W),
               xlab="F", ylab("Weight at age of maturation"))
  lines(Y$F,Y$wmat)
  hline(p$etaM*W)
  #
  # Weight at age curves
  #
  defaultpanel(xlim=c(0,3*ageMaturation(W,p)), ylim=c(0,W),
               xlab("age"), ylab="Weight")
  
  p <- paramCommunitymodel(W,100)
  p$thetaFish=cannibalism # no cannibalism
  F = c(0, 0.3, 0.6, 0.9, 1.2, 1.5)
  for (i in 1:length(F)) {
    p$F[1] <- F[i]
    resF <- runspectrum(p)
    
    weightatage = calcWeightAtAge(p, resF)
    lines(weightatage$ages, weightatage$w)
  }
  vline(ageMaturation(W,p))
  
}


#
# Base run
#
baserun <- function() {
  param <- paramTraitbasedmodel()
  res <- runspectrum(param)
  plotTraitbasedmodel(param, res)
  res
}
