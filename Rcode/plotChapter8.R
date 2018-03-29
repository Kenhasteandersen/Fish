#
# Chapter 8
#

source("Rcode/basetools.R")
source("Rcode/basefunctions.R")
source("Rcode/baseparameters.R")
source("Rcode/community.R")

require(cowplot)
require(gam)
require(numDeriv)

panelSMSspecies <- function(year=1990) {
  SMS <- read.csv("Data/SMS 2008.csv",sep=";")
  # Make a continuous time variable:
  SMS$time <- SMS$Year + 0.25*(SMS$Quarter-0.5)
  # Fix some names:
  SMS$Biomass <- SMS$Biomass..tonnes. 
  SMS$w <- SMS$Mean.Stock.weight..kg.*1000
  # Make log quantities for fitting:
  SMS$logB <- log10(SMS$Biomass)
  SMS$logw <- log10(SMS$w)
  # Fit a gam to each species:
  species <- unique(SMS$Species)
  SMSpred <- data.frame(species=NULL, w=NULL, B=NULL, W=NULL)
  grid <- makegrid(10^min(SMS$logw), 10^max(SMS$logw), .1)
  Winf <- rep(0, length(species))
  for (i in 1:length(species)) {
    # Fit:  
    speciesSMS <- subset(SMS, Species==species[i])
    fit <- gam(logB ~ s(logw) + s(time), data=speciesSMS)
    # Winf
    Winf[i] <- max(speciesSMS$w)
    #print(unique(speciesSMS$Species))
    #cat(unique(speciesSMS$Species), ", ", Winf[i],"\n")
    # Predict:
    ix <- grid$w>=10^min(speciesSMS$logw) & grid$w<=10^max(speciesSMS$logw)
    speciesDat <- data.frame(logw=log10(grid$w[ix]),
                             time=rep(year,sum(ix)))
    SMSpred <- rbind(SMSpred, data.frame(
      Species = rep(species[i], sum(ix)),
      w = 10^speciesDat$logw,
      dw = grid$dw[ix],
      B = 10^predict(fit, speciesDat),
      W = Winf[i]))
  }
  # Plot species:
  loglogpanel(ylim=c(10,1e8), xlim=c(1, 2e4),
              xlab="",ylab="Biomass spectrum (-)",
              label=TRUE)
  dat <- aggregate(SMSpred$B/SMSpred$dw, by=list(w=SMSpred$w), FUN=sum)
  lines(x=dat$w, y=dat$x, lwd=3)
  lines(x=grid$w, y=.5e7*grid$w^-1.05, lty=dotted)
  for (i in 1:length(species)) {
    ix <- SMSpred$Species == unique(SMSpred$Species)[i]
    lines(SMSpred$w[ix], SMSpred$B[ix]/SMSpred$dw[ix], col=grey(0.1+(i-1)/length(species)), lwd=2 )
  }
}


panelSMSgroups <- function(year=1990) {
  SMS <- read.csv("Data/SMS 2008.csv",sep=";")
  # Make a continuous time variable:
  SMS$time <- SMS$Year + 0.25*(SMS$Quarter-0.5)
  # Fix some names:
  SMS$Biomass <- SMS$Biomass..tonnes. 
  SMS$w <- SMS$Mean.Stock.weight..kg.*1000
  # Make log terms for fitting:
  SMS$logB <- log10(SMS$Biomass)
  SMS$logw <- log10(SMS$w)
  # Find Winf of each species:
  species <- unique(SMS$Species)
  for (i in 1:length(species)) 
    SMS$W[SMS$Species==species[i]] <- max(subset(SMS, Species==species[i])$w)
  # Fit a gam to each Winf group:
  Wgroup = c(10,100,1000,10000,100000)
  SMSpred <- data.frame(species=NULL, w=NULL, B=NULL, time=NULL, dw=NULL)
  grid <- makegrid(10^min(SMS$logw), 10^max(SMS$logw), .1)
  for (i in 2:length(Wgroup)) {
    # Fit
    groupSMS <- subset(SMS, W>Wgroup[i-1] & W<=Wgroup[i])
    fit <- gam(logB ~ s(logw) + s(time), data=groupSMS)
    # Predict:
    ix <- grid$w>=10^min(groupSMS$logw) & grid$w<=10^max(groupSMS$logw)
    years <- 1975:2008
    groupDat <- data.frame(logw = rep(log10(grid$w[ix]), length(years)),
                           time=rep(years, length.out=sum(ix)*length(years), each=sum(ix)))
    tmp <- data.frame(
      Wgrp = Wgroup[i],
      w = 10^groupDat$logw,
      time = groupDat$time,
      dw = grid$dw[ix],
      B = 10^predict(fit, groupDat)/grid$dw[ix])  # Biomass spectrum
    # Find min and max ranges
    
    tmp$Bmin <- rep(aggregate(tmp, by=list(w=tmp$w), FUN=min)$B, length(years))
    tmp$Bmax <- rep(aggregate(tmp, by=list(w=tmp$w), FUN=max)$B, length(years))
    SMSpred <- rbind(SMSpred, tmp)
  }
  
  # Plot Wgroups:
  loglogpanel(ylim=c(10,1e8), xlim=c(1, 2e4), label=TRUE, yaxis = FALSE,
              xlab="Weight (g)")
  
  for (i in 2:length(Wgroup)) {
    dat <- subset(SMSpred, Wgrp==Wgroup[i] & time==year)
    polygon(x=c(dat$w, dat$w[seq(length(dat$w),by=-1,1)]),
            y=c(dat$Bmin, dat$Bmax[seq(length(dat$w),by=-1,1)]),
            col=lightgrey, border=NA)
  }
  
  dat <- aggregate(SMSpred, by=list(w=SMSpred$w), FUN=sum)
  lines(x=unique(SMSpred$w), y=dat$B/length(years), lwd=3)
  
  # Power-law line:
  lines(grid$w, y=.5e7*grid$w^-1.05, lty=dotted)
  # Species spectra:
  for (i in unique(SMSpred$Wgrp)) {
    dat <- subset(SMSpred, time==year & Wgrp==i)
    lines(dat$w, dat$B, lwd=2)
  }
  
  #
  # Calculate extended sheldon:
  #
  # BW <- rep(0,4)
  # ix <- grid$w>10
  # for (i in 2:length(Wgroup)) {
  #   tmp <- subset( SMSpred, w>10 & Wgrp==Wgroup[i])
  #   BW[i-1] <- sum( tmp$B*tmp$dw )
  # }
  # fig2 <- ggplot() +
  #   geom_point(aes(x=Wgroup[2:5], y=BW))
  # fig2 <- loglog(fig)
  # fig2
}

panelTheory <- function() {
  param <- baseparameters()
  
  K <- .5e7
  n <- param$n
  q <- param$q
  a <- param$a
  
  w <- 10^seq(0, 5, length.out = 500)
  Nc <- K * w^(-2-q+n)
  
  W <- 10^seq(-1, 5, length.out = 700)
  dW <- 0*W
  dW[2:length(dW)] <- diff(W)
  dW[1] <- dW[2]
  Wgroup <- 10^(1:5)
  Ndat <- data.frame(w=0, W=0, N=0)
  names(Ndat) <- c("w","W","N")
  for (i in 2:length(Wgroup)) {
    Ngroup <- 0*w
    ixW <- W>Wgroup[i-1] & W<=Wgroup[i]
    for (j in 1:length(W)) 
      if (ixW[j]) {
        N <- 1/0.5383*K*W[j]^(2*n-q-3+a) * 
          w^(-n-a) * (1-(w/W[j])^(1-n))^(a/(1-n)-1)
        N[w>=W[j]] <- 0
        Ngroup <- Ngroup + N*dW[j]
      }
    Ndat <- rbind(Ndat, 
                  data.frame(w=w, W=Wgroup[i], N=Ngroup))
  }
  
  loglogpanel(ylim=c(10,1e8), xlim=c(1, 2e4), label=TRUE, yaxis = FALSE)
  lines(w, Nc*w, lwd=3)
  
  for (W in unique(Ndat$W)) {
    ix <- Ndat$W==W
    lines(Ndat$w[ix], Ndat$w[ix]*Ndat$N[ix], lwd=2)
  }
}

plotSpectra <- function() {
  defaultplothorizontal(nPanels=3) 
  # make room for title:
  oma <- par()$oma
  oma[top]<-2
  par(oma=oma) 

  panelSMSspecies()
  mtext(side=top, line=0.5, "Species")
  panelSMSgroups()
  mtext(side=top, line=0.35, TeX("$\\textit{W}_{\\infty}$ groups"))
  panelTheory()
  mtext(side=top, line=0.5, "Theory")
}


plotExtendedSheldon <- function() {
  # 
  # Data from Daan et al (2005):
  #
  length = exp( seq(from=2.75, to=5.75, by=0.5) ) # 
  Daan <- data.frame(W=0.01*length^3)
  factor = 3e7 # hand-fitted factor to calibrate between CPUE and numbers from the model
  Daan$CPUE1 = c(2.2, .95, 5.9, 5.8, 3.5, 3.2, -5.3)
  Daan$CPUE2 = c(5.5, 4.8, 7.2, 5.5, 6.7, 2, -3)
  Daan$CPUE3 = c(7, 4.2, 7.2, 6.2, 6.8, 3.8, NaN)
  #
  # Simulations with trait-based model
  #
  param <- paramTraitbasedmodel(W=10^seq(log10(4), 5, length.out = 27))
  param$F <- param$F*0
  res0 <- runspectrum(param)
  param$F <- 0.4 + param$F # Fishing with F=0.4 on all species
  resF <- runspectrum(param)#, res0)
  # calc Extended Sheldon:
  w0 <- 0.01*10^3  # Smallest fish retained by trawl
  ix <- (res0$w > w0)
  Nsummed0 <- rowSums(res0$N[res0$nSave, ,ix] * t(matrix(rep(res0$w[ix],param$nSpecies), ncol=param$nSpecies)))
  NsummedF <- rowSums(resF$N[resF$nSave, ,ix] * t(matrix(rep(resF$w[ix],param$nSpecies), ncol=param$nSpecies)))
  #
  # Theoretic results:
  #
  n <- param$n
  a <- param$a
  c <- 2*n-param$q-3+a
  W <- 10^seq(1,6,length.out=20)
  theo = data.frame(W=W, Nsummed = 1e13* (W^(c+2-n-a)/(c+2-n-a) - w0^(1-n-a)/(c+1)*W^(c+1) ))
  #
  # Plot:
  #
  
  ix <- param$W>10
  
  defaultplot()
  loglogpanel(xlim=c(10,1e6), ylim=c(1e5, 1e11),
              xlab='Asymptotic weight (g)',
              ylab="Abundance (a.u.)")
  points(Daan$W, factor*exp(Daan$CPUE1), pch=dots)
  points(Daan$W, factor*exp(Daan$CPUE2), pch=triangles)
  points(Daan$W, factor*exp(Daan$CPUE3), pch=squares)
  lines(Daan$W, factor*exp(Daan$CPUE1))
  lines(Daan$W, factor*exp(Daan$CPUE2))
  lines(Daan$W, factor*exp(Daan$CPUE3))
  
  lines(theo$W, theo$Nsummed, lty=dashed, lwd=3)
  lines(param$W[ix], y=1e14*Nsummed0[ix], lwd=3)
  lines(param$W[ix], y=1e14*NsummedF[ix], lwd=3, col=stdgrey)
}


plotGrowthF0 <- function() {
  require(deSolve)
  
  f0 <- c(0.4, 0.6, 0.8)  # The feeding levels to use for the plot
  
  ages = seq(0,50,length.out=500)
  dat <- data.frame(age=NULL, f0=NULL, w=NULL, R=NULL)
  for (i in 1:length(f0)) {
    #
    # Size at age
    #
    
    W = 2000
    p <- baseparameters()  
    p$W <- W
    etaM <- p$etaM
    p$A <- p$epsA*p$h*(f0[i]-p$fc)
    
    
    # Interval:
    out <- ode(y=p$w0, times=ages, func=function(t,w,p) list(growth(p,w)), parms=p)
    w <- out[,2]
    R <- p$epsEgg*p$A*w^p$n*psi(w/(p$etaM*p$W)) * (w/p$W)^(1-p$n)
    dat <- rbind(dat, 
                 data.frame(age=out[,1], f0=as.factor(f0[i]), w=w, R=R))
  }
  
  defaultplot()
  defaultpanel(xlim=c(0,35), ylim=c(0,2000),
               xlab="Age (years)", ylab="Weight (g)")
  # Weight:
  for (i in f0) 
    lines(dat$age[dat$f0==i], dat$w[dat$f0==i], lwd=1.5+1.5*(i==0.6))
  hline(y=p$W)
  hline(y=p$etaM*p$W)

  # Repro:
  for (i in f0) 
    lines(dat$age[dat$f0==i], dat$R[dat$f0==i], lwd=1.5+1.5*(i==0.6), col=stdgrey)
}


plotRmaxSensitivity <- function() {
  facRmax <- 2.5*c(0.03, 0.1, 0.3, 1.0)
  
  dat <- data.frame(w=NULL, muP=NULL, facRmax=NULL)
  dat2 <- data.frame(SSB=NULL, W=NULL, facRmax=NULL)
  res <- list()
  for (i in 1:length(facRmax)) {
    p <- paramTraitbasedmodel(W=10^seq(log10(4), 5, length.out = 20), facRmax=facRmax[i])
    res[[i]] <- runspectrum(p)
    dat <- rbind(dat, 
                 data.frame(w=res[[i]]$fish$w, muP=res[[i]]$fish$muP, facRmax=facRmax[i]))
    dat2 <- rbind(dat2, data.frame(
      SSB=subset(res[[i]]$R, time==p$Tend)$SSB, W=p$W, facRmax=facRmax[i]))
  }
  
  defaultplot()
  par(mar=par()$mar + c(0,0,0,5)) # space for legend
  loglogpanel(xlim=c(3,1e5), ylim=c(1e-5, 1e-2),
              xlab="Asymptotic size (g)",
              ylab = "SSB (g/$m^3$)")
  
  for (i in 1:length(facRmax)) {
    ix <- dat2$facRmax==facRmax[i]
    lines(dat2$W[ix], dat2$SSB[ix], lwd=0.75*i)
  }

  legend("right", bty="n", xpd=TRUE, 
         inset=c(-0.28,0), 
         legend=facRmax, 
         lwd=0.75*seq(1,5),
         title=TeX("$\\textit{K}_{Rmax}$"))
}


plotDynamicCommunity <- function() {
  param <- paramTraitbasedmodel()
  res <- runspectrum(param)
  wtot = sort(c(res$fish$w, res$resource$wR))  # weight classes spanning resource + fish
  xlim = c(1e-4, 1e5)
  #
  # Top frame: spectrum
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
  fig <- ggplot() +
    geom_line(aes(x=wtot, y=param$KR*wtot^(-1-param$q+param$n)), linetype="dashed") + # Theoretical spectrum
    geom_line(data=res$resource, aes(x=wR, y=NR*wR), size=thick, color="grey") + # Resource spectrum
    geom_line(data=res$fish, aes(x=w, y=Nc*w), size=thick) + # Fish community spectrum
    geom_line(data=N, aes(x=w, y=N*w, group=species)) # Asymptotic size groups
  fig <- loglog(fig, ylim=c(1e-8,10), xlim=xlim, bXaxis=FALSE, label="a") +
    ylab("Biomass spectrum (-)")
  #
  # Feeding level:
  #
  figF <- ggplot(data=res$fish) +
    geom_line(aes(x=w, y=f), size=thick) + # Feeding level
    geom_hline(yintercept = param$fc, linetype="dotted", size=thin) + # Critical feeding level
    geom_hline(yintercept = param$f0, linetype="dashed") +   # f0
    annotate("text", x=2e-4, y=0.12+param$fc, label="f[c]",parse=TRUE) +
    annotate("text", x=2e-4, y=0.12+param$f0, label="f[0]",parse=TRUE)
  figF <- semilogx(figF, xlim=xlim, ylim=c(0,1), bXaxis=FALSE, label="b") +
    ylab("feeding level")
  
  #
  # Mortality
  #
  figM <- ggplot() +
    geom_line(data=res$fish, aes(x=w, y=muP), size=thick) + # Predation mortality on fish
    geom_line(data=res$resource, aes(x=wR, y=muPR), size=thick) + # Predation mortality on resource
    geom_line(aes(x=wtot, y=param$A*param$aTheo*wtot^(param$n-1) ), linetype="dashed") +
    geom_point(aes(x=param$W, y=param$mu0Coefficient*param$W^param$mu0Exponent))
  figM <- loglog(figM, ylim=c(0.05,8), xlim=xlim, label="c") +
    xlab("Weight (g)") +
    ylab(TeX("Mortality ($yr^{-1}$)"))
  
  #
  # Recruitment
  #
  
  # Calc. "theoretical" recruitment from single-species model:
  p <- param
  W <- 10^seq(-1,6,length.out = 10)
  R0full <- 0*W
  
  for (i in 1:length(R0full)) {
    p$W <- W[i]
    
    N <- spectrum(p)
    R0full[i] <- N$R0[1]
  }
  
  figR <- ggplot() +
    geom_line(data=subset(res$R, time==param$Tend), aes(x=param$W, y=R0), size=thin) +
    geom_point(data=subset(res$R, time==param$Tend), aes(x=param$W, y=R0)) +
    geom_line(aes(x=W, y=R0full), linetype="dashed") +
    geom_hline(yintercept = 1, size=thin, linetype="dotted")
  figR <- loglog(figR, ylim=c(0.5,1000), xlim=c(1,1e5), label="d")+
    xlab(TeX("Asymptotic size, $\\textit{W}_{\\infty}$ (g)")) +
    ylab(TeX("$\\textit{R}_0$"))
  
  grid <- plot_grid(fig,figF,figM,figR, ncol=1, align="v",
                    rel_heights = c(2,1,1.3,1.3))
  ggsave("Chapter8/DynamicCommunity.pdf", grid, height=3*height, width=singlewidth)  
}

#
# Test for intermezzo:
#
plotDD <- function() {
  p <- paramCommunitymodel(W=1000, facRmax = 10000000000)
  p$rR <- 0.01
  result <- runspectrum(p)
  plotTraitbasedmodel(p,result)
}


plotAllChapter8 <- function() 
{
  pdfplot("Chapter8/spectra.pdf", plotSpectra, width=doublewidth, height = height)
  pdfplot("Chapter8/ExtendedSheldon.pdf", plotExtendedSheldon, width=singlewidth, height=height)
  pdfplot("Chapter8/GrowthF0.pdf", plotGrowthF0, width=singlewidth, height=height)
  pdfplot("Chapter8/RmaxSensitivity.pdf", plotRmaxSensitivity, width=1.5*singlewidth, height=height)

  plotDynamicCommunity()
  
  # Calc. approximation of relation between K, Linf and h:
  param <- paramTraitbasedmodel()
  cat("h =", 
      3*param$c^0.25* param$etaM^(-1/12) /(param$epsA*(param$f0-param$fc)),
      "K Linf^(3/4)\n")
  
}

