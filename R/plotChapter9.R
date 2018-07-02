#
# Chapter 9
#

source("R/basetools.R")
source("R/basefunctions.R")
source("R/baseparameters.R")
source("R/community.R")

plotDaan <- function()
{
  dat <- read.csv("Data/Daan.dat", sep=",", header=FALSE)
  names(dat) <- c("Year", "CPUE")
  dat$Year <- 0.5+round(dat$Year-0.5) # Fix years to be on the same point
  dat$CPUE[seq(9,by=120/5,length.out=5)] <- 1
  dat$Lengths <- factor(floor(seq(1,5.99,length.out=120)), 
                        labels=c("<20 cm", "20-30 cm", "30-40 cm",
                                 "40-50 cm", ">50 cm"))
  
  defaultplot()
  par(mar=par()$mar + c(0,0,0,5)) # space for legend
  semilogypanel(xlim=c(1977.5, 2000.5), ylim=c(0.3, 3),
                xlab="Year", ylab="Scaled CPUE")
  vline(1985.5)
  for (i in 1:length(unique(dat$Lengths))){
    ix <- dat$Lengths == unique(dat$Lengths)[i]
    lines(dat$Year[ix], dat$CPUE[ix], lwd=0.5)
    points(dat$Year[ix], dat$CPUE[ix], 
          pch=dots, col=gray(0.15+(1-0.15)/5*(5-i)))
  }

  legend("right", bty="n", inset=c(-0.35,0), xpd=TRUE,
         legend=unique(dat$Lengths), pch=dots, 
         col=gray(0.15+(1-0.15)/5*(c(4,3,2,1,0))),
         title="Lengths")
}




panelComparison <- function(param, res0, fish, bYaxis=TRUE) {
  xlim <- c(1, 1e5)
  #
  # Fishing pressure
  #
  p <- param
  p$F <- 0*p$F + 1
  ww = subset(fish,sim==1)$w
  muF = param$funcFishing(ww,p)
  tmp = data.frame(sp = integer(),w=double(),muF=double())
  for (i in seq(1,param$nSpecies,by=3))
    tmp <- rbind(tmp, data.frame(sp=i, w=ww, muF=muF[i,]))
  figFish <- ggplot(data=tmp) +
    geom_line(aes(w, muF,group=sp)) +
    #geom_hline(yintercept = 1, linetype="dotted", size=thin)  +
    ylab(TeX("$\\textit{\\mu_F}/\\textit{F}$ ($yr^{-1}$)")) 
  figFish <- semilogx(figFish, xlim=xlim, ylim=c(0,1.1),
                      bXaxis = FALSE, bYaxis=bYaxis)+
    scale_y_continuous(breaks=seq(0,1,0.5))
  #
  # Spectra
  #
  figS <- ggplot()
  for (i in 1:3) 
    figS <- figS +
      geom_line(data=subset(fish,sim==i), aes(x=w, y=Nc), size=i/2) +
      #geom_line(data=subset(fish,sim==i), aes(x=res0$wR, y=NR), size=thick, color=grey(0.5)) +
      ylab("Rel. change") 
  # resource:
  #figS <- figS + geom_line(data=res0$resource, aes(x=wR, y=1+0*wR), size=thick, color=grey(0.5))
  figS <- loglog(figS, ylim=c(0.05, 20), xlim=xlim, bXaxis = FALSE, bYaxis=bYaxis) +
    geom_hline(yintercept = 1, linetype="dashed") 
  #
  # Feeding level:
  #
  figF <- ggplot()
  for (i in 1:3) 
    figF <- figF +
    geom_line(data=subset(fish,sim==i), aes(x=w, y=f), size=i/2)  +
    geom_line(data=res0$fish, aes(x=w, y=f), linetype="dashed") +
    geom_hline(yintercept = param$f0, linetype="dotted", size=thin) +
    ylab("Feeding level") + xlab("Weight (g)")
  figF <- semilogx(figF, xlim=xlim, ylim=c(0,1.1),
                   bYaxis=bYaxis)+
    scale_y_continuous(breaks=seq(0,1,0.2), limits=c(0,1)) 
  
  #
  # Predation mortality:
  #
  figM <- ggplot()
  for (i in 1:3) 
    figM <- figM +
    geom_line(data=subset(fish,sim==i), aes(x=w, y=muP), size=i/2)
  figM <- figM +
    geom_line(data=res0$fish, aes(x=w, y=muP), linetype="dashed") +
    geom_line(data=res0$fish, aes(x=w, y=param$A*param$aTheo*w^(param$n-1)), linetype="dotted", size=thin) +
    ylab(TeX("Mortality ($yr^{-1}$)")) +
    xlab("Weight (g)")
  
  figM <- semilogx(figM, xlim=xlim, ylim=c(0.05,4), bXaxis = FALSE, bYaxis=bYaxis)
  
  
  
  fig <- plot_grid(figFish,figS,figM,figF, ncol=1, align="v", 
                   rel_heights = c(0.5,1,1,1.25))
  #ggsave("Comparison.pdf", fig)
  fig
}

plotCascades <- function(nSpecies=27) {
  #
  # Make base run:
  #
  param <- paramTraitbasedmodel(W=10^seq(log10(4),5,length.out = nSpecies))
  res0 <- runspectrum(param)
  F <- c(0.1, 0.3, 0.7)
  #
  # Make runs with three different fishing mortalities, while fishing on big fish:
  #
  res0$fish$sim <- 0
  fish <- res0$fish[FALSE,]
  
  param$funcFishing <- function(w, p) {
    muF <- matrix(0, ncol = length(w), nrow = length(p$W))
    for (i in 1:length(p$W))
      muF[i,] <- p$F[i] * psi(w/1000, u=3)
    muF
  }
  
  for (i in 1:length(F)) {
    param$F <- 0*param$F + F[i]
    res <- runspectrum(param, res0)
    res$fish$Nc <- res$fish$Nc / res0$fish$Nc
    #res$fish$NR <- res$resource$NR/res0$resource$NR
    res$fish$sim <- i
    fish <- rbind(fish, res$fish)
  }
  panel1 <- panelComparison(param,res0,fish, bYaxis=TRUE)
  panel1
  #
  # Make runs with three different fishing mortalities, while fishing on all fish:
  #
  param$funcFishing <- fishingTrawl
  fish <- res0$fish[FALSE,]
  
  for (i in 1:length(F)) {
    param$F <- 0*param$F + F[i]
    res <- runspectrum(param, res0)
    res$fish$Nc <- res$fish$Nc / res0$fish$Nc
    res$fish$sim <- i
    fish <- rbind(fish, res$fish)
  }
  panel2 <- panelComparison(param,res0,fish, bYaxis=FALSE)
  panel2
  
  fig <- plot_grid(panel1, panel2, ncol=2, align="h", 
                   rel_widths = c(1.3,1))
  ggsave("Chapter9/Cascades.pdf", fig, width=doublewidth, height=2.5*height)
  fig
}


plotForagefishing <- function() {
  maxSizeFF <- 150
  minSizeFF <- 5
  #
  # Base case: fishing the entire community:
  #
  param <- paramTraitbasedmodel()
  param$F <- 0*param$F + 0
  res0 <- runspectrum(param)
  #
  # Define fleets: forage fish Winf < 100 g; large consumers: Winf > 2 kg
  #
  ixForage <- param$W>=minSizeFF & param$W < maxSizeFF  # Note: also update these limits further down
  ixPelagic <- param$W >= maxSizeFF & param$W < 5000
  ixBig <- param$W >=5000
  #
  # Calc yield, SSB, R0 for a combination of forage and consumer fisheries
  #
  if (bRecalcExpensiveFunctions) {
    n <- 50
    F <- seq(0,2, length.out = n)
    R <- list()
    ix <- 1
    for (i in 1:n) { # Forage fish
      res <- res0
      for (j in 1:n) {  # Consumer fish
        param$F[ixForage] <- F[i]
        param$F[ixPelagic] <- F[j]
        param$F[ixBig] <- F[j]
        res <- runspectrum(param, res)
        R[[ix]] <- res
        ix <- ix + 1
      }
    }
    # Calc yield:
    Y <- rep(0, param$nSpecies)
    ix <- 1
    dat <- data.frame(FF=double(), FC=double(), 
                      YF=double(), YP=double(), YB=double(),
                      R0F=double(), R0P=double(), R0B=double())
    for (i in 1:n) { # Forage fish
      for (j in 1:n) {  # Consumer fish
        res <- R[[ix]]
        muF <- res$muF
        for (k in 1:param$nSpecies)
          Y[k] <- sum(res$N[res$nSave,k,]*res$fish$w*res$fish$dw*muF[k,])
        dat <- rbind(dat, data.frame(FF=F[i], FC=F[j], 
                                     YF=sum(Y[ixForage]), R0F=min(subset(res$R, param$W[species]>=minSizeFF & param$W[species]<maxSizeFF & time==param$Tend)$R0), 
                                     YP=sum(Y[ixPelagic]), R0P=min(subset(res$R, param$W[species]>=maxSizeFF &  param$W[species]<5000 & time==param$Tend)$R0), 
                                     YB=sum(Y[ixBig]), R0B=min(subset(res$R, param$W[species]>=5000 & time==param$Tend)$R0)))
        ix <- ix + 1
      }
    }
    save(dat, file="Data/Foragefishing.RData")
  }
  else
    load("Data/Foragefishing.RData")
  
  figFF <- ggplot(data=dat) +
    geom_raster(aes(x=FF, y=FC, fill=YF/max(YF)), interpolate=TRUE) +
    guides(fill=FALSE) +
    geom_contour(aes(x=FF, y=FC, z=YF/max(YF)), colour="black") +
    scale_fill_gradient(low="white", high=grey(0.3)) +
    #geom_contour(aes(x=FF,y=FC,z=R0F),breaks=1.01,colour="white",size=3) +
    geom_raster(aes(x=FF, y=FC, alpha=1-sign(R0F-1.01), fill=0), interpolate=FALSE) +
    guides(alpha=FALSE) +
    xlab(" ") +
    ylab(TeX("Consumer fishery ($yr^{-1}$)")) +
    ggtitle("Forage fish") +
    scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) +
    annotate("text", x=2, y=2, label="a", hjust=1, vjust=1)
  figFF <- mytheme(figFF)
  
  figFP <- ggplot(data=dat) +
    geom_raster(aes(x=FF, y=FC, fill=YP/max(YP)), interpolate=TRUE) +
    guides(fill=FALSE) +
    geom_contour(aes(x=FF, y=FC, z=YP/max(YP)), colour="black") +
    scale_fill_gradient(low="white", high=grey(0.3)) +
    #geom_contour(aes(x=FF,y=FC,z=R0P),breaks=1.1,colour="white",size=3) +
    geom_raster(aes(x=FF, y=FC, alpha=1-sign(R0P-1.1), fill=0), interpolate=FALSE) +
    guides(alpha=FALSE) +
    xlab(TeX("Forage fishery ($yr^{-1}$)")) +
    ylab("") +
    ggtitle("Pelagics") +
    scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) +
    annotate("text", x=0, y=2, label="b", hjust=-0.5, vjust=1)
  figFP <- mytheme(figFP, bYaxis = FALSE)
  
  figFB <- ggplot(data=dat) +
    geom_raster(aes(x=FF, y=FC, fill=YB/max(YB)), interpolate=TRUE) +
    guides(fill=FALSE) +
    geom_contour(aes(x=FF, y=FC, z=YB/max(YB)), colour="black") +
    scale_fill_gradient(low="white", high=grey(0.3)) +
    #geom_contour(aes(x=FF,y=FC,z=R0B),breaks=1.1,colour="white",size=3) +
    geom_raster(aes(x=FF, y=FC, alpha=1-sign(R0B-1.1), fill=0), interpolate=FALSE) +
    guides(alpha=FALSE) +
    xlab("") +
    ylab("") +
    ggtitle("Large fish")+
    scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) +
    annotate("text", x=2, y=2, label="c", hjust=1, vjust=1)
  figFB <- mytheme(figFB, bYaxis = FALSE)
  
  
  grid <- plot_grid(figFF, figFP, figFB, ncol=3, align="h", 
                    rel_widths = c(1.4,1,1))
  ggsave("Chapter9/Foragefishing.pdf", grid,
         width=doublewidth, height=height)
}


plotYieldvsF <- function() {
  param <- paramTraitbasedmodel(nSpecies=30)
  F <- seq(0,3, length.out = 30)
  R <- list()  # Structure to hold the results of simulations
  dat <- data.frame(F=NULL, B=NULL, Y=NULL, fracCrash=NULL, fracLowSSB=NULL)
  dat2 <- data.frame(W=NULL, F=NULL, B=NULL, Y=NULL, R=NULL)
  
  if (bRecalcExpensiveFunctions) {
    for (i in 1:length(F)) {
      param$F = 0*param$F + F[i]
      res <- runspectrum(param)
      R[[i]] <- res
      Rend <- subset(res$R, time==param$Tend)
      dat <- rbind(dat, data.frame(
        F=F[i], 
        B=sum( Rend$SSB ), 
        Y=sum( Rend$Y ), 
        fracCrash=sum(Rend$R < 0.5)/param$nSpecies,  # F > Flim
        fracLowSSB = sum( (Rend$SSB/subset(R[[1]]$R, time==param$Tend)$SSB)<0.2)/param$nSpecies)) # SSB < 0.2 SSB0
      dat2 <- rbind(dat2, data.frame(
        W=param$W,
        F=F[i],
        B=Rend$SSB/subset(R[[1]]$R, time==param$Tend)$SSB,
        Y=Rend$Y,
        R = Rend$R))
    }
    save(dat,dat2,file="Data/YieldvsF.RData")
  } else
    load(file="Data/YieldvsF.Rdata")

  fig <- ggplot(data=dat) +
    geom_line(aes(x=F, y=B/max(B)), size=thick, color=grey(0.5)) +
    geom_line(aes(x=F, y=Y/max(Y)), size=thick) +
    geom_line(aes(x=F, y=fracCrash), linetype="dashed", size=thick) +
    geom_line(aes(x=F, y=fracLowSSB), linetype="dashed") +
    xlab(TeX("Fishing mortality ($yr^{-1}$)")) +
    ylab("Scaled quantities") +
    annotate("text", x=.3, y=1, label="a", hjust=1, vjust=1)
  fig <- mytheme(fig)
  
  fig2 <- ggplot(data=dat2) +
    geom_raster(aes(x=F, y=W, fill=(Y/max(Y))), interpolate = TRUE) +
    scale_fill_gradient(low="white", high=grey(0.3)) +
    guides(fill=FALSE) +
    geom_contour(aes(x=F, y=W, z=Y/max(Y)), colour="black") +
    geom_contour(aes(x=F, y=W, z=R),breaks=0.5,colour="black",size=3) +
    geom_contour(aes(x=F, y=W, z=B), breaks=0.2, colour="black")+ 
    xlab(TeX("Fishing mortality ($yr^{-1}$)")) +
    ylab("Asymptotic weight (g)") 
  fig2 <- semilogy(fig2) 
  
  fig3 <- ggplot(data=dat2) +
    geom_raster(aes(x=F, y=W, fill=B), interpolate = TRUE) +
    scale_fill_gradient(low="white", high=grey(0.3)) +
    guides(fill=FALSE) +
    geom_contour(aes(x=F, y=W, z=B), colour="black", size=thick, breaks=1) +
    geom_contour(aes(x=F, y=W, z=B), colour="black", breaks=c(0.2,0.5)) +
    geom_contour(aes(x=F, y=W, z=B), colour="white", breaks=c(2,5,10)) +
    #geom_contour(aes(x=F, y=W, z=R),breaks=0.5,colour="black",size=3) +
    #geom_contour(aes(x=F, y=W, z=B), breaks=0.2, colour="black")+ 
    xlab(TeX("Fishing mortality ($yr^{-1}$)")) +
    ylab("Asymptotic weight (g)")+
    annotate("text", x=.5, y=1e5, label="b", hjust=1, vjust=1)
  fig3 <- semilogy(fig3)
  
  
  grid <- plot_grid(fig,fig3, ncol=2, align="h")
  ggsave("Chapter9/yield.pdf", grid , width=doublewidth, height=height)
}



plotIllustrateMSY <- function() {
  #
  # Define fleets: forage fish Winf < 100 g; large consumers: Winf > 2 kg
  #
  ixForage <- param$W < 150  # Note: also update these limits further down
  ixPelagic <- param$W >= 150 & param$W < 5000
  ixBig <- param$W >=5000
  
  # Base run:
  param <- paramTraitbasedmodel()
  res0 <- runspectrum(param)
  Fmsy0 <- calcRefpointsTraitbasedmodel(param, res0)
  
  # Fishing on large fish:
  param$F[ixForage] <- 0
  param$F[ixPelagic] <- 0.3
  param$F[ixBig] <- 0.3
  res1 <- runspectrum(param, res0)
  Fmsy1 <- calcRefpointsTraitbasedmodel(param, res1)
  
  # Fishing on everything:
  param$F[ixForage] <- 0.3
  param$F[ixPelagic] <- 0.3
  param$F[ixBig] <- 0.3
  res2 <- runspectrum(param, res0)
  Fmsy2 <- calcRefpointsTraitbasedmodel(param, res2)
  
  # Fish on everything, but lower big fish to 0.2:
  param$F[ixForage] <- 0.3
  param$F[ixPelagic] <- 0.3
  param$F[ixBig] <- 0.1
  res3 <- runspectrum(param, res0)
  Fmsy3 <- calcRefpointsTraitbasedmodel(param, res3)
  
  # Fish on everything, but increase forage fishing to 0.5:
  param$F[ixForage] <- 0.8
  param$F[ixPelagic] <- 0.3
  param$F[ixBig] <- 0.3
  res4 <- runspectrum(param, res0)
  Fmsy4 <- calcRefpointsTraitbasedmodel(param, res4)
  
  
  
  # Assemble data frame
  Fmsy <- Fmsy2
  Fmsy$case <- 1
  #Fmsy <- rbind(Fmsy, cbind(Fmsy1, data.frame(case=2)))
  #Fmsy <- rbind(Fmsy, cbind(Fmsy2, data.frame(case=3)))
  Fmsy <- rbind(Fmsy, cbind(Fmsy3, data.frame(case=4)))
  Fmsy <- rbind(Fmsy, cbind(Fmsy4, data.frame(case=5)))
  Fmsy$case <- as.factor(Fmsy$case)
  
  fig <- ggplot(data=Fmsy, aes(x=W, y=Fmsy, linetype=case)) +
    geom_line() + geom_point()
  fig <- semilogx(fig) +
    xlab("Asymptotic mass (g)") +
    ylab(TeX("$F_{msy}$ $(yr^{-1})$"))
  fig
}
#
# Plot reference points vs overall F:
#
plotRefvsF <- function() {
  nSpecies <- 18
  param <- paramTraitbasedmodel(nSpecies = nSpecies)
  F <- seq(0,1, length.out = 3)
  F <- c(0, 0.3, 0.7)
  R <- list()  # Structure to hold the results of simulations
  refs <- data.frame(F=NULL, W=NULL, Fmsy=NULL)
  iSpecies <- 1:nSpecies # Species to calc ref points for
  
  if (bRecalcExpensiveFunctions) {
    #
    # Ref points from commuinty model:
    #
    for (i in 1:length(F)) {
      print(F[i])
      param$F = 0*param$F + F[i]
      res0 <- runspectrum(param)
      R[[i]] <- res0
      refpoints <- calcRefpointsTraitbasedmodel(param,res0,iSpecies)
      refs <- rbind(refs, data.frame(
        F = rep(F[i], length(iSpecies)),
        W = refpoints$W,
        Fmsy = refpoints$Fmsy))
    }
    #
    # Ref points from single species model for reference:
    #
    MSYsingle <- data.frame(W=param$W, Fmsyminus=10+0*param$W, Fmsyplus = 0*param$W)
    paramSingle <- baseparameters()
    a0 = paramSingle$a
    a_range = a0 + seq(-0.15, 0.15, length.out = 10)
    for (i in 1:length(param$W)) {
      paramSingle$W <- param$W[i]
      paramSingle$a = a0
      tmp <- calcRefpoints(paramSingle)
      MSYsingle$W[i] <- param$W[i]
      MSYsingle$Fmsy[i] <- tmp$Fmsy
      # Calc Fmsy for a range of "a"'s and find the max and min values:
      for (a in a_range) {
        paramSingle$a = a
        tmp <- calcRefpoints(paramSingle)
        MSYsingle$Fmsyminus[i] = min(MSYsingle$Fmsyminus[i], tmp$Fmsy)
        MSYsingle$Fmsyplus[i]  = max(MSYsingle$Fmsyplus[i], tmp$Fmsy)
      }
    }
    save(MSYsingle, refs, file="Data/RefsvsF.RData")
  } 
  else {
    load("Data/RefsvsF.Rdata")
  }

  defaultplot()
  par(mar=par()$mar + c(0,0,0,6)) # space for legend
  semilogxpanel(xlim=c(4,1e5), ylim=c(0,0.8),
                xlab=("Asymptotic weight (g)"),
                ylab=("$\\textit{F}_{msy}$ (yr$^{-1}$)")
  )
  ribbon(MSYsingle$W, MSYsingle$Fmsyplus, MSYsingle$Fmsyminus, col=lightgrey)
  lines(MSYsingle$W, MSYsingle$Fmsy, lty=dashed, lwd=3, col="white")  
  for (i in 1:length(F)) {
    ix = refs$F==F[i]
    lines(refs$W[ix], refs$Fmsy[ix], lwd=i)
  }

  legend("right", bty="n", inset=c(-0.4,0), xpd=TRUE,
         legend=c("0",TeX("$0.3$ yr$^{-1}$"),TeX("$0.7$ yr$^{-1}$")),
         lwd=seq(1,3),
         title=TeX("$\\textit{F}$"))
  
}


plotCommunityMSY <- function(param=paramTraitbasedmodel(), n=20) {
  #
  # Define fleets: forage fish Winf < 100 g; large consumers: Winf > 2 kg
  #
  ixForage <- param$W < 150  # Note: also update these limits further down
  ixPelagic <- param$W >= 150 & param$W < 5000
  ixBig <- param$W >=5000
  #
  # Reference run:
  # 
  res00 <- runspectrum(param)
  bConstrainSSB <- FALSE # Flag telling algorithms whether SSB should be constrained above 0.2 SSB0
  
  F <- 10^seq(-1,log10(8),length.out=n)
  p <- param
  #F <- 10^seq(0,1,length.out=4)
  
  # Function to optimize. Returns the total yield
  funcYield <- function(fracF, F, res0) {
    p$F[ixForage] <- fracF[1]
    p$F[ixPelagic] <- fracF[2]
    p$F[ixBig] <- 1-fracF[1]-fracF[2]
    p$F <- F*p$F
    
    res <- runspectrum(p, res0)
    yield <- sum( subset(res$R, time==p$Tend)$Y )
    # Penalize if some species are below 0.2 SSB0:
    SSBratio <- subset(res$R, time==res$nSave)$SSB / subset(res00$R, time==res00$nSave)$SSB
    nCrashed <- sum(SSBratio <0.2) 
    if (bConstrainSSB & nCrashed>0)
      yield <- sum( log(SSBratio[SSBratio<0.2]/0.2) ) #-sum( 0.2-SSBratio[SSBratio<0.2] )
    #cat("frac ", fracF, " yield:", yield ," # crashed: ",nCrashed,"\n")
    return( -yield )
  }
  
  funcRent <- function(fracF, F, res0) {
    p$F[ixForage] <- fracF[1]
    p$F[ixPelagic] <- fracF[2]
    p$F[ixBig] <- 1-fracF[1]-fracF[2]
    p$F <- F*p$F
    
    res <- runspectrum(p, res0)
    price <- res$fish$w^0.41
    cost <- .5e-5 * p$F * p$W^0.4 #.5e-5 * p$F * p$W^0.4
    rent <- 0
    for (j in 1:p$nSpecies)
      rent = rent + 
      sum( res$fish$w*res$N[res$nSave,j,]*price*res$fish$dw*res$muF[j,]-cost[j])
    #cat("F ",F,"frac ", fracF, " rent:", rent,"\n")
    return( -rent )
  }
  
  optOverF <- function(funcOpt) {
    opt <- list()
    optF = data.frame(F=NULL, FF=NULL, FP=NULL, FB=NULL, Y=NULL, FFerr=NULL, FPerr=NULL, FBerr=NULL)
    f0 <- c(0.2, 0.2)
    for (i in 1:length(F)) {
      p$F[ixForage] <- f0[1]
      p$F[ixPelagic] <- f0[2]
      p$F[ixBig] <- 1-f0[1]-f0[2]
      p$F <- F[i]*p$F
      res0 <- runspectrum(p)
      ixF <- i
      opt[[i]] <- constrOptim(
        theta=f0, f=funcOpt, F=F[i], res0=res0, grad=NULL, 
        ui=matrix(c(-1,1,0,-1,0,1), ncol=2),ci=c(-1,0,0), 
        outer.iterations = 50, outer.eps=0.1)
      f0 <- opt[[i]]$par
      hess <- hessian(funcOpt, f0, F=F[i], res=res0)
      optF <- rbind(optF, data.frame(
        F = F[i], 
        FF = opt[[i]]$par[1], 
        FP = opt[[i]]$par[2], 
        FB = 1-opt[[i]]$par[1]-opt[[i]]$par[2],
        FFerr = sqrt(0.2*abs(opt[[i]]$value)/hess[1,1]), 
        FPerr = sqrt(0.2*abs(opt[[i]]$value)/hess[2,2]),
        FBerr = 0.5*( sqrt(0.2*abs(opt[[i]]$value)/hess[1,1]) + sqrt(0.2*abs(opt[[i]]$value)/hess[2,2])), 
        Y = -opt[[i]]$value ))
      
      cat(F[i], ": ",opt[[i]]$par,". Yield: ", -opt[[i]]$value, "\n-------------------\n")
    }
    return(optF)
  }
  
  if (bRecalcExpensiveFunctions) {
    #
    # Optimize yield
    #
    bConstrainSSB <- FALSE
    optYield <- optOverF(funcYield)
    #
    # Optimize rent
    #
    bConstrainSSB <- FALSE
    optRent <- optOverF(funcRent)
    #
    # Optimize yield constrained:
    #
    bConstrainSSB <- TRUE
    optYieldConstrained <- optOverF(funcYield)
    
    # Save calculations
    save(optYield, optYieldConstrained, optRent, file="Data/optimization.Rdata")
  }
  else
    load("Data/optimization.Rdata")
  
  
  # fig1 <- ggplot(dat=optYield) +
  #   #geom_ribbon(aes(x=F/3, ymin=FF-FFerr, ymax=FF+FFerr),fill=grey(0.8), alpha=0.5) +
  #   #geom_ribbon(aes(x=F/3, ymin=FP-FPerr, ymax=FP+FPerr),fill=grey(0.8), alpha=0.5) +
  #   #geom_ribbon(aes(x=F/3, ymin=FB-FBerr, ymax=FB+FBerr),fill=grey(0.8), alpha=0.5) +
  #   geom_line(aes(x=F/3, y=FF), size=thin) +
  #   geom_line(aes(x=F/3, y=FP)) +
  #   geom_line(aes(x=F/3, y=FB), size=thick) +
  #   geom_line(aes(x=F/3, y=Y/max(Y)), size=thick, colour=grey(0.8))
  # fig1 <- semilogx(fig1, ylim=c(0,1)) +
  #   xlab(TeX("Average F ($yr^{-1}$)")) + 
  #   ylab("Fraction") 
  # fig1
  # 
  # fig2 <- ggplot(dat=optRent) +
  #   #geom_ribbon(aes(x=F/3, ymin=FF-FFerr, ymax=FF+FFerr),fill=grey(0.8), alpha=0.5) +
  #   #geom_ribbon(aes(x=F/3, ymin=FP-FPerr, ymax=FP+FPerr),fill=grey(0.8), alpha=0.5) +
  #   #geom_ribbon(aes(x=F/3, ymin=FB-FBerr, ymax=FB+FBerr),fill=grey(0.8), alpha=0.5) +
  #   geom_line(aes(x=F/3, y=FF), size=thin) +
  #   geom_line(aes(x=F/3, y=FP)) +
  #   geom_line(aes(x=F/3, y=FB), size=thick) +
  #   geom_line(aes(x=F/3, y=Y/max(Y)), size=thick, colour="grey")
  # fig2 <- semilogx(fig2, ylim=c(0,1)) +
  #   xlab(TeX("Average F ($yr^{-1}$)")) + 
  #   ylab("Fraction")
  # fig2
  # 
  # fig3 <- ggplot(dat=optYieldConstrained) +
  #   #geom_ribbon(aes(x=F/3, ymin=FF-FFerr, ymax=FF+FFerr),fill=grey(0.8), alpha=0.5) +
  #   #geom_ribbon(aes(x=F/3, ymin=FP-FPerr, ymax=FP+FPerr),fill=grey(0.8), alpha=0.5) +
  #   #geom_ribbon(aes(x=F/3, ymin=FB-FBerr, ymax=FB+FBerr),fill=grey(0.8), alpha=0.5) +
  #   geom_line(aes(x=F/3, y=FF), size=thin) +
  #   geom_line(aes(x=F/3, y=FP)) +
  #   geom_line(aes(x=F/3, y=FB), size=thick) +
  #   geom_line(aes(x=F/3, y=Y/max(optYield$Y)), size=thick, colour="grey")
  # fig3 <- semilogx(fig3, ylim=c(0,1)) +
  #   xlab(TeX("Average F ($yr^{-1}$)")) + 
  #   ylab("Fraction")
  # fig3
  # 
  # grid <- plot_grid(fig1, fig2, fig3, ncol=3, align="h")
  # ggsave("Chapter9/CommunityMSY.pdf", width=doublewidth, height=height)
  # grid
  
  
  
  defaultplot()
  par(mar=par()$mar + c(0,0,0,10)) # space for legend
  defaultpanel(xlim=c(0,2.67), ylim=c(0,1),
               xlab="Average $\\textit{F}$ ($yr^{-1}$)", ylab="Fraction")
  lines(optYield$F/3, optYield$FF, col=stdgrey, lwd=1)
  lines(optYield$F/3, optYield$FP, col=stdgrey, lwd=2)
  lines(optYield$F/3, optYield$FB, col=stdgrey, lwd=3)
  lines(optYield$F/3, optYield$Y/max(optYield$Y), lwd=3)

  legend("right", bty="n", inset=c(-0.6,0), xpd=TRUE,
         legend=c(TeX("5 g < $\\textit{W}_{\\infty}$ < 150 g"),
                  TeX("150 g $\\leq$ $\\textit{W}_{\\infty}$ < 5 kg"),
                  TeX("$\\textit{W}_{\\infty}$ $\\geq$ 5 kg")),  
         col=stdgrey, lwd=seq(1,3),
         title="Fleet")
}


panelForageFishing <- function() {
  #
  # Make base run:
  #
  param <- paramTraitbasedmodel(W=10^seq(log10(4),5,length.out = 9))
  res0 <- runspectrum(param)
  F <- c(0, 0.25, 0.7)
  #
  # Make runs with three different fishing mortalities, while fishing on big fish:
  #
  res0$fish$sim <- 0
  listRes <- vector("list", length(F))
  
  for (i in 1:length(F)) {
    param$F <- 0*param$F
    param$F[param$W<250] <- 0.6
    param$F[param$W>=250] <- F[i]
    listRes[[i]] <- runspectrum(param, res0)
  }
  
  panel1 <- panelComparison(param,res0,fish)
  panel1
}




growthDynamic <- function(p, w)
{
  p$epsA*p$h*(p$f0-p$fc)*w^p$n*( 1 - psi(w/(p$etaM*p$W)) * (w/p$W)^(1-p$n)) 
}






#
# Check the MSY for various level of FF and consumer fishing
#
testForageFishingMSY <- function(param) {
  param <- paramTraitbasedmodel()
  
  param$F <- c(0,0,0,0.5,0.5,0.5,0.5,0.5,0.5)
  res0 <- runspectrum(param)
  ref0 <- calcRefpointsTraitbasedmodel(param,res0)
  
  param$F <- c(0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5)
  res1 <- runspectrum(param)
  ref1 <- calcRefpointsTraitbasedmodel(param,res1)
  
  param$F <- c(0.5,0.5,0.5,0.5,0.5,0.5,0.,0.,0.)
  res2 <- runspectrum(param)
  ref2 <- calcRefpointsTraitbasedmodel(param,res2)
}


plotDoubleYield <- function() {
  F <- seq(0,2,length.out=20);
  param <- paramTraitbasedmodel(nSpecies=30)
  
  R <- list()
  dat <- data.frame(F=NULL, B=NULL, Y=NULL, fracCrash=NULL, fracLowSSB=NULL)
  dat2 <- data.frame(W=NULL, F=NULL, B=NULL, Y=NULL, R=NULL)
  for (i in 1:length(F)) {
    param$F <- rep(F[i], param$nSpecies)
    res <- runspectrum(param)
    R[[i]] <- res
    
    Rend <- subset(res$R, time==param$Tend)
    dat <- rbind(dat, data.frame(
      F=F[i], 
      B=sum( Rend$SSB ), 
      Y=sum( Rend$Y ), 
      fracCrash=sum(Rend$R < 0.5)/param$nSpecies,  # F > Flim
      fracLowSSB = sum( (Rend$SSB/subset(R[[1]]$R, time==param$Tend)$SSB)<0.2)/param$nSpecies)) # SSB < 0.2 SSB0
    dat2 <- rbind(dat2, data.frame(
      W=param$W,
      F=F[i],
      B=Rend$SSB/subset(R[[1]]$R, time==param$Tend)$SSB,
      Y=Rend$Y,
      R = Rend$R))
  }
  
  pdf("doubleyield1.pdf", width=doublewidth, height=0.5*height)
  defaultplot()
  defaultpanel(xlim=F[2:length(F)], ylim=c(0,1), xlab="Fishing mortality", ylab="Yield")
  lines(dat$F, dat$Y/max(dat$Y), lwd=2)
  dev.off()
  
  
  #W = unique(dat2$W)
  #ix = c(6,8,11)
  #B0 = subset(dat2,F==0)$B
  #loglogpanel(xlim=W, ylim=c(0.01,10),
  #            xlab="Asymptotic weight (g)", ylab="Relative biomass")
  #for (i in 1:length(ix)) {
  #  lines(W, dat2$B[dat2$F==F[ix[i]]], lwd=i)
  #}
  #hline(1)
  
  pdf("doubleyield2.pdf", width=doublewidth, height=0.75*doublewidth)
  defaultplot()
  w <- res$w
  wR <- R[[1]]$resource$wR
  ix = c(3,6,20)
  loglogpanel(xlim=c(0.1, max(w)), ylim=c(1e-8,0.05),
              xlab="Weight (g)", ylab="Biomass density")
  lines(w, R[[1]]$fish$Nc*w, lwd=3, col="red")
  lines(wR, R[[1]]$resource$NR*wR, lwd=2, col="green")
  for (i in 1:length(ix)) {
    lines(w, R[[ix[i]]]$fish$Nc*w, lwd=3, col="orange")
    lines(wR, R[[ix[i]]]$resource$NR*wR, lwd=2, col="green")
  }
  dev.off()
}


plotAllChapter9 <- function() {
  pdfplot("Chapter9/Daan.pdf", plotDaan, width=1.3*singlewidth, height=height)
  plotCascades()
  plotForagefishing()
  plotYieldvsF()
  pdfplot("Chapter9/CommunityMSY.pdf", plotCommunityMSY, width=doublewidth, height=height)
  pdfplot("Chapter9/RefvsF.pdf", plotRefvsF, width=1.3*singlewidth, height=height)
}





