source("Rcode/basetools.R")
source("Rcode/basefunctions.R")
source("Rcode/baseparameters.R")
source("Rcode/community.R")

getOlssonData <- function() {
  data <- read.csv("Data/fecundpruned.csv",header=TRUE,sep=",")
  data$Winf <- data$aW * data$Linf^data$bW
  data$Wm  <- data$aW * data$Lm^data$bW
  data$Phylum <- data$phyl

  # Remove pike, which seems wrong (also is not in the table in her thesis)
  data <- data[which(data$species!="Esox lucius"),]
  
  ## calc W0
  ix = is.na(data$W0)
  data$W0[ix] <- data$aW[ix] * data$L0[ix]^data$bW[ix]
  data <- data[!is.na(data$W0),]
  
  #
  # Fit R, but this time correct with A
  #
  # Calc weight-specific R ; assume that it is at size at maturation
  data$cycle[is.na(data$cycle)]=1
  data$R <-  data$litter/data$cycle*data$W0 / (data$Wm)
  # no need to correct with temperature when we correct with A
  
  # Calc A
  data$A = calcA(data$K, data$Linf)
  data$k = data$A * data$Winf^(-0.25)
  
  data  
}

plotComparison <- function() {
  defaultplot(mfcol=c(1,4))
  #defaultplothorizontal(nPanels = 4)
  data <- getOlssonData()
  # -----
  # A
  # ----
  #
  # Data from Kooijman and Gislason:
  #
  data2 = read.csv("Data/KooijmanGislason.csv", header=FALSE, col.names=c("K","Linf","T"))
  data2$Phylum = "t" # all teleost data
  data2$A <- calcA(data2$K, data2$Linf)
  #
  # Combine data
  # 
  Q10 <- 1.83
  dat <- data.frame(A=data$A * Q10^((15-data$Temp)/10), W=data$Winf, Phylum=data$Phylum)
  dat2 <- data.frame(A=data2$A * Q10^((15-data2$T)/10), W=0.01*data2$Linf^3, Phylum=data2$Phylum)
  dat = rbind(dat, dat2)
  #
  # plot
  #
  boxplot(log10(A)~Phylum, data=dat, notch=TRUE, 
          ylab=TeX("Growth constant $log_{10}(\\textit{A}$) ($g^{0.25}/yr$)"),
          names=c("Teleost","Elas."), cex.lab=0.9)
  makepanellabel()
  # -----
  #  etaM
  # -----
  boxplot(log10(Wm/Winf)~Phylum, data=data, notch=TRUE, 
          ylab=TeX("Maturation ratio $log_{10}(\\textit{\\eta_m})"),
          names=c("Teleost","Elas."), cex.lab=0.9)
  makepanellabel()
  # -----
  #  R*w^0.25
  # -----
  boxplot(log10(R*Winf^0.25)~Phylum, data=data, notch=TRUE, 
          ylab=TeX("Repro. $log_{10}(\\textit{R}_{egg} W_{\\infty}$$^{0.25})$ $(yr^{-1})"),
          names=c("Teleost","Elas."), cex.lab=0.9)
  makepanellabel()
  # ----
  #  M/K
  # ----
  # Get fish data from Gislason for fish:
  A <- read.csv("Data/Gislason.csv",sep=",", header=FALSE)
  names(A) <- c("L","Linf","K","M","tau")
  A <- A[A$L > 0.5*A$Linf, ]  # Selection only adults
  p <- baseparameters()
  A$a = A$M/A$K*p$aMKconstant
  A <- A[A$a < 1.25,]  # remove outliers
  A$Phylum = "t"
  # Use Olsson for elasmobranchs
  ixElas <- which(data$Phylum=="c" & !is.na(data$M))
  adat <- rbind(data.frame(a=A$a, Phylum=A$Phylum),
                data.frame(a=0.22*(data$M/data$K)[ixElas], Phylum=data$Phylum[ixElas]))
  
  boxplot(a~Phylum, data=adat, ylab=TeX("Physiological mortality, $\\textit{a}$"),
          names=c("Teleost","Elas."), cex.lab=0.9)
  makepanellabel()
}

plotEggSize_OnlyOlssonData <- function() {
  data <- getOlssonData()
  ixTeleO <- which(data$Phylum=="t" & data$repr=="o")
  ixTeleV <- which(data$Phylum=="t" & data$repr=="v")
  ixElasO <- which(data$Phylum=="c" & data$repr=="o")
  ixElasV <- which(data$Phylum=="c" & data$repr=="v")
  
  defaultplot()
  loglogpanel(xlim=c(50,1e6), ylim=range(data$W0),
              xlab="Asymptotic weight (g)", ylab="Offspring weight (g)")
  # Points:
  points(data$Winf[ixTeleO], data$W0[ixTeleO], pch=16)
  points(data$Winf[ixTeleV], data$W0[ixTeleV], pch=15)
  points(data$Winf[ixElasO], data$W0[ixElasO], pch=16, col="grey")
  points(data$Winf[ixElasV], data$W0[ixElasV], pch=15, col="grey")
  
  # Fits:
  ixTele <- which(data$Phylum=="t")
  w0 <- mean(data$W0[ixTele])
  lines(x=range(data$Winf[ixTele]), y=w0*c(1,1), lwd=2)
  
  ixElas <- which(data$Phylum=="c")
  coef <- exp(mean(log(data$W0[ixElas]/data$Winf[ixElas])))
  lines(x=range(data$Winf[ixElas]), y=coef*range(data$Winf[ixElas]), lwd=2, col="grey")
  cat('Elasmobranch offspring mass ratio: ', coef,'\n')
  
}

plotEggSize <- function() {
  #
  # Get sharks from Olsson:
  #
  data <- getOlssonData()
  ixTeleO <- which(data$Phylum=="t" & data$repr=="o")
  ixTeleV <- which(data$Phylum=="t" & data$repr=="v")
  ixElasO <- which(data$Phylum=="c" & data$repr=="o")
  ixElasV <- which(data$Phylum=="c" & data$repr=="v")
  dat <- data.frame(Winf=data$Winf[data$Phylum=="c"], w0=data$W0[data$Phylum=="c"], 
                    repr=data$repr[data$Phylum=="c"], Phylum="c")
  #
  # Get teleosts from Fishbase:
  #
  tab = read.table("Data/Teleost offspring size.dat", header = TRUE, fill=TRUE)
  tab = tab[tab$ClassNum==6,] # Only ray-finned species
  
  tab$w0 = pi*4/3*(tab$Eggdiammod/10/2)^3
  ixVivi = tab$RepGuild1=="bearers"
  tab$w0[ixVivi] = 0.01*(tab$LhMid[ixVivi]/10)^3 
  tab$repr = "o"
  tab$repr[ixVivi] = "v"
  tab = tab[tab$w0>0,]
  dat = rbind(dat, data.frame(Winf=0.01*tab$Loo^3, w0=tab$w0, repr=tab$repr, Phylum="t"))
  
  defaultplot()
  loglogpanel(xlim=c(2,2e6), ylim=range(dat$w0),
              xlab="Asymptotic weight (g)", ylab="Offspring weight (g)")
  # Points:
  ixTeleO <- which(dat$Phylum=="t" & dat$repr=="o")
  ixTeleV <- which(dat$Phylum=="t" & dat$repr=="v")
  ixElasO <- which(dat$Phylum=="c" & dat$repr=="o")
  ixElasV <- which(dat$Phylum=="c" & dat$repr=="v")
  points(dat$Winf[ixTeleO], dat$w0[ixTeleO], pch=16, col=gray(0.4, alpha=0.75))
  points(dat$Winf[ixTeleV], dat$w0[ixTeleV], pch=1, col=gray(0, alpha=0.75))
  points(dat$Winf[ixElasO], dat$w0[ixElasO], pch=17, col=gray(0, alpha=0.75))
  points(dat$Winf[ixElasV], dat$w0[ixElasV], pch=2, col=gray(0, alpha=0.75))
  
  # Fits:
  ixTele <- which(dat$Phylum=="t")
  w0 <- exp(mean(log(dat$w0[ixTele])))
  lines(x=range(dat$Winf[ixTele]), y=w0*c(1,1), lwd=2)
  cat('Teleost egg size: ',w0,' g\n')
  
  ixElas <- which(dat$Phylum=="c")
  coef <- exp(mean(log(dat$w0[ixElas]/dat$Winf[ixElas])))
  lines(x=range(dat$Winf[ixElas]), y=coef*range(dat$Winf[ixElas]), lwd=2)
  cat('Elasmobranch offspring mass ratio: ', coef,'\n')
  
  # Pictograms:
  addEpsPicture('Chapter10/Fish_silhouette_licensed_from_123rs.eps',
    x=0.87, y=0.47, width=0.2)
  addEpsPicture('Chapter10/Shark_silhouette_licensed_from_123rf.eps',
                x=0.54, y=0.86, width=0.2)
}

panelGrowthRates <- function() {
  p <- baseparameters()
  W <- 10^seq(1,6,length.out = 25)
  Wshark <- 10^seq(2,6,length.out=25)
  
  calcr <- function(W,e,w0) {
    p <- baseparameters()
    p$epsEgg <- 1
    p$epsR <- e
    p$w0 <- w0
    calc_rana1(W,p)
  }

  semilogxpanel(xlim=W,ylim=c(0,1),
                xlab="Asymptotic size (g)",
                ylab="Population growth rate (yr$^{-1}$)")
  lines(W, calc_rana1(W, p), lty=dashed, lwd=3) # Teleosts
  lines(Wshark, calcr(Wshark, 0.5*p$epsEgg, p$ElasmobranchEggMassRatio*Wshark), lwd=1)
  lines(Wshark, calcr(Wshark, 0.2*p$epsEgg, p$ElasmobranchEggMassRatio*Wshark), lwd=3)
  lines(Wshark, calcr(Wshark, p$epsR*p$epsEgg, p$ElasmobranchEggMassRatio*Wshark), lwd=1)
  lines(W, calcr(W, 0.5*p$epsEgg, p$ElasmobranchEggMassRatio*W), lwd=1, lty=dotted)
  lines(W, calcr(W, p$epsEgg*0.2, p$ElasmobranchEggMassRatio*W), lwd=1, lty=dotted)
  lines(W, calcr(W, p$epsEgg*p$epsR, p$ElasmobranchEggMassRatio*W), lwd=0.5, lty=dotted)
  #hline(0)
}
  
panelSharkRefpoints <- function() {
  p <- baseparameters()
  W <- 10^seq(1,6,length.out = 20)
  Wshark <- 10^seq(2,6,length.out=12)
  
  calcSharkRefpoint <- function(W) {
    pp <- p
    pp$w0 <- pp$ElasmobranchEggMassRatio*W
    pp$epsR <- 0.2
    pp$W <- W
    return(calcRefpoints(pp))
  }
  
  # Sharks
  FmsyShark <- 0*Wshark
  FcrashShark <- FmsyShark
  for (i in 1:length(Wshark)) {
    ref <- calcSharkRefpoint(Wshark[i])
    FmsyShark[i] <- ref$Fmsy
    FcrashShark[i] <- ref$Fcrash
  }
  
  # Data points from Zhou et al 2011:
  dat <- read.table("Data/Zhou et al (2011) supplementary.csv", 
                    sep=";", header=TRUE)
  dat <- dat[1:103,] # extract only elasmobranchs
  
  # Teleosts
  Fmsy <- 0*W
  Fcrash <- Fmsy
  for (i in 1:length(W)) {
    p$W <- W[i]
    ref <- calcRefpoints(p)
    Fmsy[i] <- ref$Fmsy
    Fcrash[i] <- ref$Fcrash
  }
  
  #loglogpanel(xlim=W, ylim=c(0.01, 1),
  semilogxpanel(xlim=W, ylim=c(0.01, 1),
                xlab="Asymptotic size (g)",
                ylab="Fishing mortality (yr^{-1})")
  points(weight(dat$L...), dat$FBRP, pch=dots)
  lines(W, Fmsy, lty=dashed, lwd=2)
  #lines(W, Fcrash, lty=dashed, lwd=2, col="grey")
  lines(Wshark, FmsyShark, lwd=3)
  lines(Wshark, FcrashShark, lwd=3, col="grey")
}

plotComparison2 <- function() {
  defaultplot(mfcol=c(1,2))
  panelGrowthRates()
  panelSharkRefpoints()
}


plotAllChapter10 <- function() {
  pdfplot(FUN=plotComparison, "Chapter10/comparison.pdf", width=doublewidth, height=height)
  pdfplot(FUN=plotEggSize, "Chapter10/eggsize.pdf", width=1.5*singlewidth, height=1.5*height)
  pdfplot(FUN=plotComparison2, "Chapter10/comparison2.pdf", width=doublewidth, height=height)
}
