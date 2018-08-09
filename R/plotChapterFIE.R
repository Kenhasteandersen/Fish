source('R/QuantitativeGenetics.R')

dir.create("ChapterFIE")

plotConover <- function() {
  defaultplot(mfcol=c(1,2))
  lty = c(1,1,3,3,2,2)
  #
  # First panel
  #
  data <- read.csv("Data/Conover2002.csv")
  
  defaultpanel(xlab='Age (days)', ylab="Weight (g)", 
               xlim = c(85, 190), ylim=c(0,8), xaxis=FALSE, label=TRUE)
  par(xaxp=c(90,190,5))
  axis(bottom)
 # plot( x=0,y=0, xlab='Age (days)', ylab="Weight (g)", 
#        xlim = c(85, 190), ylim=c(0,8),
#        xaxp=c(90,190,5))
  
  for (i in 0:5) {
    ix <- 1:6+6*i
    lines(x=data$day[ix], data$weight[ix], lty=lty[i+1], lwd=3)
  }
  
  legend(x=80, y=8, bty='n',
         legend=c('Harvest slow growers', 'Random harvest', 'Harvest fast growers'),
         lty=c(1,3,2), y.intersp = 0.8)
  #
  # Second panel
  #
  data <- read.csv("Data/Conover2009.csv")
  data$age <- round(data$age)
  # Calc. CUMULATIVE relative selection response. 
  # Assume that a change of 9 mm corresponds to a change of 1.5 g, 
  # and that the average weight is 3.7 g:
  data$Rcum <- data$deltaL/9*1.5/3.7
  # Calc relative selection response:
  data$R <- 0
  for (i in 1:10) {
    ix <- which(data$age==i)
    data$R[ix] <- (data$deltaL[ix]-data$deltaL[ix-1])/9*1.5/3.7
  }
  
  defaultpanel(xlab="Generation", ylab="Cumulated sel. response",
               xlim=c(0,10), ylim=c(-0.6,0.74), label=TRUE)
  #  defaultpanel(xlab="Generation", ylab="Rel. sel. response ($yr^{-1}$)",
  #       xlim=c(0,10), ylim=c(-.65,.65), label=TRUE)
  lines(x=c(0,10),y=c(0,0), lty=3, lwd=1)
  lines(x=c(5,5), y=c(-15,15), lty=3,lwd=1)
  for (i in 0:5) {
    ix <- 1:11+11*i
    lines(x=data$age[ix], data$Rcum[ix], lty=lty[i+1], lwd=3)
  }
  text(x=5,y=0.6,labels='Selection      ',adj=1)
  text(x=5,y=0.6, labels='   No selection', adj=0)

}

plotQG <- function() {
  defaultplot()
  lty = c(1,1,3,3,2,2)

  wmat0 = 1000
  p <- baseparamQG(wm=wmat0)

  # Plot distributions
  wmat <- seq(0, p$wm/p$etaM,length.out = 200)
  defaultpanel(xlim=c(0.3, 2000/wmat0), ylim=c(0,1.4),
               xlab='Rel. trait values, $\\theta /\\theta$', 
               ylab='Relative fitness, $\\textit{R}_0/\\textit{R}_0$', yaxis=FALSE)
  lines(x=wmat/wmat0, y=exp(-(wmat - 1.2*p$wm)^2 / (2*(p$cv*p$wm)^2)), 
        lwd=2, col=grey(0))#, col=gray(0.5))

  polygon(x=wmat/wmat0, y=exp(-(wmat - p$wm)^2 / (2*(p$cv*p$wm)^2)  ),
          border=NA, col=grey(0.4, alpha=0.75),yaxt='n')
  axis(side=2, at=c(0, 0.5, 1))
  
  
  # Plot R0
  R00 <- R(p)
  R0 <- seq(0,0,length.out=length(wmat))
  for (i in 1:length(R0)) {
    p$wm <- wmat[i]
    R0[i] <- R(p)
  }
  lines(x=wmat/wmat0, y=R0/R00, lty='dashed', lwd=2)
  
  # Plot linear expansion
  p$wm <- wmat0*1.1
  dR0_dwm <- (R(p)-R00) / (0.1*wmat0)
  lines(x=wmat/wmat0, y=1+(wmat-wmat0)*dR0_dwm/R00, lty='dotted', lwd=1)

  # Plot annotations:
  lines(x=c(1,1), y=c(1,1.25), lty='dotted')
  lines(x=1.2*c(1,1), y=c(1,1.25), lty='dotted')
  arrows(x0=1, x1=1.2, y0=1.2, y1=1.2, lwd=2, length=0.1)
  text(x=1.1, y=1.1, labels=TeX('$\\textit{S}_\\theta$'))
  text(x=1, y=1.27, labels=TeX('$\\theta$'))
  tmp <- sqrt(2*log(2))*(p$cv)
  arrows(x0=1-tmp, x1=1+tmp, y0=0.5,y1=0.5,code=3, lwd=2, length=0.1)
  text(x=1.1,y=0.45,labels=TeX('$\\sigma_\\theta$'))
}


plotSelectionResponse <- function(p=baseparamQG(), 
                                  ylim=c(-0.0025, 0.0065)) {
  defaultplothorizontal(mfcol=c(1,2))
  
  #
  # First panel, R vs. F
  #
  F <- seq(0,1.05,length.out = 20)
  
  defaultpanel(ylim=ylim, xlim=c(0,1), 
               xlab='Fishing mortality ($yr^{-1}$)',
               ylab="Rel. sel. response ($yr^{-1})$", label=TRUE)

  lines(F, 0*F, lty='dotted')
  S <- calcSelectionResponse(p=p, F=0, W=2000)

  for (i in 1:length(F)) 
    S[i,] <- calcSelectionResponse(p=p, F=F[i], W=2000)
  lines(F, S$dwmdt, lwd=2)
  lines(F, S$dAdt, lty='dotdash', lwd=2)
  lines(F, S$dkrdt, lty='dashed', lwd=2)
  legend(x=0, y=ylim[2], bty='n', lwd=2,
         legend=c('Size at maturation', 'Growth rate', 'Reproductive investment'),
         lty=c(1,4,2), y.intersp = 0.8)
  #
  # Second panel, R vs. W
  #
  F <- 0.3
  W <- 10^seq(log10(5), log10(25000), length.out = 20)
  
  semilogxpanel(xlim=c(min(W),2e4), ylim=ylim, xlab='Asymptotic size (g)', yaxis=FALSE, label=TRUE)
  lines(x=W, y=0*W,lty='dotted')

  S <- calcSelectionResponse(p=p, F=F, W=W)
  lines(W, S$dwmdt, lwd=2)
  lines(W, S$dAdt, lwd=2, lty='dotdash')
  lines(W, S$dkrdt, lwd=2, lty='dashed')
}

plotSelectionResponseLarge <- function() {
  p = baseparamQG()
  p$etaF = 0.5
  plotSelectionResponse(p)
}


plotSelectionResponsevsSpawner <- function(p=baseparamQG()) {
  par( mfcol=c(1,1), 
       oma=c(2.1, 2.5, 2.5, 0.5), # Outer margins (bottom, left, top, right)
       mar=c(0,0,0,0), # Margins
       ps=10, # Point size
       tcl=0.2, # Tick mark length
       mgp=c(1.1,0,0))  # Margins for axis title, axis labels, and axis line

  F <- seq(0, 0.3, length.out = 20)
  ylim <- c(-0.001, 0.002)
  plot(x=F, y=0*F, type='l', lty='dotted', ylim=ylim, cex.axis=0.9)
  mtext(side=2, line=1, TeX("Rel. sel. response ($yr^{-1})$"))
  mtext(side=1, line=1, TeX("Spawner $\\textit{F}$ ($yr^{-1})$"))
  axis(side=3, at=seq(0,0.3,by=0.05), labels=seq(0.3,0,by=-0.05), cex.axis=0.9)
  mtext(side=3, line=1, TeX("Size-based $\\textit{F}$ ($yr^{-1})$"))
  
  S <- calcSelectionResponse(p=p, F=max(F), W=20000)
  for (i in 1:length(F)) {
    p$fracSpawnerFishery <- F[i]/max(F)
    S[i,] <- calcSelectionResponse(p=p, F=max(F), W=20000)
  }
  lines(F, S$dwmdt, lwd=2)
  lines(F, S$dAdt, lwd=2, lty='dotdash')
  lines(F, S$dkrdt, lwd=2, lty='dashed')
}


plotSelection100 <- function() {
  library(deSolve)
  defaultplot(mfcol=c(1,2))
  
  years <- 100
  F <- 0.3
  W <- c(20,20000)
  for (i in 1:2) {
    p <- baseparamQG(wm=W[i]*0.28)
    ages <- seq(0, 5*ageMaturation(W[i],p), length.out=500)
    defaultpanel(ylim=c(0,W[i]), xlim=range(ages), 
                 xlab='Age (yrs)',
                 ylab="Weight (g)")
    # Before selection:
    out <- ode(y=p$w0, times=ages, func=function(t,w,p) list(growth(p,w)), parms=p)
    lines(x=out[,1], y=out[,2], lty='dotted')
    # Iterate over the years:
    p <- calcIteratedSelectionResponse(p,F=F, years=years)
    # After selection
    ages <- seq(0, 5*ageMaturation(p$W,p), length.out=500)
    out <- ode(y=p$w0, times=ages, func=function(t,w,p) list(growth(p,w)), parms=p)
    lines(x=out[,1], y=out[,2])
    
    # Iterate over the years:
    p <- baseparamQG(wm=W[i]*0.28)
    p$fracSpawnerFishery <- 1
    p <- calcIteratedSelectionResponse(p,F=F, years=years)
    # After selection
    ages <- seq(0, 5*ageMaturation(p$W,p), length.out=500)
    out <- ode(y=p$w0, times=ages, func=function(t,w,p) list(growth(p,w)), parms=p)
    lines(x=out[,1], y=out[,2])
  }
}


plotSelectionSelectivity <- function() {
  library(rootSolve)
  calcYieldError <- function(F, p) {
    p$F <- F
    Y <- NULL
    for (i in 1:length(F)) {
      p$F <- F[i]
      Y[i] <- calcYield(p, spectrum(p))-yield0
    }
    #cat(F,',',Y,'\n')
    return(Y)
  }

  W <- 20000
  p <- baseparamQG(wm=0.28*W)
  p$F <- 0.1
  yield0 <- calcYield(p, spectrum(p))

  defaultplothorizontal(mfrow=c(1,2))
  ylim <- c(-0.002, 0.002)
  
  #= Calc with trawl selectivity:  
  S <- calcSelectionResponse(p, W=W)
  etaF <- 10^seq(-2,0,length.out = 20)
  for (i in 1:length(etaF)) {
    p$etaF <- etaF[i]
    F[i] <- uniroot.all(calcYieldError, interval=c(0,10), p=p)[1]
    S[i,] <- calcSelectionResponse(p, F=F[i], W=W)
  }
  semilogxpanel(xlim=etaF, ylim=ylim,
                xlab('Lower size limit, $\\eta_\\textit{F}$'),
                ylab('Rel. sel. response ($yr^{-1}$)'), label=TRUE)
  lines(etaF, S$dwmdt, lwd=2)
  lines(etaF, S$dAdt, lty='dotdash', lwd=2)
  lines(etaF, S$dkrdt, lty='dashed', lwd=2)
  lines(etaF, 0*etaF, lty='dotted')
  lines(0.05*c(1,1), ylim, lty='dotted')

  #= Calc with upper limit:  
  p <- baseparamQG(wm=0.28*W)
  etaFF <- 10^seq(log10(0.3),0,length.out = 20)
  for (i in 1:length(etaFF)) {
    p$etaFF <- etaFF[i]
    F[i] <- uniroot.all(calcYieldError, interval=c(0,10), p=p)[1]
    S[i,] <- calcSelectionResponse(p, F=F[i], W=W)
  }
  semilogxpanel(xlim=etaFF, ylim=ylim,
                xlab='Upper size limit, $\\eta_{\\textit{FF}}', yaxis = FALSE, label=TRUE)
  lines(etaFF, S$dwmdt, lwd=2)
  lines(etaFF, S$dAdt, lty='dotdash', lwd=2)
  lines(etaFF, S$dkrdt, lty='dashed', lwd=2)
  lines(etaFF, 0*etaF, lty='dotted')
  
}


plotAllChapterFIE <- function() {
  pdfplot(FUN=plotConover, "ChapterFIE/Conover.pdf", width=doublewidth, height=height)

  pdfplot(FUN=plotQG, "ChapterFIE/QG.pdf", width=1.25*singlewidth, height=1.25*height)

  pdfplot(FUN=plotSelectionResponse, "ChapterFIE/SelectionResponse.pdf", width=doublewidth, height=height)
  
  pdfplot(FUN=plotSelectionResponseLarge, "ChapterFIE/SelectionResponseLarge.pdf", width=doublewidth, height=height)
  
  pdfplot(FUN=plotSelectionResponsevsSpawner, "ChapterFIE/SpawnerFishing.pdf", height=1.2*height)
  
  pdfplot("ChapterFIE/SelectionSelectivity.pdf", plotSelectionSelectivity, height=height, width=doublewidth)
}