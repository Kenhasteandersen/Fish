source("R/basetools.R")
source("R/basefunctions.R")
source("R/baseparameters.R")
source("R/community.R")
require(jpeg)

dir.create("ChapterDynamics")

calcEigenvalue <- function(p) {
  W <- p$W
  n <- p$n
  dt <- 0.1
  # Grid:
  grid = makegrid(p$w0, W, 0.01)
  w = grid$w
  dw = grid$dw
  nGrid = grid$nGrid
  # Growth and mortality:
  g <- growth(p, w)
  mu <- p$a*p$A*w^(p$n-1) + p$funcFishing(w,p)
  #
  # Set up diagonal of matrix 
  #
  A <- matrix(data=0, nrow=nGrid, ncol=nGrid)
  A[seq(from=2, by=nGrid+1, length.out=nGrid-1)] <- g[1:(nGrid-1)]/dw[2:nGrid]*dt
  A[seq(1, by=nGrid+1, length.out = nGrid)] <- 1 - (g/dw + mu)*dt  # diagonal
  #
  # Set up reproduction part of the matrix:
  #
  R <- p$epsEgg*psi(w/(p$etaM*p$W)) * p$A * W^(n-1) * w / p$w0 * dw
  A[1,] <- A[1,] + p$epsR*R*dt / dw[1]
  #
  # Find the eigenvalue and -vector with the largest real part
  #
  ev <- eigen(A)
  ixMax <- which.max( Re(ev$values) )
  
  fig <- ggplot() +
    geom_line(aes(x=log10(w), y=log10(abs(ev$vectors[,ixMax])))) +
    geom_line(aes(x=log10(w), y=log10(w^(-p$n-p$a))))
  fig
  
  return(list(r=log(Re(ev$values[ixMax]))/dt, vec=ev$vectors[,ixMax], w=w))
}


SimulateSpectrum <- function(p, tEnd, Ninit=NULL) {
  # Set up grid
  grid = makegrid(p$w0, p$W, 0.05)
  w = grid$w
  dw = grid$dw
  wminus <- grid$w - 0.5*dw # weight at left cell boundaries
  nGrid = grid$nGrid
  
  # Growth and mortality:
  #gminus <- growth(p, wminus)
  gplus <- pmax(0,growth(p, w+0.5*dw))
  gminus <- c(0,gplus[1:(nGrid-1)])
  mu <- as.vector(mortality(p,w))
  
  # Set up matrix coefficients
  dt <- p$dt
  A <- dt*(-1/8*gplus[2:nGrid] - 6/8*gminus[2:nGrid])/dw[2:nGrid]
  A[1] <- dt*(-1/8*gplus[2] - gminus[2])/dw[2]
  B <- 1+dt*(-3/8*gminus+6/8*gplus)/dw + dt*mu
  B[1] <- 1+dt*gplus[1]/dw[1] + dt*mu[1]
  B[2] <- 1+dt*6/8*gplus[2]/dw[2]+dt*mu[2]
  C <- dt*3/8*gplus[1:(nGrid-1)]/dw[1:(nGrid-1)]
  
  g <- growth(p, w)
  A1 <- -g[1:(nGrid-1)]*dt/dw[2:nGrid]
  B1 <- 1 + dt*(g/dw + mu)
  C1 <- 0*A1
  
  # Initial conditions
  nTime <- tEnd/dt
  N <- matrix(nrow=nTime, ncol=nGrid, data=0)
  if (is.null(Ninit))
    N[1,1] <- 1e-20
  else
    N[1,] <- Ninit
  N1 <- N
  R <- rep(0, nTime)
  SSB <- R
  # Iterate forward:
  for (iTime in 1:(nTime-1)) {
    # Recruitment:
    SSB[iTime] <- trapz(w, psi(w/(p$etaM*p$W))*w*N[iTime,] )
    Rp <- p$epsR*calcKr(p)*SSB[iTime]/p$w0 
    R[iTime] <- Rp/(Rp + 1)
    
    S <- N[iTime,] - c(0,0, dt/8*N[iTime,1:(nGrid-2)]*gminus[1:(nGrid-2)])/dw
    S[1] <- N[iTime,1] + dt*R[iTime]/dw[1]
    #S[2] <- N[iTime,2]
    N[iTime+1,] <- pmax(0, t(Solve.tridiag(A, B, C, S)))
    
    
    S1 <- N1[iTime,]
    S1[1] <- R[iTime]*dt/dw[1]
    N1[iTime+1,] <- Solve.tridiag(A1,B1,C1,S1)
    #SSB[iTime] <-
  }
  
  n <- 2000
  loglogpanel(xlim=c(1e-3,p$W), ylim=c(1e-10,10)*N[n,1])
  lines(w,N[n,],col="red")
  lines(w,N1[n,],col="black")
  lines(w,1e-2*w^{-p$n-p$a})
  
  ev <- calcEigenvalue(p)
  lines(ev$w, Re(ev$vec/ev$vec[1] * N[n,1]), lwd=2)
}


SimulateSpectrum2 <- function(p, tEnd, Ninit=NULL) {
  
  growth <- function(p,w) 
    p$A*w^p$n*(1-(w/p$W)^(1-p$n))
  
  # Set up grid
  nGrid <- 1000
  halfdw <- 0.5*(log(p$W)-log(p$w0))/nGrid
  w <- seq(log(p$w0)+halfdw, log(p$W)-halfdw, length.out = nGrid)
  dw <- rep(w[2]-w[1],nGrid)
  
  # Growth and mortality:
  #gminus <- growth(p, wminus)
  gplus <- pmax(0,growth(p, exp(w+0.5*dw)))
  gminus <- c(0,gplus[1:(nGrid-1)])
  mu <- as.vector(mortality(p,exp(w)))
  
  # Set up matrix coefficients
  dt <- p$dt
  A <- dt*(-1/8*gplus[2:nGrid] - 6/8*gminus[2:nGrid])/dw[2:nGrid] * exp(-w[2:nGrid])
  g1 <- growth(p,exp(w[1]))
  A[1] <- dt*(-1/8*gplus[2] - g1)/dw[2] * exp(-w[2])
  B <- 1+dt*(-3/8*gminus+6/8*gplus)/dw * exp(-w) + dt*mu
  B[1] <- 1+dt*g1/dw[1]*exp(-w[1]) + dt*mu[1]
  B[2] <- 1+dt*6/8*gplus[2]/dw[2]*exp(-w[2])+dt*mu[2]
  C <- dt*3/8*gplus[1:(nGrid-1)]/dw[1:(nGrid-1)]*exp(-w[1:(nGrid-1)])
  
  g <- growth(p, exp(w))
  A1 <- -g[1:(nGrid-1)]*dt/dw[2:nGrid]*exp(-w[2:nGrid])
  B1 <- 1 + dt*(g/dw*exp(-w) + mu)
  C1 <- 0*A1
  
  # Initial conditions
  nTime <- tEnd/dt
  N <- matrix(nrow=nTime, ncol=nGrid, data=0)
  if (is.null(Ninit))
    N[1,1] <- 1e-20
  else
    N[1,] <- Ninit
  N1 <- N
  R <- rep(0, nTime)
  SSB <- R
  # Iterate forward:
  for (iTime in 1:(nTime-1)) {
    # Recruitment:
    SSB[iTime] <- trapz(exp(w), psi(exp(w)/(p$etaM*p$W))*exp(w)*N[iTime,] )
    Rp <- p$epsR*calcKr(p)*SSB[iTime]/p$w0 
    R[iTime] <- 1#Rp/(Rp + 1)
    
    S <- N[iTime,] - dt/8*c(0,0,N[iTime,1:(nGrid-2)])*gminus/dw*exp(-w) 
    # The factor 1.3 seems to match stuff up nicely. Wonder why?
    S[1] <- N[iTime,1] + 1.3*dt*R[iTime]/dw[1]*exp(-w[1])
    S[2] <- N[iTime,2]
    N[iTime+1,] <- pmax(0, t(Solve.tridiag(A, B, C, S)))
    
    
    S1 <- N1[iTime,]
    S1[1] <- R[iTime]*dt/dw[1]*exp(-w[1])
    N1[iTime+1,] <- Solve.tridiag(A1,B1,C1,S1)
    #SSB[iTime] <-
  }
  
  defaultplot(mfcol=c(2,1))
  n <- nTime
  loglogpanel(xlim=c(1e-3,p$W), ylim=c(1e-10,10)*N[n,1])
  lines(exp(w),N[n,],col="red")
  lines(exp(w),N1[n,],col="black")
  #lines(exp(w),1e-2*exp(w)^{-p$n-p$a})
  
  analytical <- spectrumana(p)
  lines(analytical$w, analytical$NprR, col="blue",lty="dashed")
  #ev <- calcEigenvalue(p)
  #lines(ev$w, Re(ev$vec/ev$vec[1] * N[n,1]), lwd=2)
  
  #semilogypanel(xlim=c(0,nTime*dt), ylim=range(SSB))
  plot(SSB)
  
  return(list(N,N1,SSB,R))
}


SimulateSimplePopGrowth <- function(p, tEnd) {
  
  growth <- function(p,w) 
    p$A*w^p$n
  
  # Set up grid
  nGrid <- 500
  halfdw <- 0.5*(log(p$W)-log(p$w0))/nGrid
  w <- seq(log(p$w0)+halfdw, log(p$W)-halfdw, length.out = nGrid)
  dw <- rep(w[2]-w[1],nGrid)
  
  # Growth and mortality:
  #gminus <- growth(p, wminus)
  gplus <- pmax(0,growth(p, exp(w+0.5*dw)))
  gminus <- c(0,gplus[1:(nGrid-1)])
  mu <- as.vector(mortality(p,exp(w)))
  
  # Set up matrix coefficients
  dt <- p$dt
  A <- dt*(-1/8*gplus[2:nGrid] - 6/8*gminus[2:nGrid])/dw[2:nGrid] * exp(-w[2:nGrid])
  g1 <- growth(p,exp(w[1]))
  A[1] <- dt*(-1/8*gplus[2] - g1)/dw[2] * exp(-w[2])
  A[nGrid-1] <- -dt*6/8*gminus[nGrid]/dw[nGrid]*exp(-w[nGrid])
  B <- 1+dt*(-3/8*gminus+6/8*gplus)/dw * exp(-w) + dt*mu
  B[1] <- 1+dt*g1/dw[1]*exp(-w[1]) + dt*mu[1]
  B[2] <- 1+dt*6/8*gplus[2]/dw[2]*exp(-w[2])+dt*mu[2]
  B[nGrid] <- 1+dt*(-3/8*gminus[nGrid]+6/8*gplus[nGrid])/dw[nGrid] * exp(-w[nGrid]) + dt*mu[nGrid]
  C <- dt*3/8*gplus[1:(nGrid-1)]/dw[1:(nGrid-1)]*exp(-w[1:(nGrid-1)])
  
  g <- growth(p, exp(w))
  A1 <- -g[1:(nGrid-1)]*dt/dw[2:nGrid]*exp(-w[2:nGrid])
  B1 <- 1 + dt*(g/dw*exp(-w) + mu)
  C1 <- 0*A1
  
  # Initial conditions
  nTime <- tEnd/dt
  N <- matrix(nrow=nTime, ncol=nGrid, data=0)
  N[1,nGrid] <- 1e-20
  N1 <- N
  R <- rep(0, nTime)
  R1 <- R
  SSB <- R
  # Iterate forward:
  for (iTime in 1:(nTime-1)) {
    # Recruitment:
    SSB[iTime] <- trapz(exp(w), psi(exp(w)/(p$etaM*p$W))*exp(w)*N[iTime,] )
    Rp <- p$epsR*p$epsEgg*growth(p,p$W)*N[iTime,nGrid]*p$W/p$w0
    R[iTime] <- Rp #Rp/(Rp + 1)
    
    S <- N[iTime,] - dt/8*c(0,0,N[iTime,1:(nGrid-2)])*gminus/dw*exp(-w) 
    S[1] <- N[iTime,1] + 1.3*dt*R[iTime]/dw[1]*exp(-w[1])
    S[2] <- N[iTime,2]
    N[iTime+1,] <- pmax(0, t(Solve.tridiag(A, B, C, S)))
    
    R1[iTime] <- p$epsR*p$epsEgg*growth(p,p$W)*N1[iTime,nGrid]*p$W/p$w0
    S1 <- N1[iTime,]
    S1[1] <- R1[iTime]*dt/dw[1]*exp(-w[1])
    N1[iTime+1,] <- Solve.tridiag(A1,B1,C1,S1)
    #SSB[iTime] <-
  }
  
  defaultplot(mfcol=c(2,1))
  n <- nTime
  loglogpanel(xlim=c(1e-3,p$W), ylim=c(1e-10,1))
  lines(exp(w), N[n,]/N[n,1],col="red")
  lines(exp(w), N1[n,]/N1[n,1],col="black")
  epsilon <- p$epsEgg*p$epsR
  n <- p$n
  rana1 <- (p$A*(-1 + n)*p$w0^n*p$W^n*log((p$a*p$w0^(1 - p$a)*p$W^(-1 + p$a))/(epsilon)))/
    (p$w0^n*p$W - p$w0*p$W^n)
  nn <- (exp(w)/p$w0)^(-n-p$a) * exp(-rana1*exp(w)^(1-n) / (p$A*(1-n))) / exp(-rana1*p$w0^(1-n) / (p$A*(1-n)))
  lines(exp(w), nn/nn[1], col="blue", lty='dashed')
  
  t <- p$dt:nTime*p$dt
  semilogypanel(xlim=range(t), ylim=range(exp(rana1*t)))
  lines(t, R/R[1],col="red")
  lines(t, R1/R1[1])
  lines(t, exp(rana1*t),col="blue")
  
  return(list(t=t,N=N,N1=N1,SSB=SSB,R=R, R1=R1))
}


plotEigen <- function()
{
  p <- baseparameters()
  W <- 10^seq(-.5,6,length.out = 25)
  #
  # Full numerical solution:
  #
  rfull <- 0*W
  for (i in 1:length(W)) {
    p$W <- W[i]
    rfull[i] <- calcEigenvalue(p)$r
  }
  #
  # First analytical approximation:
  #
  rana1 <- calc_rana1(W,p)
  #
  # Second analytical approximation:
  #
  rana2 <- calc_rana2(W)
  #
  # Data from Hutchings et al 2012
  #
  data <- read.table(file="Data/Hutchings2012.dat", header=TRUE)
  data$A <- (p$etaM*data$Maximum_weight_.kg.*1000)^(1-p$n) / ((1-p$n)*data$age_at_maturity_.years.)
  dat <- data[which( data$class=="Actinopterygii"),]
  #dat2 <- data[which( data$class=="Chondrichthyes"),]
  
  defaultplot(mfcol=c(1,2), mar=c(2.1,2.1,0,1.3))
  semilogxpanel(xlim=W, ylim=c(-0.1, 1.5),
                xlab('Asymptotic weight (g)'),
                ylab('Population growth rate ($yr^{-1}$)'), label=TRUE)
  points(dat$Maximum_weight_.kg.*1000, dat$rmax, pch=16)
  #points(dat2$Maximum_weight_.kg.*1000, dat2$rmax, col="grey")
  
  lines(W,rfull, lwd=2)
  lines(W, rana1, lwd=2, lty='dashed')
  lines(W, rana2, lwd=2, col='grey')
  lines(W,0*W, lwd=1, lty='dotted')
  
  W = 10000
  p$W <- W
  ev <- calcEigenvalue(p)
  w = 10^seq(-3,log10(W),length.out = 100)
  loglogpanel(xlim=w, ylim=c(1e-10,10), label=TRUE,
              xlab="Weight (g)",
              ylab="Number spectrum $\\textit{n}(\\textit{w})/\\textit{n}(\\textit{w}_0)$")

  analeigen <- function(r) {
    (w/p$w0)^(-p$n-p$a) * exp(-r*w^(1-p$n) / (p$A*(1-p$n))) /
     exp(-r*p$w0^(1-p$n) / (p$A*(1-p$n)))
  }
  lines(w, analeigen(calc_rana2(W)), lwd=2, col='grey')
  lines(w, analeigen(calc_rana1(W)), lwd=2, lty='dashed')
  lines(w, (w/p$w0)^(-p$n-p$a), lwd=1, lty='dotted')
  lines(ev$w, Re(ev$vec/ev$vec[1]), lwd=2)
}

# plotDynamics <- function(W=10000) {
#   p <- paramTraitbasedmodel(nSpecies = 1, W=W)
#   p$Rmax <- 1e-10 # simulate a species that has no food competition
#   p$mu0Coefficient <- 0
#   p$SizeBasedMortCoef <- p$a*p$A # Set up a size based mortality mimicing predation from the community
#   p$GridExpand <- 0.1
#   # Make initial conditions
#   spec0 <- runspectrum(p)
#   factorReduction <- 1e-2
#   spec0$N <- spec0$N*factorReduction
#   
#   # Simulate forward
#   p$stepSave <- 1
#   p$dt <- 0.01
#   p$Tend <- 40
#   spec1 <- runspectrum(p, prev_run = spec0)
#   
#   res <- SimulateSimplePopGrowth(p,40)
#   
#   defaultplot(mfcol=c(2,1))
#   
#   semilogypanel(xlim=c(0,p$Tend), ylim=c(factorReduction,1))
#   t <- spec1$R$time
#   lines(t, spec1$R$R)
#   lines(t, spec1$R$SSB/spec1$R$SSB[spec1$nIte])
#   lines(res$t, res$R, lty="dashed", col="red")
#   lines(res$t, res$R1, lty="dashed")
#   
#   
#   ev <- calcEigenvalue(p)
#   lines(t, factorReduction*exp(ev$r*t), col="red")
#   lines(t, spec1$N[,1,2]/spec1$N[spec1$nIte,1,2], col="blue")
#   #lines(spec1$w, spec1$N[spec1$nSave,1,],col="red")
#   
#   loglogpanel(xlim=c(p$w0, p$W), ylim=c(1e-6,1))
#   lines(ev$w, Re(ev$vec/ev$vec[1]), col="red")
#   lines(spec1$fish$w, spec1$N[spec1$nIte,1,]/spec1$N[spec1$nIte,1,1])
# }

calcDynamics <- function(W, tEnd) {
  p <- paramTraitbasedmodel(nSpecies = 1, W=W)
  p$Rmax <- 1e-10 # simulate a species that has no food competition
  p$mu0Coefficient <- 0 # No background mortality
  p$SizeBasedMortCoef <- p$a*p$A # Set up a size based mortality mimicing predation from the community
  p$GridExpand <- 0.1
  p$dt <- 0.01
  # Make initial conditions
  spec0 <- runspectrum(p)
  factorReduction <- 1e-2  # the reduction in biomass prior to the "recovery"
  spec0$N <- spec0$N*factorReduction
  
  # Simulate forward
  p$stepSave <- 1
  p$Tend <- tEnd
  spec <- runspectrum(p, prev_run = spec0)
  
  return(list(t=spec$R$time, R=spec$R$R, SSB=spec$R$SSB))
}


plotDynamics <- function() {
  W <- c(20,20000)
  tEnd <- 20
  
  defaultplotvertical(mfcol=c(2,1))
  
  for (i in 1:2) {
    # Full structured calculation:
    dyn <- calcDynamics(W[i], tEnd)
    
    # Unstructured approximation:
    #p <- baseparameters(W[i])
    #deriv <- function(t, B, p) {
    #  RoverRmax <- p$epsR*p$epsEgg*p$A*p$W^(p$n)/p$w0*B
    #  list( calc_rana1(p$W,p)* (1-B) * B )
    #}
    #out <- ode(y=0.01, times=seq(0,20,by=0.1), func=deriv, parms=p)
    
    if (i==2)
      xlab='Time (years)'
    else
      xlab=''
    semilogypanel(xlim=c(0,tEnd), ylim=c(0.01,3), xaxis=(i!=1),
                  xlab=xlab, label = TRUE)
    lines(dyn$t, dyn$SSB/max(dyn$SSB), lwd=2)
    lines(dyn$t, dyn$R, col=stdgrey, lwd=2)
    #lines(out[,1]+ageMaturation(W[i],p), out[,2], lty=dotted)
    
    if (i==1) {
      mtext('     SSB/max(SSB)', side=left, line=1.5, adj=0)
      img<-readJPEG("data/Sardinops_sagax open.jpg")
      rasterImage(img, 15, 0.02, 18, 0.05)
    }
    else {
      mtext(TeX('  R/$R_{max};   $'), side=left, line=1.5, col=stdgrey, adj=1)
      img <- readJPEG("data/Atlantic_cod open.jpg")
      rasterImage(img, 11, 0.02, 19, 0.25)
    }
    
    # Draw vertical lines
    t1 <- dyn$t[(which((dyn$SSB/dyn$SSB[1])>1.1))[1]]
    if (i==2)
      t2 <- 8.5
    else
      t2 <- 4
    vline(x=t1, lty="dotted")
    vline(x=t2, lty="dotted")
    
    # Numbers to refer to regions
    x <- c(t1/2, t1+(t2-t1)/2, t2+(tEnd-t2)/2)
    y <- 2*c(1,1,1)
    points(x=x, y=y, pch=16, col="grey", cex=2.5)
    text(x=x[1], y=2,labels='1')
    text(x=x[2], y=2,labels='2')
    text(x=x[3], y=2,labels='3')
    
    #vline(x=dyn$t[(which(dyn$R/dyn$R[length(dyn$R)]>0.95))[1]], lty="dotted")
  }  
  xlab("Time (years)")
}


plotRecoveryExample <- function() {
  p <- paramTraitbasedmodel(nSpecies = 1, W=20000)
  p$dt <- 0.025
  p$Rmax <- 1e-10 # simulate a species that has no food competition
  p$mu0Coefficient <- 0
  p$SizeBasedMortCoef <- p$a*p$A # Set up a size based mortality mimicing predation from the community
  
  # 1. Simulate an unexploited species:
  spec0 <- runspectrum(p)
  p$Tend <- 5
  p$stepSave <- 1
  spec0 <- runspectrum(p, prev_run = spec0)
  # 2. Simulate over-exploitation:
  p$F <- 1.5
  p$Tend <- 15
  spec1 <- runspectrum(p, prev_run = spec0)
  # Simulate  recovery:
  p$F <- 0
  spec2 <- runspectrum(p, prev_run=spec1)
  # Simulate Fmsy:
  p$F <- calcRefpoints(p)$Fmsy
  p$Tend <- 50
  specMSY <- runspectrum(p, prev_run=spec1)
  
  idxMSY <- 1:which( spec2$R$SSB>=specMSY$R$SSB[specMSY$nSave])
  time <- c(spec0$R$time, spec1$R$time+5, spec2$R$time[idxMSY]+20, spec2$R$time[max(idxMSY)]+20, 30)
  SSB <- c(c(spec0$R$SSB, spec1$R$SSB, spec2$R$SSB[idxMSY]) / specMSY$R$SSB[specMSY$nSave], 1,1)
  Yield <- c(c(spec0$R$Y, spec1$R$Y, spec2$R$Y[idxMSY]) / specMSY$R$Y[specMSY$nSave], 1,1)
  
  defaultplot()
  plot(time, SSB, type="l", lwd=2, ylim=c(0,10), xlab="Time (years)", ylab="")
  lines(time, Yield, col=stdgrey, lwd=2)
  hline(1)
  vline(5)
  vline(20)
  vline(time[max(idxMSY)]+20)
  text(x=2.5, y=9.5, labels=TeX("\\textit{F} = 0"))
  text(x=12.5, y=9.5, labels=TeX("\\textit{F} = 1.5 yr^{-1}"))
  text(x=23, y=9.5, labels=TeX("\\textit{F} = 0"))
  text(x=28.5, y=9.35, labels=TeX("\\textit{F} = \\textit{F}_{MSY}"))
  mtext(TeX('\\textit{Y}/\\textit{Y}_{MSY}'), left, 1, outer=FALSE, adj=0, col=stdgrey)
  mtext(TeX('\\textit{B}_{SSB}/\\textit{B}_{MSY}'), left, 1, outer=FALSE, adj=1)
}

plotRecovery <- function() {
  p <- paramTraitbasedmodel(nSpecies = 1, W=20000)
  p$dt <- 0.025
  p$tEnd <- 50
  p$Rmax <- 1e-10 # simulate a species that has no food competition
  p$mu0Coefficient <- 0
  p$SizeBasedMortCoef <- p$a*p$A # Set up a size based mortality mimicing predation from the community
  # Simulate an over-exploited species
  p$F <- 1.5
  spec0 <- runspectrum(p)
  
  #loglogpanel(xlim=c(1e-3, 1000), ylim=c(1e-20,1))
  #lines(spec0$w, spec0$N[spec0$nSave,1,])
  
  # Simulate Fmsy recovery
  Fmsy <- calcRefpoints(p)$Fmsy
  p$F <- Fmsy
  specMSY <- runspectrum(p, prev_run=spec0)
  BMSY <- specMSY$R$SSB[specMSY$nSave]
  
  defaultplot(mfcol=c(1,2))
  
  defaultpanel(xlim=c(0,20), ylim=c(0,1.2), 
               xlab="Time (years)", ylab="\\textit{B}_{SSB}/\\textit{B}_{MSY}", label = TRUE)
  # Simulate recovery scenarios
  F <- seq(0,Fmsy,length.out = 5)
  p$tEnd <- 20
  trec <- 0*F
  for (i in 1:length(F)) {
    p$F <- F[i]
    p$stepSave <- 1
    spec <- runspectrum(p, prev_run = spec0)
    SSB <- spec$R$SSB
    SSB[SSB>=BMSY]<-NaN
    trec[i] <- spec$R$time[ which( spec$R$SSB/BMSY >= .99)[1] ]
    lines(spec$R$time, SSB/BMSY)
    points(trec[i],1,pch=16)
  }
  hline(y=spec0$R$SSB[spec0$nSave]/BMSY)
  hline(y=1)
  text(trec[1], 1.1, TeX("\\textit{F} = 0"))
  text(16.5, 1.1, TeX("\\textit{F} = \\textit{F}_{MSY}"))
  
  # Do it again for two species
  W <- c(20,20000)
  defaultpanel(xlim=c(0,0.5),ylim=c(0,20),label=TRUE,
               xlab="\\textit{F} ($yr^{-1}$)", ylab="Recovery time (yr)")
  for (j in 1:2) {
    p$W <- W[j]
    # over exploited run:
    p$F <- 1.5
    spec0 <- runspectrum(p)
    # recovered state:
    Fmsy <- calcRefpoints(p)$Fmsy
    p$F <- Fmsy
    specMSY <- runspectrum(p, prev_run=spec0)
    BMSY <- specMSY$R$SSB[specMSY$nSave]
    # Simulate recovery scenarios:
    F <- seq(0,Fmsy,length.out = 25)
    p$tEnd <- 20
    trec <- 0*F
    for (i in 1:length(F)) {
      p$F <- F[i]
      p$stepSave <- 1
      spec <- runspectrum(p, prev_run = spec0)
      SSB <- spec$R$SSB
      SSB[SSB>=BMSY]<-NaN
      trec[i] <- spec$R$time[ which( spec$R$SSB/BMSY >= .99)[1] ]
    }
    lines(F, trec, lwd=j)
    vline(Fmsy)
  }
}

plotNoisyRecruitment <- function() {
  #require(gsl)  
  W <- c(20, 20000)
  defaultplothorizontal(mfcol=c(1,2))
  
  for (i in 1:length(W)) {
    p <- paramTraitbasedmodel(nSpecies = 1, W=W[i])
    p$dt <- 0.01
    p$Tend <- 1000
    p$Rmax <- 1e-10 # simulate a species that has no food competition
    p$mu0Coefficient <- 0
    p$SizeBasedMortCoef <- p$a*p$A # Set up a size based mortality mimicing predation from the community
    # noisy recruitment:
    p$thetaRecruitment <- 1
    p$sigmaRecruitment <- 2
    
    res <- runspectrum(p)
    
    t<-res$R$time
    R <- res$R$R
    SSB <- res$R$SSB
    
    semilogypanel(xlim=c(0,67), ylim=c(0.01,100), label=TRUE, 
                  yaxis = (i==1))
    idx = which((t>200) & (t<250))
    lines(t[idx]-200, R[idx], col=stdgrey)
    lines(t[idx]-200, SSB[idx]/mean(SSB))
    hline(1)
    
    idx <- which(t>50)
    # density of SSB:
    d <- density(log10(SSB[idx]/mean(SSB[idx])), adjust=2)
    polygon(x=52+5*d$y, y=10^d$x, col='black', border=NA)
    # density of R
    d <- density(log10(R[idx]))
    polygon(x=52+5*d$y, y=10^d$x, col=grey(0.4, alpha=0.5), border=NA)
    
  }
  mtext('Time (years)', bottom, 1, outer=TRUE, las=0)
  mtext(TeX('\\textit{R}/\\textit{R}_{max}'), left, 1.5, outer=TRUE, adj=0, col='grey')
  mtext(TeX('\\textit{B}_{SSB}/mean(\\textit{B}_{SSB})'), left, 1.5, outer=TRUE, adj=1)
  
}


plotNoisyFishingOld <- function() {
  F <- seq(0, 3 ,length.out = 20)
  W <- c(20, 20000)
  CV <- matrix(ncol=length(F), nrow=length(W))
  R0 <- CV
  ref <- list()
  
  if (bRecalcExpensiveFunctions) {
    for (j in 1:length(W)) {  
      res <- list()
      p <- paramTraitbasedmodel(nSpecies = 1, W=W[j])
      p$dt <- 0.01
      p$Tend <- 10000
      p$Rmax <- 1e-10 # simulate a species that has no food competition
      p$mu0Coefficient <- 0
      p$SizeBasedMortCoef <- p$a*p$A # Set up a size based mortality mimicing predation from the community
      # noisy recruitment:
      p$thetaRecruitment <- 1
      p$sigmaRecruitment <- 2
      ref[[j]] <- calcRefpoints(p)
      for (i in 1:length(F)) {
        p$F <- F[i]
        res[[i]] <- runspectrum(p)
        ix <- seq(0.75*res[[i]]$nSave, res[[i]]$nSave)
        #CV[j,i] <- std(log(res[[i]]$R$SSB[ix]))/mean(log(res[[i]]$R$SSB[ix]))
        #mn <- exp(mean(log(res[[i]]$R$SSB)[ix]))
        CV[j,i] <- std(log(res[[i]]$R$SSB[ix]))
        R0[j,i] <- mean(res[[i]]$R$R0[ix])
      }
    }
    save(CV,R0,ref, file="Data/NoisyFishing.RData")
  } else
    load("Data/NoisyFishing.RData")
  
  defaultplot()
  defaultpanel(xlim=F, ylim=c(0,6), 
               xlab="Fishing mortality $(yr^{-1})$",
               ylab="Coef. of variation")
  for (j in 1:length(W)) {
    lines(F, CV[j,], lwd=j)
    points(ref[[j]]$Fmsy, interp1(F,CV[j,],ref[[j]]$Fmsy), pch=19)
    points(ref[[j]]$Fcrash, interp1(F,CV[j,],ref[[j]]$Fcrash), pch=19, col="grey", bg="grey")
  }
}

plotNoisyFishing <- function() {
  library(parallel)

  F <- seq(0, 2.5 ,length.out = 24)
  W <- c(20, 20000)
  CV <- matrix(ncol=length(F), nrow=length(W))
  R0 <- CV
  r <- CV
  ref <- list()
  
  processInput <- function(i) {
    p$F <- F[i]
    res <- runspectrum(p)
    ix <- seq(0.75*res$nSave, res$nSave)
    list(CV=std(log(res$R$SSB[ix])), 
         R0=mean(res$R$R0[ix]),
         r=lm(formula = log(res$R$SSB) ~ res$R$time)$coefficients[[2]])
  }
  
  if (bRecalcExpensiveFunctions) {
    for (j in 1:length(W)) {  
      res <- list()
      p <- paramTraitbasedmodel(nSpecies = 1, W=W[j])
      p$dt <- 0.01
      p$Tend <- 10000
      p$Rmax <- 1e-10 # simulate a species that has no food competition
      p$mu0Coefficient <- 0
      p$SizeBasedMortCoef <- p$a*p$A # Set up a size based mortality mimicing predation from the community
      # noisy recruitment:
      p$thetaRecruitment <- 1
      p$sigmaRecruitment <- 2
      ref[[j]] <- calcRefpoints(p)
      results = mcmapply(processInput, 1:length(F), mc.cores=detectCores())
      CV[j,] <- unlist(results[1,])
      R0[j,] <- unlist(results[2,])
      r[j,] <- unlist(results[3,])
    }
    save(CV,R0,r,ref, file="ChapterDynamics/NoisyFishing.RData")
  } else
    load("ChapterDynamics/NoisyFishing.RData")
  
  defaultplot()
  defaultpanel(xlim=F, ylim=c(0,6), 
               xlab="Fishing mortality $(yr^{-1})$",
               ylab="Coef. of variation")
  for (j in 1:length(W)) {
    lines(F, CV[j,], lwd=j)
    ix <- r[j,]< -0.01
    lines(F[ix], CV[j,ix], col="white", lwd=2*j)
    points(ref[[j]]$Fmsy, interp1(F,CV[j,],ref[[j]]$Fmsy), pch=19)
    points(ref[[j]]$Fcrash, interp1(F,CV[j,],ref[[j]]$Fcrash), pch=19, col="grey", bg="grey")
  }
}

testNoisyFishing <- function() {
  library(parallel)

  F <- seq(0, 2 ,length.out = 8)
  W <- 2000
  res <- list()
  CV <- matrix(ncol=length(F), nrow=1)
  R0 <- CV
  
  p <- paramTraitbasedmodel(nSpecies = 1, W=W)
  p$dt <- 0.01
  p$Tend <- 50
  p$stepSave <- 1
  p$Rmax <- 1e-10 # simulate a species that has no food competition
  p$mu0Coefficient <- 0
  p$SizeBasedMortCoef <- p$a*p$A # Set up a size based mortality mimicing predation from the community
  # noisy recruitment:
  p$thetaRecruitment <- 1
  p$sigmaRecruitment <- 2
  #ref[[j]] <- calcRefpoints(p)
  
  cl <- makeCluster(4)

  processInput <- function(i)   {
    pp <- p
    pp$F <- F[i]
    res[[i]] <- runspectrum(pp)
    ix <- seq(0.75*res[[i]]$nSave, res[[i]]$nSave)
    #CV[j,i] <- std(log(res[[i]]$R$SSB[ix]))/mean(log(res[[i]]$R$SSB[ix]))
    mn <- exp(mean(log(res[[i]]$R$SSB)[ix]))
    CV <- std(log(res[[i]]$R$SSB[ix]/mn))
    R0 <- mean(res[[i]]$R$R0[ix])
    list(t=res[[i]]$R$time, SSB=res[[i]]$R$SSB)
  }
  
  results = mcmapply(processInput, 1:length(F), mc.cores=4)
  
  defaultplot()
  semilogypanel(xlim=c(0,500), ylim=range(c(results[2,1]$SSB,results[2,dim(results)[2]]$SSB)))
  for (i in 1:length(F)) {
    lines(results[1,i]$t, results[2,i]$SSB)
  }
  
  
}

testcalcUnstructuredDynamics <- function(W) {
  p <- baseparameters(W=W)
  p$r <- calc_rana1(W)
  p$Bmax <- 1-1/(p$epsR*p$epsEgg)* (p$w0/W)^(1-p$a)
  deriv <- function(t, B, p) {
    RoverRmax <- p$epsR*p$epsEgg*p$A*W^(p$n)/p$w0*B
    print( RoverRmax )
    #p$epsR <- p$epsR*(1-B)
    list( calc_rana1(W,p)* (1-B) * B )
    #list( calc_rana1(W,p) * B )
    #print( p$A*(p$n-1)*p$w0^(p$n-1) * log(RoverRmax) )
    #list( p$A*(p$n-1)*p$w0^(p$n-1) * log(RoverRmax) * B)
  }
  out <- ode(y=0.01, times=seq(0,20,by=0.1), func=deriv, parms=p)
  plot(out, log="y")
}



plotAllChapterDynamics <- function() {
  pdfplot(FUN=plotEigen, "ChapterDynamics/eigen.pdf", width=doublewidth, height=height)
  pdfplot(FUN=plotDynamics, "ChapterDynamics/dynamics.pdf", width=1.5*singlewidth, height=1.5*height)
  pdfplot(FUN=plotRecoveryExample, "ChapterDynamics/recoveryexample.pdf", width=doublewidth, height=height)
  pdfplot(FUN=plotRecovery, "ChapterDynamics/recovery.pdf", width=doublewidth, height=height)
  pdfplot(FUN=plotNoisyRecruitment, "ChapterDynamics/noisyrecruitment.pdf", width=doublewidth, height=height)
  #pdfplot(FUN=plotNoisyFishing, "ChapterDynamics/noisyfishing.pdf", width=singlewidth, height=height)
}
