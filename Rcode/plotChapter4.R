#
# Chapter 4
#
source("Rcode/basetools.R")
source("Rcode/basefunctions.R")
source("Rcode/baseparameters.R")

# 
# Plots
# 
plotSpectrum <- function()
{
  p <- baseparameters()
  W <- 1
  p$W <- W
  N <- spectrum(p)
  #N$NprR <- N$NprR/calcSSB(p,N)
  Nana <- spectrumana(p)
  #Nana$NprR <- Nana$NprR/calcSSB(p,Nana)
  
  fig <- ggplot() +
    geom_line(data=N, aes(x=w, y=NprR/NprR[1]), size=1) +
    geom_line(data=Nana, aes(x=w, y=NprR/NprR[1]), size=0.5) +
    geom_line(data=N, aes(x=w, y=(w/p$w0)^(-p$n-p$a)), size=0.5, linetype="dashed") +
    geom_line(aes(x=p$etaM*c(1,1), y=c(1e-5,10)), linetype="dotted", size=thin) +
    xlab(TeX("Rel. weight $\\textit{w}/\\textit{W}_\\infty$")) + 
    ylab("Number spectrum") 
  fig <- loglog(fig, ylim=c(1e-4,6), xlim=c(1e-3, 1), label="a")

  fig1 <- ggplot() +
    geom_line(data=N, aes(x=w, y=(NprR/NprR[1])*(w/p$w0)^2), size=1) +
    geom_line(data=Nana, aes(x=w, y=NprR/NprR[1]*(w/p$w0)^2), size=0.5) +
    geom_line(data=N, aes(x=w, y=(w/p$w0)^(2-p$n-p$a)), size=0.5, linetype="dashed") +
    geom_vline(xintercept = p$etaM*c(1,1), linetype="dotted", size=thin) +
    xlab(TeX("Rel. weight $\\textit{w}/\\textit{W}_\\infty$")) + 
    ylab("Sheldon spectrum")
  fig1 <- loglog(fig1, xlim=c(1e-3, 1), ylim=c(1,1.5e3), label="b")
  
  #pdf("spectra.pdf", width=doublewidth, height=height)
  #multiplot(fig, fig1, cols=2)
  #dev.off()
  grid <- plot_grid(fig, fig1, align="h", ncol=2)
  ggsave("Chapter4/spectra.pdf", grid, width=doublewidth, height=height)
  grid
}

plotCohort <- function()
{
  p <- baseparameters()
  W <- 1
  p$W <- W
  N <- spectrum(p)
  #N$NprR <- N$NprR/calcSSB(p,N)
  Nana <- spectrumana(p)
  
  fig <- ggplot() +
    geom_line(data=N, aes(x=w, y=NprR/NprR[1] * w/p$w0 * g/g[1]), size=1) +
    geom_line(data=Nana, aes(x=w, y=NprR/NprR[1] * w/p$w0 *g/g[1]), size=0.5) +
    geom_line(data=N, aes(x=w, y=(w/p$w0)^(1-p$a)), size=0.5, linetype="dashed") +
    geom_vline(xintercept=p$etaM*c(1,1), linetype="dotted", size=thin) +
    labs(x=TeX('Relative weight $\\textit{w}/\\textit{W}_\\infty$')) + 
    ylab("Cohort biomass")
  fig <- loglog(fig, ylim=c(1,100))
  ggsave("Chapter4/cohort.pdf",width=singlewidth,height=height)
  fig  
}

plotCohort_for_presentation <- function() {
  p <- baseparameters()
  W <- c(10, 1000, 100000)
  
  for (i in 1:length(W)) {
    p$W <- W[i]
    N <- spectrum(p)
    #N$NprR <- N$NprR/calcSSB(p,N)
    Nana <- spectrumana(p)
    
    if (i==1) 
      fig <- ggplot();
    fig <- fig +
      geom_line(data=N, aes(x=w, y=NprR/NprR[1] * w/p$w0 * g/g[1]), size=1, color="blue") +
      #geom_line(data=Nana, aes(x=w, y=NprR/NprR[1] * w/p$w0 *g/g[1]), size=0.5) +
      geom_line(data=N, aes(x=w, y=(w/p$w0)^(1-p$a)), size=0.5, linetype="dashed") +
      #geom_vline(xintercept=p$etaM*c(1,1), linetype="dotted", size=thin) +
      labs(x=TeX('Weight (g)')) + 
      ylab("Increase in cohort biomass")
    fig <- loglog(fig, ylim=c(1,1e6), xlim=c(1e-3, max(W)))
    ggsave(paste("Chapter4/cohort Presentation",i,".pdf"),width=doublewidth,height=2*height)
    fig  
  }
  # Sharks
  pp <- p
  pp$w0 <- 0.0027*p$W
  NN <- spectrum(pp)
  fig <- fig +
    geom_line(data=NN, aes(x=w, y=NprR/NprR[1] * w/pp$w0 * g/g[1]), size=1, color="red")
  fig <- loglog(fig, ylim=c(1,1e6), xlim=c(1e-3, max(W)))
  ggsave(paste("Chapter4/cohort Presentation 4.pdf"),width=doublewidth,height=2*height)
  
}

plotR0 <- function() {
  p <- baseparameters()
  W <- 10^seq(-1,5,length.out=10)
  R0full <- 0*W
  R0vonB <- 0*W
  R0simple <- 0*W
  
  for (i in 1:length(W)) {
    p$W <- W[i]
    
    N <- spectrum(p)
    R0full[i] <- N$R0[1]
    
    N <- spectrumana(p)
    R0vonB[i] <- N$R0[1]
    
    N <- spectrumsimple(p)
    R0simple[i] <- N$R0[1]
  }
  
  fig <- ggplot() +
      geom_line(aes(x=W, y=R0full), size=1) +
      geom_line(aes(x=W, y=R0vonB), size=0.5) +
      geom_line(aes(x=W, y=R0simple), size=0.5, linetype="dashed") +
      geom_hline(yintercept=1, size=thin, linetype="dotted") +
      xlab(TeX("Asymptotic weight, $\\textit{W}_{\\infty}$ (g)")) + 
      ylab(TeX("Eggs per recruit, $\\textit{R}_0$"))
  fig <- loglog(fig, ylim=c(0.5,1000))

  ggsave("Chapter4/R0.pdf",width=singlewidth,height=height)
  fig  
}

plotBHpanel <- function() 
{
  Rp = 10^seq(-3,4,length.out=100)  # Rp/Rmax
  fig <- ggplot()
  # Plot grey-scale background
  #fig <- fig +
  #  geom_ribbon(aes(x=Rp, ymin=0.01*Rp, ymax=1*Rp), fill=grey(0.8)) +
  #  geom_ribbon(aes(x=Rp, ymin=0.1*Rp, ymax=1*Rp), fill=grey(0.6))
  
  # Plot BH curve
  fig <- fig +
    geom_line(aes(x=Rp, y=1-(1+Rp)^(-1)), size=thick) +
    geom_line(aes(x=Rp, y=Rp), linetype="dashed", color="black", size=thin) +
    geom_hline(yintercept=1, linetype="dotted", color="black", size=thin)
  
  # Add points for species
  W <- c(10,100,1000,10000)
  p <- baseparameters()
  for (i in 1:length(W))
  {
    p$W <- W[i]
    N <- spectrum(p)
    fig <- fig +
      geom_point(data=N, aes(x=Rp[1], y=R[1]), size=i+1, colour="white") +
      geom_point(data=N, aes(x=Rp[1], y=R[1]), size=i)
  }
  
  fig <- fig +
    scale_y_continuous(limits=c(0,1.1), breaks=seq(0,1,by = 0.2), oob=rescale_none)+
    scale_x_continuous(limits=c(0,20),oob=rescale_none)+
    #    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
    #                  labels = trans_format("log10", math_format(10^.x)),oob=rescale_none,
    #                  limits=c(0.1,1000)) +
    
    xlab(TeX("Egg production, $\\textit{R}_{\\textit{p}}/\\textit{R}_{max}$")) + 
    ylab(TeX("Recruitment, $\\textit{R}/\\textit{R}_{max}$")) 
  fig <- mytheme(fig)
  
  fig
}

plotBH <- function()
{
  fig1 <- plotBHpanel()
  fig2 <- plotBHpanel()
  fig2 <- semilogx(fig2, xlim=c(0.1,1000), ylim=c(0,1.1)) +
    scale_y_continuous(breaks=seq(0,1,by = 0.2), limits = c(0,1.1))
#    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
#                  labels = trans_format("log10", math_format(10^.x)),oob=rescale_none,
#                  limits=c(0.1,1000)) +
  #pdf("BH.pdf", width=doublewidth, height=height)
  #multiplot(fig1,fig2,cols=2)
  #dev.off()
  grid <- plot_grid(fig1,fig2,align="h")
  ggsave("Chapter4/BH.pdf", grid, width=doublewidth, height=height)
  grid
}

plotBHvsSSB <- function()
{
  SSB <- seq(0,.2,length.out=100)
  W = c(1,100,10000)
  fig <- ggplot()
  p <- baseparameters()
  
  for (i in 1:length(W))
  {
    p$W <- W[i]
    N <- spectrum(p)
    kr <- calcKr(p)
    Rp <- p$epsR * kr* SSB/p$w0  # Rp/R
    R <- 1-(1+Rp)^(-1) # R/Rmax
    dat = data.frame(SSB=SSB, R=R, Rp=p$epsR*kr*SSB/p$w0 )
    fig <- fig +
      geom_line(data=dat, aes(x=SSB, y=R), size=(i-1)/3+1) +
      geom_line(data=dat, aes(x=SSB, y=Rp ), linetype="dashed", size=0.5) +
      geom_line(data=dat, aes(x=SSB, y=SSB/SSB), linetype="dotted", size=0.3)
  }
  fig <- fig +
    scale_y_continuous(limits=c(0,1.1),breaks=seq(0,1,by = 0.2),oob=rescale_none) +
    scale_x_continuous(oob=rescale_none) +
    xlab(TeX("Spawning stock biomass, $B_{SSB}/R_{max}$ (g yr/#)")) +
    ylab(TeX("Recruitment, $R/R_{max}$"))
  fig <- mytheme(fig)

  ggsave("BHvsSSB.pdf",width=singlewidth,height=height)
  fig
}

plotSteepness <- function()
{
  W <- 10^seq(0,5,length.out=100)
  steepness <- 0*W # steepness
  steepness_small <- 0*W
  
  p <- baseparameters()
  for (i in 1:length(W)) 
  {
    p$W <- W[i]
    N <- spectrum(p)
    steepness[i] <- (1 + N$RpperR) / (5 + N$RpperR)
    steepness_small[i] <- (1 + 0.1*N$RpperR) / (5 + 0.1*N$RpperR)
  }
  
  fig <- ggplot() +
    geom_ribbon(aes(x=W, ymin=steepness_small, ymax=steepness), fill=grey(0.6)) +
    geom_line(aes(x=W, y=steepness_small), size=1.5) +
    scale_y_continuous(limits=c(0,1.1),breaks=seq(0,1,by = 0.2),oob=rescale_none) +
    #scale_x_log10(breaks = bks,
    #              labels = trans_format("log10", math_format(10^.x)),oob=rescale_none) +
    xlab("Asymptotic size, Winf (g)") +
    ylab("Steepness") 
  fig <- semilogx(fig)

  fig
}

plotRecruitment <- function()
{
  W <- c(.1,.8)
  fig <- ggplot()
  p <- baseparameters()
  
  for (i in 1:length(W))
  {
    p$W <- W[i]
    N <- spectrum(p)
    ix <- N$w >= W[i]*p$etaM
    n0 <- N$NprR[1]
    fig <- fig + 
      geom_ribbon(data=N[ix,], aes(x=w, ymin=0*w, ymax=(w/p$w0)^2*NprR/n0), fill=grey(0.6)) +
      geom_line(data=N, aes(x=w, y=(w/p$w0)^2*NprR/n0),size=1.3)
  }
  
  fig <- fig +
    xlab("Weight (g)") + 
    ylab("B/R") +
    annotation_logticks(long=unit(1,"mm"),mid=unit(0,"cm"), short=unit(0,"cm"))  # l = left, b = bottom etc) 
    
  fig <- loglog(fig)
  ggsave("recruitment.pdf", width=8, height=6, units="cm")
  fig
}

plotAlpha <- function()
{
  p <- baseparameters()
  W <- 10^seq(1,5,length.out = 100)
  
  alpha <- read.csv("Chapter4/Hall et al alpha data.csv", header=TRUE)
  alpha$W <- weight(exp(alpha$log.Linf))
  alpha$alpha <- exp(alpha$log.alpha.)
  Linf <- weight2length(W)
  alphafit <- exp(11-2.3*log(Linf)) 
  
  wAgeOne <- ( p$A*(1-p$n)+p$w0^(1-p$n) )^(1/(1-p$n))
  P <- (wAgeOne/p$w0)^-p$a
  
  fig <- ggplot() +
    geom_point(data=alpha, aes(x=W, y=alpha)) +
    geom_line(aes(x=W, y=alphafit)) +
    geom_line(aes(x=W, y=P * p$epsR * p$epsEgg * p$A*W^(p$n-1)/p$w0),linetype="dashed") +
    xlab(TeX("Asymptotic size, $W_\\infty$ (g)")) +
    ylab(TeX("Recruiment $\\alpha$ (#/yr/g)"))
  fig <- loglog(fig, ylim=c(0.01,2000))

  ggsave("Chapter4/alpha.pdf", width=singlewidth, height=height)
  fig
}

plot_acrit <- function()
{
  R0minusone <- function(a,p, specfunc) {
    p$a <- a
    N <- specfunc(p)
    N$R0[1]-1
  }
  
  p <- baseparameters()
  W <- 10^seq(0.5, 7,length.out = 30)
  acrit <- 0*W
  acritvonB <- 0*W
  acritsimple <- 0*W
  for (i in 1:length(W))
  {
    p$W <- W[i]
    acrit[i] = uniroot(R0minusone, c(0,1), p, spectrum)$root
    acritvonB[i] = uniroot(R0minusone, c(0,1), p, spectrumana)$root
    acritsimple[i] = uniroot(R0minusone, c(0,1), p, spectrumsimple)$root
  }
  # Get a data from Gislason:
  A <- read.csv("Data/Gislason.csv",sep=",", header=FALSE)
  names(A) <- c("L","Linf","K","M","tau")
  A <- A[A$L > 0.5*A$Linf, ]  # Selection only adults
  A$W <- p$q*A$L^3
  p <- baseparameters()
  A$a = A$M/A$K*p$aMKconstant
  
  fig <- ggplot() +
    geom_ribbon(aes(x=W, ymin=acrit, ymax=acrit/acrit), fill=grey(0.6)) +
    geom_line(aes(x=W, y=acrit), size=thick) +
    geom_line(aes(x=W, y=acritvonB), size=thin) +
    geom_line(aes(x=W, y=acritsimple), size=thin, linetype="dashed") +
    geom_point(aes(x=A$W, y=A$a)) +
    annotate("text", x=20, y=0.75, label="R[0] < 1", parse=TRUE, colour="white") +
    annotate("text", x=20, y=0.1, label="R[0] > 1", parse=TRUE, colour="black") +
    #scale_y_continuous(limits = c(0, 1)) +
    xlab(TeX("Asymptotic size, $\\textit{W}_\\infty$ (g)")) +
    ylab(TeX("Phys. mortality, $\\textit{a}$"))
  fig <- fig +
    scale_x_log10(breaks = bks, #10^seq(floor(xrange[1]), floor(xrange[2]), by=1), #trans_breaks("log10", function(x) 10^x),
                           labels = trans_format("log10", math_format(10^.x)),
                           limits = c(min(W), max(W)), oob=rescale_none,
                  expand=c(0,0)) +
    annotation_logticks(long=unit(1,"mm"),mid=unit(0.5,"mm"), 
                        short=unit(0.5,"mm"), size=0.25, sides="b") +
    scale_y_continuous(oob=rescale_none, expand=c(0,0), limits=c(0,1))
  fig <- mytheme(fig)
  
  ggsave("Chapter4/acrit.pdf", width=singlewidth, height=height)
  fig
}

plotGislasonMK <- function()
{
  library(fitdistrplus)
  A <- read.csv("Data/Gislason.csv",sep=",", header=FALSE)
  names(A) <- c("L","Linf","K","M","tau")
  A <- A[A$L > 0.5*A$Linf, ]  # Selection only adults
  p <- baseparameters()
  A$a = A$M/A$K*p$aMKconstant
  A <- A[A$a < 1.25,]  # remove outliers
  
  amean <- mean(A$a)
  ageommean <- exp(mean(log(A$a)))
  cat("Mean value of a: ",amean ,"\n")
  #
  # Fit log-normal
  #
  mle = fitdistr(A$a, "lognormal")
  cat("Log-normal fit, meanlog: ", mle$estimate["meanlog"], " sdlog:", mle$estimate["sdlog"],"\n")
  
  fig <- ggplot(data=A, aes(x=a, y=..density..)) + 
    geom_histogram(fill=grey(0.6), size=0.2, bins=20) +
    geom_density(adjust=2, colour=NA) +
    #geom_line(stat="density", adjust=2, size=1.5) +
    geom_line(aes(x=a, y=dlnorm(x=a, meanlog=mle$estimate["meanlog"], sdlog=mle$estimate["sdlog"])),
              size=thick)+
    geom_vline(xintercept = amean, linetype="dotted", size=thin) + 
    #geom_vline(xintercept = ageommean, linetype="dotted") + 
    xlab(TeX("Physiological mortality, $\\textit{a}$")) +
    ylab("Density of observations") +
    scale_y_continuous(oob=rescale_none)+
    scale_x_continuous(limit=c(0,1.2), oob=rescale_none)
  fig <- mytheme(fig)

  ggsave("Chapter4/a.pdf", width=singlewidth, height=height)
  fig
  
}

testHerring <- function() {
  dat <- read.table("Chapter4/Herring_ICES30.csv", header=TRUE, sep=";")
  defaultplot()
  defaultpanel(xlim=c(0,1e6), ylim=dat$R)
  points(dat$SSB, dat$R, type="b")
  Rmax = 1e10
  alpha=100000
  lines(dat$SSB, Rmax*alpha*dat$SSB/(alpha*dat$SSB+Rmax))
  out = summary( nls(R ~ Rmax*(alpha*SSB/(alpha*SSB+Rmax)), data=dat, start=list(Rmax=Rmax, alpha=alpha)) )
  SSB = seq(0,6e6, length.out = 100)
  Rmax = out$coefficients[1]
  alpha = out$coefficients[2]
  lines(SSB, Rmax*alpha*SSB/(alpha*SSB+Rmax), col="red")
}

plotAllChapter4 <- function()
{
  plotSpectrum()
  plotCohort()
  plotR0()
  plotBH()
  #plotAlpha()
  plotGislasonMK()
  plot_acrit()
}


