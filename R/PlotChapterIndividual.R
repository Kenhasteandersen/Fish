#
# TeX/Chapter 3
#
library(fishsizespectrum)
source("R/basetools.R")

dir.create("TeX/ChapterIndividual")

getGrowthParameters <- function()
{
  data <- read.csv("Data/fecundpruned.csv",header=TRUE,sep=",")
  data$Winf <- data$aW * data$Linf^data$bW
  data$Wm  <- data$aW * data$Lm^data$bW
  data$Phylum <- data$phyl
  levels(data$Phylum) <- c("Elasmobranchs", "Teleosts")
  #
  # Data from Kooijman and Gislason:
  #
  data2 = read.csv("Data/KooijmanGislason.csv", header=FALSE, col.names=c("K","Linf","T"))
  #
  # Combine data for fish
  # 
  Q10 <- 1.83
  ix  <- data$Phylum == "Teleosts"
  dat <- data.frame(K=data$K[ix] * Q10^((15-data$T[ix])/10), 
                    W=data$Winf[ix])
  dat2 <- data.frame(K=data2$K * Q10^((15-data2$T)/10), W=0.01*data2$Linf^3)
  dat = rbind(dat, dat2)
  
  p <- baseparameters()
  dat$L <- (dat$W/p$c)^(1/3)
  
  dat$A = calcA(dat$K, dat$L)
  
  dat
}

panelvonB <- function()
{
  #
  # Sketch of von B curve
  # 
  Linf <- 100
  K <- 0.18
  
  ages = seq(0,30,length.out=100)
  
  defaultpanel(xlim=ages, ylim=c(0,c(0,110)),
               xlab="Age (years)",
               ylab="Length (cm)",
               label=TRUE)
  lines(ages, Linf*(1-exp(-K*ages)), lwd=3)
  lines(ages, y=K*Linf*ages, lty="dashed", lwd=1)
  lines(x=ages, y=Linf*ages/ages, lty="dashed", lwd=1)
}

panelSketch <- function()
{
  #
  # Sketch of Winf
  # 
  xlim = c(1e-11,1e5)
  ylim = c(5e-10,1e6)
  p <- baseparameters()
  
  n <- p$n
  W1 = 10
  W2 = 10000
  A <- p$A
  k = function(W) A*W^(n-1)
  
  w=10^seq(from=-3, to=6, length.out=100)
  Ea = A

  loglogpanel(xlim=w, ylim=A*w^n, 
              xlab="Body weight (g)", 
              ylab="Available energy, $E_a$ (g/yr)",
              label=TRUE)
  lines(w, A*w^n, lwd=3)
  
  ix = w<10*W1
  lines(w[ix], k(W1)*w[ix])
  lines(W1*c(1,1), c(1e-4,6*W1^n), lty=dashed)

  ix = w<10*W2
  lines(w[ix], k(W2)*w[ix])
  lines(W2*c(1,1), c(1e-4,6*W2^n), lty=dashed)
}


plotSketch_and_vonB <- function()
{
  defaultplot(mfcol=c(1,2))
  panelvonB()
  panelSketch()
}


plotKvsLinf <- function()
{
  dat <- getGrowthParameters()
  
  fig <- ggplot() +
    geom_point(data=dat, aes(x=L, y=K), colour=grey(0.3))
  # fit
  fit <- lm(log(K) ~ log(L), data=dat)
  tmp <- data.frame(x=dat$L, y=exp(predict(fit)))
  cat('Fit of K: C = ', exp(fit$coefficients[1]), 
      ' and exponent ', fit$coefficients[2], "\n")
  
  # fit with exponent -0.75:
  C <- exp(mean(log(dat$K * dat$L^0.75)))
  cat('K fitted wit n=0.75: ', C, "\n")
  
  # fit with exponent -1:
  C2 <- exp(mean(log(dat$K * dat$L)))
  
  fig <- fig +
    geom_line(data=tmp, aes(x=x,y=y),linetype="dashed") +
    geom_line(data=tmp, aes(x=x, y=C*x^-0.75)) #+
  #geom_line(data=tmp, aes(x=x, y=C2/x), linetype="dotted")
  
  fig <- fig +
    xlab(TeX("Asymptotic length, $\\textit{L}_{\\infty}$ (cm)")) + 
    ylab(TeX("von Bertalanffy coefficient, $\\textit{K}$ (1/yr)"))
  fig <- loglog(fig, ylim=c(0.03, 3))
  
  ggsave("TeX/ChapterIndividual/KvsLinf.pdf",width=singlewidth,height=singlewidth)
  fig
}

plotA <- function() 
{
  dat <- getGrowthParameters()
  
  fig <- ggplot() +
    geom_point(data=dat, aes(x=W, y=A), colour=grey(0.3))
  
  # fit
  fit <- lm(log(A) ~ log(W), data=dat)
  tmp <- data.frame(x=dat$W, y=exp(predict(fit)))
  fig <- fig + geom_line(data=tmp, aes(x=x,y=y),linetype="dashed")
  cat("Fit to A, coef.: ", exp(fit$coefficients[1]), ", exponent: ", fit$coefficients[2],'\n')
  
  # Mean value:
  Amean <- exp(mean(log(dat$A)))
  cat("Amean ",Amean,"\n")
  fig <- fig + 
    geom_line(data=data.frame(x=c(min(dat$W), max(dat$W))), 
              aes(x=x, y=Amean)) #+
  #
  # Finalize figure
  #
  fig <- fig +
    theme(legend.position=c(1,0), legend.justification=c(1,0)) +
    scale_colour_grey()+ 
    xlab(TeX("Asymptotic weight, $\\textit{W}_{\\infty}$ (g)")) + 
    ylab(TeX("Growth coef. $\\textit{A}$ ($g^{1/4}$/yr)")) +
    labs(colour="",shape="")
  fig <- loglog(fig)  
  
  ggsave("TeX/ChapterIndividual/A.pdf", width=singlewidth, height=height)
  fig
}


panelEtam <- function()
{
  p <- baseparameters()
  # 
  # etam
  # 
  data <- read.csv("Data/fecundpruned.csv",header=TRUE,sep=",")
  data$Winf <- data$aW * data$Linf^data$bW
  data$Wm  <- data$aW * data$Lm^data$bW
  data$etam <- data$Wm/data$Winf
  data$Phylum <- data$phyl
  levels(data$Phylum) <- c("Elasmobranchs", "Teleosts")
  
  loglogpanel(xlim=data$Winf, ylim=c(0.01,2),
              xlab="Asymptotic weight $\\textit{W}_{\\infty}$ (g)", 
              ylab="Rel. maturation size, $\\textit{\\eta}_\\textit{m}$",
              label=TRUE)
  ixElasmo = data$Phylum=="Elasmobranchs"
  ixTele = data$Phylum=="Teleosts"
  points(data$Winf[ixElasmo], data$etam[ixElasmo], pch=dots, col=stdgrey)
  points(data$Winf[ixTele], data$etam[ixTele], pch=dots)
  #
  # Fits
  # 
  
  # all:
  fit <- lm(log(etam) ~ log(Winf), data=data)
  #tmp <- data.frame(x=data$Winf, y=exp(predict(fit)))
  #fig <- fig + geom_line(data=tmp, aes(x=x,y=y),linetype="dashed")
  lines(range(data$Winf), exp(fit$coefficients[1])*range(data$Winf)^fit$coefficients[2],
        lty=dashed, lwd=3)  
  # teleosts:
  ix <- data$Phylum=="Teleosts"
  fitfish <- lm(log(Wm[ix]/Winf[ix]) ~ log(Winf[ix]), data=data)
  tmp <- data.frame(x=data$Winf[ix], y=exp(predict(fitfish)))
  #fig <- fig + geom_line(data=tmp, aes(x=x,y=y), linetype="dotted")
  # teleosts with exponent 1:
  etam <- exp(mean(log(data$Wm/data$Winf)))
  cat('Average eta_m',etam,'\n')
  lines(range(data$Winf), etam*c(1,1), lwd=3)
  
  legend("bottomright", bty="n", 
         legend=c("Teleosts", "Elasmobranchs"),
         pch=c(dots,dots),
         col=c("black", stdgrey))
}

panelOgive <- function() {
  ogive <- as.matrix(read.table(file="Data/maturity_saithe.dat", header=FALSE, skip=5))
  weights <- as.matrix(read.table(file="Data/weightatage_saithe.dat", header=FALSE, skip=5))
  
  weights[weights==0] <- max(weights)
  nages <- dim(weights)[1]
  nyears <- dim(weights)[2]
  dat <- data.frame(weight=log(as.vector(weights)), ogive=as.vector(ogive), 
                    year=rep(seq(1,nyears),nages))
  dat$ogive[ogive==0] <- NaN
  dat$ogive[ogive==1] <- NaN
  dat$ogive.logit <- qlogis(dat$ogive)
  
  dat$year = factor(dat$year)
  
  m = glm( ogive.logit ~ 0 + year + weight, data=dat )
  
  semilogxpanel(xlim=c(0.1,10), ylim=c(0.01,1),
                xlab="\\textit{w}/\\textit{w}_{maturation}", 
                ylab="Maturity",
                label=TRUE)
  
  for(i in 1:nyears) {
    ix = dat$year==i
    points(exp(dat$weight[ix]+coef(m)[i]/coef(m)[nyears+1]), dat$ogive[ix], pch=dots, col=stdgrey)
  }
  w = 10^seq(-2,1, length.out = 100)
  lines( w, psi(w, coef(m)[nyears+1]) , lwd=3)  
  vline(1)
  
  print(coef(m)[nyears+1])
}

plotMaturity <- function() {
  defaultplot(mfcol=c(1,2))
  panelEtam()
  panelOgive()
  
}

plotR <- function()
{
  #
  # Use gunderson data
  #
  p <- baseparameters()
  A <- read.csv("Data/Gunderson.csv", header = TRUE, sep = ",")
  A$Winf <- weight(A$Linf)
  
  # Estimate A from age at maturation + 0.5:
  A$A <- p$etaM^0.25/((1-p$n)*(A$alpha+0.5))*A$Winf^(1-p$n) #log(1-p$etaM^0.25*p$etaA) / (-0.25*p$etaA*(A$alpha+0.5) ) * A$Winf^(1-p$n)
  
  fig <- ggplot(data=A) +
    geom_point(aes(x=Winf, y=GSI/A), colour=grey(0.3))
  
  # all:
  fit <- lm(log(GSI/A) ~ log(Winf), data=A)
  tmp <- data.frame(x=A$Winf, y=exp(predict(fit)))
  fig <- fig + geom_line(data=tmp, aes(x=x,y=y),linetype="dashed")
  cat('Fit to R', fit$coefficients,'\n')
  
  epsA <- exp(mean(log(A$GSI/A$A*A$Winf^0.25)))
  cat('eps_repro = ', epsA, '.\n')
  
  fig <- fig + 
    geom_line(aes(x=Winf, y=Winf^(p$n-1))) +
    geom_line(aes(x=Winf, y=Winf^(p$n-1)*epsA)) +
    xlab(TeX("Asymptotic weight, $\\textit{W}_{\\infty}$ (g)")) + 
    ylab(TeX("Repro. output  ($g^{-0.25}$)"))
  
  fig <- loglog(fig)
  ggsave("TeX/ChapterIndividual/R.pdf", width=singlewidth, height=height)
  fig
}  


plotGrowth <- function()
{
  #
  # Plot some growth functions
  #
  library(deSolve)
  width = 10
  height = 6
  
  W = 2000
  fig <- ggplot()
  p <- baseparameters()  
  p$W <- W
  # Add my growth curve:
  #n <- p$n
  #epsA <- p$etaA
  etaM <- p$etaM
  A <- p$A
  
  # Interval:
  CV = 1.95 # Estimated as exp(sqrt(var(log(dat$A)))) in plotA.R
  ages = seq(0,50,length.out=500)
  p1 <- p
  p1$W <- p$W*0.8
  p1$A <- p$A/CV
  p2 <- p
  p2$W <- p$W/0.8
  p2$A <- p2$A*CV
  out1 <- ode(y=p$w0, times=ages, func=function(t,w,p) list(growth(p,w)), parms=p1)
  out2 <- ode(y=p$w0, times=ages, func=function(t,w,p) list(growth(p,w)), parms=p2)
  
  fig <- fig + 
    geom_ribbon(aes(x=out1[,1], ymin=out1[,2], ymax=out2[,2]),
                fill=grey(0.5))
  
  # Add observed von Bertalanffy curves:
  dat <- getGrowthParameters()
  
  ix <- dat$W>0.8*W & dat$W<(1/0.8)*W
  
  ages <- seq(0,50, length.out = 150)
  for (i in 1:sum(ix)) {
    w <- weight( (dat[ix,]$L[i]*(1-exp(-dat[ix,]$K[i]*ages))) )
    fig <- fig + 
      geom_line(data=data.frame(x=ages,y=w),
                aes(x=x,y=y), size=thin, colour=grey(0.75))
  }
  # "Median" line
  out <- ode(y=p$w0, times=ages, func=function(ages,w,p) list(growth(p,w)), parms=p)
  fig <- fig + 
    geom_line(aes(x=out[,1],y=out[,2]), size=1) +
    geom_hline(yintercept=etaM*W, linetype="dotted",size=thin) +
    geom_vline(xintercept = ageMaturation(W=W), linetype="dotted",size=thin)
  
  fig <- fig +
    xlab("Age (yr)") + 
    ylab("Weight (g)") +
    scale_x_continuous(limits=c(0,0.5*A*W^0.25), oob=rescale_none) +
    scale_y_continuous(oob=rescale_none)
  fig <- mytheme(fig)
  
  ggsave("TeX/ChapterIndividual/Growth.pdf",width=width, height=height, units="cm")
  fig
}


plotBudget <- function()
{
  p <- baseparameters()
  W <- 1
  eta <- 0.001
  A <- p$A
  h <- p$h
  n <- p$n
  epsEgg <- p$epsEgg
  etaM <- p$etaM
  
  k  = A*W^(n-1)
  maxrate = A*W^n
  dat = data.frame(w = 10^seq(from=log10(eta*W), to=log10(W), length.out=100))
  dat$psi = (1 + (dat$w/(etaM*W))^(-10))^(-1)
  
  fig <- ggplot() +
    # consumption:
    geom_line(data=dat, aes(x=w, y=p$f0*h*w^n/maxrate), size=thick) +
    # consumption minus egestion and excretion:
    geom_ribbon(data=dat, aes(x=w, ymax=(1-0.15-0.1)*p$f0*h*w^n/maxrate, 
                              ymin=0.01/maxrate), fill=grey(0.9)) +
    # available energy
    geom_ribbon(data=dat, aes(x=w, ymax=A*w^n/maxrate, ymin=0.01/maxrate), 
                fill=grey(0.8)) +
    geom_line(data=dat, aes(x=w, y=A*w^n/maxrate)) +
    #geom_ribbon(data=dat, aes(x=w, ymax=(A*w^n-epsR*k*w)/maxrate, 
    #                          ymin=0.01/maxrate), fill=grey(0.6)) +
    # Reproductive out:
    geom_ribbon(data=dat, 
                aes(x=w, ymax=(A*w^n-(1-epsEgg)*psi*k*w)/maxrate, ymin=0.01/maxrate), 
                fill=grey(0.6)) +
    # Reproductive output:
    #geom_ribbon(data=dat, 
    #            aes(x=w, ymax=(A*w^n-epsR*psi*k*w)/maxrate, ymin=0.01/maxrate), 
    #            fill=grey(0.4)) +
    # Growth
    #geom_line(data=dat, aes(x=w, y=(A*w^n-psi*k*w)/maxrate), linetype="dashed") +
    geom_ribbon(data=dat, 
                aes(x=w, ymax=(A*w^n-psi*k*w)/maxrate, ymin=0.01/maxrate), 
                fill=grey(0.2)) +
    #geom_line(data=dat, aes(x=w, y=(A*w^n-(1-epsEgg)*psi*k*w)/maxrate), 
    #          linetype="dashed", size=thick) +
    geom_vline(xintercept = etaM, linetype="dotted", size=thin) +
    xlab(TeX("Relative body weight $(\\textit{w}/\\textit{W}_{\\infty})$")) + 
    ylab("Normalized rates")
  fig <- loglog(fig, xlim=c(0.01,1), ylim=c(2e-2,3)) 
  
  # W <- 1
  # fig2 <- ggplot() +
  #   geom_line(data=dat, aes(x=w, 
  #       y=((0.15*p$f0 + p$fc*p$epsA)*p$h*w^p$n + psi*(1-p$epsEgg)*A*W^(p$n-1)*w) / w^p$n)) +
  #   geom_line(data=dat, aes(x=w, A*W^(p$n-1)*w / w^p$n), size=0.5)
  # fig2 <- semilogx(fig2)
  # fig2
  
  ggsave("TeX/ChapterIndividual/budget.pdf",width=singlewidth, height=height)
  fig
}

plotTestGrowth <- function(K, Linf)
{
  p <- baseparameters()
  A <- calcA(K,Linf)
  p$A <- A
  W <- weight(Linf)
  p$W <- W
  t <- seq(0,30,length.out = 100)
  cat(A,'\n')
  
  out <- ode(y=p$w0, times=t, func=function(t,w,p) list(growth(p,w)), parms=p)
  
  fig <- ggplot() +
    geom_line(aes(x=t, y=weight(Linf*(1-exp(-K*t))))) +
    geom_line(aes(x=t, y=(0.25*p$A*t + p$w0^0.25)^4)) +
    geom_hline(yintercept = p$etaM*W) +
    geom_line(aes(x=out[,1], y=out[,2])) +
    ylim(c(0, W))
  fig
}


plotAllChapterIndividual <- function()
{
  pdfplot("TeX/ChapterIndividual/vonB2.pdf", plotSketch_and_vonB, width=doublewidth, height=height)
  plotKvsLinf()
  plotA()
  pdfplot("TeX/ChapterIndividual/maturity.pdf", plotMaturity, width=doublewidth, height=height)
  plotR()
  plotGrowth()
  plotBudget()
  
}
