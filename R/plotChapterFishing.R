#
# TeX/Chapter 5
#
library(fishsizespectrum)
source("R/basetools.R")

dir.create("TeX/ChapterFishing")

plotSelectivity <- function()
{
  getData <- function(filename)
  {
    data <- read.csv(filename, header=FALSE)
    zeros <- 0*(1:(dim(data)[1]/3))
    dat =data.frame(L=zeros, min=zeros, sel=zeros, max=zeros)
    for (i in 1:(dim(data)[1]/3))
    {
      ix <- (i-1)*3+1
      dat$L[i] <- mean(data[ix:ix+2,1])
      dat$min[i] <- data[ix,2]
      dat$sel[i] <- data[ix+1,2]
      dat$max[i] <- data[ix+2,2]
    }
    dat
  }
  
  LL <- seq(20,130,length.out = 100)
  breaks <- c(30,40,50,60,70,80,90,100,120)
  
  # Gillnet
  dat <- getData("Data/Selectivity_gillnet.csv")
  fit1 <- nls(sel ~ exp(-(log(L/L0))^2/sigma^2), data=dat, start=list(L0=60,sigma=.1))
  fig1 <- ggplot() +
    geom_line(data=data.frame(L=LL, y=predict(fit1,list(L=LL))), aes(x=L, y=y), size=thick, color=grey(0.4)) +
    geom_point(data=dat,aes(x=L, y=sel)) +
    geom_errorbar(data=dat,aes(x=L, ymin=min, ymax=max), width=0) +
    geom_hline(yintercept = 1, linetype="dotted", size=thin) +
    scale_x_log10(breaks=breaks, labels=breaks, 
                  limit=c(30,120)) +
    scale_y_continuous(limit=c(0,1.2), oob = rescale_none) +
    ylab("Selectivity") +
    xlab("") +
    annotate("text", x=30, y=1.12,label="Gill net", hjust=0)
  fig1 <- mytheme(fig1)
  
  # Trawl:
  dat <- getData("Data/Selectivity_trawl.csv")
  fit2 <- nls(sel ~ 1/(1+(L/L0)^-u), data=dat, start=list(L0=30,u=4))
  fig2 <- ggplot(data=dat) +
    geom_line(data=data.frame(L=LL, y=predict(fit2,list(L=LL))), aes(x=L, y=y), size=thick, color=grey(0.4)) +
    geom_point(aes(x=L, y=sel)) +
    geom_errorbar(aes(x=L, ymin=min, ymax=max), width=0) +
    geom_hline(yintercept = 1, linetype="dotted", size=thin) +
    scale_x_log10(breaks=breaks, labels=breaks, 
                  limit=c(30,120)) +
    scale_y_continuous(limit=c(0,1.2), oob = rescale_none) +
    xlab("Length (cm)") +
    ylab("Selectivity")+
    annotate("text", x=30, y=1.12,label="Otter trawl", hjust=0)
  fig2 <- mytheme(fig2)
  
  # Longline:
  dat <- getData("Data/Selectivity_longline.csv")
  fit3 <- nls(sel ~ 1/(1+(L/L0)^-u), data=dat, start=list(L0=70,u=1))
  fig3 <- ggplot(data=dat) +
    geom_line(data=data.frame(L=LL, y=predict(fit3,list(L=LL))), aes(x=L, y=y), size=thick, color=grey(0.4)) +
    geom_point(aes(x=L, y=sel)) +
    geom_errorbar(aes(x=L, ymin=min, ymax=max), width=0) +
    geom_hline(yintercept = 1, linetype="dotted", size=thin) +
    scale_x_log10(breaks=breaks, labels=breaks, 
                  limit=c(30,120)) +
    scale_y_continuous(limit=c(0,1.2), oob = rescale_none) +
    xlab("") +
    ylab("") +
    annotate("text", x=30, y=1.12,label="Long line", hjust=0)
  fig3 <- mytheme(fig3)
  
  # Trap:
  dat <- getData("Data/Selectivity_trap.csv")
  fit4 <- nls(sel ~ 1/(1+(L/L0)^-u), data=dat, start=list(L0=70,u=1))
  fig4 <- ggplot(data=dat) +
    geom_line(data=data.frame(L=LL, y=predict(fit4,list(L=LL))), aes(x=L, y=y), size=thick, color=grey(0.4)) +
    geom_point(aes(x=L, y=sel)) +
    geom_errorbar(aes(x=L, ymin=min, ymax=max), width=0) +
    geom_hline(yintercept = 1, linetype="dotted", size=thin) +
    scale_x_log10(breaks=breaks, labels=breaks, 
                  limit=c(30,120)) +
    scale_y_continuous(limit=c(0,1.2), oob = rescale_none) +
    xlab("Length (cm)") +
    ylab("") +
    annotate("text", x=30, y=1.12,label="Trap", hjust=0)
  fig4 <- mytheme(fig4)
  
  pdf("TeX/ChapterFishing/selectivity.pdf", width=doublewidth, height=2*height, useDingbats=FALSE)
  multiplot(fig1,fig2,fig3,fig4, cols=2)
  dev.off()
}


plotTrawl <- function()
{
  p <- baseparameters()
  
  etaF <- p$etaF
  F <- 0.3
  fig <- ggplot(data=data.frame(w=10^seq(-7,1,length.out = 100))) +
    geom_line(aes(x=100000*w, y=p$a*p$A*(100000*w)^(p$n-1)), color=grey(0.5), size=thick) +
    geom_line(aes(x=10*w, y=F*psi(w/etaF,3)),size=thin) +
    geom_line(aes(x=333*w, y=F*psi(w/etaF,3)),size=0.5*(thin+thick)) +
    geom_line(aes(x=10000*w, y=F*psi(w/etaF,3)),size=thick) +
    geom_line(aes(x=p$etaM*10*w/w, y=seq(-1,0.35,length.out=100)),linetype='dotted',size=thin) +
    geom_line(aes(x=p$etaM*333*w/w, y=seq(-1,0.35,length.out=100)),linetype='dotted',size=thin) +
    geom_line(aes(x=p$etaM*10000*w/w, y=seq(-1,0.35,length.out=100)),linetype='dotted',size=thin) +
    xlab(TeX("Body weight (g)")) +
    ylab(TeX("Fishing mort., $\\mu_{\\textit{F}}$ (yr$^{-1}$)")) 
  fig <- semilogx(fig, xlim=c(1e-1,1e5), ylim=c(0,1))
  #fig <- mytheme(fig)
  ggsave("TeX/ChapterFishing/trawl.pdf",width=singlewidth, height=height)
  fig
}

# plotGillnet <- function()
# {
#   p <- baseparameters()
#   p$W <- 10000
#   p$F <- 0.3
#   p2 <- p
#   p2$etaF <- 0.1*p$etaF
#   fig <- ggplot(data=data.frame(w=10^seq(-3,log10(p$W) ,length.out = 100))) +
#     geom_line(aes(x=w, y=fishingGillnet(w,p))) +
#     geom_line(aes(x=w, y=fishingGillnet(w,p2))) +
#     xlab(TeX("Weight, w (g)")) +
#     ylab(TeX("Fishing mortality, $\\mu_F$ (yr$^{-1}$)")) +
#     scale_y_continuous(limit=c(0,1), oob=rescale_none)
#   fig <- semilogx(fig, xlim=c(1e-1,p$W))
#   
#   #ggsave("trawl.pdf",width=singlewidth, height=height)
#   fig
# }

# plotSpectrum<-function(p)
# {
#   spec <- spectrum(p)
#   fig <- ggplot(data=spec, aes(x=w)) +
#     geom_line(aes(y=NprR*R*w))
#   cat("yield: ", calcYield(p,spec),"\n")
#   fig <- loglog(fig, ylim=c(1e-3,0.1))
#   fig
#   
# }

plotSpectraFishing <- function()
{
  F <- 0.3
  p <- baseparameters()
  W <- c(10, 333, 10000)
  
  defaultplot()
  par(mar=par()$mar + c(0,0,0,10)) # space for legend
  loglogpanel(xlim=c(5e-3, 1),
              ylim=c(5e-3, 0.05),
              xlab="Relative weight $\\textit{w}/\\textit{W}_\\infty$",
              ylab="Biomass spectrum, $\\textit{wN}(\\textit{w})  ")

  for (i in 1:length(W))
  {
    p$F <- F
    p$W <- W[i]
    spec <- spectrum(p)
    lines(spec$w/p$W, spec$w*spec$NprR*p$W^(p$n+p$a-1), lwd=i)
  }
  p$F <- 0
  spec <- spectrum(p)
  lines(spec$w/p$W, spec$w*spec$NprR*p$W^(p$n+p$a-1), lwd=3, col=stdgrey)
  vline(p$etaM)
  legend("right", bty="n", inset=c(-0.5,0), xpd=TRUE,
         legend=c("Unfished",
                  TeX("$\\textit{W}_{\\infty}$ = 10 g"),
                  TeX("$\\textit{W}_{\\infty}$ = 333 g"),
                  TeX("$\\textit{W}_{\\infty}$ = 10 kg")), 
         col=c(stdgrey, black, black, black),
         lwd=c(3,1,3,4),
         title="")
}

plotFishing <- function()
{
  F <- seq(0,1.6,length.out=150)
  W <- c(10, 333, 10000)
  p <- baseparameters()
  
  dat <- data.frame(F=1:(3*length(F)), R=0, SSBprR=0)
  
  ix <- 1
  fig <- ggplot()
  for (j in 1:length(W))
  {
    p$W <- W[j]
    for (i in 1:length(F)) {
      p$F <- F[i]
      spec <- spectrum(p)
      dat$F[ix] <- F[i]
      dat$R[ix] <- spec$R[1]
      dat$SSBprR[ix] <- spec$SSBperR[1]
      dat$W[ix] <- W[j]
      ix <- ix + 1
    }
    #dat$R[dat$R<0] <- 0
    fig <- fig +
      geom_line(data=split(dat, f=dat$W)[[j]], 
                aes(x=F, y=SSBprR*R / max(SSBprR*R)), color=grey(0.5), size=j/2) +
      geom_line(data=split(dat, f=dat$W)[[j]], aes(x=F, y=R/max(R)),color=grey(0.75),size=j/2) 
  }
  
  fig <- fig +
    geom_hline(yintercept = 0.2, linetype="dotted", size=thin) +
    ylim(c(0,1)) +
    xlim(c(0,1.6)) +
    xlab(TeX("Fishing mortality \\textit{F} (yr$^{-1}$)")) +
    ylab("Scaled quantities")
  fig <- mytheme(fig)
  fig
  #ggsave("fishing.pdf", width=singlewidth, height=height)
  #fig
}

plotRecruitmentFishing <- function() 
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
    geom_line(aes(x=Rp, y=Rp/Rp), linetype="dotted", color="black", size=thin)
  
  # Add points for species
  W <- c(10,10000)
  p <- baseparameters()
  for (i in 1:length(W))
  {
    p$W <- W[i]
    p$F <- 0
    N <- spectrum(p)
    fig <- fig +
      geom_point(data=N, aes(x=Rp[1], y=R[1]), size=i*2+1, colour="white") +
      geom_point(data=N, aes(x=Rp[1], y=R[1]), size=i*2)
    
    p$F <- 0.3
    NF <- spectrum(p)
    fig <- fig +
      geom_point(data=NF, aes(x=Rp[1], y=R[1]), size=i*2+1, colour="white") +
      geom_point(data=NF, aes(x=Rp[1], y=R[1]), size=i*2, colour=grey(0.5))
    
    fig <- fig +
      annotate("segment", 
               x=N$Rp[1], xend=NF$Rp[1], 
               y=N$R[1]-0.1, yend=NF$R[1]-0.1,
               arrow=arrow(length=unit(2,"mm")))
  }
  
  fig <- fig +
    scale_y_continuous(limits=c(0,1.2), breaks=seq(0,1,by = 0.2), oob=rescale_none)+
    xlab(TeX("Egg production, $\\textit{R_p}/\\textit{R}_{max}$")) + 
    ylab(TeX("Recruitment, $\\textit{R}/\\textit{R}_{max}$")) 
  
  fig <- semilogx(fig, xlim=c(0.1,1e3))
  
  fig
}



plotFishingComplete <- function()
{
  F <- seq(0, 1.6, length.out=50)
  p <- baseparameters()
  
  panel <- function(W, bAnnotate=FALSE)
  {
    p$W <- W
    dat <- data.frame(F=F, R=F, YprR=F, SSBprR=F)
    for (i in 1:length(F)) {
      p$F <- F[i]
      spec <- spectrum(p)
      dat$R[i] <- spec$R
      dat$YprR[i] <- spec$YprR
      dat$SSBprR[i] <- spec$SSBperR
    }
    
    refs <- calcRefpoints(p)
    Bmsy <- refs$Bmsy / max(dat$SSBprR*dat$R)
    Blim <- refs$Blim / max(dat$SSBprR*dat$R)
    
    ymax <- 1.25
    if (bAnnotate==FALSE)
      ymax <- 1.05
    
    fig <- ggplot(dat) +
      geom_line(aes(x=F, y=SSBprR*R / max(SSBprR*R)), color=grey(0.5), size=thick) +
      geom_line(aes(x=F, y=R/max(R)),color=grey(0.75), size=thick) +
      geom_line(aes(x=F, y=YprR / max(YprR)), linetype="dashed") +
      geom_line(aes(x=F, y=YprR*R / max(YprR*R)),  size=thick) +
      scale_y_continuous(limit=c(0,ymax), oob=rescale_none, breaks=c(0,0.25,0.5,0.75,1)) +
      xlim(c(0,1.6)) +
      xlab(TeX("Fishing mortality, \\textit{F} (yr$^{-1}$)")) +
      ylab("Scaled quantities")
    
    if (bAnnotate==TRUE)
    {
      fig <- fig + geom_vline(xintercept = refs$Fmsy, linetype="dotted", size=thin) +
        #geom_vline(xintercept = refs$Fmax, linetype="dotted", size=thin) +
        geom_vline(xintercept = refs$Flim, linetype="dotted", size=thin) +
        annotate("segment", x=refs$Fmsy, xend=refs$Fmsy, y=1.2, yend=1,arrow=arrow(length=unit(2,"mm"))) +
        annotate("text", x=refs$Fmsy, y=1.2, label=("F[msy]"), parse=TRUE,vjust=0) +
        annotate("segment", x=refs$Fmax, xend=refs$Fmax, y=1.20, yend=1,arrow=arrow(length=unit(2,"mm"))) +
        annotate("text", x=refs$Fmax, y=1.21, label="F[max]", parse=TRUE,vjust=0) +
        annotate("segment", x=refs$Fcrash, xend=refs$Fcrash, y=.5, yend=0,arrow=arrow(length=unit(2,"mm"))) +
        annotate("text", x=refs$Fcrash, y=.51, label="F[crash       ]",parse=TRUE, vjust=0) +
        annotate("segment", x=refs$Flim, xend=refs$Flim, y=0.75, yend=0.5,arrow=arrow(length=unit(2,"mm"))) +
        annotate("text", x=refs$Flim, y=0.76, label="F[lim]", parse=TRUE,vjust=0) +
        annotate("segment", x=refs$Fmsy+0.1, xend=refs$Fmsy, y=Bmsy, yend=Bmsy,arrow=arrow(length=unit(2,"mm"))) +
        annotate("text", x=refs$Fmsy+0.1, y=Bmsy, label=" B[msy]", parse=TRUE,hjust=0) +
        annotate("segment", x=refs$Flim+0.1, xend=refs$Flim, y=Blim+0.05, yend=Blim,arrow=arrow(length=unit(2,"mm"))) +
        annotate("text", x=refs$Flim+0.1, y=Blim+0.05, label=" B[lim]", parse=TRUE,hjust=0) 
    }
    fig <- mytheme(fig)
    fig
  }
  
  fig <- panel(333, TRUE) + xlab("")
  ggsave("TeX/ChapterFishing/fishingmain.pdf",width=0.75*doublewidth, height=0.75*2*height)
  
  
  pdf("TeX/ChapterFishing/fishingsubs.pdf", width=doublewidth, height=height)
  multiplot(panel(10),panel(10000), cols=2)
  dev.off()
}


plotRefPoints <- function() 
{
  W <- 10^seq(0.1,5.5, length.out = 20)
  p <- baseparameters()
  #
  # Calc the ref points for all values of W:
  #
  p$W <- W[1]
  refs <- as.data.frame( calcRefpoints(p) )
  for (i in 2:length(W))
  {
    p$W <- W[i]
    refs <- rbind(refs, calcRefpoints(p))
    #Fmsy[i] <- refs$Fmsy
    #Fcrash[i] <- refs$Fcrash
  }
  refs$W <- W
  #
  # Plot the F-figure
  #
  fig <- ggplot(refs) +
    geom_ribbon(aes(x=W, ymin=Fcrash, ymax=2), fill=grey(0.6)) +
    geom_line(aes(x=W, y=Fmsy), size=thick, color=grey(0.5)) +
    geom_line(aes(x=W, y=Fmax), size=thick, color=grey(0.5), linetype="dashed") +
    geom_line(aes(x=W, y=Fcrash), size=thick) +
    geom_line(aes(x=W, y=Flim), size=thick, linetype="dashed") +
    #scale_y_continuous(limit=c(0,1.5), oob=rescale_none) +
    xlab(TeX("Asymptotic size, $\\textit{W}_{\\infty}$ (g)")) +
    ylab(TeX("Fishing mortality, \\textit{F} (yr$^{-1}$)"))
  fig <- semilogx(fig, ylim=c(0,1.7), xlim=c(4,1e5), label="a")
  #
  # Add the ICES data points:
  #
  dat <- read.csv("Data/ICESrefpoints.csv", header=FALSE)
  names(dat) <- c("W","K","Flim","Fpa","Fmsy")
  fig <- fig +
    geom_point(data=dat, aes(x=W, y=Fmsy), color=grey(0.6)) +
    geom_point(data=dat, aes(x=W, y=Flim))
  #
  # Plot the B figure:
  #
  ix <- refs$Blim==0
  refs$Blim[ix] <- refs$B0[ix]
  figB <- ggplot(refs) +
    geom_line(aes(x=W, y=Bmsy/B0), size=thick, color=grey(0.5)) +
    geom_line(aes(x=W, y=Bmax/B0), size=thick, color=grey(0.5), linetype="dashed") +
    geom_line(aes(x=W, y=Blim/B0), size=thick, linetype="dashed") +
    geom_hline(yintercept = 0.2, linetype="dotted", size=thin) +
    scale_y_continuous(limit=c(0,1.01), oob=rescale_none)+
    xlab(TeX("Asymptotic size, $\\textit{W}_{\\infty}$ (g)")) +
    ylab(TeX("Scaled biomass"))
  
  figB <- semilogx(figB, xlim=c(4,1e5), label=" b") 
  figB
  
  pdf("TeX/ChapterFishing/refpoints.pdf", width=doublewidth, height=height, useDingbats = FALSE)
  multiplot(fig, figB, cols=2)
  dev.off()
}


plotRefvs_a <- function()
{
  p <- baseparameters()
  a <- seq(0,0.7,length.out=20)
  
  panel <- function(W, bYaxis=TRUE)
  {
    p$W <- W
    a0 <- p$a
    dat <- data.frame(a=a,Fmsy=a,Flim=a,Fmax=a,Fcrash=a)
    for (i in 1:length(a))
    {
      p$a <- a[i]
      ref <- calcRefpoints(p)
      dat$Fmsy[i] <- ref$Fmsy
      dat$Flim[i] <- ref$Flim
      dat$Fmax[i] <- ref$Fmax
      dat$Bmsy[i] <- ref$Bmsy
      dat$Blim[i] <- ref$Blim
      dat$Bmax[i] <- ref$Bmax
      dat$B0[i] <- ref$B0
      dat$Fcrash[i] <- ref$Fcrash
    }
    
    fig <- ggplot(dat) +
      geom_ribbon(aes(x=a, ymin=Fcrash, ymax=2), fill=grey(0.6)) +
      geom_line(aes(x=a, y=Fmsy), size=thick, color=grey(0.5)) +
      geom_line(aes(x=a, y=Fmax), size=thick, color=grey(0.5), linetype="dashed") +
      geom_line(aes(x=a, y=Fcrash), size=thick) +
      geom_line(aes(x=a, y=Flim), size=thick, linetype="dashed") +
      geom_vline(xintercept = a0, linetype="dotted", size=thin) +
      scale_y_continuous(limit=c(0,1.5), oob=rescale_none) +
      ggtitle(TeX(paste("$\\textit{W}_{\\infty}$ = ",W," g")))
    fig <- mytheme(fig, bYaxis = bYaxis)
    
    # ix <- dat$Blim==0
    # dat$Blim[ix] <- dat$B0[ix]
    # figB <- ggplot(dat, aes(x=a)) +
    #   geom_line(aes(y=Bmsy/B0), size=thick, color=grey(0.5)) +
    #   geom_line(aes(y=Blim/B0), size=thick, linetype="dashed") +
    #   geom_line(aes(y=Bmax/B0), size=thick, color=grey(0.5), linetype="dashed") +
    #   scale_y_continuous(limit=c(0,1.01), oob=rescale_none) +
    #   geom_vline(xintercept = a0, linetype="dashed", size=thin)
    # figB <- mytheme(figB, bYaxis=bYaxis)
    
    fig
  }
  
  fig1 <- panel(10) + ylab(TeX("$\\textit{F}_{msy}$ (yr$^{-1}$)")) + xlab("")
  fig2 <- panel(333, bYaxis=FALSE) + xlab("Physiological mortality") + ylab("")
  fig3 <- panel(10000, bYaxis=FALSE) + xlab("") + ylab("")
  
  grid = plot_grid(fig1,fig2,fig3, ncol=3, align="h",
                   rel_widths = c(1.26,1,1))
  ggsave("TeX/ChapterFishing/refpointsvsa.pdf", grid, width=doublewidth, height=height)
}


plotRefvs_A <- function()
{
  p <- baseparameters()
  A <- p$A*10.^seq(log10(0.5), log10(2),length.out=20)
  
  panel <- function(W)
  {
    p$W <- W
    dat <- data.frame(A=A,Fmsy=A,Flim=A,Fmax=A,Fcrash=A)
    for (i in 1:length(A))
    {
      p$A <- A[i]
      ref <- calcRefpoints(p)
      dat$Fmsy[i] <- ref$Fmsy
      dat$Flim[i] <- ref$Flim
      dat$Fmax[i] <- ref$Fmax
      dat$Fcrash[i] <- ref$Fcrash
    }
    
    fig <- ggplot(dat) +
      geom_ribbon(aes(x=A, ymin=Fcrash, ymax=2), fill=grey(0.6)) +
      geom_line(aes(x=A, y=Fmsy), size=thick, color=grey(0.5)) +
      geom_line(aes(x=A, y=Fmax), size=thick, color=grey(0.5), linetype="dashed") +
      geom_line(aes(x=A, y=Fcrash), size=thick) +
      geom_line(aes(x=A, y=Flim), size=thick, linetype="dashed") +
      scale_y_continuous(limit=c(0,1.5), oob=rescale_none) +
      ggtitle(TeX(paste("$W_{\\infty} = ",W," g")))
    fig <- mytheme(fig)  
    fig
  }
  
  fig1 <- panel(10) + ylab(TeX("F (yr$^{-1}$)")) + xlab("")
  fig2 <- panel(333) + xlab("Growth coef, A") + ylab("")
  fig3 <- panel(10000) + xlab("") + ylab("")
  
  pdf("refpointsvsA.pdf", width=doublewidth, height=height)
  multiplot(fig1,fig2,fig3, cols=3)
  dev.off()
}


# plotRefvs_eps <- function()
# {
#   p <- baseparameters()
#   epsR0 <- p$epsR
#   epsR <- p$epsR *10.^seq(-1,1,length.out=20)
#   
#   panel <- function(W)
#   {
#     p$W <- W
#     dat <- data.frame(epsR=epsR,Fmsy=epsR,Flim=epsR,Fmax=epsR,Fcrash=epsR)
#     for (i in 1:length(epsR))
#     {
#       p$epsR <- epsR[i]
#       ref <- calcRefpoints(p)
#       dat$Fmsy[i] <- ref$Fmsy
#       dat$Flim[i] <- ref$Flim
#       dat$Fmax[i] <- ref$Fmax
#       dat$Fcrash[i] <- ref$Fcrash
#     }
#     
#     fig <- ggplot(dat) +
#       geom_ribbon(aes(x=epsR, ymin=Fcrash, ymax=2), fill=grey(0.6)) +
#       geom_line(aes(x=epsR, y=Fmsy), size=thick, color=grey(0.5)) +
#       geom_line(aes(x=epsR, y=Fmax), size=thick, color=grey(0.5), linetype="dashed") +
#       geom_line(aes(x=epsR, y=Fcrash), size=thick) +
#       geom_line(aes(x=epsR, y=Flim), size=thick, linetype="dashed") +
#       geom_vline(xintercept = epsR0, linetype="dotted", size=thin) +
#       scale_y_continuous(limit=c(0,1.5), oob=rescale_none) +
#       ggtitle(TeX(paste("$W_{\\infty} = ",W,"\ g")))
#     fig <- semilogx(fig)  
#     fig
#   }
#   
#   fig1 <- panel(10) + ylab(TeX("F (yr$^{-1}$)")) + xlab("")
#   fig2 <- panel(333) + xlab(TeX("Repro. eff. $\\epsilon_r$")) + ylab("")
#   fig3 <- panel(10000) + xlab("") + ylab("")
#   
#   pdf("refpointsvsepsr.pdf", width=doublewidth, height=height)
#   multiplot(fig1,fig2,fig3, cols=3)
#   dev.off()
# }


plotRefvs_eps2 <- function()
{
  p <- baseparameters()
  epsR0 <- p$epsR
  epsR <- 10.^seq(-1,1,length.out=20)
  W <- c(10, 333, 10000)
  
  dat <- data.frame(epsR=1:(length(epsR)*length(W)), W=0, Fmsy=0)
  ix <- 1
  fig <- ggplot()
  for (j in 1:length(W))
  {
    p$W <- W[j]
    for (i in 1:length(epsR))
    {
      p$epsR <- epsR[i] * epsR0
      ref <- calcRefpoints(p)
      dat$Fmsy[ix] <- ref$Fmsy
      dat$W[ix] <- p$W
      dat$epsR[ix] <- p$epsR
      ix <- ix + 1
    }
    idx <- (ix-length(epsR)):(ix-1)
    fig <- fig +
      geom_line(data=dat[idx,], aes(x=epsR, y=Fmsy), size=j/2)
  }
  fig <- semilogx(fig, xlim=c(0.003,0.15))  +
    scale_y_continuous(limit=c(0,1), oob=rescale_none) +
    geom_vline(xintercept = epsR0, linetype="dotted", size=thin) +
    xlab(TeX("Recruitment efficiency $\\epsilon_\\textit{R}$")) +
    ylab(TeX("$\\textit{F}_{msy}$ (yr$^{-1}$)"))
    
  ggsave("TeX/ChapterFishing/refsvs_eps.pdf", width=singlewidth, height=height)
  fig

}



plotRefUncertainty <- function()
{
  nPoints <- 100  # no. of calculations to make
  W <- 10^seq(0.5,5, length.out = 10)
  p <- baseparameters()
  
  p$W <- 100
  refs <- as.data.frame( calcRefpoints(p) )
  for (i in 1:nPoints) 
  {
    # 
    # Calc random ref points:
    #
    p$a <- exp(rnorm(n=1,mean=-1.36,sd=0.76))
    p$epsR <- exp(rnorm(n=1, mean=-3.91, sd=2.3))
    #
    # Calc the ref points for all values of W:
    #
    p$W <- W[1]
    for (i in 2:length(W))
    {
      p$W <- W[i]
      refs <- rbind(refs, calcRefpoints(p))
    }
  }
  
  fig <- ggplot(refs,aes(x=W, y=Fmsy)) +
    geom_point() + geom_density2d() +
    geom_line(aes(x=W, y=W^-0.25))
  #geom_point(aes(x=W, y=Flim))
  fig <- semilogx(fig)
  fig
}

panelYieldvsEtaF <- function(p, n=30, bYaxis=TRUE)
{
  etaF <- 10^seq(-2,0,length.out=n)
  F <- seq(0,2,length.out=n)
  
  dat = data.frame(F=1:(n*n), Y=0, etaF=0, R=0)
  
  ix <- 1
  for (j in 1:length(F))
  {
    for (i in 1:length(etaF))
    {
      p$etaF <- etaF[i]
      p$F <- F[j]
      spec <- spectrum(p)
      if (p$etaF == 1)
        dat$Y[ix] <- 0
      else
        dat$Y[ix] <- max(0, calcYield(p,spec))
      dat$R[ix] <- max(0, spec$R[1])
      dat$etaF[ix] <- etaF[i]
      dat$F[ix] <- F[j]
      ix <- ix + 1
    }
  }
  fig <- ggplot(data=dat, aes(x=etaF, y=F, z=Y, fill=R)) +
    #geom_raster() + #guides(fill=FALSE)+
    geom_contour(colour="black") +
    geom_vline(xintercept = p$etaM, linetype="dotted", size=thin) +
    xlab("") + ylab("")
  fig <- semilogx(fig, xlim=c(min(etaF),1), bYaxis=bYaxis, bXaxis = FALSE)
  fig
}

panelMSYvsetaF <- function(p, n=30, bYaxis=TRUE)
{
  W <- c(10,333,10000)
  etaF <- 10^seq(-2,0,length.out=n+1)
  etaF <- etaF[1:(length(etaF)-1)]
  
  ix <- 1
  dat <- data.frame(etaF=1:(length(W)*n), Ymsy=0, W=0,Fmsy=0)
  fig <- ggplot()
  for (j in 1:length(W))
  {
    p$W <- W[j]
    # Find Fmsy for each etaF
    for (i in 1:length(etaF))
    {
      p$etaF <- etaF[i]
      #cat(etaF[i],", ")
      refs <- calcRefpoints(p)
      dat$etaF[ix] <- etaF[i]
      dat$Ymsy[ix] <- refs$Ymsy
      dat$Fmsy[ix] <- refs$Fmsy
      dat$W[ix] <- p$W
      ix <- ix + 1
    }
    idx <- (n*(j-1)+1):(n*j)
    dat$Ymsy[idx] <- dat$Ymsy[idx]/max(dat$Ymsy[idx])
    fig <- fig +
      geom_line(data=dat[idx,], aes(x=etaF, y=Ymsy), size=j/2)
  }  
  fig <- semilogx(fig, xlim=c(min(etaF),1), bYaxis=bYaxis) +
    geom_vline(xintercept=p$etaM, linetype="dotted", size=thin) +
    xlab("") + ylab("") +
    ylim(c(0,1))
  fig
}

panelSelectivity <- function(p)
{
  #
  # Plot selectivity function
  #
  w <- logseq(0.0005,1)
  p$W <- 1
  p$F <- 1
  panel <- ggplot() +
    geom_line(aes(x=w, y=as.vector(p$funcFishing(w,p))), size=thick) +
    geom_vline(xintercept = p$etaF,linetype="dotted", size=thin) +
    scale_y_continuous(limit=c(0,1.4), breaks=c(0, 0.5, 1)) +
    ylab("") + xlab("") +
    annotate("text", x=p$etaF, y=1, 
             parse=TRUE, label="eta[F]", hjust=1.2, vjust=-0.4 ) +
    coord_fixed(ratio=1)
  panel <- semilogx(panel, ylim=c(0,1.2), bXaxis = FALSE, bYaxis = FALSE)
  panel
}

plotYield <- function(n=50)
{
  p<-baseparameters()
  p$W <- 10000
  p$F <- 1
  
  p$funcFishing <- fishingKnifeedge
  fig1 <- list(panelSelectivity(p), panelYieldvsEtaF(p,n), panelMSYvsetaF(p,n))
  fig1[[1]] <- fig1[[1]] + 
    ggtitle("Knife edge") #+
    #ylab("Selectivity")
  fig1[[2]] <- fig1[[2]] + ggtitle("Knife edge") +
    ylab(TeX("Fishing, F (yr$^{-1}$)"))
  fig1[[3]] <- fig1[[3]] +
    ylab(TeX("$\\textit{Y}/\\textit{Y}_{msy}$"))
  
  
  p$funcFishing <- fishingTrawl
  fig2 <- list(panelSelectivity(p), 
               panelYieldvsEtaF(p,n, bYaxis=FALSE), panelMSYvsetaF(p,n, bYaxis=FALSE))
  fig2[[1]] <- fig2[[1]] + 
    ggtitle("Trawl") +
    xlab(TeX("Rel. weight $\\textit{w}/\\textit{W}_{\\infty}$ (g)"))
  fig2[[2]] <- fig2[[2]] + ggtitle("Trawl")
  #  xlab(TeX("$\\eta_F$"))
  fig2[[3]] <- fig2[[3]] +
    xlab(TeX("Selectivity, $\\eta_\\textit{F}$"))
  
  p$funcFishing <- fishingGillnet
  fig3 <- list(panelSelectivity(p), 
               panelYieldvsEtaF(p,n,  bYaxis=FALSE), panelMSYvsetaF(p,n,bYaxis=FALSE))
  fig3[[2]] <- fig3[[2]] + 
    ggtitle("Gill net")
    
    
    pdf("yield.pdf", width=doublewidth, height=3*height)
  multiplot(fig1[[1]],fig1[[2]],fig1[[3]],
            fig2[[1]],fig2[[2]],fig2[[3]],
            fig3[[1]],fig3[[2]],fig3[[3]],
            cols=3)
  dev.off()
  
  grid <- plot_grid(
                    fig1[[2]],fig2[[2]],fig3[[2]],
                    fig1[[3]],fig2[[3]],fig3[[3]],
                    nrow=2, align="hv")
  ggsave("TeX/ChapterFishing/yield.pdf", grid, width=doublewidth, height=2*height)
}


# x <- cbind(etaF = 10^seq(-4,0,length.out=10),  F = seq(0,3,length.out=10))
# YY <- function(x) {
#   cat(x)
#   p$etaF <- x[1]
#   p$F <- x[2]
#   spec <- spectrum(p)
#   max(0, calcYield(p, spec))
# }
# apply(x, 1, YY)

plotAllChapterFishing <- function()
{
  plotSelectivity()
  pdfplot(FUN=plotSpectraFishing, "TeX/ChapterFishing/spectrumfishing.pdf", width=doublewidth, height=height)
  plotTrawl()
  plotSpectraFishing()
  
  grid <- plot_grid(plotFishing(), plotRecruitmentFishing(), ncol=2, align="h")
  ggsave("TeX/ChapterFishing/fishing.pdf", grid, width=doublewidth, height=height)
  
  plotFishingComplete()
  plotRefPoints()
  plotRefvs_a()
  #plotRefvs_A()
  plotRefvs_eps2()
  plotYield(n=50)
  #plotRefUncertainty()
}

