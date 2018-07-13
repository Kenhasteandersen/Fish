source("R/basetools.R")
source("R/basefunctions.R")
source("R/baseparameters.R")
source("R/community.R")

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



plotSpectra_ConsumerResourceModel <- function(W=1000) {
  p = paramConsumerResourcemodel(W, 1e-8) # std. model driven by SR relation
  r_std = runspectrum(p)
  
  p = paramConsumerResourcemodel(W, 1000) # Fully driven by physiological recruitment
  r = runspectrum(p)
  
  defaultplot(mfcol=c(2,2))
  #
  # Spectra
  #
  wrange = c(0.1,W)
  loglogpanel(xlim = wrange, 
              ylim=c(5, 1e6),
              ylab = "Spectrum, \\textit{Nw^2/(N(w_0)w_0^2)} (g)",
              xlab = "Weight (g)", label=TRUE)
  # Resource:
  # wR = r$resource$wR
  # NR = r$resource$NR/r$N[r$nSave,1,1]*(wR/p$w0)^2
  # KR = p$KR * wR^(p$lambdaR)/r$N[r$nSave,1,1]*(wR/p$w0)^2
  # ribbon(wR, NR, KR)
  # lines(wR, NR, lwd=2, col=stdgrey)
  # lines(wR, KR, lwd=2, col=stdgrey, lty=dashed)

  w = r$fish$w
  lines(w, r_std$N[r$nSave,1,]/r_std$N[r$nSave,1,1]*(w/p$w0)^2, lty=dashed, lwd=2)
  lines(w, r$N[r$nSave,1,]/r$N[r$nSave,1,1]*(w/p$w0)^2, lwd=2)
  
  legend("bottomright", 
         legend=c("Early-life density dep.", "Emergent density dep."), 
         lwd=2, lty=c(dashed,1),
         bty="n")
  #
  # Feeding and losses:
  #
  semilogxpanel(xlim = wrange, 
                ylim=c(0,1),
                ylab = "Feeding and loss",
                xlab = "Weight (g)",
                label = TRUE)
  # feeding
  f0 = r_std$fish$f
  f = r$fish$f
  ribbon(w, ymax=f0, ymin=f)
  lines(w, f0, lty=dashed, lwd=2, col=darkgrey)
  lines(w, f, lwd=2, col=darkgrey)
  
  text(x=0.1, y=0.64, "Feeding", pos=4)
  fmin = min(f)
  arrows(x1=w[which(f==fmin)], y1=fmin,
         x=w[which(f==fmin)], y=fmin-0.1,
         length=0.1)
  text(x=60, y=fmin-0.1,
       labels="Resource-driven\nbottleneck",
       pos=1)
  
  # loss
  ref = p$h*w^(p$n-1)
  mort = p$SizeBasedMortCoef*w^(p$n-1)
  loss_std = mort/ref
  loss = (mort+r$fish$muP)/ref
  ribbon(w, ymax=loss, ymin=loss_std)
  lines(w, loss_std, lty=dashed, lwd=2)
  lines(w, loss, lwd=2)
  
  text(x=0.1, y=0.05, "Loss", pos = 4)
  lossmax = max(loss)
  w1 = w[which(loss==lossmax)]
  arrows(x1=w1, y1=lossmax,
         x=w1, y=lossmax+0.1,
         length=0.1)
  text(x=w1, y=lossmax+0.075,
       labels="Cannibalistic\nbottleneck",
       pos=3, adj=0.5)
  #
  # Mortality
  #
  semilogxpanel(xlim = wrange,
                ylim = c(0,7),
                ylab = "Mortality (yr$^{-1}$)", 
                xlab = "Weight (g)",
                label = TRUE)
  
  ribbon(w, ymin=mort, ymax=mort+r$fish$muP)
  lines(w, mort, lty=dashed, lwd=2)
  lines(w, mort+r$fish$muP, lwd=2)
  #
  # Weight at age
  #
  defaultpanel(xlim = c(0,30),
               ylim = c(0,W),
               xlab = "Age",
               ylab = "Weight (g)",
               label = TRUE)
  
  wa_std = calcWeightAtAge(p, r_std) # calc without DD
  wa = calcWeightAtAge(p, r)
  
  ribbon(wa$ages, wa_std$w, wa$w)
  lines(wa_std$ages, wa_std$w, lwd=2, lty=dashed)
  lines(wa$ages, wa$w, lwd=2)
}

plotFishingvsWf <- function(W=1000, n=20) {
  
  panelFishingvsWf <- function(p, F=seq(0,1, length.out = 9), etaF=10^seq(-3,0,length.out = 10)) {
    Y = matrix(nrow=length(etaF), ncol=length(F))
    for (i in 1:length(F)) {
      
      for (j in 1:length(etaF)) {
        p$F = F[i]
        p$etaF = etaF[j]
        res = runspectrum(p)
        Y[j,i] <- mean(res$R$Y[floor(res$nSave/4):res$nSave])
      }
    }
    return(Y/max(Y))
    
    #semilogxpanel(xlim=range(etaF), ylim=range(F) )
    #contour(x=etaF, y=F, z=t(Y),
    #        drawlabels=FALSE, frame.plot=FALSE, add=TRUE,
    #        levels=seq(0,1,length.out = 10), axes=FALSE)
  }
  
  #
  # Do the calculations
  #
  if (bRecalcExpensiveFunctions==TRUE) {
    F = seq(0,1, length.out = n)
    etaF=10^seq(-2.9,0,length.out = n)
    Y = array(dim=c(3, length(etaF), length(F)))
    
    # Standard model driven purely by SR relation:
    p = paramConsumerResourcemodel(W=W, facRmax = 1e-8)
    Y[1,,] = panelFishingvsWf(p, F=F, etaF=etaF)
    
    # Emergent functional response without cannibalism:
    p = paramCommunitymodel(W=W, facRmax=1e8)
    p$thetaFish = 0
    Y[2,,] = panelFishingvsWf(p, F=F, etaF=etaF)
    
    # Emergent functional response with cannibalism:
    p = paramCommunitymodel(W=W, facRmax=1e8)
    Y[3,,] = panelFishingvsWf(p, F=F, etaF=etaF)
    
    save(file="Data/ConsumerResourceYield.Rsave", list=c("Y","F","etaF"))
  }
  else
    load('Data/ConsumerResourceYield.Rsave')
  #
  # Plots
  #
  defaultplothorizontal(nPanel=3)
  # make room for titles:
  oma <- par()$oma
  oma[top]<-1.5
  par(oma=oma) 
  
  tightaxes()
  semilogxpanel(xlim=etaF, ylim=F, ylab="Fishing mortality (yr$^{-1})")
  polygon(c(p$etaM,1,1,p$etaM), c(0,0,1,1), col=lightgrey, border=NA)
  contour(x=etaF, y=F, z=Y[1,,],
          drawlabels=FALSE, frame.plot=FALSE, add=TRUE,
          levels=seq(0,1,length.out = 10), axes=FALSE)
  semilogxpanel(xlim=etaF, ylim=F, ylab="Fishing mortality (yr$^{-1})", new=TRUE)
  mtext(side=top, line=0.5, "Early density dep.")
  
  semilogxpanel(xlim=etaF, ylim=F, yaxis=FALSE,
                xlab="Start of fishing relative to asymptotic size")
  polygon(c(p$etaM,1,1,p$etaM), c(0,0,1,1), col=lightgrey, border=NA)
  tightaxes()
  contour(x=etaF, y=F, z=Y[2,,],
          drawlabels=FALSE, frame.plot=FALSE, add=TRUE,
          levels=seq(0,1,length.out = 10), axes=FALSE)
  semilogxpanel(xlim=etaF, ylim=F, yaxis=FALSE,
                xlab="Start of fishing relative to asymptotic size", new=TRUE)
  mtext(side=top, line=0.5, "Competition")
  
  semilogxpanel(xlim=etaF, ylim=F, yaxis=FALSE)
  polygon(c(p$etaM,1,1,p$etaM), c(0,0,1,1), col=lightgrey, border=NA)
  tightaxes()
  contour(x=etaF, y=F, z=Y[3,,],
          drawlabels=FALSE, frame.plot=FALSE, add=TRUE,
          levels=seq(0,1,length.out = 10), axes=FALSE)
  semilogxpanel(xlim=etaF, ylim=F, yaxis=FALSE, new=TRUE)
  mtext(side=top, line=0.5, "Competition + cannibalism")
}

plotOptimalEtaF <- function(W=1000, nPoints = 10) {
  
  opt = function(theta) {
    p$F=theta[1]
    p$etaF=10^theta[2]
    res = runspectrum(p)
    #cat(theta,":",-mean(res$R$Y[floor(res$nSave/4):res$nSave]),"\n")
    return(-mean(res$R$Y[floor(res$nSave/4):res$nSave])/1e-12)
  }
  #
  # Calculations
  #
  if (bRecalcExpensiveFunctions==TRUE) {
    Rmaxscaled = 10^seq(-1, 3, length.out = nPoints)
    etaFoptRmax = matrix(nrow=2, ncol=length(Rmaxscaled))
    for (thetaFish in c(0,1)) {
      theta = c(0.5,-1)
      for (i in 1:length(Rmaxscaled)) {
        p = paramConsumerResourcemodel(W=W, facRmax = Rmaxscaled[i])
        p$thetaFish = thetaFish
        mm = optim(fn=opt, par=theta, method="L-BFGS-B", 
                   lower=c(0,-3), upper=c(3,0))
        theta = mm$par
        cat(theta,"\n")
        etaFoptRmax[thetaFish+1, i] = 10^mm$par[2]
      }  
    }
    
    W = 10^seq(1,5, length.out = nPoints)
    etaFoptW = matrix(nrow=3, ncol=length(W))
    theta = list(c(3,-0.3), c(3,-0.3), c(3,-0.3))
    for (i in 1:length(W)) {
      # SR relation
      p = paramConsumerResourcemodel(W=W[i], facRmax = 1e-8)
      mm = optim(fn=opt, par=theta[[1]], method="L-BFGS-B", 
                 lower=c(0,-3), upper=c(3,0))
      theta[[1]] = mm$par
      cat(theta[[1]],"\n")
      etaFoptW[1, i] = 10^mm$par[2]
      
      # No cannibalism
      p = paramConsumerResourcemodel(W=W[i], facRmax = 1e8)
      p$thetaFish = 0
      mm = optim(fn=opt, par=theta[[2]], method="L-BFGS-B", 
                 lower=c(0,-3), upper=c(3,0))
      theta[[2]] = mm$par
      cat(theta[[2]],"\n")
      etaFoptW[2, i] = 10^mm$par[2]
      
      # With cannibalism
      p = paramConsumerResourcemodel(W=W[i], facRmax = 1e8)
      mm = optim(fn=opt, par=theta[[3]], method="L-BFGS-B", 
                 lower=c(0,-3), upper=c(3,0))
      theta[[3]] = mm$par
      cat(theta[[3]],"\n")
      etaFoptW[3, i] = 10^mm$par[2]
    }
    save(file="Data/OptimalEtaF.Rdata", list=c("Rmaxscaled","etaFoptRmax",
                                        "W", "etaFoptW"))
  }
  else
    load("Data/OptimalEtaF.Rdata")
  #
  # Plots
  #
  defaultplothorizontal(nPanel=2)
  
  loglogpanel(xlim=Rmaxscaled, ylim=c(0.001,1), 
              xlab="\\textit{R}_{max}/\\tilde{\\textit{R}}_{max}", 
              ylab="Optimal selection, $\\eta^*_\\textit{F}",
              label=TRUE)
  
  lines(Rmaxscaled, etaFoptRmax[1,], lwd=2)
  lines(Rmaxscaled, etaFoptRmax[2,], lwd=2, lty=dashed)
  lines(Rmaxscaled, etaFoptRmax[1,1]*Rmaxscaled/Rmaxscaled, lty=dotted)
  
  legend("bottomleft", 
         legend=c("Stock-recruitment", "Emergent", "Emergent with cannibalism"), 
         lwd=c(1,2,2), lty=c(dotted,1,dashed),
         bty="n")
  
  loglogpanel(xlim=W, ylim=c(0.001,1),
              xlab="Asymptotic size (g)", yaxis = FALSE)
  makepanellabel(line = -2)
  
  lines(W, etaFoptW[1,], lty=dotted)
  lines(W, etaFoptW[2,], lwd=2)
  lines(W, etaFoptW[3,], lwd=2, lty=dashed)
}

plotSpatialPlaice = function() {
  library('ggplot2')
  library(rworldmap)
  library(rworldxtra)
  library(grid)
  library(gridExtra)
  library(rgdal)
  library(mapplots)
  library(scales)
  library(mapproj)
  
  
  wmap <- readOGR(dsn=paste('Data/shapefiles', sep =''), layer="ne_110m_land")
  wmap_df <- fortify(wmap)
  
  wmap_df <- fortify(wmap)
  
  IcesMap <- readOGR(dsn=paste('Data/shapefiles', sep =''), 
                     layer = "ices_rectangles")
  
  Icesmap_df <- fortify(IcesMap)
  
  newmap <- getMap(resolution = 'low')
  newmap_df <- fortify(newmap)
  
  xL <- c(-6,12)
  yL <- c(50,60)
  
  # Get data 
  # Plot age 0 
  theme_opts <- list(theme(panel.grid.minor = element_blank(),
                           panel.grid.major = element_blank(),
                           panel.background = element_rect(fill = 'white'),
                           plot.background = element_rect(fill="white"),
                           panel.border = element_blank(),
                           plot.title = element_text(size=12)))
  
  df <- read.csv('Data/MapData_NS_Demersals.csv')
  pop_n <- 0.9
  # Subset plaice
  plaice <- df[which(df$Species == "Pleuronectes platessa"),]
  # Lat/long for plaice
  coords <- ices.rect(plaice$Subarea)
  plaice <- data.frame(plaice, coords)
  
  # Read the mean weights 
  weights <- read.csv('Data/NS_IBTS_MeanWeights.csv', sep = ',')
  head(weights)
  
  plaice_weight <- weights[weights$Species == "Pleuronectes platessa",]
  plaice$weight <-  0.0093*(plaice$LngtClas/10)^3.03
  
  plaice1 <- plaice[which(plaice$weight < 30),] # Only less than 30g
  
  # Normalize Cpue for plotting
  plaice1$lCPUE = (plaice1$CPUE_number_per_hour)/(max(plaice1$CPUE_number_per_hour))
  
  plaice1 <- plaice1[with(plaice1, order(lCPUE)),]
  Idx <- which.min((pop_n*sum(plaice1$CPUE_number_per_hour)-cumsum(sort(plaice1$CPUE_number_per_hour, decreasing = F)))^2)
  
  plaice1_plot <- plaice1[1:Idx,]
  
  p1plaice <- ggplot(data =plaice1_plot,aes(x =lon,y=lat))+
    geom_tile(aes(fill=lCPUE))+
    geom_polygon(data=newmap_df, aes(x =long,y=lat,group=group), fill = 'gray', color = 'black') + 
    scale_fill_continuous(low = "gray90",high = "gray0", guide = FALSE,limits = c(0,1))+
    theme_opts+
    #theme(legend.position="topright")+
    xlab('Longitude')+theme(axis.text.y = element_blank())+ylab('Latitude')+
    ggtitle('Plaice < 30 g')+
    coord_map(xlim = xL, ylim = yL)
  p1plaice
  
  # For larger plaice
  plaice2 <- plaice[which(plaice$weight < 100),] 
  plaice2$lCPUE = (plaice2$CPUE_number_per_hour)/(max(plaice2$CPUE_number_per_hour))
  plaice2 <- plaice2[with(plaice2, order(lCPUE)),]
  Idx <- which.min((pop_n*sum(plaice2$CPUE_number_per_hour)-cumsum(sort(plaice2$CPUE_number_per_hour, decreasing = F)))^2)
  
  plaice2_plot <- plaice2[1:Idx,]
  
  p2plaice <- ggplot(data =plaice2_plot,aes(x =lon,y=lat))+
    geom_tile(aes(fill=lCPUE))+
    geom_polygon(data=newmap_df, aes(x =long,y=lat,group=group), fill = 'gray', color = 'black') + 
    scale_fill_continuous(low = "gray90",high = "gray0", guide = FALSE,limits = c(0,1))+
    theme_opts+
    #theme(legend.position="none")+
    ggtitle('Plaice < 100 g')+
    xlab('Longitude')+
    ylab('')+theme(axis.text.y = element_blank())+
    coord_map(xlim = xL, ylim = yL)
  
  p2plaice
  
  # Only the largest plaice
  plaice3 <- plaice[which(plaice$weight > 200),]
  plaice3$lCPUE = (plaice3$CPUE_number_per_hour)/(max(plaice3$CPUE_number_per_hour))
  plaice3 <- plaice3[with(plaice3, order(lCPUE)),]
  Idx <- which.min((pop_n*sum(plaice3$CPUE_number_per_hour)-cumsum(sort(plaice3$CPUE_number_per_hour, decreasing = F)))^2)
  
  plaice3_plot <- plaice3[1:Idx,]
  
  
  p3plaice <- ggplot(data =plaice3_plot,aes(x =lon,y=lat))+
    # stat_contour(geom="polygon", aes(fill=..level..))+
    geom_tile(aes(fill=lCPUE))+
    geom_polygon(data=newmap_df, aes(x =long,y=lat, group=group), fill = 'gray', color = 'black') + 
    scale_fill_continuous(low = "gray90",high = "gray0", guide = FALSE,limits = c(0,1))+
    theme_opts+
    #theme(legend.position="none")+
    xlab('Longitude')+ylab('')+theme(axis.text.y = element_blank())+
    ggtitle('Plaice > 200 g')+
    coord_map(xlim = xL, ylim = yL)
  p3plaice
  
  ggsave("Chapter7b/SpatialPlaice.pdf",
         grid.arrange(p1plaice, p2plaice,p3plaice ,ncol=3),
         width=doublewidth, height=height)
}

plotConsumerResourceVariation <- function(W=1000) {
  
  p <- paramConsumerResourcemodel(W,1)
  
  r0scaling = p$rR / (p$f0*p$h)
  r0scaled = c(0.01, 1, 100)
  
  Rmaxscaling = p$A*p$KR*p$w0^(-p$a)*W^(-2-p$q+2*p$n+p$a)
  Rmaxscaled = c(10, 1, 0.1)
  
  p$Rmax = 1e-8
  r_std = runspectrum(p)
  
  defaultplot(mfcol=c(3,3))
  
  for (i in 1:length(r0scaled)) {
    
    for (j in 1:length(Rmaxscaled)) {
      
      p$Rmax = Rmaxscaled[j] * Rmaxscaling
      p$rR = r0scaled[i] * r0scaling
      r = runspectrum(p)
      
      semilogxpanel(xlim=c(0.1,W),
                    ylim=c(0,1))
      w = r$fish$w
      f0 = r_std$fish$f
      f = r$fish$f
      ribbon(w, ymax=f0, ymin=f)
      lines(w, f0, lty=dashed, lwd=2)
      lines(w, f, lwd=2)
      
      # loss
      ref = p$h*w^(p$n-1)
      mort = p$SizeBasedMortCoef*w^(p$n-1)
      loss_std = mort/ref
      loss = (mort+r$fish$muP)/ref
      ribbon(w, ymax=loss, ymin=loss_std)
      lines(w, loss_std, lty=dashed, lwd=2)
      lines(w, loss, lwd=2)
      
      # wrange = c(0.1,W)
      # loglogpanel(xlim = wrange, 
      #             ylim=c(10, 1e6),
      #             ylab = "Spectrum, \\textit{Nw^2/(N(w_0)w_0^2)} (g)",
      #             xlab = "Weight (g)")
      # 
      # w = r$fish$w
      # lines(w, r_std$N[r$nSave,1,]/r_std$N[r$nSave,1,1]*(w/p$w0)^2, lty=dashed, lwd=2)
      # lines(w, r$N[r$nSave,1,]/r$N[r$nSave,1,1]*(w/p$w0)^2, lwd=2)
      
    }
  }
}

plotConsumerResourceVariation2 <- function(W=1000) {
  
  facRmax = 10^seq(-1, 3, length.out = 6)
  
  f = matrix(nrow=length(facRmax), ncol = 154)
  loss = f
  B = f
  BR= matrix(nrow=length(facRmax), ncol=273)
  
  for (i in 1:length(facRmax)) {
    # with cannibalism
    p = paramConsumerResourcemodel(W, facRmax = facRmax[i])
    r = runspectrum(p)
    
    ix = floor(r$nSave/4):r$nSave
    w = r$fish$w
    f[i,] =  r$fish$f 
    loss[i,] = (r$fish$muP + p$SizeBasedMortCoef*w^(p$n-1)) /
                       (p$h*w^(p$n-1) )
    B[i,] =  colMeans( r$N[ix,1,], dims=1) *w^2
    BR[i,] = r$resource$NR*r$resource$wR^2
  }
  
  defaultplot(mfcol = c(1,2))
  #
  # Spectra
  #
  loglogpanel(xlim=w, ylim=c(1e-9,1e-1), 
              xlab="Weight (g)", ylab="Biomass spectrum",
              label=TRUE)
  lines(r$w, p$KR * w^(p$lambdaR+2), lwd=3, col=stdgrey)
  for (i in 1:length(facRmax)) {
    lines(w, B[i,])
    lines(r$resource$wR, BR[i,], col=stdgrey)
  }
  arrows(x0=c(0.01,1), y0=c(5e-9, 2e-2),
         x1=c(0.01,1), y1=c(5e-5, 0.6e-3),
         length=0.1)
  text(x=1e-2, y=5e-9, pos=4,
       labels = TeX("Increasing $\\textit{R}_\\mathrm{max}"))
  #text(x=7e-3, y=2e-2, "Resource", col=stdgrey, pos=4)
  #
  # Feeding and loss:
  #
  semilogxpanel(xlim=w, ylim = c(0,1),
                xlab="Weight (g)", ylab="feeding and loss",
                label=TRUE)
  for (i in 1:length(facRmax))  {
    lines(w, f[i,], col=stdgrey)
    lines(w, loss[i,])
  }
  arrows(x0=c(1.4, 3e2), y0=c(0.05, 0.63), 
         y1=c(0.4, 0.4), 
         length=0.1) 
  text(1e-3, 0.65, "Feeding level", col=stdgrey, pos=4)
  text(1e-3, 0.15, "Loss", pos=4)
}

plotPlaiceGrowth = function() {
  dat = read.csv("Data/Plaice growth.csv", header=TRUE)
  load("Data/plaice_timeseries.Rda") # in  "plaice"
  
  defaultplot(mfcol=c(1,2))
  #
  # Cohorts
  #
  defaultpanel(xlim=c(1980,2016), ylim=c(5,45),
               xlab="Year", "Length (cm)")
  ribbon(x=plaice$Year, ymax=10*plaice$StockSize/max(plaice$StockSize)+5, ymin=plaice$StockSize*0)
  defaultpanel(xlim=c(1980,2016), ylim=c(5,45),
               xlab="Year", "Length (cm)", new=TRUE)
  
  for (age in 0:6) {
    ix = (dat$AGE==age)
    lines(dat$YEAR[ix], dat$mean[ix], lwd=0.5+age/2)
  }
  text(1986,13,"Age 0")
  text(1986,42, "Age 6")
  #
  # Weigth at age
  #
  w1985 = weight( dat$mean[dat$YEAR==1985] )
  w2016 = weight( dat$mean[dat$YEAR==2016] )
  age = 0:6
  defaultpanel(xlim=age, range(c(w1985, w2016)), 
               xlab="Age", ylab="Weight (g)")
  lines(age, w1985, lwd=3)
  lines(age, w2016, lwd=3)
  text(4,380,"1985")
  text(4,190,"2016")
  
}

compareModels <- function(W=1000) {
  F = seq(0,1, length.out = 10)
  p = baseparameters(W=1000)
  pc = paramCommunitymodel(W=1000, facRmax = 1e-8)
  
  Y = array(length(F))
  YC = Y
  
  for (i in 1:length(F)) {
    p$F = F[i]
    r = spectrum(p)
    Y[i] = calcYield(p,r)
    
    pc$F=F[i]
    rc = runspectrum(pc)
    YC[i] = rc$R$Y[100]
  }
  
  plot(F, Y/max(Y))
  lines(F, YC/max(YC), lty=dashed)
}

plotAllChapter7b = function() {
  pdfplot("Chapter7b/GrowthF0.pdf", plotGrowthF0, width=singlewidth, height=height)
  pdfplot("Chapter7b/Spectra_ConsumerResourceModel.pdf", 
          plotSpectra_ConsumerResourceModel, width=doublewidth, height = 2*height)
  pdfplot("Chapter7b/ConsumerResourceVariation.pdf",
          plotConsumerResourceVariation2, width=doublewidth,
          height=height)
  pdfplot("Chapter7b/FishingvsWf.pdf",
          plotFishingvsWf, width=doublewidth, height=height)
  pdfplot("Chapter7b/OptimalEtaF.pdf",
          plotOptimalEtaF, width=doublewidth, height=height)
  pdfplot("Chapter7b/PlaiceGrowth.pdf",
          plotPlaiceGrowth, width=doublewidth, height=height)
  plotSpatialPlaice()
}
