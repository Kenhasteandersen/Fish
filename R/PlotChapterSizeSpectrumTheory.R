#
# Plots for TeX/Chapter 2
# 

source("R/basetools.R")
source("R/basefunctions.R")
source("R/baseparameters.R")
source("R/PlotChapterIndividual.R")

dir.create("TeX/ChapterSizeSpectrumTheory")

plotSheldon <- function() {
  data <- read.csv("Data/sheldon.dat")
  data <- 10^data
  # Change the conversion from cal to g with a factor 1.3 for fish
  # (for other sizes the conversion is 1-to-1, see Boudrey & Dickie appendix)
  ixFish <- data[,1] > 0.1
  data[ixFish,] <- data[ixFish,]*1.3
  #
  # 
  d <- data.frame(x=data[,1], y=data[,2], Ecosystem=1)
  d$Ecosystem[1:8] <- 2
  d$Ecosystem[9:24] <- 1
  d$Ecosystem[25:41] <- 3
  d$Ecosystem[42:56] <- 4
  d$Ecosystem <- factor(d$Ecosystem, 
                        labels=c("Inland lakes","Lake Superior","Lake Michigan","Scotian shelf"))
  
  # -----
  # Make fits:
  #
  
  fit = lm(log(y) ~ log(x) + Ecosystem, data=d)
  exponent = fit$coefficients[2]
  
  DD = data.frame(x=d$x,y=exp(predict(fit)),Ecosystem=d$Ecosystem)
  
  fig <- ggplot()
  for (Ecosystem in unique(d$Ecosystem)) {
    fig <- fig + 
      geom_line(data=DD[d$Ecosystem==Ecosystem,], aes(x=x,y=y,colour=Ecosystem),
                linetype="dashed")
  }
  # -----
  # Plot points
  #
  d1 <- d[1:4,] #data.frame(x=data[1:4,1], y=data[1:4,2])
  d2 <- d[5:6,]
  d3 <- d[7:8,]
  
  d4 <- d[9:18,]
  d5 <- d[19:24,]
  
  d6 <- d[25:41,]
  
  d7 <- d[42:48,]
  d8 <- d[49:52,]
  d9 <- d[53:56,]
  
  fig <- fig +
    geom_line(data=d1, aes(x=x,y=y,colour=Ecosystem)) + 
    geom_point(data=d1,aes(x=x,y=y,colour=Ecosystem,shape=Ecosystem),size=2) + 
    geom_line(data=d2, aes(x=x,y=y,colour=Ecosystem)) +
    geom_point(data=d2,aes(x=x,y=y,colour=Ecosystem,shape=Ecosystem),size=2) +
    geom_line(data=d3, aes(x=x,y=y,colour=Ecosystem)) +
    geom_point(data=d3,aes(x=x,y=y,colour=Ecosystem,shape=Ecosystem),size=2)
  
  fig <- fig + 
    geom_line(data=d4,aes(x=x,y=y,colour=Ecosystem)) +
    geom_point(data=d4,aes(x=x,y=y,colour=Ecosystem,shape=Ecosystem),size=2) +
    geom_line(data=d5,aes(x=x,y=y,colour=Ecosystem)) +
    geom_point(data=d5,aes(x=x,y=y,colour=Ecosystem,shape=Ecosystem),size=2)
  
  fig <- fig + 
    geom_line(data=d6, aes(x=x,y=y,colour=Ecosystem)) + 
    geom_point(data=d6,aes(x=x,y=y,colour=Ecosystem,shape=Ecosystem),size=2)  
  
  fig <- fig + 
    geom_line(data=d7, aes(x=x,y=y,colour=Ecosystem)) + 
    geom_point(data=d7,aes(x=x,y=y,colour=Ecosystem,shape=Ecosystem),size=2) +
    geom_line(data=d8, aes(x=x,y=y,colour=Ecosystem)) + 
    geom_point(data=d8,aes(x=x,y=y,colour=Ecosystem,shape=Ecosystem),size=2) +
    geom_line(data=d9, aes(x=x,y=y,colour=Ecosystem)) + 
    geom_point(data=d9,aes(x=x,y=y,colour=Ecosystem,shape=Ecosystem),size=2)
  #
  # Annotations
  #
  fig <- fig +
    annotate("segment", x=1e-13, xend=5e-7, y=2e2, yend=2e2,
             arrow=arrow(ends="last", angle=90, length=unit(.2,"cm"))) +
    annotate("text", x=3e-10, y=4e2, label="unicellular plankton ",size=2.6) +
    #
    annotate("segment", x=1e-6, xend=1e-2, y=0.7e2, yend=0.7e2,
             arrow=arrow(ends="both", angle=90, length=unit(.2,"cm"))) +
    annotate("text", x=1e-4, y=1.3e2, label="copepods",size=2.6) +
    #
    annotate("segment", x=2e-2, xend=1e6, y=0.4e2, yend=0.4e2,
             arrow=arrow(ends="first", angle=90, length=unit(.2,"cm"))) +
    annotate("text", x=100, y=0.8e2, label="fish",size=2.6) 
  
  
  #
  # Finalize figure
  #
  fig <- loglog(fig, xlim=c(1e-13,1e6), ylim=c(0.01,500), ixStep=2) +
    scale_colour_grey()+
    scale_fill_grey()+
    xlab("Body weight (g)") + 
    ylab(expression(Biomass~(g/m^2))) +
    coord_fixed(ratio=2) #+
  #scale_x_log10(breaks = 10^seq(from=-14,by=2,to=10),
  #              labels = trans_format("log10", math_format(10^.x)))
  
  ggsave("TeX/ChapterSizeSpectrumTheory/sheldon.pdf", width=14, height=5, units="cm")
  fig
}

plotClearance <- function() {
  fig <- ggplot() 
  
  #
  #  data from Acuna et al (2011)
  #
  #datA <- read.csv('Clearance_rate Acuna.csv', header=TRUE)
  #datA$type <- "0"
  #datA$type[1:418] <- "Acuna fish"
  #datA$type[783:1449] <- "Acuna crustaceans"
  #datA <- na.omit(datA[datA$type != "0",])
  #
  # Data from TK
  #
  #datTK <- na.omit( read.csv("Clearance_rate Thomas.csv", header=TRUE, sep=";"))
  #datTK$type = 0
  #datTK$type[1:79] <- "TKunicellular"
  #datTK$type[79:127] <- "TKfish"
  #
  # Data from TK & Hirst
  # Ki??rboe, T; Hirst, AG (2013): Compilation of maximum ingestion and maximum clearance rate data for marine pelagic organisms. doi:10.1594/PANGAEA.819856,
  # In Supplement to: Ki??rboe, Thomas; Hirst, Andrew G (2014): Shifts in mass-scaling of respiration, feeding, and growth rates across life-form transitions in marine pelagic organisms. American Naturalist, 183(4), E118-E130, doi:10.1086/675241
  #datTK2 <- read.csv("feeding_rates_compilation_Kiorboe_Hirst.tab.tsv", header=TRUE, sep="\t")
  # Remove Acuna's data
  #datTK2 <- datTK2[grep("Acuna", dd$Reference, invert=TRUE, ignore.case=TRUE),]
  #datTK2$type="Kiorboe & Hirst"
  #
  # Data from TK & hirst again
  #
  datTK3 <- read.csv("Data/TK Appendix feeding rates - revised.csv",header=TRUE,sep=";")
  #The factor 8 is a good rough estimate to get "wet weight" from carbon weight
  tmp4 <- data.frame(w=8e-3*datTK3$Body.mass, beta=365*24*0.001*datTK3$Fmax.1, type=datTK3$Group)
  tmp4 <- na.omit(tmp4)
  tmp4$type[tmp4$type=="Tunicates"]="Tunicata"
  tmp4$type[tmp4$type=="Ciliate"]="Ciliates"
  tmp4$type[tmp4$type=="Dinoflagellate"]="Dinoflagellates"
  
  
  #
  # Merge data sets
  #
  #tmp1 <- data.frame(w=datA$WW, beta=datA$beta, type=datA$type)
  #tmp2 <- data.frame(w=datTK$Body.carbon*8, beta=datTK$beta/1000, type=datTK$type)
  #
  #tmp3 <- data.frame(w=8e-6*datTK2$Biom.C.ind....g..., beta=24*0.001*datTK2$CR..ml...h., type=datTK2$type)
  #tmp3 <- na.omit(tmp3)
  #tmp3 <- tmp3[tmp3$w>0,] # get rid of zero values
  #
  # Focus on fish and "others"
  #
  tmp = 0*tmp4$w
  tmp[tmp4$type=="Pisces"] = 1
  tmp4$type2 <- factor(tmp,labels=c("Other","Fish"))
  
  dat <- tmp4 #rbind(tmp4,tmp1)
  
  fig <- fig + geom_point(data=dat, aes(x=w, y=beta, colour=type2))
  #
  # fit all
  # 
  fit <- lm(log(beta) ~ log(w), data=dat)
  fit_q <- fit$coefficients[2]
  fit_gamma <- exp(fit$coefficients[1])
  tmp <- data.frame(x=dat$w, y=exp(predict(fit)))
  fig <- fig + geom_line(data=tmp, aes(x=x,y=y),linetype="dashed")
  #fig <- fig + annotate("text",x=1e-8,y=1e4,parse=TRUE,label="'V(w) = ' *900*w^{0.95}")
  #fig <- fig + annotate("text", x=2e-4, y=5e5, label="beta== asdf", parse=TRUE)
  # 
  # fit  fish
  # 
  ix <- dat$type=="Pisces"
  fitfish = lm(log(beta[ix]) ~ log(w[ix]), data=dat)
  fitfish_q <- fitfish$coefficients[2]
  fitfish_gamma <- exp(fitfish$coefficients[1])
  tmp <- data.frame(x=dat$w[ix], y=exp(predict(fitfish)))
  fig <- fig + geom_line(data=tmp, aes(x=x,y=y),linetype="dotted")
  
  canonical_fish_gamma <- exp(mean(log(dat$beta[ix]/ (dat$w[ix]^0.8))))
  gammastd <- sqrt(exp(var(log(dat$beta[ix]/ (dat$w[ix]^0.8)))))
  fig <- fig + geom_line(data=dat[ix,], aes(x=w, y=canonical_fish_gamma*w^0.8))
  
  canonical_gamma <- exp(mean(log(dat$beta/ (dat$w^0.8))))
  
  cat("All gamma", fit_gamma, "\n")
  cat("All q", fit_q, "\n")
  
  cat("Fish gamma", fitfish_gamma, "\n")
  cat("Fish q", fitfish_q, "\n")
  
  cat("Canonical gamma", canonical_gamma, "\n")
  #-------------------
  # Finalize figure
  #
  fig <- loglog(fig, iyStep=2) +
    theme(legend.position=c(1,0), legend.justification=c(1.1,-0.1)) +
    scale_colour_grey()+
    #                    labels=c("","","","")) +
    #scale_fill_brewer(palette="PuBu") +
    xlab("Body weight (g)") + 
    ylab("Clearance rate (liter/yr)") +
    labs(colour="",shape="")
  
  # fig <- fig +
  #   theme_classic(base_size=10) +
  #   theme(legend.position=c(1,0), legend.justification=c(1,0)) +
  #   scale_colour_grey()+
  #   #                    labels=c("","","","")) +
  #   #scale_fill_brewer(palette="PuBu") +
  #   xlab("Body weight (g)") + 
  #   ylab("Clearance rate (l/yr)") +
  #   labs(colour="",shape="") +
  #   scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #                 labels = trans_format("log10", math_format(10^.x)),
  #                 oob=squish) +
  #   scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #                 labels = trans_format("log10", math_format(10^.x)),
  #                 oob=squish) +
  #   #guides(colour=FALSE) +
  #   annotation_logticks(long=unit(1,"mm"),mid=unit(0,"cm"), short=unit(0,"cm")) + # l = left, b = bottom etc) +
  #   coord_fixed(ratio=1)
  
  fig
  
  ggsave("TeX/ChapterSizeSpectrumTheory/Clearance_rate.pdf", width=1.5*singlewidth, height=height)
}


plotImax <- function() {
  xlim = c(1e-11,1e6)
  ylim = c(1e-10,1e6)
  quartz(width = doublewidth, height=height)
  # =============================================================
  #  Respiration
  #  
  
  #------------------------
  # Respiration Data from TK & hirst again
  #
  dat <- read.csv("Data/Appendix respiration rates -revised.csv",header=TRUE,sep=";")
  #
  # Conversion is: 1 ml O2 corresponds to 0.5 mg C:
  # 1 ml O2 = 1.33 mg O2.
  # Weight O2: 32 g/mol
  # Weight C:  12 g/mol
  # 1 mg O2 is therefore 1.33/32 mol O2.  This corresponds to 1 mol C (to produce CO2), so 
  # 1 ml O2 = 1.33/32*12 = 0.5 mgC.
  # Note: do not convert this to ww, because water is supposed to be lost at the same rate
  dat <- data.frame(w=8e-3*dat$Body.mass, resp=1e-6*0.5*24*365*dat$Spec.Resp.15C*dat$Body.mass, type=dat$Group)
  dat <- na.omit(dat)
  dat <- dat[dat$resp>0,]
  #
  # Focus on fish and "others"
  #
  tmp = 0*dat$w
  tmp[dat$type=="Pisces"] = 1
  dat$type2 <- factor(tmp,labels=c("Other","Fish"))
  
  fig2 <- ggplot()
  fig2 <- fig2 + geom_point(data=dat, aes(x=w, y=resp, colour=type2))
  
  fit <- lm(log(resp) ~ log(w), data=dat)
  fit_kr <- exp(fit$coefficients[1])
  fit_kr_p <- fit$coefficients[2]
  tmp <- data.frame(x=dat$w, y=exp(predict(fit)))
  fig2 <- fig2 + geom_line(data=tmp, aes(x=x,y=y), linetype="dashed", color="white")
  
  cat("kr ", fit_kr, "p ", fit_kr_p, "\n")
  
  ix <- dat$type=="Pisces"
  fitfish_kr <- lm(log(resp[ix]) ~ log(w[ix]), data=dat)
  tmp <- data.frame(x=dat$w[ix], y=exp(predict(fitfish_kr)))
  fig2 <- fig2 + geom_line(data=tmp, aes(x=x,y=y),linetype="dotted", color="white")
  canonical_kr <-  exp(mean(log(dat$resp/ (dat$w^0.75))))
  
  canonical_fish_kr <- exp(mean(log(dat$resp[ix]/ (dat$w[ix]^0.75))))
  krstd <- sqrt(exp(var(log(dat$resp[ix]/ (dat$w[ix]^0.75)))))
  fig2 <- fig2 + geom_line(data=dat[ix,], aes(x=w, y=canonical_fish_kr*w^0.75), color="white")
  
  cat("kr canonical ", canonical_kr, "\n")
  cat("kr fish ", exp(fitfish_kr$coefficients[1]), "p ",fitfish_kr$coefficients[2],"\n")
  cat("kr fish canonical ", canonical_fish_kr, "\n")
  
  #-------------------
  # Finalize respiration panel
  #
  fig2 <- loglog(fig2, ixStep=3, iyStep=2, xlim=xlim, ylim=ylim) +
    theme(legend.position=c(.95,0.05), legend.justification=c(1,0)) +
    scale_colour_grey() +
    xlab("Body weight (g)") + 
    ylab("Respiration (g/yr)") +
    labs(colour="",shape="") +
    coord_fixed(ratio=1)
  
  # ==========================================================
  # Imax
  # 
  tmpImax <- read.csv("Data/TK Appendix feeding rates - revised.csv",header=TRUE,sep=";")
  datImax <- data.frame(w=8e-3*tmpImax$Body.mass, Imax=8*0.000001*24*365*tmpImax$Specific.Imax*tmpImax$Body.mass, type=tmpImax$Group)
  datImax <- na.omit(datImax)
  datImax <- datImax[datImax$Imax>0,]
  #
  # Focus on fish and "others"
  #
  tmp = 0*datImax$w
  tmp[datImax$type=="Pisces"] = 1
  datImax$type2 <- factor(tmp,labels=c("Other","Fish"))
  
  #source("R/PlotAllTeX/Chapter3.R")
  dat3 <- getGrowthParameters()
  epsa <- 0.6
  f0 <- 0.6
  dat3$h <- dat3$A/(0.6*0.4) # (dat3$A+canonical_kr)/(epsa*f0)
  # add to Imax data:
  w = 0.25*dat3$W  # use the data point around size at maturation
  tmp <- data.frame(w=w, Imax=dat3$h * (0.25*dat3$W)^0.75, 
                    type="Fish",type2="Fish")
  datImax <- rbind(datImax, tmp)

  
  fig <- ggplot() 
  fig <- fig + geom_point(data=datImax, aes(x=w, y=Imax, colour=type2))
  #
  # Fits
  # 
  fit <- lm(log(Imax) ~ log(w), data=datImax)
  #fit_n <- fit$coefficients[2]
  fit_h <- exp(fit$coefficients[1])
  tmp <- data.frame(x=datImax$w, y=exp(predict(fit)))
  fig <- fig + geom_line(data=tmp, aes(x=x,y=y),linetype="dashed", color="white")
  
  cat('h ', fit_h, ' n ',fit$coefficients[2],"\n")
  canonical_h <- exp(mean(log(datImax$Imax / datImax$w^0.75)))
  cat("Canonical h ", canonical_h, "\n")
  
  ix <- datImax$type2=="Fish"
  fitfish_h <- lm(log(Imax[ix]) ~ log(w[ix]), data=datImax)
  tmp <- data.frame(x=datImax$w[ix], y=exp(predict(fitfish_h)))
  fig <- fig + geom_line(data=tmp, aes(x=x,y=y), linetype="dotted", color="white")
  
  canonical_fish_h <- exp(mean(log(datImax$Imax[ix]/ (datImax$w[ix]^0.75))))
  hstd <- sqrt(exp(var(log(datImax$Imax[ix]/ (datImax$w[ix]^0.75)))))
  fig <- fig + geom_line(data=data.frame(x=c(min(datImax$w[ix]), max(datImax$w))), 
                         aes(x=x, y=canonical_fish_h*x^0.75), color="white")
  fc_fish <- canonical_fish_kr/(0.6*canonical_fish_h)
  
  cat("Fish h ", exp(fitfish_h$coefficients[1]), " p ", fitfish_h$coefficients[2], "\n")
  cat("Canonical fish h ", canonical_fish_h, "\n")
  cat("Canonical fish kr ", canonical_fish_kr, "\n")
  cat("fc fish", fc_fish,"\n")
  cat("fc canonical", canonical_kr/(0.6*canonical_h), "\n")
  
  #-------------------
  # Finalize Imax panel
  #
  fig <- loglog(fig, ixStep=3, iyStep=2, xlim=xlim, ylim=ylim) +
    theme(legend.position=c(0.95,0.05), legend.justification=c(1,0)) +
    scale_colour_grey()+ 
    xlab("Body weight (g)") + 
    ylab("Max. ingestion rate (g/yr)") +
    labs(colour="",shape="") +
    guides(colour=FALSE) +
    coord_fixed(ratio=1)
  
  
  # =======================================
  #  Assemble final figure
  # 
  
  grid <- plot_grid(fig2,fig, ncol=2, align="h")
  ggsave("TeX/ChapterSizeSpectrumTheory/RespirationAndImax.pdf", grid, width=doublewidth, height=doublewidth/2)
}

plotRes2<-function() {
  defaultplothorizontal(2)
  #
  # Respiration
  #
  dat <- read.csv("TeX/ChapterSizeSpectrumTheory/Appendix respiration rates -revised.csv",header=TRUE,sep=";")
  dat$w <- 8e-3*dat$Body.mass
  dat$resp <- 1e-6*0.5*24*365*dat$Spec.Resp.15C*dat$Body.mass
  dat <- data.frame(dat)
  #dat <- na.omit(dat)
  dat <- dat[dat$resp>0 & !is.na(dat$resp),]
  #
  # Focus on fish and "others"
  #
  ixfish <- dat$Group=="Pisces" & !is.na(dat$resp)
  
  # Acuna data:
  #datA <- read.csv("TeX/ChapterSizeSpectrumTheory/Acuna_et_al_2011_respiration.csv", header=TRUE, sep=";")
  #datA$w <- datA$WW
  # conv: Cww/C * gC/mol * C/O * milli * days/year
  #Q10 = 1.83
  #datA$resp <- datA$Respiration * 8 * 12 * 0.5 * 1e-3 *365 *Q10^((15-datA$Temp)/10)
  #ixfishA <- datA$group=="F"
  #dat <- datA
  
  loglogpanel(xlim=c(1e-10, 1e5), ylim=c(1e-10,1e6), label = TRUE)
  points(dat$w[!ixfish], dat$resp[!ixfish], pch=dots, cex=0.3)
  points(dat$w[ixfish], dat$resp[ixfish], pch=dots, col=stdgrey, cex=0.3)
  
  fit_all = 
  
  #
  # Cmax
  #
  tmpImax <- read.csv("Data/TK Appendix feeding rates - revised.csv",header=TRUE,sep=";")
  datImax <- data.frame(w=8e-3*tmpImax$Body.mass, Imax=8e-6*24*365*tmpImax$Specific.Imax*tmpImax$Body.mass, type=tmpImax$Group)
  datImax <- na.omit(datImax)
  datImax <- datImax[datImax$Imax>0,]
  #
  # Focus on fish and "others"
  #
  tmp = 0*datImax$w
  tmp[datImax$type=="Pisces"] = 1
  datImax$type2 <- factor(tmp,labels=c("Other","Fish"))
  
  # fish:
  dat2 <- getGrowthParameters()
  dat2$h <- dat2$A/(0.6*0.4)
  datImax <- rbind(datImax, 
                   data.frame(w=0.25*dat2$W, Imax=dat2$h*(0.25*dat2$W)^0.75, 
                    type="Fish",type2="Fish"))
  ixfish <- datImax$type2=="Fish"
  
  points(datImax$w[!ixfish], datImax$Imax[!ixfish], pch=dots, col=stdgrey)
  points(datImax$w[ ixfish], datImax$Imax[ ixfish], pch=squares, col=stdgrey)
  
  canonical_fish_h <- exp(mean(log(dat2$h)))
  canonical_fish_kr <- exp(mean(log(dat$resp[ixfish]/ (dat$w[ixfish]^0.75))))
  
  fc <- canonical_fish_kr/(0.6*canonical_fish_h)
  
  cat("Canonical h ", canonical_fish_h, "\n")
  cat("Canonical kr ", canonical_fish_kr, "\n")
  cat("fc ", fc, "\n")
}


plotUrsin <- function() {
  p <- baseparameters()
  beta = p$beta
  sigma = p$sigma
  
  phi  <- function(xp) exp( -(log(1/(beta * xp)))^2/(2 * sigma^2))
  fig <- ggplot(data.frame(x=c(1e-5,2)), aes(x=x)) +
    stat_function(fun = phi ) 
  
  phi2  <- function(xp) exp( -(log(1/(1.7*beta * xp)))^2/(2 * sigma^2))
  fig <- fig + stat_function(fun = phi2, linetype="dashed") 
  #
  # Add annotations:
  # 
  fig <- fig +
    annotate("segment", x=1/beta, xend=1, y=1.05, yend=1.05,
             arrow=arrow(ends="both", angle=90, length=unit(.2,"cm"))) +
    annotate("text", x=beta^(-0.5), y=1.1, label="beta",size=3,parse=TRUE) 
  
  fig <- fig +
    annotate("segment", x=1/(1.7*beta), xend=1, y=1.23, yend=1.23,
             arrow=arrow(ends="both", angle=90, length=unit(.2,"cm"))) +
    annotate("text", x=beta^(-0.5), y=1.3, label="beta[PPMR]",size=3,parse=TRUE) 
  
  
  
  #
  # Finalize:
  # 
  fig <- semilogx(fig)+
    xlab("Prey:predator weight ratio") + 
    ylab("Preference") 
  # fig <- fig +
  #   theme_classic(base_size=10) +
  #   theme(legend.position=c(1,0), legend.justification=c(1,0)) +
  #   xlab("Prey:predator weight ratio") + 
  #   ylab("Preference") +
  #   labs(colour="",shape="") +
  #   scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #                 labels = trans_format("log10", math_format(10^.x)),
  #                 oob=squish) +
  #   annotation_logticks(long=unit(1,"mm"),mid=unit(0,"cm"), short=unit(0,"cm"))
  
  fig
  
  ggsave("TeX/ChapterSizeSpectrumTheory/Ursin.pdf", width=1.5*singlewidth, height=height)
  
}

plotPredPrey <- function() {
  dat <- read.csv("Data/Predator_and_prey_body_sizes_in_marine_food_webs_vsn4.dat",header=TRUE,sep="\t")
  # Get only fish:
  ixFish <- dat$Predator..taxon=="ectotherm vertebrate"
  dat <- dat[ixFish,]
  # Fix prey measured in mg:
  ix = dat$Prey.mass.unit=="mg"
  dat$Prey.mass[ix] = dat$Prey.mass[ix]/1000
  
  # --------------
  fig <- ggplot()
  
  fig <- fig + 
    geom_point(data=dat[seq(1, dim(dat)[1],by=4),], aes(x=Predator.mass, y=Prey.mass),
               size=.5, colour="grey",fill="grey", alpha=0.2)
  fig
  #
  # Fit
  # 
  fit <- lm(log(Prey.mass) ~ log(Predator.mass), data=dat)
  xx = c(1e-4, 1e6)
  tmp <- data.frame(x=xx, y=exp(fit$coefficients[1])*xx^fit$coefficients[2])
  fig <- fig + geom_line(data=tmp, aes(x=x,y=y), linetype="dashed")
  
  ix = which(dat$Predator.mass>1 & dat$Predator.mass<1e5)
  fitLarge <- lm(log(Prey.mass[ix]) ~ log(Predator.mass[ix]), data=dat)
  tmp <- data.frame(x=xx, y=exp(fitLarge$coefficients[1])*xx^fitLarge$coefficients[2])
  fig <- fig + geom_line(data=tmp, aes(x=x,y=y), linetype="dotted")
  
  PPMR <- exp(mean(log(dat$Predator.mass[ix]/ dat$Prey.mass[ix])))
  PPMRstd <- sqrt(exp(var(log(dat$Predator.mass[ix]/ dat$Prey.mass[ix]))))
  cat('PPMR fish:', PPMR, 'std:', PPMRstd,'\n')
  fig <- fig + geom_line(data=data.frame(x=c(min(dat$Predator.mass[ix]), max(dat$Predator.mass))), 
                         aes(x=x, y=x/PPMR))
  
  
  
  tmp <- data.frame(x=c(1e-5,1e4,1e-5), y=c(1e-5,1e4,1e4))
  fig <- fig + 
    geom_polygon(data=tmp, aes(x=x,y=y),fill=grey(0.2), colour=grey(0.2))
  
  #
  # Finalize:
  # 
  fig <- loglog(fig, ixStep = 2, iyStep = 2, bExpand=TRUE) +
    xlab("Predator weight (g)") + 
    ylab("Prey weight (g)") +
    labs(colour="",shape="") 
  
  # fig <- fig +
  #   theme_classic(base_size=10) +
  #   theme(legend.position=c(1,0), legend.justification=c(1,0)) +
  #   xlab("Predator weight (g)") + 
  #   ylab("Prey weight (g)") +
  #   labs(colour="",shape="") +
  #   scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #                 labels = trans_format("log10", math_format(10^.x)),
  #                 oob=squish,limits=c(1e-4,1e6)) +
  #   scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #                 labels = trans_format("log10", math_format(10^.x)),
  #                 oob=squish,limits=c(1e-9,1e4)) +
  #   annotation_logticks(long=unit(1,"mm"),mid=unit(0,"cm"), short=unit(0,"cm")) + # l = left, b = bottom etc) +
  #   coord_fixed(ratio=1)
  
  
  ggsave("TeX/ChapterSizeSpectrumTheory/PredPrey.pdf", width=singlewidth, height=singlewidth)

  fig
}

plotCalculations <- function() {
  #
  # Calculate size spectra from first principles
  #
  sigma = 1
  #
  # canonical
  # 
  q = 0.8
  gamma = 120556.4  # l/yr; from plotClearance()
  n = 0.75
  h = 16 # from plotImax()
  f0 = 0.6
  alpha = 0.6
  fc <- 0.12 # from plotImax()
  lambda = 2-n+q
  PPMRratio = exp((lambda-3/2)*sigma^2)
  beta = 708/PPMRratio
  Phi = beta^(lambda-2)*sigma*exp((lambda-2)^2*sigma^2/2)*sqrt(2*pi)
  kappa = f0 * h/(gamma*Phi) 
  corr = (10^(2-lambda)-1) / (2-lambda)  # For comparison with B&D Sheldon spectra
  Phi_p = beta^(2*n-q-1) * exp((2*n-q-1)*(q-1)*sigma^2/2)
  a = Phi_p * f0/(alpha*(f0-fc))
  eT = beta^(2*n-q-1)
  cat("\n-------------\nCanonical","\n")
  cat("gamma",gamma,"\n")
  cat("h",h,"\n")
  cat("fc",fc,"\n")
  cat("beta",beta,"\n")
  cat("Kappa", kappa,"\n")
  cat("Phi_a", Phi,"\n")
  cat("Phi_p", Phi_p, "\n")
  cat("eT", eT, "\n")
  #
  # Fit fish
  #
  cat("\nFit fish","\n")
  cat("q",0.76,"\n")  # from plotClearance
  cat("gamma",607009,"\n")  # from plotClearance
  cat("n",0.76,"\n") # from plotImax
  cat("h",17.6,"\n")# from plotImax
  cat("fc",0.15,"\n")# from plotImax
  cat("kr", 1.67,"\n")# from plotImax
  #
  # All:
  # 
  q = 0.96 # 0.96
  gammaFit = 330629 # 3.3e5
  n = 0.77 #0.8527
  hFit = 17.23 # 42.6 # 24.5
  krFit = 1.65
  lambdaFit = 2-n+q
  PPMRratio = exp((lambdaFit-3/2)*sigma^2)
  beta = 1224/PPMRratio
  fc <- krFit/(alpha*hFit)
  PhiFit = beta^(lambdaFit-2)*sigma*exp((lambdaFit-2)^2*sigma^2*sqrt(2*pi)/2)
  kappaFit = f0*hFit/(gammaFit*PhiFit) 
  corrFit = (10^(2-lambdaFit)-1) / (2-lambdaFit)  # For comparison with B&D Sheldon spectra
  Phi_pFit = beta^(2*n-q-1) * exp((2*n-q-1)*(q-1)*sigma^2/2)
  aFit = Phi_p * f0/(alpha*(f0-fc))
  eTfit = beta^(2*n-q-1)
  cat("\n-------------\nAll","\n")
  cat("q",q,"\n")
  cat("gamma",gammaFit,"\n")
  cat("n",n,"\n")
  cat("h",hFit,"\n")
  cat("lambda",lambdaFit,"\n")
  cat("kappa",kappaFit,"\n")
  cat("beta",beta,"\n")
  cat("fc",fc,"\n")
  cat("Phi_a", PhiFit,"\n")
  cat("Phi_p", Phi_pFit, "\n")
  cat("eT", eTfit, "\n")
  
  #
  # Plot
  # 

    #
  # Load the B & D data to plot on top
  #
  data <- read.csv("Data/sheldon.dat")
  data <- 10^data
  # Change the conversion from cal to g with a factor 1.3 for fish
  # (for other sizes the conversion is 1-to-1, see Boudrey & Dickie appendix)
  ixFish <- data[,1] > 0.1
  data[ixFish,] <- data[ixFish,]*1.3
  d <- data.frame(x=data[,1], y=data[,2], series=1)
  d$series[1:8] <- 1
  d$series[9:24] <- 0
  d$series[25:41] <- 2
  d$series[42:56] <- 3
  #d$series <- factor(d$series,labels=c("Inland lakes","Lake Superior","Lake Michigan","Scotian shelf"))
  
  defaultplot()
  par(mar=par()$mar + c(0,0,0,10)) # space for legend
  loglogpanel(xlim=d$x, ylim=c(1e-7, 1e-2),
              xlab="Weight (g)", ylab="Sheldon spectrum (g/L)",
              powx = seq(-12, 4, by = 2))
  for (series in unique(d$series))
    points(d$x[d$series==series], d$y[d$series==series]/100/300, # first 100 converts from m^2 to dm^2. 
           # The last 300 converts to l, assuming a depth of 30 m
           pch=c(dots, triangles, squares, 3)[series+1],
           col=c("black", darkgrey, stdgrey, lightgrey)[series+1])
  lines(sort(d$x), kappaFit*sort(d$x)^(2-lambdaFit), lty=dashed, lwd=3)
  lines(d$x, kappa*d$x^(2-lambda), lwd=3)
  
  legend("right", bty="n", inset=c(-0.5,0), xpd=TRUE,
         legend=c("Inland lakes","Lake Superior","Lake Michigan","Scotian shelf"),
         pch=c(dots, triangles, squares, 3), 
         col=c("black", darkgrey, stdgrey, lightgrey),
         title="")
  

    # Figure for presentation:
  #d$series <- factor(d$series,labels=c("Inland lakes","Lake Superior","Lake Michigan","Scotian shelf"))
  #fig2 <- ggplot(data.frame(x=c(1e-13, 1e6)), aes(x)) + 
  #  #geom_point(data=d, aes(x=x,y=y/100/1000,shape=series,colour=series)) +
  #  geom_point(data=subset(d, series=="Scotian shelf"),aes(x=x,y=y/100/1000),colour="grey",size=2) +
  #  stat_function(fun=function(x) kappa*x^(2-lambda), geom="line", size=2) 
  
  
  #fig2 <- loglog(fig2, ixStep=5, iyStep=5, ylim=c(5e-8, 1e-3)) + 
  #  theme(legend.position="none") +
  #  scale_colour_grey() +
  #  xlab("Weight (g)") + 
  #  ylab("Sheldon spectrum (g/L)") +
  #  labs(colour="",shape="")
  #fig2
  #ggsave("TeX/ChapterSizeSpectrumTheory/PresentationSheldon.pdf", width=10, height=7, units="cm")
  
}


plotMortality <- function() {
  #
  # Read mortality of fish:
  # 
  A=read.csv("Data/Gislason.csv")
  names(A) <- c("L","Linf","K","M","T")
  A$w <- 0.01*A$L^3
  #
  # Correct to 15 degrees C using Q10 = 2
  #
  Q10 = 1.83
  A$M <- A$M * Q10^((15-A$T)/10)
  fig <- ggplot() + 
    geom_point(data=A, aes(x=w, y=M), colour="grey")
  #
  # Read mortality of copepods:
  # 
  B = read.csv("Data/HirstMortality.csv",header=TRUE, sep=";")
  fig <- fig + geom_point(data=B, aes(x=8*1e-6*Weight.mugDW, y=365*0.5*(MinMortality+MaxMortality)))
  #
  # Fit
  # 
  fit <- lm(log(M) ~ log(w), data=A)
  tmp <- data.frame(x=A$w, y=exp(predict(fit)))
  fig <- fig + geom_line(data=tmp, aes(x=x,y=y), linetype="dashed")
  
  tmp <- data.frame(x=8*c(1e-7,1e-4))
  tmp$y = 365*exp(0.047*15 -0.154* (log(1e6*tmp$x/8))-2.532)
  fig <- fig + geom_line(data=tmp, aes(x=x,y=y), linetype="dotted")
  
  # add a line with metabolic scaling:
  tmp <- data.frame(x=8*c(1e-8,1e5))
  tmp$y = 2.66*tmp$x^-0.25
  fig <- fig + geom_line(data=tmp, aes(x=x,y=y))
  
  
  #-------------------
  # Finalize figure
  #
  fig <- loglog(fig, ixStep = 1) +
    scale_colour_grey() +
    xlab("Body weight (g)") + 
    ylab("Mortality rate (1/yr)") +
    labs(colour="",shape="")
  
  ggsave("TeX/ChapterSizeSpectrumTheory/Mortality.pdf", width=1.5*singlewidth, height=1.5*height)
  fig
}

SheldonSketch <- function() {
  w = c(1,1e4)
  defaultplothorizontal(3)
  loglogpanel(xlim=w, ylim=c(1,1e-8),
              ylab="Biomass (biomass)", label=TRUE)
  lines(w, w^-0.05, lwd=3)
  
  loglogpanel(xlim=w, ylim=c(1,1e-8),
              ylab="Biomass spectrum (biomass/weight)", label=TRUE)
  lines(w, w^-1.05, lwd=3)
  
  loglogpanel(xlim=w, ylim=c(1,1e-8),
              ylab="Number spectrum (numbers/weight)", label=TRUE)
  lines(w, w^-2.05, lwd=3)
}


plotAllChapterSizeSpectrumTheory <- function() 
{
  plotSheldon()
  plotClearance()
  plotImax()
  #pdfplot("TeX/ChapterSizeSpectrumTheory/PredPrey.pdf", plotPredPrey, width=singlewidth, height=singlewidth)
  plotPredPrey()
  plotUrsin()
  pdfplot("TeX/ChapterSizeSpectrumTheory/Calculations.pdf", plotCalculations, width=5, height=height)
  plotMortality()
}
