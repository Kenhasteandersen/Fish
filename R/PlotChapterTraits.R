library(fishsizespectrum)
source("R/basetools.R")

dir.create("TeX/ChapterTraits")

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


plotTraits <- function(bWithNames = FALSE) {
  data <- getOlssonData()
  # -----
  # A
  # ----
  #
  # Data from Kooijman and Gislason:
  #
  #data2 = read.csv("TeX/Chapter3/KooijmanGislasonPruned.csv", header=FALSE, col.names=c("K","Linf","T"))
  data2 = read.csv("Data/KooijmanGislasonPruned.csv", sep=',')
  data2$Phylum = "t" # all teleost data
  data2$A <- calcA(data2$K, data2$Linf)
  #
  # Combine data
  # 
  Q10 <- 1.83
  dat <- data.frame(A=data$A * Q10^((15-data$Temp)/10), W=data$Winf, Phylum=data$Phylum, Name=data$species)
  dat2 <- data.frame(A=data2$A * Q10^((15-data2$T)/10), W=0.01*data2$Linf^3, Phylum=data2$Phylum, Name=data2$Name)
  dat = rbind(dat, dat2)
  #
  # Population growth rates:
  #
  param <- baseparameters()
  calcr <- function(A,W,epsR,w0) {
    param$A <- A
    param$epsR <- epsR
    param$w0 <- w0
    calc_rana2(W,param)
  }

  dat$e <- 0
  dat$w0 <- 0
  ixT <- dat$Phylum == "t";
  ixE <- dat$Phylum == "c";
  dat$e[ixT] <- param$epsR
  dat$w0[ixT] <- param$w0
  dat$e[ixE] <- 0.2
  dat$w0[ixE] <- param$ElasmobranchEggMassRatio*dat$W[ixE]
  dat$r <- calcr(dat$A, dat$W, dat$e, dat$w0)
  dat$pch <- 21 + 3*ixE #[ixT] <- 0
  #dat$pch[ixE] <- 2
  #
  # Plot
  #
  loglogpanel(xlim=dat$W, ylim=c(1, 18), xlab="Asymptotic weight (g)", ylab="Growth coefficient (g$^{0.25}$yr$^{-1}$)")
  for (i in 1:length(dat$A)) {
    points(dat$W[i], dat$A[i], cex=2*sqrt(dat$r[i]), pch=dat$pch[i], col=gray(0.4), bg=gray(0.5, alpha=0.5) )
    if (bWithNames) {
      points(dat$W[i], dat$A[i], cex=2*sqrt(dat$r[i]), pch=dat$pch[i], col="black")#, bg="Grey")
      text(dat$W[i], dat$A[i], labels=dat$Name[i], cex=0.5)
    }
  }
  
}

plotTraitsFishbase <- function() {
  #
  # Download fishbase
  #
  require(rfishbase)
  require(taxize)
  if (bRecalcExpensiveFunctions) {
    fish_list <- species_list()
    flatfish_list <- species_list(Order="Pleuronectiformes")
    elasmobranch_list <- species_list(Class="Elasmobranchii")
    acti_list <- species_list(Class="Actinopterygii")
    clupeid_list <- species_list(Order="Clupeiformes")
    fish <- species(fish_list)
    growth <- popgrowth(fish_list)
    
    save(
         fish_list, flatfish_list, clupeid_list, elasmobranch_list, acti_list, fish, growth, 
         file='Data/fishbasedata.Rdata')
  }
  #
  # .. or load it:
  #
  load('Data/fishbasedata.Rdata')
  # Keep on good quality data
  growth <- growth[which(growth$to>-1 & growth$to<1),]
  
  # create indices for: teleosts and elasmobranchs
  ixActi <- ( growth$sciname %in% acti_list )
  ixElasmo <- which( growth$sciname %in% elasmobranch_list )
  
  Loo <- growth$Loo
  Woo <- growth$Winfinity
  A <- calcA(growth$K, Loo)
  cat("Geom mean A: ", exp(mean(log(A[ixActi]))),".\n")
  
  popgrowth <- function(W,A,epsR,w0) {
    n = 3/4
    a = baseparameters()$a
    epsEgg = baseparameters()$epsEgg
    A*(1-n) * (W^(1-n) - w0^(1-n))^(-1) * ( (1-a)*log(W/w0)+log(epsEgg*epsR))
  }
  
  # Growth for fish assuming Woffspring = 1 mg
  r <- popgrowth(weight(Loo), A, baseparameters()$epsR, 0.001)
  # Growth for Elasmobranchs assuming Woffspring = 0.0044*Winf and eR=0.2:
  r[ixElasmo] <- popgrowth(weight(Loo[ixElasmo]), A[ixElasmo], 0.2, 0.0044*weight(Loo[ixElasmo]))
  # Set those with negative growth to a small number, just to have the points in:
  r[r<0] <- 0.1
  #
  # Trait space
  #
  factor = 1
  #loglogpanel(xlim=c(4,500), ylim=c(0.5, 150), bExponential=FALSE,
  Woo = weight(Loo)
  loglogpanel(xlim=c(.7,1e6), ylim=c(0.5, 150), bExponential=TRUE,
                          xlab='Asymptotic weight $\\textit{W}_{\\infty}$ (g)', 
              ylab='Growth rate, $\\textit{A}$ (g$^{0.25}yr^{-1})')
  points(Woo[ixActi], A[ixActi], pch=dots, cex=factor*sqrt(r[ixActi]), col=gray(0.3, alpha=0.4))
  points(Woo[ixElasmo], A[ixElasmo], pch=triangles, cex=2*factor*sqrt(r[ixElasmo]), col="white" )
  points(Woo[ixElasmo], A[ixElasmo], pch=triangles, cex=factor*sqrt(r[ixElasmo]), col="black" )
  
    
  # Highlight Coryphaena hippurus:
  ixCoryphaena <- growth$sciname == "Coryphaena hippurus"
  points(Woo[ixCoryphaena], A[ixCoryphaena], pch=dots, cex=factor*sqrt(r[ixCoryphaena]), col="blue")
  
  # Highlight Cod
  ixCod <- growth$sciname == "Gadus morhua"
  points(Woo[ixCod], A[ixCod], pch=dots, cex=sqrt(r[ixCod]), col="orange")
  
  # Highlight tiger shark
  ixTiger <- growth$sciname == "Galeocerdo cuvier"
  points(Woo[ixTiger], A[ixTiger], pch=dots, cex=sqrt(r[ixTiger]), col="yellow")
  
  # Highlight herring
  #ixHerring <- growth$sciname == "Clupea harengus"
  #points(Loo[ixHerring], A[ixHerring], pch=21, cex=sqrt(r[ixHerring]), col="green", bg=rgb(0,1,0) )
  
  # Highlight Gobies
  goby_list <- species_list(Family="Gobiidae")
  ixGoby <- which( growth$sciname %in% goby_list )
  points(Woo[ixGoby], A[ixGoby], pch=dots, cex=factor*sqrt(r[ixGoby]), col="green" )
  
  # Highlight Gasterosteidae (sticklbacks)
  stickle_list <- species_list(Family="Gasterosteidae")
  ixStickle <- which( growth$sciname %in% stickle_list )
  points(Woo[ixStickle], A[ixStickle], pch=dots, cex=factor*sqrt(r[ixStickle]), col="magenta" )
  
  # highlight rockfish
  rockfish_list <- species_list(Family="Sebastidae")
  ixRockfish <- which( growth$sciname %in% rockfish_list )
  points(Woo[ixRockfish], A[ixRockfish], pch=dots, cex=factor*sqrt(r[ixGoby]), col="red" )
  

  # Distribution of growth:
  out=density(log(A[ixActi]))
  defaultpanel(log(c(0.7,1e6)), ylim=log(c(0.5, 150)), new=TRUE, xaxis = FALSE, yaxis=FALSE)
  polygon(x=log(1.75e6)-1.5*(c(out$y, 0*out$y)), y=c(out$x, out$x[seq(length(out$x),1,by=-1)]),
          border=NA, col=stdgrey)
  
  # Distribution of Winf
  out=density(log(Woo[ixActi]))
  polygon(x=(c(out$x,out$x[1])), y=2*c(out$y, out$y[1])-0.92, border=NA, col=stdgrey)

  loglogpanel(xlim=c(.7, 1e6), ylim=c(0.5, 150), bExponential=TRUE, new=TRUE)
}

plotAllChapterTraits <- function() {
  #pdfplot(FUN=plotTraits, "TeX/ChapterTraits/traits.pdf", width=doublewidth, height=2*height)
  pdfplot(FUN=plotTraitsFishbase, "TeX/ChapterTraits/traits.pdf", width=doublewidth, height=2*height)
}
  
