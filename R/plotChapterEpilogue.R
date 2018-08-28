dir.create("TeX/ChapterEpilogue")

plotOptimalForaging = function() {
  p = baseparameters()
  defaultplothorizontal(nPanels = 2)
  # Expand right side of clipping rect to make room for the legend
  #par(mar=par()$mar+c(0,0,0,9),
  #    oma=par()$oma+c(0,0,0,9))
  
  #par(xaxp=c(90,190,5))
  defaultpanel(xlim=c(0,2), ylim=c(0,1), 
               xlab="Encountered food relative to $\\textit{C}_{max}$", ylab="Time in arena $\\tau^*")##), xaxis = TRUE, yaxis=TRUE)
  x = seq(0,4, length.out = 1000)
  tau = pmin(1, 1/x*sqrt(p$fc)/(1-sqrt(p$fc)))
  lines(x, tau, lwd=2)
  hline(1, lwd=2, lty=dashed)
  
  lines(x, x/(x+1), lwd=2, col=stdgrey, lty=dashed)
  lines(x, tau*x/(tau*x+1), lwd = 2, col=stdgrey)
  hline(y=p$fc)
  
  defaultpanel(xaxis = FALSE, yaxis = FALSE, 
               xlim=c(0,1), ylim=c(0,1), bty="n")
  legend("topleft", 
         bty="n",
         legend = c("Without behaviour", "With behaviour", "Time in arena", 
                    "Feeding level","Critical feeding level"),
         lty=c(dashed, 1, 1, 1, dotted),
         lwd=c(2,2,2,2,1),
         col=c(black, black, black, stdgrey, black))
}

plotAllChapterEpilogue = function() {
  pdfplot("TeX/ChapterEpilogue/OptimalForaging.pdf",
          plotOptimalForaging, width=doublewidth,
          height=height)
}