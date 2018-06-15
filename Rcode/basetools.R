#
# Background tools, mainly for plotting
#
require(ggplot2)
require(scales)
require(caTools)
require(latex2exp)
require(cowplot)
require(grImport)

singlewidth <- 8/2.54 # panel width in cm
doublewidth <- 13/2.54
height <- 8/1.6/2.54 + 2.5*0.138889

# Base R plotting tools =====================================================

bottom <- 1
left <- 2
top <- 3
right <- 4
ticklength <- 0.2 # tick mark length
omargin <- 0.7 # outer margins
cex <- 0.9
cexaxis <- 0.8 # Scaling of size of numbers on axes
iPlot <- 1 # Static variable used for labels on plots
dashed <- 2
dotted <- 3
dots <- 16
triangles <- 17
squares <- 15
stdgrey <- grey(0.4)
darkgrey <- grey(0.15)
lightgrey <- grey(0.7)
axis.lwd <- 0.35

defaultplot <- function(
  mfcol=c(1,1), 
  oma=c(0, 0.4, omargin, omargin), # Outer margins (bottom, left, top, right)
  mar=c(2.1,2.3,0,0), # Margins
  ps=10, # Point size
  tcl=ticklength, # Tick mark length
  mgp=c(1.1,0,0), ...)  # Margins for axis title, axis labels, and axis line
{
  par(ps=ps)
  par(mfcol=mfcol, oma=oma, mar=mar, tcl=tcl, mgp=mgp, cex=cex, cex.axis=cexaxis, ...) 
  assign("iPlot", 1, envir = .GlobalEnv)
}

defaultplothorizontal <- function(
  nPanels=1,
  mfcol=c(1,nPanels), 
  oma=c(2.1, 2.7, omargin, omargin), # Outer margins (bottom, left, top, right)
  mar=c(0,0,0,0.3), # Margins
#  oma=c(0,0.4, omargin, omargin), # Outer margins (bottom, left, top, right)
#  mar=c(2.7,2.5,0,0.3), # Margins
  ps=10, # Point size
  tcl=ticklength, # Tick mark length
  mgp=c(1.1,0,0), ...)  # Margins for axis title, axis labels, and axis line
{
  par(mfcol=mfcol, oma=oma, mar=mar, ps=ps, tcl=tcl, mgp=mgp, cex=cex, cex.axis=cexaxis, ...) 
  assign("iPlot", 1, envir = .GlobalEnv)
}

defaultplotvertical <- function(
  nPanels=1,
  mfcol=c(nPanels,1), 
  oma=c(2.1, 2.5, omargin, omargin), # Outer margins (bottom, left, top, right)
  mar=c(0.3,0,0,0), # Margins
  ps=10, # Point size
  tcl=ticklength, # Tick mark length
  mgp=c(2,0,0), ...)  # Margins for axis title, axis labels, and axis line
{
  par(mfcol=mfcol, oma=oma, mar=mar, ps=ps, tcl=tcl, mgp=mgp, cex=cex, cex.axis=cexaxis, ...) 
  assign("iPlot", 1, envir = .GlobalEnv)
}

test <- function() {
  defaultplot()
  defaultpanel(xlab="x", ylab="y",xlim=c(0.1,10), ylim=c(0.1,10))
  points(1,1)
  
  defaultplot()
  loglogpanel(xlab="x", ylab="y",xlim=c(0.1,10), ylim=c(0.1,10))
  points(1,1)

  defaultplothorizontal(nPanels = 3)
  loglogpanel(xlab="x", ylab="y",xlim=c(0.1,10), ylim=c(0.1,10))
  points(1,1)
  loglogpanel(xlab="x", ylab="",xlim=c(0.1,10), ylim=c(0.1,10), yaxis = FALSE)
  points(1,1)
  loglogpanel(xlab="x", ylab="",xlim=c(0.1,10), ylim=c(0.1,10), yaxis = FALSE)
  points(1,1)
  
  defaultplothorizontal(nPanels = 2)
  defaultpanel(xlab="x", ylab="y",xlim=c(0.1,10), ylim=c(0.1,10))
  points(1,1)
  defaultpanel(xlab="x", ylab="",xlim=c(0.1,10), ylim=c(0.1,10), yaxis = FALSE)
  points(1,1)
  
  defaultplotvertical(nPanels = 2)
  loglogpanel(xlab="", ylab="y",xlim=c(0.1,10), ylim=c(0.1,10), xaxis=FALSE)
  points(1,1)
  loglogpanel(xlab="x", ylab="y",xlim=c(0.1,10), ylim=c(0.1,10))
  points(1,1)
  
}

defaultpanel <- function(xlim, ylim, 
                         xlab='', ylab='', 
                         xaxis=TRUE, yaxis=TRUE, label=FALSE, new=FALSE) {
  plot(1, type='n', 
       ylim=range(ylim[!is.na(ylim)]), 
       xlim=range(xlim[!is.na(xlim)]), axes=FALSE, xlab='', ylab='', par(new=new))
  mtext(side=bottom, line=1, TeX(xlab), cex=par()$cex)
  mtext(side=left, line=1, TeX(ylab), cex=par()$cex)
  if (label) 
    makepanellabel()
  if (xaxis)
    axis(bottom, labels=xaxis, cex.axis=par()$cex.axis, lwd=axis.lwd, lwd.ticks=axis.lwd)
  if (yaxis)
    axis(left, labels=yaxis, cex.axis=par()$cex.axis, lwd=axis.lwd, lwd.ticks=axis.lwd)
  box()

}

semilogxpanel <- function(xlim, ylim, xlab='', ylab='', 
                          xaxis=TRUE, yaxis=TRUE, label=FALSE, new=FALSE) {
  ylim <- range(na.omit(ylim))
  xlim <- range(na.omit(xlim))
  plot(1, type='n', log='x',
       ylim=ylim, 
       xlim=xlim, axes=FALSE, xlab='',ylab='', par(new=new))

  mtext(side=bottom, line=1, TeX(xlab))
  mtext(side=left, line=1, TeX(ylab))
  if (label) 
    makepanellabel()
#  if (xaxis)
  logaxes(bottom, lim=xlim, labels=xaxis)
#  if (yaxis)
  axis(left, labels=yaxis, lwd=axis.lwd, lwd.ticks=axis.lwd)
  box(lwd=axis.lwd)
}

semilogypanel <- function(xlim, ylim, xlab='', ylab='', 
                          xaxis=TRUE, yaxis=TRUE, label=FALSE) {
  plot(1, type='n', log='y',
       ylim=range(ylim[!is.na(ylim)]), 
       xlim=range(xlim[!is.na(xlim)]), axes=FALSE, xlab='',ylab='')
  mtext(side=bottom, line=1, TeX(xlab))
  mtext(side=left, line=1.5, TeX(ylab))
  if (label) 
    makepanellabel()
  axis(bottom, labels=xaxis, lwd=axis.lwd, lwd.ticks=axis.lwd)
  logaxes(left, lim=ylim, labels=yaxis)
  box(lwd=axis.lwd)
}

loglogpanel <- function(xlim, ylim, xlab='', ylab='', 
                        xaxis=TRUE, yaxis=TRUE, label=FALSE,
                        bExponential=TRUE, new=FALSE, powx=NA, powy=NA) {
  ylim <- range(na.omit(ylim))
  xlim <- range(na.omit(xlim))
  plot(1, type='n', log='xy',
       ylim=ylim, 
       xlim=xlim, axes=FALSE, xlab='',ylab='', par(new=new))
  mtext(side=bottom, line=1.1, TeX(xlab))
  mtext(side=left, line=1.5, TeX(ylab))
  if (label) 
    makepanellabel()
  logaxes(bottom, lim=xlim, bExponential = bExponential, labels=xaxis, pow=powx)
  logaxes(left, lim=ylim, bExponential = bExponential, labels=yaxis, pow=powy)
  box(lwd=axis.lwd)
}


## side:  1 for x-axis, 2 for y-axis (3, 4 secondary x and y axes)
## labels: logical, if TRUE, adds them
## las: integer in [0,1,2,3] see ?par
logaxes <- function(side = bottom, 
                    lim, pow = NA,
                    labels=TRUE, las = 1, col = 1, bExponential=TRUE) {
  poww <- pow
  if (is.na(pow[1]))
    poww = ceiling(min(log10(lim))):floor(max(log10(lim)))
  
  axis(side, at = 10^poww, lwd=0, labels = FALSE, tcl=ticklength, lwd.ticks = axis.lwd, col = col)

  if (labels) 
    axis(side, at = 10^poww, lwd=0, labels = FALSE, tcl=-0., lwd.ticks = axis.lwd)

  if (side==bottom)
    line <- 0.2/(par()$oma[bottom] + par()$mar[bottom])
  else
    line <- 0.6/(par()$oma[left] + par()$mar[left])
  
  if (labels) {
    for(i in poww) {
      if (bExponential)
        mtext(side = side, at = 10^i, text=bquote(10^.(i)), line = line, 
              las = las, cex=cexaxis)
      else
        mtext(side = side, at = 10^i, text=10^i, line = 0.2, 
              las = las, cex=par()$cex.axis)
    }
  }
  # Minor ticks:
  if (is.na(pow[1])) {
    at <- as.vector(sapply(seq(range(poww)[1]-1, range(poww)[2]), function(p) (2:9)*10^p))
    axis(side, at = at, labels = FALSE, tcl=0.5*ticklength, lwd=0, lwd.ticks=axis.lwd, col = col)
  }
} 

makepanellabel <- function() {
  mtext(letters[iPlot], side=top, line=-1.1, adj=0.05)
  assign("iPlot", iPlot+1, envir = .GlobalEnv)
}

hline <- function(y=0, lty='dotted') {
  if (par("xlog"))
    lines(x=10^par("usr")[1:2], y=y*c(1,1), lty=lty)
  else
    lines(x=par("usr")[1:2], y=y*c(1,1), lty=lty)
}

vline <- function(x=0, lty='dotted', col="black") {
  if (par("ylog"))
    lines(x=x*c(1,1), y=10^par("usr")[3:4], lty=lty, col=col)
  else
    lines(x=x*c(1,1), y=par("usr")[3:4], lty=lty, col=col)
}

pdfplot <- function(filename, FUN, ..., width=singlewidth, height=height) {
  pdf(filename, width=width, height=height, useDingbats=FALSE)
  FUN(...)
  dev.off()  
}

addEpsPicture <- function(sName, x, y, width=1) {
  # Convert the picture
  sOutName = paste(sName,'.xml',sep='')
  PostScriptTrace(sName, sOutName)
  # Add it
  grid.picture(readPicture(sOutName), x=x, y=y, width=width)
}

ribbon <- function(x,ymin=NA,ymax,col=lightgrey) {
  x <- c(x, x[seq(length(x),1,by = -1)])
  polygon(x, c(ymin, ymax[seq(length(ymax),1,by = -1)]), col=col, border=NA)
}

tightaxes <- function()
  par(xaxs="i", yaxs="i")

#
# ggplot plotting tools ==============================================
#

logseq <- function(xmin, xmax, n=100)
  10^seq(log10(xmin), log10(xmax), length.out = n)

thick = 1 # thick line
thin = 0.25 # thin line

# Breaks for a log axis:
bks <- function(lims) {
  10^seq(floor(log10(lims[1])), floor(log10(lims[2])))
}

# Standard for a double-log plot:
loglog <- function(fig, ylim=NULL, xlim=NULL, label="", 
                   bXaxis=TRUE, bYaxis=TRUE,
                   ixStep=1, iyStep=1, bExpand=FALSE) {
  
    expand <- c(0.05,0)
    if (bExpand)
      expand <- c(0,0)

  fig <- fig + 
    scale_x_log10(breaks = function(lims) 10^seq(floor(log10(lims[1])), floor(log10(lims[2])), by=ixStep),
                  labels = trans_format("log10", math_format(10^.x)),
                  limits = xlim, oob=rescale_none,
                  expand=expand) +
    scale_y_log10(breaks = function(lims) 10^seq(floor(log10(lims[1])), floor(log10(lims[2])), by=iyStep),
                  labels = trans_format("log10", math_format(10^.x)),
                  oob=rescale_none, limits=ylim,
                  expand=expand) +
    annotation_logticks(long=unit(1,"mm"),mid=unit(0.5,"mm"), 
                        short=unit(0.5,"mm"), size=0.25)  
  # Add label
  #ylimit <- ggplot_build(fig)$panel$ranges[[1]]$y.range
  #xlimit <- ggplot_build(fig)$panel$ranges[[1]]$x.range
  if (label != "")
    fig <- fig + 
      annotate("text", x=xlim[1], y=ylim[2], label=label, hjust=0, vjust=1)
  
  fig <- mytheme(fig, bXaxis=bXaxis, bYaxis=bYaxis)
  fig
}

# standard for semilogx plot:
semilogx <- function(fig, xlim=NULL, ylim=NULL, label="", 
                     bXaxis=TRUE, bYaxis=TRUE) {
  fig <- fig + 
    scale_x_log10(breaks = bks, #10^seq(floor(xrange[1]), floor(xrange[2]), by=1), #trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)),
                  limits = xlim, oob=rescale_none) +
    annotation_logticks(long=unit(1,"mm"),mid=unit(0.5,"mm"), 
                        short=unit(0.5,"mm"), size=0.25, sides="b") +
    scale_y_continuous(oob=rescale_none, limits=ylim)
  # Add label
  #xlimit <- ggplot_build(fig)$layout$panel_ranges[[1]]$x.range
  if (label != "")
    fig <- fig + 
      annotate("text", x=xlim[1], y=Inf, label=label, hjust=0, vjust=1.5)
  
  fig <- mytheme(fig, bXaxis=bXaxis, bYaxis=bYaxis)
  fig
}

# standard for semilogy plot:
semilogy <- function(fig, ylim=NULL, label="") {
  fig <- fig + 
    scale_y_log10(breaks = bks,
                  labels = trans_format("log10", math_format(10^.x)),
                  oob=rescale_none, limits=ylim) +
    annotation_logticks(long=unit(1,"mm"),mid=unit(0.5,"mm"), 
                        short=unit(0.5,"mm"), size=0.25, sides="l") 
  fig <- mytheme(fig)
  fig
}

# My standard theme
mytheme <- function(fig, bXaxis=TRUE, bYaxis=TRUE) {
  fig <- fig + 
    theme_classic(base_size=9) +
    theme(plot.title=element_text(size=9, hjust=0.5)) +
    theme(panel.border = element_rect(fill=NA, size=0.35)) +
    theme(axis.line = element_blank()) +
    # Ticks facing inwards:
    theme(axis.ticks = element_line(size=thin),
          axis.ticks.length = unit(-0.1, "cm"),
          axis.text.x = element_text(margin=unit(c(0.2,0.0,0.0,0.0), "cm")), 
          axis.text.y = element_text(margin=unit(c(0.0,0.15,0.0,0.1), "cm")) )
  # Remove axis tick labels?
  if (!bXaxis)
    fig <- fig +
      labs(x=NULL) +
      theme(axis.text.x=element_blank()) # Remove ticks
  # Remove y axis tick lables?
  if (!bYaxis)
    fig <- fig +
      labs(y=NULL) +
      theme(axis.text.y=element_blank()) # Remove ticks

  fig
}

plotlabel <- function(fig, label) {
  fig <- fig + 
    annotate("text", x=-Inf, y=Inf, label=label, hjust=-.2, vjust=2)
  fig
}


# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}