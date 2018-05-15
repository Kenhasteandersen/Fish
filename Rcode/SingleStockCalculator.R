#
# Shiny web application to calculate demographic quantities for a fish stock.
#
# Run the application with:
#  shinyApp(ui = SingleStockCalculatorUI, server = SingleStockCalculatorServer)
#
setwd("~/Documents/Projects/Fish")
library(shiny)
source("Rcode/basetools.R")
source("Rcode/basefunctions.R")
source("Rcode/baseparameters.R")
source("Rcode/QuantitativeGenetics.R")

# Define the user interface
SingleStockCalculatorUI <- fluidPage(
  # Make rules widers:
  tags$head(
    tags$style(HTML("hr {border-top: 1px solid #444444;}"))
  ),
  #
  # Main text
  #
  h1('Single stock size spectrum calculator'),
  p('Calculate all demographic quantities and the response to fishing of single stock',
    'using the procedures described in ',em('Fish: Ecology, Evolution, and Exploitation.')),
  
  sidebarLayout(
    #
    # The side bar where parameters are defined:
    #
    sidebarPanel(
      h3('Fishing'),
      
      checkboxInput(inputId="bMSY", 
                    label="Fish at MSY", 
                    value = FALSE, width = NULL),
      
      conditionalPanel(
        condition = "input.bMSY==false",
        sliderInput(inputId = "F",
                    "Fishing mortality",
                    min = 0,
                    max = 3,
                    step = 0.05,
                    value = 0.3)
      ),
      
      sliderInput(inputId = "etaF10",
                  "Log10(Minimum size for fishing rel. to W)",
                  min = -3,
                  max = 0,
                  step = 0.05,
                  value = log10(baseparameters()$etaF)),
      
      hr(),
      h3('Parameters for the stock'),
      
      # Select type of parameters:
      checkboxInput("bClassicParameters", 
                    label = "Show classic parameters", 
                    value = FALSE,
                    width = NULL),
      
      conditionalPanel(
        condition = "input.bClassicParameters==false",
        #
        # Physiological parameters:
        #
        sliderInput(inputId = "W10",
                    "Log10(Asymptotic weight)",
                    min = 0,
                    max = 6,
                    step = 0.1,
                    value = 4),
        
        sliderInput(inputId = "A10",
                    "Log10(Growth rate coefficient)",
                    min = 0,
                    max = 2,
                    step = 0.1,
                    value = log10(baseparameters()$A)),
        
        sliderInput(inputId = "a",
                    "Physiological mortality",
                    min = 0.1,
                    max = 1,
                    step = 0.05,
                    value = baseparameters()$a),
        
        sliderInput(inputId = "epsR10",
                    "Log10(Reproductive efficiency)",
                    min = -4,
                    max = 0,
                    step = 0.05,
                    value = log10(baseparameters()$epsR)),
        
        sliderInput(inputId = "etaM",
                    "Maturation rel. to asymp. weight",
                    min=0.01,
                    max = 1,
                    step = 0.01,
                    value = baseparameters()$etaM)
      ),
      
      conditionalPanel(
        condition = "input.bClassicParameters==true",
        #
        # Classic parameters:
        #
        #hr(),
        #h4('Classic parameters'),
        
        sliderInput(inputId = "L",
                    "Asymptotic length",
                    min = 1,
                    max = 200,
                    step = 1,
                    value = -1)
        ,
        sliderInput(inputId = "K10",
                    "log10(von B. growth coef.)",
                    min = -2,
                    max = log10(20),
                    step = .1,
                    value = 1.2)
        ,
        sliderInput(inputId = "M",
                    "Adult natural mortality",
                    min = 0,
                    max = 2,
                    step = .1,
                    value = .2)
        ,
        sliderInput(inputId = "alpha10",
                    "log10(Recruitment coefficient)",
                    min = -1,
                    max = 4,
                    step = 1,
                    value = log10(baseparameters()$epsR*baseparameters()$epsEgg/baseparameters()$w0 * baseparameters()$A * (10000)^(1-baseparameters()$n)))
        ,
        sliderInput(inputId = "tmat",
                    "Age at maturation",
                    min = 0,
                    max = 10,
                    step = .1,
                    value = ageMaturation(10000))
      )
      
    ),
    # Tabs with results:
    mainPanel(
      tabsetPanel(
        tabPanel('Parameters',
                 uiOutput("physparameters")),
        tabPanel('Size spectrum', 
                 plotOutput("plotSpectrum"), 
                 plotOutput("plotSizeAtAge")),
        tabPanel('Overview', 
                 uiOutput("results"),
                 plotOutput("plotState", width="30%")),
        tabPanel('Fisheries reference points', 
                 plotOutput("plotRefpoints", width="30%"),
                 plotOutput("plotBiomassRefpoints", width="30%")),
        tabPanel('Fisheries induced evolution',
                 plotOutput("plotFIE", width="30%")),
        selected='Size spectrum'
      )
    )
  )
)

updatePhysParam <- function(input, session) {
  W = weight(input$L)
  p <- baseparameters()
  K = 10^input$K10
  tmat = input$tmat
  L = input$L
  alpha = 10^input$alpha10
  M = input$M
  
  etaM = (27*K^3*tmat^3)/64
  A = sqrt(2)*3^(3/4)*((p$c*K^3*L^3)/tmat)^(1/4)
  a = (M*tmat)/4
  epsR = ((tmat/K^3)^(1/4)*p$w0*alpha)/(sqrt(2)*3^(3/4)*p$epsEgg)
  
  updateSliderInput(session=session,
                    inputId = "W10",
                    value = log10(W))
  updateSliderInput(session=session,
                    inputId = "A10",
                    value = log10(A))
  updateSliderInput(session=session,
                    inputId = "a",
                    value = a)
  updateSliderInput(session=session,
                    inputId = "etaM",
                    value = etaM)
  updateSliderInput(session=session,
                    inputId = "epsR10",
                    value = log10(epsR))
  
}

updateClassicParam <- function(input, session) {
  L = weight2length(10^input$W10)
  p = baseparameters()
  K = 10^input$A10 * L^(-3/4) / (3*p$c^(1/4) * input$etaM^(-1/12) )
  M = 3*input$a * input$etaM^(-1/3)*K
  n = baseparameters()$n
  tmat = input$etaM^(1-n)/(10^input$A10*(1-n)) * (10^input$W10)^(1-n)
  alpha = 10^input$epsR10*p$epsEgg/p$w0 * 10^input$A10 * (10^input$W10)^(n-1)
  
  updateSliderInput(session=session, 
                    inputId = "L", 
                    value = L)
  updateSliderInput(session=session, 
                    inputId = "K10", 
                    value = log10(K))
  updateSliderInput(session=session, 
                    inputId = "M", 
                    value = M)
  updateSliderInput(session=session, 
                    inputId = "tmat", 
                    value = tmat)
  updateSliderInput(session=session, 
                    inputId = "alpha10", 
                    value = log10(alpha))
}


# Define server logic
SingleStockCalculatorServer <- function(input, output, session) {
  #
  # Updated connected sliders between physilogical and classic parameters
  #
  observeEvent({ 
    input$bClassicParameters
  },
  {
    if ((input$bClassicParameters==TRUE) | (input$L==1)) {
      updateClassicParam(input, session)
    }
  })
  observeEvent({
    input$bClassicParameters
  }, {
    if ((input$bClassicParameters==FALSE) & (input$L!=1)) {
      updatePhysParam(input, session)
    }
  } )
  #
  # Simulate the stock whenever a parameter is updated
  #
  simResults <- eventReactive({
    input$W10
    input$A10
    input$a
    input$etaF10
    input$epsR10
    input$etaM
    input$bMSY
    input$F
    input$L
    input$K10
    input$M
    input$bClassicParameters
    
  }, {
    param <- baseparameters()
    
    #
    # Extract physiological parameters:
    #
    if (input$bClassicParameters == TRUE) 
      updatePhysParam(input, session)
    
    W = 10^input$W10
    param$W = W
    param$A <- 10^input$A10
    param$a <- input$a
    param$epsR <- 10^input$epsR10
    param$etaM <- input$etaM
    #
    # Extract parameters related to fishing
    #
    param$etaF <- 10^input$etaF10
    #
    # Calculate stock without fishing:
    #
    param$F=0
    spec0 = spectrum(param)
    #
    # Calculate reference points:
    #
    refs <- calcRefpoints(param)
    #
    # Calculate stock with fishing
    #
    param$F <- refs$Fmsy
    if (!input$bMSY)
      param$F = input$F
    spec <- spectrum(param)
    #
    # Setup physioligical parameters data frame:
    #
    physparameters <- data.frame(
      Parameter = c("<i>Physiological parameters</i>",
                    "&nbsp;&nbsp;Asymptotic weight, <i>W</i><sub>&#x221e;</sub>",
                    "&nbsp;&nbsp;Growth coefficient, <i>A</i>",
                    "&nbsp;&nbsp;Physiological mortality, <i>a</i>",
                    "&nbsp;&nbsp;Reproductive efficiency, &epsilon;<sub><i>r</i></sub>",
                    "&nbsp;&nbsp;Maturation relative to asymptotic weight, &eta;<sub><i>m</i></sub>",
                    "",
                    "<i>Classic parameters</i>",
                    "&nbsp;&nbsp;Asymptotic length, <i>L</i><sub>&#x221e;</sub>",
                    "&nbsp;&nbsp;Von Bertalanffy growth coefficient, <i>K</i>",
                    "&nbsp;&nbsp;Adult mortality, <i>M</i>",
                    "&nbsp;&nbsp;Stock-recruitment parameter, <i>&alpha;</i>",
                    "&nbsp;&nbsp;Age at maturation, <i>t</i><sub>mat</sub>"),
      Value = c("",
                sprintf("%6.0f g", param$W),
                sprintf("%2.1f g<sup>1/4</sup>year<sup>-1</sup>", param$A),
                sprintf("%1.2f", param$a),
                sprintf("%1.3f", param$epsR),
                sprintf("%1.2f", param$etaM),
                "",
                "",
                sprintf("%3.1f cm", input$L),
                sprintf("%2.2f yr<sup>-1</sup>", 10^input$K10),
                sprintf("%1.2f yr<sup>-1</sup>", input$M),
                sprintf("%3.0f g<sup>-1</sup>yr<sup>-1</sup>", 10^input$alpha10),
                sprintf("%2.1f yr", input$tmat)
      )
    )
    #
    # Setup classic parameters data frame:
    #
    
    #
    # Calculate selection respones:
    #
    response <- calcSelectionResponse(p=baseparamQG(wm=param$etaM*param$W,p=param), F=param$F, W=param$W)
    #
    # Set the results data frame:
    #
    results <- data.frame(
      Quantity = c("<i>Recruitment</i>",
                   "&nbsp;&nbsp;Without fishing",
                   "&nbsp;&nbsp;When fishing",
                   "<i>Spawning stock biomass</i>",
                   "&nbsp;&nbsp;When fishing",
                   "<i>Reference points</i>",
                   "&nbsp;&nbsp;Maximum sustainable yield",
                   "&nbsp;&nbsp;Limit",
                   "&nbsp;&nbsp;Crash",
                   "<i>Selection responses",
                   "&nbsp;&nbsp;Size at maturation",
                   "&nbsp;&nbsp;Growth rate",
                   "&nbsp;&nbsp;Reproductive investment"
      ),
      Symbol = c("", 
                 "<i>R</i>/<i>R</i><sub>max</sub>",
                 "<i>R</i>/<i>R</i><sub>max</sub>",
                 "", 
                 "<i>B<sub>SSB</sub>/<i>B</i><sub>max</sub>",
                 "",
                 "<i>F</i><sub>msy</sub>",
                 "<i>F</i><sub>lim</sub>",
                 "<i>F</i><sub>crash</sub>",
                 "",
                 "<i>G</i><sub>rl.wm</sub>",
                 "<i>G</i><sub>rl.A</sub>",
                 "<i>G</i><sub>rl.R</sub>"
      ),
      Value=c("",
              sprintf("%0.2f",spec0$R[1]),
              sprintf("%0.2f",spec$R[1]),
              "",
              sprintf("%0.2f",spec$SSB[1]/spec0$SSB[1]),
              "",
              sprintf("%0.2f year<sup>-1</sup>",refs$Fmsy),
              sprintf("%0.2f year<sup>-1</sup>",refs$Flim),
              sprintf("%0.2f year<sup>-1</sup>",refs$Fcrash),
              "",
              sprintf("%0.3f &percnt;/year",100*response$dwmdt),
              sprintf("%0.3f &percnt;/year",100*response$dAdt),
              sprintf("%0.3f &percnt;/year",100*response$dkrdt)))
    
    return(list(param=param, spec0=spec0, spec=spec, refs=refs, response=response, physparameters=physparameters, table=results))#, parameters=parameters))
  })
  #
  # Table with physiological parameters
  #
  output$physparameters <- renderTable({
    simResults()$physparameters
  }, sanitize.text.function = function(x) x)
  #
  # Table with classic parameters:
  #
  output$classicparameters <- renderTable({
    simResults()$classicparameters
  }, sanitize.text.function = function(x) x)
  
  #
  # Table with results:
  #
  output$results <- renderTable({
    simResults()$table
  }, sanitize.text.function = function(x) x)
  #
  # Plot the size spectrum
  # 
  output$plotSpectrum <- renderPlot({
    spec0 = simResults()$spec0
    spec = simResults()$spec
    param = simResults()$param
    W = param$W
    
    defaultplot(mar=c(5.1,2.3,0,0), cex=2, cex.axis=1.5, ps=16)
    B0 = spec0$N * spec0$w^2
    B0 = B0/B0[1]
    yrange = c(1e-4,1)*max(B0)
    loglogpanel(xlim = c(1e-3, W), ylim=yrange,
                xlab='', ylab=('Biomass (g)'))
    mtext("Weight (g)     ", side=bottom, line=0, adj=1, at=1e-3)
    # Axis with ages:
    wa = weightatage(param, ages=seq(0,60))
    axis(side=bottom, line=2,
         at=wa$w,
         labels=FALSE, col="blue")
    mtext("Age      ", side=bottom, line=1.5, at=1e-3, adj=1, col="blue")
    # Remove superfluous age-labels:
    wa$w10 = log10(wa$w)
    mindist = (max(wa$w10)-min(wa$w10))/30
    ix = 2
    while (ix<dim(wa)[1]) {
      if ( (wa$w10[ix]-wa$w10[ix-1]) < mindist) {
        wa = rbind( wa[1:ix-1,], wa[(ix+1):dim(wa)[1],] )
      } else
        ix = ix + 1
    }
    wa = wa[1:dim(wa)[1]-1,]
    mtext(side=bottom, 
          at = wa$w,
          line=2,
          text=wa$age, col="blue")
    lines(spec0$w, B0, lty='dashed', lwd=2)
    # Size spectrum
    B = spec$N * spec$w^2
    B = B/B[1]
    lines(spec$w, B, lwd=3)
    # Vertical line at maturation
    vline(param$etaM*param$W)
    #
    # Add fishing mortality:
    #
    F = matrix( param$funcFishing(spec$w,param) )
    ymin = 10^par("usr")[3]
    w <- spec$w
    w[length(F)] = 10^par("usr")[2]
    ribbon(w, ymin=0*F+ymin, ymax=10^(F)*ymin)
    #print(10^(F)* min(yrange))
  })
  #
  # Plot size-at-age:
  #
  output$plotSizeAtAge <- renderPlot({
    param = simResults()$param
    W = param$W
    
    defaultplot(cex=2, cex.axis=.8)#, ps=16)
    ages = seq(0,5*ageMaturation(W, param),length.out=500)
    out <- ode(y=param$w0, times=ages, func=function(ages,w,p) list(growth(param,w)), parms=param)
    plot(x=out[,1],y=out[,2], lwd=3, type = "l",
         xlim=c(0.01, 5*ageMaturation(W, param)), ylim=c(0.001,W),
         xlab="Age", ylab="Weight (g)")
    vline(ageMaturation(W,param))
    hline(param$etaM*W)
  })
  #
  # Stock state.
  #
  output$plotState <- renderPlot({
    par(mar=c(5,12,4,2))
    barplot(c(simResults()$spec0$R[1], 
              simResults()$spec$R[1], 
              simResults()$spec$SSB[1]/simResults()$spec0$SSB[1])
            ,
            names.arg = c(TeX("$\\frac{\\textit{R}}{\\textit{R}_{max}}$   without fishing"),
                          TeX("$\\frac{\\textit{R}}{\\textit{R}_{max}}$   "),
                          TeX("$\\frac{\\textit{B}}{\\textit{B}_{max}}$   "))
            ,
            main="Stock state",
            xlim=c(0,1),
            horiz=TRUE, las=1,
            border=NA)
  }, width=500)
  #
  # Reference points:
  #
  output$plotRefpoints <- renderPlot({
    par(mar=c(5,12,4,2))
    barplot(c(simResults()$refs$Fcrash,
              simResults()$refs$Flim, 
              simResults()$refs$Fmax,
              simResults()$refs$Fmsy),
            names.arg = c(TeX("$\\textit{F}_{crash}$"),
                          TeX("$\\textit{F}_{lim}$"),
                          TeX("$\\textit{F}_{max}$"),
                          TeX("$\\textit{F}_{msy}$")),
            main="Moratlity reference points",
            xlab=TeX("Fishing mortality (yr$^{-1}$)"),
            horiz=TRUE, las=1,
            border=NA)
    if (!input$bMSY)
      vline(x=input$F, col="blue")
  }, width=500)
  
  output$plotBiomassRefpoints <- renderPlot({
    par(mar=c(5,12,4,2))
    barplot(c(simResults()$refs$Bmsy/simResults()$spec0$SSB[1],
              simResults()$refs$Blim/simResults()$spec0$SSB[1]), 
            names.arg = c(
              TeX("$\\frac{\\textit{B}_{msy}}{\\textit{B}_{max}}$"),
              TeX("$\\frac{\\textit{B}_{lim}}{\\textit{B}_{max}}$")),
            main="Biomass reference points",
            xlab=TeX(""),
            xlim=c(0,1),
            horiz=TRUE, las=1,
            border=NA)
    if (!input$bMSY)
      vline(x=simResults()$spec$SSB[1]/simResults()$spec0$SSB[1], col="blue")
  }, width=500)
  #
  # Fisheries induced evolution (selection responses):
  # 
  output$plotFIE <- renderPlot({
    par(mar=c(5,12,4,2))
    barplot(c(100*simResults()$response$dkrdt,
              100*simResults()$response$dAdt, 
              100*simResults()$response$dwmdt),
            names.arg = c("Investment in reproduction",
                          "Growth rate",
                          "Size at maturation"),
            xlab="Relative selection response (%/yr)",
            horiz=TRUE, las=1,
            border=NA)
    vline(0)
  }, width=500)
  
  
}

# Run the application 
shinyApp(ui = SingleStockCalculatorUI, server = SingleStockCalculatorServer)

