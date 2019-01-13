#
# Plot all figures
#
devtools::install_github("https://github.com/Kenhasteandersen/FishSizeSpectrum.git")

bRecalcExpensiveFunctions = FALSE

plotAllFigures = function(bRecalcExpensiveFunctionsFlag=FALSE) {
  # 
  # Note: update the flag "bRecalcExpensiveFunctions <- FALSE" in basefunction.R
  #
  source("R/PlotChapterSizeSpectrumTheory.R")
  bRecalcExpensiveFunctions = bRecalcExpensiveFunctionsFlag
  plotAllChapterSizeSpectrumTheory()
  
  source("R/PlotChapterIndividual.R")
  bRecalcExpensiveFunctions = bRecalcExpensiveFunctionsFlag
  plotAllChapterIndividual()
  
  source("R/PlotChapterDemography.R")
  bRecalcExpensiveFunctions = bRecalcExpensiveFunctionsFlag
  plotAllChapterDemography()
  
  source("R/PlotChapterFishing.R")
  bRecalcExpensiveFunctions = bRecalcExpensiveFunctionsFlag
  plotAllChapterFishing()
  
  source("R/PlotChapterFIE.R")
  bRecalcExpensiveFunctions = bRecalcExpensiveFunctionsFlag
  plotAllChapterFIE()
  
  source("R/PlotChapterDynamics.R")
  bRecalcExpensiveFunctions = bRecalcExpensiveFunctionsFlag
  plotAllChapterDynamics()
  
  source("R/plotChapterSharks_vs_teleosts.R")
  bRecalcExpensiveFunctions = bRecalcExpensiveFunctionsFlag
  plotAllChapterSharks_vs_teleosts()
  
  source("R/PlotChapterTraits.R")
  bRecalcExpensiveFunctions = bRecalcExpensiveFunctionsFlag
  plotAllChapterTraits()
  
  source("R/plotChapterCommunity.R")
  bRecalcExpensiveFunctions = bRecalcExpensiveFunctionsFlag
  plotAllChapterCommunity()
  
  source("R/plotChapterCommunityFishing.R")
  bRecalcExpensiveFunctions = bRecalcExpensiveFunctionsFlag
  plotAllChapterCommunityFishing()
  
  source("R/plotChapterEpilogue.R")
  bRecalcExpensiveFunctions = bRecalcExpensiveFunctionsFlag
  plotAllChapterEpilogue()
}
