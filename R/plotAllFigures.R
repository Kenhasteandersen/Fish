#
# Plot all figures
#
rm(list=ls())

# 
# Note: update the flag "bRecalcExpensiveFunctions <- FALSE" in basefunction.R
#
source("R/PlotChapterSizeSpectrumTheory.R")
plotAllChapterSizeSpectrumTheory()

source("R/PlotChapterIndividual.R")
plotAllChapterIndividual()

source("R/PlotChapterDemography.R")
plotAllChapterDemography()

source("R/PlotChapterFishing.R")
plotAllChapterFishing()

source("R/PlotChapterFIE.R")
plotAllChapterFIE()

source("R/PlotChapterDynamics.R")
plotAllChapterDynamics()

source("R/plotChapterSharks_vs_teleosts.R")
plotAllChapterSharks_vs_teleosts()

source("R/PlotChapterTraits.R")
plotAllChapterTraits()

source("R/plotChapterCommunity.R")
plotAllChapterCommunity()

source("R/plotChapterCommunityFishing.R")
plotAllChapterCommunityFishing()

source("R/plotChapterEpilogue.R")
plotAllChapterEpilogue()

