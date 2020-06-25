## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE, echo = TRUE,
  comment = "#>", fig.align = "center"
)
require(BIOMASS)
require(knitr)

## -----------------------------------------------------------------------------
trees <- read.csv(system.file("external", "NouraguesPlot.csv",
  package = "BIOMASS", mustWork = T
))

## ----echo=FALSE---------------------------------------------------------------
kable(head(trees), digits = 3, row.names = F, caption = "Head of the table trees")

## ---- fig.cap="Plot the coordinate long lat"----------------------------------
coord <- read.csv(system.file("external", "Coord.csv",
  package = "BIOMASS", mustWork = T
))

plot(coord[, c("Long", "Lat")], asp = 1)

## ----echo=FALSE---------------------------------------------------------------
kable(head(coord), digits = 3, row.names = F, caption = "Head of the table coord")

## ---- cache=FALSE-------------------------------------------------------------
correct_plot <- correctCoordGPS(
  longlat = coord[, c("Long", "Lat")],
  coordRel = coord[, c("xRel", "yRel")],
  rangeX = c(0, 100), rangeY = c(0, 100), drawPlot = T,
  maxDist = 10, rmOutliers = T
)

str(correct_plot, max.level = 1)

## -----------------------------------------------------------------------------
coord_num <- numberCorner(
  projCoord = correct_plot$cornerCoords,
  plot = rep("NB1", 4),
  origin = c(F, F, F, T),
  clockWise = T
)


plot(coord_num[, c("X", "Y")], asp = 1)
text(coord_num[, c("X", "Y")], labels = coord_num[, "corner"], pos = 2, offset = 0.2)

## -----------------------------------------------------------------------------
subplot <- cutPlot(
  projCoord = coord_num[, c("X", "Y")],
  plot = coord_num[, c("plot")],
  corner = coord_num[, c("corner")],
  gridsize = 25, dimX = 100, dimY = 100
)

## ----echo=FALSE---------------------------------------------------------------
kable(head(subplot))

## -----------------------------------------------------------------------------
trees$subplot <- attributeTree(trees[, c("xRel", "yRel")], rep("NB1", nrow(trees)), subplot)

## -----------------------------------------------------------------------------
trees$AGB <- computeAGB(trees$D, trees$WD, H = trees$H)

AGB <- summaryByPlot(trees$AGB, trees$subplot, drawPlot = T, subplot = subplot)

print(AGB)

## -----------------------------------------------------------------------------
TreeCoord <- attributeTreeCoord(
  xy = trees[, c("xRel", "yRel")],
  plot = trees$plot,
  coordAbs = subplot,
  dim = c(100, 100)
)

## ----echo=FALSE---------------------------------------------------------------
kable(head(TreeCoord), digits = 3, row.names = F, caption = "Head of the table TreeCoord")

## -----------------------------------------------------------------------------
#TreeCoord <- as.data.frame( proj4::project(TreeCoord, proj = correct_plot$codeUTM, inverse = T) )

## ----echo=FALSE---------------------------------------------------------------
kable(head(TreeCoord), digits = 3, row.names = F, caption = "Head of the table TreeCoord")

## -----------------------------------------------------------------------------
coordAbs = data.frame(X = c(4.066923, 4.067865, 4.067842, 4.066905), Y = c(52.68883, 52.68877, 52.68793, 52.68783))

ncoordAbs = numberCorner(projCoord = coordAbs, 
                                plot = rep("NB1", 4), 
                                origin = c(T, F, F, F), 
                                clockWise = T)

TreeCoord <- attributeTreeCoord(
  xy = trees[, c("xRel", "yRel")],
  plot = trees$plot,
  coordAbs = ncoordAbs,
  dim = c(100, 100)
)

## ----echo=FALSE---------------------------------------------------------------
kable(head(TreeCoord), digits = 3, row.names = F, caption = "Head of the table TreeCoord")

