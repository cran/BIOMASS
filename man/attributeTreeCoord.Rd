% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/attributeTreeCoord.R
\name{attributeTreeCoord}
\alias{attributeTreeCoord}
\title{Attribute GPS coordinates to trees}
\usage{
attributeTreeCoord(xy, plot, dim, coordAbs)
}
\arguments{
\item{xy}{The relative coordinates of the trees within each plot}

\item{plot}{The label of the plot (same length as the number of rows of \code{xy} or length of 1)}

\item{dim}{The dimension of the plot (either one value if the plot is a square or a vector if a rectangle)}

\item{coordAbs}{The result of the function \code{\link[=cutPlot]{cutPlot()}} or \code{\link[=numberCorner]{numberCorner()}}}
}
\value{
A data frame with two columns:
- \code{Xproj}: The \code{X} coordinates in the absolute coordinate system
- \code{Yproj}: The \code{Y} coordinates in the absolute coordinate system
}
\description{
Attribute GPS coordinates to trees
}
\examples{

# Trees relative coordinates
xy <- data.frame(x = runif(200, min = 0, max = 200), y = runif(200, min = 0, max = 200))


# cut the plot in multiple part
coord <- data.frame(X = rep(c(0, 200, 0, 200), 2), Y = rep(c(0, 0, 200, 200), 2))
coord[1:4, ] <- coord[1:4, ] + 5000
coord[5:8, ] <- coord[5:8, ] + 6000
corner <- rep(c(1, 2, 4, 3), 2)
Forestplot <- rep(c("plot1", "plot2"), each = 4)

Outcut <- cutPlot(coord, Forestplot, corner, gridsize = 100, dimX = 200, dimY = 200)


# Assign a plot to 200 trees
Forestplot <- rep(c("plot1", "plot2"), 100)

# attribute trees to subplots
attributeTreeCoord(xy, Forestplot, dim =100,coordAbs = Outcut)
}
