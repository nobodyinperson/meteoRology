% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meteorology.functions.R
\name{wind.direction}
\alias{wind.direction}
\title{Meteorological wind direction}
\usage{
wind.direction(wind = c(u, v), u = wind[1], v = wind[2])
}
\arguments{
\item{wind}{numerical vector of wind components c(u,v)}

\item{u}{lateral wind component}

\item{v}{longitudinal wind component}
}
\value{
The meteorological wind direction of wind vector c(u,v), that is the angle from y-axis in mathematical negative direction to the inverted vector c(-u,-v).
}
\description{
Meteorological wind direction
}

