% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/twoDxc.R
\name{subtractIon}
\alias{subtractIon}
\title{Subtract ion}
\usage{
subtractIon(object, ion, tol)
}
\arguments{
\item{object}{An xsAnnotate object returned by CAMERA that needs 2D peak
grouping}

\item{ion}{The ion to subtract}

\item{tol}{The tolerance window for the ion in ppm}
}
\value{
An xsAnnotate2D object with the ion intensity subtracted
}
\description{
This function subtracts an ion from the MSnbase object by removing its
intensity from the TIC.
}
