% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/p-variation.R
\name{IsEqualPvar}
\alias{IsEqualPvar}
\title{Test if two `pvar` objects are equivalent.}
\usage{
IsEqualPvar(pv1, pv2)
}
\arguments{
\item{pv1}{an object of the class \code{pvar}.}

\item{pv2}{an object of the class \code{pvar}.}
}
\description{
Two \code{pvar} objects are considered to be equal 
if they have the same \code{x}, \code{p}, \code{value} and the same value of \code{x} 
in the points of \code{partition} (the index of partitions are not necessary the same).
All other tributes like \code{dname} or \code{TimeLabel} are not important.
}
\examples{
x <- rwiener(100)
pv1 <- pvar(x, 2)
pv2 <- pvar(x[1:50], 2) + pvar(x[50:101], 2)
IsEqualPvar(pv1, pv2)
}
