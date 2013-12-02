\name{DataSets}

\alias{PvarQuantileDF}
\alias{MeanCoef}
\alias{SdCoef}


\title{
  Data sets of Monte-Carlo simulations results
}
\description{  
  The test \code{\link{PvarBreakTest}} uses quantiles from Monte-Carlo simulations. 
  The results of the simulations are saved in these data sets.
}
\usage{
	PvarQuantileDF
	
	MeanCoef
	
	SdCoef
}



\details{
The distribution of p-variation of \code{BridgeT(x)} are unknown. 
The distribution of p-variance is approximated form Monte-Carlo simulation based on 141 millions iterations. 
The data frame \code{PvarQuantile} summarise the distribution of normalised statistics.
Meanwhile, \code{MeanCoef} and \code{SdCoef} defines the coefficients of functional form of \code{mean} and \code{sd} statistics.

}


\seealso{
 \code{\link{PvarBreakTest}},  
 \code{\link{PvarPvalue}},  
 \code{\link{PvarQuantile}},   
 \code{\link{getMean}}, 
 \code{\link{getSd}},
}



