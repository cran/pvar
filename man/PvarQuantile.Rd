\name{PvarQuantile}


\alias{PvarQuantile}
\alias{PvarPvalue}

\alias{getMean}
\alias{getSd}
\alias{NormalisePvar}



\title{
  Quantiles and probabilities of p-variation
}
\description{
  The distribution of p-variation of \code{BridgeT(x)} depends on \code{n=length(x)}. 
  This fact is important for getting appropriate quantiles (or p-value).
  These functions helps to deal with it. 
}
\usage{
  PvarQuantile(n, prob = c(0.9, 0.95, 0.99), DF = PvarQuantileDF)
  PvarPvalue(n, stat, DF = PvarQuantileDF)
  
  getMean(n, bMean = MeanCoef) 
  getSd(n, bSd = SdCoef)
  NormalisePvar(x, n, bMean = MeanCoef, bSd = SdCoef) 
}

\arguments{
  \item{n}{a positive integer indicating the length of data vector.}
  \item{prob}{a cumulative probabilities of p-variation distribution.}
  \item{stat}{a vector of p-variation statistics.}
  \item{x}{a numeric vector of data values.}
  \item{DF}{a data frame with the information of p-variance distribution.}
  \item{bMean}{a coefficient vector that defines a function of the mean of p-variation.}
  \item{bSd}{a coefficient vector that defines a function of the standard deviation of p-variation.}
 
}

\details{
The distribution of p-variance is form Monte-Carlo simulation based on 140 millions iterations. 
The data frame \code{\link{PvarQuantileDF}} saves the results of Monte-Carlo simulation.

Meanwhile, \code{MeanCoef} and \code{SdCoef} defines the coefficients of functional 
form (according to \code{n}) of \code{mean} and \code{sd} statistics.

A functianal form of \code{mean} and \code{sd} statistics are the same, namely
\deqn{
	f(n) = b_1 + b_2 n^b_2 .
}{
	f(n) = b_1 + b_2 * n^b_2 .
}

And the coeficients \eqn{(b_1, b_2, b_3)} are saved in vectors \code{MeanCoef} and \code{SdCoef}.
Those vectors are estimated with \code{nls} function form Monte-Carlo simulation.


}

\value{
Functions \code{PvarQuantile} and \code{PvarPvalue} returns a corresponding value quantile or the probability.

Functions \code{getMean} and \code{getSd} returns a corresponding value of \code{mean} and \code{sd} statistics.
Function \code{NormalisePvar} returns normalise values

}

\note{
  Arguments \code{n}, \code{stat} and \code{prob} might be vectors,
  but they can't be vectors simultaneously (at least one of then must be a number).
}

\seealso{
 \code{\link{PvarBreakTest}}, \code{\link{PvarQuantileDF}},  \code{\link{NormalisePvar}},  \code{\link{getMean}}, \code{\link{getSd}}
}
\examples{
### Critical value of 0.05 significant level
CV = PvarQuantile(n=1000, prob = 0.95)
CV

### The p-value of the same significant level.
1-PvarPvalue(n=1000, stat=CV)

\dontrun{
# this may take for a while:
x100 = replicate(2000, BridgeT(rnorm(100)), simplify=FALSE)
pvar100 = sapply(x100, pvar, p=4, info=FALSE)

# sample mean:
mean(pvar100)
# a mean according to Monte-Carlo simulation:
getMean(100) 

# sample 0.95 quantile
quantile(pvar100, probs=0.95, type=8)

# a mean according to Monte-Carlo simulation:
PvarQuantile(n=100, prob = 0.95)
}



}


