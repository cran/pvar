\name{ProcessGenerators}
\alias{rbridge}
\alias{rwiener}

\title{
Random process generators
}
\description{
Generate a trajectory of random processes.  
}
\usage{
rwiener(end = 1, frequency = 1000)

rbridge(end = 1, frequency = 1000)
}

\arguments{
  \item{end}{a number. The end point of the process in the 'time' scale. }
  \item{frequency}{a number specifying the size of trajectory vector. The trajectory will start at point 0 
	and will have \code{frequency} more observations. The length of the results will be \code{frequency+1} . }
}
\details{
\code{rwiener} generate Wiener process via partial sums process. 

\code{rbridge} generate Brownian bridge via \code{rwiener}.
}
\value{
A numerics vector of the length \code{frequency+1}.
}

\examples{

x = rwiener(1, 10000)
plot(x, type="l")

x = rbridge(1, 10000)
plot(x, type="l")

}

\keyword{ts}

