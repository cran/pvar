\name{Sum_p}
\alias{Sum_p}

\title{
p-variation summation function.
}
\description{
It is the sum of absolute differences in the power of p.
}

\usage{
Sum_p(x, p, lag = 1)
}

\arguments{
  \item{x}{A numeric vector.}
  \item{p}{A number indicating the power in summing function.}
  \item{lag}{A number, indicating the lag.}
}
\details{
It is the the function that should be maximised in p-variation calculations.
}
\value{
The number equal to \code{sum((abs(diff(x, lag)))^p)}
}

\seealso{
\code{\link{pvar}}
}


\examples{

set.seed(10)
x = rbridge(1, 1000)
# the p-variation
pv = pvar(x, 2); pv			
# the same value the used the the p-varation partitnio.
Sum_p(x[pv$Partition], 2)	
# any other partition must give not greater result
Sum_p(x, 2)

}


