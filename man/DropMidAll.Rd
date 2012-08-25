\name{DropMidAll}

\alias{DropMidAll}
\alias{DropMidAllB}

\title{
Drop all mid points that appier to be insignificant is 'small' intervals.
}

\description{
Checks all inner points of small intervals. If a point appears to be insignificant in at least one case it will be drooped out.
This is an internal function and it has no practical use for users.
}
\usage{
DropMidAll(x, p, dn)

DropMidAllB(x, p, dn)
}

\arguments{
  \item{x}{a numeric vector.}
  \item{p}{a positive number indicating the power p in p-variation.}
  \item{dn}{a positive add number determining the length of 'small' interval.}
}
\details{
The requarements for the \code{x} is that the monotanic points should be removed.

}
\value{
Function \code{DropMidAll} returns values from \code{x} that have been left after checking.

Function \code{DropMidAllB} returns a \code{logical} vector indicating which values from \code{x} should be left.

}

\note{
This is inner function and shoud not be used by users. This function is loaded for illustration purport.
}


\seealso{
\code{\link{pvar}}, \code{\link{SplitByExtremum}},  \code{\link{ChangePointsId}}
}

\examples{

set.seed(6)
xx = rwiener(end = 1, frequency = 1000)
x = xx[ChangePointsId(xx)]

### finds the points that remains
id = which(DropMidAllB(x, p=4, dn=5))	
### all p-variation partition points are form id
pvar(x, p=4)$Partition %in% id 

### we can find the x points directly by DropMidAll
all(x[id] == DropMidAll(x, p=4, dn=5))
}

\keyword{ internal }

