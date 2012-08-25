\name{ChangePointsId}
\alias{ChangePointsId}

\title{Change Points of a \code{numeric} vector.}

\description{
Finds changes points in the \code{numeric} vector.
}
\usage{
ChangePointsId(x)
}

\arguments{
  \item{x}{a \code{numeric} vector.}
}
\details{

The end points of the vector will be always included in the results. 
}
\value{
The vector of index of change points.
}



\examples{
# as long as 1:10 is monotonic vector, 
# the results are only the end points of the vector.
ChangePointsId(1:10)

# there is one change point in the middle that can by find by	
ChangePointsId(c(1:10,9:1))	
which.max(c(1:10,9:1))	

# To illustrate that it relay does, 
# let look at random processes with illustration:
x = rwiener(1, 30)
plot(1:31, x, type="l")
CId = ChangePointsId(x)
lines(CId, x[CId], col=3)
points(CId, x[CId], col=2, pch=19)

}

