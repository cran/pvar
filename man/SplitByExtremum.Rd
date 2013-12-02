\name{SplitByExtremum}

\alias{SplitByExtremum}
\alias{Split_MinMax}
\alias{Split_MinMax_req}


\title{Split vector according extremum (minimum and maximum).}

\description{
Splits vectors in pseudo-monotonic parts.
}
\usage{
SplitByExtremum(x)
Split_MinMax(x)
Split_MinMax_req(x, add = 1)
}

\arguments{
  \item{x}{a numeric vector.}
  \item{add}{the inner parameter for accounting in recursion.}
}

\details{
According to the properties of p-variation, 
the minimum and maximum are the points where vector might be spitted into separate subsets and analysed independent of each other. 
These functions does that.

Only the function \code{SplitByExtremum} are now used in the program. 
The functions \code{Split_MinMax_req} and \code{Split_MinMax} are older versions.
They work fine, but function \code{Split_MinMax_req} uses recursion, therefor it has some limitations.

Functions \code{Split_MinMax_req} and \code{Split_MinMax} are left for illustration purpose - 
they are much more intuitive then \code{SplitByExtremum}.
\code{Split_MinMax} splits vector only once, while \code{Split_MinMax_req} splits vector recurrently.

}
\value{
The index of the points where vector might be spitted.
}


\seealso{
\code{\link{pvar}}
}


\examples{

##### randomised sample #####
set.seed(16)
x = BridgeT(rnorm(30))

# when there are no repeated walues, the results of two functionas are identical
s1 = SplitByExtremum(x); s1
s2 = Split_MinMax_req(x); s2

# only one iteration, one min and on max.
s3 = Split_MinMax(x); s3	
c(1, which.max(x), which.min(x), length(x))

# grafical iliustration:
plot(x, type="l")
points(s1, x[s1], col=2, cex=1, pch=19)
points(s3, x[s3], col=4, cex=2)

##### special sample #####

x = c(0, -0.5, -1, 2, -1, 3, -1, 0)

### in this case functions are not identical. SplitByExtremum makes a better job.
s1 = SplitByExtremum(x); s1
s2 = Split_MinMax_req(x); s2



}



