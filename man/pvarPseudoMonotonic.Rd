\name{pvarPseudoMonotonic}

\alias{pvarPseudoMonotonic}

\title{
p-variation for a special form of \code{x}
}

\description{
Calculates p-variation of special form of function (for a part of vector \code{x}, that satisfy the requirements).
It is an inner function and should not be applied by users.
}
\usage{
pvarPseudoMonotonic(x, p, a = 1, b = length(x), sizeN = 1)
}

\arguments{
  \item{x}{a numeric vector.}
  \item{p}{a number.}
  \item{a}{a number indicating the index the first member of subset of  \code{x} that should be analysed.}
  \item{b}{a number indicating the index the last member of subset of  \code{x} that should be analysed.}
  \item{sizeN}{a number indicating the length of 'small' interval. See \code{\link{pvar}}.}
}
\details{

This function is applicable only for special form of sample. The requirements for the subset:
1. it should be pseudo-monotonic. i.e. \code{min} and \code{max} should be in the end points of the vector. 
2. removed all monotonic points.

}
\value{
The list:
  \item{ChB }{indicator of p-variation partition.}
  \item{pvarVec }{The vector of p-variations between partition points. 
	It has meaning only for those points that \code{ChB=TRUE}, elsewier the numbers are meaningless.}

}

\note{
This is inner function and shoud not be used by users. 
On the other hand this is the main function of \code{\link{pvar}}
that does the final step in p-variation calcualtion.
}


\seealso{
\code{\link{pvar}}, \code{\link{ChangePointsId}}, \code{\link{SplitByExtremum}}, \code{\link{FindAnalysingPoints}}.
}

\examples{
############################################################
xx = c(0.9552822,  0.8317665,  0.9347037,  0.3861130,  0.6546701,  0.3600638,
	0.6767309,  0.1090384,  0.4675789,  0.3223871,  0.4312971,  0.3514759,
	0.4097339,  0.2601273,  0.4978533,  0.3255348,  0.7301583,  0.5317979,
	0.6366720,  0.5096362,  0.8123783, -0.8443095)
pvM = pvarPseudoMonotonic(x = xx, p = 2, a = 1, b = 22, sizeN = 7)
sum(pvM$pvarVec[pvM$ChB]) 
pvar(xx, 2)		#gives the same result.
		
############################################################
set.seed(51)
### get x
x = rbridge(1, 1000)		
### remuves monotonic points
xx = x[ChangePointsId(x)]
### splits into pseudo-monotonic parts
s = SplitByExtremum(xx)
### for parts we can calculate p-variaiotn the partitinio
pvM = pvarPseudoMonotonic(x = xx, p = 2, a = s[4], b = s[5], sizeN = 1)
sum(pvM$pvarVec[pvM$ChB]) 
pvar(xx[s[4]:s[5]], 2)		# the same result		
 		
		
}

\keyword{ internal }

