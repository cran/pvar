\name{FindAnalysingPoints}

\alias{FindAnalysingPoints}
\alias{FindAnalysingPointsOP}

\title{
Finds the points that should by analysed.
}

\description{
Finds the points that should by analysed as potentially meaningful partition points.
This is an internal function and it has no practical use for users.
}
\usage{
FindAnalysingPoints(x)
FindAnalysingPointsOP(x)
}

\arguments{
  \item{x}{a numeric vector.}
}
\details{
This function is applicable only for special form of sample. The requirements for the subset:
1. it should be quasi-monotonic. i.e. min and max should be in the end points. 
2. removed all monotonic points.

}
\value{
The index of points from \code{x} that should by analysed as potentially meaningful partition points.
}

\note{
This is inner function and should not be used by users. This function is loaded for illustration purport.
}


\seealso{
\code{\link{pvar}}, \code{\link{ChangePointsId}}, \code{\link{SplitByExtremum}}, \code{\link{pvarMon}}.
}

\examples{

xx = c(0.9552822,  0.8317665,  0.9347037,  0.3861130,  0.6546701,  0.4000638,
	0.6767309,  0.1090384,  0.4675789,  0.3223871,  0.4312971,  0.3514759,
	0.4097339,  0.2601273,  0.4978533,  0.3255348,  0.7301583,  0.5317979,
	0.6366720,  0.5096362,  0.8123783, -0.4,-0.3, -0.8443095)

### the analysing points are below all points.	
plot(xx, type="l")
AP = FindAnalysingPoints(xx)	
points(seq_along(xx)[AP], xx[AP], pch=19, col=2)
lines(seq_along(xx)[AP], xx[AP], col=2)	

### the analysing points are above all points.
plot(xx, type="l", lwd=2)
APOP = FindAnalysingPointsOP(xx)		
points(seq_along(xx)[APOP], xx[APOP], pch=19, col=3)	
lines(seq_along(xx)[APOP], xx[APOP], col=3)		

}

\keyword{ internal }

