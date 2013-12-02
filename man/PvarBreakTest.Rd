\name{PvarBreakTest}

\alias{PvarBreakTest}
\alias{plot.PvarBreakTest}
\alias{summary.PvarBreakTest}
\alias{print.summary.PvarBreakTest}

\alias{BreakPoints}

\title{
	Structural break test
}
\description{
	This function performs structural break test based on p-variation. 
}
\usage{
PvarBreakTest(x, TimeLabel = as.vector(time(x)), alpha = 0.05, FullInfo = TRUE) 

\method{summary}{PvarBreakTest}(object, ...)

\method{plot}{PvarBreakTest}(x, main1 = "Data", main2 = "Bridge transformation", 
    ylab1 = x$dname, ylab2 = "BridgeT(" \%.\% x$dname \%.\% ")", sub2=NULL, 
	col.PP = 3, cex.PP = 0.5, col.BP = 2, cex.BP = 1, cex.DP = 0.5, ...)
	
BreakPoints(object) 	
}

\arguments{
  \item{x}{a numeric vector of data values or an object of class \code{pvar}.}
  \item{TimeLabel}{numeric, a time index of \code{x}. Used only for plotting.}
  \item{alpha}{a small number greater then 0. It indicates the significant level of the test.}
  \item{FullInfo}{\code{logical}. If \code{TRUE} (the default) the function will return an object of the class \code{PvarBreakTest} 
		that saves all useful information.	Otherwise only the statistics will by returned. }
  \item{object}{the object of the class \code{PvarBreakTest}.}
  \item{main1}{the \code{main} parameter of the data graph.}
  \item{main2}{the \code{main} parameter of the Bridge transformation graph.}
  \item{ylab1}{the \code{ylab} parameter of the data graph.}
  \item{ylab2}{the \code{ylab} parameter of the Bridge transformation graph.}
  \item{sub2}{the \code{sub} parameter of the Bridge transformation graph. By default it reports the number of break points.}
  \item{cex.DP}{the cex of data points.}
  \item{col.PP}{the color of partition points.}
  \item{cex.PP}{the cex of partition points.}
  \item{col.BP}{the color of break points.}
  \item{cex.BP}{the cex of break points.}
  \item{\dots}{further arguments.}	
}

\details{
Lets \code{x} be a data that should be tested of structial breaks. 
Then the p-variation of the \code{BridgeT(x)} with \code{p=4} is the test's statistics.

The quantiles of H0 distribution is based on Monte-Carlo simulation of 140 millions iterations. 
The test is reliable then \code{length(x)} is between 100 and 10000.
The test might work with other lengths too, but it is not tested well.

The test will not compute then \code{length(x)<20}.

}

\value{
 If \code{FullInfo=TRUE} then function returns an object of the class \code{PvarBreakTest}. 
 It is the \code{list} that contains:
 \item{Stat}{a value of statistics (p-variation of transformed data).}
 \item{CriticalValue}{the critical value of the test according to significant level.} 
 \item{alpha}{the significant level.} 
 \item{p.value}{approximate p-value.} 
 \item{reject}{\code{logical}. If \code{TRUE}, the H0 was rejected.} 
 \item{dname}{the name of data vector.} 
 \item{p}{the power in p-variation calculus. The test performs only with the \code{p=4}.} 
 \item{x}{a vector of original data.}
 \item{y}{a vector of transformed data (\code{y=BridgeT(x)}).}
 \item{Timelabel}{time label of \code{x}. Used only for ploting.}
 \item{BreakPoints}{the indexes of break points suggestion.}
 \item{Partition}{a vector of indexes that indicates the partition of \code{y} that achieves the p-variation maximum.}
 \item{SplitPoints}{the indexes of \code{y} that indicates the points where \code{y} might be spited and analysed separately.}
 \item{info}{the list of extra information about the calculation.}
 
 The \code{info} is a list that contains:
 \item{PartitionN}{the length of the partition.}
 \item{AnalysedPointsN}{the length of the prepared \code{y}, i.e. after removing monotonic points and meaningless pints in small intervals.}
 \item{SegmentsN}{the number of segments the \code{y} was splinted.}
 \item{sizeN}{the inner \code{sizeN} parameter}
 \item{TakesTime}{the total time taken for calculation.}
 
}

\references{
[1] R. Norvaisa, A. Rackauskas. Convergence in law of partial sum processes in p-variation norm. 
   Lth. Math. J., 2008., Vol. 48, No. 2, 212-227.   
}

\author{
Vygantas Butkus
}


\seealso{
  Tests statistics is  \code{\link{pvar}} of the data \code{BridgeT(x)}(see  \code{\link{BridgeT}}) with (p=4).
  
  The critical value and the approximate  p-value of the test might by found by functions
  \code{\link{PvarQuantile}} and  \code{\link{PvarPvalue}}.	
 
}
\examples{



set.seed(1)
# x contains 4 blocks with diferent mean.
x = c(rnorm(250,0),rnorm(250,0.3), rnorm(250,0), rnorm(250,0.3))

# run test
PVBT = PvarBreakTest(x); PVBT
PVBTSummary = summary(PVBT); PVBTSummary

# the plot of the test
plot(PVBT)

# the sugestion of break points. True points are 250, 500 and 750
BreakPoints(PVBT)

# if we only interesting in main numbers.
PvarBreakTest(x, FullInfo = FALSE)

##################################
# the statistic maigth by find with pvar
Y = BridgeT(x)
pv = pvar(Y, p = 4)
pv

# critical value
PvarQuantile(n=length(x), prob = c(0.95))
# apricimal p-value
PvarPvalue(n=length(x), stat=pv$value)


}


