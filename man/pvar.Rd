\name{pvar}

\alias{pvar}
\alias{plot.pvar}
\alias{summary.pvar}
\alias{print.summary.pvar}

\title{
	p-variation calculation
}
\description{
	Calculates p-variation of the sample.
}
\usage{

pvar(x, p, TimeLabel=as.vector(time(x)), info = TRUE, sizeN = 7)

\method{summary}{pvar}(object, ...)

\method{plot}{pvar}(x, main = "p-variation", ylab = x$dname, 
	sub="p="\%.\%round(x$p,5)\%.\%", p-variation: "\%.\% formatC(x$value, 5, format = "f"), 
	col.PP = 3, cex.PP = 0.5, col.SP = 2, cex.SP = 1,	...)

}

\arguments{

  \item{x}{a (non-empty) numeric vector of data values or an object of the class \code{pvar}.}
  \item{p}{a positive number indicating the power p in p-variation.}
  \item{TimeLabel}{numeric, a time index of \code{x}. Used only for plotting.}
  \item{info}{\code{logical}. If \code{TRUE} (the default) the results will be an object of the class \code{pvar}
	with all the information included. Otherwise it will be only the p-variation value.}
  \item{sizeN}{a positive and odd number. This is inner parameter that may significantly influence the speed of calculation. 
    7 is recommended for randomized processes (for example Wiener process). }
  \item{object}{an object of class \code{pvar}.}  
  \item{main}{a \code{main} parameter in \code{plot} function.}
  \item{ylab}{a \code{ylab} parameter in \code{plot} function.}
  \item{sub}{a \code{sub} parameter in \code{plot} function.}
  \item{col.PP}{the color of partition points.}
  \item{cex.PP}{the cex of partition points.}
  \item{col.SP}{the color of split points.}
  \item{cex.SP}{the cex of split points.}
  \item{\dots}{further arguments.}	
}

\details{
To make it clear lets star with the definitions. Originally p-variation is defined for a functions.
For a function \eqn{f:[0,1] \rightarrow R}{f:[0,1] -> R} and \eqn{0 < p < \infty}{0 < p <  \infty} p-variation is defined as

\deqn{
	v_p(f) = \sup \left\{ \sum_{i=1}^m |f(t_i) - f(t_{i-1})|^p : 0=t_0<t_1<\dots<t_m=1 \right\}
}{
	v_p(f) = sup { \sum |f(t_i) - f(t_{i-1})|^p : 0=t_0<t_1<\dots<t_m=1}
}

Analogically for a sequences of values \eqn{X_0, X_1,..., X_n}, the p-variation is defined as
\deqn{
	v_p(\{X_i\}_{i=0}^n) = \max\left\{ \sum_{i=1}^k |X_{j_i}-X_{j_{i-1}}|^p: 0=j_0<j_1<\dots<j_k=n, \; k=1,2,...,n \right\} 
}{
	v_p(\{X_i\}_{i=0}^n) = max { \sum |X_{j_i}-X_{j_{i-1}}|^p :0=j_0<j_1<\dots<j_k=n, \; k=1,2,\dots,n }
}

The points \eqn{0=t_0<t_1<\dots<t_m=1} (or \eqn{0=j_0<j_1<\dots<j_k=n}) that achieves the maximums is called p-variation partition.

The functions\code{pvar} calculates p-variation for a sequence of values \eqn{X_0, X_1,\dots, X_n} 
that are written in the \code{numeric} vector \code{x}. 

This function also finds extra information that concerns the calculation of the p-variation. 
For example it also indicates the partition (\eqn{j_0, j_1,\dots, j_k} that achieves the maximum.

The variable \code{sizeN} is inner parameter that may influence the speed of calculation. Before the real maximisation
of the function \eqn{|X_{j_i}-X_{j_{i-1}}|^p} starts some of the points may be drooped out by analysing small intervals
(see  \code{\link{DropMidAll}}). 
The \code{sizeN} defines the length of 'small' interval. For randomised samples it is recommended to set \code{sizeN=7}.
But if \code{x} has special form the optimal \code{sizeN} value might be different.


}


\value{
 If \code{info=TRUE} then function returns the object of the class \code{pvar}. 
 It contains the \code{list} with information, namely:
 \item{value}{a value of p-variation.}
 \item{x}{a vector of original data \code{x}.}
 \item{p}{the value of p.}
 \item{dname}{a name of data vector.}
 \item{TimeLabel}{a time label of \code{x}.}
 \item{Partition}{a vector of indexes that indicates the partition that achieves the maximum.}
 \item{SplitPoints}{the indexes of \code{x} that indicates the points where \code{x} might be spited and analysed separately.}
 \item{info}{the list of extra information about the calculation.}
 
 The \code{info} is a vector that contains:
 \item{PartitionN}{the length of the partition.}
 \item{AnalysedPointsN}{the length of the prepared \code{x}, i.e. after removing monotonic points and meaningless points in small intervals.}
 \item{SegmentsN}{the number of segments the \code{x} was splinted.}
 \item{sizeN}{the inner \code{sizeN} parameter.}
 \item{TakesTime}{the total time taken for p-variation calculation. If it has negative value \code{-1}, 
	this means, that this information is not available.}
 
}

\references{
~~~ Some articles about p-variation.~~~
}
\author{
Vygantas Butkus
}

\note{
  The are additional fuctians that are now under development. 
  Anfortunately those functions are not relaible. But still migth be usefull.
  \code{\link{AddPvar}}
  \code{\link{PvarSubset}}
}

\section{Warning}{
	Unfortunately \code{pvar} is recursive function. Therefore, for large data sets this might by real limitation.
	For a randomised samples this function works just fine. The problem arise then \code{x} has very special kind of form.
}




\seealso{
 \code{pvar} maximises function \code{\link{Sum_p}}. 
 Firstly, it removes all monotonic points with function \code{\link{ChangePointsId}}.
 After that, some points are drop out by function \code{\link{DropMidAll}}.
 Moreover, vector is spited by \code{\link{SplitByExtremum}} and each part is calculated independently by \code{\link{pvarMon}}.

 In some points \code{pvar} objects migth be merged with function \code{\link{MergePvar}}.

 One of the p-variations application is realized in \code{\link{PvarBreakTest}} function. 
}
\examples{

set.seed(26)
### randomised data:
x = rbridge(1, 1000)

### the main functions:
pv = pvar(x, 2)
print(pv)
summary(pv)
plot(pv)


### By taking only change point, we get the same results:
CP = ChangePointsId(x)
x.CP = x[CP]
pv.CP = pvar(x.CP, TimeLabel=time(x)[CP], 2)
pv.CP == pv


### We can drop out some points withs any lost as well:
DM = DropMidAllB(x.CP, 2, dn=7)
x.DM = x.CP[DM] 
pv.DM = pvar(x.DM, TimeLabel=time(x)[CP][DM], 2)
pv.DM == pv

### We can witness the efectivness of each trick by the length of x
c(length(x), length(x.CP), length(x.DM))

### we olso can splti vector in independent samples by finding the extremums
SP = which.max(x)
pv.SP = MergePvar(pvar(x[1:SP],2), pvar(x[SP:length(x)],2))
pv.SP == pv


### Finly, all we must be concern is partition ponts.
pv; Sum_p(x[pv$Partition], 2)	#direct calculus

pv.PP = pvar(x[pv$Partition], TimeLabel=time(x)[pv$Partition], 2)
pv == pv.PP

op <- par(mfrow = c(2, 1))
	plot(pv, main="pvar with original data")
	plot(pv.PP, main="the same pvar without meaningless points")
par(op)


### the illustration of bad form for p-variation calculation
N = 2000
x = c(-N-1000, (1001:(1000+N-2))*((-1)^(0:(N-3))), N+1000)
plot(x, type = "p", pch=19)
\dontrun{
  # the recursion is too deep:
  var(x, p = 2)
}

### the problem might be solved by reversing x
rx = rev(x)
plot(rx, type = "p", pch=19)
pvar(rx, p = 2) #now it is OK



}


