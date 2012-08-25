\name{PvarExtraFun}

\alias{MergePvar}
\alias{PvarSubset}



\title{
	Extra functions for \code{pvar} objects.
}
\description{
	Some extra functions for \code{pvar} objects that might by useful.
}
\usage{

PvarSubset(pv, a=1, b=length(pv$x), dname=NULL)

MergePvar(..., JoinIfPossible = TRUE, dname = NULL)


}

\arguments{
  \item{\dots}{\code{pvar} objects that should be merged.}
  \item{JoinIfPossible}{\code{logical}. If \code{TRUE} (the default) the ending points of 
	diferent vectors \code{pvar$x} will be treated as one point in result vector. Points will be joint only if they has the same value.}
  \item{dname}{the name of data vector of new \code{pvar} object.}
  \item{pv}{a object of \code{pvar}.}
  \item{a}{a number indicating the start of the subset.}
  \item{b}{a number indicating the end of the subset.}
}

\details{
The function \code{PvarSubset} finds the p-variation of the subset of \code{pv$x}. 

The function \code{MergePvar} merges \code{pvar} objects. It does not do any recalculations, it just combine the information into one object.
Therefore p-variation value form \code{MergePvar} function not necessary real p-variation. 
The user must use this function with careful and some knowledge: by p-variation properties you can merge p-variation straightforward
if joint point remains in the partition of joint set. The extremum points are such points for sure. See examples for illustration.

}

\value{
a object of the class \code{pvar}.
}


\seealso{
 \code{\link{pvar}}
}

\examples{

#### pvSubset ####
set.seed(12)
x = BBT(rnorm(1000))

pv = pvar(x, 2)	

pvSubset = PvarSubset(pv, 200, 800)
pvSmall = pvar(x[200:800], 2)
pvSubset == pvSmall


op = par(mfrow = c(2, 1))
	plot(pvSubset)
	plot(pvSmall, xlim=range(pvSmall$TimeLabel), ylim=range(pvSmall$x))
par(op)

#### MergePvar ####
set.seed(13)
x = BBT(rnorm(1000))
S = Split_MinMax(x)

pv = pvar(x, 2)	

pv1 = pvar(x[1:S[2]], 2)	 
pv2 = pvar(x[S[2]:length(x)], 2)	 

pvJ = MergePvar(pv1, pv2)
pvJ == pv	# since joint point is extremum, we get the correct value

op = par(mfrow = c(2, 1))
	plot(pv)
	plot(pvJ, xlim=range(pv$TimeLabel), ylim=range(pv$x))
par(op)

# on the other hand if we take other points we usually do not get correct value
pvJ = MergePvar(pvar(x[1:500], 2), pvar(x[500:length(x)], 2))
pvJ == pv




}

