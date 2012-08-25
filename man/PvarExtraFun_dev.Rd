\name{PvarExtraFun_dev}


\alias{Add2Pvar}
\alias{AddPvar}


\title{
	Extra functions for p-variation that are under development.
}
\description{
	Some useful functions for objects of class \code{pvar}. Those functions now under development. 
	It is known that those functions not correct therefor it should not be used by users (at least up until now).
}
\usage{

Add2Pvar(pv1, pv2, JoinIfPossible = TRUE, dname = NULL)

AddPvar(..., JoinIfPossible = TRUE, dname = "data")

}

\arguments{
  \item{\dots}{\code{pvar} objects that should be added.}
  \item{JoinIfPossible}{\code{logical}. If \code{TRUE} (the default) the ending points of 
	diferent vectors \code{pvar$x} will be treated as one point in result vector. Points will be joint only if they has the same value.}
  \item{dname}{The name of data vector of new \code{pvar} object.}
  \item{pv1, pv2}{a object of \code{pvar}.}
}

\details{
The initial idea was to write the function that adds \code{pvar} object (with recalculations). 
The purpose is to avoid complete recalculation. It is meant to use the information already available.
Unfortunately, the way functions was tritten is no correct, no optimal. It should be reconsider.

}

\value{
a object of the class \code{pvar}.
}


\seealso{
 \code{\link{pvar}}, \code{\link{MergePvar}}.
}
\examples{

#### in this case it is OK:
	set.seed(20)
	XX = BBT(rnorm(1000))
	pv = pvar(XX, p = 2)
	pv1 = pvar(XX[1:round(length(XX)/2)], p = 2)
	pv2 = pvar(XX[round(length(XX)/2):length(XX)], p = 2)
	
	pvJ = AddPvar(pv1, pv2);pvJ
	Sum_p(pvJ$x[pvJ$Partition],pvJ$p)
	pv
	pvJ == pv	

	
	op = par(mfrow = c(2, 1))
		plot(pv)
		plot(pvJ, xlim=range(pv$TimeLabel), ylim=range(pv$x))
	par(op)
	
#### in this case it is NOT OK:
	set.seed(30)
	XX = BBT(rnorm(1000))
	pv = pvar(XX, p = 2)
	pv1 = pvar(XX[1:round(length(XX)/2)], p = 2)
	pv2 = pvar(XX[round(length(XX)/2):length(XX)], p = 2)
	
	pvJ = AddPvar(pv1, pv2);pvJ
	Sum_p(pvJ$x[pvJ$Partition],pvJ$p)
	pv
	pvJ == pv	
	
	op = par(mfrow = c(2, 1))
		plot(pv)
		plot(pvJ, xlim=range(pv$TimeLabel), ylim=range(pv$x))
	par(op)
	
}

\keyword{ internal }
