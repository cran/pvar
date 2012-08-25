\name{pvar-package}
\alias{pvar-package}

\docType{package}
\title{
p-variation calculation and it's application
}
\description{
This is a package about the p-variation.
}
\details{
\tabular{ll}{
Package: \tab pvar\cr
Type: \tab Package\cr
Version: \tab 1.0.6\cr
Stage: \tab development \cr
Date: \tab 2012-09-28\cr
License: \tab GPL-2\cr
}

This package is about p-variation. It deals only of p-variation of a finite sample data values.
To make it clear lets star with the definitions. Originally p-variation is defined for a functions.

For a function \eqn{f:[0,1] \rightarrow R}{f:[0,1] -> R} and \eqn{0 < p < \infty}{0 < p <  \infty} 
p-variation is defined as

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

The points \eqn{0=t_0<t_1<\dots<t_m=1}(or \eqn{0=j_0<j_1<\dots<j_k=n}) that achieves the maximums is called p-variation partition.

There are two main functions that this package is all about, namely it is \code{\link{pvar}} and \code{\link{PvarBreakTest}}.
The main function in this package is \code{\link{pvar}}. 
It calculates the p-variation and the partition.
And the function \code{\link{PvarBreakTest}} is one of the examples of p-variation applications.
It performs structural break test of vector \code{x} that exams whether there are multiple
shifts in mean inside vector \code{x}.
The test statistics distribution if drown from Monte-Carlo simulations of 140 millions iterations. 
Therefore test's critical value and p-value are not exact, but onyly approximate .

 
All other functions are loaded only for supporting and illustrating purposes.

}

\note{
Up till now this package is still in the development stage. 
The written functions are well tested and I have strong believe that this package is already  reliable.

But the manual stays behind, because the lack of time and motivation.
Sorry for poor English and lack of information. 
I hope this imperfection will be eliminated in some future time.

}



\author{
Author and Maintainer: Vygantas Butkus <Vygantas.Butkus@gmail.com>.

Special thanks to Rimas Norvaisa - the superviser of my studies.
}
\references{
~~ Literature or other references for background information ~~
}




\seealso{
The main function is \code{\link{pvar}} - it finds p-variation and the partition that maximises \code{\link{Sum_p}} function.

Other important functions is \code{\link{PvarBreakTest}} it performs structural break test of vector \code{x}
by calculating p-variations of \code{BBT(x)} (see \code{\link{BBT}}). 
The critical value and the approximate  p-value of the test might by found by functions
  \code{\link{PvarQuantile}} and  \code{\link{PvarPvalue}}.
}

\examples{

set.seed(1)

######### pvar ########
x = rbridge(1, 1000)

### simple run:
pv = pvar(x, 2)
print(pv)
summary(pv)
plot(pv)

### partition points are those that maximises Sum_p function
Sum_p(x[pv$Partition], 2)

### lets take only the partition points:
pv2 = pvar(x[pv$Partition], TimeLabel=pv$TimeLabel[pv$Partition], 2)

op <- par(mfrow = c(2, 1))
	plot(pv, main="pvar with original data")
	plot(pv2, main="the same pvar without meaningless points")
par(op)

######### PvarBreakTest ########
x = c(rnorm(250,0),rnorm(250,0.3), rnorm(250,0), rnorm(250,0.3))

# Run test
PVBT = PvarBreakTest(x); PVBT
PVBTSummary = summary(PVBT); PVBTSummary

# the plot of the test
plot(PVBT)

# the suggestion of break points. True points are 250, 500 and 750
BreakPoints(PVBT)


##################################
# the statistic maigth by find with pvar
Y = BBT(x)
pv = pvar(Y, p = 4)
pv

# critical value
PvarQuantile(n=length(x), prob = c(0.95))
# apricimal p-value
PvarPvalue(n=length(x), stat=pv$value)

PVBT


}
