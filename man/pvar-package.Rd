\name{pvar-package}
\alias{pvar-package}


\encoding{utf8} 
 

\docType{package}
\title{
p-variation calculation and application
}
\description{
This package deals with p-variation for the sequence of data values. It gives opportunity to calculate the p-variation -- this is the main purpose of this package.
Moreover, it includes one example of practical application of the p-variation.

}
\details{
\tabular{ll}{
Package: \tab pvar\cr
Type: \tab Package\cr
Version: \tab 1.0.7\cr
Date: \tab 2013-12-02\cr
License: \tab GPL-2\cr
Institution: \tab Vilnius University Faculty of Mathematics and Informatics \cr
}

This package is about p-variation. It deals only of p-variation of a finite sample data values.
To be precise, lets star with the definitions. Originally p-variation is defined for a functions.

For a function \eqn{f:[0,1] \rightarrow R}{f:[0,1] -> R} and \eqn{0 < p < \infty}{0 < p <  \infty} 
p-variation is defined as

\deqn{
	v_p(f) = \sup \left\{ \sum_{i=1}^m |f(t_i) - f(t_{i-1})|^p : 0=t_0<t_1<\dots<t_m=1 \right\}
}{
	v_p(f) = sup { \sum |f(t_i) - f(t_{i-1})|^p : 0=t_0<t_1<\dots<t_m=1}
}

Analogically, for a sequences of values \eqn{X_0, X_1,..., X_n}, the p-variation is defined as
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
 
All other functions are loaded only for supporting and illustrating purposes.

}



\author{
Author and Maintainer: Vygantas Butkus <Vygantas.Butkus@gmail.com>.

Special thanks to Rimas Norvaiša the superviser of my studies.
}
\references{
[1] R. M. Dudley, R. Norvaisa. An Introduction to
  p-variation and Young Integrals, Cambridge, Mass., 1998.
  
[2] R. M. Dudley, R. Norvaisa. Differentiability of 
  Six Operators on Nonsmooth Functions and p-Variation, Springer Berlin Heidelberg, 
  Print ISBN 978-3-540-65975-4, Lecture Notes in Mathematics Vol. 1703, 1999.
  
[3] R. Norvaisa, A. Rackauskas. Convergence in law of partial sum processes in p-variation norm. 
   Lth. Math. J., 2008., Vol. 48, No. 2, 212-227.   
   
[4] J. Qian. The p-variation of Partial Sum Processes and the Empirical Process. 
   The Annals of Probability, 1998, Vol. 26, No. 3, 1370-1383.
}


\seealso{
The main function is \code{\link{pvar}} - it finds p-variation and the partition that maximises \code{\link{Sum_p}} function.

Other important functions is \code{\link{PvarBreakTest}} it performs structural break test of vector \code{x}
by calculating p-variations of \code{BridgeT(x)} (see \code{\link{BridgeT}}). 

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

# run test
PVBT = PvarBreakTest(x); PVBT
PVBTSummary = summary(PVBT); PVBTSummary

# the plot of the test
plot(PVBT)

# the suggestion of break points. True points are 250, 500 and 750
BreakPoints(PVBT)


##################################
# the statistic maigth by find with pvar
Y = BridgeT(x)
pv = pvar(Y, p = 4)
pv

# critical value
PvarQuantile(n=length(x), prob = c(0.95))
# apricimal p-value
PvarPvalue(n=length(x), stat=pv$value)

PVBT


}
