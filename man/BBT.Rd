\name{BBT}

\alias{BBT}


\title{
	Brownian Bridge transformation
}
\description{
	Transforms data according to empirical Brownian Bridge precess.
}
\usage{
BBT(x)
}

\arguments{
  \item{x}{ a numeric vector of data values. The length must be greater then 1 and variance greater then 0.}
}

\details{
Let \code{n} denotes the length ox \code{x}. For  each \eqn{m \in [1,n]} transformations \code{BBT} difined as  

\deqn{
	BBT(m, x) = \frac{1}{\sqrt{n var(x)}} \left\{ \sum_{i=1}^m x_i - \frac{m}{n} \sum_{i=1}^n x_i  \right\} .
}{
	BBT(m, x) = ( \sum^m x_i - m/n \sum^n x_i  ) / (n var(x))^0.5.
}


}

\value{
A numeric vector. 
}


\seealso{
 \code{\link{PvarBreakTest}},  \code{\link{rbridge}}
}
\examples{

x = c(rnorm(250,0),rnorm(250,0.3), rnorm(250,0), rnorm(250,0.3))
Y = BBT(x)
Pv = pvar(Y, p = 4)
plot(Pv)

PvarBreakTest(x, FullInfo = FALSE)
}


