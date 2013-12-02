\name{BridgeT}

\alias{BridgeT}


\title{
	Bridge transformation
}
\description{
	Transforms data by Bridge trasnformation.
}
\usage{
BridgeT(x, normalise=TRUE)
}

\arguments{
  \item{x}{ a numeric vector of data values. The length must be greater then 1 and variance greater then 0.}
  \item{normalise}{\code{logical}, indicating whether the vector should be normalised. }
}

\details{
Let \code{n} denotes the length ox \code{x}. For  each \eqn{m \in [1,n]} bridge transformations \code{BridgeT}
is difined as  

\deqn{
	BridgeT(m, x) = \left\{ \sum_{i=1}^m x_i - \frac{m}{n} \sum_{i=1}^n x_i  \right\} .
}{
	BridgeT(m, x) = ( \sum^m x_i - m/n \sum^n x_i ).
}

Meanwhile, the transformation with normalisation is 

\deqn{
	BridgeT(m, x) = \frac{1}{\sqrt{n var(x)}} \left\{ \sum_{i=1}^m x_i - \frac{m}{n} \sum_{i=1}^n x_i  \right\} .
}{
	BridgeT(m, x) = ( \sum^m x_i - m/n \sum^n x_i  ) / (n var(x))^0.5.
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
Y = BridgeT(x)
Pv = pvar(Y, p = 4)
plot(Pv)

PvarBreakTest(x, FullInfo = FALSE)
}


