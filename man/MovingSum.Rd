\name{MovingSum}

\alias{MovingSum}

\title{
Moving sum
}
\description{
The filter of moving sum.
}
\usage{
MovingSum(x, n = 2)
}

\arguments{
  \item{x}{a numeric vector.}
  \item{n}{a number indicating the size of moving window.}
}

\value{
The numeric vector.
}

\seealso{
\code{\link{filter}}
}
\examples{

MovingSum(1:5)	# the results: 3  5  7  9


}



