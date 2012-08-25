\name{ConcatenateStrings}
\alias{\%.\%}

\title{
Concatenate strings
}
\description{
Concatenate Strings
}
\usage{
x \%.\% y
}

\arguments{
  \item{x, y}{strings that should be concatenated. }

}
\details{
The same result may be achieved with \code{paste}, but this function is more user friendly then you need to add to strings.
}
\value{
A character string of the concatenated values.
}


\seealso{
 \code{\link{paste}}
}
\examples{
paste("I ", "love ", "R.", sep="")
"I " \%.\% "love " \%.\% "R."

x = c(2,1,6,7,9)
paste("The length of vector (", paste(x , sep="", collapse =","), ") is ", length(x) , sep="")
"The length of vector (" \%.\% paste(x , sep="", collapse =",") \%.\% ") is " \%.\% length(x)


}

\keyword{ print }

