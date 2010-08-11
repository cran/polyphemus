\name{lowessNorm}
\alias{lowessNorm}
\title{
lowess normalization}
\description{
\code{lowessNorm}, perform the lowess normalisation. Called by the \code{getNorm} function.

}
\usage{
lowessNorm(listInt, ff = 2/3, plotMVA = NULL)
}
\arguments{
  \item{listInt}{
	\code{numeric}, vector of intensity (SWI) values.
}
  \item{ff}{
	\code{numeric}, the span value for lowess. see \code{lowess} help.
}
  \item{plotMVA}{
	\code{Boolean}, whether to display the MVA plot (on the X).
}
}
\details{
subfunction of \code{getNorm}.
}


\author{
Sankar Martial
}

\seealso{
See Also \code{\link{lowess}}, \code{\link{getNorm}}
}
\keyword{ normalization }
