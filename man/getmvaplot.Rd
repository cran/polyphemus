\name{getmvaplot}
\alias{getmvaplot}
\title{
Display the MVA plot}
\description{
\code{getmvaplot}, call by lowessNorm & quantileNorm function, display the mva plot}
\usage{
getmvaplot(matRaw, matNorm, normMeth = NULL)
}
\arguments{
  \item{matRaw}{
	\code{matrix}, initial matrix, contains the raw intensity
}
  \item{matNorm}{
	\code{matrix}, normalized matrix
}
  \item{normMeth}{
	\code{character}, normalization methods
}
}
\details{
subfunction of \code{quantileNorm} and \code{lowessNorm}.}

\author{
Sankar Martial}
\seealso{
See also \code{\link{quantileNorm}} and \code{\link{lowessNorm}}.
}
\keyword{ normalization }
