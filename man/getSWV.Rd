\name{getSWV}
\alias{getSWV}
\title{
get the Sliding Windows values by genes
}
\description{
\code{getSWV}, get the Sliding Windows Values by genes. Values are obtained using the median of the tags counts.
}

\usage{
getSWV(gID, wiggTable, commonTable, dbcurr, window = 250, span)
}

\arguments{
  \item{gID}{
  \code{numeric}, the Gene ID.
}
  \item{wiggTable}{
  \code{data.frame},the entire wiggTable.
}
  \item{commonTable}{
  \code{data.frame}, the results matrix of 'overlap' function.
}
  \item{dbcurr}{
  \code{data.frame}, the entire database.
}
  \item{window}{
  \code{numeric}, the sliding window size.
}
  \item{span}{
  \code{numeric}, if NULL take the span of the WIGG. 
}
}
\details{
subfunction of \code{getCodReg}.
}
\value{
return average inside the sliding window intensity.
}
\author{
Sankar Martial}

\seealso{
 See Also as \code{\link{getCodReg}}
}
