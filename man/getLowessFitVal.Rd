\name{getLowessFitVal}
\alias{getLowessFitVal}
\title{
get fitted values from a lowess regression.
}
\description{
\code{getLowessFitVal}, return a list of table where it was extracted the fit values from a lowess regression,
}
\usage{
getLowessFitVal(id, list, toPlot = FALSE)
}
\arguments{
  \item{id}{
	\code{character}, geneId from the final output of polyphemus
}
  \item{list}{
	\code{character}, splitted polyphemus output (by ID)
}
  \item{toPlot}{
	\code{Boolean}, whether the the fitted values should be added to the mva plot
}
}
\details{
subfunction of getOutput}
\author{
Sankar Martial}
\seealso{
See also \code{\link{getOutput}}.}
