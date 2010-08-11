\name{stretchGene}
\alias{stretchGene}
\title{
"stretch" genes.
}
\description{
return a table with stretched body of the profile
}
\usage{
stretchGene(genetab, tssCut, lengthRef)
}
\arguments{
  \item{genetab}{
\code{data.frame}, the gene table (id, logratio, fittedval).
}
  \item{tssCut}{
\code{numeric},number of line considered as part of the TSS region.
}
  \item{lengthRef}{
  \code{numeric}, the length of the reference gene.
}
}
\details{
subfunction of \code{getOutput}.
}
\value{
data frame of "streched" or "compacted" body region of the gene.
}
\author{
Sankar Martial
}
\seealso{
See Also \code{\link{getOutput}}
}
