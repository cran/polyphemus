\name{medichiMerge}
\alias{medichiMerge}
\title{
Merge MeDiChi files.
}
\description{
\code{medichiMerge}, merge all medichi files for one chromosome & return the list of peaks and their position.
}
\usage{
medichiMerge(files, windo, cutoff = 0.01, write = NA)
}
\arguments{
  \item{files}{
	\code{vector}, vector of String containing the path and name of each files.
}
  \item{windo}{
	\code{numeric}, the length of the window range
}
  \item{cutoff}{
	\code{numeric}, filers on the cutoff values of MediChi[1]
}
  \item{write}{
	\code{String}, path and file name for the output merged file. if NA, not output.
}
}
\details{
subfunction in \code{getInt}}

\author{
Sankar Martial}

\references{
[1] Reiss, D. J., M. T. Facciotti, et al. (2008). "Model-based deconvolution of genome-wide DNA binding." Bioinformatics 24(3): 396-403.
}
\seealso{
See Also \code{\link{getInt}}
}
