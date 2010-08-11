\name{getNorm}
\alias{getNorm}
\title{
	Normalization
}
\description{
	The \code{getNorm} permits the normalization of the intensity (count) values. 
}
\usage{
getNorm(objInt, norm = c("linear", "lowess", "quantile"), lowessSpan = NULL, mva = NULL)
}
\arguments{
  \item{objInt}{
		\code{list}, the intensity object. Output of the \code{\link{getInt}} function.
}
  \item{norm}{
		\code{character}, the normalization method : "linear" for linear normalization, "lowess" and "quantile" for local normalization. See details section.
}
  \item{lowessSpan}{
		\code{numeric}, the smoother span. This gives the proportion of points in the plot which influence the smooth at each value. Larger values give more smoothness but requires more computing time. if NULL, lowess span = 2/3.
}
  \item{mva}{
  		\code{logical}, if TRUE, the MVA figures are plotted.
}
}
\details{
	\code{Polyphemus} provides three normalization procedure. Choosing a "linear" normalization, a correction factor is applied such that the sum of the intensities for a track is the same for all the tracks. Choosing a local normalisation, "lowess" or "quantile", a correction factor is locally calculating and applied to the track intensities. This turns to be more reliable than classic "linear" normalisation. \emph{Taslim et al.} already mentionned it to normalize PolII ChIP-seq data. However when dealing with more than two samples, the "quatile" procedure turns to be more suitable as it requires less computing time. The "lowess" normalization used the \code{\link{lowess}} implementation of the R base package. The "quantile" procedure used the \code{normalizeQuantiles} of the limma package. \code{getNorm} is called in the polyphemus wrapper \emph{gwComp}.  
}
\value{
	returns an object (list) containing three elements.
	\item{gene}{
		\code{character}, the gene annotation.
	}
	\item{chr}{
		\code{chr}, the chromosome name. 
	}
	\item{matrix}{
		\code{matrix},the matrix of normalized intensity values.
	}
}
\references{
\emph{Taslim, C., J. Wu, et al. (2009).} "Comparative study on ChIP-seq data: normalization and binding pattern characterization." Bioinformatics 25(18): 2334-2340.
}
\author{
Martial Sankar
}

\seealso{
 See Also \code{\link{gwComp}}.
}

\keyword{normalization}
\keyword{normalisation}
