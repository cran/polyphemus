\name{getResult}
\alias{getResult}

\title{
final result output
}
\description{
the \code{getResult} function permits to get the final table output suitable for further clustering analysis. 
}
\usage{
getResult(objInt, iter = "", outIntTable = FALSE, out, norm = "quantile", stretchVal, tssLength, denoString)
}
\arguments{
  \item{objInt}{
		\code{list}, is the object output of the \code{\link{getNorm}} function. 
}
  \item{iter}{
		\code{character}. not used.
}
  \item{outIntTable}{
  		\code{logical}, if TRUE save the intermediate intensity tables.
  
}
  \item{out}{
		\code{character}, is the path of the output folder.
}
  \item{norm}{
		\code{character}, is the used normalization method.
}
  \item{stretchVal}{
		\code{numeric}, is the length value to use to uniformize the coding region length. See 'Details' section (units: number of span). 
}
  \item{tssLength}{
		\code{numeric}, is the approximate length of the TSS to used to separate TSS and body regions (units : number of span).
}
  \item{denoString}{
		\code{character}, is the name of the denominator track. 
}
}
\details{
	The \code{stretchVal} and \code{tssLength} arguments takes a numerical value in -number of span- units. The \code{span} length (in bp) is an argument of the \code{\link{gwComp}} and \code{\link{getInt}} functions. Users can easily determine the number of span length to define the \code{TSSlength}. if one assumes the TSS length is 1000 bp and he uses a span of 10 bp. The \code{TSSLength} arguments will be 100. The same method can be used to calculate the \code{stretchVal} argument.    
}
\value{
	returns one file by ratio between tracks. These output files can then reload in R or used in external program for clustering analysis.
}

\author{
Martial Sankar
}

\seealso{
See also \code{\link{gwComp}}, \code{\link{getInt}}, \code{\link{getNorm}}.
}

\keyword{ratio}
