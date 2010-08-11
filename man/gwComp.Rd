\name{gwComp}
\alias{gwComp}

\title{
Polyphemus Wrapper
}
\description{
The \code{gwComp} function is the polyphemus wrapper. This calls the sub-fonctions requires to performs the comparative analysis of the ChIP-seq profiles (multi-states analysis) in a genome wide manner (whole genome analysis). 
}
\usage{
gwComp(InitDir, posProcParam = c(NULL, NULL, NULL), prepost = c(500, 2000), window = 250, denomTrack, dataBase, putTSS = FALSE, span = NULL, normMethod = "quantile", stretched = "median", TSSLength = NULL, outIntTable = FALSE, plotMVA = FALSE, logFile = NULL)
}
\arguments{
  \item{InitDir}{
 	 	\code{character} the path of the folder containing the pooled files. If the folder already contained hierarchized folder. nothing is done.
}
  \item{posProcParam}{
  	     \code{vector}, are the post-processing parameters c(window, pvalue filters, String path and name for output (NULL default)).	
}
  \item{prepost}{
  		\code{vector}, contains two values, the number of bp (base pair) before and the number of bp after the coding region.
}
  \item{window}{
  		\code{numeric}, size of the sliding window (in bp).
}
  \item{denomTrack}{
  		\code{character}, the denominator track name (string contained in the filename)
}
  \item{dataBase}{
  \code{character}, path and name of the database file.
}
  \item{putTSS}{
  	\code{logical}, if TRUE, detect the putative TSS inside the coding region. (not yet implemented).
}
  \item{span}{
  	\code{numeric}, if NULL take the span of the wig file else take the specified span.
}
  \item{normMethod}{
  		\code{character}, the normalization method : "linear" for linear normalization, "lowess" and "quantile" for local normalization.
}
  \item{stretched}{
  		\code{numeric}, is the length value to use to uniformize the coding region length. See 'Details' section (units: number of span). 
}
  \item{TSSLength}{
  		\code{numeric}, is the approximate length of the TSS to used to separate TSS and body regions (units : number of span).
}
  \item{outIntTable}{
  		\code{logical}, if TRUE save the intermediate intensity tables.
}
  \item{plotMVA}{
  		\code{character}, path and name to save the MVA plots. if NULL, the MVA figures are not plotted.
}
  \item{logFile}{
  		\code{logical}, if TRUE a log file is save in the current folder.
}
}
\details{
	The \code{denomTrack} argument is the exact String of the denominator track contained in the filename. Its intensity values will be used as steady denominator when calculating the ratio between tracks. The \code{posProcParam} argument is a vector that contains three elements. The first is the window size where the TSS is centered. The second is a filter on the pvalue. The third is the file path and name of the output table. if NULL, the output file is not provided.
	The \code{Polyphemus} package provides three normalization procedure. Choosing a "linear" normalization, a correction factor is applied such that the sum of the intensities for a track is the same for all the tracks. Choosing a local normalisation, "lowess" or "quantile", a correction factor is locally calculating and applied to the track intensities. This turns to be more reliable than classic "linear" normalisation. \emph{Taslim et al.} already mentionned it to normalize PolII ChIP-seq data. However when dealing with more than two samples, the "quatile" procedure turns to be more suitable as it requires less computing time. The "lowess" normalization used the \code{\link{lowess}} implementation of the R base package. The "quantile" procedure used the \code{normalizeQuantiles} of the limma package.  
	The \code{stretched} and \code{tssLength} arguments takes a numerical value in -number of span- units. The \code{span} length (in bp) is an argument of the \code{\link{gwComp}} and \code{\link{getInt}} functions. Users can easily determine the number of span length to define the \code{TSSlength}. if one assumes the TSS length is 1000 bp and he uses a span of 10 bp. The \code{TSSLength} arguments will be 100. The same method can be used to calculate the \code{stretched} argument.
}
\value{
	returns one file by ratio between tracks. These output files can then reload in R or used in external program for clustering analysis.
}
\references{
\emph{Taslim, C., J. Wu, et al. (2009).} "Comparative study on ChIP-seq data: normalization and binding pattern characterization." Bioinformatics 25(18): 2334-2340.
}
\author{
Martial Sankar
}
\seealso{
See also the wrapped functions : \code{\link{getInt}}, \code{\link{getNorm}}, \code{\link{getResult}}. 
}
\examples{
fpath1<-system.file("extdata", package="polyphemus")
fpath2<-file.path(fpath1,"ATRA")
ex <- gwComp(InitDir=fpath2,
		posProcParam=c(300, 0.1,NA), denomTrack="EtOH", 
		dataBase=file.path(fpath1,"mouserefseq-genenames.txt"),
		normMethod="quantile", stretched =50, TSSLength = 10, window=1000,
		 span=500,prepost=c(0,0), outIntTable=TRUE, plotMVA=TRUE) 
}
\keyword{wrapper}
