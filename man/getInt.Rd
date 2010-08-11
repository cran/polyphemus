\name{getInt}
\alias{getInt}
\title{
Extract intensities and start the comparative analysis.
}
\description{
\code{getInt} permits to extract intensities from files in one chromosome folder and start the comparative analysis.
}
\usage{
getInt(folderPath, medichiProc = c(NULL, NULL, NULL), prepost = c(500, 2000), window = 250, denomTrack, dataBaseFile, putTSS = FALSE, span = NULL)
}
\arguments{
  \item{folderPath}{
	\code{character}, path to the wig folder.
}
  \item{medichiProc}{
     \code{vector}, c(window, pvalue filters, String path and name for output (NULL default)). see the details section.
}
  \item{prepost}{
  		\code{vector}, contains two values, the number of bp (base pair) before and the number of bp after the coding region.
}
  \item{window}{
		\code{numeric}, size of the sliding window (in bp).
}
  \item{denomTrack}{
		\code{character}, the denominator track name (string contained in the filename).
}
  \item{dataBaseFile}{
		\code{character}, path and name of the database file.
}
  \item{putTSS}{
		\code{logical}, if TRUE, detect the putative TSS inside the coding region. (not yet implemented).
}
  \item{span}{
		\code{numeric}, if NULL take the span of the wig file else take the specified span.
		
}
}
\details{
  The \code{denomTrack} argument is the exact String of the denominator track contained in the filename. Its intensity values will be used as steady denominator when calculating the ratio between tracks.
  The \code{medichiProc} argument is a vector that contains three elements. The first is the window size where the TSS is centered. The second is a filter on the pvalue. The third is the file path and name of the output table. if NULL, the output file is not provided.
  The \code{getInt} function is also called by the wrapper \code{\link{gwComp}}.
}
\value{
  The \code{getInt} function returns an object (list) of, at least, four elements. It can be more depending on the number of tracks in the experiments.
  \item{gene}{
  	\code{character}, the gene annotation.
  }
  \item{chr}{
  	\code{character}, the chromosome name.
   }
  \item{denominator}{
  	\code{numeric}, intensities value for the denominator track.
   }
  \item{numerator}{
  	\code{numeric}, intensities value for the first numerator track.
   }
  \item{...}{ 
  	\code{numeric}, intensities value for the other track.
  }
}

\author{
Martial Sankar
}

\seealso{
See also \code{\link{gwComp}} for genome wide analysis (in all chromosomes). See also \code{\link{getCodReg}}.
}
\examples{
fpath1<-system.file("extdata", package="polyphemus")
fpath2<-file.path(file.path(fpath1,"ATRA"),"chr19")
x <- getInt(folderPath=fpath2, # path of the parent folder of the wig directory for the chr19
		medichiProc = c(300, 0.1, NULL), prepost = c(0, 0),
		window = 1000, denomTrack="EtOH", putTSS = FALSE,
		 dataBaseFile=file.path(fpath1,"mouserefseq-genenames.txt"), span = 500)
}
\keyword{intensity}
\keyword{count}

