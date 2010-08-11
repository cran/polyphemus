\name{generateWig}
\alias{generateWig}

\title{
Generate Wig files
}
\description{
\code{generateWig} Generate Wig files from a count files [deprecated]
}
\usage{
generateWig(countdf, col, chr, dbori, window, span, prespost = c(1000, 1000), out)
}
\arguments{
  \item{countdf}{
	\code{data.frame}, two columns dataset (Start, End position)
}
  \item{col}{
	\code{numeric}, column index of the count column 
}
  \item{chr}{
	\code{numeric}, chromosome index

}
  \item{dbori}{
	\code{data.frame}, intial database (refseq)	
}
  \item{window}{
 	 \code{numeric}, size of the sliding window
}
  \item{span}{
 	 \code{numeric}, span of the sliding window
}
  \item{prespost}{
  \code{vector}, the number of pb before and after the gene, \code{c(int, int)}
}
  \item{out}{
	  \code{character}, output path and file name of the generated wwig file
}
}
\details{
deprecated function.
}

\author{
Sankar Martial
}
\keyword{ Wig }
