\name{getIncrByG}
\alias{getIncrByG}
\title{
sub-function of generateWig
}
\description{
\code{getIncrByG}, sub-function of \code{generateWig}.
}
\usage{
getIncrByG(gCurr, col, countdf, dbori, window, span, prespost)
}
\arguments{
  \item{gCurr}{
	\code{String}, gene name.
}
  \item{col}{
	\code{numeric}, column index.
}
  \item{countdf}{
	\code{data.frame}, two columns dataset (Start, End position).
}
  \item{dbori}{
	\code{data.frame}, intial database (refseq).	
}
  \item{window}{
 	 \code{numeric}, size of the sliding window.
}
  \item{span}{
 	 \code{numeric}, span of the sliding window.
}
  \item{prespost}{
  \code{vector}, the number of pb before and after the gene, \code{c(int, int)}.
}
}
\details{
\code{generateWig}
}

\author{
Sankar Martial
}
\keyword{ wig }
