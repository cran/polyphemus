\name{getCodReg}
\alias{getCodReg}
\title{
extract intensity (count) values from the wig file.
}
\description{
\code{getCodReg} permits to define coding regions and extract intensity (count) values from the wig file.
}
\usage{
getCodReg(wiggFile, prepost, window, dbCurr, chrCurr, span = NULL)
}
\arguments{
  \item{wiggFile}{ 
  		\code{character}, path and filename to a wig file.
}
  \item{prepost}{ 
  		\code{vector}, contains two values, the number of bp (base pairs) before and the number of bp after the coding region.
}
  \item{window}{
		\code{numeric}, size of the sliding window (in bp).
}
  \item{dbCurr}{
		\code{data.frame}, corresponding to a subset of the database file for the processing chromosome.
}
  \item{chrCurr}{
        \code{character}, name of the processing chromosome.
}
  \item{span}{
		\code{numeric}, if NULL take the span of the wig file else take the specified span.
}
}

\value{
returns a list. each element corresponds to the vector of intensities for the specific gene. 
}

\author{
Martial Sankar
}



