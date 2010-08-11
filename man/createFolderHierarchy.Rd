\name{createFolderHierarchy}
\alias{createFolderHierarchy}
\title{
Create hierarchy of folder.
}
\description{
\code{createFolderHierarchy} create a hierarchy of folders. The functions takes in input a folder of pooled files (bed and wig format). One folder is then created for each chromosome. the folder name corresponds to the name of the chromosome (i.e: "chr1", "chrX"\dots). The subfolders "pos" and "wig" are then created. They contain respectively the bed and wig files.  
}
\usage{
createFolderHierarchy(initialDir)
}
\arguments{
  \item{initialDir}{ \code{character} the path of the folder containing the pooled files.
}
}
\details{
If the folders and subfolders are already created nothing is done.
}

\author{
Martial Sankar
}

\examples{
fpath<-system.file("extdata", package="polyphemus")
createFolderHierarchy(file.path(fpath,"ATRA"))
}

\keyword{hierarchy}
