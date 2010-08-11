\name{polyphemus-package}
\alias{polyphemus-package}
\alias{polyphemus}
\docType{package}
\title{
R package for comparative analysis of genome-wide RNA Polymerase II ChIP-seq data using non-linear normalisation approaches.
}
\description{
The \code{polyphemus} package was designed to extract comparative information about PolII action at different gene loci from ChIP-seq profiling. This permits to (i) integrates TSS annotation with PolII enrichment over coding regions, (ii) normalizes signal intensity profiles using non-linear approaches to correct for technical/experimental variations and iii) generates outputs for gene classification according to differential transcription characteristics.
}
\details{
\tabular{ll}{
Package: \tab polyphemus\cr
Type: \tab Package\cr
Version: \tab 0.3.3\cr
Date: \tab 2010-08-11\cr
License: \tab GPL (>=2)\cr
LazyLoad: \tab yes\cr
}

}
\author{
Martial Sankar directed by Marco Mendoza and Hinrich Gronemeyer
Maintainer: Martial Sankar <martial.sankar@unil.ch>
}
\references{
Reiss, D. J., M. T. Facciotti, et al. (2008). "Model-based deconvolution of genome-wide DNA binding." Bioinformatics 24(3): 396-403.
Rozowsky, J., G. Euskirchen, et al. (2009). "PeakSeq enables systematic scoring of ChIP-seq experiments relative to controls." Nat Biotechnol 27(1): 66-75.
Taslim, C., J. Wu, et al. (2009). "Comparative study on ChIP-seq data: normalization and binding pattern characterization." Bioinformatics 25(18): 2334-2340.
}

\keyword{ package }

\examples{
fpath1<-system.file("extdata", package="polyphemus")
fpath2<-file.path(fpath1,"ATRA")
ex <- gwComp(InitDir=fpath2,
		posProcParam=c(300, 0.1,NA), denomTrack="EtOH", 
		dataBase=file.path(fpath1,"mouserefseq-genenames.txt"),
		normMethod="quantile", stretched =50, TSSLength = 10, window=1000,
		 span=500,prepost=c(0,0), outIntTable=TRUE, plotMVA=TRUE) 
}






