# Polyphemus version 3.3 dev
#
#- Contain a collection of Polyphemus functions.
#- It is the sample comparison step. It will process a WIGG file
#- containing the count.
#- impovements of the version .2: - make a folder hierarchy a genome wide.


#-------------- Local --------------------#
library(IRanges)
library(zoo)
library(RUnit)
library(parallel)
library(R.utils)
library(limma)

#-----------------------------------------#

#library(Cairo)

cat("Polyphemus 0.3.4 \n")


# create the wrapper
gwComp <- function (InitDir, posProcParam=c(NULL, NULL, NULL), prepost=c(500,2000), window=250, denomTrack, dataBase, putTSS=FALSE, span=NULL, normMethod="quantile", stretched="median",TSSLength=NULL, outIntTable=FALSE, plotMVA=FALSE, logFile=NULL){
	if(	is.null(logFile)){
		toLog <- paste(InitDir,"/","log.txt",sep="")
		}
	#create folder structure
	createFolderHierarchy(InitDir)

	#Iter on each Chr folder
	chrFolderList <- dir(InitDir, full.names=TRUE)

	#Global
	InitDir <<- InitDir
	normMethod<<-normMethod
        browser()
	cat("******* extract  SWI ******\n")
	intsbyChr <- mclapply(chrFolderList, getInt, medichiProc= posProcParam, prepost= prepost, window= window, span=span, denomTrack=denomTrack, dataBaseFile= dataBase, putTSS= putTSS)

	#assign names (chr)
	cat("******* get Chr names ******\n")
	nm <- NULL
	index2rm <- NULL
	for (i in 1:length(intsbyChr)){
		if(all(unlist(lapply(intsbyChr[[i]][3:length(intsbyChr[[i]])],is.numeric)))){
			nm <- c(nm,intsbyChr[[i]][2])
		}else {
			log <- which(unlist(lapply(intsbyChr[[i]][3:length(intsbyChr[[i]])],is.numeric))==FALSE)
			write(paste(intsbyChr[[i]][2+log],"\n"), logFile, append=FALSE, ncolumns = 1)
			cat(logFile,"\n")
			index2rm <- i
			nm <- c(nm,intsbyChr[[i]][2])
			}
		}
	nm <- as.character(unlist(nm))
	names(intsbyChr) <- nm
	#norm
	if(length(index2rm)>0){
   		intsbyChr <- intsbyChr[-index2rm]
     }
     if(length(intsbyChr)==0){
     	stop(paste("[gwComp] Error: No sample processed, see Log file :", logFile,"\n"))
     	}

   cat("******* normalize *******\n")
	if (!is.null(normMethod)){
	normIntsbyChr <- mclapply(intsbyChr, getNorm, norm=normMethod, mva=plotMVA)
	}else {
		normIntsbyChr <- intsbyChr
		}

	#normVal <<- normIntsbyChr
	cat("******* get output *******\n")
	resbyChr <- mclapply(normIntsbyChr, getResult, outIntTable= outIntTable,
	out= InitDir, norm=normMethod, denoString= denomTrack, stretchVal=stretched, tssLength= TSSLength)

	}



# getInt: (by Chr) Extract intensities from files in the directory and start the comparative analysis.
#           Takes in parameters :
#           - wiggDir, String,a path to a directory cotaining the wigg files to analyse,
#           - prepost,vector, containing the number of pb before and after the coding region
#           - denomTrack, String, the denominator track name (string contained in the filename)
#           - dataBaseFile, String,pathname to  the organism data base file
#           - putTSS, Logical, detect putative TSS
#           - medichiDir, String, folder of medichi coeficients
#           - medichiProc, vector, c(windows, pvalue filters, String path&name for output (NULL default))
#           - span, Integer, if NULL take the WIGG span.
getInt <- function (folderPath, medichiProc=c(NULL,NULL,NULL),
                      prepost=c(500,2000),window=250, denomTrack, dataBaseFile,
                       putTSS=FALSE,span=NULL){


  wiggDir <- paste(folderPath,"/wig", sep="")
  medichiDir <- paste(folderPath, "/pos", sep="")
  wiggs <- dir(wiggDir, ".wig", full.names=TRUE)
  wiggstoSearch <- dir(wiggDir, ".wig", full.names=FALSE)
  indexDenom <- which(wiggstoSearch== grep(denomTrack, wiggstoSearch,value=TRUE))
  #checkEquals(length(wiggs), length(dir(medichiDir))) # test
if(length(indexDenom)==0){
	stop(paste("no denom Track in ", wiggDir, "\n",sep=""))
	}
  #- extract Current Chr:
  chr <- strsplit(folderPath, "/")[[1]]
  chr <- chr[length(chr)]
  cat("start", chr, "\n")
  dbori <- read.table(dataBaseFile,as.is=FALSE)
  dbi <- dbori[dbori[,1]==chr,,drop=F] # db for the current chr
  dbmatrix <- data.matrix(dbi)# 2 <=> +; 1<=>-
  db <- dbmatrix
  rm(dbori)
  if(nrow(dbi)==0){
  	stop(paste("chr name *", chr, "* not referenced in ", dataBaseFile, "\n",sep=""))
  	}
  #- cross with medichi identified peaks

  db[,2] <- db[,2]-prepost[1]
  db[,3] <- db[,3]+prepost[2]


  irDb <- IRanges(db[,2], db[,3])


  #- merge MeDiChI files
  if (!is.null(medichiDir)){
    medichis <- dir(medichiDir,full.names=TRUE)
    mergedMedichi <- medichiMerge(medichis, windo =as.numeric(medichiProc[1]) ,
                                  cutoff=as.numeric(medichiProc[2]),
                                  write=medichiProc[3])#merge site
    #- reduce the dbCurr,
    irMM <- IRanges(mergedMedichi[,1], mergedMedichi[,2])
    ov <- findOverlaps(irDb, irMM)
    #possibleG <- unique(as.matrix(ov)[,2])
    possibleG <- unique(as.matrix(ov)[,1])
    db <- db[sort(possibleG),]
    rm(ov,irMM,mergedMedichi,medichis,irDb)
  }


  #- get codingwindow=250 region & BW intensity
  listOfGbyTrack <- lapply(wiggs, getCodReg, prepost=prepost,window=window,
                           dbCurr=db, chrCurr=chr,span=span)

  #checkEquals(length(listOfGbyTrack[[1]]), length(listOfGbyTrack[[2]])) #test


  denomTrack <- listOfGbyTrack[[indexDenom]]#get the denom Track
  vectIndex <- 1:length(wiggs)
  vectIndex <- vectIndex[-indexDenom]
  expName <- wiggstoSearch[vectIndex]

  obj <- list()

  #- take shorter length track to intersect with the denom track
  i <- vectIndex[iii <- which.min(unlist(lapply(vectIndex,function(x){length(names(listOfGbyTrack[[x]]))})))]
  #i <- vectIndex[1]
  cId <- intersect(as.integer(names(denomTrack)), as.integer(names(listOfGbyTrack[[i]])))
  sortedCId <- as.character(sort(cId)) #-- order gene
  refCId <- sortedCId
  denomTrackA <- denomTrack[refCId]
  trackB <- listOfGbyTrack[[i]][refCId]
  # annotate
  nmID<- db[as.numeric(names(denomTrackA)),4]
  lines<-which(dbmatrix[,4]%in%nmID)# get line
  nm<- dbi[lines,4]                 #get gene name

	nnn<-1:length(denomTrackA)

    colGID <- unlist(lapply(nnn, function(x) rep(as.character(nm[x]),length(denomTrackA[[x]]))))
    colG<-colGID
    chrIntTrackAdenom <- as.numeric(unlist(denomTrackA))+1
    chrIntTrackB <- as.numeric(unlist(trackB))+1
    checkEquals(length(chrIntTrackAdenom), length(chrIntTrackB))

  obj$name <- colG
  obj$chr <- chr
  obj$denominator <- as.numeric(chrIntTrackAdenom)
  obj$numerator <- as.numeric(chrIntTrackB)
  iter <- length(obj)+1
  vectIndex1 <- vectIndex[-iii]

  for (i in vectIndex1) {
  	g <- as.integer(names(listOfGbyTrack[[i]]))
  	cId <- g[g%in%cId]
    sortedCId <- as.character(sort(cId)) #-- order gene
    if(length(which(is.na(match(sortedCId, refCId))))>0) {
    	stop("[getInts] error  in the length Id extraction\n")
    	}
	trackB <- listOfGbyTrack[[i]][sortedCId]

    chrIntTrackB <- as.numeric(unlist(trackB))+1
    #checkEquals(length(chrIntTrackAdenom), length(chrIntTrackB), tolerance = 0.001 )


    # test
    if(length(chrIntTrackAdenom)!=length(chrIntTrackB)){
    	obj[[iter]] <- paste("[getInt] length of SW differs in ",chr,"\n" ,sep="")
    	cat(obj[[iter]],"\n")
    	break
    	#return(obj)
    	}else {
#checkEquals(length(chrIntTrackAdenom), length(chrIntTrackB))
		obj[[iter]] <- chrIntTrackB
	}
    iter <- iter+1
  }


  if(length(names(obj))==length(c("gene","chr", wiggstoSearch[indexDenom], wiggstoSearch[vectIndex]))){
  	names(obj) <- c("gene","chr", wiggstoSearch[indexDenom], wiggstoSearch[vectIndex])
  }

  cat("end", chr, "\n")
  return(obj)

}




#- getNorm, return the intensity num, denom and num normalized.
#           takes in argument,
#           - objInt, List, object containing the intensity num and denominator
#           - the  method of normalization.
#           - lowessSpan, integer, the smoother span to lowess
#           - mva, String, path and name to plot the mva graph ("should finish by .png)")
getNorm <- function(objInt, norm=c("linear","lowess", "quantile"),
                    lowessSpan=NULL, mva=NULL){


  #denom <- objInt$denominator
  #num <- objInt$numerator
  cat("start",objInt[[2]],"\n")
  intList <- objInt[3:length(objInt)]
  if(norm=="linear"){
  	denom <- intList[[1]]
    num <- intList[[2]]
    norm <- getNormRatio(num, denom)
    intensityNorm <- norm$intNorm # non log2
  }else if (norm=="lowess"){
    object <- lowessNorm(intList, ff=lowessSpan, plotMVA=mva)

    objInt <- objInt[1:2]
    objInt$lowessMat <- object


   }else if (norm=="quantile"){
    intensityNorm=NULL
    object <- quantileNorm(intList, plotMVA=mva)
    objInt <- objInt[1:2]
    objInt$quantMat <- object
  }else {
    cat("norm should be quantile, lowess, linear\n")
  }
  cat("end",objInt[[2]],"\n")
  return(objInt)
}




#- getResult, return a tsv file containing the current intensities value and ratio
#             takes in argument:
#             - objInt, List, the object containing a list of intensity
#             - norm, normalzation methods
#             - outIntTable, boolean, if TRUE save
#             - initial Directory
getResult <- function(objInt, iter="", outIntTable=FALSE, out, norm="quantile",                         stretchVal, tssLength, denoString){

  cat("start",objInt$chr,"\n")
  chrCurr <- objInt$chr
  name <- objInt$gene
  if(norm=="quantile"){
  	intdf <- objInt$quantMat
  	}else if(norm=="lowess") {
  		intdf <- objInt$lowessMat
  		}else {

  denom <- objInt$denominator
  num <- objInt$numerator
  numNorm <- objInt$intNumNorm
  checkEquals(length(name), length(denom))
  intdf <- base::cbind(as.numeric(denom),  as.numeric(num), as.numeric(numNorm))
  #if(length(objInt)==6){
  #  intdf <- base::cbind(intdf, objInt$intDenoNorm)
  #}
  }

  #- generate SW id
  unic <- unique(name)
  l <- split(name, name)
  checkEquals(length(unic), length(l))#test

  vectorSW <- NULL
  for ( i in 1:length(unic)){
    ntime <- length(l[[unic[i]]])
    vectorSW <- c(vectorSW, 1:ntime)
  }

  #checkEquals(length(denom), length(vectorSW))#test

  outputFrame <- data.frame(base::cbind(name, as.numeric(vectorSW),intdf),
                            row.names=NULL,
                            stringsAsFactors=FALSE)
  #if (length(objInt)==6){
  #  names(outputFrame) <- c("RefSeq Gene", "SW", "Denominator","Nominator","Nominator.Norm", "Denominator.Norm")
  #}else if (length(objInt)==5){
  #  names(outputFrame) <- c("RefSeq Gene", "SW", "Denominator","Nominator","Nominator.Norm")
  #}else {
  #	names(outputFrame) <- c("RefSeq Gene", "SW", "Denominator","Nominator")
  #	}

  if(norm=="quantile"){
  	  names(outputFrame) <- c("RefSeq Gene", "SW", names(intdf))
  	}else if (length(objInt)==5) {
  		names(outputFrame) <- c("RefSeq Gene", "SW", "Denominator","Nominator","Nominator.Norm")
  		}
  		else {
  			names(outputFrame) <- c("RefSeq Gene", "SW", "Denominator","Nominator") 	}


  outChr<-paste(out,"/",chrCurr,"/",sep="")
  timecurr <- gsub(" " ,"",Sys.time())
  timecurr <- gsub("[-:]","",timecurr)
  if (outIntTable==TRUE){
  outfile <- paste(outChr,"IntTable","_", denoString,"_",timecurr,".txt",sep="")
  write.table(outputFrame, outfile,
                sep="\t",row.names=FALSE,quote=FALSE)
    cat(outfile,"\n")
  }

  #return(outputFrame)

  objRes <- list()
  objRes$chr <- chrCurr
  objRes$matrix <- getOutput(inputTable=outputFrame, norm=norm, stretch=stretchVal, tssL=tssLength, chr=chrCurr, outputdir=outChr, denoString=denoString, time=timecurr)

  return(matrix)
  cat("end",chrCurr,"\n")
}











# getCodReg : extract coding region from Wigg and extract intensity values (moyenne).
#             Takes in input:
#             - pathname to a wigg file
#             - vector, containing the number of pb before and after the coding region
#             - pathname to athe organism data base file
#             - span, integer, if NULL take the span of the WIGG else take the specified span
getCodReg <- function(wiggFile,prepost,window, dbCurr,chrCurr,span=NULL){

  #get chr Curr


  #tryCatch({
    #browser()
    irDbCurr <- IRanges(dbCurr[,2], dbCurr[,3])
    wigg <- read.table(wiggFile, header=FALSE,skip=2,colClasses=c("numeric","integer"))
    wigg <- as.matrix(wigg)
    spanWigg <- diff(wigg[1:2,1])
    cat("span in wigg =", spanWigg, "span =", span, "window =", window, "preStamp =",prepost[1], "postStamp =",prepost[2],"\n")

    span <- ifelse(is.null(span), spanWigg,span)

    irWigg <- IRanges(wigg[,1], wigg[,1]+span)

    #ov <- findOverlaps(irDbCurr, irWigg)
    ov <- findOverlaps( irWigg,irDbCurr)
    ov <- as.matrix(ov)

    #- split wigg by G
    uniqueG <- unique(ov[,2])
    swv <- lapply(uniqueG, getSWV, wiggTable=wigg,
                  commonTable=ov,
                  dbcurr=dbCurr,
                  window=window,
                  span=span)

                                        #swvlg <- sort(unlist(lapply(swv,length)))
    names(swv) <- uniqueG
  #}, error = function(ex) {
  #    cat("An error was detected in *getCodReg* :\n",wiggFile,"\n");
  #    print(ex)
  #  }
  # )
  return(swv)
}



#getSWV : get the Sliding Windows values by G, values are obtained using the median of the tags counts
#         takes in input:
#         - gID, numeric, the Gene ID
#         - wiggTable, the entire wiggTable
#         - commonTable, the results matrix of 'overlap' function
#         - dbCurr, the entire database
#         - window,
#         - span, integer, if NULL take the span of the WIGG else take the specified span
getSWV<- function(gID, wiggTable, commonTable,
                  dbcurr, window=250,span){
  #cat("span=",span,"\n")
  tabCurr <- wiggTable[commonTable[which(commonTable[,2]==gID),1],,drop=F]
  gSpec <- dbcurr[gID,,drop=F] # specific gene
  gSize <- as.numeric(gSpec[3])-as.numeric(gSpec[2]) # gene width
  makeGridSum <- function(x){
    if (length(x) == 0 ){
      x <- 0
    }else if (length(x)>1){
      x <- mean(x)
    }
    return(x)
  }




  if(round(gSize/span) > 1){
    gridded <- split(tabCurr[,2],cut(tabCurr[,1],breaks=round(gSize/span),#modif
                                     include.lowest=FALSE))
  }else {
    gridded <- list(tabCurr[,2])
    tmp <-makeGridSum(gridded[[1]])
    vectInt <-tmp
    rm(tmp)
    return(vectInt)
  }


 tmp <-lapply(gridded, makeGridSum)

  #- TSS Study
  #int <- rollapply(as.zoo(unlist(tmp))[1:6],(window/span),mean) # TSS 3000pb; span =500
  #int <- rollapply(zoo::as.zoo(unlist(tmp)),(window/span),mean)
  int <- rollapply(unlist(tmp),(window/span),mean)
  rm(tmp)

  vectInt <- as.vector(int)
  if(as.numeric(gSpec[5])==2){ #gene orientation
    vectInt <- rev(as.vector(int))
  }
  return(vectInt)#vector
}

#-- Normalisation methods


#- Linear Methods
#getNormRatio <- function(trackCurr, trackDenom){
#  ratioList <- as.numeric(trackCurr)/as.numeric(trackDenom)
#  medfact <- median(ratioList,na.rm=TRUE)
#  cat(medfact)
#  ratioListNorm <- ratioList/medfact
#  return(ratioListNorm)
#}

## getNormSumRatio <- function(trackCurr, trackDenom){
##   ## obj <- list()
##   sumfact <- sum(as.numeric(trackCurr))/ sum(as.numeric(trackDenom))
##   ratioList <- as.numeric(trackCurr)/as.numeric(trackDenom)
##   ratioListNorm <- ratioList/sumfact
##   cat(sumfact)
##   obj$ratioNorm <- ratioListNorm
##   obj$intNorm <- as.numeric(trackCurr)/sumfact
##   return(obj)
## }


#- Local Methods : Loess (M/A plot)
lowessNorm <- function (listInt, ff=2/3, plotMVA=NULL){
  out <- list()
  #na.point <- (1:length(xDeno))[!is.na(xDeno) & !is.na(yNum)]
  #xDeno <- xDeno[na.point]
  #yNum <- yNum[na.point]

  xDeno <- as.numeric(listInt[[1]])

  mat <- matrix(nrow= length(listInt[[1]]), ncol =length(listInt))
  for (i in 1:length(listInt)){
  	mat[,i] <- as.numeric(listInt[[i]])
  }
  #names(mat) <- names(listInt)



  #create a norm matrix
  matNorm <- matrix(nrow= length(listInt[[1]]), ncol =length(listInt))
  matNorm[,1] <- xDeno
  #posDeno <- 1
  #for (i in 1:(ncol(mat)-1)){
  for (i in 2:ncol(mat)){
  	yNum <- as.numeric(mat[, i])
  	Mi <-  log2(yNum/xDeno)
    Ai <- log2(sqrt(yNum* xDeno))
  	fit <- lowess(Ai,Mi,f=2/3)
    diff.fit <- approx(fit$x, fit$y, Ai, ties = mean)
    Mnorm <- Mi - diff.fit$y
	ynorm <- yNum*2^(-diff.fit$y)# no log
  	matNorm[,i] <- ynorm
  	}


	#plot
  if(plotMVA==TRUE){
      #InitDir <- getwd()
  	#dirToPrint <- paste(InitDir,"/",grep("chr", strsplit(names(listInt)[1], c("[_\\.]"))[[1]], value=TRUE),sep="")
	#filename <- paste(dirToPrint,"/mvaPlot_lowess.png", sep="")
  	#png(filename)
      X11()
	getmvaplot( as.data.frame(mat), as.data.frame(matNorm), normMeth="lowess")
	#cat("[MVA] ",filename,"\n")
	dev.off()
  	}


  nmdfNorm <- paste("Norm_",names(listInt),sep="")
  matTot <- base::cbind(as.data.frame(mat),as.data.frame(matNorm))
  names(matTot) <- c(names(listInt), nmdfNorm)

  return(matTot)

}




#- quantileNorm, return numerator and denominator normalized according to quantile norm
quantileNorm <- function (listInt, plotMVA=FALSE){
  out <- list()
  #- create matrix
  mat <- matrix(nrow= length(listInt[[1]]), ncol =length(listInt))
  for (i in 1:length(listInt)){
  	mat[,i] <- as.numeric(listInt[[i]])
  	}
  	names(mat) <- names(listInt)
  #mat[,1] <- as.numeric(yNum)
  #mat[,2] <- as.numeric(xDeno)
  #- apply quantile norm
  matnorm <- normalizeQuantiles(mat)
  dfmatnorm <- as.data.frame(matnorm)

  # plot
  if(plotMVA==TRUE){
        #InitDir <- getwd()
	#dirToPrint <- paste(InitDir,"/",grep("chr", strsplit(names(listInt)[1], c("[_\\.]"))[[1]], value=TRUE),sep="")
	#filename <- paste(dirToPrint,"/mvaPlot_quantile.png", sep="")
	#png(filename)
      X11()
	getmvaplot( mat, dfmatnorm, normMeth="quantile")
	#cat("[MVA] ",filename,"\n")
	dev.off()
	}

  nmdfNorm <- paste("Norm_",names(listInt),sep="")
  matTot <- base::cbind(as.data.frame(mat),dfmatnorm)
  names(matTot) <- c(names(listInt), nmdfNorm)
  #out$intNumNorm <- matnorm[,1]
  #out$intDenoNorm <- matnorm[,2]
  return(matTot)
}




#- generateWig: generate a wigg file, takes as input :
#     - dbmat, matrix,the datamatrix
#     - countdf, data.frame, output data frame,
#     - window, the window used
#     - span, int, the span
#     - prepost, c(int, int), the number of pb before and after the gene
#     - col, int, index of the count column
generateWig <- function (countdf ,col,chr, dbori, window, span,prespost=c(1000,1000),out) {
  header <- paste("track type=wiggle_0 name=\"Polyphemus_counts_after_smoothing\" description=\"smoothing window=250 for every 50 bp\"
variableStep chrom=",chr," span=",span,sep="")
  write(header,out)
  #- get the vector of start position
  uniG <- unique(countdf[,1])
  #browser()

  #aa <- unlist(lapply(uniG,getRowNumber,countdf=countdf,dbori=dbori,col=col,
  #               window=window,span=span,prespost=prespost))
  lout <- lapply(uniG, getIncrByG,
         countdf=countdf,dbori=dbori,col=col,
                 window=window,span=span,prespost=prespost)
  resWigg <- do.call("rbind", lout)
  write.table(resWigg,out,append=TRUE,sep="\t",quote=FALSE, row.names=FALSE,col.names=FALSE )
  cat(out, "\n")
  #return(resWigg)
}

## getRowNumber <- function(gCurr,col, countdf, dbori, window,span,prespost){
##   gCurr <- as.character(gCurr)
##   voriG <- as.character(dbori[,4])
##   gRefCurr <- dbori[which(voriG==gCurr),]
##   subV <- as.numeric(countdf[countdf[,1]==gCurr,col]) #matrix
##   return(length(subV))
## }

getIncrByG <- function(gCurr,col, countdf, dbori, window,span,prespost){
  gCurr <- as.character(gCurr)
  voriG <- as.character(dbori[,4])
  gRefCurr <- dbori[which(voriG==gCurr),]
  subV <- as.numeric(countdf[countdf[,1]==gCurr,col]) #matrix
  tmpStartPos <-gRefCurr[,2]
  if (length(tmpStartPos )>1){
    startPos <- min(tmpStartPos)
  }else {startPos <- tmpStartPos}
  startPos <- startPos-prespost[1]
  tmpEndPos <- gRefCurr[,3]
  if (length(tmpEndPos)>1){
    endPos <- max(tmpEndPos)
  }else {endPos <- tmpEndPos}
  endPos <- endPos+prespost[2]
  rm(tmpStartPos,tmpEndPos)
  width <- endPos-startPos
  #browser()
  #endPosReel <- endPos-window+span
  #seqPos <- seq(startPos,endPosReel, span)
  seqq <- seq(startPos,endPos, span)
  seqPos <- seqq[1:(length(seqq)-4)]
  if (length(seqPos)==length(subV)){
    out <- cbind(format(seqPos,scientific=FALSE),subV)
  }else if (length(seqPos)>length(subV)){
    out <- cbind(format(seqPos[1:length(subV)],scientific=FALSE),subV)
    cat("Pos Inf Val:",length(subV)-length(seqPos) ,"\n")

  }else if (length(seqPos)<length(subV)) {
    #browser()
    cat("Pos Sup Val:",length(subV)-length(seqPos) ,"\n")
    out <- cbind(seqPos,subV[1:length(seqPos)])
  }

  return(out)
}



#medichiMerge, merge all medichi files for one chromosome.
#              return the list of peaks and their position. takes in argument:
#              files: vector of String, contains the path and name of each files
#
medichiMerge <- function(files, windo,
                         cutoff=0.01, write=NA){
  tabList <- lapply(files, read.table, header=FALSE,skip=1)
  t1 <- tabList[[1]]
  #- filter
  t1 <- t1[t1[,4]<=cutoff,]
  #nrow(t1)#ok
  ir <- IRanges((t1[,2]-(windo/2)),(t1[,2]+(windo/2)))

  for (i in 2:length(tabList)){
    tcurr<- tabList[[i]]
    tcurr <- tcurr[tcurr[,4]<=cutoff,]

    irCurr <- IRanges((tcurr[,2]-(windo/2)), (tcurr[,2]+(windo/2)))
    ir <- union(ir,irCurr)
  }

  tFinal <- as.data.frame(ir)
  if(!is.na(write)){
    write(paste("windows=",windo,", cutoff=",cutoff,", mean width after merging=",mean(tFinal[,3]),"\n"),
          write)
    write.table(tFinal,write,sep="\t",append=TRUE,row.names=FALSE,quote=FALSE,col.names=F)
    cat(write,"\n")
  }
  return(tFinal)

}

# createFolderHierarchy, create the folder heierarchy for input files

createFolderHierarchy <- function(initialDir){

	#get the chr.
	files <- dir(initialDir, full.names=FALSE)
	filesfull <- dir(initialDir, full.names=TRUE)
	filesplit <- strsplit(files, c("[_\\.]"))
	#store chr num
	chr <- grep("chr", unlist(filesplit), value = TRUE)
	chr <- unique(chr)

	#create hierarchy
	dirH <- function(x) {

		match1 <- which(unlist(lapply(filesplit,function(t){ifelse(length(which(grepl(paste("^",x,"$",sep=""),t))), return(TRUE),return(FALSE))})))

		mainfolder <- paste(initialDir,"/",x, sep ="")
		subfolder1 <- paste(mainfolder,"/","wig", sep = "")
		subfolder2 <- paste(mainfolder,"/","pos", sep = "")
		lapply(c(mainfolder, subfolder1, subfolder2 ), dir.create)

		for (i in match1){

			if(filesplit[[i]][length(filesplit[[i]])]=="wig"){
				system(paste("mv ", filesfull[i], " ",subfolder1,sep="" ))
				}else {
					system(paste("mv ", filesfull[i], " ", subfolder2,sep=""))
					}
				}
	}#end function

	logg <- lapply (chr, dirH)

	if(is.null(unlist(logg))){
		cat("Successfuly create folder hierarchy.\n")
		}
	}
#createFolderHierarchy("/Users/admin/Documents/polyphemus_Paper/polyphemus_Code/workingSpace")





#- getLowessFitVal, return a list of table where it was extracted the fit values from a lowess regression, takes in input
#           - id, geneId from the final output of polyphemus
#           - splited by id the final output <=> list of id table
getLowessFitVal <- function(id,list, toPlot = FALSE){
  tabCurr <- list[[id]]
  y <- as.numeric(tabCurr[,2])
  x <- 1:length(y)
  fit <- lowess(x,y, f=0.6)
  diffFit <- approx(fit$x, fit$y, x, ties=mean)
  if(!is.null(toPlot)){
    lines(diffFit$x, diffFit$y, col = rgb(255,0,0,50,maxColorValue=255))
  }
  tab <- base::cbind(tabCurr, as.numeric(diffFit$y))
  return(tab)
}



#- getOutput, data Preprocessing function, return a list of table with
#               3 columns (ID, ratio, approx).It takes in input
#          - filename, String, the path and name of polyphemus output
#          - interpol, String, function use for interpolation (lowess (default), bspline etc..)
#          - colnum, integer, index of the column to use as numerator of the ratio
#          - coldeno, integer, index of the column to use as denominator of the ratio.
#          - plot, Logical, whether the results should be plot.
#          - chr, String, which Chr?
#          - inputTable, results table used (NULL if filename is used)
#          - tssLength, Integer, length of the tss (=number of span to define the tss)
	getOutput <- function (inputTable=NULL, filename=NULL, norm, stretch, tssL, chr, outputdir, denoString, time, toPlot=NULL){

if(!is.null(filename)){
  df <- read.table(filename, header=TRUE,as.is=TRUE,sep="\t" )
  }else {
  	df <- inputTable
  	}
  gid <- unique(df[,1])

  #browser()
  #folder curr
 folderCurr <- paste(InitDir,"/",chr,"/",sep="")

 if(!is.null(norm)){
	if(norm=="quantile"){
		# find out colDeno
		coldeno  <- max(grep( denoString, names(df),value=FALSE))
		#colnum <- 5
		#coldeno <- 6
	}else if (norm=="lowess"){
		coldeno  <- max(grep( denoString, names(df),value=FALSE))
		#colnum <- 5
		#coldeno <- 3
		}
	}else {
		coldeno  <- max(grep( denoString, names(df),value=FALSE))

		#colnum <- 4
		#coldeno <- 3
		}


for (i in (coldeno+1):ncol(df)){


  colnum <- i

  #- get string num

  numString <- strsplit(names(df)[i], "_")[[1]][2]


  logratio <- log2(as.numeric(df[,colnum])/as.numeric(df[,coldeno]))
  ratiodf <- R.utils::dataFrame(colClasses=c(ID="character", LogRatio="numeric"),
                     nrow=length(logratio))
  #max(unlist(lapply(split(as.numeric(df[,2]),df[,1]),length)))# 224357
  ratiodf[,1] <-  df[,1]
  ratiodf[,2] <-  logratio
  ratiolist <- split(ratiodf, ratiodf[,1])

  xinitval <- 400
  #- filter gene if length is inferior to tsslength
  del <- which(as.numeric(unlist(lapply(ratiolist, nrow))) < (tssL +1))
  if(length(del)>0){
    ratiolist <- ratiolist[-del]
  }
  gid <- names(ratiolist)

    lFittedTab <- lapply(gid, getLowessFitVal,list= ratiolist, toPlot=toPlot)

  #ratiodf2 <- do.call("rbind", lFittedTab)
  if(!is.null(toPlot)){
    dev.off()
  }
  names(lFittedTab) <- gid

  if(!is.null(stretch)){
    mid <- median(unlist(lapply(lFittedTab, nrow)))
    if(stretch != "median"){
      cat("[stretch] default stretch value : ",ceiling(mid)," your stretch value :", ceiling(stretch),".\n")
      mid<-stretch
    } else {
	cat("[stretch] default stretch value : ",ceiling(mid),".\n")
    }
    lFittedTab <- lapply(lFittedTab, stretchGene, tssCut=tssL,
                         lengthRef=ceiling(mid))
 }
  #return(lFittedTab)

  #ratiotable <- data.frame(do.call("rbind", lFittedTab), row.names=NULL)
ratiotable <- NULL
  for(i in 1:length(lFittedTab)){
  	ratiotable <-base::rbind(ratiotable, lFittedTab[[i]])
  	}
  gid <- unique(ratiotable[,1])

  widthmax <- max(unlist(lapply(split(as.numeric(ratiotable[,2]),ratiotable[,1]),length))) # in term of SW
  # initialize
  fakedf <- matrix(nrow=length(gid),
                   ncol=widthmax)

  for (i in 1:length(gid)){

    ratios <- as.numeric(ratiotable[which(ratiotable[,1]==gid[i]),2]) #ou 3
    fakedf[i,1:length(ratios)] <- ratios

  }

  #save output
  #get output directory

  row.names(fakedf) <- gid
  filename <- paste(outputdir,"matrix_", denoString,"VS",numString,"_",time,".txt",sep="")
  write.table(fakedf, filename,
                sep="\t",quote=FALSE)
  cat(filename,"\n")

  }

  return(fakedf)
}





#- stretchGene, return a table with stretched body of the profile, take in input
#               - genetab, the gene table (id, logratio, fittedval)
#               - tssCut, Integer,number of line considered as part of the TSS region
#               - lengthRef, Numeric, the length of the reference gene
stretchGene <- function(genetab, tssCut, lengthRef){

  tab <- genetab[-c(1:tssCut),]

  x <- 0:(nrow(tab)) # normally -1 but to obtain an interpolate value for the last element we put +1 more.
  y <- c(tab[,2], tab[nrow(tab),2]) # concatenate the last element

  lRefReal <- lengthRef-tssCut

  corf <- lRefReal/ (nrow(tab))
  X <- corf*x

  fit <- lowess(X, y, f=2/3)
  diff.fit <- approx(fit$x, fit$y, 1:lRefReal, ties=mean)

  newVal <-  diff.fit$y
  newtab <- data.frame(rep(as.character(tab[1,1]), lengthRef), c(genetab[1:tssCut,3],diff.fit$y))
  return(newtab)
}





#getmvaplot: Call by by quantilenorm function display mva plot

getmvaplot <- function(matRaw, matNorm, normMeth=NULL){

    if(!is.null(normMeth)){
        normMethod <- normMeth
    }
	checkEquals(ncol(matRaw), ncol(matNorm))
	posDeno <- 1
	par(mfrow=c((ncol(matRaw)-1),2))
	for (i in 1:(ncol(matRaw)-1)){
		Mbrut <- log2(matRaw[,posDeno+i]/matRaw[,posDeno])
		Abrut <- log2(sqrt(matRaw[,posDeno+i]*matRaw[,posDeno]))
		fitbrut <- lowess(Abrut, Mbrut, f=2/3)
		M <- log2(matNorm[,posDeno+i]/matNorm[,posDeno])
		A <- log2(sqrt(matNorm[,posDeno+i]* matNorm[,posDeno]))
		fit <- lowess(A, M, f=2/3)

		ymax <- ceiling(max(c(max(Mbrut), max(M))))
		xmax <- ceiling(max(c(max(Abrut), max(A))))

		plot(Abrut,Mbrut, col="green", main="raw", ylim=c(-ymax,ymax), xlim=c(0,xmax))
		lines(fitbrut, col = "red")
		abline(h=0,col="grey", lty=2, lwd = 2)

		plot(A,M, col="green", main=paste(normMethod,"normalization"),ylim=c(-ymax,ymax), xlim=c(0,xmax))
		lines(fit, col = "red")
		abline(h=0,col="grey", lty=2, lwd=2)

		}

	}
