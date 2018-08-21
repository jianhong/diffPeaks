#' @title Counting reads in features
#' @description Counting reads in giving features.
#' @param features an object of \link[GenomicRanges]{GRanges}
#' @param bamfiles filenames of aligned reads
#' @param samples names of samples, could be same length with bamfiles
#' @param windowSize the size of windows for counts
#' @param colData a DataFrame or data.frame with at least a single column. 
#' Rows of colData correspond to bamfiles
#' @param mode mode of counts. See \link[GenomicAlignments]{summarizeOverlaps} for details. 
#' Default is \link{IntersectionNotStrict}.
#' @param ... parameters could be passed to \link[GenomicAlignments]{summarizeOverlaps}
#' @return a list of \link[SummarizedExperiment]{RangedSummarizedExperiment}
#' @import SummarizedExperiment
#' @import GenomicAlignments
#' @import GenomicRanges
#' @import S4Vectors
#' @export
#' @author Jianhong Ou
#' @examples 
#' path <- system.file("extdata", package = "diffPeaks", mustWork = TRUE)
#' bamfiles <- dir(path, "bam$")
#' peaks <- dir(path, "bed$")
#' p <- mergePeaks(file.path(path, peaks))
#' colData <- DataFrame(samples=bamfiles, condition=sub(".rep..bam", "", bamfiles))
#' cnt <- countTable(p, file.path(path, bamfiles), colData)

countTable <- function(features, bamfiles, colData,
                       samples=bamfiles, 
                       windowSize=200L, 
                       mode=IntersectionNotStrict,
                       ...){
  stopifnot(inherits(features, "GRanges"))
  stopifnot(is.character(samples))
  if(missing(colData)){
    stop("colData is required")
  }
  stopifnot(nrow(colData)==length(bamfiles))
  stopifnot(inherits(colData, "DataFrame"))
  ## count for DESeq2
	so <- summarizeOverlaps(features, bamfiles, mode=mode, ...)
	## count for splitted two features.
	## split each feature into two parts, in next step, will test it using fisher's exact test.
	features2 <- features
	features2 <- tile(features2, n = 2)
	features2 <- unlist(features2)
	features2$feature_oid <- rep(seq_along(features), each=2)
	so2 <- summarizeOverlaps(features2, bamfiles, mode = mode, ...)
	## count for splited features
	features$feature_oid <- seq_along(features)
	fw <- ceiling(width(features)/windowSize) + 2 ## set minimal width == 3*wind
	## recenter the features
	f.center <- floor((start(features) + end(features))/2)
	f.start <- f.center - ceiling(windowSize * fw/2)
	f.start[f.start<1] <- 1
	start(features) <- f.start
	width(features) <- fw * windowSize
	features.tile <- tile(features, width = windowSize)
	features.tile.l <- lengths(features.tile)
	features.tile <- unlist(features.tile)
	mcols(features.tile) <- mcols(features)[rep(seq_along(features), features.tile.l), ]
	so3 <- summarizeOverlaps(features.tile, bamfiles, mode=mode, ...)
	colData(so) <- colData(so2) <- colData(so3) <- colData
	list(feature=so, two.feature=so2, tile.feature=so3, signature="countTable")
}

#' is output of countTable
#' @description check an object is output of \link{countTable} or not.
#' @param counts the object to be checked
#' @param error the error message
#' 
isCountTable <- function(counts, error="counts must be output of countTable!"){
  if(length(counts$signature)!=1){
    stop(error)
  }
  if(counts$signature!="countTable"){
    stop(error)
  }
  if(any(names(counts)!=c("feature", "two.feature", "tile.feature", "signature"))){
    stop(error)
  }
  invisible(return(TRUE))
}
