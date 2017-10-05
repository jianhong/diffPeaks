#' @title Counting reads in features
#' @description Counting reads in giving features.
#' @param features an object of \link[GenomicRanges]{GRanges}
#' @param bamfiles filenames of aligned reads
#' @param samples names of samples, could be same length with bamfiles
#' @param windowSize the size of windows for counts
#' @param colData a DataFrame or data.frame with at least a single column. 
#' Rows of colData correspond to bamfiles
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
                       windowSize=200L, ...){
  stopifnot(inherits(features, "GRanges"))
  stopifnot(is.character(samples))
  if(missing(colData)){
    stop("colData is required")
  }
  stopifnot(nrow(colData)==length(bamfiles))
  stopifnot(inherits(colData, "DataFrame"))
  ## count for DESeq2
	so <- summarizeOverlaps(features, bamfiles, ...)
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
	so2 <- summarizeOverlaps(features.tile, bamfiles, ...)
	colData(so) <- colData(so2) <- colData
	list(feature=so, tile.feature=so2, signature="countTable")
}

