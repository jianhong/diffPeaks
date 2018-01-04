#' @title Merge BED files
#' @description Merge multiple peak files
#' @param peaks BED filenames which indicate the peaks or an GRangesList 
#' or GRanges object could be used as input peaks.
#' @param maxPeakWidth maximal peak width. If greater than maxPeakWidth,
#'        the peak will be divided in half to fit the parameter.
#' @param ... parameters could be passed to \link[rtracklayer]{import}
#' @return an object of \link[GenomicRanges]{GRanges}
#' @import GenomicRanges
#' @import rtracklayer
#' @export
#' @author Jianhong Ou
#' @examples 
#' path <- system.file("extdata", package = "diffPeaks", mustWork = TRUE)
#' peaks <- dir(path, "bed$")
#' p <- mergePeaks(file.path(path, peaks))

mergePeaks <- function(peaks, maxPeakWidth=5000, ...){
  if(inherits(peaks, c("GRanges", "GRangesList", "list"))){
    if(is(peaks, "list")){
      peaks <- GRangesList(peaks)
    }
    d <- sort(reduce(unlist(peaks)))
  }else{
    stopifnot(length(peaks)>1)
    d <- lapply(peaks, function(.ele) import(.ele, ...))
    d <- sort(reduce(unlist(GRangesList(d))))
  }
	w <- width(d)
	if(any(w>maxPeakWidth)){
	  bins <- ceiling(w/maxPeakWidth)
	  d <- tile(d, n=bins)
	  d <- unlist(d)
	}
	d
}
