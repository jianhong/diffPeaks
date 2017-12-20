#' @title Merge BED files
#' @description Merge multiple peak files
#' @param beds BED filenames
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

mergePeaks <- function(beds, maxPeakWidth=5000, ...){
  stopifnot(length(beds)>1)
	d <- lapply(beds, function(.ele) import(.ele, ...))
	d <- sort(reduce(unlist(GRangesList(d))))
	w <- width(d)
	if(any(w>maxPeakWidth)){
	  bins <- ceiling(w/maxPeakWidth)
	  d <- tile(d, n=bins)
	}
	d
}
