#' @title Merge BED files
#' @description Merge multiple peak files
#' @param beds BED filenames
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

mergePeaks <- function(beds, ...){
  stopifnot(length(beds)>1)
	d <- lapply(beds, function(.ele) import(.ele, ...))
	d <- sort(reduce(unlist(GRangesList(d))))
	d
}
