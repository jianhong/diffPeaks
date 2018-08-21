#' Count reads overlapping genomic ranges
#'
#' Count reads overlapping a set of genimc features represented as
#' genomic ranges. This function does not work for parallel.
#' @param features A object of \link[GenomicRanges:GRanges-class]{GRanges} representing the
#' feature regions to be counted.
#' @param reads An object that represents the data to be counted. See
#' \link[GenomicAlignments:summarizeOverlaps-methods]{summarizeOverlaps}. If reads are more than 1 bam files,
#' it should be a vector of character with full path, otherwise current working directory 
#' is the default directory. For paired end reads, 
#' @param ignore.strand logical(1). ignore strand?
#' @param inter.feature not used. This parameter is required by
#' \link[GenomicAlignments:summarizeOverlaps-methods]{summarizeOverlaps}.
#' @export
#' @return return a summarized experiment object with chromosome-level depth
#' information for each input sample as metadata.
#'
IntersectionNotStrict <-function(features,
                                 reads,
                                 ignore.strand = TRUE,
                                 inter.feature = FALSE) 
{
  ## NOT work for parallel
  ov <- findOverlaps(reads,
                     features,
                     type = "within",
                     ignore.strand = ignore.strand)
  countSubjectHits(ov)
}
