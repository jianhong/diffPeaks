#' @title differential binding analysis by Fisher's exact test
#' @description use Fishter's exact test to test differential peaks
#' @param counts output of \link{countTable}
#' @param conditionA,conditionB condition A and B. Comparison will be A-B.
#' @param ... not used.
#' @return an object of \link[GenomicRanges]{GRanges}
#' @import GenomicAlignments
#' @importFrom stats fisher.test
#' @export
#' @author Jianhong Ou
#' @examples 
#' path <- system.file("extdata", package = "diffPeaks", mustWork = TRUE)
#' bamfiles <- dir(path, "bam$")
#' peaks <- dir(path, "bed$")
#' p <- mergePeaks(file.path(path, peaks))
#' colData <- DataFrame(samples=bamfiles, condition=sub(".rep..bam", "", bamfiles))
#' cnt <- countTable(p, file.path(path, bamfiles), colData)
#' FTpeaks(cnt, conditionA="inj", conditionB="uni")

FTpeaks <- function(counts, conditionA, conditionB, ...){
  isCountTable(counts)
  rr <- rowRanges(counts$two.feature)
  count <- assays(counts$two.feature)$counts
  colData <- colData(counts$two.feature)
  stopifnot(conditionA %in% colData$condition)
  stopifnot(conditionB %in% colData$condition)
  cntA <- count[, colData$condition %in% conditionA, drop=FALSE]
  cntB <- count[, colData$condition %in% conditionB, drop=FALSE]
  count <- data.frame(A=rowSums(cntA), B=rowSums(cntB))
  count <- split(count, as.character(rr$feature_oid))
  pval <- sapply(count, function(.ele) fisher.test(.ele)$p.value)
  rr1 <- GRanges(seqnames=as.character(rr$feature_oid), 
                 ranges=ranges(rr),
                 strand=strand(rr))
  rr1 <- reduce(rr1)
  rr1 <- GRanges(seqnames = seqnames(rr[match(as.character(seqnames(rr1)), as.character(rr$feature_oid))]),
                 ranges = ranges(rr1), strand = strand(rr1), pvalue=pval[as.character(seqnames(rr1))])
  seqinfo(rr1) <- seqinfo(rr)
  rr1$padj <- p.adjust(rr1$pvalue, method = "BH")
  rr1$type <- "Fisher's exact test"
  rr1
}