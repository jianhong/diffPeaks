#' @title differential binding analysis by Cochran-Mantel-Haenszel Test
#' @description use Cochran-Mantel-Haenszel Test to test differential peaks
#' @param counts output of \link{countTable}
#' @param conditionA,conditionB condition A and B. Comparison will be A-B.
#' @param ... not used.
#' @return an object of \link[GenomicRanges]{GRanges}
#' @import GenomicAlignments
#' @importFrom stats mantelhaen.test
#' @export
#' @author Jianhong Ou
#' @examples 
#' path <- system.file("extdata", package = "diffPeaks", mustWork = TRUE)
#' bamfiles <- dir(path, "bam$")
#' peaks <- dir(path, "bed$")
#' p <- mergePeaks(file.path(path, peaks))
#' colData <- DataFrame(samples=bamfiles, condition=sub(".rep..bam", "", bamfiles), 
#'                      pairs=rep(c(1,2,3), 2))
#' cnt <- countTable(p, file.path(path, bamfiles), colData)
#' CMHpeaks(cnt, conditionA="inj", conditionB="uni")

CMHpeaks <- function(counts, conditionA, conditionB, ...){
  isCountTable(counts)
  rr <- rowRanges(counts$two.feature)
  count <- assays(counts$two.feature)$counts
  colData <- colData(counts$two.feature)
  if(!"pairs" %in% colnames(colData)){
    stop("pairs must be a column of colData")
  }
  stopifnot(conditionA %in% colData$condition)
  stopifnot(conditionB %in% colData$condition)
  cntA <- count[, colData$condition %in% conditionA, drop=FALSE]
  cntB <- count[, colData$condition %in% conditionB, drop=FALSE]
  if(ncol(cntA)!=ncol(cntB)){
    stop("sample number in each condition are not identical.")
  }
  colData <- rbind(colData[colData$condition %in% conditionA, ], 
                   colData[colData$condition %in% conditionB, ])
  count <- data.frame(cbind(cntA, cntB))
  count <- split(count, as.character(rr$feature_oid))
  pval <- sapply(count, function(.ele) {
    .ele <- split(data.frame(t(.ele)), colData$pairs)
    .e <- array(0, dim = c(2, 2, length(.ele)))
    for(i in seq_along(.ele)){
      .e[, , i] <- as.matrix(.ele[[i]])
    }
    .e <- .e + 1
    mantelhaen.test(.e)$p.value
  })
  rr1 <- GRanges(seqnames=as.character(rr$feature_oid), 
                 ranges=ranges(rr),
                 strand=strand(rr))
  rr1 <- reduce(rr1)
  rr1 <- GRanges(seqnames = seqnames(rr[match(as.character(seqnames(rr1)), as.character(rr$feature_oid))]),
                 ranges = ranges(rr1), strand = strand(rr1), pvalue=pval[as.character(seqnames(rr1))])
  seqinfo(rr1) <- seqinfo(rr)
  rr1$padj <- p.adjust(rr1$pvalue, method = "BH")
  rr1$type <- "Cochran-Mantel-Haenszel Test"
  rr1
}