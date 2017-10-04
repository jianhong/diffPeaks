#' @title differential binding analysis
#' @description use DEseq2 to test differential peaks
#' @param counts output of \link{countTable}
#' @param ... parameters could be passed to \link[DESeq2]{DESeqDataSet}
#' @return an object of \link[GenomicRanges]{GRanges}
#' @import DESeq2
#' @export
#' @author Jianhong Ou
#' @examples 
#' path <- system.file("extdata", package = "diffPeaks", mustWork = TRUE)
#' bamfiles <- dir(path, "bam$")
#' peaks <- dir(path, "bed$")
#' p <- mergePeaks(file.path(path, peaks))
#' colData <- DataFrame(samples=bamfiles, condition=sub(".rep..bam", "", bamfiles))
#' cnt <- countTable(p, file.path(path, bamfiles), colData)
#' DEpeaks(cnt, design= ~condition)

DEpeaks <- function(counts, ...){
  if(length(counts$signature)!=1){
    stop("counts must be output of countTable!")
  }
  if(counts$signature!="countTable"){
    stop("counts must be output of countTable!")
  }
  if(any(names(counts)!=c("feature", "tile.feature", "signature"))){
    stop("counts must be output of countTable!")
  }
  dds <- DESeqDataSet(counts$feature, ...)
  dds <- DESeq(dds, betaPrior = FALSE)
  if(length(resultsNames(dds))==2 && resultsNames(dds)[1]=="Intercept"){
    res <- results(dds)
    resLFC <- lfcShrink(dds, coef = 2, res = res)
    gr <- rowRanges(counts$feature)
    mcols(gr) <- resLFC
    gr
  }else{
    stop("Only simple comparison is supported.")
  }
}